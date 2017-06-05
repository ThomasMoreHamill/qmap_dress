import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import sys
import pygrib
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import os
import time as timey

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonth = sys.argv[1]
imonth = cmonths.index(cmonth) + 1
if imonth < 10:
    cmono = '0' + str(imonth)
else:
    cmono = str(imonth)

print 'cmonth, cmono = ',cmonth, cmono

# ---- determine lat/lon information and masks from a sample file
#      conusmask is the strict mask for inside CONUS, practical_mask is the 
#      set of points with CCPA data not missing (1=available, 0 = no)

filename= '/data/thamill/Rf2_tests/ccpa_v1/precip_ccpav1_2002010200_to_2016123100.nc'
print filename
nc = Dataset(filename)
lats_anal = nc.variables['lats_anal'][:]
lons_anal = nc.variables['lons_anal'][:]
xavals = nc.variables['xa'][:]
yavals = nc.variables['ya'][:]
nia = len(xavals)
nja = len(yavals)
conusmask_in = nc.variables['conusmask'][:,:]
iyyyymmddhh_list = nc.variables['yyyymmddhh'][:]
cyyyymmddhh_list = str(iyyyymmddhh_list)

# ---- make a list of 1/0 for whether or not to read in this date in iyyyymmddhh_list
#      for a given month we want data from that month's central date +/- 31 days

ndates = len(iyyyymmddhh_list)
usethisdate = np.zeros((ndates),dtype=np.int16)

for iyear in range(2002, 2017):
    cyyyymmddhh_center = str(iyear) + cmono + '1500'
    cyyyymmddhh_begin = dateshift(cyyyymmddhh_center,-24*31)
    cyyyymmddhh_end = dateshift(cyyyymmddhh_center,24*31)
    print 'flagging dates from ', cyyyymmddhh_begin, ' to ',cyyyymmddhh_end,' to load'
    date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
    print date_list
    for cdate in date_list:
        idate = int(cdate)
        idx = np.where(iyyyymmddhh_list == idate)
        if idx >= 0:
            usethisdate[idx] = 1
            
print 'number of good dates = ',np.sum(usethisdate)

# ---- estimate the climatology from relative frequency

ktr = 0    
climo_POP = np.zeros((nja,nia), dtype=np.float32) 
climo_1mm = np.zeros((nja,nia), dtype=np.float32) 
climo_2p5mm = np.zeros((nja,nia), dtype=np.float32) 
climo_5mm = np.zeros((nja,nia), dtype=np.float32) 
climo_10mm = np.zeros((nja,nia), dtype=np.float32) 
climo_25mm = np.zeros((nja,nia), dtype=np.float32)
climo_50mm = np.zeros((nja,nia), dtype=np.float32) 
ones  = np.ones((nja,nia), dtype=np.float32)
zeros = np.zeros((nja,nia), dtype=np.float32)       
for idate in range(ndates):
    if usethisdate[idate] == 1:
        print 'loading date = ',iyyyymmddhh_list[idate]
        apcp_anal = nc.variables['apcp_anal'][idate,:,:]
        if apcp_anal[nja/2,nia/2] >= 0. :
            work = np.where(apcp_anal > 0.254, ones, zeros)
            climo_POP[:,:] = climo_POP[:,:] + work[:,:]
            work = np.where(apcp_anal > 1.0, ones, zeros)
            climo_1mm[:,:] = climo_1mm[:,:] + work[:,:]
            work = np.where(apcp_anal > 2.5, ones, zeros)
            climo_2p5mm[:,:] = climo_2p5mm[:,:] + work[:,:]
            work = np.where(apcp_anal > 5.0, ones, zeros)
            climo_5mm[:,:] = climo_5mm[:,:] + work[:,:]
            work = np.where(apcp_anal > 10.0, ones, zeros)
            climo_10mm[:,:] = climo_10mm[:,:] + work[:,:]
            work = np.where(apcp_anal > 25.0, ones, zeros)
            climo_25mm[:,:] = climo_25mm[:,:] + work[:,:]
            work = np.where(apcp_anal > 50.0, ones, zeros)
            climo_50mm[:,:] = climo_50mm[:,:] + work[:,:]
            ktr = ktr+1

nc.close()

print 'number of samples: ', ktr

climo_POP = climo_POP / np.float(ktr)
climo_1mm = climo_1mm / np.float(ktr)
climo_2p5mm = climo_2p5mm / np.float(ktr)
climo_5mm = climo_5mm / np.float(ktr)
climo_10mm = climo_10mm / np.float(ktr)
climo_25mm = climo_25mm / np.float(ktr)
climo_50mm = climo_50mm / np.float(ktr)

# ---- set points outside CONUS to -99.99

mninetynine = -99.99*np.ones((nja,nia), dtype=np.float32)
climo_POP = np.where(conusmask_in == 0, mninetynine, climo_POP)
climo_1mm = np.where(conusmask_in == 0, mninetynine, climo_1mm)
climo_2p5mm = np.where(conusmask_in == 0, mninetynine, climo_2p5mm)
climo_5mm = np.where(conusmask_in == 0, mninetynine, climo_5mm)
climo_10mm = np.where(conusmask_in == 0, mninetynine, climo_10mm)
climo_25mm = np.where(conusmask_in == 0, mninetynine, climo_25mm)
climo_50mm = np.where(conusmask_in == 0, mninetynine, climo_50mm)

# ---- open and initialize the netCDF file we'll be writing to.

outfilename = \
    '/data/thamill/Rf2_tests/ccpa_v1/0.125d/apcp_climatologies_12_to_00UTC_'+\
    cmonth+'.nc'
print outfilename
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
    
xa = rootgrp.createDimension('xa',nia)
xav = rootgrp.createVariable('xa','f4',('xa',))
xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
xav.units = "grid index (dimensionless)" 

ya = rootgrp.createDimension('ya',nja)
yav = rootgrp.createVariable('ya','f4',('ya',))
yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
yav.units = "grid index (dimensionless)"
 
lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
lonsa.long_name = "longitude" 
lonsa.units = "degrees_east" 

latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
latsa.long_name = "latitude" 
latsa.units = "degrees_north" 

conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
conusmask.units=""

climo_prob_POP = rootgrp.createVariable('climo_prob_POP','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_POP.units = "" 
climo_prob_POP.long_name = \
    "Climatological probability of exceeding 0.254 mm in accumulation period" 
climo_prob_POP.valid_range = [0.0,1.0]
climo_prob_POP.missing_value = -99.99

climo_prob_1mm = rootgrp.createVariable('climo_prob_1mm','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_1mm.units = "" 
climo_prob_1mm.long_name = \
    "Climatological probability of exceeding 1.0 mm in accumulation period" 
climo_prob_1mm.valid_range = [0.0,1.0]
climo_prob_1mm.missing_value = -99.99

climo_prob_2p5mm = rootgrp.createVariable('climo_prob_2p5mm','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_2p5mm.units = "" 
climo_prob_2p5mm.long_name = \
    "Climatological probability of exceeding 2.5 mm in accumulation period" 
climo_prob_2p5mm.valid_range = [0.0,1.0]
climo_prob_2p5mm.missing_value = -99.99

climo_prob_5mm = rootgrp.createVariable('climo_prob_5mm','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_5mm.units = "" 
climo_prob_5mm.long_name = \
    "Climatological probability of exceeding 5 mm mm in accumulation period" 
climo_prob_5mm.valid_range = [0.0,1.0]
climo_prob_5mm.missing_value = -99.99

climo_prob_10mm = rootgrp.createVariable('climo_prob_10mm','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_10mm.units = "" 
climo_prob_10mm.long_name = \
    "Climatological probability of exceeding 10.0 mm in accumulation period" 
climo_prob_10mm.valid_range = [0.0,1.0]
climo_prob_10mm.missing_value = -99.99

climo_prob_25mm = rootgrp.createVariable('climo_prob_25mm','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_25mm.units = "" 
climo_prob_25mm.long_name = \
    "Climatological probability of exceeding 25.0 mm in accumulation period" 
climo_prob_25mm.valid_range = [0.0,1.0]
climo_prob_25mm.missing_value = -99.99

climo_prob_50mm = rootgrp.createVariable('climo_prob_50mm','f4',('ya','xa',),
    zlib=True,least_significant_digit=4)  
climo_prob_50mm.units = "" 
climo_prob_50mm.long_name = \
    "Climatological probability of exceeding 50.0 mm in accumulation period" 
climo_prob_50mm.valid_range = [0.0,1.0]
climo_prob_50mm.missing_value = -99.99

rootgrp.latcorners = \
    [lats_anal[0,0], lats_anal[0,-1], lats_anal[-1,0], lats_anal[-1,-1]]
rootgrp.loncorners = \
    [lons_anal[0,0], lons_anal[0,-1], lons_anal[-1,0], lons_anal[-1,-1]]

rootgrp.stream = "s4" # ????
rootgrp.title = \
    "12-h climatological probabilities of exceeding various precip thresholds"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created 5 May 2017 by Hamill" 
rootgrp.institution = "ESRL/PSD using CCPA data from NCEP/EMC"
rootgrp.platform = "Model" 
rootgrp.references = "n/a"

# --- now write data to records

climo_prob_POP[:] = climo_POP
climo_prob_1mm[:] = climo_1mm
climo_prob_2p5mm[:] = climo_2p5mm
climo_prob_5mm[:] = climo_5mm
climo_prob_10mm[:] = climo_10mm
climo_prob_25mm[:] = climo_25mm
climo_prob_50mm[:]= climo_50mm

print 'climo_POP[nja/2,0:-1:10] = ',climo_POP[nja/2,0:-1:10]
print 'climo_10mm[nja/2,0:-1:10] = ',climo_10mm[nja/2,0:-1:10]

xav[:]   = np.arange(nia)
yav[:]   = np.arange(nja)
lonsa[:]  = lons_anal
latsa[:]  = lats_anal

conusmask[:] = conusmask_in
rootgrp.close()

