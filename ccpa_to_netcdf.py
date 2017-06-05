import numpy as np
import sys
#import /Users/jwhitaker/python/pygrib.git/pygrib
import pygrib
import os
import time as timey
from netCDF4 import Dataset
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import cPickle

def read_gribdata(idxfilename, gribfilename, validityDate, validityTime):
    grb=[]
    istat = -1
    for iter in range(10):
        fexist_grib = False
        fexist_grib = os.path.exists(gribfilename)
        fexist_idx  = False
        fexist_idx  = os.path.exists(idxfilename)
        if fexist_grib and fexist_idx:
            try:
                fcstfile = pygrib.index(idxfilename)
                grb = fcstfile.select(validityDate=validityDate, \
                    validityTime=validityTime)[0]
                istat = 0
                fcstfile.close()
                return istat,grb
            except (IOError, ValueError, RuntimeError):
                print 'Error reading ', gribfilename
                istat = -1
    return istat,grb



# =====================================================================================

# ---- read in sample precip analysis file in order to get the lats/lons

infile = '/Users/thamill/refcst2/ccpa_conus_0.125d_t12z_06h'
lonsa = None; latsa = None
fexist_apcp = os.path.exists(infile)
if fexist_apcp:
    print infile
    ffcst = pygrib.open(infile)
    grb = ffcst.select(shortName='tp')[0]
    #print 'latsa = ',latsa
    if latsa is None:
        lats_anal, lons_anal = grb.latlons()
        latsa1d = lats_anal[:,0]
        lonsa1d = lons_anal[0,:]
        nlatsa, nlonsa = lons_anal.shape
        coslatsa = np.cos((np.pi/180.)*lats_anal)
        ffcst.close()
        zeros = np.zeros((nlatsa,nlonsa),dtype=np.float)
        ones = np.ones((nlatsa,nlonsa),dtype=np.float)

# ---- initialize

date_begin = '2002010200'
date_end   = '2016123100'
date_list  = daterange(date_begin, date_end, 12)
ndates = len(date_list)

# ---- read in CONUS mask

infile = '/Users/thamill/precip/conusmask_ccpa.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conusmask'][:,:]
nc.close()

# ---- set up netCDF file particulars

outfilename = '/data/thamill/Rf2_tests/ccpa_v1/precip_ccpav1_'+\
    date_begin+'_to_'+date_end+'.nc'

print outfilename
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
xa = rootgrp.createDimension('xa',nlonsa)
xva = rootgrp.createVariable('xa','f4',('xa',))
xva.long_name = "eastward grid point number, precip analysis" 
xva.units = "n/a" 

ya = rootgrp.createDimension('ya',nlatsa)
yva = rootgrp.createVariable('ya','f4',('ya',))
yva.long_name = "northward grid point number, precip analysis" 
yva.units = "n/a" 

time = rootgrp.createDimension('time',None)
timev = rootgrp.createVariable('time','f4',('time',))
timev.units = "hours since 1-1-1 00:00:0.0" 

lonsa = rootgrp.createVariable('lons_anal','f4',('ya','xa',))
lonsa.long_name = "longitude" 
lonsa.units = "degrees_east" 

latsa = rootgrp.createVariable('lats_anal','f4',('ya','xa',))
latsa.long_name = "latitude" 
latsa.units = "degrees_north"  

conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
conusmask.long_name = "CONUS mask" 
conusmask.units = "1 or 0 for inside/out of the CONUS" 

yyyymmddhh = rootgrp.createVariable('yyyymmddhh','i4',('time',))
yyyymmddhh.longname = "End of period of 12-h precipitation accumulation"

apcp_anal = rootgrp.createVariable('apcp_anal','f4',('time','ya','xa',), \
   zlib=True,least_significant_digit=2)
apcp_anal.units = "mm" 
apcp_anal.long_name = "Analyzed Accumulated Precipitation from CCPA v1" 
apcp_anal.valid_range = [0.,1000.]
apcp_anal.missing_value = -99.99

xva[:] = np.arange(nlonsa)
yva[:] = np.arange(nlatsa)

lonsa[:] = lons_anal
latsa[:] = lats_anal
conusmask[:] = conusmask_in

rootgrp.stream = "s4" # ????
rootgrp.title = \
   "CCPA 1/8 degree precip analysis over CONUS"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created May 2017 by Tom Hamill" 
rootgrp.institution =  "CCPA from NCEP/EMC"
rootgrp.platform = "Model" 
rootgrp.references = "http://www.esrl.noaa.gov/psd/forecasts/reforecast2/" 

ktr = 0
missing_values_anal = -99.99 * np.ones((nlatsa,nlonsa),dtype=np.float)

# ---- loop thru all dates ...
ktr = 0

for idate, date in zip(range(ndates),date_list):

    print idate, date

    # ---- read in CCPA 6-hourly precipitation analysis data

    precipver = False
    date_pplate  = date
    date_ppearly = dateshift(date, -6)
    date_lista   = [date_ppearly, date_pplate]
    print 'loading analysis data for ',date_lista
    tp_anal = np.zeros((nlatsa, nlonsa),dtype=np.float)
    nogood = False
    for datea in date_lista:
        cyear = datea[0:4]
        chour = datea[8:10]
        cyyyymmdd = datea[0:8]
        cyyyymmddhh = datea
        infile = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+cyyyymmdd+\
            '/'+chour+'/ccpa.t'+chour+'z.06h.0p125.conus'
        print infile
        fexist_apcp = os.path.exists(infile)
        if fexist_apcp and nogood == False:
            ffcst = pygrib.open(infile)
            grb = ffcst.select(shortName='tp')[0]
            tp_anal = tp_anal + grb.values
            ffcst.close()
        else: 
            nogood = True
            tp_anal[:,:] = -99.99
            print 'missing precip analysis data for ', date
    
    tp_anal = np.where(tp_anal < 80000., tp_anal, missing_values_anal)
    print 'ktr, max(tp_anal) = ',ktr, np.max(tp_anal)

    timev[ktr] = datetohrs(date)  # in our file, since 0 AD
    yyyymmddhh[ktr] = int(date)
    apcp_anal[ktr] = tp_anal
    ktr = ktr + 1 

print 'writing to ',outfilename
rootgrp.close()

