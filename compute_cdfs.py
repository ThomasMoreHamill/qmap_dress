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
from fortran_routine2 import get_ceedeef_precip, get_quantiles_linear


nmembers = 11
npct = 99 + 8 # quantiles, from 0.01 to 0.99 by 0.01, also (.0001, .005, .001, .0005), 
              # and 1 - (.0001, .005, .001, .0005)

# ---- determine lat/lon information and masks from a sample file
#      conusmask is the strict mask for inside CONUS, practical_mask is the 
#      set of points with CCPA data not missing (1=available, 0 = no)

filename= 'precip_12h_anal_ccpa_47km_2polarstereo.nc'
print filename
nc = Dataset(filename)
lats_anal = nc.variables['lats_anal'][:]
lons_anal = nc.variables['lons_anal'][:]
xavals = nc.variables['xa'][:]
yavals = nc.variables['ya'][:]
nia = len(xavals)
nja = len(yavals)
conusmask_in = nc.variables['conusmask'][:,:]
practical_mask_in = nc.variables['practical_mask'][:,:]
zeros = np.zeros((nja,nia),dtype=np.int16)
conusmask_in = np.where(conusmask_in < 0, zeros, conusmask_in)

nc.close()
print 'nia, nja = ', nia, nja
print 'lons_anal[0,:] = ',lons_anal[0,:]
print 'lats_anal[:,0] = ',lats_anal[:,0]

# ---- define the corner lat/lons.

llcrnrlat = lats_anal[0,0]
llcrnrlon = lons_anal[0,0]
urcrnrlat = lats_anal[-1,-1]
urcrnrlon = lons_anal[-1,-1]

# ---- define the precipitation amount thresholds that we will calculate CDF at.

xc = range(90)
thresh = [.001,.003,.005,.01,.03, .05,.07,.1,.2,.3,  .4,.5,.6,.7,.8,  \
      .9,1.0,1.2, 1.4, 1.6,    1.8, 2.0, 2.25, 2.5, 2.75,   3.0, 3.5, 4.0, 4.5, 5.0, \
      6.0,7.0,8.0,9.0,10.0,  11.0,12.0,13.0,14.0,15.0,  16.0,17.0,18.0,19.0,19.5,  \
      20.0,22.5,25.,27.5,30.0,   32.5,35.0,37.5,40.0,42.5,  45.0,50.0,55.0,60.0,65.0,  \
      70.0,75.0,80.0,85.0,90.0,   95.0,100.0,105.0,110.0,120.0,  130.0,140.0,150.0,160.0,170.0,  \
      180.0,190.0,200.0,220.0,240.0,   260.0,280.0,300.0,325.0,350.0,   400.0,500.0,600.0,700.0,1000.]
nthresh = len(xc)
print 'len xc, thresh = ',len(xc), len(thresh)

# --- loop thru and process each month, reading that data plus the surrounding two months.

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
icenter = [14,45,73,104,134,165,195,226,257,287,318,348] # 15th of each month, in Julian Days



infilename = 'precip_12h_anal_ccpa_47km_2polarstereo.nc'
print 'reading ',infilename
cfile1 = Dataset(infilename,"r")
yyyymmddhh_date_end = cfile1.variables['yyyymmddhh_date_end'][:]

for imonth, cmonth in zip(range(12), cmonths): 

    print 'processing month = ',cmonth

    # ---- initialize CDFs and counter

    CDFa = np.zeros((nthresh,nja,nia), dtype=np.float64) # analyzed
    CDFworka = np.zeros((nthresh,nja,nia), dtype=np.float64)
    icdf = 0

    # ---- open and initialize the netCDF file we'll be writing to.

    outfilename = '/Volumes/Drobo/ccpa_v1/0.125d/apcp_CCPA_CDF_24h_accum_'+cmonth+'.nc'
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
    pct = rootgrp.createDimension('pct',npct)
    pctv = rootgrp.createVariable('pct','f4',('pct'))
    pctv.long_name = "quantiles of the distribution"
    pctv.units = "fraction"
    thrnum = rootgrp.createDimension('thrnum',nthresh)
    thrnumv = rootgrp.createVariable('thrnum','i4',('thrnum',))
    thrnumv.long_name = "Threshold iterator (0:nthresh)" 
    thrnumv.units = " "  
    thrval = rootgrp.createDimension('thrval',nthresh)
    thrvalv = rootgrp.createVariable('thrval','f4',('thrval',))
    thrvalv.long_name = "Precip thresholds (mm) that precip_CDFs evaluated at" 
    thrvalv.units = "K"  
    time = rootgrp.createDimension('time',None)
    timev = rootgrp.createVariable('time','f4',('time',))
    timev.units = "hours since 1-1-1 00:00:0.0" 
    lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
    lonsa.long_name = "longitude" 
    lonsa.units = "degrees_east" 
    latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
    latsa.long_name = "latitude" 
    latsa.units = "degrees_north" 

    conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
    conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
    conusmask.units=""

    practical_mask = rootgrp.createVariable('practical_mask','i2',('ya','xa',))
    practical_mask.long_name = "mask for grid points with CCPA data (1=yes,0=no)"
    practical_mask.units=""

    panal_CDF = rootgrp.createVariable('panal_CDF','f4',('thrnum','ya','xa',),
        zlib=True,least_significant_digit=3)  
    panal_CDF.units = "" 
    panal_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    panal_CDF.valid_range = [0.0,1.0]
    panal_CDF.missing_value = -99.99

    panal_quantiles = rootgrp.createVariable('panal_quantiles','f4',('pct','ya','xa',),
        zlib=True,least_significant_digit=2)  
    panal_quantiles.units = "mm" 
    panal_quantiles.long_name = "Precip amount associated with quantiles of analyzed distribution" 
    panal_quantiles.valid_range = [0.0,10000.]
    panal_quantiles.missing_value = -99.99

    pctv[:] = [.0001, .0005, .001, .005, \
            .01, .02, .03, .04, .05, .06, .07, .08, .09, .10, \
            .11, .12, .13, .14, .15, .16, .17, .18, .19, .20, \
            .21, .22, .23, .24, .25, .26, .27, .28, .29, .30, \
            .31, .32, .33, .34, .35, .36, .37, .38, .39, .40, \
            .41, .42, .43, .44, .45, .46, .47, .48, .49, .50, \
            .51, .52, .53, .54, .55, .56, .57, .58, .59, .60, \
            .61, .62, .63, .64, .65, .66, .67, .68, .69, .70, \
            .71, .72, .73, .74, .75, .76, .77, .78, .79, .80, \
            .81, .82, .83, .84, .85, .86, .87, .88, .89, .90, \
            .91, .92, .93, .94, .95, .96, .97, .98, .99, \
            .995, .999, .9995, .9999]

    xav[:]   = np.arange(nia)
    yav[:]   = np.arange(nja)
    lonsa[:]  = lons_anal
    latsa[:]  = lats_anal

    conusmask[:] = conusmask_in
    practical_mask[:] = practical_mask_in

    thrnumv[:] = xc[:]
    thrvalv[:] = thresh[:]
 
    rootgrp.latcorners = [lats_anal[0,0], lats_anal[0,-1], lats_anal[-1,0], lats_anal[-1,-1]]
    rootgrp.loncorners = [lons_anal[0,0], lons_anal[0,-1], lons_anal[-1,0], lons_anal[-1,-1]]

    rootgrp.stream = "s4" # ????
    rootgrp.title = "anal CDF"
    rootgrp.Conventions = "CF-1.0"  # ????
    rootgrp.history = "Created 23 Nov 2015 by Hamill" 
    rootgrp.institution = "ESRL/PSD using CCPA data from NCEP/EMC"
    rootgrp.platform = "Model" 
    rootgrp.references = "http://www.esrl.noaa.gov/psd/forecasts/reforecast2/" 
    
    # ---- read in data, and augment cdf 
    #      information for that year if the sample is within the month of interest
    #      or the neighboring month

    nyears = len(range(2002,2016))
    print 'nsamps = ',92*nyears
    precipa = ma.zeros((nja,nia),dtype=np.float)
    precipa3d = ma.zeros((92*nyears,nja,nia),dtype=np.float)
    ipktr = 0

    for iyyyymmddhh, iday in zip(yyyymmddhh_date_end, range(len(yyyymmddhh_date_end))):

        if iday%2 == 0 and iday !=0: # files are 12-hourly data, we need to process 24-h, so every other

            if iday % 90 == 0: print iyyyymmddhh, imonth
            cyyyymmddhh = str(iyyyymmddhh)
            cmonth2 = cyyyymmddhh[4:6]
            imonth2 = int(cmonth2) - 1
            if imonth2 == 11 and imonth == 0: imonth2 = -1
            if imonth2 == 0 and imonth == 11: imonth2 = 12
            if imonth2 == imonth-1 or imonth2 == imonth or imonth2 == imonth+1:
                 precipa[:,:] = cfile1.variables['apcp_anal'][iday,:,:] + cfile1.variables['apcp_anal'][iday-1,:,:]
                 precipa3d[ipktr,:,:] = precipa[:,:]
                 ipktr = ipktr + 1
                 CDFworka, istata = get_ceedeef_precip(nthresh,nja,nia,thresh,precipa)
                 if istata == 1:
                     CDFa = CDFa + CDFworka
                     icdf = icdf + 1

    flicdf = np.float(icdf)
    for it in range(nthresh):
        CDFa[it,:,:] = CDFa[it,:,:] / flicdf


    # --- get quantiles

    precip_qret = np.zeros((npct,nja,nia),dtype=np.float)
    precip_qret,istat = get_quantiles_linear(nthresh,npct,nja,nia,pctv,thresh,CDFa)

    # --- Yan Luo has many points where the precip is outside CONUS, inside practical_mask, but 
    #     values are set to zero.  Evaluate and reset practical_mask and set CDFs back to 
    #     missing in such a case.

    nthinned = 0
    for jya in range (nja):
        for ixa in range(nia):
            if conusmask[jya,ixa] == 0 and practical_mask[jya,ixa] == 1 \
                    and CDFa[0,jya,ixa] == CDFa[-1,jya,ixa] :
                practical_mask[jya,ixa] = 0
                CDFa[:,jya,ixa] = -99.99
                precip_qret[:,jya,ixa] = -99.99
                nthinned = nthinned + 1
            if practical_mask[jya,ixa] == 0:
                CDFa[:,jya,ixa] = -99.99
                precip_qret[:,jya,ixa] = -99.99

    # --- write out the CDF record

    panal_CDF[:] = CDFa
    print 'sample panal_CDF     = ',CDFa[:,nja/2,nia/2]

    # --- write out the quantiles record

    panal_quantiles[:] = precip_qret
    print 'sample panal_quantiles     = ',precip_qret[:,nja/2,nia/2]

    print 'nthinned  = ',nthinned,' of ',nja*nia
    rootgrp.close()

cfile1.close()
