import numpy as np
from netCDF4 import Dataset
import sys
import pygrib
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import os
import scipy.stats as stats

cleadb = sys.argv[1]  # enter begin lead time in hours
cleade = sys.argv[2]  # enter end lead time in hours
clon = sys.argv[3]
clat = sys.argv[4]
cmmdd = sys.argv[5] # the month / day of target
ileadb = int(cleadb)
ileade = int(cleade)
rlon = float(clon)
rlat = float(clat)
nmembers = 11
npct = 99 + 8 # quantiles, from 0.01 to 0.99 by 0.01, also (.0001, .005, .001, .0005), 
              # and 1 - (.0001, .005, .001, .0005)
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonth = cmonths[int(cmmdd[0:2])-1]

# ---- determine lat/lon information from a sample file

filename=\
'/data/thamill/Rf2_tests/ccpa/netcdf/refcstv2_precip_ccpav3_'+\
    cleadb+'_to_'+cleade+'.nc'
print filename
nc = Dataset(filename)
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
yyyymmddhh_fcstb = nc.variables['yyyymmddhh_fcstb'][:]
yyyymmddhh_fcste = nc.variables['yyyymmddhh_fcste'][:]
xfvals = nc.variables['xf'][:]
yfvals = nc.variables['yf'][:]
lons_fcst = nc.variables['lons_fcst'][:]
lats_fcst = nc.variables['lats_fcst'][:]
lons_fcst_1d = lons_fcst[0,:]
lats_fcst_1d = lats_fcst[:,0]
lons_anal = nc.variables['lons_anal'][:]
lats_anal = nc.variables['lats_anal'][:]
lons_anal_1d = lons_anal[0,:]
lats_anal_1d = lats_anal[:,0]
nya, nxa = np.shape(lons_anal)

# ---- read in the supplemental locations

xoffset = 21
yoffset = 22
infile = '/Users/thamill/precip/'+\
    'supplemental_locations_eighth_degree_'+cmonth+'_v9.nc'
print infile
nc2 = Dataset(infile)
nsupplemental = nc2.variables['nsupplemental'][:]
xlocations = nc2.variables['xlocations'][:]
ylocations = nc2.variables['ylocations'][:]
lons_supp = nc2.variables['longitudes'][:]
lats_supp = nc2.variables['longitudes'][:]
lons_supp_1d = lons_anal[0,:]
lats_supp_1d = lats_anal[:,0]
nc2.close()

# --- find indices of closest points

difflons_a = np.abs(rlon-lons_anal_1d)
difflats_a = np.abs(rlat-lats_anal_1d)
inearest_a = np.argmin(difflons_a)
jnearest_a = np.argmin(difflats_a)
print 'jnearest_a, inearest_a = ',jnearest_a, inearest_a

difflons_s = np.abs(rlon-lons_supp_1d)
difflats_s = np.abs(rlat-lats_supp_1d)
inearest_s = np.argmin(difflons_s)
jnearest_s = np.argmin(difflats_s)

# ---- define the precipitation amount thresholds that we will calculate CDF at.

xc = range(90)
thresh = [.001,.003,.005,.01,.03, .05,.07,.1,.2,.3,  .4,.5,.6,.7,.8,  \
      .9,1.0,1.2, 1.4, 1.6,    1.8, 2.0, 2.25, 2.5, 2.75,   3.0, 3.5, 4.0, 4.5, 5.0, \
      6.0,7.0,8.0,9.0,10.0,  11.0,12.0,13.0,14.0,15.0,  16.0,17.0,18.0,19.0,19.5,  \
      20.0,22.5,25.,27.5,30.0,   32.5,35.0,37.5,40.0,42.5,  45.0,50.0,55.0,60.0,65.0,  \
      70.0,75.0,80.0,85.0,90.0,   95.0,100.0,105.0,110.0,120.0,  130.0,140.0,150.0,160.0,170.0,  \
      180.0,190.0,200.0,220.0,240.0,   260.0,280.0,300.0,325.0,350.0,   400.0,500.0,600.0,700.0,1000.]
      
npct = 99 + 8 # quantiles, from 0.01 to 0.99 by 0.01, also (.0001, .005, .001, .0005),
              # and 1 - (.0001, .005, .001, .0005)      
pctvalues = [.0001, .0005, .001, .005, \
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
nthresh = len(xc)
print 'len xc, thresh = ',len(xc), len(thresh)

# ---- open and initialize the netCDF file we'll be writing to.

outfilename = '/Users/thamill/precip/ecmwf_data/CDF_refcst_byyear_lon'+\
    clon+'_lat'+clat+'_'+cleadb+'_to_'+cleade+'_'+cmmdd+'.nc'
print outfilename
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

thrnum = rootgrp.createDimension('thrnum',nthresh)
thrnumv = rootgrp.createVariable('thrnum','i4',('thrnum',))
thrnumv.long_name = "Threshold iterator (0:nthresh)" 
thrnumv.units = " "  
    
thrval = rootgrp.createDimension('thrval',nthresh)
thrvalv = rootgrp.createVariable('thrval','f4',('thrval',))
thrvalv.long_name = "Precip thresholds (mm) that precip_CDFs evaluated at" 
thrvalv.units = "K" 
     
pct = rootgrp.createDimension('pct',npct)
pctv = rootgrp.createVariable('pct','f4',('pct'))
pctv.long_name = "quantiles of the distribution"
pctv.units = "fraction"

time = rootgrp.createDimension('time',12) # 2002 to 2013
timev = rootgrp.createVariable('time','i4',('time',))
timev.units = "year" 

pforecast_CDF = rootgrp.createVariable('pforecast_CDF','f4',('time','thrnum'),
    zlib=True,least_significant_digit=3)  
pforecast_CDF.units = "" 
pforecast_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
pforecast_CDF.valid_range = [0.0,1.0]
pforecast_CDF.missing_value = -99.99

pfcst_quantiles = rootgrp.createVariable('pfcst_quantiles','f4',('time','pct'),
    zlib=True,least_significant_digit=2)
pfcst_quantiles.units = "mm"
pfcst_quantiles.long_name = "Precip amount associated with quantiles of forecast distribution"
pfcst_quantiles.valid_range = [0.0,10000.]
pfcst_quantiles.missing_value = -99.99

panal_CDF = rootgrp.createVariable('panal_CDF','f4',('time','thrnum'),
    zlib=True,least_significant_digit=3)  
panal_CDF.units = "" 
panal_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
panal_CDF.valid_range = [0.0,1.0]
panal_CDF.missing_value = -99.99

panal_quantiles = rootgrp.createVariable('panal_quantiles','f4',('time','pct'),
    zlib=True,least_significant_digit=2)
panal_quantiles.units = "mm"
panal_quantiles.long_name = "Precip amount associated with quantiles of forecast distribution"
panal_quantiles.valid_range = [0.0,10000.]
panal_quantiles.missing_value = -99.99

thrnumv[:] = xc[:]
thrvalv[:] = thresh[:]
pctv[:] = pctvalues[:]
timev[:] = range(2002,2014)
 
rootgrp.stream = "s4" # ????
rootgrp.title = \
    "Reforecast V2 forecast and CCPA anal CDF at chosen (lon,lat), one for each year + suppl locns"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created May 2017 by Hamill" 
rootgrp.institution = \
    "Reforecast from ERSL/PSD using NCEP/EMC GEFS, circa 2012"
rootgrp.platform = "Model" 
rootgrp.references = "http://www.esrl.noaa.gov/psd/forecasts/reforecast2/" 
  
# ---- cycle thru years, producing a CDF for each year including data from supplemental locations.

for iyear in range(2002,2014):
    cyear = str(iyear)
    centerdate = cyear + cmmdd + '00'
        
    idx_center = np.where(yyyymmddhh_init == int(centerdate))
    idx_tostart = idx_center[0] - 30
    idx_toend = idx_center[0] + 30

    print 'year = ',cyear,'idx range ', idx_tostart, idx_toend
        
    # ---- loop over dates
        
    CDFf = np.zeros((nthresh), dtype=np.float64)
    CDFa = np.zeros((nthresh), dtype=np.float64)
    nsamps_a = 0
    nsamps_f = 0
    for idx in range(idx_tostart, idx_toend):
            
        # ---- read in forecast and analyzed 
  
        if idx%10 == 0: print 'idx = ',idx
        apcp_fcst_ens = nc.variables['apcp_fcst_ens'][idx,:,:,:]
        apcp_anal = nc.variables['apcp_anal'][idx,:,:]
        #print 'apcp_anal[nya/2,0:nxa:10] = ', apcp_anal[nya/2,0:nxa:10]
          
        for isupp in range(nsupplemental[jnearest_a,inearest_a]):
            iloc = xlocations[isupp,jnearest_s-yoffset-1,inearest_s-xoffset-1]
            jloc = ylocations[isupp,jnearest_s-yoffset-1,inearest_s-xoffset-1]
            #print 'jloc,iloc = ',jloc,iloc            
            for ithr in range(nthresh):
                if apcp_anal[jloc,iloc] <= thresh[ithr] and \
                apcp_anal[jloc,iloc] > -98.:
                    CDFa[ithr] = CDFa[ithr] + 1
                if ithr == 0: nsamps_a = nsamps_a + 1

        # ---- interpolate forecast to analyzed grid
            
        apcp_fcst_analgrid = np.zeros((nya,nxa), dtype=np.float32)
        for imem in range(11):
            apcp_fcst_analgrid = interp(apcp_fcst_ens[imem,:,:], \
                lons_fcst_1d, lats_fcst_1d, lons_anal, lats_anal, \
                checkbounds=False, masked=False, order=1)
            #if imem == 0 : print 'apcp_fcst_analgrid[nya/2,0:nxa:10] = ',\
            #    apcp_fcst_analgrid[nya/2,0:nxa:10]
            for isupp in range(nsupplemental[jnearest_a,inearest_a]):
                iloc = xlocations[isupp,jnearest_s-yoffset-1,inearest_s-xoffset-1]
                jloc = ylocations[isupp,jnearest_s-yoffset-1,inearest_s-xoffset-1]
                for ithr in range(nthresh):
                    if apcp_fcst_analgrid[jloc,iloc] <= thresh[ithr] and \
                    apcp_fcst_analgrid[jloc,iloc] > -98.:
                        CDFf[ithr] = CDFf[ithr] + 1
                    if ithr == 0: nsamps_f = nsamps_f + 1

    CDFf[:] = CDFf[:] / float(nsamps_f)
    CDFa[:] = CDFa[:] / float(nsamps_a)

    # ---- get precip amounts associated with quantiles

    aquants = np.zeros((npct), dtype=np.float32)
    fquants = np.zeros((npct), dtype=np.float32)
    for iquant in range(npct):
        t = pctv[iquant]
        for ithr in range(nthresh-1):
            if CDFf[ithr] <= t and CDFf[ithr+1] > t:
                fac = (t-CDFf[ithr]) / (CDFf[ithr+1] - CDFf[ithr])
                value = (1.-fac)*thresh[ithr] + fac*thresh[ithr+1]
                fquants[iquant] = value
            if CDFa[ithr] <= t and CDFa[ithr+1] > t:
                fac = (t-CDFa[ithr]) / (CDFa[ithr+1] - CDFa[ithr])
                value = (1.-fac)*thresh[ithr] + fac*thresh[ithr+1]
                aquants[iquant] = value
      
    # --- write out the CDF record for this year
    

    

    pforecast_CDF[iyear-2002] = CDFf[:] 
    panal_CDF[iyear-2002] = CDFa[:]  
    
    pfcst_quantiles[iyear-2002] = fquants[:] 
    panal_quantiles[iyear-2002] = aquants[:]

    print 'sample pforecast_CDF = ',CDFf[:]
    print 'sample panal_CDF     = ',CDFa[:]
    print 'pct[0:60] = ',pctv[0:60]
    print 'fquants[0:60] = ',fquants[0:60]
    print 'aquants[0:60] = ',aquants[0:60]

rootgrp.close()
nc.close()
