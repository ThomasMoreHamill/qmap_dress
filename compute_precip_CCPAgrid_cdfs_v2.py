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

from get_cdf_precip_anal_v2 import get_cdf_precip_anal_v2
from get_cdf_cmc_v2 import get_cdf_cmc_v2
from get_cdf_ecmwf_ncep_v2 import get_cdf_ecmwf_ncep_v2


# --- queries from command line

center = sys.argv[1]
cleade = sys.argv[2]
date_end = sys.argv[3]
idaysbefore = -1 -int(int(cleade)/24)
ihoursbefore = idaysbefore*24
print 'idaysbefore, ihoursbefore = ',idaysbefore, ihoursbefore
date_end_shift = dateshift(date_end, ihoursbefore)
date_begin_shift = dateshift(date_end_shift, -24*61)
print 'date_begin_shift, date_end_shift = ', date_begin_shift, date_end_shift

# --- initialize stuff

date_list = daterange(date_begin_shift, date_end_shift, 24)
print 'date_list = ', date_list
ndates = len(date_list)
date_middle = date_list[ndates/2]
cmonthno = date_middle[4:6]
imo = int(cmonthno)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonth = cmonths[imo]
if center == 'ECMWF':
    nmembers = 50
else:
    nmembers = 20

# ---- define the precipitation amount thresholds that we will calculate CDF at.

npct = 99 + 8 # quantiles, from 0.01 to 0.99 by 0.01, also (.0001, .005, .001, .0005), 
              # and 1 - (.0001, .005, .001, .0005)

xc = range(90)
thresh = [.001,.003,.005,.01,.03, .05,.07,.1,.2,.3,  .4,.5,.6,.7,.8,  \
      .9,1.0,1.2, 1.4, 1.6,    1.8, 2.0, 2.25, 2.5, 2.75,   3.0, 3.5, 4.0, 4.5, 5.0, \
      6.0,7.0,8.0,9.0,10.0,  11.0,12.0,13.0,14.0,15.0,  16.0,17.0,18.0,19.0,19.5,  \
      20.0,22.5,25.,27.5,30.0,   32.5,35.0,37.5,40.0,42.5,  45.0,50.0,55.0,60.0,65.0,  \
      70.0,75.0,80.0,85.0,90.0,   95.0,100.0,105.0,110.0,120.0,  130.0,140.0,150.0,160.0,170.0,  \
      180.0,190.0,200.0,220.0,240.0,   260.0,280.0,300.0,325.0,350.0,   400.0,500.0,600.0,700.0,1000.]
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


# ---- read in the supplemental locations information for the CCPA grid.

infile = '/Users/thamill/precip/supplemental_locations_eighth_degree_'+cmonth+'_v9.nc'
print infile
nc = Dataset(infile)
conusmask_lg_in = nc.variables['conusmask'][:,:]
lons_lg = nc.variables['longitudes'][:,:]
lats_lg = nc.variables['latitudes'][:,:]
nc.close()

# ---- extract subset on smaller CCPA grid.

nxa = 464
nya = 224
xoffset=21  
yoffset=22  
conusmask_in = np.zeros((nya,nxa), dtype=np.int16)
lons_anal = np.zeros((nya,nxa), dtype=np.float32)
lats_anal = np.zeros((nya,nxa), dtype=np.float32)

conusmask_in[:,:] = conusmask_lg_in[yoffset:yoffset+nya,xoffset:xoffset+nxa]
lons_anal[:,:] = lons_lg[yoffset:yoffset+nya,xoffset:xoffset+nxa]
lats_anal[:,:] = lats_lg[yoffset:yoffset+nya,xoffset:xoffset+nxa]


# ---- open and initialize the netCDF file we'll be writing to.

outfilename = '/data/thamill/ecmwf_data/'+center+'_CDF_flead'+cleade+\
    '_'+date_begin_shift+'_to_'+date_end_shift+'.nc'
print outfilename
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
    
xa = rootgrp.createDimension('xa',nxa)
xav = rootgrp.createVariable('xa','f4',('xa',))
xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
xav.units = "grid index (dimensionless)" 
    
ya = rootgrp.createDimension('ya',nya)
yav = rootgrp.createVariable('ya','f4',('ya',))
yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
yav.units = "grid index (dimensionless)"
    
pct = rootgrp.createDimension('pct',npct)
pctv = rootgrp.createVariable('pct','f4',('pct'))
pctv.long_name = "quantiles of the distribution"
pctv.units = "fraction"

time = rootgrp.createDimension('time',None)
timev = rootgrp.createVariable('time','i4',('time',))
timev.units = "index into the file for time dimension"

thrnum = rootgrp.createDimension('thrnum',nthresh)
thrnumv = rootgrp.createVariable('thrnum','i4',('thrnum',))
thrnumv.long_name = "Threshold iterator (0:nthresh)" 
thrnumv.units = " "  
    
thrval = rootgrp.createDimension('thrval',nthresh)
thrvalv = rootgrp.createVariable('thrval','f4',('thrval',))
thrvalv.long_name = "Precip thresholds (mm) that precip_CDFs evaluated at" 
thrvalv.units = "K" 
    
lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
lonsa.long_name = "longitude" 
lonsa.units = "degrees_east" 
    
latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
latsa.long_name = "latitude" 
latsa.units = "degrees_north" 

conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
conusmask.units=""

panal_CDF = rootgrp.createVariable('panal_CDF','f4',('time','thrnum','ya','xa',),
    zlib=True,least_significant_digit=3)  
panal_CDF.units = "" 
panal_CDF.long_name = "Cumulative distribution function of analyzed precip, defined at thrval" 
panal_CDF.valid_range = [0.0,1.0]
panal_CDF.missing_value = -99.99
  
# ----  set for CMC differently
  
if center == 'CMC': 
    
    ens = rootgrp.createDimension('ens',nmembers)
    ensv = rootgrp.createVariable('ensv','i4',('ens',))
    ensv.long_name = "Ensemble member number" 
    ensv.units = " " 

    pfcst_CDF = rootgrp.createVariable('pfcst_CDF','f4',('time','thrnum','ens','ya','xa',),
        zlib=True,least_significant_digit=3)  
    pfcst_CDF.units = "" 
    pfcst_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    pfcst_CDF.valid_range = [0.0,1.0]
    pfcst_CDF.missing_value = -99.99
    
else:
    
    pfcst_CDF = rootgrp.createVariable('pfcst_CDF','f4',('time','thrnum','ya','xa',),
        zlib=True,least_significant_digit=3)  
    pfcst_CDF.units = "" 
    pfcst_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    pfcst_CDF.valid_range = [0.0,1.0]
    pfcst_CDF.missing_value = -99.99
    
rootgrp.latcorners = [lats_anal[0,0], lats_anal[0,-1], lats_anal[-1,0], lats_anal[-1,-1]]
rootgrp.loncorners = [lons_anal[0,0], lons_anal[0,-1], lons_anal[-1,0], lons_anal[-1,-1]]

rootgrp.stream = "s4" # ????
rootgrp.title = "anal CDF"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created 8 May 2017 by Tom Hamill" 
rootgrp.institution = "ESRL/PSD using CCPA data from NCEP/EMC"
rootgrp.platform = "Model" 
rootgrp.references = "" 

pctv[:] = pctvalues[:] 
thrnumv[:] = xc[:]
thrvalv[:] = thresh[:]
if center == 'CMC':
    ens[:] = range(nmembers)




# ---- define the corner lat/lons for output netcdf.  Also x and y coordinates,
#      lat, lon, conusmask

llcrnrlat = lats_anal[0,0]
llcrnrlon = lons_anal[0,0]
urcrnrlat = lats_anal[-1,-1]
urcrnrlon = lons_anal[-1,-1]

xav[:]   = np.arange(nxa)
yav[:]   = np.arange(nya)
lonsa[:]  = lons_anal
latsa[:]  = lats_anal

conusmask[:] = conusmask_in

# --- loop thru forecasts and analyses for each date in date_list    

for idate, date in zip(range(ndates), date_list):
    
    # ---- set up the CDF arrays

    CDFa = np.zeros((nthresh,nya,nxa), dtype=np.float64) # analyzed
    CDFworka = np.zeros((nthresh,nya,nxa), dtype=np.float64)
    if center == 'CMC':
        CDFf = np.zeros((nthresh,nmembers,nya,nxa), dtype=np.float64) # forecast
        CDFworkf = np.zeros((nthresh,nmembers,nya,nxa), dtype=np.float64)
    else:
        CDFf = np.zeros((nthresh,nya,nxa), dtype=np.float64) # forecast
        CDFworkf = np.zeros((nthresh,nya,nxa), dtype=np.float64)
        
    icdf_a = np.zeros((nya,nxa), dtype=np.float32)
    icdf_f = np.zeros((nya,nxa), dtype=np.float32)
    
    # ---- read in ccpa precip analysis data for this date (12-h accum)

    print 'processing date = ', date
    if idate == 0: 
        infilename = '/data/thamill/Rf2_tests/ccpa_v1/precip_ccpav1_'+\
            '2002010200_to_2016123100.nc'            
        print infilename
        afile1 = Dataset(infilename,"r")
        yyyymmddhh = afile1.variables['yyyymmddhh'][:]
    
    fdate_today = int(dateshift(date, int(cleade)))
    idx_today = int(np.where(yyyymmddhh == fdate_today)[0])
        
    if idx_today >= 0:
        apcp_anal = afile1.variables['apcp_anal'][idx_today,:,:]
    else:
        apcp_anal = -99.99*np.ones((nya, nxa), dtype=np.float32)
            
    # ---- read in the forecast data for this date, and for the date 12 h previous
    #      then subtract to get accumulated precip during this period
    
    infile = '/Users/thamill/precip/ecmwf_data/'+center+'_'+date+'_leadtime'+cleade+'h.nc'
    print infile
    nc = Dataset(infile)
    apcp_fcst_ens_late = nc.variables['apcp_fcst_ens'][:,:,:] # mbr,x,y
    
    cleadb = str(int(cleade)-12)
    infile = '/Users/thamill/precip/ecmwf_data/'+center+'_'+date+'_leadtime'+cleadb+'h.nc'
    print infile
    nc = Dataset(infile)
    apcp_fcst_ens_early = nc.variables['apcp_fcst_ens'][:,:,:] # mbr,x,y
    
    apcp_fcst_ens = apcp_fcst_ens_late - apcp_fcst_ens_early

    # ---- populate work arrays with CDF data for this date.

    print 'get_cdf_precip_anal '
    CDFworka, icounta, istata = get_cdf_precip_anal_v2(nthresh, nya, nxa, \
        thresh, apcp_anal, conusmask_in)
    print 'CDFworka[0:15,nya/2,nxa/2] = ',CDFworka[0:15,nya/2,nxa/2]
    print 'icounta[nya/2,nxa/2]  = ',icounta[nya/2,nxa/2]
    
    if center == 'CMC':
        print 'get_cdf_precip_cmc '
        CDFworkf, icountf, istatf = get_cdf_cmc(nthresh, nmembers, nya, nxa, nsuppmax, \
            yoffset, xoffset, thresh, apcp_fcst_ens, xlocations, ylocations, \
            nsupplemental, conusmask_in)
        print 'CDFworkf[0:15,0,nya/2,nxa/2] = ',CDFworkf[0:15,0,nya/2,nxa/2]
        print 'icountf[nya/2,nxa/2]  = ',icountf[nya/2,nxa/2]
        
    else:
        print 'get_cdf_precip_ecmwf_ncep '
        CDFworkf, icountf, istatf = get_cdf_ecmwf_ncep(nthresh, nmembers, nya, nxa, nsuppmax, \
            yoffset, xoffset, thresh, apcp_fcst_ens, xlocations, ylocations, \
            nsupplemental, conusmask_in)     
        print 'CDFworkf[0:15,nya/2,nxa/2] = ',CDFworkf[0:15,nya/2,nxa/2]
        print 'icountf[nya/2,nxa/2]  = ',icountf[nya/2,nxa/2]
    
    # ---- if both forecast and analysis data exist for this date, increment CDF arrays

# ---- divide thru by the number of samples.  keep track of cdfs for each member
#      in the Canadian system, since that doesn't have members with exchangeable
#      error statistics.

for it in range(nthresh):
    CDFa[it,:,:] = CDFa[it,:,:] / icdf_a[:,:]
if center == 'CMC':
    for imem in range(nmembers):
        CDFf[it,imem,:,:] = CDFf[it,imem,:,:] / icdf_f[:,:]
else:
    CDFf[it,:,:] = CDFf[it,:,:] / icdf_f[:,:]
    
# --- get quantiles associated with the precipitation amounts, both for the forecast and analyzed

precip_qret_a = np.zeros((npct,nya,nxa),dtype=np.float)
precip_qret_a, istat = get_quantiles_linear(nthresh,npct,nya,nxa,pctvalues,thresh,CDFa)

if center != 'CMC':
    precip_qret_f = np.zeros((npct,nya,nxa),dtype=np.float)
    precip_qret_f, istat = get_quantiles_linear(nthresh,npct,nya,nxa,pctvalues,thresh,CDFf)
else:
    precip_qret_f = np.zeros((npct,nmembers,nja,nia),dtype=np.float)
    precip_qret_f, istat = get_quantiles_linear_cmc(nthresh,npct,nmembers,nya,nxa,pctvalues,thresh,CDFf)



# --- write out the CDF record

panal_CDF[:] = CDFa
print 'sample panal_CDF = ',CDFa[:,nya/2,nxa/2]

pfcst_CDF[:] = CDFf
if center == 'CMC':
    print 'sample pfcst_CDF = ',CDFf[:,0,nya/2,nxa/2]
else:
    print 'sample pfcst_CDF = ',CDFf[:,nya/2,nxa/2]
    
# --- write out the quantiles record

panal_quantiles[:] = precip_qret_a
print 'sample panal_quantiles = ',precip_qret_a[:,nya/2,nxa/2]

pfcst_quantiles[:] = precip_qret_f
if center == 'CMC':
    print 'sample pfcst_quantiles = ',precip_qret_f[:,0,nya/2,nxa/2]
else:
    print 'sample pfcst_quantiles = ',precip_qret_f[:,nya/2,nxa/2]

rootgrp.close()
