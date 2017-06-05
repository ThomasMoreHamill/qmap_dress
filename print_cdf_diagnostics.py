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

from get_cdf_precip_anal import get_cdf_precip_anal
from get_cdf_cmc import get_cdf_cmc
from get_cdf_ecmwf_ncep import get_cdf_ecmwf_ncep
from get_quantiles_linear import get_quantiles_linear
from get_quantiles_linear_cmc import get_quantiles_linear_cmc

# --- queries from command line

center = sys.argv[1]
cleade = sys.argv[2]
date_begin_shift = sys.argv[3]
date_end_shift = sys.argv[4]

# ---- open and initialize the netCDF file we'll be writing to.

infilename = '/Users/thamill/precip/ecmwf_data/'+center+'_CDF_flead'+cleade+\
    '_'+date_begin_shift+'_to_'+date_end_shift+'.nc'
print infilename
nc = Dataset(infilename,'r')
    
lons = nc.variables['lonsa'][:,:]
lats = nc.variables['latsa'][:,:]
conusmask = nc.variables['conusmask'][:,:]
panal_CDF = nc.variables['panal_CDF'][:,:,:]
pctiles = nc.variables['pct'][:]
threshes = nc.variables['thrval'][:]
panal_quantiles = nc.variables['panal_quantiles'][:,:,:]
if center == 'CMC': 
    pfcst_CDF = nc.variables['pfcst_CDF'][:,:,:,:]
    pfcst_quantiles = nc.variables['pfcst_quantiles'][:,:,:,:]
else:
    pfcst_CDF = nc.variables['pfcst_CDF'][:,:,:]
    pfcst_quantiles = nc.variables['pfcst_quantiles'][:,:,:]
nya, nxa = np.shape(lons)
nc.close()

print 'lats[:,0] = ',lats[:,0]
print 'lons[0,:] = ',lons[0,:]
print 'conusmask[nya/2,0:-1:10] = ',conusmask[nya/2,0:-1:10]
print 'pctiles = ', pctiles
print 'threshes = ', threshes 
print 'panal_CDF[0:20,nya/4,nxa/2] = ', panal_CDF[0:20,nya/4,nxa/2]
if center == 'CMC':
    print 'pfcst_CDF[0:20,0,nya/4,nxa/2] = ', pfcst_CDF[0:20,0,nya/4,nxa/2]
else:
    print 'pfcst_CDF[0:20,nya/4,nxa/2] = ', pfcst_CDF[0:20,nya/4,nxa/2]
    
print 'panal_CDF[0:20,3*nya/4,nxa/2] = ', panal_CDF[0:20,3*nya/4,nxa/2]
if center == 'CMC':
    print 'pfcst_CDF[0:20,0,3*nya/4,nxa/2] = ', pfcst_CDF[0:20,0,3*nya/4,nxa/2]
else:
    print 'pfcst_CDF[0:20,3*nya/4,nxa/2] = ', pfcst_CDF[0:20,3*nya/4,nxa/2]
    

