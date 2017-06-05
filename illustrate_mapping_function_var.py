import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os, sys

rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'

def bound(value, arr):
    return arr[arr <= value].max(), arr[arr >= value].min()

cleadb = '108' #sys.argv[1]  # enter begin lead time in hours
cleade = '120' #sys.argv[2]  # enter end lead time in hours
clon = '90 W' #sys.argv[3]
clat = '40 N' #sys.argv[4]
cmmdd = '0515' #sys.argv[5] # the month / day of target

# ---- read in CDFs

infile =  '/Users/thamill/precip/ecmwf_data/'\
    +'CDF_refcst_byyear_lon-90.0_lat40.0_108_to_120_0515.nc'
nc = Dataset(infile)
pct = nc.variables['pct'][:]
pfcst_quantiles = nc.variables['pfcst_quantiles'][:,:]
panal_quantiles = nc.variables['panal_quantiles'][:,:]
year = nc.variables['time'][:]       
nc.close()
print 'years = ',year

# --- now plot mapping functions

rcParams['legend.fontsize']='medium'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'

fig = plt.figure(figsize=(6.5,7.))


a1 = fig.add_axes([.09,.09,.88,.83])
title = r'May quantile mapping functions for $\lambda$ = '+clon+', $\phi$='+clat
a1.set_title(title,fontsize=15)

for ix in range(12):
    cletters = ['2002', '2003','2004','2005','2006',\
        '2007','2008','2009','2010','2011', '2012','2013']
    colors = ['LightGray','Red','Cyan', 'RoyalBlue','LightGreen',
        'Orange','Maroon','Orchid','DarkGreen','Black','Yellow','Purple']
    
    a1.plot(pfcst_quantiles[ix,:],panal_quantiles[ix,:],\
        '-',color=colors[ix],lw=2,label=cletters[ix])
    
    if ix == 0: 
        a1.set_xlabel('Forecast amount (mm)', fontsize=15)
        a1.set_ylabel('Analyzed amount (mm)', fontsize=15)
    a1.set_ylim(-0.01,25.1)
    a1.set_xlim(-0.01,25.1)
    a1.set_xticks(range(0,26,5))
    a1.set_yticks(range(0,26,5))
    a1.grid(color='LightGray',linestyle='--',linewidth=1)
    a1.grid (True)
    a1.legend(loc=0)


plot_title = 'mapping_functions.pdf'
print 'saving plot to file = ',plot_title
plt.savefig(plot_title)
print 'Plot done'


