""" 
"""

from mpl_toolkits.basemap import Basemap
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from numpy import ma
import os, sys
from netCDF4 import Dataset
from verify_relia_bss_boots import verify_relia_bss_boots
import pygrib
from dateutils import daterange, dateshift
rcParams['legend.fontsize']='small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'

cleade = sys.argv[1]
cmodelcombo = sys.argv[2] # E, C, N, EC, EN, NC, or ENC permitted
cthresh = sys.argv[3]  # right now set for POP or 10mm

cyyyymmddhh_start = '2016040100'
cyyyymmddhh_end = '20160701000'
date_list = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)

ileade = int(cleade)
ileadb = ileade - 12
cleadb = str(ileadb)
if ileadb < 100: cleadb='0'+cleadb
if cthresh == 'POP':
    rthresh = 0.254
else:
    rthresh = 10.0

nclasses = 21 # 0 to 100% probability by 5%
nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 #

# ---- read in precipitation analyses

date_list_anal = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
print 'date_list_anal = ',date_list_anal
ndates = len(date_list_anal)
apcp_anal_t = np.zeros((nxa,nya,ndates), dtype=np.float32)
mninetynine = -99.99*np.ones((nya,nxa), dtype=np.float32)

for idate, date in zip(range(ndates), date_list_anal):

    print '------------------------- getting precipitation for idate = ',idate,date
    date_fearly = dateshift(date, ileade-6)
    date_flate = dateshift(date, ileade)
    print 'date_fearly, date_flate = ', date_fearly, date_flate

    # --- read in the first of the two precip analysis files, 6 hourly
    
    infile = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+date_fearly[0:8]+'/18/ccpa.t18z.06h.0p125.conus.gb2'
    #infile = '/Users/thamill/ccpa_v1/ccpa.'+date_fearly[0:8]+'/18/ccpa.t18z.06h.0p125.conus.gb2'
    print infile
    fexist1 = os.path.exists(infile)
    if fexist1:
        afile1 = pygrib.open(infile)
        grb1 = afile1.select()[0]    # --- read in the first of the two precip analysis files, 6 hourly
    
    infile2 = '/data/thamill/Rf2_tests/ccpa_v1/0.125d/ccpa.'+date_flate[0:8]+'/00/ccpa.t00z.06h.0p125.conus.gb2'
    #infile2 = '/Users/thamill/ccpa_v1/ccpa.'+date_fearly[0:8]+'/18/ccpa.t18z.06h.0p125.conus.gb2'
    print infile2
    fexist2 = os.path.exists(infile2)
    if fexist2:
        afile2 = pygrib.open(infile2)
        grb2 = afile2.select()[0]
    
    if fexist1 and fexist2:
        apcp_anal = grb1.values + grb2.values
        apcp_anal = np.where(apcp_anal > 500., mninetynine, apcp_anal)
    else:
        print 'no analysis data for this date'
        apcp_anal = -99.99*np.ones((nya,nxa),dtype=np.float32)
    print np.shape(apcp_anal)

    apcp_anal_t[:,:,idate] = np.transpose(apcp_anal[:,:])
    print 'min, max apcp_anal = ', np.min(apcp_anal), np.max(apcp_anal)

    afile1.close()
    afile2.close()

# ---- call fortran routine to do the reading in of forecasts and the generation
#      of reliability, frequency of use, and bss information.

relia = np.zeros(21, dtype=np.float32)
relia_raw = np.zeros(21, dtype=np.float32)
relia_raw_CMC = np.zeros(21, dtype=np.float32)
relia_raw_NCEP = np.zeros(21, dtype=np.float32)
relia_raw_ECMWF = np.zeros(21, dtype=np.float32)
relia_cdf = np.zeros(21, dtype=np.float32)
relia_cdf_CMC = np.zeros(21, dtype=np.float32)
relia_cdf_NCEP = np.zeros(21, dtype=np.float32)
relia_cdf_ECMWF = np.zeros(21, dtype=np.float32)

relia_05 = np.zeros(21, dtype=np.float32)
relia_raw_05 = np.zeros(21, dtype=np.float32)
relia_raw_05_NCEP = np.zeros(21, dtype=np.float32)
relia_raw_05_CMC = np.zeros(21, dtype=np.float32)
relia_raw_05_ECMWF = np.zeros(21, dtype=np.float32)
relia_cdf_05 = np.zeros(21, dtype=np.float32)
relia_cdf_05_NCEP = np.zeros(21, dtype=np.float32)
relia_cdf_05_CMC = np.zeros(21, dtype=np.float32)
relia_cdf_05_ECMWF = np.zeros(21, dtype=np.float32)

relia_95 = np.zeros(21, dtype=np.float32)
relia_raw_95 = np.zeros(21, dtype=np.float32)
relia_raw_95_NCEP = np.zeros(21, dtype=np.float32)
relia_raw_95_CMC = np.zeros(21, dtype=np.float32)
relia_raw_95_ECMWF = np.zeros(21, dtype=np.float32)
relia_cdf_95 = np.zeros(21, dtype=np.float32)
relia_cdf_95_NCEP = np.zeros(21, dtype=np.float32)
relia_cdf_95_CMC = np.zeros(21, dtype=np.float32)
relia_cdf_95_ECMWF = np.zeros(21, dtype=np.float32)

frequse_95 = np.zeros(21, dtype=np.float32)
frequse_raw_95 = np.zeros(21, dtype=np.float32)
frequse_raw_95_NCEP = np.zeros(21, dtype=np.float32)
frequse_raw_95_CMC = np.zeros(21, dtype=np.float32)
frequse_raw_95_ECMWF = np.zeros(21, dtype=np.float32)
frequse_cdf_95 = np.zeros(21, dtype=np.float32)
frequse_cdf_95_NCEP = np.zeros(21, dtype=np.float32)
frequse_cdf_95_CMC = np.zeros(21, dtype=np.float32)
frequse_cdf_95_ECMWF = np.zeros(21, dtype=np.float32)

bss = 0.
bss_raw = 0.
bss_raw_CMC = 0.
bss_raw_NCEP = 0.
bss_raw_ECMWF = 0.
bss_cdf = 0.
bss_cdf_CMC = 0.
bss_cdf_NCEP = 0.
bss_cdf_ECMWF = 0.

nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 #

print 'shape(apcp_anal_t)= ',np.shape(apcp_anal_t)
print verify_relia_bss_boots.__doc__

relia_raw, relia_raw_05, relia_raw_95, frequse_raw, bss_raw, \
    relia_raw_NCEP, relia_raw_05_NCEP, relia_raw_95_NCEP, frequse_raw_NCEP, bss_raw_NCEP, \
    relia_raw_CMC, relia_raw_05_CMC, relia_raw_95_CMC, frequse_raw_CMC, bss_raw_CMC, \
    relia_raw_ECMWF, relia_raw_05_ECMWF, relia_raw_95_ECMWF, frequse_raw_ECMWF, bss_raw_ECMWF, \
    relia_cdf, relia_cdf_05, relia_cdf_95, frequse_cdf, bss_cdf, \
    relia_cdf_NCEP, relia_cdf_05_NCEP, relia_cdf_95_NCEP, frequse_cdf_NCEP, bss_cdf_NCEP, \
    relia_cdf_CMC, relia_cdf_05_CMC, relia_cdf_95_CMC, frequse_cdf_CMC, bss_cdf_CMC, \
    relia_cdf_ECMWF, relia_cdf_05_ECMWF, relia_cdf_95_ECMWF, frequse_cdf_ECMWF, bss_cdf_ECMWF = \
    verify_relia_bss_boots(cleade, cthresh, nclasses, \
    rthresh, cmodelcombo, date_list_anal, apcp_anal_t, nxa, nya, ndates)

    
    
# ---- now make multi-panel reliability diagrams for each of the forecasts

fig = plt.figure(figsize=(12.,6.5))
plt.suptitle(cthresh+' reliability diagrams for +'+cleadb+' to +'+cleade+' hour forecasts',fontsize=18)
rcParams['xtick.labelsize']='xx-small'
rcParams['ytick.labelsize']='xx-small'

for itype in range(8):
    if itype == 0:
        ov_ax = [.04,.55,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_raw
        relia05 = relia_raw_05
        relia95 = relia_raw_95
        frequse_out = frequse_raw
        bss_out = bss_raw
        extra = False
        ctitle = '(a) Raw Multi-Model Ensemble'
    elif itype == 1:
        ov_ax = [.29,.55,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_raw_NCEP
        relia05 = relia_raw_05_NCEP
        relia95 = relia_raw_95_NCEP
        frequse_out = frequse_raw_NCEP
        bss_out = bss_raw_NCEP
        extra = False
        ctitle = '(b) Raw NCEP ensemble '
    elif itype == 2:
        ov_ax = [.54,.55,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_raw_CMC
        relia05 = relia_raw_05_CMC
        relia95 = relia_raw_95_CMC
        frequse_out = frequse_raw_CMC
        bss_out = bss_raw_CMC
        extra = False
        ctitle = '(c) Raw CMC ensemble'
    elif itype == 3:
        ov_ax = [.79,.55,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_raw_ECMWF
        relia05 = relia_raw_05_ECMWF
        relia95 = relia_raw_95_ECMWF
        frequse_out = frequse_raw_ECMWF
        bss_out = bss_raw_ECMWF
        extra = False
        ctitle = '(d) Raw ECMWF ensemble'
    elif itype == 4:
        ov_ax = [.04,.06,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_cdf
        relia05 = relia_cdf_05
        relia95 = relia_cdf_95
        frequse_out = frequse_cdf
        bss_out = bss_cdf
        extra = False
        ctitle = '(e) Postprocessed MME'
    elif itype == 5:
        ov_ax = [.29,.06,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_cdf_NCEP
        relia05 = relia_cdf_05_NCEP
        relia95 = relia_cdf_95_NCEP
        frequse_out = frequse_cdf_NCEP
        bss_out = bss_cdf_NCEP
        extra = False
        ctitle = '(f) Postprocessed NCEP'
    elif itype == 6:
        ov_ax = [.54,.06,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_cdf_CMC
        relia05 = relia_cdf_05_CMC
        relia95 = relia_cdf_95_CMC
        frequse_out = frequse_cdf_CMC
        bss_out = bss_cdf_CMC
        extra = False
        ctitle = '(g) Postprocessed CMC'
    elif itype == 7:
        ov_ax = [.79,.06,.2,.35]
        a1 = fig.add_axes(ov_ax)
        relia_out = relia_cdf_ECMWF
        relia05 = relia_cdf_05_ECMWF
        relia95 = relia_cdf_95_ECMWF
        frequse_out = frequse_cdf_ECMWF
        bss_out = bss_cdf_ECMWF
        extra = False
        ctitle = '(h) Postprocessed ECMWF'

    # -- make reliability diagram for statistically downscaled and RAW

    a1.set_title(ctitle,fontsize=11)

    # --- add basic reliability diagram 

    yerrs = np.squeeze(np.array([[100.*(relia_out-relia05)],[100.*(relia95-relia_out)]]))
    
    #rcParams['xtick.labelsize']='medium'
    #rcParams['ytick.labelsize']='medium'
    strbss = 'BSS = %0.3f' %(bss_out)
    relia_m = ma.array(relia_out)
    relia_m = ma.masked_where(relia_m < 0.0, relia_m)
    probs = np.arange(nclasses) * 100./np.real(nclasses-1)
    a1.plot(probs,100.*relia_m,'o-',color='r',markersize=3)
    a1.errorbar(probs,100.*relia_m,yerr=yerrs,fmt='-',color='r')
    if extra == True: a1.plot(probs,100.*relia_mo,'o-',color='b')
    a1.plot([0,100],[0,100],'--',color='k')
    a1.set_ylabel('Observed Relative Frequency (%)',fontsize=10)
    a1.set_xlabel('Forecast Probability (%)',fontsize=10)
    a1.set_ylim(-1,101)
    a1.set_xlim(-1,101)

    # -- BSS inserted here

    a1.text(51,7,strbss,fontsize=10)
    
    # --- Frequency of usage inset diagram

    h0 = ov_ax[0] + ov_ax[2]*0.2
    h1 = ov_ax[1] + ov_ax[3]*0.6
    h2 = ov_ax[2]*0.3
    h3 = ov_ax[3]*0.3
    a2 = fig.add_axes([h0,h1,h2,h3])

    a2.bar(probs,frequse_out,width=5,bottom=0.001,log=True,color='red',edgecolor='black',align='center')
    if extra == True:
        for i in range(len(frequse_out_o)):
            a2.plot([probs[i]-2.5, probs[i]+2.5], [frequse_out_o[i],frequse_out_o[i]], color='Blue', linewidth=2)
    a2.set_xticks([0,25,50,75,100])
    a2.set_xlim(-3,103)
    a2.set_ylim(0.0,1.)
    a2.set_title('Frequency of usage',fontsize=7)
    a2.set_xlabel('Probability',fontsize=7)
    a2.set_ylabel('Frequency',fontsize=7)
    a2.hlines([.01,.1],0,100,linestyles='dashed',colors='black')


plot_title = 'relia_multipanel_'+cthresh+'_hour'+cleade+'.pdf'
print 'saving plot to file = ',plot_title
plt.savefig(plot_title)
print 'Plot done'


sys.exit()



    
# ---- now make reliability diagrams for each of the forecasts

for itype in range(8):
    if itype == 0:
        relia_out = relia_raw
        relia05 = relia_raw_05
        relia95 = relia_raw_95
        frequse_out = frequse_raw
        bss_out = bss_raw
        extra = False
        ctitle = 'Raw MME, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_raw_'+cthresh+'_hour'+cleade+'.pdf'
    elif itype == 1:
        relia_out = relia_raw_NCEP
        relia05 = relia_raw_05_NCEP
        relia95 = relia_raw_95_NCEP
        frequse_out = frequse_raw_NCEP
        bss_out = bss_raw_NCEP
        extra = False
        ctitle = 'Raw NCEP, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_raw_NCEP_'+cthresh+'_hour'+cleade+'.pdf'
    elif itype == 2:
        relia_out = relia_raw_CMC
        relia05 = relia_raw_05_CMC
        relia95 = relia_raw_95_CMC
        frequse_out = frequse_raw_CMC
        bss_out = bss_raw_CMC
        extra = False
        ctitle = 'Raw CMC, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_raw_CMC_'+cthresh+'_hour'+cleade+'.pdf' 
    elif itype == 3:
        relia_out = relia_raw_ECMWF
        relia05 = relia_raw_05_ECMWF
        relia95 = relia_raw_95_ECMWF
        frequse_out = frequse_raw_ECMWF
        bss_out = bss_raw_ECMWF
        extra = False
        ctitle = 'Raw ECMWF, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_raw_ECMWF_'+cthresh+'_hour'+cleade+'.pdf'
    elif itype == 4:
        relia_out = relia_cdf
        relia05 = relia_cdf_05
        relia95 = relia_cdf_95
        frequse_out = frequse_cdf
        bss_out = bss_cdf
        extra = False
        ctitle = 'Postprocessed MME, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_cdf_'+cthresh+'_hour'+cleade+'.pdf'
    elif itype == 5:
        relia_out = relia_cdf_NCEP
        relia05 = relia_cdf_05_NCEP
        relia95 = relia_cdf_95_NCEP
        frequse_out = frequse_cdf_NCEP
        bss_out = bss_cdf_NCEP
        extra = False
        ctitle = 'Postprocessed NCEP, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_cdf_NCEP_'+cthresh+'_hour'+cleade+'.pdf'
    elif itype == 6:
        relia_out = relia_cdf_CMC
        relia05 = relia_cdf_05_CMC
        relia95 = relia_cdf_95_CMC
        frequse_out = frequse_cdf_CMC
        bss_out = bss_cdf_CMC
        extra = False
        ctitle = 'Postprocessed CMC, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_cdf_CMC_'+cthresh+'_hour'+cleade+'.pdf' 
    elif itype == 7:
        relia_out = relia_cdf_ECMWF
        relia05 = relia_cdf_05_ECMWF
        relia95 = relia_cdf_95_ECMWF
        frequse_out = frequse_cdf_ECMWF
        bss_out = bss_cdf_ECMWF
        extra = False
        ctitle = 'Postprocessed ECMWF, lead ='+cleadb+' to '+cleade+' h'
        plot_title = 'relia_cdf_ECMWF_'+cthresh+'_hour'+cleade+'.pdf'        

    # -- make reliability diagram and frequency of usage

    fig = plt.figure(figsize=(6.,6.5))
    a1 = fig.add_axes([.12,.1,.83,.8])
    a1.set_title(ctitle,fontsize=16)

    # --- add basic reliability diagram 

    strbss = 'BSS = %0.2f' %(bss_out)
    relia_m = ma.array(relia_out)
    relia_m = ma.masked_where(relia_m < 0.0, relia_m)

    probs = np.arange(nclasses) * 100./np.real(nclasses-1)
    a1.plot(probs,100.*relia_m,'o-',color='r')
    
    if extra == True: a1.plot(probs,100.*relia_mo,'o-',color='b')
    a1.plot([0,100],[0,100],'--',color='k')
    a1.set_ylabel('Observed Relative Frequency (%)',fontsize=14)
    a1.set_xlabel('Forecast Probability (%)',fontsize=14)
    a1.set_ylim(-1,101)
    a1.set_xlim(-1,101)

    # -- BSS inserted here

    a1.text(51,7,strbss,fontsize=16,color='Red')
    
    # --- Frequency of usage inset diagram

    rcParams['xtick.labelsize']='small'
    rcParams['ytick.labelsize']='small'
    a2 = fig.add_axes([.25,.65,.34,.21])
    a2.bar(probs,frequse_out,width=5,bottom=0.001,log=True,color='red',edgecolor='black',align='center')
    if extra == True:
        for i in range(len(frequse_out_o)):
            a2.plot([probs[i]-2.5, probs[i]+2.5], [frequse_out_o[i],frequse_out_o[i]], color='Blue', linewidth=2)
        
    a2.set_xlim(-3,103)
    a2.set_ylim(0.0,1.)
    a2.set_title('Frequency of usage',fontsize=11)
    a2.set_xlabel('Probability',fontsize=10)
    a2.set_ylabel('Frequency',fontsize=10)
    a2.hlines([.01,.1],0,100,linestyles='dashed',colors='black')

    print 'saving plot to file = ',plot_title
    plt.savefig(plot_title)
    print 'Plot done'

