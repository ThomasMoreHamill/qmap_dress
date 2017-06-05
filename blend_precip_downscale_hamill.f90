PROGRAM blend_precip_downscale

! --- this program controls the execution of multi-model ensemble adjustment of
!     probabilities of precipitation and deterministic precipitation.
!
! POP.  The basic algorithm uses "quantile mapping" (aka CDF-based bias correction)
!    to adjust each ensemble forecast.  The process determines which percentile of
!    the CDF is associated with the current precipitation amount and then re-sets
!    that ensemble member to the analyzed CDF value associated with that same 
!    quantile.
!          Or at least that's the simple, earlier method.  Now, in order to deal
!    with some lack of spread in the ensemble, there is a more complicated method
!    whereby the quantile mapping also uses surrounding grid points and adds 
!    some stochasticism.
!          A last step is a location-dependent smoothing of the probability field
!    where the smoothing is less aggressive in areas with variabile terrain, e.g.,
!    the western US.
!
! QPF.   Simpler.  A quantile mapping of the ensemble mean to the analyzed distribution.
!    This tends to sharpen up the forecasts.
!
! coded by: Tom Hamill, NOAA, (303) 497-3060, tom.hamill@noaa.gov
!               August-September 2014
!           MDL adaptations (2014-early 2016) by Chrissy Finan and Eric Engle
!           Additional changes by Tom Hamill, May 2016
USE netcdf

!integer, parameter :: nxa = 2145 ! number of 1/8-deg analysis grid pts in the x-dir
!integer, parameter :: nya = 1597 ! number of 1/8-deg analysis grid pts in the y-dir
integer, parameter :: nxa = 515 ! number of expanded 1/8-deg analysis grid pts in the x-dir
integer, parameter :: nya = 262 ! number of expanded 1/8-deg analysis grid pts in the y-dir
integer, parameter :: nx47 = 297 ! number of grid pts in x-dir for 47-km polar stereo grid
integer, parameter :: ny47 = 169 !      "                y-dir 
integer, parameter :: nx95 = 137 ! number of grid pts in x-dir for 95-km polar stereo grid
integer, parameter :: ny95 = 81  !      "                y-dir 
integer, parameter :: nxs = 464  ! number of grid pts in x-dir for 1/8-deg analysis grid
integer, parameter :: nys = 224  ! number of grid pts in y-dir for 1/8-deg analysis grid 
!integer, parameter :: ncount = 625 ! max number of 1/8-deg grid pts that could be associated with fcst grid pt
integer, parameter :: ncount = 70 ! max number of 1/8-deg grid pts that could be associated with fcst grid pt
integer, parameter :: nens_ecmwf = 50 ! number of ECMWF perturbed ensemble members
integer, parameter :: nens_cmc   = 20 ! number of CMC perturbed ensemble members
integer, parameter :: nens_ncep  = 20 ! number of NCEP perturbed ensemble members
integer, parameter :: ndays = 31 ! number of forecast days in data set (00 UTC Jan 2014 forecasts)
integer, parameter :: ndays3mo = 91*14 ! +/- 45 days from center of month x 12 years (2002-2016)
integer, parameter :: npct = 91 ! number of thresholds for CDFs

integer, parameter :: window_size = 9 ! used in Savitzky-Golay smoothing of prob field
integer, parameter :: order = 3 ! fit 3rd-order polynomial in Savitzky Golay smoother
!integer, parameter :: nstride = 7 
!real, parameter :: stdran = 0.22 ! magnitude of noise for perturbing forecast quantile
integer :: nstride
real :: stdran ! magnitude of noise for perturbing forecast quantile

REAL, DIMENSION(npct) :: thresh ! the precip amount thresholds for CDFs

CHARACTER*2 chh
CHARACTER*256 panal_infile ! name of file with coarse and fine-res precip analyses, 2002-2013
CHARACTER*256 pfcst_infile ! name of file with netCDFized forecast data, deterministic and ensemble
CHARACTER*256 iccpa_infile ! has the lists of grid pts on CCPA grid assoc'd with polar-stereo grid
!CHARACTER*256 cdf_infile   ! name of CDF file to read
CHARACTER*256 outfile      ! name of flat fortran file with output prob forecasts
CHARACTER*256 outfile_d    ! name of deterministic forecast
CHARACTER*256 ptopo_infile ! topo netcdf file
CHARACTER*256 pmask_infile ! ccpa mask
CHARACTER*256 griboutfile  ! name of grib2 output file with output prob fore
CHARACTER*256 griboutfile_ecmwf  ! name of grib2 output file with output prob fore
CHARACTER*256 griboutfile_cmc  ! name of grib2 output file with output prob fore
CHARACTER*256 griboutfile_ncep  ! name of grib2 output file with output prob fore
CHARACTER*10 cyyyymmddhh   ! year,month,day,hour of initial time of forecast
CHARACTER*3 cleade         ! ending hour of precip forecast accumulation, 3 digits, e.g., '024'
CHARACTER*3 cpcpvar        !'qpf' or 'pop'
CHARACTER*1 cdfbcorr       ! Flag for applying CDF-based bias corrections (y|Y or n|N).
CHARACTER*5 cthresh        ! 'POP','1mm','2p5mm','5mm','10mm','25mm','50mm'
CHARACTER*20 cfield
INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! inherited from CCPA data set
INTEGER*2, DIMENSION(nxa,nya) :: ihavedata  ! a broader mask including coastal waters.
REAL*4, DIMENSION(window_size,window_size) :: weights ! weightings for Savitzky-Golay smooth

! ---- 47km arrays
!REAL, ALLOCATABLE, DIMENSION(:,:) :: ecmwf_deterministic  !deterministic precip forecast
!REAL, ALLOCATABLE, DIMENSION(:,:) :: cmc_control ! cmc control precip forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: climo_prob ! climatological prob
!REAL, ALLOCATABLE, DIMENSION(:,:) :: ncep_control ! ncep control precip forecast
!!REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_ensemble ! ecmwf ens precip forecast
!REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_ensemble ! cmc ens precip forecast
!REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble ! ncep ens precip forecast

! ---- 1/8 deg. Lat/Lon arrays (i.e. CCPA grid)
!REAL, ALLOCATABLE, DIMENSION(:,:) :: ecmwf_deterministic_ccpa  !deterministic precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:) :: cmc_control_ccpa ! cmc control precip forecast  on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:) :: climo_prob_ccpa ! climatological prob  on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:) :: ncep_control_ccpa ! ncep control precip forecast on 1/8-deg ccpa grid 
REAL, ALLOCATABLE, DIMENSION(:,:) :: gem_today_before ! grand ens mean for today's forecast before bias corr
REAL, ALLOCATABLE, DIMENSION(:,:) :: gem_today_after ! grand ensemble mean for today's after bias corr
!REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_ensemble_ccpa ! ecmwf ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_ensemble_ccpa ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble_ccpa ! ncep ens precip forecast on 1/8-deg ccpa grid

! ---- x9 arrays
!REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_deterministic_ccpa_x9  !deterministic precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_control_ccpa_x9 ! cmc control precip forecast  on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_control_ccpa_x9 ! ncep control precip forecast on 1/8-deg ccpa grid 
!REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ecmwf_ensemble_ccpa_x9 ! ecmwf ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_ccpa_x9 ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ncep_ensemble_ccpa_x9 ! ncep ens precip forecast on 1/8-deg ccpa grid


REAL, ALLOCATABLE, DIMENSION(:,:) :: rlonsa ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlatsa ! precip analysis grid lat/lons

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: precip_anal_fine ! precipitation analysis on fine grid

REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast  ! output downscaled prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw  ! output raw ensemble prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_ECMWF ! output raw NCEP ensemble prob 
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_NCEP ! output raw NCEP ensemble prob 
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_unsmoothed ! output raw NCEP ensemble prob 
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_CMC ! output raw NCEP ensemble prob 
REAL, ALLOCATABLE, DIMENSION(:,:) :: topo_eighth
REAL, ALLOCATABLE, DIMENSION(:,:) :: raw_weight 
REAL, ALLOCATABLE, DIMENSION(:,:) :: determ_smoothed

REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: precip_anal_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_deterministic_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_control_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_control_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_ensemble_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gem_cdf ! grand ensemble mean CDF

REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf  ! output cdf bcorr ensemble prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_NCEP ! output cdf NCEP ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_CMC

integer :: itot_c     ! flag for total model datasets loaded
integer :: icm_c_ccpa ! flag for cmc control data on 1/8-deg ccpa grid
integer :: icm_e_ccpa ! flag for cmc ensemble data on 1/8-deg ccpa grid
integer :: iec_d_ccpa ! flag for ecmwf deterministic forecast data on 1/8-deg ccpa grid
integer :: iec_e_ccpa ! flag for ecmwf ensemble on 1/8-deg ccpa grid
integer :: inc_c_ccpa ! flag for ncep control on 1/8-deg ccpa grid
integer :: inc_e_ccpa ! flag for ncep ensemble on 1/8-deg ccpa grid

integer :: ierr      ! return variable for BAOPEN
integer :: ios       ! return variable for Fortran I/O, Allocation statements
integer :: ifcstint  ! forecast interval. Set depending on cpcpvar.
integer :: ipcpvar   ! 0 = pop; 1 = qpf
integer :: icdfbcorr ! 0 = no; 1 = yes

integer :: iyyyymmddhh,jyyyymmddhh
integer :: iyear,imo,iday,ihour,idoy ! Parsed date variables from iyyyymmddhh
integer :: jyear,jmo,jday,jhour,jdoy ! Parsed date variables for valid date

! ---- Initialize

! ---- Print some parameters
write(6,*)' Parameters for blend_precip_downscale :'
write(6,100)nxa,nya,nxs,nys,ndays3mo
100 format('nxa: ',i0.1/,'nya: ',i0.1/,'nxs: ',i0.1/,&
           'nys: ',i0.1/,'ndays3mo: ',i0.1/)

! --- Via command line, read in the input year/mo/day/hr and the forecast resolution 
!     we're working with.  Process date to determine the day of the month as an integer
CALL getarg(1,cyyyymmddhh)  ! input year month day hour of initial condition, 'yyyymmddhh' format
CALL getarg(2,cpcpvar)      ! Precipitation Variable to run for ("pop" or "qpf")
CALL getarg(3,cleade)       ! forecast lead time for beginning of precip accum period, hours, e.g.'060'
CALL getarg(4,cdfbcorr)     ! "Y"/"N" for performing CDF bias correction
CALL getarg(5,cthresh)      ! Precip threshold

write(6,*)' Command line arguments:'
write(6,110)cyyyymmddhh,cpcpvar,cleade,cdfbcorr,cthresh
110 format(1x,'cyyyymmddhh: ',A/1x,'cpcpvar: ',A/1x,'cleade: ',A/1x,'cdfbcorr: ',A/1x,'cthresh: ',A/)

IF(TRIM(cthresh).eq.'POP')THEN
   rthresh = 0.254 ! mm
ELSE IF (TRIM(cthresh) .eq. '1mm') THEN
   rthresh = 1.0
ELSE IF (TRIM(cthresh) .eq. '2p5mm') THEN
   rthresh = 2.5
ELSE IF (TRIM(cthresh) .eq. '5mm') THEN
   rthresh = 5.0
ELSE IF (TRIM(cthresh) .eq. '10mm') THEN
   rthresh = 10.0
ELSE IF (TRIM(cthresh) .eq. '25mm') THEN
   rthresh = 25.0
ELSE IF (TRIM(cthresh) .eq. '50mm') THEN
   rthresh = 50.0
ELSE
   write(6,*)' **** Invalid threshold value: ',cthresh
   stop
ENDIF

! ---- Convert character based variables from command line to integers
READ(cyyyymmddhh,'(i10)')iyyyymmddhh
READ(cleade,'(i3)')ileade

! ---- Set icdfbcorr according cdfbcorr
if(cdfbcorr.eq."y".or.cdfbcorr.eq."Y")then
   icdfbcorr=1
elseif(cdfbcorr.eq."n".or.cdfbcorr.eq."N")then
   icdfbcorr=0
else
   write(6,*)' **** Invalid CDF Bias Correction value: ',cdfbcorr
   stop
endif

! ---- Set ifcstint and ipcpvar according to cpcpvar
if(cpcpvar.eq."pop".or.cpcpvar.eq."POP")then
   ifcstint=12
   ipcpvar=0
elseif(cpcpvar.eq."qpf".or.cpcpvar.eq."QPF")then
   ifcstint=6
   ipcpvar=1
endif
ileadb=ileade-ifcstint

! ---- Set stdran and nstride as a function of ileade
nstride=nint(3.+8.*ileade/168.)
stdran=0.5   ! 0.15+(real(ileade)/1680.)
write(6,*)'nstride = ',nstride
write(6,*)'stdran = ',stdran

! ---- Parse the initializtion date; determine the valid hour.  This is dependent on the
!      precip variable (cpcpvar), the model initialization (iyyyymmddhh), and forecast 
!      ending hour (ileade).
iendhour=0
call doy(iyyyymmddhh,iyear,imo,iday,ihour,idoy)
call updat(iyyyymmddhh,ileade,jyyyymmddhh)
call doy(jyyyymmddhh,jyear,jmo,jday,jhour,jdoy)
iendhour=jhour
write(6,*)'Model Initialization: ',iyyyymmddhh
write(6,*)'Forecast Projection Ending: ',ileade
write(6,*)'Forecast Valid Date/Hour: ',jyyyymmddhh,iendhour
write(6,*)'Precipitation Element Flag: ',ipcpvar

! ---- Get filenames from environment
call get_environment_variable('FORT20',panal_infile)
call get_environment_variable('FORT25',pmask_infile)
write(6,fmt='(A,A)')'CCPA NetCDF Archive file: ',trim(panal_infile)
write(6,fmt='(A,A)')'CCPA Mask file: ',trim(pmask_infile)

! ---- Allocate dynamic arrays
write(6,*)'Allocating dynamic arrays...'
!ALLOCATE (ecmwf_deterministic_ccpa(nxa,nya))
ALLOCATE (cmc_control_ccpa(nxa,nya))
ALLOCATE (ncep_control_ccpa(nxa,nya))
ALLOCATE (prob_forecast_unsmoothed(nxa,nya))
!ALLOCATE (ecmwf_ensemble_ccpa(nxa,nya,nens_ecmwf))
ALLOCATE (cmc_ensemble_ccpa(nxa,nya,nens_cmc))
ALLOCATE (ncep_ensemble_ccpa(nxa,nya,nens_ncep))
ALLOCATE (rlonsa(nxa,nya),rlatsa(nxa,nya))
ALLOCATE (precip_anal_fine(nxa,nya,ndays3mo))
ALLOCATE (prob_forecast(nxa,nya))
ALLOCATE (prob_forecast_raw(nxa,nya))
ALLOCATE (prob_forecast_raw_ECMWF(nxa,nya))
ALLOCATE (prob_forecast_raw_CMC(nxa,nya))
ALLOCATE (prob_forecast_raw_NCEP(nxa,nya))
ALLOCATE (prob_forecast_cdf(nxa,nya))
ALLOCATE (prob_forecast_cdf_CMC(nxa,nya))
ALLOCATE (prob_forecast_cdf_NCEP(nxa,nya))
ALLOCATE (climo_prob(nxa,nya))
ALLOCATE (topo_eighth(nxa,nya))
ALLOCATE (gem_today_before(nxa,nya))
ALLOCATE (gem_today_after(nxa,nya))
ALLOCATE (determ_smoothed(nxa,nya))
ALLOCATE (raw_weight(nxa,nya))
ALLOCATE (cmc_control_ccpa_x9(9,nxa,nya))
!ALLOCATE (ecmwf_deterministic_ccpa_x9(9,nxa,nya))
ALLOCATE (ncep_control_ccpa_x9(9,nxa,nya))
ALLOCATE (cmc_ensemble_ccpa_x9(9,nxa,nya,nens_cmc))
!ALLOCATE (ecmwf_ensemble_ccpa_x9(9,nxa,nya,nens_ecmwf))
ALLOCATE (ncep_ensemble_ccpa_x9(9,nxa,nya,nens_ncep))

! ---- Read in the 1/8-degree precipitation analyses.  For a given month, 
!      load in only the precip analyses for the surrounding months, 
!      e.g., for January load up Dec-Jan-Feb.
write(6,*)'Calling read_precip_analyses_x9'
CALL read_precip_analyses_x9(nxa,nya,nxs,nys,ndays3mo,ipcpvar,&
                             iyyyymmddhh,iendhour,panal_infile,&
                             pmask_infile,precip_anal_fine,conusmask)
ihavedata = conusmask
write(6,*)'Max, Min ihavedata = ',maxval(ihavedata),minval(ihavedata)
write(6,*)'Max, Min conusmask = ',maxval(conusmask),minval(conusmask)

! ---- Read precipitation foreasts for all models/ensemble members from TDLPACK Vector files.
!      These files contains the 1/8 deg. grid.
write(6,*)'Calling read_forecasts_ccpa_tdlp'
CALL read_forecasts_tdlp_vect(nxa,nya,nens_cmc,nens_ncep,ipcpvar,iyyyymmddhh,ileade,&
                              cmc_control_ccpa,ncep_control_ccpa,&
                              cmc_ensemble_ccpa,ncep_ensemble_ccpa,&
                              icm_c_ccpa,inc_c_ccpa,&
                              icm_e_ccpa,inc_e_ccpa,ier)

! ---- Check to make sure data is avaiable
itot_c=icm_c_ccpa+inc_c_ccpa+icm_e_ccpa+inc_e_ccpa
if(itot_c.eq.0)then
   write(6,*)'**** NO MODEL DATA AVAILABLE FOR DOWNSCALE. PROGRAM STOPPING.'
   call w3tage('BLEND_PRECIP_DOWNSCALE')
   stop
else 
   if(icm_c_ccpa.eq.0)then
      write(6,*)'**** CMC CONTROL DATA MISSING. ARRAY SET TO MISSING.'
      cmc_control_ccpa = -99.99
   else 
      write(6,*)'CMC dataset loaded'
   end if

   if(inc_c_ccpa.eq.0)then
      write(6,*)'**** NCEP CONTROL DATA MISSING. ARRAY SET TO MISSING.'
      ncep_control_ccpa = -99.99
   else
      write(6,*)'GFS dataset loaded'
   end if

   if(icm_e_ccpa.eq.0)then
      write(6,*)'**** CMC ENSEMBLE DATA MISSING. ARRAY SET TO MISSING.'
      cmc_ensemble_ccpa = -99.99
   else
      write(6,*)'CMCE dataset loaded'
   end if

   if(inc_e_ccpa.eq.0)then
      write(6,*)'**** NCEP ENSEMBLE DATA MISSING. ARRAY SET TO MISSING.'
      ncep_ensemble_ccpa = -99.99
   else
      write(6,*)'GEFS dataset loaded'
   end if
end if

write(6,120)
120 format('MODEL AVAILABILITY:'/'     1/8 DEG.')
write(6,130)'CMC',icm_c_ccpa
write(6,130)'CMCE',icm_e_ccpa
write(6,130)'GFS',inc_c_ccpa
write(6,130)'GEFS',inc_e_ccpa
130 format(A4,1X,I2,4X,I2)

! ---- Read in the 1/8 degree topography data (created by create_eighth_degree_conus_elev.f90)
!      This will be used to determine how much smoothing of the prob forecast to do, 
!      with more smoothing applied in the flat regions of the CONUS than in mountainous 
!      regions, since we want to preserve prob detail related to terrain features.
call get_environment_variable('FORT24',ptopo_infile)
write(6,fmt='(A,A)')'1/8 Deg. Terrain File: ',trim(ptopo_infile)
CALL check(nf90_open(ptopo_infile,NF90_NOWRITE,netid))
CALL check(nf90_inq_varid(netid,"MTERH_surface",ivar))
CALL check(nf90_get_var(netid,ivar,topo_eighth,start=(/1,1,1/),count=(/nxa,nya,1/)))
CALL check(nf90_close(netid))

! ---- For purposes of having a baseline for comparison, generate an ensemble probability
!      simply from the relative frequency.
if(ipcpvar.eq.0)then

   write(6,*)'Calling raw_ensemble_probs'
   CALL raw_ensemble_probs(nxa,nya,nens_cmc,nens_ncep,rthresh,cmc_control_ccpa,&
                        ncep_control_ccpa,cmc_ensemble_ccpa,ncep_ensemble_ccpa,&
                        icm_c_ccpa,inc_c_ccpa,inc_e_ccpa,icm_e_ccpa,prob_forecast_raw,&
                        prob_forecast_raw_CMC,prob_forecast_raw_NCEP)
end if

! ---- Read in and apply CDF-based bias corrections if indicated on command line

PRINT *,'icdfbcorr = ', icdfbcorr
IF(icdfbcorr.eq.1)THEN

   ! Allocate dynamic CDF arrays
   ALLOCATE(precip_anal_cdf(nxa,nya,npct))
   ALLOCATE(ecmwf_deterministic_cdf(nxa,nya,npct))
   ALLOCATE(ncep_control_cdf(nxa,nya,npct))
   ALLOCATE(cmc_control_cdf(nxa,nya,npct))
   ALLOCATE(ecmwf_ensemble_cdf(nxa,nya,npct))
   ALLOCATE(ncep_ensemble_cdf(nxa,nya,npct))
   ALLOCATE(cmc_ensemble_cdf(nxa,nya,nens_cmc,npct))
   ALLOCATE(gem_cdf(nxa,nya,npct))

   !write(6,*)'Calling read_cdfs'
   !CALL read_cdfs(nxa,nya,npct,ipcpvar,cdf_infile,thresh,precip_anal_cdf,&
   !            ecmwf_deterministic_cdf,ncep_control_cdf,cmc_control_cdf,&
   !            ecmwf_ensemble_cdf,ncep_ensemble_cdf,cmc_ensemble_cdf,gem_cdf)
   write(6,*)'Calling read_cdf_netcdf'
   CALL read_cdf_netcdf(nxa,nya,npct,nens_cmc,ipcpvar,thresh,precip_anal_cdf,&
                        ecmwf_deterministic_cdf,ncep_control_cdf,cmc_control_cdf,&
                        ecmwf_ensemble_cdf,ncep_ensemble_cdf,cmc_ensemble_cdf,&
                        gem_cdf, rlatsa, rlonsa)

   PRINT *,'ipcpvar = ',ipcpvar
   if(ipcpvar.eq.0)then ! POP

      ! Perform CDF-based bias correction.  Here is where quantile mapping occurs.
      ! New (May 2016) version introduces use of more forecast grid points (original + 8 surrounding ones)
      ! and adds stochasticism to the quantile mapping to account for under-spread nature of ensemble.
      write(6,*)'Calling control_cdf_biascorrection_x9'

      print *,'precip_anal_cdf(nxa/2,nya/2,:) = ', precip_anal_cdf(nxa/2,nya/2,:)
      print *,'ncep_ensemble_cdf(nxa/2,nya/2,:) = ', ncep_ensemble_cdf(nxa/2,nya/2,:)
      print *,'before: ncep_ensemble_ccpa(nxa/2,nya/2,:) = ',ncep_ensemble_ccpa(nxa/2,nya/2,:)
      CALL control_cdf_biascorrection_x9(nxa,nya,npct,nstride,nens_cmc,&
           nens_ncep,stdran,icm_c_ccpa,inc_c_ccpa,inc_e_ccpa,icm_e_ccpa,&
           thresh,conusmask,precip_anal_cdf,ncep_control_cdf,&
           cmc_control_cdf,ncep_ensemble_cdf,cmc_ensemble_cdf,&
           cmc_control_ccpa,ncep_control_ccpa,&
           cmc_ensemble_ccpa,ncep_ensemble_ccpa, &
           cmc_control_ccpa_x9,ncep_control_ccpa_x9,&
           cmc_ensemble_ccpa_x9,ncep_ensemble_ccpa_x9)

      ! Get probabilities from quantile-mapped (CDF bias corrected) ensemble
      write(6,*)'Calling raw_ensemble_probs_x9'
      CALL raw_ensemble_probs_x9(nxa,nya,nens_cmc,nens_ncep,&
           rthresh,cmc_control_ccpa_x9,ncep_control_ccpa_x9,&
           cmc_ensemble_ccpa_x9,ncep_ensemble_ccpa_x9,icm_c_ccpa,inc_c_ccpa,&
           inc_e_ccpa,icm_e_ccpa,prob_forecast_cdf,&
           prob_forecast_cdf_CMC,prob_forecast_cdf_NCEP)

   elseif (ipcpvar.eq.1)then ! QPF

      ! ---- Generate grand ensemble mean forecast for current model forecast
      write(6,*)'Calling grand_ensemble_mean'
      CALL grand_ensemble_mean(nxa, nya, nens_cpc, nens_ncep, icm_c_ccpa, inc_c_ccpa,&
           inc_e_ccpa, icm_e_ccpa, ncep_ensemble_ccpa, ncep_control_ccpa,&
           cmc_ensemble_ccpa, cmc_control_ccpa, gem_today_before)

      ! ---- CDF bias correct (quantile map) the ensemble-mean forecast
      gem_today_after = gem_today_before
      write(6,*)'Calling cdf_correct'
      !CALL cdf_correct(nxa,nya,npct,nstride,thresh,conusmask,gem_cdf,precip_anal_cdf,&
      CALL cdf_correct(nxa,nya,npct,thresh,conusmask,gem_cdf,precip_anal_cdf,&
            gem_today_after)

      ! ---- Now Savitzky-Golay smooth the deterministic quantile mapped determinstic forecast, 
      !      first determining the weight array to multiply by the forecast.
      write(6,*)'Calling sgolay_2d_weights'
      CALL sgolay_2d_weights(window_size,order, istat, weights)

      ! ---- Determine how much we will weight the Savitzky-Golay smoothed fields vs.
      !      the raw input.  It makes sense to weight the former more in regions
      !      where there is not much topographic variation, and the latter more where
      !      there is.
      istat = 0
      write(6,*)'Calling raw_vs_smoothed_weight_x9'
      CALL raw_vs_smoothed_weight_x9(nxa,nya,topo_eighth,conusmask,raw_weight)

      ! ---- Now do the Savitzky-Golay smoothing, and the blending of raw and S-G smoothed
      !      based on terrain variation
      determ_smoothed = gem_today_after
      write(6,*)'Calling sgolay_smooth_x9'
      CALL sgolay_smooth_x9(nxa,nya,ipcpvar,determ_smoothed,gem_today_before,&
           conusmask,weights,window_size,order,istat)
      determ_smoothed = raw_weight*gem_today_after + (1.-raw_weight)*determ_smoothed

   else
      write(6,*)' **** ipcpvar = ',ipcpvar,' not a valid option.  Stopping.'
      stop
   endif

   deallocate(precip_anal_cdf,ecmwf_deterministic_cdf,ncep_control_cdf,&
              cmc_control_cdf,ecmwf_ensemble_cdf,ncep_ensemble_cdf,&
              cmc_ensemble_cdf)

end if ! icdfbcorr=1

write(6,fmt='(A)')' Model Stats (1/8 deg. grid):'
write(6,fmt='(4(A10,1X))')'MODEL','MIN','MAX','MEAN'
write(6,fmt='(A10,1X,3(F10.5,1X))')'NCEP',minval(ncep_control_ccpa),maxval(ncep_control_ccpa),sum(ncep_control_ccpa)/(nxa*nya)
write(6,fmt='(A10,1X,3(F10.5,1X))')'NCEP ENS',minval(ncep_ensemble_ccpa),maxval(ncep_ensemble_ccpa),sum(ncep_ensemble_ccpa)/(nxa*nya)
write(6,fmt='(A10,1X,3(F10.5,1X))')'CMC',minval(cmc_control_ccpa),maxval(cmc_control_ccpa),sum(cmc_control_ccpa)/(nxa*nya)
write(6,fmt='(A10,1X,3(F10.5,1X)/)')'CMC ENS',minval(cmc_ensemble_ccpa),maxval(cmc_ensemble_ccpa),sum(cmc_ensemble_ccpa)/(nxa*nya)

!
if(ipcpvar.eq.0)then ! POP

   ! ---- Now Savitzky-Golay smooth the deterministic quantile mapped determinstic forecast, 
   !      first determining the weight array to multiply by the forecast.
   write(6,*)'Calling sgolay_2d_weights'
   CALL sgolay_2d_weights(window_size,order,istat,weights)

   ! ---- Determine how much we will weight the Savitzky-Golay smoothed fields vs.
   !      the raw input.  It makes sense to weight the former more in regions
   !      where there is not much topographic variation, and the latter more where
   !      there is.
   istat = 0
   write(6,*)'Calling raw_vs_smoothed_weight_x9'
   CALL raw_vs_smoothed_weight_x9(nxa,nya,topo_eighth,conusmask,raw_weight)

   ! ---- Smooth via Savitzky-Golay 
   prob_forecast=prob_forecast_cdf  ! prob_forecast_cdf is raw/jittered MME outside of conusmask
   write(6,*)'Calling sgolay_smooth_x9'
   CALL sgolay_smooth_x9(nxa,nya,ipcpvar,prob_forecast,prob_forecast_cdf,&
                         conusmask,weights,window_size,order,istat)
   prob_forecast=raw_weight*prob_forecast_cdf+(1.-raw_weight)*prob_forecast

   ! ---- Final pass through the final blended probability forecast to check for
   !      "bad" values.
   !      1) Set negative values to 0.0.
   !      2) Check for NaN? Necessary?
   !      
   do j=1,nya
      do i=1,nxa
         if(prob_forecast(i,j).lt.0.0) prob_forecast(i,j)=0.0
         if(isnan(prob_forecast(i,j))) prob_forecast(i,j)=9999.0
      end do
   end do

   write(6,*)'Final prob_forecast max,min = ', maxval(prob_forecast), minval(prob_forecast)

   ! ---- Compute the climatological probability
   PRINT *,'Calling compute_prob_climatology'
   CALL compute_prob_climatology(nxa,nya,ndays3mo,precip_anal_fine,ihavedata,rthresh,climo_prob)
   print *,'climo_prob(1:nxa:5,nya/2) = ',climo_prob(1:nxa:5,nya/2)

   CALL GET_ENVIRONMENT_VARIABLE('FORT43',outfile)
   write(6,*)'Writing data to ',TRIM(outfile)
   OPEN(unit=43,file=outfile,status='unknown',form='unformatted')
   WRITE(43)nxa,nya
   WRITE(43)prob_forecast
   WRITE(43)prob_forecast_raw
   WRITE(43)climo_prob
   WRITE(43)rlonsa
   WRITE(43)rlatsa
   WRITE(43)conusmask
   WRITE(43)prob_forecast_cdf
   WRITE(43)prob_forecast_raw_CMC
   WRITE(43)prob_forecast_raw_NCEP
   WRITE(43)prob_forecast_cdf
   WRITE(43)prob_forecast_cdf_CMC
   WRITE(43)prob_forecast_cdf_NCEP
   WRITE(43)ncep_control_ccpa
   WRITE(43)ncep_ensemble_ccpa
   WRITE(43)cmc_control_ccpa
   WRITE(43)cmc_ensemble_ccpa
   CLOSE(43)

elseif(ipcpvar.eq.1)then

   CALL GET_ENVIRONMENT_VARIABLE('FORT43',outfile)
   write(6,*)'Writing data to ',TRIM(outfile)
   OPEN(unit=43,file=outfile,status='unknown',form='unformatted')
   WRITE(43)nxa,nya
   WRITE(43)gem_today_before
   WRITE(43)gem_today_after
   WRITE(43)determ_smoothed
   WRITE(43)conusmask
   WRITE(43)rlonsa
   WRITE(43)rlatsa
   CLOSE(43) 

endif

! ---- Write Final Blended forecasts to GRIB2 file
CALL GET_ENVIRONMENT_VARIABLE('FORT40',griboutfile)
CALL BAOPEN(1,TRIM(griboutfile),ierr)
if(ipcpvar.eq.0)CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,prob_forecast,ipcpvar,1)
if(ipcpvar.eq.1)CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,determ_smoothed,ipcpvar,1)
CALL BACLOSE(1,ierr)
write(6,*)'GRIB2 File written. ierr = ',ierr

! ---- For POP, we also want the model-specific blended forecasts
if(ipcpvar.eq.0)then

   ! Write raw NCEP field
   CALL GET_ENVIRONMENT_VARIABLE('FORT41',griboutfile_ncep)
   CALL BAOPEN(1, TRIM(griboutfile_ncep), ierr)
   CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,prob_forecast_raw_NCEP,ipcpvar,1)
   CALL BACLOSE(1, ierr)

   ! Write raw CMC field
   CALL GET_ENVIRONMENT_VARIABLE('FORT42',griboutfile_cmc)
   CALL BAOPEN(1, TRIM(griboutfile_cmc), ierr)
   CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,prob_forecast_raw_CMC,ipcpvar,1)
   CALL BACLOSE(1, ierr)
endif

DEALLOCATE(rlonsa,rlatsa,&
           precip_anal_fine,prob_forecast,prob_forecast_raw,climo_prob,&
           cmc_control_ccpa,ncep_control_ccpa,&
           cmc_ensemble_ccpa,ncep_ensemble_ccpa,topo_eighth,&
           prob_forecast_raw_ecmwf,prob_forecast_raw_CMC,prob_forecast_raw_NCEP,&
           gem_today_before,gem_today_after,determ_smoothed,raw_weight,&
           cmc_control_ccpa_x9,ncep_control_ccpa_x9,cmc_ensemble_ccpa_x9,&
           ncep_ensemble_ccpa_x9,stat=ios)

! ---- Removed ECMWF arrays from above DEALLOCATE
!   ecmwf_deterministic,ecmwf_ensemble,ecmwf_deterministic_ccpa,ecmwf_ensemble_ccpa,

write(6,*)'Deallocation Status = ',ios
write(6,*)'Done!'

END PROGRAM blend_precip_downscale
