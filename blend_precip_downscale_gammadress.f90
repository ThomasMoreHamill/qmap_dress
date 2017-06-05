PROGRAM blend_precip_downscale_gammadress

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
!    some stochasticism, i.e., best-member dressing, with fitted Gamma distributions
!          A last step is a location-dependent smoothing of the probability field
!    where the smoothing is less aggressive in areas with variabile terrain, e.g.,
!    the western US.
!
! QPF.   Simpler.  A quantile mapping of the ensemble mean to the analyzed distribution.
!    This tends to sharpen up the forecasts. Smoothing follows. 
!
! coded by: Tom Hamill, NOAA, (303) 497-3060, tom.hamill@noaa.gov
!               August-September 2014
!           MDL adaptations (2014-early 2016) by Chrissy Finan and Eric Engle
!           Additional changes by Tom Hamill, May 2016, Dec 2016
USE netcdf

!INTEGER, PARAMETER  :: nxa = 2145 ! number of 1/8-deg analysis grid pts in the x-dir
!INTEGER, PARAMETER  :: nya = 1597 ! number of 1/8-deg analysis grid pts in the y-dir

INTEGER, PARAMETER  :: nxa = 515 ! number of expanded 1/8-deg analysis grid pts in the x-dir
INTEGER, PARAMETER  :: nya = 262 ! number of expanded 1/8-deg analysis grid pts in the y-dir
INTEGER, PARAMETER  :: nx47 = 297 ! number of grid pts in x-dir for 47-km polar stereo grid
INTEGER, PARAMETER  :: ny47 = 169 !      "                y-dir
INTEGER, PARAMETER  :: nx95 = 137 ! number of grid pts in x-dir for 95-km polar stereo grid
INTEGER, PARAMETER  :: ny95 = 81  !      "                y-dir
INTEGER, PARAMETER  :: nxs = 464  ! number of grid pts in x-dir for 1/8-deg analysis grid
INTEGER, PARAMETER  :: nys = 224  ! number of grid pts in y-dir for 1/8-deg analysis grid
!INTEGER, PARAMETER  :: ncount = 625 ! max number of 1/8-deg grid pts that could be associated with fcst grid pt
INTEGER, PARAMETER  :: ncount = 70 ! max number of 1/8-deg grid pts that could be associated with fcst grid pt
INTEGER, PARAMETER  :: nens_cmc   = 20 ! number of CMC perturbed ensemble members
INTEGER, PARAMETER  :: nens_ncep  = 20 ! number of NCEP perturbed ensemble members
INTEGER, PARAMETER  :: ndays = 31 ! number of forecast days in data set (00 UTC Jan 2014 forecasts)
INTEGER, PARAMETER  :: ndays3mo = 91*14 ! ± 45 days from center of month x 12 years (2002-2016)
INTEGER, PARAMETER  :: npct = 91 ! number of thresholds for CDFs
INTEGER, PARAMETER  :: nmult_other = 1 ! number of replicates of CMC and NCEP deterministic
INTEGER, PARAMETER  :: nmembers = nens_cmc + nens_ncep + nmult_other*2 ! total # mbrs, NCEP+CMC ens+determ.
INTEGER, PARAMETER  :: n_amounts = 68 ! the number of precip values where best-member dressing stats tallied.
INTEGER, PARAMETER  :: n_climocats = 8 ! the categories in the fraction_zero_fclimpop array, used for
                       ! setting the fraction of zero observed in situations where the forecast is zero.
                       ! This now varies by region according to the climatological POP

INTEGER, PARAMETER  :: window_size = 9 ! used in Savitzky-Golay smoothing of prob field
INTEGER, PARAMETER  :: order = 3 ! fit 3rd-order polynomial in Savitzky Golay smoother
REAL, PARAMETER :: rthresh_pop = 0.4  ! POP threshold in mm

INTEGER :: nstride ! stride length between grid pts used in 3x3 array later

REAL, DIMENSION(npct) :: thresh ! the precip amount thresholds for CDFs

CHARACTER*2 chh
CHARACTER*256 panal_infile ! name of file with coarse and fine-res precip analyses, 2002-2013
CHARACTER*256 closest_infile ! name of netcdf file with histogram of closest sorted member
CHARACTER*256 dressing_infile ! name of netcdf file with dressing gamma dist fittted stats
CHARACTER*256 pfcst_infile ! name of file with netCDFized forecast data, deterministic and ensemble
CHARACTER*256 iccpa_infile ! has the lists of grid pts on CCPA grid assoc'd with polar-stereo grid
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
CHARACTER*1 cdiagnostics   ! Flag for whether extra diagnostics to be written out (y|Y or n|N).
CHARACTER*1 cdressing      ! Flag for whether do dress with gamma-distributed noise (y|Y or n|N).
CHARACTER*5 cthresh        ! 'POP','1mm','2p5mm','5mm','10mm','25mm','50mm'
CHARACTER*20 cfield

INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! inherited from CCPA data set
REAL*4, DIMENSION(window_size,window_size) :: weights ! weightings for Savitzky-Golay smooth

REAL, DIMENSION(nmembers) :: closest_weight ! weight to apply to sorted ensemble member
REAL, DIMENSION(n_climocats) :: fraczero_fclimpop ! fraction with zero precip as f(climo POP) 
REAL, DIMENSION(n_climocats-1) :: climo_pop_thresholds ! associated boundaries between fraczero_fclimpop elements

! ---- 1/8 deg. Lat/Lon arrays (i.e. CCPA grid)

REAL, ALLOCATABLE, DIMENSION(:,:) :: cmc_control_ccpa ! cmc control precip forecast  on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:) :: climo_prob ! climatological prob of exceeeding threshold from analyses
REAL, ALLOCATABLE, DIMENSION(:,:) :: climo_pop ! climatological prob of nonzero precip
REAL, ALLOCATABLE, DIMENSION(:,:) :: ncep_control_ccpa ! ncep control precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:) :: gem_today_before ! grand ens mean for today's forecast before bias corr
REAL, ALLOCATABLE, DIMENSION(:,:) :: gem_today_after ! grand ensemble mean for today's after bias corr
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_ensemble_ccpa ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble_ccpa ! ncep ens precip forecast on 1/8-deg ccpa grid

! ---- x9 arrays

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_control_ccpa_x9 ! cmc control precip forecast  on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_control_ccpa_x9 ! ncep control precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_ccpa_x9 ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ncep_ensemble_ccpa_x9 ! ncep ens precip forecast on 1/8-deg ccpa grid

REAL, ALLOCATABLE, DIMENSION(:,:) :: rlonsa ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlatsa ! precip analysis grid lat/lons

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: precip_anal_fine ! precipitation analysis on fine grid

REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast  ! output downscaled probabilistic forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw  ! output raw ensemble probabilisitic forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_ECMWF ! output raw ECMWF ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_NCEP ! output raw NCEP ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_unsmoothed ! MME unsmoothed ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_CMC ! output raw CMC ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_smoothed ! output of S-G smoother
REAL, ALLOCATABLE, DIMENSION(:,:) :: topo_eighth
REAL, ALLOCATABLE, DIMENSION(:,:) :: raw_weight
REAL, ALLOCATABLE, DIMENSION(:,:) :: determ_smoothed

REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: precip_anal_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_deterministic_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_control_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_control_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_ensemble_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble_cdf
REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_cdf ! ensemble not exchangeable, track CDF for each mbr.
REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: gem_cdf ! grand ensemble mean CDF

REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf  ! output cdf bcorr ensemble prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_NCEP ! output cdf NCEP ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_CMC

REAL, ALLOCATABLE, DIMENSION(:) :: fraczero ! fraction of dressed dist'n with zero precip
REAL, ALLOCATABLE, DIMENSION(:) :: gamma_shape ! gamma dressing shape parameter (alpha)
REAL, ALLOCATABLE, DIMENSION(:) :: gamma_scale ! gamma dressing scale parameter (beta)
REAL, ALLOCATABLE, DIMENSION(:) :: ramt ! precipitation amount gamma dist parameters determined at
REAL, ALLOCATABLE, DIMENSION(:) :: closest_histogram

INTEGER :: itot_c     ! flag for total model datasets loaded
INTEGER :: icm_c_ccpa ! flag for cmc control data on 1/8-deg ccpa grid
INTEGER :: icm_e_ccpa ! flag for cmc ensemble data on 1/8-deg ccpa grid
INTEGER :: iec_d_ccpa ! flag for ecmwf deterministic forecast data on 1/8-deg ccpa grid
INTEGER :: iec_e_ccpa ! flag for ecmwf ensemble on 1/8-deg ccpa grid
INTEGER :: inc_c_ccpa ! flag for ncep control on 1/8-deg ccpa grid
INTEGER :: inc_e_ccpa ! flag for ncep ensemble on 1/8-deg ccpa grid

INTEGER :: ierr      ! return variable for BAOPEN
INTEGER :: ios       ! return variable for Fortran I/O, Allocation statements
INTEGER :: ifcstint  ! forecast interval. Set depending on cpcpvar.
INTEGER :: ipcpvar   ! 0 = pop; 1 = qpf
INTEGER :: icdfbcorr ! 0 = no; 1 = yes

INTEGER :: iyyyymmddhh,jyyyymmddhh
INTEGER :: iyear,imo,iday,ihour,idoy ! Parsed date variables from iyyyymmddhh
INTEGER :: jyear,jmo,jday,jhour,jdoy ! Parsed date variables for valid date

! ---- Print some parameters

WRITE (6,*) ' Parameters for blend_precip_downscale :'
WRITE (6,100) nxa,nya,nxs,nys,ndays3mo
100 FORMAT('nxa: ',i0.1/,'nya: ',i0.1/,'nxs: ',&
    i0.1/,'nys: ',i0.1/,'ndays3mo: ',i0.1/)

! --- Via command line, read in the input year/mo/day/hr and the forecast resolution
!     we're working with.  Process date to determine the day of the month as an INTEGER

! **** ERIC, note getting rid of cdfbcorr input on command line;
!      No need for option, we will always do this.  This getarg replaced with
!      a variable for setting whether diagnostic data for me is output  

CALL getarg(1,cyyyymmddhh) ! input year month day hour of initial condition, 'yyyymmddhh' format
CALL getarg(2,cpcpvar) ! Precipitation Variable to run for ("pop" or "qpf")
CALL getarg(3,cleade) ! forecast lead time for beginning of precip accum period, hours, e.g.'060'
CALL getarg(4,cdiagnostics) ! "Y"/"N" for writing extra output diagnostics
CALL getarg(5,cthresh) ! Precip threshold for which to estimate exceedance probabilities
CALL getarg(6,cdressing) ! "Y" or "y" to dress with Gamma-distributed random noise; when

!  doing a run to determine best member statistics (make sure cdiagnostics = 'y' then
!  then we want to set this to "N")
!   
WRITE (6,*) ' Command line arguments:'
WRITE (6,110) cyyyymmddhh,cpcpvar,cleade,cdiagnostics,cthresh, cdressing
110 FORMAT (1x,'cyyyymmddhh: ',A/1x,'cpcpvar: ',A/1x,&
    'cleade: ',A/1x,'cdiagnostics: ',A/1x,'cthresh: ',A/1x,'cdressing :',A/)

! --- **************************************************************************
!     ERIC: you'll probably want to associate this filename in the script that 
!     runs the blend with a unit number of your choice
! --- **************************************************************************

IF (TRIM(cthresh) .eq. 'POP' .or. TRIM(cthresh) .eq. 'pop') THEN
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
    WRITE (6,*) ' **** Invalid threshold value: ',cthresh
    STOP
END IF

! ---- Convert character based variables from command line to INTEGERs

READ (cyyyymmddhh,'(i10)') iyyyymmddhh
READ (cleade,'(i3)') ileade

! ---- Set ifcstint (forecast accumulatio interval) and 
!      ipcpvar according to cpcpvar

IF (cpcpvar .eq. "pop" .or. cpcpvar .eq. "POP") THEN
    ifcstint = 12
    ipcpvar = 0
ELSE IF (cpcpvar .eq. "qpf" .or. cpcpvar .eq. "QPF") THEN
    ifcstint = 6
    ipcpvar = 1
ENDIF
ileadb = ileade - ifcstint

! ---- Set nstride as a function of forecast lead time

nstride = NINT(3. + 8.*ileade/168.)
WRITE (6,*)'nstride = ',nstride

! ---- Parse the initializtion date; determine the valid hour.  This is dependent on the
!      precip variable (cpcpvar), the model initialization (iyyyymmddhh), and forecast
!      ending hour (ileade).

iendhour  = 0
CALL doy(iyyyymmddhh, iyear, imo, iday, ihour, idoy)
CALL updat(iyyyymmddhh, ileade, jyyyymmddhh)
CALL doy(jyyyymmddhh, jyear, jmo, jday, jhour, jdoy)
iendhour = jhour
WRITE (6,*) 'Model Initialization: ', iyyyymmddhh
WRITE (6,*) 'Forecast Projection Ending: ', ileade
WRITE (6,*) 'Forecast Valid Date/Hour: ', jyyyymmddhh, iendhour
WRITE (6,*) 'Precipitation Element Flag: ', ipcpvar

! ---- Get filenames from environment

CALL get_environment_variable('FORT20', panal_infile)
CALL get_environment_variable('FORT25', pmask_infile)
WRITE (6,fmt='(A,A)')'CCPA NetCDF Archive file: ',trim(panal_infile)
WRITE (6,fmt='(A,A)')'CCPA Mask file: ',trim(pmask_infile)

! ---- Allocate dynamic arrays

WRITE (6,*) 'Allocating dynamic arrays...'

ALLOCATE(cmc_control_ccpa(nxa,nya))
ALLOCATE(ncep_control_ccpa(nxa,nya))
ALLOCATE(cmc_ensemble_ccpa(nxa,nya,nens_cmc))
ALLOCATE(ncep_ensemble_ccpa(nxa,nya,nens_ncep))
ALLOCATE(rlonsa(nxa,nya),rlatsa(nxa,nya))
ALLOCATE(precip_anal_fine(nxa,nya,ndays3mo))
ALLOCATE(prob_forecast(nxa,nya))
ALLOCATE(prob_forecast_smoothed(nxa,nya))
ALLOCATE(prob_forecast_unsmoothed(nxa,nya))
ALLOCATE(prob_forecast_raw(nxa,nya))
ALLOCATE(prob_forecast_raw_ECMWF(nxa,nya))
ALLOCATE(prob_forecast_raw_CMC(nxa,nya))
ALLOCATE(prob_forecast_raw_NCEP(nxa,nya))
ALLOCATE(prob_forecast_cdf(nxa,nya))
ALLOCATE(prob_forecast_cdf_CMC(nxa,nya))
ALLOCATE(prob_forecast_cdf_NCEP(nxa,nya))
ALLOCATE(climo_prob(nxa,nya))
ALLOCATE(climo_pop(nxa,nya))
ALLOCATE(topo_eighth(nxa,nya))
ALLOCATE(gem_today_before(nxa,nya))
ALLOCATE(gem_today_after(nxa,nya))
ALLOCATE(determ_smoothed(nxa,nya))
ALLOCATE(raw_weight(nxa,nya))
ALLOCATE(cmc_control_ccpa_x9(9,nxa,nya))
ALLOCATE(ncep_control_ccpa_x9(9,nxa,nya))
ALLOCATE(cmc_ensemble_ccpa_x9(9,nxa,nya,nens_cmc))
ALLOCATE(ncep_ensemble_ccpa_x9(9,nxa,nya,nens_ncep))
ALLOCATE(precip_anal_cdf(nxa,nya,npct))
ALLOCATE(ecmwf_deterministic_cdf(nxa,nya,npct))
ALLOCATE(ncep_control_cdf(nxa,nya,npct))
ALLOCATE(cmc_control_cdf(nxa,nya,npct))
ALLOCATE(ecmwf_ensemble_cdf(nxa,nya,npct))
ALLOCATE(ncep_ensemble_cdf(nxa,nya,npct))
ALLOCATE(cmc_ensemble_cdf(nxa,nya,nens_cmc,npct))
ALLOCATE(gem_cdf(nxa,nya,npct))
ALLOCATE(fraczero(n_amounts))
ALLOCATE(gamma_shape(n_amounts))
ALLOCATE(gamma_scale(n_amounts))
ALLOCATE(ramt(n_amounts))
ALLOCATE(closest_histogram(nmembers))

! ---- Read in the 1/8-degree precipitation analyses.  For a given month,
!      load in only the precip analyses for the surrounding months,
!      e.g., for January load up Dec-Jan-Feb.

! ****  Eric, all this precip_anal_file data is read in to calculate the
! climatological probability.  We could precalculate this for each month
! and save a lot of I/O time and memory.  Let me know if you wish to pursue this.

WRITE (6,*) 'Calling read_precip_analyses_x9'
CALL read_precip_analyses_x9(nxa, nya, nxs, nys, ndays3mo, ipcpvar, &
    iyyyymmddhh, iendhour, panal_infile, pmask_infile, precip_anal_fine, &
    conusmask)

! ---- Read precipitation foreasts for all models/ensemble 
!      members from TDLPACK Vector files.
!      These files contains the 1/8 deg. grid.

WRITE (6,*) 'Calling read_forecasts_ccpa_tdlp'
CALL read_forecasts_tdlp_vect(nxa, nya, nens_cmc, nens_ncep, ipcpvar,&
    iyyyymmddhh, ileade, cmc_control_ccpa, ncep_control_ccpa, &
    cmc_ensemble_ccpa, ncep_ensemble_ccpa, icm_c_ccpa, inc_c_ccpa,&
    icm_e_ccpa, inc_e_ccpa, ier)

! ---- Check to make sure data is avaiable

itot_c = icm_c_ccpa + inc_c_ccpa + icm_e_ccpa + inc_e_ccpa

IF (itot_c .eq. 0) THEN
    WRITE (6,*) '**** NO MODEL DATA AVAILABLE FOR DOWNSCALE. PROGRAM STOPPING.'
    CALL w3tage('BLEND_PRECIP_DOWNSCALE')
    STOP
ELSE
    IF (icm_c_ccpa .eq. 0) THEN 
        WRITE (6,*) '**** CMC CONTROL DATA MISSING. ARRAY SET TO MISSING.'
        cmc_control_ccpa = -99.99
    ELSE
        WRITE (6,*) 'CMC dataset loaded'
    ENDIF
    
    IF (inc_c_ccpa .eq. 0) THEN
        WRITE (6,*) '**** NCEP CONTROL DATA MISSING. ARRAY SET TO MISSING.'
        ncep_control_ccpa = -99.99
    ELSE
        WRITE (6,*)'GFS dataset loaded'
    END IF

    IF (icm_e_ccpa.eq.0) THEN
        WRITE (6,*) '**** CMC ENSEMBLE DATA MISSING. ARRAY SET TO MISSING.'
        cmc_ensemble_ccpa = -99.99
    ELSE
        WRITE (6,*)'CMCE dataset loaded'
    ENDIF

    IF (inc_e_ccpa .eq. 0) THEN
        WRITE (6,*) '**** NCEP ENSEMBLE DATA MISSING. ARRAY SET TO MISSING.'
        ncep_ensemble_ccpa = -99.99
    ELSE
        WRITE (6,*)'GEFS dataset loaded'
    ENDIF
END IF

WRITE (6,120)
120 format('MODEL AVAILABILITY:'/'     1/8 DEG.')
WRITE (6,130) 'CMC',icm_c_ccpa
WRITE (6,130) 'CMCE',icm_e_ccpa
WRITE (6,130) 'GFS',inc_c_ccpa
WRITE (6,130) 'GEFS',inc_e_ccpa
130 FORMAT(A4,1X,I2,4X,I2)

! ---- Read in the 1/8 degree topography data (created by 
!      create_eighth_degree_conus_elev.f90).   This will be used to 
!      determine how much smoothing of the prob forecast to do, with more
!      smoothing applied in the flat regions of the CONUS than in mountainous
!      regions, since we want to preserve prob detail related to terrain features.

CALL get_environment_variable('FORT24', ptopo_infile)
WRITE (6,fmt='(A,A)')'1/8 Deg. Terrain File: ',trim(ptopo_infile)
CALL check(nf90_open(ptopo_infile, NF90_NOWRITE, netid))
CALL check(nf90_inq_varid(netid, "MTERH_surface", ivar))
CALL check(nf90_get_var(netid, ivar, topo_eighth, &
    start=(/1,1,1/), count=(/nxa,nya,1/)))
CALL check(nf90_close(netid))

! ---- read in cumulative distribution function data

WRITE (6,*)'Calling read_cdf_netcdf'
CALL read_cdf_netcdf(nxa,nya,npct,nens_cmc,ipcpvar,thresh,precip_anal_cdf,&
     ecmwf_deterministic_cdf,ncep_control_cdf,cmc_control_cdf,&
     ecmwf_ensemble_cdf,ncep_ensemble_cdf,cmc_ensemble_cdf,gem_cdf, &
     rlatsa, rlonsa)

IF (ipcpvar .eq. 0) THEN ! probabilistic forecast
    
    ! ---- For purposes of having a baseline for comparison, generate an ensemble probability
    !      simply from the relative frequency.

    WRITE (6,*)'Calling raw_ensemble_probs'
    CALL raw_ensemble_probs(nxa, nya,nens_cmc,nens_ncep,rthresh,cmc_control_ccpa,&
        ncep_control_ccpa,cmc_ensemble_ccpa,ncep_ensemble_ccpa,&
        icm_c_ccpa,inc_c_ccpa,inc_e_ccpa,icm_e_ccpa,prob_forecast_raw,&
        prob_forecast_raw_CMC,prob_forecast_raw_NCEP)

    ! ---- Compute the climatological probability
        
    PRINT *, 'Calling compute_prob_climatology'
    CALL compute_prob_climatology(nxa,nya,ndays3mo,&
        precip_anal_fine,conusmask,rthresh,climo_prob)
    IF (rthresh .ne. rthresh_pop) THEN
       CALL compute_prob_climatology(nxa,nya,ndays3mo,&
           precip_anal_fine,conusmask,rthresh_pop,climo_pop)
    ELSE
       climo_pop = climo_prob
    ENDIF

    ! ---- Here is where deterministic quantile mapping occurs, using 
    !      original + 8 surrounding grid points
    
    WRITE (6,*) 'Calling control_quantile_mapping_x9'
    CALL control_quantile_mapping_x9(nxa, nya, npct, nstride, nens_cmc, &
        nens_ncep, icm_c_ccpa, inc_c_ccpa, inc_e_ccpa, icm_e_ccpa, &
        thresh, conusmask, precip_anal_cdf, ncep_control_cdf, &
        cmc_control_cdf, ncep_ensemble_cdf, cmc_ensemble_cdf, &
        cmc_control_ccpa, ncep_control_ccpa, cmc_ensemble_ccpa, &
        ncep_ensemble_ccpa, cmc_control_ccpa_x9, ncep_control_ccpa_x9, &
        cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9)

    ! ---- open the file that has the best member dressing parameters 
    !      for the Gamma distribution

    CALL get_environment_variable('FORT86', dressing_infile)
    PRINT *,'dressing file = ',TRIM(dressing_infile)
    CALL read_dressing_parameters_v2 (n_amounts, n_climocats, &
        dressing_infile, ramt, fraczero, fraczero_fclimpop, &
        climo_pop_thresholds, gamma_shape, gamma_scale)

    ! ---- Read in the histogram that describes the weights to apply to the sorted ensemble

    CALL get_environment_variable('FORT87', closest_infile)
    WRITE (6,fmt='(A,A)')' Closest best-member histogram file ',trim(closest_infile)
    CALL read_closest_histogram(nmembers, closest_infile, closest_histogram)

    ! ---- Apply dressing and weighting of ensemble members based on closest_histogram.
    !      Get probabilities from quantile-mapped, dressed, weighted  ensemble.
      
    WRITE (6,*) 'Calling ensemble_probs_x9_dressweight'
    !CALL ensemble_probs_x9_weighted (nxa, nya, nens_cmc, nens_ncep, nmult_other, &
    !    nmembers, n_amounts, cdressing, rthresh, ramt, fraczero, gamma_shape, gamma_scale, &
    !    closest_histogram, conusmask, cmc_control_ccpa_x9, ncep_control_ccpa_x9, &
    !    cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, icm_c_ccpa, inc_c_ccpa, &
    !    inc_e_ccpa, icm_e_ccpa, prob_forecast_cdf, prob_forecast_cdf_CMC, &
    !    prob_forecast_cdf_NCEP, prob_forecast_unsmoothed)
    CALL ensemble_probs_x9_dressweight (nxa, nya, nens_cmc, nens_ncep, &
        nmult_other, nmembers, n_amounts, n_climocats, cdressing, rthresh, &
        ramt, fraczero, fraczero_fclimpop, climo_pop_thresholds, &
        gamma_shape, gamma_scale, closest_histogram, climo_pop, &
        conusmask, cmc_control_ccpa_x9, ncep_control_ccpa_x9, &
        cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, icm_c_ccpa, inc_c_ccpa, &
        inc_e_ccpa, icm_e_ccpa, prob_forecast_cdf, prob_forecast_cdf_CMC, &
        prob_forecast_cdf_NCEP, prob_forecast_unsmoothed)

    ! ---- Determine how much we will weight the Savitzky-Golay smoothed fields vs.
    !      the raw input.  It makes sense to weight the former more in regions
    !      where there is not much topographic variation, and the latter more where
    !      there is.
            
    WRITE (6,*) 'Calling sgolay_2d_weights'
    CALL sgolay_2d_weights(window_size, order, istat, weights)
    istat = 0
    WRITE (6,*) 'Calling raw_vs_smoothed_weight_x9'
    CALL raw_vs_smoothed_weight_x9(nxa, nya, topo_eighth, conusmask, raw_weight)

    ! ---- Smooth via Savitzky-Golay
        
    WRITE (6,*)'Calling sgolay_smooth_x9'
    prob_forecast_smoothed =  prob_forecast_unsmoothed
    PRINT *,'before sgolay_smooth_x9'
    CALL sgolay_smooth_x9(nxa, nya, ipcpvar, prob_forecast_smoothed, &
         prob_forecast_unsmoothed, conusmask, weights, window_size, order, istat)
    prob_forecast = raw_weight*prob_forecast_unsmoothed + &
         (1.-raw_weight)*prob_forecast_smoothed

    PRINT *,'final probabilities'
    DO jya = nya/3, nya/3
       DO ixa = nxa/5, nxa/2,5
          print 317,'ixa,jya,prob = ',ixa, jya, &
               prob_forecast(ixa,jya)
317       format(a33,2(i3,1x),4(f8.5,1x))
       END DO
    END DO


    ! ---- Final pass through the final blended probability forecast to check for
    !      "bad" values.
    !      1) Set negative values to 0.0.
    !      2) Check for NaN? Necessary?
    !
    DO j = 1, nya
        DO i = 1,nxa
            IF (prob_forecast(i,j) .lt. 0.0) prob_forecast(i,j) = 0.0
            IF (isnan(prob_forecast(i,j))) prob_forecast(i,j) = 9999.0
        END DO
    END DO

    ! ---- WRITE out the resulting data if desired

    PRINT *,'cdiagnostics = ',cdiagnostics
    IF (cdiagnostics .eq. 'Y' .or. cdiagnostics .eq. 'y') THEN
        CALL GET_ENVIRONMENT_VARIABLE('FORT43',outfile)
        WRITE (6,*)'Writing data to ',TRIM(outfile)
        OPEN(unit=43,file=TRIM(outfile),status='unknown',form='unformatted')
        WRITE (43) nxa,nya
        WRITE (43) prob_forecast
        WRITE (43) prob_forecast_raw
        WRITE (43) climo_prob
        WRITE (43) rlonsa
        WRITE (43) rlatsa
        WRITE (43) conusmask
        WRITE (43) prob_forecast_cdf
        WRITE (43) prob_forecast_raw_CMC
        WRITE (43) prob_forecast_raw_NCEP
        WRITE (43) prob_forecast_unsmoothed
        WRITE (43) prob_forecast_cdf_CMC
        WRITE (43) prob_forecast_cdf_NCEP
        WRITE (43) prob_forecast_smoothed
        WRITE (43) raw_weight
        !IF (cdiagnostics .eq. 'Y' .or. cdiagnostics .eq. 'y') THEN
        !    WRITE(43)ncep_control_ccpa_x9(5,:,:)
        !    WRITE(43)ncep_ensemble_ccpa_x9(5,:,:,:)
        !    WRITE(43)cmc_control_ccpa_x9(5,:,:)
        !    WRITE(43)cmc_ensemble_ccpa_x9(5,:,:,:)
        !ENDIF
        CLOSE (43)
    ENDIF
    
    ! ---- Write Final Blended forecasts to GRIB2 files

    CALL GET_ENVIRONMENT_VARIABLE('FORT40',griboutfile)
    CALL BAOPEN(1,TRIM(griboutfile),ierr)
    CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,&
        prob_forecast,ipcpvar,1)
    CALL BACLOSE(1, ierr)
    WRITE (6,*) 'GRIB2 File written. ierr = ', ierr
    
    ! ---- For prob. precip, we also want the raw model probabilities

    !      Write raw NCEP field
   
    CALL GET_ENVIRONMENT_VARIABLE('FORT41',griboutfile_ncep)
    CALL BAOPEN(1, TRIM(griboutfile_ncep), ierr)
    CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,&
        prob_forecast_raw_NCEP,ipcpvar,1)
    CALL BACLOSE(1, ierr)

    !      Write raw CMC field
   
    CALL GET_ENVIRONMENT_VARIABLE('FORT42',griboutfile_cmc)
    CALL BAOPEN(1, TRIM(griboutfile_cmc), ierr)
    CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,&
            prob_forecast_raw_CMC,ipcpvar,1)
    CALL BACLOSE(1, ierr)
   
ELSE IF (ipcpvar .eq. 1) THEN ! QPF
        
    ! ---- Generate grand ensemble mean forecast for current model forecast
      
    WRITE (6,*) 'Calling grand_ensemble_mean'
    CALL grand_ensemble_mean(nxa, nya, nens_cpc, nens_ncep, icm_c_ccpa, inc_c_ccpa,&
        inc_e_ccpa, icm_e_ccpa, ncep_ensemble_ccpa, ncep_control_ccpa,&
        cmc_ensemble_ccpa, cmc_control_ccpa, gem_today_before)
    PRINT *,'max, min gem_today_before = ',maxval(gem_today_before), minval(gem_today_before)
    PRINT *,'gem_today_before(:,nya/2) = ',gem_today_before(:,nya/2)

    ! ---- Quantile map the ensemble-mean forecast
      

    PRINT *, 'before_cdf_correct, precip_anal_cdf(479,199,:) = ',&
         precip_anal_cdf(479,199,:)
    gem_today_after = gem_today_before
    WRITE (6,*) 'Calling cdf_correct'
    CALL cdf_correct(nxa, nya, npct, thresh, conusmask, gem_cdf, precip_anal_cdf, &
        gem_today_after)
    PRINT *,'max, min gem_today_after = ',maxval(gem_today_after), minval(gem_today_after)
    PRINT *,'gem_today_after(:,nya/2) = ',gem_today_after(:,nya/2)


    ! ---- Now Savitzky-Golay smooth the deterministic quantile mapped determinstic forecast,
    !      first determining the weight array to multiply by the forecast.
        
    WRITE (6,*)'Calling sgolay_2d_weights'
    CALL sgolay_2d_weights(window_size,order, istat, weights)

    ! ---- Determine how much we will weight the Savitzky-Golay smoothed fields vs.
    !      the raw input.  It makes sense to weight the former more in regions
    !      where there is not much topographic variation, and the latter more where
    !      there is.
        
    istat = 0
    WRITE (6,*)'Calling raw_vs_smoothed_weight_x9'
    CALL raw_vs_smoothed_weight_x9(nxa, nya, topo_eighth, conusmask, raw_weight)

    ! ---- Now do the Savitzky-Golay smoothing, and the blending of raw and S-G smoothed
    !      based on terrain variation
        
    determ_smoothed = gem_today_after
    WRITE (6,*)'Calling sgolay_smooth_x9'
    CALL sgolay_smooth_x9(nxa, nya, ipcpvar, determ_smoothed, gem_today_before,&
        conusmask, weights, window_size, order, istat)

    PRINT *,'max, min determ_smoothed before = ',maxval(determ_smoothed), minval(determ_smoothed)
    PRINT *,'determ_smoothed(:,nya/2) = ',determ_smoothed(:,nya/2)

    determ_smoothed = raw_weight*gem_today_after + (1.-raw_weight)*determ_smoothed
    PRINT *,'max, min determ_smoothed = ',maxval(determ_smoothed), minval(determ_smoothed)
    PRINT *,'determ_smoothed(:,nya/2) = ',determ_smoothed(:,nya/2)

    ! ---- write the blended forecasts and associated diagnostics to a f77_unformatted file
    
    IF (cdiagnostics .eq. 'Y' .or. cdiagnostics .eq. 'y') THEN
        CALL GET_ENVIRONMENT_VARIABLE('FORT43',outfile)
        WRITE (6,*)'Writing data to ',TRIM(outfile)
        OPEN (unit=43, file=outfile, status='unknown', form='unformatted')
        WRITE (43) nxa,nya
        WRITE (43) gem_today_before
        WRITE (43) gem_today_after
        WRITE (43) determ_smoothed
        WRITE (43) conusmask
        WRITE (43) rlonsa
        WRITE (43) rlatsa
        CLOSE (43)
    ENDIF
    
    ! ---- Write Final Blended forecasts to GRIB2 files

    CALL GET_ENVIRONMENT_VARIABLE('FORT40',griboutfile)
    CALL BAOPEN(1,TRIM(griboutfile),ierr)
    CALL WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,determ_smoothed,ipcpvar,1)
    CALL BACLOSE(1, ierr)
    WRITE (6,*) 'GRIB2 File written. ierr = ', ierr

ELSE
    WRITE (6,*)' **** ipcpvar = ',ipcpvar,' not a valid option. Stopping.'
    STOP
END IF

! ---- some final diagnostic print statements

WRITE (6,fmt='(A)')' Model Stats (1/8 deg. grid):'
WRITE (6,fmt='(4(A10,1X))')'MODEL','MIN','MAX','MEAN'
WRITE (6,fmt='(A10,1X,3(F10.5,1X))')'NCEP',minval(ncep_control_ccpa),&
    maxval(ncep_control_ccpa),sum(ncep_control_ccpa)/(nxa*nya)
WRITE (6,fmt='(A10,1X,3(F10.5,1X))')'NCEP ENS',&
    minval(ncep_ensemble_ccpa),maxval(ncep_ensemble_ccpa),&
    sum(ncep_ensemble_ccpa)/(nxa*nya*nens_ncep)
WRITE (6,fmt='(A10,1X,3(F10.5,1X))')'CMC',minval(cmc_control_ccpa),&
    maxval(cmc_control_ccpa),sum(cmc_control_ccpa)/(nxa*nya)
WRITE (6,fmt='(A10,1X,3(F10.5,1X)/)')'CMC ENS',&
    minval(cmc_ensemble_ccpa),maxval(cmc_ensemble_ccpa),&
    sum(cmc_ensemble_ccpa)/(nxa*nya*nens_cmc)

DEALLOCATE(rlonsa, rlatsa, precip_anal_fine, prob_forecast, prob_forecast_raw, &
    climo_prob, cmc_control_ccpa, ncep_control_ccpa, cmc_ensemble_ccpa, &
    climo_pop, ncep_ensemble_ccpa, topo_eighth, prob_forecast_raw_ecmwf, &
    prob_forecast_raw_CMC,prob_forecast_raw_NCEP, gem_today_before, &
    gem_today_after, determ_smoothed, raw_weight, cmc_control_ccpa_x9, &
    ncep_control_ccpa_x9, cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, &
    precip_anal_cdf, ecmwf_deterministic_cdf, ncep_control_cdf, &
    cmc_control_cdf, ecmwf_ensemble_cdf, ncep_ensemble_cdf, &
    cmc_ensemble_cdf, gem_cdf, fraczero, gamma_shape, gamma_scale, &
    ramt, closest_histogram, prob_forecast_smoothed, &
    prob_forecast_unsmoothed, stat=ios)

WRITE (6,*)'Deallocation Status = ',ios
WRITE (6,*)'Done!'

END PROGRAM blend_precip_downscale_gammadress
