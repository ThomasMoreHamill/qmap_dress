PROGRAM generate_gammadressing_stats

USE netcdf

integer, parameter :: nxa = 464  ! number of grid pts in x-dir for 1/8-deg analysis grid
integer, parameter :: nya = 224  ! number of grid pts in y-dir for 1/8-deg analysis grid 
integer, parameter :: nens_ecmwf = 50 ! number of ECMWF perturbed ensemble members
integer, parameter :: nens_cmc   = 20 ! number of CMC perturbed ensemble members
integer, parameter :: nens_ncep  = 20 ! number of NCEP perturbed ensemble members
integer, parameter :: npct = 90 ! number of thresholds for CDFs
INTEGER, PARAMETER :: n_climocats = 8 ! the categories in the fraction_zero_fclimpop 
	! array, used for setting the fraction of zero observed in situations where the forecast
	! zero. This now varies by region according to the climatological POP
INTEGER, PARAMETER :: nthreshes = 68 ! number precip threshold amts to calculate gamma parameters
REAL, PARAMETER :: thresh_low = 2.5
REAL, PARAMETER :: thresh_high = 10.0

integer :: nstride
real :: stdran ! magnitude of noise for perturbing forecast quantile

REAL, DIMENSION(npct) :: thresh ! the precip amount thresholds for CDFs
REAL, DIMENSION(nthreshes) :: gamma_threshes

REAL, DIMENSION(n_climocats-1) :: climo_pop_thresholds ! associated boundaries between fraczero_fclimpop elements


CHARACTER*2 chh, cmm
CHARACTER*3, DIMENSION(12) :: cmonths 
CHARACTER*5 :: center
CHARACTER*256 pclimo_infile ! name of file with climo probs
CHARACTER*256 infile_ecmwf_early ! ecmwf forecast data for today
CHARACTER*256 infile_cmc_early ! cmc forecast data for today
CHARACTER*256 infile_ncep_early ! ncep forecast data for today
CHARACTER*256 infile_ecmwf_late ! ecmwf forecast data for today
CHARACTER*256 infile_cmc_late ! cmc forecast data for today
CHARACTER*256 infile_ncep_late ! ncep forecast data for today
CHARACTER*256 infile
CHARACTER*256 outfile      ! name of flat fortran file with output prob forecasts
CHARACTER*10 cyyyymmddhh   ! year,month,day,hour of initial time of forecast
CHARACTER*3 cleade         ! ending hour of precip forecast accumulation, 3 digits, e.g., '024'
CHARACTER*3 cleadb         ! beginning hour of precip forecast accumulation, 3 digits, e.g., '012'
CHARACTER*5 cthresh        ! 'POP','1mm','2p5mm','5mm','10mm','25mm','50mm'
CHARACTER*3 cmodelcombo  ! E, C, N, EC, EN, NC, or ENC permitted

INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! inherited from CCPA data set

! ---- 1/8 deg. Lat/Lon arrays (i.e. CCPA grid)

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_ensemble_ccpa ! ecmwf ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: cmc_ensemble_ccpa ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble_ccpa ! ncep ens precip forecast on 1/8-deg ccpa grid

! ---- x9 arrays

REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ensemble_ccpa_x9 ! ecmwf ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ecmwf_ensemble_ccpa_x9 ! ecmwf ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_ccpa_x9 ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ncep_ensemble_ccpa_x9 ! ncep ens precip forecast on 1/8-deg ccpa grid

REAL, ALLOCATABLE, DIMENSION(:,:) :: analysis ! precipitation analysis
REAL, ALLOCATABLE, DIMENSION(:,:) :: climo_prob ! climatological event probability
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlonsa ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: rlatsa ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: ensemble_mean ! precip analysis grid lat/lons
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast  ! output downscaled prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw  ! output raw ensemble prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_ECMWF ! output raw NCEP ensemble prob 
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_NCEP ! output raw NCEP ensemble prob 
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_raw_CMC ! output raw NCEP ensemble prob 

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: precip_anal_cdf
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ecmwf_ensemble_cdf
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ncep_ensemble_cdf
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_cdf

REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf  ! output cdf bcorr ensemble prob of precip (prob) forecast
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_NCEP ! output cdf NCEP ensemble prob
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_CMC
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_forecast_cdf_ECMWF

integer :: ierr      ! return variable for BAOPEN
integer :: ios       ! return variable for Fortran I/O, Allocation statements
integer :: ifcstint  ! forecast interval. Set depending on cpcpvar.

integer :: iyyyymmddhh,jyyyymmddhh
integer :: iyear,imo,iday,ihour,idoy ! Parsed date variables from iyyyymmddhh
integer :: jyear,jmo,jday,jhour,jdoy ! Parsed date variables for valid date

! ---- Initialize

DATA cmonths /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

DATA gamma_threshes / &  ! (1:68)
    0.0, 0.05, 0.1, 0.2, 0.3, &
    0.4, 0.5, 0.6, 0.7, 0.8, &
    0.9, 1.0, 1.2, 1.4, 1.6, &
    1.8, 2.0, 2.3, 2.6, 3.0, &
    3.3, 3.6, 4.0, 4.5, 5.0, &
    5.5, 6.0, 6.5, 7.0, 7.5, &
    8.0, 8.5, 9.0, 10.0, 11.0, &
	12.0, 13.0, 14.0, 15.0, 16.0, &
	18.0, 20.0, 22.5, 25.0, 27.5, &
	30.0, 33.0, 36.0, 40.0, 45.0, &
	50.0, 55.0, 60.0, 65.0, 70.0, &
	75.0, 80.0, 85.0, 90.0, 95.0, &
    100.0, 120.0, 140.0, 160.0, 180.0, &
	200.0, 250.0, 300.0 /
	
DATA climo_pop_thresholds /0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21/

! --- Via command line, read in the input year/mo/day/hr and the forecast resolution 
!     we're working with.  Process date to determine the day of the month as an integer

CALL getarg(1,cyyyymmddhh)  ! input year month day hour of initial condition, 'yyyymmddhh' format
CALL getarg(2,cleade)       ! forecast lead time for beginning of precip accum period, hours, e.g.'060'
cmodelcombo = 'ENC'
!CALL getarg(3,cmodelcombo)  ! ! E, C, N, EC, EN, NC, or ENC permitted
cthresh = 'POP'
!CALL getarg(4,cthresh)      ! Precip threshold

cmm = cyyyymmddhh(5:6)
READ (cmm,'(i2)') imonth

write(6,*)' Command line arguments:'
write(6,110)  cyyyymmddhh, cleade, cmodelcombo, cthresh
110 format(1x,'cyyyymmddhh: ',A/1x,'cleade: ',A/1x,'cmodelcombo: ',A/1x,'cthresh: ',A/)

IF (TRIM(cthresh).eq.'POP')THEN
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

READ (cyyyymmddhh,'(i10)') iyyyymmddhh
READ (cleade,'(i3)') ileade
!PRINT *,'iyyyymmddhh, ileade = ',iyyyymmddhh, ileade

! ---- Set ifcstint and ipcpvar according to cpcpvar

ifcstint = 12
ileadb = ileade - ifcstint
!PRINT *,'ileadb = ',ileadb
WRITE (cleadb,'(i3)') ileadb

! ---- Set stdran and nstride as a function of ileade

nstride = nint(3.+8.*ileade/168.)
stdran = 0.25 + real(ileade)/1680.
!write(6,*)'nstride = ',nstride
!write(6,*)'stdran = ',stdran

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

! ---- Allocate dynamic arrays

write(6,*)'Allocating dynamic arrays...'
ALLOCATE (analysis(nxa,nya))
ALLOCATE (ecmwf_ensemble_ccpa(nxa,nya,nens_ecmwf))
ALLOCATE (cmc_ensemble_ccpa(nxa,nya,nens_cmc))
ALLOCATE (ncep_ensemble_ccpa(nxa,nya,nens_ncep))
ALLOCATE (ensemble_mean(nxa,nya))
ALLOCATE (rlonsa(nxa,nya),rlatsa(nxa,nya))
ALLOCATE (prob_forecast(nxa,nya))
ALLOCATE (prob_forecast_raw(nxa,nya))
ALLOCATE (prob_forecast_raw_ECMWF(nxa,nya))
ALLOCATE (prob_forecast_raw_CMC(nxa,nya))
ALLOCATE (prob_forecast_raw_NCEP(nxa,nya))
ALLOCATE (prob_forecast_cdf(nxa,nya))
ALLOCATE (prob_forecast_cdf_CMC(nxa,nya))
ALLOCATE (prob_forecast_cdf_NCEP(nxa,nya))
ALLOCATE (prob_forecast_cdf_ECMWF(nxa,nya))
ALLOCATE (climo_prob(nxa,nya))
ALLOCATE (ensemble_ccpa_x9(9,nxa,nya,nens_cmc+nens_ncep+nens_ecmwf))
ALLOCATE (ecmwf_ensemble_ccpa_x9(9,nxa,nya,nens_ecmwf))
ALLOCATE (ncep_ensemble_ccpa_x9(9,nxa,nya,nens_ncep))
ALLOCATE (cmc_ensemble_ccpa_x9(9,nxa,nya,nens_cmc))

! ---- read in the precipitation climatology appropriate to this threshold

write(6,*) 'Calling read_precip_climatology_local'
pclimo_infile = &
    '/data/thamill/Rf2_tests/ccpa_v1/0.125d/apcp_climatologies_12_to_00UTC_'//&
    cmonths(imonth)//'.nc'
CALL read_precip_climatology_local(nxa, nya, pclimo_infile, cthresh, climo_prob, &
    rlonsa, rlatsa, conusmask)
write(6,*) 'Max, Min conusmask = ',maxval(conusmask), minval(conusmask)

! ---- Read precipitation forecasts for all models/ensemble members from TDLPACK Vector files.
!      These files contains the 1/8 deg. grid.

! E, C, N, EC, EN, NC, or ENC permitted

infile_cmc_late = &
    '/Users/thamill/precip/ecmwf_data/CMC_' // cyyyymmddhh // &
    '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'
infile_cmc_early = &
    '/Users/thamill/precip/ecmwf_data/CMC_' // cyyyymmddhh // &
    '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'
center = 'CMC'
PRINT *,'reading from ', TRIM(infile_cmc_late)
CALL read_forecasts_local (nxa, nya, nens_cmc, center, &
    infile_cmc_early, infile_cmc_late, cmc_ensemble_ccpa)
PRINT *,'maxval(cmc_ensemble_ccpa) in driver = ',&
	maxval(cmc_ensemble_ccpa)

infile_ncep_late = &
    '/Users/thamill/precip/ecmwf_data/NCEP_' // cyyyymmddhh // &
    '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'   
infile_ncep_early = &
    '/Users/thamill/precip/ecmwf_data/NCEP_' // cyyyymmddhh // &
    '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'  
center = 'NCEP'
PRINT *,'reading from ', TRIM(infile_ncep_late)
CALL read_forecasts_local (nxa, nya, nens_ncep, center, &
    infile_ncep_early, infile_ncep_late, ncep_ensemble_ccpa)
PRINT *,'maxval(ncep_ensemble_ccpa) in driver = ',&
	maxval(ncep_ensemble_ccpa)

infile_ecmwf_late = &
    '/Users/thamill/precip/ecmwf_data/ECMWF_' // cyyyymmddhh // &
    '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'
infile_ecmwf_early = &
    '/Users/thamill/precip/ecmwf_data/ECMWF_' // cyyyymmddhh // &
    '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'
center = 'ECMWF'
PRINT *,'reading from ', TRIM(infile_ecmwf_late)
CALL read_forecasts_local (nxa, nya, nens_ecmwf, center, &
    infile_ecmwf_early, infile_ecmwf_late, ecmwf_ensemble_ccpa)
PRINT *,'maxval(ecmwf_ensemble_ccpa) in driver = ',&
	maxval(ecmwf_ensemble_ccpa)

! ---- read the precipitation analysis valid for this lead time.

infile = '/data/thamill/Rf2_tests/ccpa_v1/'//&
	'precip_ccpav1_2002010200_to_2016123100.nc'
CALL read_precipitation_analysis(nxa, nya, jyyyymmddhh,&
	infile, analysis, istat)

! ---- For purposes of having a baseline for comparison, generate an ensemble probability
!      simply from the relative frequency.

PRINT *, 'Calling raw_ensemble_probs_local'
CALL raw_ensemble_probs_local(nxa, nya, nens_cmc, nens_ncep, nens_ecmwf, &
    rthresh, cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    prob_forecast_raw, prob_forecast_raw_CMC, prob_forecast_raw_NCEP, &
    prob_forecast_raw_ECMWF, ensemble_mean)

! ---- read in the CDFs for each forecast and for the analyzed precip

ALLOCATE(precip_anal_cdf(nxa,nya,npct))
ALLOCATE(ecmwf_ensemble_cdf(nxa,nya,npct))
ALLOCATE(ncep_ensemble_cdf(nxa,nya,npct))
ALLOCATE(cmc_ensemble_cdf(nxa,nya,nens_cmc,npct))

write(6,*)'Calling read_cdf_netcdf_local'
CALL read_cdf_netcdf_local(nxa, nya, npct, nens_cmc, iyyyymmddhh, &
    cleade, precip_anal_cdf, ecmwf_ensemble_cdf, ncep_ensemble_cdf, &
    cmc_ensemble_cdf, thresh)
    
! ---- perform the quantile mapping and fill the _x9 arrays

CALL control_quantile_mapping_x9_local(nxa, nya, npct, nstride, nens_cmc, &
	nens_ncep, nens_ecmwf, thresh, conusmask, precip_anal_cdf, &
    ncep_ensemble_cdf, cmc_ensemble_cdf, ecmwf_ensemble_cdf, &
    cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, &
	ecmwf_ensemble_ccpa_x9, ensemble_ccpa_x9)	
	
! ---- tally up the information needed to generate dressing statistics
!      Do this first for the MME and then individually for each
!      prediction system in the MME

nine = 9
nmembers = nens_cmc + nens_ncep + nens_ecmwf
outfile = '/Users/thamill/precip/ecmwf_data/gamma_stats_'//&
	cyyyymmddhh//'_fhour'//TRIM(cleade)//'.dat'
PRINT *,'will write data to ', TRIM(outfile)
CALL tally_gamma_stats(nine, nxa, nya, nmembers, n_climocats, nthreshes, &
	thresh_low, thresh_high, ensemble_ccpa_x9, analysis, conusmask, &
	outfile, gamma_threshes, climo_prob, climo_pop_thresholds, &
	istat)
	
nmembers = nens_ncep 
outfile = '/Users/thamill/precip/ecmwf_data/gamma_stats_NCEPonly_'//&
	cyyyymmddhh//'_fhour'//TRIM(cleade)//'.dat'
PRINT *,'will write data to ', TRIM(outfile)
CALL tally_gamma_stats(nine, nxa, nya, nmembers, n_climocats, nthreshes, &
	thresh_low, thresh_high, ncep_ensemble_ccpa_x9, analysis, conusmask, &
	outfile, gamma_threshes, climo_prob, climo_pop_thresholds, &
	istat)	
	
nmembers = nens_cmc
outfile = '/Users/thamill/precip/ecmwf_data/gamma_stats_CMConly_'//&
	cyyyymmddhh//'_fhour'//TRIM(cleade)//'.dat'
PRINT *,'will write data to ', TRIM(outfile)
CALL tally_gamma_stats(nine, nxa, nya, nmembers, n_climocats, nthreshes, &
	thresh_low, thresh_high, cmc_ensemble_ccpa_x9, analysis, conusmask, &
	outfile, gamma_threshes, climo_prob, climo_pop_thresholds, &
	istat)	
		
nmembers = nens_ecmwf
outfile = '/Users/thamill/precip/ecmwf_data/gamma_stats_EConly_'//&
	cyyyymmddhh//'_fhour'//TRIM(cleade)//'.dat'
PRINT *,'will write data to ', TRIM(outfile)
CALL tally_gamma_stats(nine, nxa, nya, nmembers, n_climocats, nthreshes, &
	thresh_low, thresh_high, ecmwf_ensemble_ccpa_x9, analysis, conusmask, &
	outfile, gamma_threshes, climo_prob, climo_pop_thresholds, &
	istat)	
		
DEALLOCATE(precip_anal_cdf,ecmwf_ensemble_cdf,ncep_ensemble_cdf,&
    cmc_ensemble_cdf)

DEALLOCATE(rlonsa, rlatsa, prob_forecast, prob_forecast_raw, climo_prob,&
    cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    prob_forecast_raw_ecmwf, prob_forecast_raw_CMC, prob_forecast_raw_NCEP,&
    cmc_ensemble_ccpa_x9,ncep_ensemble_ccpa_x9,ecmwf_ensemble_ccpa_x9,&
    ensemble_mean, analysis, stat=ios)

write(6,*)'Deallocation Status = ',ios
write(6,*)'Done!'

END PROGRAM generate_gammadressing_stats
