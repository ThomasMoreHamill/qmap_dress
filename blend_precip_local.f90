PROGRAM blend_precip_local

USE netcdf

integer, parameter :: nxa = 464  ! number of grid pts in x-dir for 1/8-deg analysis grid
integer, parameter :: nya = 224  ! number of grid pts in y-dir for 1/8-deg analysis grid 
integer, parameter :: nens_ecmwf = 50 ! number of ECMWF perturbed ensemble members
integer, parameter :: nens_cmc   = 20 ! number of CMC perturbed ensemble members
integer, parameter :: nens_ncep  = 20 ! number of NCEP perturbed ensemble members
integer, parameter :: npct = 90 ! number of thresholds for CDFs

integer :: nstride
real :: stdran ! magnitude of noise for perturbing forecast quantile

REAL, DIMENSION(npct) :: thresh ! the precip amount thresholds for CDFs

CHARACTER*2 chh, cmm
CHARACTER*3, DIMENSION(12) :: cmonths 
CHARACTER*256 pclimo_infile ! name of file with climo probs
CHARACTER*256 infile_ecmwf_early ! ecmwf forecast data for today
CHARACTER*256 infile_cmc_early ! cmc forecast data for today
CHARACTER*256 infile_ncep_early ! ncep forecast data for today
CHARACTER*256 infile_ecmwf_late ! ecmwf forecast data for today
CHARACTER*256 infile_cmc_late ! cmc forecast data for today
CHARACTER*256 infile_ncep_late ! ncep forecast data for today
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

REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ecmwf_ensemble_ccpa_x9 ! ecmwf ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: cmc_ensemble_ccpa_x9 ! cmc ens precip forecast on 1/8-deg ccpa grid
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ncep_ensemble_ccpa_x9 ! ncep ens precip forecast on 1/8-deg ccpa grid

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

! --- Via command line, read in the input year/mo/day/hr and the forecast resolution 
!     we're working with.  Process date to determine the day of the month as an integer

CALL getarg(1,cyyyymmddhh)  ! input year month day hour of initial condition, 'yyyymmddhh' format
CALL getarg(2,cleade)       ! forecast lead time for beginning of precip accum period, hours, e.g.'060'
CALL getarg(3,cmodelcombo)  ! ! E, C, N, EC, EN, NC, or ENC permitted
CALL getarg(4,cthresh)      ! Precip threshold

cmm = cyyyymmddhh(5:6)
READ (cmm,'(i2)') imonth

write(6,*)' Command line arguments:'
write(6,110)  cyyyymmddhh, cleade, cmodelcombo, cthresh
110 format(1x,'cyyyymmddhh: ',A/1x,'cleade: ',A/1x,'cmodelcombo: ',A/1x,'cthresh: ',A/)

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
ALLOCATE (cmc_ensemble_ccpa_x9(9,nxa,nya,nens_cmc))
ALLOCATE (ecmwf_ensemble_ccpa_x9(9,nxa,nya,nens_ecmwf))
ALLOCATE (ncep_ensemble_ccpa_x9(9,nxa,nya,nens_ncep))

! ---- read in the precipitation climatology appropriate to this threshold

write(6,*) 'Calling read_precip_climatology_local'
pclimo_infile = &
    '/data/thamill/Rf2_tests/ccpa_v1/0.125d/apcp_climatologies_12_to_00UTC_'//&
    cmonths(imonth)//'.nc'
CALL read_precip_climatology_local(nxa, nya, pclimo_infile, cthresh, climo_prob, &
    rlonsa, rlatsa, conusmask)
write(6,*) 'Max, Min conusmask = ',maxval(conusmask), minval(conusmask)
!print *,'climo_prob(1:nxa:10,nya/2) = ', climo_prob(1:nxa:10,nya/2)
!print *,'min, max rlonsa = ', minval(rlonsa), maxval(rlonsa)
!print *,'rlatsa(1,1), rlatsa(nxa,nya) = ', rlatsa(1,1), rlatsa(nxa,nya)

! ---- Read precipitation forecasts for all models/ensemble members from TDLPACK Vector files.
!      These files contains the 1/8 deg. grid.

! E, C, N, EC, EN, NC, or ENC permitted

IF (cmodelcombo .eq. 'C' .or. cmodelcombo .eq. 'EC' .or. &
    cmodelcombo .eq. 'NC' .or. cmodelcombo .eq. 'ENC') THEN    
    infile_cmc_late = &
        '/Users/thamill/precip/ecmwf_data/CMC_' // cyyyymmddhh // &
        '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'
    infile_cmc_early = &
        '/Users/thamill/precip/ecmwf_data/CMC_' // cyyyymmddhh // &
        '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'
    CALL read_forecasts_local (nxa, nya, nens_cmc, &
        infile_cmc_early, infile_cmc_late, cmc_ensemble_ccpa)
ELSE
    cmc_ensemble_ccpa = -99.99
ENDIF

IF (cmodelcombo .eq. 'N' .or. cmodelcombo .eq. 'EN' .or. &
    cmodelcombo .eq. 'NC' .or.  cmodelcombo .eq. 'ENC') THEN   
    infile_ncep_late = &
        '/Users/thamill/precip/ecmwf_data/NCEP_' // cyyyymmddhh // &
        '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'   
    infile_ncep_early = &
        '/Users/thamill/precip/ecmwf_data/NCEP_' // cyyyymmddhh // &
        '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'   
    CALL read_forecasts_local (nxa, nya, nens_ncep, &
        infile_ncep_early, infile_ncep_late, ncep_ensemble_ccpa)
ELSE
    ncep_ensemble_ccpa = -99.99
ENDIF    

IF (cmodelcombo .eq. 'E' .or. cmodelcombo .eq. 'EN' .or. &
    cmodelcombo .eq. 'EC' .or.  cmodelcombo .eq. 'ENC') THEN  
    infile_ecmwf_late = &
        '/Users/thamill/precip/ecmwf_data/ECMWF_' // cyyyymmddhh // &
        '_leadtime' // TRIM(ADJUSTL(cleade)) // 'h.nc'
    infile_ecmwf_early = &
        '/Users/thamill/precip/ecmwf_data/ECMWF_' // cyyyymmddhh // &
        '_leadtime' // TRIM(ADJUSTL(cleadb)) // 'h.nc'
    CALL read_forecasts_local (nxa, nya, nens_ecmwf, &
        infile_ecmwf_early, infile_ecmwf_late, ecmwf_ensemble_ccpa)
ELSE
    ecmwf_ensemble_ccpa = -99.99
ENDIF  

! ---- For purposes of having a baseline for comparison, generate an ensemble probability
!      simply from the relative frequency.

PRINT *, 'Calling raw_ensemble_probs_local'
CALL raw_ensemble_probs_local(nxa, nya, nens_cmc, nens_ncep, nens_ecmwf, &
    rthresh, cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    prob_forecast_raw, prob_forecast_raw_CMC, prob_forecast_raw_NCEP, &
    prob_forecast_raw_ECMWF, ensemble_mean)
    
!PRINT *, 'prob_forecast_raw(1:nxa:10,nya/2) = ', prob_forecast_raw(1:nxa:10,nya/2)
!PRINT *, 'prob_forecast_raw_CMC(1:nxa:10,nya/2) = ', prob_forecast_raw_CMC(1:nxa:10,nya/2)
!PRINT *, 'prob_forecast_raw_NCEP(1:nxa:10,nya/2) = ', prob_forecast_raw_NCEP(1:nxa:10,nya/2)
!PRINT *, 'prob_forecast_raw_ECMWF(1:nxa:10,nya/2) = ', prob_forecast_raw_ECMWF(1:nxa:10,nya/2)

! ---- read in the CDFs for each forecast and for the analyzed precip

ALLOCATE(precip_anal_cdf(nxa,nya,npct))
ALLOCATE(ecmwf_ensemble_cdf(nxa,nya,npct))
ALLOCATE(ncep_ensemble_cdf(nxa,nya,npct))
ALLOCATE(cmc_ensemble_cdf(nxa,nya,nens_cmc,npct))

write(6,*)'Calling read_cdf_netcdf_local'
CALL read_cdf_netcdf_local(nxa, nya, npct, nens_cmc, iyyyymmddhh, &
    cleade, precip_anal_cdf, ecmwf_ensemble_cdf, ncep_ensemble_cdf, &
    cmc_ensemble_cdf, thresh)
    
! ---- compute and apply the quantile mapping bias correction, 
!      including the use of surrounding grid points and dressing

CALL control_quantile_mapping_x9_local(nxa, nya, npct, nstride, nens_cmc, & 
    nens_ncep, nens_ecmwf, stdran, thresh, conusmask, precip_anal_cdf, &
    ncep_ensemble_cdf, cmc_ensemble_cdf, ecmwf_ensemble_cdf, &
    cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, ecmwf_ensemble_ccpa_x9)
    
    
    SUBROUTINE  control_quantile_mapping_x9_local(nxa, nya, npct, nstride, nens_cmc, & 
            nens_ncep, nens_ecmwf, stdran, thresh, conusmask, precip_anal_cdf, &
            ncep_ensemble_cdf, cmc_ensemble_cdf, ecmwf_ensemble_cdf, &
            cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
            cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, ecmwf_ensemble_ccpa_x9)
    
! ---- Get probabilities from quantile-mapped (CDF bias corrected) ensemble

write(6,*)'Calling raw_ensemble_probs_x9'
CALL ensemble_probs_x9_local(nxa, nya, nens_cmc, nens_ncep, nens_ecmwf, &
    rthresh, cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, &
    ecmwf_ensemble_ccpa_x9, prob_forecast_cdf, prob_forecast_cdf_CMC, &
    prob_forecast_cdf_NCEP, prob_forecast_cdf_ECMWF)

DEALLOCATE(precip_anal_cdf,ecmwf_ensemble_cdf,ncep_ensemble_cdf,&
    cmc_ensemble_cdf)

write(6,fmt='(A)')' Model Precipitation Amount Stats (1/8 deg. grid):'
write(6,fmt='(4(A10,1X))')'MODEL','MIN','MAX','MEAN'
write(6,fmt='(A10,1X,3(F10.5,1X))')'NCEP ENS',minval(ncep_ensemble_ccpa),&
    maxval(ncep_ensemble_ccpa),sum(ncep_ensemble_ccpa)/(nxa*nya*nens_ncep)
write(6,fmt='(A10,1X,3(F10.5,1X))')'CMC ENS',minval(cmc_ensemble_ccpa),&
    maxval(cmc_ensemble_ccpa),sum(cmc_ensemble_ccpa)/(nxa*nya*nens_cmc)
write(6,fmt='(A10,1X,3(F10.5,1X))')'ECMWF ENS',minval(ecmwf_ensemble_ccpa),&
    maxval(ecmwf_ensemble_ccpa),sum(ecmwf_ensemble_ccpa)/(nxa*nya*nens_ecmwf)
    
! ---- Final pass through the final blended probability forecast to check for
!      "bad" values.
!      1) Set negative values to 0.0.
!      2) Check for NaN? Necessary?
!      
DO j=1,nya
    DO i=1,nxa
        IF (prob_forecast(i,j).lt.0.0) prob_forecast(i,j)=0.0
        IF (isnan(prob_forecast(i,j))) prob_forecast(i,j)=9999.0
    END DO
END DO

write (6,*) 'Final prob_forecast max, min = ', &
    maxval(prob_forecast_cdf), minval(prob_forecast_cdf)

IF (cmodelcombo .eq. 'E') THEN
    outfile = '/Users/thamill/precip/ecmwf_data/ECMWF_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ELSEIF (cmodelcombo .eq. 'N') THEN
    outfile = '/Users/thamill/precip/ecmwf_data/NCEP_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ELSEIF (cmodelcombo .eq. 'C') THEN
    outfile = '/Users/thamill/precip/ecmwf_data/CMC_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ELSEIF (cmodelcombo .eq. 'EC') THEN
    outfile = '/Users/thamill/precip/ecmwf_data/ECMWF_CMC_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ELSEIF (cmodelcombo .eq. 'EN') THEN
    outfile = '/Users/thamill/precip/ecmwf_data/ECMWF_NCEP_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ELSEIF (cmodelcombo .eq. 'NC') THEN
    outfile = '/Users/thamill/precip/ecmwf_data/NCEP_CMC_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ELSE
    outfile = '/Users/thamill/precip/ecmwf_data/ECMWF_NCEP_CMC_'//TRIM(cleade)//'h_IC'//cyyyymmddhh//&
        '_thresh'//TRIM(cthresh)//'.dat'
ENDIF

WRITE (6,*)'Writing data to ',TRIM(outfile)
OPEN (unit=43, file=outfile, status='unknown', form='unformatted')
WRITE(43) nxa,nya
WRITE(43) prob_forecast_raw
WRITE(43) prob_forecast_raw_CMC
WRITE(43) prob_forecast_raw_NCEP
WRITE(43) prob_forecast_raw_ECMWF
WRITE(43) prob_forecast_cdf
WRITE(43) prob_forecast_cdf_CMC
WRITE(43) prob_forecast_cdf_NCEP
WRITE(43) prob_forecast_cdf_ECMWF
WRITE(43) climo_prob
WRITE(43) rlonsa
WRITE(43) rlatsa
WRITE(43) conusmask
WRITE(43) ensemble_mean
CLOSE(43)

DEALLOCATE(rlonsa, rlatsa, prob_forecast, prob_forecast_raw, climo_prob,&
    cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    prob_forecast_raw_ecmwf, prob_forecast_raw_CMC, prob_forecast_raw_NCEP,&
    cmc_ensemble_ccpa_x9,ncep_ensemble_ccpa_x9,ecmwf_ensemble_ccpa_x9,&
    ensemble_mean,stat=ios)

write(6,*)'Deallocation Status = ',ios
write(6,*)'Done!'

END PROGRAM blend_precip_local
