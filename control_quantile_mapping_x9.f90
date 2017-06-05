SUBROUTINE control_quantile_mapping_x9(nxa, nya, npct, nstride, nens_cmc, &
    nens_ncep, icm_c_ccpa, inc_c_ccpa, inc_e_ccpa, icm_e_ccpa, &
    thresh, ihavedata, precip_anal_cdf, ncep_control_cdf, &
    cmc_control_cdf, ncep_ensemble_cdf, cmc_ensemble_cdf, &
    cmc_control_ccpa, ncep_control_ccpa, cmc_ensemble_ccpa, &
    ncep_ensemble_ccpa, cmc_control_ccpa_x9, ncep_control_ccpa_x9, &
    cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9)

! purpose: control the quantile mapping of forecast precipitation amounts.  
!
! this version of the routine applies the quantile mapping 
! not only from forecast (i,j) but from 8 other surrounding grid points.  
! In this way we attempt to account for position errors.  There
! is an input "nstride" that controls how afar we go looking for other grid points;
! possibly there are larger position errors at longer leads and we want to permit
! a larger nstride for those.

INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions 
INTEGER, INTENT(IN) :: npct ! number of precip thresholds where CDF tallied
INTEGER, INTENT(IN) :: nstride ! stride length when skipping grid pts
INTEGER, INTENT(IN) :: nens_cmc, nens_ncep ! ensemble sizes
INTEGER, INTENT(IN) :: icm_c_ccpa, inc_c_ccpa ! control data availability flags
INTEGER, INTENT(IN) :: inc_e_ccpa, icm_e_ccpa ! ensemble data availability flags

REAL, INTENT(IN), DIMENSION(npct) :: thresh ! threshold amts for CDF
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata ! conus mask
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct):: precip_anal_cdf
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct):: ncep_control_cdf
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct):: cmc_control_cdf
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct):: ncep_ensemble_cdf
REAL*8, INTENT(IN), DIMENSION(nxa,nya,nens_cmc,npct) :: cmc_ensemble_cdf
! note that Canadian system has different forecast biases assoc'd with diff mbrs

REAL, INTENT(IN), DIMENSION(nxa,nya) :: cmc_control_ccpa !fcst interp'd to CCPA grid
REAL, INTENT(IN), DIMENSION(nxa,nya) :: ncep_control_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc) :: cmc_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa

REAL, INTENT(OUT), DIMENSION(9,nxa,nya) :: cmc_control_ccpa_x9 ! quantile mapped
! forecasts for CMC control including 8 surrounding grid points
REAL, INTENT(OUT), DIMENSION(9,nxa,nya) :: ncep_control_ccpa_x9 ! same but NCEP
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_cmc)  :: cmc_ensemble_ccpa_x9
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_ncep) :: ncep_ensemble_ccpa_x9

REAL, DIMENSION(9,nxa,nya) :: forecast_x9 ! work array
REAL*8, DIMENSION(nxa,nya,npct):: work_cdf

! ----- determine the quantile mapping for the NCEP control forecast

IF (inc_c_ccpa .eq. 1) THEN
   PRINT *,'Quantile mapping of NCEP control'
   CALL cdf_correct_nonoise_x9(nxa, nya, npct, nstride, &
       thresh, ihavedata, ncep_control_cdf, precip_anal_cdf, &
       ncep_control_ccpa, ncep_control_ccpa_x9)
ENDIF

! ----- determine the quantile mapping for the CMC control forecast

IF (icm_c_ccpa .eq. 1) THEN
   PRINT *,'Quantile mapping of CMC control'
   CALL cdf_correct_nonoise_x9(nxa, nya, npct, nstride, &
       thresh, ihavedata, cmc_control_cdf, precip_anal_cdf, &
       cmc_control_ccpa, cmc_control_ccpa_x9)
ENDIF

! --- and now the quantile mappings for the ensemble forecasts

IF (icm_e_ccpa .eq. 1) THEN
   PRINT *,'Quantile mapping of CMC ensemble, # members = ',nens_cmc
   DO imem = 1, nens_cmc
      work_cdf(:,:,:) = cmc_ensemble_cdf(:,:,imem,:)
      CALL cdf_correct_nonoise_x9(nxa, nya, npct, nstride, &
          thresh, ihavedata, work_cdf, &
          precip_anal_cdf, cmc_ensemble_ccpa(1,1,imem), &
          forecast_x9)
      cmc_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
   END DO
ENDIF

IF (inc_e_ccpa .eq. 1) THEN
   PRINT *,'Quantile mapping of NCEP ensemble, # members = ',nens_ncep
   DO imem = 1, nens_ncep
      CALL cdf_correct_nonoise_x9(nxa, nya, npct, nstride, &
          thresh, ihavedata, ncep_ensemble_cdf, &
          precip_anal_cdf, ncep_ensemble_ccpa(1,1,imem), &
          forecast_x9)      
      ncep_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
   END DO
ENDIF

!PRINT *, 'ncep_ensemble_ccpa_x9(:,479,199,:) = ',ncep_ensemble_ccpa_x9(:,479,199,:)
!PRINT *, 'ncep_ensemble_ccpa(479,199,:) = ',ncep_ensemble_ccpa(479,199,:)
!PRINT *, 'nstride = ',nstride
!stop

RETURN
END SUBROUTINE control_quantile_mapping_x9
