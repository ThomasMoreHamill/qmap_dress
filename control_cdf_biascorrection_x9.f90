SUBROUTINE control_cdf_biascorrection_x9(nxa, nya, npct, nstride, nens_cmc, &
   nens_ncep, stdran, icm_c_ccpa, inc_c_ccpa, inc_e_ccpa, icm_e_ccpa, &
   thresh, ihavedata, precip_anal_cdf, &
   ncep_control_cdf, cmc_control_cdf, ncep_ensemble_cdf, cmc_ensemble_cdf, &
   cmc_control_ccpa, ncep_control_ccpa, cmc_ensemble_ccpa, ncep_ensemble_ccpa, &
   cmc_control_ccpa_x9, ncep_control_ccpa_x9, cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9)

! purpose: control the CDF bias correct forecast precipitation amounts.  Note that the CDF-based 
!    correction is done on the 1/8-degree CCPA scale, but we also need a CDF-based correction
!    to the coarser-resolution forecasts.  This is also done here.
!
! this new _x9 version of the routine is going to do two new things:
!   (1) the quantile mapping will happen not only from forecast (i,j) but from 8 other
!       surrounding grid points.  In this way we attempt to account for position errors.  There
!       is an input "nstride" that controls how afar we go looking for other grid points; 
!       presumably there are larger position errors at longer leads and we want to permit
!       a larger nstride for those.
!   (2) the routine now permits a stochastic perturbation of the forecast quantile 
!       via the "stdran" passed in

INTEGER, INTENT(IN) :: nxa, nya, npct, nstride, nens_cmc, &
     nens_ncep, icm_c_ccpa, inc_c_ccpa, &
     inc_e_ccpa, icm_e_ccpa
REAL, INTENT(IN) :: stdran
REAL, INTENT(IN), DIMENSION(npct) :: thresh ! threshold amts for CDF
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct)::  precip_anal_cdf, &
     ncep_control_cdf, cmc_control_cdf, ncep_ensemble_cdf
REAL*8, INTENT(IN), DIMENSION(nxa,nya,nens_cmc,npct)::  cmc_ensemble_cdf
REAL*8, DIMENSION(nxa,nya,npct) :: work_cdf

REAL, INTENT(IN), DIMENSION(nxa,nya) :: cmc_control_ccpa, ncep_control_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc)  :: cmc_ensemble_ccpa 
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa

REAL, INTENT(OUT), DIMENSION(9,nxa,nya) :: cmc_control_ccpa_x9, ncep_control_ccpa_x9
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_cmc)  :: cmc_ensemble_ccpa_x9
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_ncep) :: ncep_ensemble_ccpa_x9

REAL, DIMENSION(9,nxa,nya) :: forecast_x9

! ----- determine the CDF-based bias correction for the NCEP control 
!       forecast

IF (inc_c_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to NCEP control'
   CALL cdf_correct_x9(nxa, nya, npct, nstride, stdran, thresh, ihavedata, ncep_control_cdf, &
        precip_anal_cdf, ncep_control_ccpa, ncep_control_ccpa_x9)
ENDIF

! ----- determine the CDF-based bias correction for the ECMWF deterministic 
!       forecast

IF (icm_c_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to CMC control'
   CALL cdf_correct_x9(nxa, nya, npct, nstride, stdran, thresh, ihavedata, cmc_control_cdf, &
        precip_anal_cdf, cmc_control_ccpa, cmc_control_ccpa_x9)
ENDIF

! --- and now the bias corrections for the ensemble forecasts

IF (icm_e_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to CMC ensemble, members = ',nens_cmc
   DO imem = 1, nens_cmc
      PRINT *,'calling cdf_correct_x9 for member = ',imem
      work_cdf(:,:,:) = cmc_ensemble_cdf(:,:,imem,:)
      CALL cdf_correct_x9(nxa, nya, npct, nstride, stdran, thresh, ihavedata, work_cdf, &
           precip_anal_cdf, cmc_ensemble_ccpa(1,1,imem), forecast_x9)
      PRINT *,'imem, cmc_ensemble_ccpa(nxa/2,nya/2,imem) before = ',&
           imem, cmc_ensemble_ccpa(nxa/2,nya/2,imem)
      cmc_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
      PRINT *,'imem, cmc_ensemble_ccpa_x9(5,nxa/2,nya/2,imem) after = ',&
           imem, cmc_ensemble_ccpa_x9(5,nxa/2,nya/2,imem)
   END DO
ENDIF

IF (inc_e_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to NCEP ensemble, members = ',nens_ncep
   DO imem = 1, nens_ncep
      PRINT *,'calling cdf_correct_x9 for member = ',imem
      CALL cdf_correct_x9(nxa, nya, npct, nstride, stdran, thresh, ihavedata, ncep_ensemble_cdf, &
           precip_anal_cdf, ncep_ensemble_ccpa(1,1,imem), forecast_x9)
      ncep_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
   END DO
ENDIF
print *,'done in control_cdf_biascorrection_x9'

RETURN
END SUBROUTINE control_cdf_biascorrection_x9
