SUBROUTINE control_cdf_biascorrection(nxa, nya, nxf, nyf, npct, nens_cmc, &
   nens_ncep, nnearest, icm_c_ccpa, inc_c_ccpa, inc_e_ccpa, icm_e_ccpa, &
   ncount, ilist, jlist, thresh, ihavedata, &
   precip_anal_cdf, ncep_control_cdf, cmc_control_cdf, &
   ncep_ensemble_cdf, cmc_ensemble_cdf, &
   cmc_control, ncep_control, cmc_ensemble, ncep_ensemble, &
   cmc_control_ccpa, ncep_control_ccpa, &
   cmc_ensemble_ccpa, ncep_ensemble_ccpa)

! purpose: control the CDF bias correct forecast precipitation amounts.  Note that the CDF-based 
!    correction is done on the 1/8-degree CCPA scale, but we also need a CDF-based correction
!    to the coarser-resolution forecasts.  This is also done here.

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf, npct, nens_cmc, &
     nens_ncep, icm_c_ccpa, inc_c_ccpa, &
     inc_e_ccpa, icm_e_ccpa
INTEGER, INTENT(IN), DIMENSION(nxf,nyf) :: nnearest
INTEGER, INTENT(IN), DIMENSION(nxf,nyf,ncount) :: ilist, jlist
REAL, INTENT(IN), DIMENSION(npct) :: thresh ! threshold amts for CDF
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct)::  precip_anal_cdf, &
     ncep_control_cdf, cmc_control_cdf, &
     ncep_ensemble_cdf, cmc_ensemble_cdf
REAL, INTENT(INOUT), DIMENSION(nxf,nyf) :: cmc_control, ncep_control
REAL, INTENT(INOUT), DIMENSION(nxf,nyf,nens_cmc)  :: cmc_ensemble 
REAL, INTENT(INOUT), DIMENSION(nxf,nyf,nens_ncep) :: ncep_ensemble

REAL, INTENT(INOUT), DIMENSION(nxa,nya) :: cmc_control_ccpa, ncep_control_ccpa
REAL, INTENT(INOUT), DIMENSION(nxa,nya,nens_cmc)  :: cmc_ensemble_ccpa 
REAL, INTENT(INOUT), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa
!REAL, INTENT(INOUT), DIMENSION(nxa,nya,nens_ukmo) :: ukmo_ensemble_ccpa

REAL*8, DIMENSION(nxa,nya,npct) :: cdfwork
REAL work(nxa,nya)

! ----- determine the CDF-based bias correction for the NCEP control 
!       forecast

IF (inc_c_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to NCEP control'
   work = ncep_control_ccpa
   CALL cdf_correct(nxa, nya, npct, thresh, ihavedata, ncep_control_cdf, &
        precip_anal_cdf, ncep_control_ccpa)
   CALL upscale_adjusted(nxa, nya, nxf, nyf, nnearest, ncount, ilist, jlist, &
        ncep_control_ccpa, ncep_control)
ENDIF

! ----- determine the CDF-based bias correction for the ECMWF deterministic 
!       forecast

IF (icm_c_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to CMC control'
   CALL cdf_correct(nxa, nya, npct, thresh, ihavedata, cmc_control_cdf, &
        precip_anal_cdf, cmc_control_ccpa)
   CALL upscale_adjusted(nxa, nya, nxf, nyf, nnearest, ncount, ilist, jlist, &
        cmc_control_ccpa, cmc_control)
ENDIF

! --- and now the bias corrections for the ensemble forecasts

IF (icm_e_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to CMC ensemble'
   DO imem = 1, nens_cmc
      CALL cdf_correct(nxa, nya, npct, thresh, ihavedata, cmc_ensemble_cdf, &
           precip_anal_cdf, cmc_ensemble_ccpa(1,1,imem))
      CALL upscale_adjusted(nxa, nya, nxf, nyf, nnearest, ncount, ilist, jlist, &
           cmc_ensemble_ccpa(1,1,imem), cmc_ensemble(1,1,imem))
   END DO
ENDIF

IF (inc_e_ccpa .eq. 1) THEN
   PRINT *,'CDF adjustment to NCEP ensemble'
   DO imem = 1, nens_ncep
      cdfwork(:,:,:) = ncep_ensemble_cdf(:,:,:)
      CALL cdf_correct(nxa, nya, npct, thresh, ihavedata, ncep_ensemble_cdf, &
           precip_anal_cdf, ncep_ensemble_ccpa(1,1,imem))
      CALL upscale_adjusted(nxa, nya, nxf, nyf, nnearest, ncount, ilist, jlist, &
           ncep_ensemble_ccpa(1,1,imem), ncep_ensemble(1,1,imem))
   END DO
ENDIF


RETURN
END SUBROUTINE control_cdf_biascorrection
