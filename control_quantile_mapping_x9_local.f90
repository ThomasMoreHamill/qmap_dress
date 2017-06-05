    
SUBROUTINE  control_quantile_mapping_x9_local (&
	nxa, nya, npct, nstride, nens_cmc, & 
	nens_ncep, nens_ecmwf, thresh, conusmask, precip_anal_cdf, &
    ncep_ensemble_cdf, cmc_ensemble_cdf, ecmwf_ensemble_cdf, &
    cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, &
	ecmwf_ensemble_ccpa_x9, ensemble_ccpa_x9)

INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions 
INTEGER, INTENT(IN) :: npct ! number of precip thresholds where CDF tallied
INTEGER, INTENT(IN) :: nstride ! stride length when skipping grid pts
INTEGER, INTENT(IN) :: nens_cmc, nens_ncep, nens_ecmwf ! ensemble sizes


REAL, INTENT(IN), DIMENSION(npct) :: thresh ! threshold amts for CDF
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask ! conus mask
REAL, INTENT(IN), DIMENSION(nxa,nya,npct):: precip_anal_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya,npct):: ecmwf_ensemble_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya,npct):: ncep_ensemble_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc,npct) :: cmc_ensemble_cdf
! note that Canadian system has different forecast biases assoc'd with diff mbrs

REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc) :: cmc_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ecmwf) :: ecmwf_ensemble_ccpa

REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_cmc)  :: cmc_ensemble_ccpa_x9
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_ncep) :: ncep_ensemble_ccpa_x9
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_ecmwf) :: ecmwf_ensemble_ccpa_x9
REAL, INTENT(OUT), DIMENSION(9,nxa,nya,nens_ecmwf+nens_ncep+nens_cmc) :: &
	ensemble_ccpa_x9

REAL, DIMENSION(9,nxa,nya) :: forecast_x9 ! work array
REAL, DIMENSION(nxa,nya,npct):: work_cdf

PRINT *,'inside control_quantile_mapping_x9_local'
!print *,'cmc_ensemble_ccpa(1:nxa:10,nya/2,1) = ', cmc_ensemble_ccpa(1:nxa:10,nya/2,1)

! --- and now the quantile mappings for the ensemble forecasts

PRINT *,'Quantile mapping and dressing of CMC ensemble, # members = ',nens_cmc
imemtot = 1
DO imem = 1, nens_cmc
	
	!PRINT *,'********** processing CMC member = ', imem
    work_cdf(:,:,:) = cmc_ensemble_cdf(:,:,imem,:)
    !print *,'work_cdf(3*nxa/4,nya/2,:) = ', work_cdf(3*nxa/4,nya/2,:) 
    !print *,'thresh = ', thresh(:)
    CALL cdf_correct_x9_local(nxa, nya, npct, nstride, &
        thresh, conusmask, work_cdf, &
        precip_anal_cdf, cmc_ensemble_ccpa(1,1,imem), &
        forecast_x9)
    cmc_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
    ensemble_ccpa_x9(:,:,:,imemtot) = forecast_x9(:,:,:)
	imemtot = imemtot + 1
END DO

!print *,'thresh = ', thresh
!stop

!print *,'after: cmc_ensemble_ccpa_x9(1,1:nxa:10,nya/2,1) = ', cmc_ensemble_ccpa_x9(1,1:nxa:10,nya/2,1)
!print *,'after: cmc_ensemble_ccpa_x9(5,1:nxa:10,nya/2,1) = ', cmc_ensemble_ccpa_x9(5,1:nxa:10,nya/2,1)
!print *,'after: cmc_ensemble_ccpa_x9(9,1:nxa:10,nya/2,1) = ', cmc_ensemble_ccpa_x9(9,1:nxa:10,nya/2,1)

PRINT *,'Quantile mapping and dressing of NCEP ensemble, # members = ',nens_ncep
DO imem = 1, nens_ncep
    CALL cdf_correct_x9_local(nxa, nya, npct, nstride, &
        thresh, conusmask, ncep_ensemble_cdf, &
        precip_anal_cdf, ncep_ensemble_ccpa(1,1,imem), &
        forecast_x9)      
    ncep_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
    ensemble_ccpa_x9(:,:,:,imemtot) = forecast_x9(:,:,:)
	imemtot = imemtot + 1
END DO

PRINT *,'Quantile mapping and dressing of ECMWF ensemble, # members = ',nens_ecmwf
DO imem = 1, nens_ecmwf
    CALL cdf_correct_x9_local(nxa, nya, npct, nstride, &
        thresh, conusmask, ecmwf_ensemble_cdf, &
        precip_anal_cdf, ecmwf_ensemble_ccpa(1,1,imem), &
        forecast_x9)      
    ecmwf_ensemble_ccpa_x9(:,:,:,imem) = forecast_x9(:,:,:)
    ensemble_ccpa_x9(:,:,:,imemtot) = forecast_x9(:,:,:)
	imemtot = imemtot + 1
END DO

!print *,'cmc_ensemble_ccpa_x9(:,3*nxa/4, nya/2, :) = ', cmc_ensemble_ccpa_x9(:,3*nxa/4, nya/2, :)
!print *,'ncep_ensemble_ccpa_x9(:,3*nxa/4, nya/2, :) = ', ncep_ensemble_ccpa_x9(:,3*nxa/4, nya/2, :)
!print *,'ecmwf_ensemble_ccpa_x9(:,3*nxa/4, nya/2, :) = ', ecmwf_ensemble_ccpa_x9(:,3*nxa/4, nya/2, :)

RETURN
END SUBROUTINE control_quantile_mapping_x9_local
