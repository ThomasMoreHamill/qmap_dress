SUBROUTINE raw_ensemble_probs (nxa, nya, nens_cmc, nens_ncep, &
     rthresh, cmc_control_ccpa, ncep_control_ccpa, &
     cmc_ensemble_ccpa, ncep_ensemble_ccpa, icm_c_ccpa, inc_c_ccpa, &
     inc_e_ccpa, icm_e_ccpa, prob_forecast_raw, &
     prob_forecast_raw_CMC, prob_forecast_raw_NCEP)

! --- purpose:  conduct the statistical downscaling.  For each forecast grid point and 
!     each forecast (deterministic, control, or ensemble member), find the set of dates
!     that have past coarse-resolution analyses (at the grid scale of the forecast)
!     closest to that of the forecast.  Use that set of dates to determine a statistical
!     downscaling using the deviations between the fine-scale precip analyses and the
!     coarser-resolution precip analyses.

INTEGER, INTENT(IN) :: nxa, nya, nens_cmc, nens_ncep, &
     icm_c_ccpa, inc_c_ccpa, &
     inc_e_ccpa, icm_e_ccpa

REAL, INTENT(IN) :: rthresh  ! event threshold amount
REAL, INTENT(IN), DIMENSION(nxa,nya) :: cmc_control_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya) :: ncep_control_ccpa
!REAL, INTENT(IN), DIMENSION(nxa,nya) :: ukmo_control_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc) :: cmc_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa
!REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ukmo) :: ukmo_ensemble_ccpa

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_NCEP
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_CMC
!REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_UKMO
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: raw_ensemble

idum = -23434 ! seed for random number generator



! ---- determine the number of members in the downscaled multi-model ensemble, and then allocate
!      arrays to hold that data.  In order to give the deterministic/control forecasts from the 
!      various centers more weight, we'll let them be downscaled multiple times (nmult_ecmwf and
!      nmult_other) with different analysis data used for the downscaling each time.  Also allocate
!      array that will hold a list of indices of the analysis data that is closest to this particular
!      forecast (idate_list)

PRINT *,'icm_c_ccpa, inc_c_ccpa, icm_e_ccpa,&
     inc_e_ccpa = ', icm_c_ccpa, inc_c_ccpa, &
     icm_e_ccpa, inc_e_ccpa

nmult_other = 1
nens = icm_c_ccpa*nmult_other + inc_c_ccpa*nmult_other + &
     + icm_e_ccpa*nens_cmc + inc_e_ccpa*nens_ncep
ALLOCATE (raw_ensemble(nxa,nya,nens))

! ---- copy the input deterministic and ensemble data into the work array

iktr = 1
IF (icm_c_ccpa .eq. 1) THEN  ! copy in CMC control
   DO i = 1, nmult_other
      raw_ensemble(:,:,iktr) = cmc_control_ccpa(:,:)
      iktr = iktr+1
   END DO
ENDIF

IF (inc_c_ccpa .eq. 1) THEN ! copy in NCEP control
   DO i = 1, nmult_other
      raw_ensemble(:,:,iktr) = ncep_control_ccpa(:,:)
      iktr = iktr+1
   END DO
ENDIF

IF (icm_e_ccpa .eq. 1) THEN ! copy in CMC ensemble
   DO i = 1, nens_cmc
      raw_ensemble(:,:,iktr) = cmc_ensemble_ccpa(:,:,i)
      iktr = iktr+1
   END DO
ENDIF

IF (inc_e_ccpa .eq. 1) THEN ! copy in NCEP ensemble
   DO i = 1, nens_ncep
      raw_ensemble(:,:,iktr) = ncep_ensemble_ccpa(:,:,i)
      iktr = iktr+1
   END DO
ENDIF

!IF (iuk_e_ccpa .eq. 1) THEN ! copy in UK Met ensemble
!   DO i = 1, nens_ukmo
!      raw_ensemble(:,:,iktr) = ukmo_ensemble_ccpa(:,:,i)
!      iktr = iktr+1
!   END DO
!ENDIF
iktr = iktr - 1

IF (iktr .ne. nens) THEN
   PRINT *,'something screwy with the number of ensemble members in raw_ensemble'
   PRINT *,'iktr = ',iktr,' nens = ',nens,' ... so stopping.'
   STOP
ENDIF

PRINT *, 'rthresh = ', rthresh

! ---- loop through all the forecast grid points, but bother to process only the 
!      ones that actually have fine-grid analysis points over the CONUS associated with them.

prob_forecast_raw(:,:) = 9999 ! initialize to a missing value flag.
DO jya = 1, nya
   DO ixa = 1, nxa

      ! ---- estimate prob from ensemble relative frequency.  

      ncount2 = 0
      ncount_ncep = 0
      ncount_cmc = 0

      DO imem = 1, nens
         IF (raw_ensemble(ixa,jya,imem) .ge. rthresh) ncount2 = ncount2+1
      END DO
      prob_forecast_raw(ixa,jya) = REAL(ncount2) / REAL(nens)

      IF (inc_e_ccpa .eq. 1) THEN
         DO imem = 1, nens_ncep
            IF (NCEP_ensemble_ccpa(ixa,jya,imem) .ge. rthresh) ncount_ncep = ncount_ncep + 1
         END DO
         prob_forecast_raw_NCEP(ixa,jya) = REAL(ncount_ncep) / REAL(nens_ncep)
      ELSE
         prob_forecast_raw_NCEP(ixa,jya) = 9999
      ENDIF

      IF (icm_e_ccpa .eq. 1) THEN
         DO imem = 1, nens_cmc
            IF (CMC_ensemble_ccpa(ixa,jya,imem) .ge. rthresh) ncount_cmc = ncount_cmc + 1
         END DO
         prob_forecast_raw_CMC(ixa,jya) = REAL(ncount_cmc) / REAL(nens_cmc)
      ELSE
         prob_forecast_raw_CMC(ixa,jya) = 9999
      ENDIF

   END DO ! ixa
END DO ! jya

PRINT *,'RAW=', minval(prob_forecast_raw), maxval(prob_forecast_raw),sum(prob_forecast_raw)/(ixa*jya)
PRINT *,'NCEP=', minval(prob_forecast_raw_NCEP), maxval(prob_forecast_raw_NCEP),sum(prob_forecast_raw_NCEP)/(ixa*jya)
PRINT *,'CMC=', minval(prob_forecast_raw_CMC), maxval(prob_forecast_raw_CMC),sum(prob_forecast_raw_CMC)/(ixa*jya)

! ---- deallocate the arrays

DEALLOCATE (raw_ensemble)

RETURN
END SUBROUTINE raw_ensemble_probs
