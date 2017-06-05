SUBROUTINE raw_ensemble_probs_x9 (nxa, nya, nens_cmc, nens_ncep, &
     rthresh, cmc_control_ccpa_x9, ncep_control_ccpa_x9, &
     cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, icm_c_ccpa, inc_c_ccpa, &
     inc_e_ccpa, icm_e_ccpa, prob_forecast_raw, &
     prob_forecast_raw_CMC, prob_forecast_raw_NCEP)

! --- purpose:  generate probabilities from the quantile-mapped ensemble, here using 
!     also data from not only the current grid point but also 8 surrounding grid
!     points  (Hamill May 2016 modification)

INTEGER, INTENT(IN) :: nxa, nya, nens_cmc, nens_ncep, &
     icm_c_ccpa, inc_c_ccpa, &
     inc_e_ccpa, icm_e_ccpa

REAL, INTENT(IN) :: rthresh  ! event threshold amount
REAL, INTENT(IN), DIMENSION(9,nxa,nya) :: cmc_control_ccpa_x9  ! super-sized arrays holding also 8 surrounding points
REAL, INTENT(IN), DIMENSION(9,nxa,nya) :: ncep_control_ccpa_x9
REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_cmc) :: cmc_ensemble_ccpa_x9
REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_ncep) :: ncep_ensemble_ccpa_x9

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_NCEP
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_CMC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: raw_ensemble
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: raw_ensemble_ncep
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: raw_ensemble_cmc

INTEGER, DIMENSION(nxa,nya) :: iktr
INTEGER, DIMENSION(nxa,nya) :: iktr_cmc
INTEGER, DIMENSION(nxa,nya) :: iktr_ncep

INTEGER :: ncount2, ncount_ncep, ncount_cmc

ncount2 = 0
ncount_ncep = 0
ncount_cmc = 0

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
ALLOCATE (raw_ensemble(nxa,nya,nens*9),raw_ensemble_ncep(nxa,nya,(nens_ncep+nmult_other)*9), &
   raw_ensemble_cmc(nxa,nya,(nens_cmc+nmult_other)*9))

! ---- copy the input deterministic and ensemble data into the work array.  Now modified to
!      include the use of 8 surrounding points.

iktr(:,:) = 0
iktr_cmc(:,:) = 0
iktr_ncep(:,:) = 0

PRINT *,'raw_ensemble populating with CMC control...'
IF (icm_c_ccpa .eq. 1) THEN  ! copy in CMC control
   DO imult = 1, nmult_other
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (cmc_control_ccpa_x9(i9,ixa,jya) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_cmc(ixa,jya) = iktr_cmc(ixa,jya) + 1
                  raw_ensemble(ixa,jya,iktr(ixa,jya)) = cmc_control_ccpa_x9(i9,ixa,jya)
                  raw_ensemble_cmc(ixa,jya,iktr_cmc(ixa,jya)) = cmc_control_ccpa_x9(i9,ixa,jya)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO           ! imult
ENDIF

PRINT *,'raw_ensemble populating with NCEP control...'
IF (inc_c_ccpa .eq. 1) THEN  ! copy in NCEP control
   DO imult = 1, nmult_other
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (ncep_control_ccpa_x9(i9,ixa,jya) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_ncep(ixa,jya) = iktr_ncep(ixa,jya) + 1
                  raw_ensemble(ixa,jya,iktr(ixa,jya)) = ncep_control_ccpa_x9(i9,ixa,jya)
                  raw_ensemble_ncep(ixa,jya,iktr_ncep(ixa,jya)) = ncep_control_ccpa_x9(i9,ixa,jya)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO           ! imult
ENDIF

PRINT *,'raw_ensemble populating with CMC ensemble...'
IF (icm_e_ccpa .eq. 1) THEN ! copy in CMC ensemble
   DO imem = 1, nens_cmc
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (cmc_ensemble_ccpa_x9(i9,ixa,jya,imem) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_cmc(ixa,jya) = iktr_cmc(ixa,jya) + 1
                  raw_ensemble(ixa,jya,iktr(ixa,jya)) = cmc_ensemble_ccpa_x9(i9,ixa,jya,imem)
                  raw_ensemble_cmc(ixa,jya,iktr_cmc(ixa,jya)) = cmc_ensemble_ccpa_x9(i9,ixa,jya,imem)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO           ! imem
ENDIF

PRINT *,'raw_ensemble populating with NCEP ensemble...'
IF (inc_e_ccpa .eq. 1) THEN ! copy in NCEP ensemble
   DO imem = 1, nens_ncep
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (ncep_ensemble_ccpa_x9(i9,ixa,jya,imem) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_ncep(ixa,jya) = iktr_ncep(ixa,jya) + 1
                  raw_ensemble(ixa,jya,iktr(ixa,jya)) = ncep_ensemble_ccpa_x9(i9,ixa,jya,imem)
                  raw_ensemble_ncep(ixa,jya,iktr_ncep(ixa,jya)) = ncep_ensemble_ccpa_x9(i9,ixa,jya,imem)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO
ENDIF

PRINT *, 'rthresh = ', rthresh

! ---- loop through all the forecast grid points, but bother to process only the 
!      ones that actually have fine-grid analysis points over the CONUS associated with them.

PRINT *,'estimating probabilities'
prob_forecast_raw(:,:) = 9999.0 ! initialize to a missing value flag.
DO jya = 1, nya
   DO ixa = 1, nxa

      ! ---- estimate prob from ensemble relative frequency.  

      ncount2 = 0
      ncount_ncep = 0
      ncount_cmc = 0

      DO imem = 1, iktr(ixa,jya)
         IF (raw_ensemble(ixa,jya,imem) .ge. rthresh) ncount2 = ncount2+1
      END DO
      prob_forecast_raw(ixa,jya) = REAL(ncount2) / REAL(iktr(ixa,jya))

      IF (inc_e_ccpa .eq. 1) THEN
         DO imem = 1, iktr_ncep(ixa,jya)
            IF (raw_ensemble_ncep(ixa,jya,imem) .ge. rthresh) ncount_ncep = ncount_ncep + 1
         END DO
         prob_forecast_raw_NCEP(ixa,jya) = REAL(ncount_ncep) / REAL(iktr_ncep(ixa,jya))
      ELSE
         prob_forecast_raw_NCEP(ixa,jya) = 9999.0
      ENDIF

      IF (icm_e_ccpa .eq. 1) THEN
         DO imem = 1, iktr_cmc(ixa,jya)
            IF (raw_ensemble_cmc(ixa,jya,imem) .ge. rthresh) ncount_cmc = ncount_cmc + 1
         END DO
         prob_forecast_raw_CMC(ixa,jya) = REAL(ncount_cmc) / REAL(iktr_cmc(ixa,jya))
      ELSE
         prob_forecast_raw_CMC(ixa,jya) = 9999.0
      ENDIF

   END DO ! ixa
END DO ! jya

write(6,fmt='(A)')' Raw Model Probabilities:'
write(6,fmt='(4(A10,1X))')'MODEL','MIN','MAX','MEAN'
write(6,fmt='(A10,1X,3(F10.5,1X))')'NCEP',minval(prob_forecast_raw_NCEP),maxval(prob_forecast_raw_NCEP),sum(prob_forecast_raw_NCEP)/(nxa*nya)
write(6,fmt='(A10,1X,3(F10.5,1X))')'CMC',minval(prob_forecast_raw_CMC),maxval(prob_forecast_raw_CMC),sum(prob_forecast_raw_CMC)/(nxa*nya)
write(6,fmt='(A10,1X,3(F10.5,1X)/)')'COMBINED',minval(prob_forecast_raw),maxval(prob_forecast_raw),sum(prob_forecast_raw)/(nxa*nya)

! ---- Deallocate the arrays
DEALLOCATE (raw_ensemble, raw_ensemble_cmc, raw_ensemble_ncep)

RETURN
END SUBROUTINE raw_ensemble_probs_x9
