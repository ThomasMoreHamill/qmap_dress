SUBROUTINE ensemble_probs_x9_local(nxa, nya, nens_cmc, nens_ncep, nens_ecmwf, &
    rthresh, cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, &
    ecmwf_ensemble_ccpa_x9, prob_forecast, prob_forecast_CMC, &
    prob_forecast_NCEP, prob_forecast_ECMWF)

! --- purpose:  generate probabilities from the quantile-mapped ensemble, here using 
!     also data from not only the current grid point but also 8 surrounding grid
!     points  (Hamill May 2016 modification)

INTEGER, INTENT(IN) :: nxa, nya, nens_cmc, nens_ncep, nens_ecmwf 
REAL, INTENT(IN) :: rthresh  ! event threshold amount

REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_cmc) :: cmc_ensemble_ccpa_x9
REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_ncep) :: ncep_ensemble_ccpa_x9
REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_ecmwf) :: ecmwf_ensemble_ccpa_x9

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_CMC
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_NCEP
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_ECMWF
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble_ncep
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble_cmc
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble_ecmwf

INTEGER, DIMENSION(nxa,nya) :: iktr
INTEGER, DIMENSION(nxa,nya) :: iktr_cmc
INTEGER, DIMENSION(nxa,nya) :: iktr_ncep
INTEGER, DIMENSION(nxa,nya) :: iktr_ecmwf

INTEGER :: ncount2, ncount_ncep, ncount_cmc, ncount_ecmwf

ncount2 = 0
ncount_ncep = 0
ncount_cmc = 0
ncount_ecmwf = 0

idum = -23434 ! seed for random number generator

! ---- determine the number of members in the downscaled multi-model ensemble, and then allocate
!      arrays to hold that data.  In order to give the deterministic/control forecasts from the 
!      various centers more weight, we'll let them be downscaled multiple times (nmult_ecmwf and
!      nmult_other) with different analysis data used for the downscaling each time.  Also allocate
!      array that will hold a list of indices of the analysis data that is closest to this particular
!      forecast (idate_list)

nens = nens_cmc + nens_ncep + nens_ecmwf
!print *,'nens, nens_cmc, nens_ncep, nens_ecmwf = ', nens, nens_cmc, nens_ncep, nens_ecmwf

ALLOCATE(ensemble(nxa,nya,nens*9), ensemble_ncep(nxa,nya,nens_ncep*9), &
    ensemble_cmc(nxa,nya,nens_cmc*9), ensemble_ecmwf(nxa,nya,nens_ecmwf*9) )

iktr(:,:) = 0
iktr_cmc(:,:) = 0
iktr_ncep(:,:) = 0
iktr_ecmwf(:,:) = 0

PRINT *,'raw_ensemble populating with CMC ensemble...'
DO imem = 1, nens_cmc
    IF (cmc_ensemble_ccpa_x9(5,nxa/2,nya/2,imem) .ge. -98.) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                DO i9 = 1,9
                    iktr(ixa,jya) = iktr(ixa,jya) + 1
                    iktr_cmc(ixa,jya) = iktr_cmc(ixa,jya) + 1
                    ensemble(ixa,jya,iktr(ixa,jya)) = cmc_ensemble_ccpa_x9(i9,ixa,jya,imem)
                    ensemble_cmc(ixa,jya,iktr_cmc(ixa,jya)) = cmc_ensemble_ccpa_x9(i9,ixa,jya,imem)
                END DO  ! i9
            END DO     ! ixa
        END DO        ! jya
    ENDIF
END DO  ! imem

PRINT *,'raw_ensemble populating with NCEP ensemble...'
DO imem = 1, nens_ncep
    IF (ncep_ensemble_ccpa_x9(5,nxa/2,nya/2,imem) .ge. -98.) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                DO i9 = 1,9
                    iktr(ixa,jya) = iktr(ixa,jya) + 1
                    iktr_ncep(ixa,jya) = iktr_ncep(ixa,jya) + 1
                    ensemble(ixa,jya,iktr(ixa,jya)) = ncep_ensemble_ccpa_x9(i9,ixa,jya,imem)
                    ensemble_ncep(ixa,jya,iktr_ncep(ixa,jya)) = ncep_ensemble_ccpa_x9(i9,ixa,jya,imem)
                END DO  ! i9
            END DO     ! ixa
        END DO        ! jya
    ENDIF
END DO  ! imem

PRINT *,'raw_ensemble populating with ECMWF ensemble...'
DO imem = 1, nens_ecmwf
    IF (ecmwf_ensemble_ccpa_x9(5,nxa/2,nya/2,imem) .ge. -98.) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                DO i9 = 1,9
                    iktr(ixa,jya) = iktr(ixa,jya) + 1
                    iktr_ecmwf(ixa,jya) = iktr_ecmwf(ixa,jya) + 1
                    ensemble(ixa,jya,iktr(ixa,jya)) = ecmwf_ensemble_ccpa_x9(i9,ixa,jya,imem)
                    ensemble_ecmwf(ixa,jya,iktr_ecmwf(ixa,jya)) = ecmwf_ensemble_ccpa_x9(i9,ixa,jya,imem)
                END DO  ! i9
            END DO     ! ixa
        END DO        ! jya
    ENDIF
END DO ! imem

!PRINT *, 'rthresh = ', rthresh

! ---- loop through all the forecast grid points, but bother to process only the 
!      ones that actually have fine-grid analysis points over the CONUS associated with them.

PRINT *,'estimating probabilities'
prob_forecast(:,:) = -99.99 ! initialize to a missing value flag.
DO jya = 1, nya
    DO ixa = 1, nxa

        ! ---- estimate prob from ensemble relative frequency.  

        ncount2 = 0
        ncount_ncep = 0
        ncount_cmc = 0
        ncount_ecmwf = 0

        IF (iktr(ixa,jya) .gt. 0) THEN
            DO imem = 1, iktr(ixa,jya)
                IF (ensemble(ixa,jya,imem) .ge. rthresh) ncount2 = ncount2+1
            END DO
            prob_forecast(ixa,jya) = REAL(ncount2) / REAL(iktr(ixa,jya))
        ELSE
            prob_forecast(ixa,jya) = -99.99
        ENDIF

        IF (iktr_ncep(ixa,jya) .gt. 0) THEN
            DO imem = 1, iktr_ncep(ixa,jya)
                IF (ensemble_ncep(ixa,jya,imem) .ge. rthresh) ncount_ncep = ncount_ncep + 1
            END DO
            prob_forecast_NCEP(ixa,jya) = REAL(ncount_ncep) / REAL(iktr_ncep(ixa,jya))
        ELSE
            prob_forecast_NCEP(ixa,jya) = -99.99
        ENDIF
        
        IF (iktr_ecmwf(ixa,jya) .gt. 0) THEN
            DO imem = 1, iktr_ecmwf(ixa,jya)
                IF (ensemble_ecmwf(ixa,jya,imem) .ge. rthresh) ncount_ecmwf = ncount_ecmwf + 1
            END DO
            prob_forecast_ECMWF(ixa,jya) = REAL(ncount_ecmwf) / REAL(iktr_ecmwf(ixa,jya))
        ELSE
            prob_forecast_ECMWF(ixa,jya) = -99.99
        ENDIF
        
        IF (iktr_cmc(ixa,jya) .gt. 0) THEN
            DO imem = 1, iktr_cmc(ixa,jya)
                IF (ensemble_cmc(ixa,jya,imem) .ge. rthresh) ncount_cmc = ncount_cmc + 1
            END DO
            prob_forecast_CMC(ixa,jya) = REAL(ncount_cmc) / REAL(iktr_cmc(ixa,jya))
        ELSE
            prob_forecast_CMC(ixa,jya) = -99.99
        ENDIF

   END DO ! ixa
END DO ! jya

write (6,fmt='(A)')'Postprocessed Model Probabilities:'
write (6,fmt='(4(A10,1X))')'MODEL','MIN','MAX','MEAN'
write (6,fmt='(A10,1X,3(F10.5,1X))')'NCEP', minval(prob_forecast_NCEP), &
    maxval(prob_forecast_NCEP), sum(prob_forecast_NCEP)/(nxa*nya)
write (6,fmt='(A10,1X,3(F10.5,1X))')'CMC',minval(prob_forecast_CMC), &
    maxval(prob_forecast_CMC), sum(prob_forecast_CMC)/(nxa*nya)
write (6,fmt='(A10,1X,3(F10.5,1X))')'ECMWF',minval(prob_forecast_ECMWF), &
    maxval(prob_forecast_ECMWF), sum(prob_forecast_ECMWF)/(nxa*nya)
write (6,fmt='(A10,1X,3(F10.5,1X)/)')'COMBINED',minval(prob_forecast), &
    maxval(prob_forecast), sum(prob_forecast)/(nxa*nya)

! ---- Deallocate the arrays

DEALLOCATE(ensemble, ensemble_cmc, ensemble_ncep, ensemble_ecmwf)

RETURN
END SUBROUTINE ensemble_probs_x9_local
