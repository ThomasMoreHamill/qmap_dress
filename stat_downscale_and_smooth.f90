SUBROUTINE stat_downscale_and_smooth (nxa, nya, nxf, nyf, nens_cmc, nens_ncep, &
     ndays3mo, ncount, ideterministic, rnoise_percent, cmc_control_ccpa, &
     ncep_control_ccpa, cmc_ensemble_ccpa, ncep_ensemble_ccpa, &
     ihavedata, precip_anal_fine, rthresh, rlonsa, rlatsa, rlonsf, rlatsf, &
     window_size, order, topo_eighth, icm_c, inc_c, inc_e, icm_e, &
     prob_forecast, prob_forecast_raw, prob_forecast_unsmoothed)

! --- purpose:  
!     (1) add noise to each member of a certain percentage of the ensemble mean
!         to simulate the effects of sub-gridscale uncertainty
!     (2) compute a POP from ensemble relative frequency
!     (3) smooth with a Savitzky-Golay smoother

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf, nens_cmc, &
     nens_ncep, ndays3mo, ncount, icm_c, inc_c,  inc_e, &
     icm_e, window_size, order, ideterministic
REAL, INTENT(IN) :: rnoise_percent 
REAL, INTENT(IN) :: rthresh ! forecast event threshold

REAL, INTENT(IN), DIMENSION(nxa,nya) :: cmc_control_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya) :: ncep_control_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc) :: cmc_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,ndays3mo) :: precip_anal_fine
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlonsa
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlatsa
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlonsf
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlatsf
REAL, INTENT(IN), DIMENSION(nxa,nya) :: topo_eighth

INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast
REAL, INTENT(IN), DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, DIMENSION(nxa,nya) :: raw_weight
REAL, DIMENSION(nxa,nya) :: prob_forecast_smoothed
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_unsmoothed

! ---- local variables

REAL, DIMENSION(window_size, window_size) :: weights ! for Savitzky-Golay smoother
REAL, DIMENSION(ndays3mo) :: difference  ! used to tally difference between F and O

! ---- local, allocatable arrays

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: downscaled_ensemble, raw_ensemble_ccpa
REAL, ALLOCATABLE, DIMENSION(:,:) :: prob_smoothed

idum = -23434 ! seed for random number generator
    
! ---- determine the number of members in the downscaled multi-model ensemble, and then allocate
!      arrays to hold that data.  In order to give the deterministic/control forecasts from the
!      various centers more weight, we'll let them be downscaled multiple times (nmult_ecmwf and
!      nmult_other) with different analysis data used for the downscaling each time.  Also allocate
!      array that will hold a list of indices of the analysis data that is closest to this particular
!      forecast (idate_list)

!nmult_ecmwf = 5
nmult_other = 3
nens =  icm_c*nmult_other + inc_c*nmult_other + &
      icm_e*nens_cmc + inc_e*nens_ncep
ALLOCATE (downscaled_ensemble(nxa,nya,nens), raw_ensemble_ccpa(nxa,nya,nens), prob_smoothed(nxa,nya))

! ---- copy the input deterministic and ensemble data into the work array

iktr = 1
!IF (iec_d .eq. 1) THEN  ! copy in ECMWF deterministic
!   DO i = 1,nmult_ecmwf
!      raw_ensemble(:,:,iktr) = ecmwf_deterministic(:,:)
!      raw_ensemble_ccpa(:,:,iktr) = ecmwf_deterministic_ccpa(:,:)
!      iktr = iktr+1
!   END DO
!ENDIF

IF (icm_c .eq. 1) THEN  ! copy in CMC control
   DO i = 1, nmult_other
      raw_ensemble_ccpa(:,:,iktr) = cmc_control_ccpa(:,:)
      iktr = iktr+1
   END DO
ENDIF

IF (inc_c .eq. 1) THEN ! copy in NCEP control
   DO i = 1, nmult_other
      raw_ensemble_ccpa(:,:,iktr) = ncep_control_ccpa(:,:)
      iktr = iktr+1
   END DO
ENDIF

!IF (iuk_c .eq. 1) THEN ! copy in UK Met Office control
!   DO i = 1, nmult_other
!      raw_ensemble_ccpa(:,:,iktr) = ukmo_control_ccpa(:,:)
!      iktr = iktr+1
!   END DO
!ENDIF

!IF (iec_e .eq. 1) THEN ! copy in ECMWF ensemble
!   DO i = 1, nens_ecmwf
!      raw_ensemble(:,:,iktr) = ecmwf_ensemble(:,:,i)
!      raw_ensemble_ccpa(:,:,iktr) = ecmwf_ensemble_ccpa(:,:,i)
!      iktr = iktr+1
!   END DO
!ENDIF

IF (icm_e .eq. 1) THEN ! copy in CMC ensemble
   DO i = 1, nens_cmc
      raw_ensemble_ccpa(:,:,iktr) = cmc_ensemble_ccpa(:,:,i)
      iktr = iktr+1
   END DO
ENDIF

IF (inc_e .eq. 1) THEN ! copy in NCEP ensemble
   DO i = 1, nens_ncep
      raw_ensemble_ccpa(:,:,iktr) = ncep_ensemble_ccpa(:,:,i)
      iktr = iktr+1
   END DO
ENDIF

!IF (iuk_e .eq. 1) THEN ! copy in UK Met ensemble
!   DO i = 1, nens_ukmo
!      raw_ensemble_ccpa(:,:,iktr) = ukmo_ensemble_ccpa(:,:,i)
!      iktr = iktr+1
!   END DO
!ENDIF
iktr = iktr - 1


IF (iktr .ne. nens) THEN
   PRINT *,'something screwy with the number of ensemble members in statistical_downscale'
   PRINT *,'iktr = ',iktr,' nens = ',nens,' ... so stopping.'
   STOP
ENDIF

! ---- loop through all the members and then all of the training sample days.

prob_forecast(:,:) = 9999. ! initialize to a missing value flag.

downscaled_ensemble(:,:,:) = raw_ensemble_ccpa(:,:,:)

! ---- for each point, add random uniform noise to account for statistical downscaling,
!      i.e., uncertainty contributed by coarse-resolution ensemble's inability to
!      model the uncertainty at sub-grid scales

DO jya = 1, nya
   DO ixa = 1, nxa
      ensmean_point = SUM(downscaled_ensemble(ixa,jya,:)) / REAL(nens)
      DO imem = 1, nens
         rnoise_sample = (ran3(idum)-0.5)*(rnoise_percent/100.)*ensmean_point
         downscaled_ensemble(ixa,jya,imem) = downscaled_ensemble(ixa,jya,imem) + rnoise_sample
         IF (downscaled_ensemble(ixa,jya,imem) .lt. 0.0) downscaled_ensemble(ixa,jya,imem) = 0.0
      END DO ! imem
   END DO ! ixa
END DO ! jya

! ---- estimate event probability from ensemble relative frequency.

DO jya = 1, nya
   DO ixa = 1, nxa
      ncount2 = 0
      IF (ihavedata(ixa,jya) .eq. 1) THEN
         DO imem = 1, nens
            IF (downscaled_ensemble(ixa,jya,imem) .ge. rthresh) ncount2 = ncount2+1
         END DO
         prob_forecast_unsmoothed(ixa,jya) = REAL(ncount2) / REAL(nens)
      ELSE
         prob_forecast_unsmoothed(ixa,jya) =  prob_forecast_raw(ixa,jya)
      ENDIF

   END DO  ! ixa
END DO ! jya

! ---- precalculate the weights to be used in subsequent Savitzky-Golay smoother

PRINT *,'calling sgolay_2d_weights'
CALL sgolay_2d_weights(window_size, order, istat, weights)

! ---- determine how much we will weight the Savitzky-Golay smoothed fields vs.
!      the raw input.  It makes sense to weight the former more in regions
!      where there is not much topographic variation, and the latter more where
!      there is.

istat = 0
PRINT *,'calling raw_vs_smoothed_weight'
CALL raw_vs_smoothed_weight(nxa, nya, nxf, nyf, rlonsf, rlatsf, rlonsa, &
   rlatsa, topo_eighth, ihavedata, raw_weight)

! ---- now do the Savitzky-Golay smoothing, and the blending of raw and S-G smoothed
!      based on terrain variation

prob_smoothed = prob_forecast_unsmoothed
PRINT *, 'calling sgolay_smooth'
CALL sgolay_smooth(nxa, nya, ideterministic, prob_smoothed, ihavedata, weights, &
     window_size, order, istat)
prob_forecast = raw_weight*prob_forecast_unsmoothed + (1.-raw_weight)*prob_smoothed
PRINT *,' prob_smoothed min,max = ', minval(prob_smoothed), maxval(prob_smoothed)
PRINT *, 'prob_unsmoothed min,max = ', minval(prob_forecast_unsmoothed), maxval(prob_forecast_unsmoothed)
PRINT *, 'prob_forecast min,max = ', minval(prob_forecast),maxval(prob_forecast)
PRINT *, 'prob_forecast_raw min,max = ', minval(prob_forecast_raw),maxval(prob_forecast_raw)

! ---- deallocate the arrays

DEALLOCATE (raw_ensemble_ccpa, downscaled_ensemble, prob_smoothed)

PRINT *,'end stat_downscale_and_smooth'
RETURN
END SUBROUTINE stat_downscale_and_smooth
