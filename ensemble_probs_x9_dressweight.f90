SUBROUTINE ensemble_probs_x9_dressweight (nxa, nya, nens_cmc, nens_ncep, &
    nmult_other, nmembers, n_amounts, n_climocats, cdressing, rthresh, &
    ramt, fraczero, fraczero_fclimpop, climo_pop_thresholds, &
    gamma_shape, gamma_scale, closest_histogram, climo_pop, &
    ihavedata, cmc_control_ccpa_x9, ncep_control_ccpa_x9, &
    cmc_ensemble_ccpa_x9, ncep_ensemble_ccpa_x9, icm_c_ccpa, inc_c_ccpa, &
    inc_e_ccpa, icm_e_ccpa, prob_forecast_qmap, prob_forecast_qmap_CMC, &
    prob_forecast_qmap_NCEP, prob_forecast_unsmoothed)
     
! --- purpose:  generate probabilities from the NCEP, CMC, and 
!     combined quantile-mapped ensemble, here using
!     also data from not only the current grid point but also 8 surrounding grid
!     points  (Hamill May 2016 modification).

!     Modification Jan 2017: in generating the probabilities, use the 
!     weights specified in the input closest_histogram vector and calculate
!     probabilities from weighted sum of exceedance probabilities of
!     the ensemble dressing distributions associated with every quantile-mapped
!     ensemble member.  Another modification is that the dressed probabilities 
!     are weighted by "closest_histogram" information which reflects how near
!     the truth is to a given sorted ensemble member (based on past stats).
! 
!     Gamma CDF probabilities taken from fortran code cumgam.f borrowed from 
!     netlib.org/random (http://www.netlib.org/random/dcdflib.f.tar.gz)

INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions
INTEGER, INTENT(IN) :: nens_cmc, nens_ncep ! # of CMC, NCEP members
INTEGER, INTENT(IN) :: nmembers ! expected number of members when all systems available
INTEGER, INTENT(IN) :: n_amounts ! dimension of gamma distribution arrays below
INTEGER, INTENT(IN) :: n_climocats
INTEGER, INTENT(IN) :: icm_c_ccpa, inc_c_ccpa, inc_e_ccpa, icm_e_ccpa ! data avail flag
CHARACTER*(*), INTENT(IN) :: cdressing ! 'y' or 'Y' if we wish to best-member dress with gamma
!  distributed noise

INTEGER*2, INTENT(IN), DIMENSION(nxa, nya) :: ihavedata ! conus mask

REAL*8, external :: random_gamma   ! generates random gamma deviate

REAL, INTENT(IN) :: rthresh  ! event threshold amount
REAL, INTENT(IN), DIMENSION(9,nxa,nya) :: cmc_control_ccpa_x9  
    ! super-sized arrays holding also 8 surrounding points
REAL, INTENT(IN), DIMENSION(9,nxa,nya) :: ncep_control_ccpa_x9
REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_cmc) :: cmc_ensemble_ccpa_x9
REAL, INTENT(IN), DIMENSION(9,nxa,nya,nens_ncep) :: ncep_ensemble_ccpa_x9
REAL, INTENT(IN), DIMENSION(nmembers) :: closest_histogram ! likelihood this sorted mbr closest
REAL, INTENT(IN), DIMENSION(nxa,nya) :: climo_pop ! climatological probability of precip
REAL, INTENT(IN), DIMENSION(n_amounts) :: ramt ! list of precip amts for quantifying best
    ! member dressing parameters.
REAL, INTENT(IN), DIMENSION(n_amounts) :: fraczero ! fraction of best-mbr with zero precip
REAL, INTENT(IN), DIMENSION(n_climocats) :: fraczero_fclimpop ! fraction zero as f(clim POP)
REAL, INTENT(IN), DIMENSION(n_climocats-1) :: climo_pop_thresholds ! climatology POP dividers

REAL, INTENT(IN), DIMENSION(n_amounts) :: gamma_shape ! alpha parameter for Gamma dist
REAL, INTENT(IN), DIMENSION(n_amounts) :: gamma_scale ! beta parameter for Gamma dist

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_qmap ! final output probability
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_qmap_NCEP ! output from NCEP w/o stoch noise
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_qmap_CMC ! output from NCEP w/o stoch noise
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_unsmoothed ! quantile mapped and dressed

REAL*8 alpha ! gamma shape parameter
REAL*8 beta  ! gamma scale parameter
REAL*8 e8    ! ensemble value in double precision
REAL*8 ksi   ! threshold / beta
REAL*8 cum   ! cumulative incomplete gamma distribution
REAL*8 ccum  ! complement of the cumulative incomplete gamma distribution
REAL*8 pz
REAL*8 pzero
REAL*8 probsum_belowthresh
REAL*8 pgamma_lt_thresh
REAL*8 probsum_pzero
REAL*8 probsum_gamma

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: qmap_ensemble
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: qmap_ensemble_ncep
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: qmap_ensemble_cmc
REAL, ALLOCATABLE, DIMENSION(:) :: member_weights
REAL, ALLOCATABLE, DIMENSION(:) :: ensemble_today

INTEGER, DIMENSION(nxa,nya) :: iktr
INTEGER, DIMENSION(nxa,nya) :: iktr_cmc
INTEGER, DIMENSION(nxa,nya) :: iktr_ncep

LOGICAL first

first = .True.
idum = -1*INT(SUM(prob_forecast_qmap)/nxa)

! ---- determine the number of members in the downscaled multi-model ensemble, 
!      and then allocate arrays to hold that data.  In order to give the 
!      deterministic/control forecasts from the various centers more weight, 
!      we'll let them be downscaled multiple times (nmult_ecmwf and nmult_other) 
!      with different analysis data used for the downscaling each time. 

PRINT *,'icm_c_ccpa, inc_c_ccpa, icm_e_ccpa, inc_e_ccpa = ', &
    icm_c_ccpa, inc_c_ccpa, icm_e_ccpa, inc_e_ccpa

nens = icm_c_ccpa*nmult_other + inc_c_ccpa*nmult_other &
     + icm_e_ccpa*nens_cmc + inc_e_ccpa*nens_ncep     
     
ALLOCATE (qmap_ensemble(nxa,nya,nens*9),&
    qmap_ensemble_ncep(nxa,nya,(nens_ncep+nmult_other)*9), &
    qmap_ensemble_cmc(nxa,nya,(nens_cmc+nmult_other)*9))

! ---- copy the input deterministic and ensemble data into the work array. 
!      Now modified to include the use of 8 surrounding points.

iktr(:,:) = 0
iktr_cmc(:,:) = 0
iktr_ncep(:,:) = 0

PRINT *,'quantile-mapped ensemble populating with CMC control...'
IF (icm_c_ccpa .eq. 1) THEN  ! copy in CMC control
   DO imult = 1, nmult_other
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (cmc_control_ccpa_x9(i9,ixa,jya) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_cmc(ixa,jya) = iktr_cmc(ixa,jya) + 1
                  qmap_ensemble(ixa,jya,iktr(ixa,jya)) = &
                       cmc_control_ccpa_x9(i9,ixa,jya)
                  qmap_ensemble_cmc(ixa,jya,iktr_cmc(ixa,jya)) = &
                       cmc_control_ccpa_x9(i9,ixa,jya)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO           ! imult
ENDIF

PRINT *,'quantile-mapped ensemble populating with NCEP control...'
IF (inc_c_ccpa .eq. 1) THEN  ! copy in NCEP control
   DO imult = 1, nmult_other
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (ncep_control_ccpa_x9(i9,ixa,jya) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_ncep(ixa,jya) = iktr_ncep(ixa,jya) + 1
                  qmap_ensemble(ixa,jya,iktr(ixa,jya)) = &
                       ncep_control_ccpa_x9(i9,ixa,jya)
                  qmap_ensemble_ncep(ixa,jya,iktr_ncep(ixa,jya)) = &
                       ncep_control_ccpa_x9(i9,ixa,jya)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO           ! imult
ENDIF

PRINT *,'quantile-mapped ensemble populating with CMC ensemble...'
IF (icm_e_ccpa .eq. 1) THEN ! copy in CMC ensemble
   DO imem = 1, nens_cmc
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (cmc_ensemble_ccpa_x9(i9,ixa,jya,imem) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_cmc(ixa,jya) = iktr_cmc(ixa,jya) + 1
                  qmap_ensemble(ixa,jya,iktr(ixa,jya)) = &
                       cmc_ensemble_ccpa_x9(i9,ixa,jya,imem)
                  qmap_ensemble_cmc(ixa,jya,iktr_cmc(ixa,jya)) = &
                       cmc_ensemble_ccpa_x9(i9,ixa,jya,imem)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO           ! imem
ENDIF

PRINT *,'quantile-mapped ensemble populating with NCEP ensemble...'
IF (inc_e_ccpa .eq. 1) THEN ! copy in NCEP ensemble
   DO imem = 1, nens_ncep
      DO jya = 1, nya
         DO ixa = 1, nxa
            DO i9 = 1,9
               IF (ncep_ensemble_ccpa_x9(i9,ixa,jya,imem) .ge. 0.0) THEN
                  iktr(ixa,jya) = iktr(ixa,jya) + 1
                  iktr_ncep(ixa,jya) = iktr_ncep(ixa,jya) + 1
                  qmap_ensemble(ixa,jya,iktr(ixa,jya)) = &
                       ncep_ensemble_ccpa_x9(i9,ixa,jya,imem)
                  qmap_ensemble_ncep(ixa,jya,iktr_ncep(ixa,jya)) = &
                       ncep_ensemble_ccpa_x9(i9,ixa,jya,imem)
               ENDIF
            END DO  ! i9
         END DO     ! ixa
      END DO        ! jya
   END DO
ENDIF

! ---- Now, dress the ensemble members over areas with valid training data with 
!      Gamma-distributed noise.  

!      Procedure is : 
!      (1) take care of NCEP and CMC ensemble probabilities alone, w/o noise added 
!      and with equal weighting.  
!      (2) sort the multi-model ensemble at each grid point from lowest to highest. 
!      (3) determine nonexceedance probability associated with each member's dressing info
!      (4) tally up the overall probability, not using equal weights for each member, but 
!      rather weighting them according to "closest_histogram" which describes the 
!      probability that the sorted (quantile-mapped, w/o noise added) member is the
!      best member.

PRINT *,'estimating MME probabilities with Gamma dressing'
prob_forecast_qmap(:,:) = 9999.0 ! initialize to a missing value flag.
DO jya = 1, nya
    DO ixa = 1, nxa

        nmembers_today = iktr(ixa,jya)
        ALLOCATE (member_weights(nmembers_today), ensemble_today(nmembers_today))

        ! ---- (1) first let's get the quantile-mapped probabilities from the qmap ensemble,
        !      with equal weighting

        ncount_ncep = 0
        ncount_cmc = 0
        ncount = 0
        
        IF (inc_e_ccpa .eq. 1) THEN
           DO imem = 1, iktr_ncep(ixa,jya)
              IF (qmap_ensemble_ncep(ixa,jya,imem) .ge. rthresh) &
                   ncount_ncep = ncount_ncep + 1
           END DO
           prob_forecast_qmap_NCEP(ixa,jya) = REAL(ncount_ncep) / REAL(iktr_ncep(ixa,jya))
        ELSE
           prob_forecast_qmap_NCEP(ixa,jya) = 9999.0
        ENDIF

        IF (icm_e_ccpa .eq. 1) THEN
           DO imem = 1, iktr_cmc(ixa,jya)
              IF (qmap_ensemble_cmc(ixa,jya,imem) .ge. rthresh) ncount_cmc = ncount_cmc + 1
           END DO
           prob_forecast_qmap_CMC(ixa,jya) = REAL(ncount_cmc) / REAL(iktr_cmc(ixa,jya))
        ELSE
           prob_forecast_qmap_CMC(ixa,jya) = 9999.0
        ENDIF

        IF (iktr(ixa,jya) .gt. 0) THEN
           DO imem = 1, iktr(ixa,jya)
              IF (qmap_ensemble(ixa,jya,imem) .ge. rthresh) ncount = ncount + 1
           END DO
           prob_forecast_qmap(ixa,jya) = REAL(ncount) / REAL(iktr(ixa,jya))
        ELSE
           prob_forecast_qmap(ixa,jya) = 9999.0
        ENDIF
        
        IF (ihavedata(ixa, jya) .eq. 1 .and. &
            (cdressing .eq. 'Y' .or. cdressing .eq. 'y'))  THEN  

            ! inside CONUS, where data shd have been qmapped
            ! and with dressing flag set to yes

            ! ---- now let's start the process of determining the weights to apply to each
            !      member.  We have previously determined the weights that would be used
            !      with an ensemble of CMC + NCEP deterministic plus ensemble member's data
            !      at just one grid point.  We are now potentially using up to 9x more
            !      data by using data from surrounding grid points, and it's also possible
            !      that sample size has been decremented should CMC or NCEP data be missing.
            !      Given the number of samples available at this grid point, let's readjust
            !      the weights, spreading them out over the available number of members.
        
            nmembers_today = iktr(ixa,jya)
            ensemble_today(1:nmembers_today) = qmap_ensemble(ixa,jya,1:nmembers_today)
            CALL rejigger_weights(nmembers, nmembers_today, closest_histogram, member_weights) 

            ! ---- sort the ensemble, lowest to highest
        
            CALL sort(nmembers_today, ensemble_today)
            !PRINT *, 'ensemble_today = ', ensemble_today(:)
        
            ! ---- tally the nonexceedance probability from the weighted sum of 
            !      nonexceedance probabilities generated from the Gamma dressing statistics
            !      parameters

            probsum_belowthresh = 0.0
            probsum_pzero = 0.0
            probsum_gamma = 0.0

            DO imem = 1, nmembers_today

                !  ---- find bounding values and interpolate to get pz, alpha, beta

                DO i = 1, n_amounts - 1
                    IF (ensemble_today(imem) .GE. ramt(i) .AND. &
                    ensemble_today(imem) .LT. ramt(i+1)) THEN
                        IF (gamma_shape(i+1) .GE. 0.0  .AND. gamma_shape(i) .GE. 0.0) THEN
                           weight = 1. - (ensemble_today(imem) - ramt(i)) / (ramt(i+1) - ramt(i))
                           IF (ensemble_today(imem) .eq. 0.0) THEN
                              cpop = climo_pop(ixa,jya)
                              DO icat = 1, n_climocats-1
                                 IF (icat .eq. 1) THEN
                                    clower = 0.0
                                 ELSE
                                    clower = climo_pop_thresholds(icat-1)
                                 ENDIF
                                 IF (cpop .ge. clower .and. cpop .lt. climo_pop_thresholds(icat)) THEN
                                    iclim = icat
                                    GOTO 2349
                                 ENDIF
                              END DO
                              iclim = n_climocats
2349                          CONTINUE
                              pz = fraczero_fclimpop(iclim)
                              !pz = weight*fraczero(i) + (1.-weight)*fraczero(i+1)
                           ELSE
                              pz = weight*fraczero(i) + (1.-weight)*fraczero(i+1)
                           ENDIF
                           alpha = weight*gamma_shape(i) + (1.-weight)*gamma_shape(i+1)
                           beta = weight*gamma_scale(i) + (1.-weight)*gamma_scale(i+1)
                        ELSE  ! precip so high that no reliable gamma params found; use last value
                           pz = fraczero(n_amounts)
                           alpha = gamma_shape(n_amounts)
                           beta = gamma_scale(n_amounts)
                        ENDIF
                    END IF
                END DO
                !PRINT *,'fraczero_fclimpop = ',fraczero_fclimpop
                !PRINT *,'climo_pop_thresholds = ',climo_pop_thresholds
                !PRINT *,'pz, alpha, beta, cpop = ',pz, alpha, beta, cpop

                ! ---- determine the total event probability from a weighted combination of
                !      probabilities of event exceedance from the dressed distributions for
                !      each ensemble member

                ksi = rthresh / beta
                CALL cumgam(ksi, alpha, cum, ccum)                      
                pgamma_lt_thresh = member_weights(imem)*(1.-pz)*cum
                pzero = member_weights(imem)*pz
                probsum_belowthresh = probsum_belowthresh + pzero + pgamma_lt_thresh
                probsum_pzero = probsum_pzero + pzero
                probsum_gamma = probsum_gamma + pgamma_lt_thresh

            END DO ! IMEM
            prob_forecast_unsmoothed(ixa,jya) = 1.0 - probsum_belowthresh
        ELSE 

            ! ihavedata(ixa,jya) = 0; generate probs from the input ens.

            ncount = 0
            DO imem = 1, iktr(ixa,jya)
                IF (qmap_ensemble(ixa,jya,imem) .GE. rthresh) ncount = ncount + 1
            END DO
            prob_forecast_unsmoothed(ixa,jya) = REAL(ncount) / REAL(iktr(ixa,jya))
            
        ENDIF

        DEALLOCATE (member_weights, ensemble_today)

   END DO ! IXA
END DO ! JY

DO jya = nya/3, nya/3
    DO ixa = nxa/5, nxa/2,5
       print 317,'ixa,jya,prob = ',ixa, jya, &
          prob_forecast_unsmoothed(ixa,jya)
317    format(a33,2(i3,1x),4(f8.5,1x))
    END DO
END DO

WRITE (6,fmt='(A)')' Quantile-mapped probabilities:'
WRITE (6,fmt='(4(A10,1X))')'MODEL','MIN','MAX','MEAN'
WRITE (6,fmt='(A10,1X,3(F10.5,1X))')'NCEP',minval(prob_forecast_qmap_NCEP),&
    maxval(prob_forecast_qmap_NCEP),sum(prob_forecast_qmap_NCEP)/(nxa*nya)
WRITE (6,fmt='(A10,1X,3(F10.5,1X))')'CMC',minval(prob_forecast_qmap_CMC),&
    maxval(prob_forecast_qmap_CMC),sum(prob_forecast_qmap_CMC)/(nxa*nya)
WRITE (6,fmt='(A)')' Quantile-mapped and Gamma-distribution dressed probabilities:'
WRITE (6,fmt='(A10,1X,3(F10.5,1X)/)')'COMBINED',minval(prob_forecast_qmap),&
    maxval(prob_forecast_qmap),sum(prob_forecast_qmap*REAL(ihavedata))/&
    SUM(REAL(ihavedata))

! ---- Deallocate the arrays
DEALLOCATE (qmap_ensemble, qmap_ensemble_cmc, qmap_ensemble_ncep)

RETURN
END SUBROUTINE ensemble_probs_x9_dressweight
