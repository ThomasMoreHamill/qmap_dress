SUBROUTINE tally_gamma_stats(nine, nxa, nya, nmembers, n_climocats, nthreshes, &
	thresh_low, thresh_high, ensemble_x9, analysis, conusmask, &
	outfilename, gamma_threshes, climo_prob, climo_pop_thresholds, &
	istat)

INTEGER, INTENT(IN) :: nine, nxa, nya, nmembers
INTEGER, INTENT(IN) :: n_climocats ! thresholds for prob zero precip climo
INTEGER, INTENT(IN) :: nthreshes ! discretization of precip amt
REAL, INTENT(IN) :: thresh_low, thresh_high ! for breaking up closest histogram and such
REAL, INTENT(IN), DIMENSION (nine, nxa, nya, nmembers) :: ensemble_x9
REAL, INTENT(IN), DIMENSION (nxa, nya) :: analysis
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
CHARACTER*(*), INTENT(IN) :: outfilename
REAL, INTENT(IN), DIMENSION (nxa, nya) :: climo_prob
REAL, INTENT(IN), DIMENSION(n_climocats-1) :: climo_pop_thresholds
REAL, INTENT(IN), DIMENSION(nthreshes) :: gamma_threshes

REAL, INTENT(OUT) :: istat

REAL*8, DIMENSION(nthreshes+1) :: gamma_sum, gamma_ln_sum
REAL*8, DIMENSION(nthreshes+1) :: gamma_sum_low, gamma_ln_sum_low
REAL*8, DIMENSION(nthreshes+1) :: gamma_sum_mod, gamma_ln_sum_mod
REAL*8, DIMENSION(nthreshes+1) :: gamma_sum_high, gamma_ln_sum_high
REAL, DIMENSION(nmembers*9) :: closest_histogram
REAL, DIMENSION(nmembers*9) :: closest_histogram_low
REAL, DIMENSION(nmembers*9) :: closest_histogram_mod
REAL, DIMENSION(nmembers*9) :: closest_histogram_high
REAL, DIMENSION(nthreshes+1) :: nzeros
REAL, DIMENSION(nthreshes+1) :: nzeros_low
REAL, DIMENSION(nthreshes+1) :: nzeros_mod
REAL, DIMENSION(nthreshes+1) :: nzeros_high
REAL, DIMENSION(nthreshes+1) :: npositive
REAL, DIMENSION(nthreshes+1) :: npositive_low
REAL, DIMENSION(nthreshes+1) :: npositive_mod
REAL, DIMENSION(nthreshes+1) :: npositive_high

REAL, DIMENSION(n_climocats+1) :: nzeros_fclim, npositive_fclim
REAL, DIMENSION(n_climocats+1) :: fraction_zero_fclim

REAL rdum
INTEGER idum

gamma_sum(:) = 0.0 
gamma_ln_sum(:) = 0.0
gamma_sum_low(:) = 0.0 
gamma_ln_sum_low(:) = 0.0
gamma_sum_mod(:) = 0.0
gamma_ln_sum_mod(:) = 0.0
gamma_sum_high(:) = 0.0
gamma_ln_sum_high(:) = 0.0
closest_histogram(:) = 0.0
closest_histogram_low(:) = 0.0
closest_histogram_mod(:) = 0.0
closest_histogram_high(:) = 0.0
nzeros(:) = 0.0
nzeros_low(:) = 0.0
nzeros_mod(:) = 0.0
nzeros_high(:) = 0.0
npositive(:) = 0.0
npositive_low(:) = 0.0
npositive_mod(:) = 0.0
npositive_high(:) = 0.0
fraction_zero_fclim(:) = 0.0
nzeros_fclim(:) = 0.0
npositive_fclim(:) = 0.0

! --- loop thru and process all points in the CONUS

!print *, 'tally_gamma_stats'
!print *, 'nine, nxa, nya, nmembers = ',nine, nxa, nya, nmembers
!print *, 'nthreshes, gamma_threshes = ', nthreshes, gamma_threshes
istat = 1

rdum = 0.0
DO ixa = 1, nxa
	rdum = rdum + ensemble_x9(1,ixa,nya/2,1)
END DO
idum = NINT(rdum)
!print *,'idum = ',idum
rm = MINVAL(ensemble_x9)
rma = MAXVAL(analysis)
!print *,'rm, rma = ',rm,rma
IF (rm .lt. -98. .and. rma .lt. -98.) THEN
	istat = -1
	PRINT *,'identified bad ensemble or analysis data.  rm, rma = ',rm, rma
	PRINT *,'setting istat to -1 and exiting.'
ELSE
	DO jya = 1, nya
		DO ixa = 1, nxa
	!DO jya = nya/2, nya/2
	!	DO ixa = 251, 251  ! 1, nxa, 5
			
			!PRINT *,'***** ixa, jya = ', ixa, jya
			IF (conusmask(ixa, jya) .eq. 1 .and. analysis(ixa,jya) .ge. 0.0) THEN
			
				! --- determine which member is the closest to the
				!     analyzed and how many members have values lower 
				!     than or equal to the analyzed value
			
				rclosest = 9999.
				rsum = 0.0
				iktr = 0
				DO inine = 1, 9
					DO imem = 1, nmembers
						!PRINT *,'inine, imem, ens = ', inine, imem, ensemble_x9(inine,ixa,jya,imem)
						IF (ensemble_x9(inine,ixa,jya,imem) .gt. 0.0) THEN
							rsum = rsum + ensemble_x9(inine,ixa,jya,imem)
							iktr = iktr+1
						ENDIF
					END DO
				END DO
				emean = rsum / REAL(iktr)
				!PRINT *,'emean = ', emean
								
				a = analysis(ixa, jya)

				! ---- find the member that is the closest

				inineclosest = 1
				imemclosest = 1
				eclosest = 9999.
				DO inine = 1, nine
					DO imem = 1, nmembers
						e = ensemble_x9(inine,ixa,jya,imem)
						diff = ABS(a - e)
						IF (diff .lt. rclosest .and. e .gt. -99) THEN
							rclosest = diff
							eclosest = e
							inineclosest = inine
							imemclosest = imem
						ENDIF
					END DO
				END DO

				! ---- determine how many other members are lower than the
				!      closest member, and how many are equal
				
				ibelow = 0
				iequal = 0
				DO inine = 1, nine
					DO imem = 1, nmembers
						e = ensemble_x9(inine,ixa,jya,imem)
						IF (inine .ne. inineclosest .and. imem .ne. imemclosest) THEN
							IF (e .lt. a) ibelow = ibelow + 1
							IF (e .eq. a) iequal = iequal + 1
						ENDIF
					END DO
				END DO

				! --- determine the closest_histogram rank 
				
				IF (iequal .eq. 0) THEN
					iclosest = ibelow + 1			
				ELSE
					r = ran3(idum) * iequal
					ir = INT(r)
					IF (ir .gt. iequal) ir = iequal
					iclosest = ibelow + iequal + 1
				ENDIF
				
				!PRINT *,'histogram rank iclosest = ', iclosest
				
				closest_histogram(iclosest) = closest_histogram(iclosest) + 1
				IF (emean .lt. thresh_low) THEN
					closest_histogram_low(iclosest) = &
						closest_histogram_low(iclosest) + 1
				ELSE IF (emean .ge. thresh_low .and. emean .lt. thresh_high) THEN
					closest_histogram_mod(iclosest) = &
						closest_histogram_mod(iclosest) + 1
				ELSE IF (emean .ge. thresh_high) THEN
					closest_histogram_high(iclosest) = &
						closest_histogram_high(iclosest) + 1
				ENDIF
				
				! ---- determine the gamma amount threshold that is closest
				!      to this ensemble member
				
				ipclosest = 1
				dmin = 9999.
				DO icat = 1, nthreshes
					diff = ABS(gamma_threshes(icat) - eclosest)
					!PRINT *,'diff, gamma_threshes(icat), icat, eclosest = ',diff, gamma_threshes(icat), icat, eclosest 
					IF (diff .lt. dmin) THEN
						dmin = diff
						ipclosest = icat
					ENDIF
				END DO
				!PRINT *,'ipclosest, gamma_threshes(ipclosest), eclosest = ',ipclosest, gamma_threshes(ipclosest), eclosest

				IF (a .gt. 0) THEN
					
					! ---- tally up the sum of the positive values of the analyzed
					!      for this closest precipitation category and of the
					!      ln of the analyzed value
				
					gamma_sum(ipclosest) = gamma_sum(ipclosest) + a
					gamma_ln_sum(ipclosest) = gamma_ln_sum(ipclosest) + &
						alog(a)
					npositive(ipclosest) = npositive(ipclosest) + 1
					!PRINT *,'incrementing npositive ',ipclosest
					IF (emean .lt. thresh_low) THEN
						gamma_sum_low(ipclosest) = gamma_sum_low(ipclosest) + a
						gamma_ln_sum_low(ipclosest) = gamma_ln_sum_low(ipclosest) + &
							alog(a)
						npositive_low(ipclosest) = npositive_low(ipclosest) + 1
					ELSE IF (emean .ge. thresh_low .and. emean .lt. thresh_high) THEN
						gamma_sum_mod(ipclosest) = gamma_sum_mod(ipclosest) + a
						gamma_ln_sum_mod(ipclosest) = gamma_ln_sum_mod(ipclosest) + &
							alog(a)
						npositive_mod(ipclosest) = npositive_mod(ipclosest) + 1
					ELSE IF (emean .ge. thresh_high) THEN
						gamma_sum_high(ipclosest) = gamma_sum_high(ipclosest) + a
						gamma_ln_sum_high(ipclosest) = gamma_ln_sum_high(ipclosest) + &
							alog(a)
						npositive_high(ipclosest) = npositive_high(ipclosest) + 1
					ENDIF
					
				ELSE
					nzeros(ipclosest) = nzeros(ipclosest) + 1
					!PRINT *,'incrementing zeros',ipclosest
					IF (emean .lt. thresh_low) THEN
						nzeros_low(ipclosest) = nzeros_low(ipclosest) 
					ELSE IF (emean .ge. thresh_low .and. emean .lt. thresh_high) THEN
						nzeros_mod(ipclosest) = nzeros_mod(ipclosest)
					ELSE IF (emean .ge. thresh_high) THEN
						nzeros_high(ipclosest) = nzeros_high(ipclosest)
					ENDIF
				ENDIF	
				
				!PRINT *, 'gamma_sum(ipclosest) = ',gamma_sum(ipclosest) 
				!PRINT *, 'gamma_ln_sum(ipclosest) = ',gamma_ln_sum(ipclosest) 
				!PRINT *, 'npositive(ipclosest), nzeros(ipclosest) = ',npositive(ipclosest) , nzeros(ipclosest)
				
				CALL find_climo_category (n_climocats-1, &
					climo_pop_thresholds, climo_prob(ixa,jya), iclim)
					
				!PRINT *, 'climo_prob(ixa,jya) = ', climo_prob(ixa,jya)
				!PRINT *, 'climo_pop_thresholds = ', climo_pop_thresholds
				!PRINT *, 'iclim = ', iclim
					
				IF (eclosest .eq. 0.0) THEN
					IF (a .eq. 0.0) THEN
					    nzeros_fclim(iclim) = nzeros_fclim(iclim) + 1
					ELSE 
					    npositive_fclim(iclim) = npositive_fclim(iclim) + 1
					ENDIF
				ENDIF
				
			ENDIF ! conusmask
		END DO ! ixa
	END DO ! jya
ENDIF ! rm .ge. -98.

IF (istat .ne. -1) THEN

	! ---- finally set the probability mass at zero in the situation where the forecast

	DO iclim = 1, n_climocats+1
    	ntot = nzeros_fclim(iclim) + npositive_fclim(iclim)
    	IF (ntot .gt. 100) THEN
        	fraction_zero_fclim(iclim) = REAL(nzeros_fclim(iclim)) / REAL(ntot)
    	ENDIF   
    	!PRINT *,'iclim, nzeros, npositive, fraction_zero_fclim, thresh = ',&
    	!   iclim, nzeros_fclim(iclim), npositive_fclim(iclim), &
    	!   fraction_zero_fclim(iclim), climo_pop_thresholds(iclim)       
	END DO

	!PRINT *, 'gamma_sum = ',gamma_sum 
	!PRINT *, 'gamma_ln_sum = ', gamma_ln_sum
	!PRINT *, 'closest_histogram = ',closest_histogram
	!PRINT *, 'nzeros = ', nzeros
	!PRINT *, 'npositive = ',npositive

	! ---- write the output to file

	PRINT *,'writing to ',outfilename
	OPEN(unit=15, file=outfilename, status='replace', form='unformatted')
	!PRINT *,'writing nmembers, nthreshes+1'
	WRITE (15) nmembers, nthreshes+1, n_climocats

	!PRINT *,'writing overall'
	WRITE (15) gamma_sum
	WRITE (15) gamma_ln_sum
	WRITE (15) closest_histogram
	WRITE (15) nzeros
	WRITE (15) npositive

	!PRINT *,'writing low'
	WRITE (15) gamma_sum_low
	WRITE (15) gamma_ln_sum_low
	WRITE (15) closest_histogram_low
	WRITE (15) nzeros_low
	WRITE (15) npositive_low

	!PRINT *,'writing mod'
	WRITE (15) gamma_sum_mod
	WRITE (15) gamma_ln_sum_mod
	WRITE (15) closest_histogram_mod
	WRITE (15) nzeros_mod
	WRITE (15) npositive_mod

	!PRINT *,'writing high'
	WRITE (15) gamma_sum_high
	WRITE (15) gamma_ln_sum_high
	WRITE (15) closest_histogram_high
	WRITE (15) nzeros_high
	WRITE (15) npositive_high

	WRITE (15) nzeros_fclim
	WRITE (15) npositive_fclim
	
	WRITE (15) climo_pop_thresholds
	WRITE (15) gamma_threshes

	CLOSE (15)
	PRINT *,'Done writing'
	
ENDIF
				
RETURN
END SUBROUTINE tally_gamma_stats
				
				