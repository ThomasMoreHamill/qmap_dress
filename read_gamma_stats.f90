! f2py -c -m read_gamma_stats read_gamma_stats.f90

SUBROUTINE read_gamma_stats(nmembers, n_climocats, &
    nthreshes, infile, gamma_sum, gamma_ln_sum, gamma_sum_low, &
    gamma_ln_sum_low, gamma_sum_mod, gamma_ln_sum_mod, gamma_sum_high, &
    gamma_ln_sum_high, closest_histogram, closest_histogram_low, &
	closest_histogram_mod, closest_histogram_high, nzeros, nzeros_low, &
	nzeros_mod, nzeros_high, npositive, npositive_low, npositive_mod, &
	npositive_high, nzeros_fclim, npositive_fclim, fraction_zero_fclim, & 	
	istat)

INTEGER, INTENT(IN) :: nmembers
INTEGER, INTENT(IN) :: n_climocats ! thresholds for prob zero precip climo
INTEGER, INTENT(IN) :: nthreshes ! discretization of precip amt
CHARACTER*(*) :: infile

REAL, INTENT(OUT) :: istat

REAL*8, INTENT(OUT), DIMENSION(nthreshes+1) :: gamma_sum, gamma_ln_sum
REAL*8, INTENT(OUT), DIMENSION(nthreshes+1) :: gamma_sum_low, gamma_ln_sum_low
REAL*8, INTENT(OUT), DIMENSION(nthreshes+1) :: gamma_sum_mod, gamma_ln_sum_mod
REAL*8, INTENT(OUT), DIMENSION(nthreshes+1) :: gamma_sum_high, gamma_ln_sum_high
REAL, INTENT(OUT), DIMENSION(nmembers*9) :: closest_histogram
REAL, INTENT(OUT), DIMENSION(nmembers*9) :: closest_histogram_low
REAL, INTENT(OUT), DIMENSION(nmembers*9) :: closest_histogram_mod
REAL, INTENT(OUT), DIMENSION(nmembers*9) :: closest_histogram_high
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: nzeros
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: nzeros_low
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: nzeros_mod
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: nzeros_high
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: npositive
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: npositive_low
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: npositive_mod
REAL, INTENT(OUT), DIMENSION(nthreshes+1) :: npositive_high

REAL, INTENT(OUT), DIMENSION(n_climocats+1) :: nzeros_fclim, npositive_fclim
REAL, INTENT(OUT), DIMENSION(n_climocats+1) :: fraction_zero_fclim

LOGICAL exist

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

INQUIRE (file = infile, exist=exist)


IF (exist) THEN

	istat = 1
	PRINT *,'reading from ',infile
	OPEN(unit=15, file=infile, status='old', form='unformatted')
	!PRINT *,'writing nmembers, nthreshes+1'
	READ (15) nmembers_in, nthreshes_in_plus1, n_climocats_in

	READ (15) gamma_sum
	READ (15) gamma_ln_sum
	READ (15) closest_histogram
	READ (15) nzeros
	READ (15) npositive

	READ (15) gamma_sum_low
	READ (15) gamma_ln_sum_low
	READ (15) closest_histogram_low
	READ (15) nzeros_low
	READ (15) npositive_low

	READ (15) gamma_sum_mod
	READ (15) gamma_ln_sum_mod
	READ (15) closest_histogram_mod
	READ (15) nzeros_mod
	READ (15) npositive_mod

	READ (15) gamma_sum_high
	READ (15) gamma_ln_sum_high
	READ (15) closest_histogram_high
	READ (15) nzeros_high
	READ (15) npositive_high

	READ (15) nzeros_fclim
	READ (15) npositive_fclim
	
	READ (15) climo_pop_thresholds
	READ (15) gamma_threshes

	CLOSE (15)
	PRINT *,'Done reading'
	
ELSE
    istat = -1
ENDIF
				
RETURN
END SUBROUTINE read_gamma_stats
				
				