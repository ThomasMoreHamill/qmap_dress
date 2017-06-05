SUBROUTINE raw_vs_smoothed_weight(nxa, nya, nxf, nyf, rlons_fcst, rlats_fcst, &
   rlons_anal, rlats_anal, topo_eighth, ihavedata, raw_weight)

! ---- we want a linear combination of the analog ensemble (raw) probability and the 
!      Savitzky-Golay smoothed probabilities, weighting more toward the former
!      in complex terrain and latter in smooth terrain.  Determine the weight
!      to apply to the raw analog ensemble.

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf

INTEGER*2, INTENT(IN), DIMENSION(nxa,nya):: ihavedata
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlons_fcst, rlats_fcst
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal, topo_eighth
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: raw_weight

REAL work(9), wmean, workdev(9)
REAL stdev2d(nxa,nya)

! ---- determine the local standard deviation of the terrain about its mean and weight
DO i = 1, nxa
   imin = MAX(1,i-1)
   imax = MIN(nxa,i+1)
   DO j = 1, nya
      IF (ihavedata(i,j) .eq. 1) THEN
         jmin = MAX(1,j-1)
         jmax = MIN(nya,j+1)
         ktr = 0
         DO i2 = imin,imax
            DO j2 = jmin, jmax
               ktr = ktr+1
               work(ktr) = topo_eighth(i2,j2)
            END DO
         END DO
         wmean = SUM(work(1:ktr))/REAL(ktr)
         workdev(1:ktr) = work(1:ktr) - wmean
         var = 0.
         DO k = 1, ktr
            var = var + workdev(k)**2
         END DO
         stdev = SQRT(var/REAL(ktr-1))
         sqrt_stdev = SQRT(stdev)
         stdev2d(i,j) = sqrt_stdev

         !weight = 0.8*(sqrt_stdev - 10.) 
         !IF (weight .lt. 0.0) weight = 0.0
         !IF (weight .ne. weight) weight = 0.0
         !IF (weight .gt. 0.8) weight = 0.8

         weight = 0.6*(sqrt_stdev - 13.) 
         IF (weight .lt. 0.0) weight = 0.0
         IF (weight .ne. weight) weight = 0.0
         IF (weight .gt. 0.6) weight = 0.6

         raw_weight(i,j) = weight
      ELSE
         raw_weight(i,j) = 1.0
      ENDIF
   END DO
END DO
PRINT *,'finished with raw_vs_smoothed_weight'

RETURN
END SUBROUTINE raw_vs_smoothed_weight
