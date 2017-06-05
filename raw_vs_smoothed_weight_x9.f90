SUBROUTINE raw_vs_smoothed_weight_x9(nxa,nya,topo_eighth,ihavedata,raw_weight)

! ---- we want a linear combination of the analog ensemble (raw) probability and the 
!      Savitzky-Golay smoothed probabilities, weighting more toward the former
!      in complex terrain and latter in smooth terrain.  Determine the weight
!      to apply to the raw analog ensemble.

INTEGER, INTENT(IN) :: nxa, nya
REAL, INTENT(IN), DIMENSION(nxa,nya) :: topo_eighth
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya):: ihavedata
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

         weight = 0.6*(sqrt_stdev - 13.) 
         IF (weight .lt. 0.0) weight = 0.0
         IF (weight .ne. weight) weight = 0.0
         IF (weight .gt. 0.6) weight = 0.6
         raw_weight(i,j) = weight
      ELSE   ! May 2016 modification, set weight to 0 outside CONUS.
         raw_weight(i,j) = 0.0     ! intended to force used of smoothed data along edges later in sgolay_smooth_x9
      ENDIF
   END DO
END DO

WRITE(6,*)'Finished with raw_vs_smoothed_weight_x9'

RETURN
END SUBROUTINE raw_vs_smoothed_weight_x9
