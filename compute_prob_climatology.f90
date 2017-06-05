SUBROUTINE compute_prob_climatology(nxa, nya, ndays3mo, precip_anal_fine, &
   ihavedata, rthresh, climo_prob)

! compute the climatological probability.  
!
! coded Aug 2014 by Tom Hamill, tom.hamill@noaa.gov, (303) 497-3060

INTEGER, INTENT(IN) :: nxa, nya,ndays3mo
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata
REAL, INTENT(IN), DIMENSION(nxa,nya,ndays3mo) :: precip_anal_fine
REAL, INTENT(IN) :: rthresh
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: climo_prob
INTEGER :: ktr, ktr_prob
REAL :: rktr

DO ix = 1, nxa
   DO jy = 1, nya
      IF (ihavedata(ix,jy) .eq. 1) THEN
         ktr_prob = 0
         ktr = 0
         DO isamp = 1, ndays3mo
            IF (precip_anal_fine(ix,jy,isamp) .ge. rthresh) ktr_prob = ktr_prob + 1
            IF (precip_anal_fine(ix,jy,isamp) .ge. 0.0) ktr = ktr+1
         END DO
         rktr = REAL(ktr)
         IF (rktr .ge. 1) THEN
            climo_prob(ix,jy) = REAL(ktr_prob) / rktr
         ELSE
            climo_prob(ix,jy) = -99.99
         ENDIF
      ELSE
         climo_prob(ix,jy) = -99.99
      END IF
   END DO ! jy
END DO    ! ix

RETURN
END SUBROUTINE compute_prob_climatology
