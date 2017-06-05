SUBROUTINE upscale_adjusted (nxa, nya, nxf, nyf, nnearest, ncount, ilist, jlist, &
     forecast_ccpa, forecast)

! --- having CDF bias-corrected the forecasts on the hi-res CCPA grid, now upscale them
!     to the forecast grid.  We need to do this to preserve the property that the 
!     forecast grid's precipitation is the area integral of the CCPA grid's precip.

INTEGER, INTENT(IN) :: nxa,nya, nxf, nyf
!INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask

INTEGER, DIMENSION(nxf,nyf) :: nnearest
INTEGER, DIMENSION(nxf,nyf,ncount) :: ilist, jlist
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast_ccpa
REAL, INTENT(INOUT), DIMENSION(nxf,nyf) :: forecast

PRINT *, 'max forecast value before = ', maxval(forecast)
DO jyf = 1, nyf
   DO ixf = 1, nxf
     
      ! ---- only process this forecast point if we know there are fine-scale
      !      data points nearest to it.
      IF (nnearest(ixf,jyf) .ge. 1) THEN
         fsum = 0.
         ktr = 0
         ! ---- find the average of all of the fine-scale points associated
         !      with this forecast point, and reset the forecast point to the 
         !      mean of these.
!         IF (nnearest(ixf,jyf) .gt. 25) THEN
!            nnearest(ixf,jyf)=25
!         ENDIF
         DO ipt = 1, nnearest(ixf,jyf)
            ixa = ilist(ixf,jyf,ipt)
            jya = jlist(ixf,jyf,ipt)
            IF (forecast_ccpa(ixa,jya) .ge. 0.0) THEN
               fsum = fsum + forecast_ccpa(ixa,jya)
               ktr = ktr+1
            ENDIF
         END DO
         IF (ktr .gt. 0) forecast(ixf,jyf) = fsum / REAL(ktr)
      ENDIF
   END DO
END DO
PRINT *, 'max forecast value after = ', maxval(forecast)
RETURN
END SUBROUTINE upscale_adjusted
