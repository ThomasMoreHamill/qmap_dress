SUBROUTINE get_85_to_99(npct, acdf, fcdf, thresh, a85_to_a99, f85_to_f99)

! purpose: find the precipitation amount, forecast and analyzed, associated
!    with the 95th to the 99th percentiles of the CDF distributions

INTEGER, INTENT(IN) :: npct
REAL*8, INTENT(IN), DIMENSION(npct) :: acdf, fcdf
REAL, INTENT(IN), DIMENSION(npct) :: thresh
REAL, INTENT(OUT), DIMENSION(15) :: a85_to_a99, f85_to_f99

REAL, DIMENSION(npct) :: cdf
REAL, DIMENSION(15) :: aorf

REAL, DIMENSION(15) :: tval
DATA tval /0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, &
   0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99/

DO itype = 1,2
   IF (itype .eq. 1) THEN
        cdf = acdf
   ELSE
        cdf = fcdf
   ENDIF
        
   DO ival = 1, 15
                
       ! --- handle case of dry CDFs, where > 85% are zero
        
       IF (tval(ival) .lt. cdf(1)) THEN
           aorf(ival) = 0.0
           GOTO 3245
       END IF
                
       ! ---- normal case
                
       DO ipct = 1, npct-1
          IF (tval(ival) .gt. cdf(ipct)  .and. tval(ival) .le. cdf(ipct+1))  THEN
             weight = (tval(ival)-cdf(ipct)) / (cdf(ipct+1)-cdf(ipct))
             aorf(ival) = (1.-weight)*thresh(ipct) + weight*thresh(ipct+1)
             GOTO 3245
          ENDIF
       END DO
                
       3245 CONTINUE
   END DO
        
   IF (itype .eq. 1) THEN
        a95_to_a99 = aorf
   ELSE
        f95_to_f99 = aorf
   ENDIF
        
END DO
RETURN

END SUBROUTINE get_95_to_99
