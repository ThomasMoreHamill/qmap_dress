! FILE: get_cdf_precip_anal_v2.f90
!  f2py -c -m get_cdf_precip_anal_v2 get_cdf_precip_anal_v2.f90
! ======================

SUBROUTINE get_cdf_precip_anal_v2(nthresh, nya, nxa, &
    thresh, apcp_anal, conusmask, CDFwork, icount, istat)
    
INTEGER, INTENT(IN) :: nthresh, nya, nxa
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
REAL, INTENT(IN), DIMENSION(nya,nxa) :: apcp_anal
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: conusmask

INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount
REAL*8, INTENT(OUT), DIMENSION(nthresh,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nya, nxa, thresh
!f2py intent(in) apcp_anal, conusmask
!f2py depend(nthresh) thresh
!f2py depend(nya,nxa) apcp_anal
!f2py depend(nya,nxa) conusmask
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

CDFwork = 0.
icount(:,:) = 0
rmaxprecip = maxval(apcp_anal)

!print *,'apcp_anal(nya/2,:)= ', apcp_anal(nya/2,:)
!print *,'conusmask(nya/2,:) = ', conusmask(nya/2,:)
!print *,'thresh(1:30) = ', thresh(1:30)

DO jya = 1, nya
    DO ixa = 1, nxa
        IF (conusmask(jya,ixa) .eq. 1 .and. rmaxprecip .ge. 0.0) THEN
            icount(jya,ixa) = icount(jya,ixa) + 1
            DO ithr = 1, nthresh
                IF (apcp_anal(jya,ixa) .LE. thresh(ithr)) &
                    CDFwork(ithr,jya,ixa) = CDFwork(ithr,jya,ixa) + 1.0
            END DO ! ithr
        END IF ! conusmask
    END DO ! ixa
END DO ! jya

IF (rmaxprecip .lt. 0.0) THEN
   istat = -1
ELSE
   istat = 1
ENDIF

RETURN
END SUBROUTINE get_cdf_precip_anal_v2



