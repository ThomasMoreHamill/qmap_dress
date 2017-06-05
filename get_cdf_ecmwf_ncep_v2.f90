! FILE: get_cdf_ecmwf_ncep_v2.f90
!  f2py -c -m get_cdf_ecmwf_ncep_v2 get_cdf_ecmwf_ncep_v2.f90
! ================================================================================

SUBROUTINE get_cdf_ecmwf_ncep_v2(nthresh, nmembers, nya, nxa, &
    thresh, apcp_fcst_ens, conusmask, CDFwork, icount, istat)
    
INTEGER, INTENT(IN) :: nthresh, nmembers, nya, nxa
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
REAL, INTENT(IN), DIMENSION(nmembers, nya,nxa) :: apcp_fcst_ens
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: conusmask

INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount
REAL*8, INTENT(OUT), DIMENSION(nthresh,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nmembers, nya, nxa, thresh
!f2py intent(in) apcp_fcst_ens, conusmask
!f2py depend(nthresh) thresh
!f2py depend(nmembers,nya,nxa) apcp_fcst_ens
!f2py depend(nya,nxa) conusmask
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

CDFwork = 0.
icount(:,:) = 0
rminv = minval(apcp_fcst_ens)
rmaxv = maxval(apcp_fcst_ens)
!print *,'min,max value of forecast = ',rminv, rmaxv
!print *,'conusmask(nya/2,1:nxa:5) = ',conusmask(nya/2,1:nxa:5)
!print *,'thresh(1:16) =',thresh(1:16) 
!print *,'apcp_fcst_ens(:,nya/2,nxa/2) = ',apcp_fcst_ens(:,nya/2,nxa/2)
DO jya = 1, nya
    !PRINT *,'jya = ',jya,' of ',nya
    DO ixa = 1, nxa
        IF (conusmask(jya,ixa) .eq. 1 .and. rminv .ge. -98.) THEN
            DO imem = 1, nmembers
                icount(jya,ixa) = icount(jya,ixa) + 1
                !IF (jya .eq. nya/2 .and. ixa .eq. nxa/2) &
                !    print *,'imem, icount(nya/2,nxa/2) = ',imem, icount(nya/2,nxa/2) 
                DO ithr = 1, nthresh
                    IF (apcp_fcst_ens(imem,jya,ixa) .LE. thresh(ithr)) &
                        CDFwork(ithr,jya,ixa) = CDFwork(ithr,jya,ixa) + 1.0
                END DO ! ithr
            END DO ! imem
        END IF ! conusmask, good fcst data
    END DO ! ixa
END DO ! jya

IF (rminv .ge. 0.0) THEN
    istat = 1
ELSE
    istat = -1
ENDIF

RETURN
END SUBROUTINE get_cdf_ecmwf_ncep_v2



