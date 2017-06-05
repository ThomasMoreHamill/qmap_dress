! FILE: get_cdf_cmc_v2.f90
!  f2py -c -m get_cdf_cmc_v2 get_cdf_cmc_v2.f90
! ================================================================================

SUBROUTINE get_cdf_cmc_v2(nthresh, nmembers, nya, nxa, &
    thresh, apcp_fcst_ens, conusmask, CDFwork, icount, istat)
    
INTEGER, INTENT(IN) :: nthresh, nmembers, nya, nxa
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
REAL, INTENT(IN), DIMENSION(nmembers, nya,nxa) :: apcp_fcst_ens
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: conusmask

INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount
REAL*8, INTENT(OUT), DIMENSION(nthresh,nmembers,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nmembers, nya, nxa, thresh
!f2py intent(in) apcp_fcst_ens,conusmask
!f2py depend(nthresh) thresh
!f2py depend(nmembers,nya,nxa) apcp_fcst_ens
!f2py depend(nya,nxa) conusmask
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

PRINT *,'get_cdf_cmc: nthresh, nmembers, nya, nxa, nsuppmax = ',&
    nthresh, nmembers, nya, nxa, nsuppmax
CDFwork = 0.
icount(:,:) = 0
rmv = minval(apcp_fcst_ens)

DO jya = 1, nya
    !PRINT *,'jya = ',jya,' of ',nya
    DO ixa = 1, nxa
        IF (conusmask(jya,ixa) .eq. 1 .and. rmv .ge. -98.) THEN
            icount(jya,ixa) = icount(jya,ixa) + 1
            DO imem = 1, nmembers
                DO ithr = 1, nthresh
                    IF (apcp_fcst_ens(imem,jya,ixa) .LE. thresh(ithr)) &
                        CDFwork(ithr,imem,jya,ixa) = CDFwork(ithr,imem,jya,ixa) + 1.0
                END DO ! ithr
            END DO ! imem
        END IF ! conusmask, good fcst data
    END DO ! ixa
END DO ! jya

istat = 1
IF (rmv .lt. 0.0) istat = -1
RETURN
END SUBROUTINE get_cdf_cmc_v2



