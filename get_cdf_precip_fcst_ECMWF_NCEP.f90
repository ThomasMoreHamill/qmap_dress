! FILE: get_cdf_precip_fcst_ECMWF_NCEP.f90
!  f2py -c -m get_cdf_precip_fcst_ECMWF_NCEP get_cdf_precip_fcst_ECMWF_NCEP.f90
! ================================================================================

SUBROUTINE get_cdf_precip_fcst_ECMWF_NCEP(nthresh, nmembers, nya, nxa, nsuppmax, &
    thresh, apcp_fcst_ens, xlocations, ylocations, nsupplemental, conusmask, &
    CDFwork, icount, istat)
    
INTEGER, INTENT(IN) :: nthresh, nmembers, nya, nxa, nsuppmax
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
REAL, INTENT(IN), DIMENSION(nmembers, nya,nxa) :: apcp_fcst_ens
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: xlocations
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: ylocations
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: nsupplemental, conusmask
INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount

REAL*8, INTENT(OUT), DIMENSION(nthresh,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nmembers, nya, nxa, nsuppmax, thresh
!f2py intent(in) apcp_fcst_ens, xlocations, ylocations, nsupplemental, conusmask
!f2py depend(nthresh) thresh
!f2py depend(nmembers,nya,nxa) apcp_fcst_ens
!f2py depend(nsuppmax,nya,nxa) xlocations, ylocations
!f2py depend(nya,nxa) nsupplemental, conusmask
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

CDFwork = 0.
icount(:,:) = 0
DO jya = 1, nya
    DO ixa = 1, nxa
        IF (conusmask(jya,ixa) .eq. 1 .and. minval(apcp_fcst_ens) .ge. 0.0) THEN
            istat = 1
            DO imem = 1, nmembers
                icount(jya,ixa) = icount(jya,ixa) + nsupplemental(jya,ixa)
                DO isupp = 1, nsupplemental(jya,ixa)
                    DO ithr = 1, nthresh
                        jy2 = ylocations(isupp,jya,ixa)
                        ix2 = xlocations(isupp,jya,ixa)
                        IF (apcp_fcst_ens(imem,jy2,ix2) .LE. thresh(ithr)) &
                            CDFwork(ithr,jya,ixa) = CDFwork(ithr,jya,ixa) + 1.0
                    END DO ! ithr
                END DO ! isupp
            END DO ! imem
        ELSE
            istat = -1
        END IF ! conusmask, good fcst data
    END DO ! ixa
END DO ! jya

RETURN
END SUBROUTINE get_cdf_precip_fcst_ECMWF_NCEP



