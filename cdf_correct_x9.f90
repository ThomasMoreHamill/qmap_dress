SUBROUTINE cdf_correct_x9(nxa,nya,npct,nstride,stdran,thresh,ihavedata,forecast_cdf,&
    analysis_cdf,forecast,forecast_x9)

! --- Perform the CDF bias correction; here we use 8 surrounding grid points forecasts as well
!     to increase sample size and account for position error

INTEGER, INTENT(IN) :: nxa, nya, npct, nstride
REAL, INTENT(IN), DIMENSION(npct) :: thresh
REAL, INTENT(IN) :: stdran
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata
REAL, INTENT(IN), DIMENSION(nxa,nya,npct) :: forecast_cdf, analysis_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast
REAL, INTENT(OUT), DIMENSION(9,nxa,nya) :: forecast_x9

REAL :: weight, weight2, rmean
INTEGER :: idum,itest,jtest

REAL*8, DIMENSION(npct) :: fcdf, acdf ! ---- set random seed from mean forecast
REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99

itest=0
jtest=0

!PRINT *,'inside cdf_correct_x9: nxa,nya,npct,nstride,stdran = ',nxa,nya,npct,nstride,stdran
!PRINT *,'  thresh = ', thresh
!PRINT *,'  ihavedata(1:nxa:10,nya/2) = ', ihavedata(1:nxa:10,nya/2) 
!PRINT *,'  forecast_cdf(3*nxa/4, nya/2,:) = ',forecast_cdf(3*nxa/4, nya/2,:)
!PRINT *,'  analysis_cdf(3*nxa/4, nya/2,:) = ',analysis_cdf(3*nxa/4, nya/2,:)
!print *,'  forecast(3*nxa/4, nya/2) = ',forecast(3*nxa/4, nya/2)


rmean=1000.*SUM(forecast)/REAL(nxa*nya)
idum=-1*INT(rmean)

forecast_x9(:,:,:) = -99.99 ! set to missing

idum = -99345 + INT(SUM(forecast_x9)/(nxa*nya*9))  ! set random seed for gaussian number generation

! ---- Not for points inside conus mask, apply the procedure of doing quantile
!      mapping using the forecast at (i,j), but also surrounding locations.
!      also add random number to the input forecast quantile to also account for the
!      tendency of the ensemble to be over-certain of its amount.
weight = 0.0
DO jya = 1, nya
   DO ixa = 1, nxa

      IF(ihavedata(ixa,jya).le.0) THEN

         ! If outside the conus, then make the output forecast array simply replicates
         ! of the input forecast array
         forecast_x9(:,ixa,jya) = forecast(ixa,jya)

      ELSE

         forecast_x9(:,ixa,jya) = -99.99
         ktr = 0

         ! Loop thru the 9 grid points with  (ixa,jya) in center
         DO jyn=jya-nstride,jya+nstride,nstride
            DO ixn=ixa-nstride,ixa+nstride,nstride

               ktr=ktr+1
            
               ! Test point
               
               IF(itest.ne.0.and.jtest.ne.0)then
                  write(6,*)'ihavedata(itest,jtest), forecast = ',ihavedata(itest,jtest),ixn,jyn,forecast(ixn,jyn)
               ENDIF

               ! Only process if this point is legitimately within the domain
               IF((jyn.ge.1.and.jyn.le.nya).and.(ixn.ge.1.and.ixn.lt.nxa))THEN

                  IF(ihavedata(ixn,jyn).eq.1.and.forecast(ixn,jyn).ge.0.0)THEN

                     ! Find loop thru the precipitation thresholds till we find ones
                     ! that bound today's forecast precip value

                     DO ipct=1,npct-1

                        IF(forecast(ixn,jyn).ge.thresh(ipct).and.&
                             forecast(ixn,jyn).lt.thresh(ipct+1))THEN

                           ! Determine the percentile in the CDF associated with this forecast amount
                           weight=(forecast(ixn,jyn)-thresh(ipct))/&
                                  (thresh(ipct+1)-thresh(ipct))
                           cdf_interpolated=forecast_cdf(ixn,jyn,ipct)*(1.-weight)+&
                                            forecast_cdf(ixn,jyn,ipct+1)*weight

                           ! For non-extreme values, set the new forecast to be the analysis 
                           ! value associated with the same quantile of the cdf.  
                           
                           IF (cdf_interpolated .lt. 0.85)THEN
                              weight2=0.0
                              IF(cdf_interpolated.lt.analysis_cdf(ixa,jya,1))THEN
                                 !ktr = ktr + 1
                                 forecast_x9(ktr,ixa,jya)=0.0
                                 GOTO 3000
                              ENDIF
                       
                              DO icdf=1,npct-1
                                 IF(cdf_interpolated.ge.analysis_cdf(ixa,jya,icdf).and.&
                                    cdf_interpolated.lt.analysis_cdf(ixa,jya,icdf+1))THEN
                                    weight2=(cdf_interpolated-analysis_cdf(ixa,jya,icdf))/&
                                            (analysis_cdf(ixa,jya,icdf+1)-analysis_cdf(ixa,jya,icdf))
                                    !ktr = ktr+1
                                    forecast_x9(ktr,ixa,jya)=thresh(icdf)*(1.-weight2)+thresh(icdf+1)*weight2
                                    GOTO 3000
                                 END IF
                              END DO
                       
                           ELSE

                              ! cdf_interpolated >= 0.95; apply Scheuerer regression analysis approach
                              ! from appendix A of Scheuerer and Hamill, MWR, 143, 4578-4596. The 
                              ! underlying rationale is that quantile mapping produces potentially
                              ! especially unrealistic values at the extreme high percentiles, so the
                              ! regression approach should diminish this tendency.
                        
                              fcdf(:) = forecast_cdf(ixn,jyn,:)
                              acdf(:) = analysis_cdf(ixa,jya,:)

                              ! Find the forecast and analyzed values associated with the 95th
                              ! thru 99th percentiles of the distribution.
                              
                              CALL get_95_to_99(npct, acdf, fcdf, thresh, a95_to_a99, f95_to_f99)

                              ! Determine the regression slope associated with the correction
                              ! following Schuerer's method.  Apply that assuming regression intercept
                              ! is set by the analyzed value at 95th percentile.
                       
                              IF (SUM((f95_to_f99(2:5) - f95_to_f99(1))**2).gt.0.0) THEN
                                 slope=SUM((a95_to_a99(2:5)-a95_to_a99(1))*(f95_to_f99(2:5)-f95_to_f99(1)))/&
                                       SUM((f95_to_f99(2:5)-f95_to_f99(1))**2)
                                 forecast_x9(ktr,ixa,jya)=a95_to_a99(1)+(forecast(ixn,jyn)-f95_to_f99(1))*slope
                              ELSE
                                 forecast_x9(ktr,ixa,jya)=forecast(ixn,jyn)
                              ENDIF 

                              GOTO 3000

                           ENDIF ! cdf_interpoloated < 0.95
                        ENDIF  ! forecast >= thresh(ipct), < thresh(ipct)+1
                     END DO ! ipct

3000                 CONTINUE

                     IF(itest.ne.0.and.jtest.ne.0)THEN
                        write(6,*)'Before:  ktr, forecast, qmapped f = ',& 
                                   ktr,forecast(ixn,jyn),forecast_x9(ktr,itest,jtest)
                     ENDIF

                     ! Add a scaled  ~N(0,1) random number to the quantile-mapped precipitation 
                     ! amount, somewhat akin to a best-member dressing of the ensemble 
                     ! (though the statistics of the best member are guessed at rather than 
                     ! formally estimated here).
                           
                     !stddev=stdran*forecast_x9(ktr,ixa,jya)

                     IF (forecast_x9(ktr,ixa,jya) .eq. 0.0) THEN
                        stddev = 0.
                     ELSE
                        stddev = 0.25 + stdran*forecast_x9(ktr,ixa,jya)
                        !stddev = 0. 
                     END IF

                     randev=gasdev(idum)
                     forecast_x9(ktr,ixa,jya)=forecast_x9(ktr,ixa,jya)+stddev*randev

                     IF(forecast_x9(ktr,ixa,jya).lt.0.0)forecast_x9(ktr,ixa,jya)=0.0

                     IF(itest.ne.0.and.jtest.ne.0)THEN
                        write(6,*)'After: stddev, randev, forecast  = ',& 
                                  stddev,randev,forecast_x9(ktr,itest,jtest)
                     ENDIF

                  END IF ! ihavedata, forecast>0
               END IF ! jyn>1, jyn < nya, etc.
            END DO ! ixn
         END DO ! jyn
      END IF ! ihavedata
   END DO ! ixa
END DO ! jya

!print *,'  forecast_x9(:, 3*nxa/4, nya/2) = ',forecast_x9(:,3*nxa/4, nya/2)

RETURN
END SUBROUTINE cdf_correct_x9
