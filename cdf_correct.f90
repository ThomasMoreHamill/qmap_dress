SUBROUTINE cdf_correct(nxa, nya, npct, thresh, ihavedata, forecast_cdf, analysis_cdf, forecast)

! ---- Perform the CDF bias correction

! ---- Input/Output Vaiables
integer, intent(in) :: nxa, nya, npct
real, intent(in), dimension(npct) :: thresh
integer*2, intent(in), dimension(nxa,nya) :: ihavedata
real*8, intent(in), dimension(nxa,nya,npct) :: forecast_cdf, analysis_cdf
real, intent(inout), dimension(nxa,nya) :: forecast

! ---- Local Vaiables
integer :: itest,jtest
real :: weight, weight2
real*8, dimension(npct) :: fcdf, acdf ! ---- set random seed from mean forecast
real, dimension(5) :: a95_to_a99, f95_to_f99

! ---- Initialize
itest=0
jtest=0
weight=0.0
weight2=0.0

print *,'at beginning of cdf_correct, analysis_cdf(479,199,:) = ',analysis_cdf(479,199,:)

! ---- Iterate over the gridpoints.  The CDF bias correction will occur for gridpoints
!      where ihavedata = 1 and the forecast >= 0.
do jya=1,nya
   do ixa=1,nxa

      if(ihavedata(ixa,jya).eq.1.and.forecast(ixa,jya).ge.0.0)then
         do ipct=1,npct-1

            if(forecast(ixa,jya).ge.thresh(ipct).and.&
               forecast(ixa,jya).lt.thresh(ipct+1))then

               ! ---- Determine the forecast CDF associated with this forecast amount.
               forecast_old = forecast(ixa,jya)
               weight=(forecast(ixa,jya)-thresh(ipct))/&
                      (thresh(ipct+1)-thresh(ipct))
               cdf_interpolated=forecast_cdf(ixa,jya,ipct)*(1.-weight)+&
                                forecast_cdf(ixa,jya,ipct+1)*weight

               ! ---- Set the new forecast to be the analysis value associated with 
               !      the same quantile of the cdf.  Only bother to do so for non-extreme values
               if(cdf_interpolated.lt.0.95)then
                  weight2=0.0
                  if(cdf_interpolated.lt.analysis_cdf(ixa,jya,1))then
                     forecast(ixa,jya)=0.0
                     goto 3000
                  endif
                  do icdf=1,npct-1
                     if(cdf_interpolated.ge.analysis_cdf(ixa,jya,icdf).and.&
                        cdf_interpolated.lt.analysis_cdf(ixa,jya,icdf+1))then
                        weight2=(cdf_interpolated-analysis_cdf(ixa,jya,icdf))/&
                                (analysis_cdf(ixa,jya,icdf+1)-analysis_cdf(ixa,jya,icdf))
                        forecast(ixa,jya)=thresh(icdf)*(1.-weight2)+thresh(icdf+1)*weight2

                        if (forecast(ixa,jya) .lt. 0.0) then
                           print *,'forecast before = ',forecast_old
                           print *,'cdf_interpolated = ',cdf_interpolated
                           print *,'icdf = ',icdf
                           print *,'bounding analysis cdf values = ',&
                                analysis_cdf(ixa,jya,icdf), analysis_cdf(ixa,jya,icdf+1)
                           print *,'bounding thresholds = ',thresh(icdf), thresh(icdf+1)
                           print *,'weight lower, upper = ', 1-weight2, weight2
                           stop
                        endif

                        goto 3000
                     end if
                  end do
               else 
                  ! ---- cdf_interpolated >= 0.95; apply Scheuerer regression analysis approach
                  !      from appendix A of Scheuerer and Hamill, MWR, 143, 4578-4596. The
                  !      underlying rationale is that quantile mapping produces potentially
                  !      especially unrealistic values at the extreme high percentiles, so the
                  !      regression approach should diminish this tendency.
                  fcdf(:) = forecast_cdf(ixa,jya,:)
                  acdf(:) = analysis_cdf(ixa,jya,:)
                  IF (ixa .eq. 479 .and. jya .eq. 199) THEN 
                     PRINT *,'acdf = ',acdf
                     PRINT *,'fcdf = ',fcdf
                     PRINT *,'thresh = ', thresh
                  ENDIF

                  ! ---- Find the forecast and analyzed values associated with the
                  !      95th thru 99th percentiles of the distribution.
                  call get_95_to_99(npct, acdf, fcdf, thresh, a95_to_a99, f95_to_f99)

                  IF (ixa .eq.479 .and. jya .eq. 199) THEN
                     PRINT *,'a95_to_a99 = ',a95_to_a99
                     PRINT *,'f95_to_f99 = ',f95_to_f99

                  ENDIF




                  ! ---- Determine the regression slope associated with the correction
                  !      following Schuerer's method.  Apply that assuming regression intercept
                  !      is set by the analyzed value at 95th percentile.
                  if(sum((f95_to_f99(2:5)-f95_to_f99(1))**2).gt.0.0)then
                     slope=SUM((a95_to_a99(2:5)-a95_to_a99(1))*(f95_to_f99(2:5)-f95_to_f99(1)))/&
                           sum((f95_to_f99(2:5)-f95_to_f99(1))**2)
                     forecast(ixa,jya)=a95_to_a99(1)+(forecast(ixa,jya)-f95_to_f99(1))*slope
                  else
                     forecast(ixa,jya)=forecast(ixa,jya)
                  endif


                        if (forecast(ixa,jya) .lt. 0.0) then
                           print *,'ixa, jya, mask = ',ixa, jya, ihavedata(ixa,jya)
                           print *,'forecast before = ',forecast_old
                           print *,'fcdf = ',fcdf
                           print *,'acdf = ',acdf
                           print *,'a95_to_a99 = ',a95_to_a99
                           print *,'f95_to_f99 = ',f95_to_f99
                           print *,'slope = ',slope
                           stop
                        endif



                  goto 3000
               ENDIF
            endif
         end do
      endif
      3000 continue

      ! Write information for a particular gridpoint (itest,jtest).
      if(itest.eq.ixa.and.jtest.eq.jya)then
          write(6,*),'inside cdf_correct, for point ',ixa,', ',jya
          write(6,*),'ipct, bounding thresh = ',icdf, thresh(ipct), thresh(ipct+1)
          write(6,*),'weight = ',weight
          write(6,*),'cdf_interpolated = ',cdf_interpolated
          write(6,*),'forecast = ',forecast(ixa,jya)
      endif

   end do
end do

return
end subroutine cdf_correct
