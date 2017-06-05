SUBROUTINE cdf_correct_nonoise_x9(nxa, nya, npct, nstride, &
    thresh, ihavedata, forecast_cdf, analysis_cdf, &
    forecast, forecast_x9)

! --- Perform the CDF bias correction; here we use 8 surrounding 
!     grid points forecasts as well as the central point to increase 
!     sample size and account for position error.  
!
!     This version (Dec 2016) removes the adding of random Gaussian
!     noise to each member; the new procedure will, at a later 
!     point in the code, apply Gamma-distributed noise and 
!     weighting by ensemble member.  


INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions
INTEGER, INTENT(IN) :: npct ! number of percentiles (quantiles) where CDF is stored
INTEGER, INTENT(IN) :: nstride ! how many grid points to skip with the 
    ! 3x3 array of ensemble data used here
REAL, INTENT(IN), DIMENSION(npct) :: thresh ! list of precipitation thresholds
    ! where CDF calculated
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata ! conus mask
REAL*8, INTENT(IN), DIMENSION(nxa,nya,npct) :: forecast_cdf, analysis_cdf
    ! forecast and analyzed CDFs for each model grid point
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast
    ! input forecast

REAL, INTENT(OUT), DIMENSION(9,nxa,nya) :: forecast_x9

REAL :: weight, weight2, rmean
LOGICAL first
REAL*8, DIMENSION(npct) :: fcdf, acdf ! ---- set random seed from mean forecast
REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99

first = .True.
!forecast_x9(:,:,:) = -99.99 ! set to missing
idum = -99345 + INT(SUM(forecast))  ! set random seed for gaussian number generation

! ---- Not for points inside conus mask, apply the procedure of doing quantile
!      mapping using the forecast at (i,j), but also surrounding locations.
!      also add random number to the input forecast quantile to also account for the
!      tendency of the ensemble to be over-certain of its amount.


weight = 0.0
DO jya = 1, nya
    DO ixa = 1, nxa
        IF (ihavedata(ixa,jya) .le. 0) THEN
            
            ! ---- If outside the conus, then make the output forecast array simply replicates
            !      of the input forecast array
            
            forecast_x9(:,ixa,jya) = forecast(ixa,jya)
        ELSE
            
            ! ---- Loop thru the 9 grid points with  (ixa,jya) in center
            
            forecast_x9(:,ixa,jya) = -99.99
            ktr = 0            
            DO jyn = jya-nstride, jya+nstride, nstride
                DO ixn = ixa-nstride, ixa+nstride, nstride
                    
                    ibadcdf = 0
                    ktr = ktr+1
                    
                    ! ---- Only process if this point is legitimately within the domain and
                    !      has reasonable value.
                    
                    !IF (ixa .eq. 479 .and. jya .eq. 199) &
                    !     PRINT *,'ixn, jyn, ixa, jya = ',ixn, jyn, ixa, jya

                    IF ((jyn .ge. 1 .and. jyn .le. nya) .and. (ixn .ge.1 .and. ixn .lt. nxa)) THEN

                        !IF (ixa .eq. 479 .and. jya .eq. 199) &
                        !     PRINT *,'ihavedata(ixn,jyn), forecast(ixn,jyn) = ',&
                        !     ihavedata(ixn,jyn), forecast(ixn,jyn)

                        IF (ihavedata(ixn,jyn) .eq. 1 .and. forecast(ixn,jyn) .ge. 0.0) THEN

                            ! ---- Find loop thru the precipitation thresholds till we find ones
                            !      that bound today's forecast precip value

                            DO ipct = 1, npct-1
                                IF (forecast(ixn,jyn) .ge. thresh(ipct) .and. &
                                forecast(ixn,jyn).lt.thresh(ipct+1))THEN
                                
                                    ! ---- Determine the percentile in the CDF associated 
                                    !      with this forecast amount
                                    
                                    weight = (forecast(ixn,jyn) - thresh(ipct)) / &
                                        (thresh(ipct+1) - thresh(ipct))
                                    cdf_interpolated = &
                                        forecast_cdf(ixn,jyn,ipct)*(1.-weight) + &
                                        forecast_cdf(ixn,jyn,ipct+1)*weight

                                    !IF (ixa .eq. 479 .and. jya .eq. 199) &
                                    !     PRINT *,'weight, cdf_interpolated = ',weight, cdf_interpolated 
            
                                    ! ---- For non-extreme values, set the new forecast to be 
                                    !      the analysis value associated with the same quantile 
                                    !      of the cdf.
                                    
                                    IF (cdf_interpolated .lt. 0.95) THEN
                                        weight2 = 0.0
                                        IF (cdf_interpolated .le. analysis_cdf(ixa,jya,1)) THEN
                                            forecast_x9(ktr,ixa,jya)=0.0
                                            GOTO 3000
                                        ENDIF

                                        DO icdf = 1, npct-1
                                            IF (cdf_interpolated .ge. analysis_cdf(ixa,jya,icdf) .and. &       
                                            cdf_interpolated .lt. analysis_cdf(ixa,jya,icdf+1)) THEN
                                                weight2 = &
                                                    (cdf_interpolated - analysis_cdf(ixa,jya,icdf)) / &
                                                    (analysis_cdf(ixa,jya,icdf+1) - &
                                                    analysis_cdf(ixa,jya,icdf))
                                                forecast_x9(ktr,ixa,jya) = &
                                                    thresh(icdf)*(1.-weight2) + thresh(icdf+1)*weight2
                                                GOTO 3000
                                            END IF
                                        END DO
                                        forecast_x9(ktr,ixa,jya) = forecast(ixn,jyn)
                                    ELSE

                                        ! -- cdf_interpolated >= 0.95; apply Scheuerer regression 
                                        !   analysis approach from appendix A of Scheuerer and Hamill, 
                                        !   MWR, 143, 4578-4596. The underlying rationale is that 
                                        !   quantile mapping produces potentially especially unrealistic 
                                        !   values at the extreme high percentiles, so the
                                        !   regression approach should diminish this tendency.

                                        fcdf(:) = forecast_cdf(ixn,jyn,:)
                                        acdf(:) = analysis_cdf(ixa,jya,:)

                                        ! -- Find the forecast and analyzed values associated with
                                        !   the 95th thru 99th percentiles of the distribution.
                                        
                                        CALL get_95_to_99(npct, acdf, fcdf, thresh, &
                                            a95_to_a99, f95_to_f99)
                                        IF (SUM(a95_to_a99) .eq. 0.0 .or. SUM(f95_to_f99) .eq. 0.0 .or. &
                                        MAXVAL(a95_to_a99) .gt. 1.0 .or. MINVAL(a95_to_a99) .lt. 0.0 &
                                        .or. MAXVAL(f95_to_f99) .gt. 1.0 .or. MINVAL(f95_to_f99) .lt. 0.0) &
                                           ibadcdf = 1

                                        ! --  Determine the regression slope associated with the 
                                        !   correction following Schuerer's method.  Apply that 
                                        !   assuming regression intercept is set by the analyzed 
                                        !   value at 95th percentile.

                                        IF (SUM((f95_to_f99(2:5) - f95_to_f99(1))**2) .gt. 0.0) THEN
                                            slope = SUM((a95_to_a99(2:5)-a95_to_a99(1))* &
                                                (f95_to_f99(2:5)-f95_to_f99(1))) / &
                                                SUM((f95_to_f99(2:5)-f95_to_f99(1))**2)
                                            IF (ibadcdf .eq. 0) THEN
                                               forecast_x9(ktr,ixa,jya) = &
                                                   a95_to_a99(1)+(forecast(ixn,jyn)- &
                                                   f95_to_f99(1))*slope
                                            ELSE
                                               forecast_x9(ktr,ixa,jya) = forecast(ixn,jyn)
                                            ENDIF

                                            !IF (ixa .eq. 479 .and. jya .eq. 199) THEN 
                                            !   print *,'ibadcdf = ',ibadcdf
                                            !   print *,'a95_to_a99 = ',a95_to_a99
                                            !   print *,'f95_to_f99 = ',f95_to_f99
                                            !   print *,'slope = ',slope
                                            !   print *,'forecast_x9(ktr,ixa,jya) = ',forecast_x9(ktr,ixa,jya)
                                            !ENDIF
                                        ELSE
                                            forecast_x9(ktr,ixa,jya) = forecast(ixn,jyn)
                                        ENDIF

                                        GOTO 3000

                                    ENDIF ! cdf_interpoloated < 0.95
                                ENDIF  ! forecast >= thresh(ipct), < thresh(ipct)+1
                            END DO ! ipct
3000                        CONTINUE
                        END IF ! ihavedata, forecast>0
                    END IF ! jyn>1, jyn < nya, etc.
                END DO ! ixn
            END DO ! jyn
        END IF ! ihavedata
    ENDDO ! ixa
END DO ! jya

RETURN
END SUBROUTINE cdf_correct_nonoise_x9
