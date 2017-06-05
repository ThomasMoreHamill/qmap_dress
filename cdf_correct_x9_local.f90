SUBROUTINE cdf_correct_x9_local(nxa, nya, npct, nstride, &
	thresh, conusmask, forecast_cdf, analysis_cdf,&
	forecast, forecast_x9)

! --- Perform the CDF bias correction; here we use 8 surrounding grid points forecasts as well
!     to increase sample size and account for position error

INTEGER, INTENT(IN) :: nxa, nya, npct, nstride
REAL, INTENT(IN), DIMENSION(npct) :: thresh
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa,nya,npct) :: forecast_cdf, analysis_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast
REAL, INTENT(OUT), DIMENSION(9,nxa,nya) :: forecast_x9

REAL :: weight, weight2, rmean
INTEGER :: idum,itest,jtest, jyn, ixn

REAL*8, DIMENSION(npct) :: fcdf, acdf ! ---- set random seed from mean forecast
REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99

itest=0
jtest=0


!print *,'cdf_correct_x9_local'
!print *,'nxa, nya, nstride = ', nxa, nya, nstride
!print *,'forecast_cdf(nxa/2, nya/2, :) = ',forecast_cdf(nxa/2, nya/2, :)
!print *,'analysis_cdf(nxa/2, nya/2, :) = ',analysis_cdf(nxa/2, nya/2, :)
!print *,'thresh = ', thresh
!PRINT *,'maxval(forecast) = ', maxval(forecast)
!PRINT *,'forecast(nxa/2-10, nxa/2+10, nya/2) = ',forecast(nxa/2-10:nxa/2+10, nya/2)

rmean = 1000.*SUM(forecast)/REAL(nxa*nya)
idum = -1*INT(rmean)

forecast_x9(:,:,:) = -99.99 ! set to missing

idum = -99345 + INT(SUM(forecast_x9)/(nxa*nya*9))  ! set random seed for gaussian number generation

! ---- Not for points inside conus mask, apply the procedure of doing quantile
!      mapping using the forecast at (i,j), but also surrounding locations.
!      also add random number to the input forecast quantile to also account for the
!      tendency of the ensemble to be over-certain of its amount.
weight = 0.0
DO jya = 1, nya
	DO ixa = 1, nxa
!DO jya = nya/2, nya/2
	!DO ixa = 3*nxa/4, 3*nxa/4
	!DO ixa = nxa/2, nxa/2
		
		IF (conusmask(ixa,jya) .le. 0) THEN
         	! If outside the conus, then make the output forecast array simply replicates
         	! of the input forecast array
         	forecast_x9(:,ixa,jya) = forecast(ixa,jya)
      	ELSE
         	forecast_x9(:,ixa,jya) = -99.97
         	ktr = 0

         	! Loop thru the 9 grid points with  (ixa,jya) in center

         	DO jyn = jya-nstride, jya+nstride, nstride
            	DO ixn = ixa-nstride, ixa+nstride, nstride
					
					ktr = ktr + 1
					!print *,'ktr, jyn, ixn = ',ktr,jyn,ixn

               	 	! Only process if this point is legitimately within the domain
               	 	IF ( (jyn.ge.1 .and. jyn .le. nya) .and. &
			   	 	(ixn .ge. 1 .and. ixn .lt. nxa)) THEN

						!PRINT *, 'conusmask(ixn,jyn), forecast(ixn,jyn) = ',&
						!	conusmask(ixn,jyn), forecast(ixn,jyn)
						IF (forecast(ixn,jyn) .eq. 0.0) THEN
							forecast_x9(ktr,ixa,jya) = forecast(ixn,jyn)
							GOTO 3000
						ENDIF
						
                  		IF (conusmask(ixn,jyn) .eq. 1 .and. forecast(ixn,jyn) .ge. 0.0) THEN

                     		! ---- Find loop thru the precipitation thresholds till we find ones
                     	    !      that bound today's forecast precip value

                     	   	DO ipct = 1, npct-1

                        		IF (forecast(ixn,jyn) .ge. thresh(ipct) .and. &
                             	forecast(ixn,jyn) .lt. thresh(ipct+1)) THEN

                           	 		!  ---- Determine the percentile in the 
									!       CDF associated with this forecast amount
									
                           		 	weight = (forecast(ixn,jyn)-thresh(ipct)) / &
                                  		(thresh(ipct+1)-thresh(ipct))
                           			cdf_interpolated = forecast_cdf(ixn,jyn,ipct)*(1.-weight)+&
                                    	forecast_cdf(ixn,jyn,ipct+1)*weight
									!print *,'weight, cdf_interpolated = ',weight, cdf_interpolated

                           			! ---- For non-extreme values, set the new forecast to be the analysis 
                           		 	! 	   value associated with the same quantile of the cdf.  
                           
                           		 	IF (cdf_interpolated .lt. 0.95)THEN
                              			weight2 = 0.0
                              			IF (cdf_interpolated .lt. analysis_cdf(ixa,jya,1)) THEN
                                 			!ktr = ktr + 1
                                 		    forecast_x9(ktr,ixa,jya) = 0.0
											!print *,'cdf_interpolated, analysis_cdf(1) = ',&
											!	cdf_interpolated, analysis_cdf(ixa,jya,1)
                                 		    GOTO 3000
                              		  	END IF
                       
                              		    DO icdf = 1,npct-1
                                 			IF (cdf_interpolated .ge. analysis_cdf(ixa,jya,icdf) .and. &
                                    		cdf_interpolated .lt. analysis_cdf(ixa,jya,icdf+1)) THEN
                                    			weight2 = (cdf_interpolated-analysis_cdf(ixa,jya,icdf)) / &
                                            	(analysis_cdf(ixa,jya,icdf+1)-analysis_cdf(ixa,jya,icdf))
                                    			!ktr = ktr+1
                                    			forecast_x9(ktr,ixa,jya) = &
													thresh(icdf)*(1.-weight2) + thresh(icdf+1)*weight2
												!print *,'weight2, analysis_cdf(icdf), analysis_cdf(icdf+1) = ',&
												!	weight2, analysis_cdf(ixa,jya,icdf), &
												!	analysis_cdf(ixa,jya,icdf+1)
												!print *,'forecast_x9(ktr,ixa,jya) = ', forecast_x9(ktr,ixa,jya)
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
                                 		   	slope = SUM((a95_to_a99(2:5)-a95_to_a99(1)) * &
												(f95_to_f99(2:5)-f95_to_f99(1))) / &
                                       		SUM((f95_to_f99(2:5)-f95_to_f99(1))**2)
                                 		    forecast_x9(ktr,ixa,jya) = a95_to_a99(1) + &
												(forecast(ixn,jyn)-f95_to_f99(1))*slope
                              			ELSE
                                 		    forecast_x9(ktr,ixa,jya) = forecast(ixn,jyn)
                              		  	ENDIF 

                              		  	GOTO 3000

                           		 	ENDIF ! cdf_interpoloated < 0.95
                        		ENDIF  ! forecast >= thresh(ipct), < thresh(ipct)+1
                     	   	END DO ! ipct

                  	 	ELSE ! conusmask, forecast>0
							forecast_x9(ktr,ixa,jya) = forecast(ixn,jyn)
						ENDIF
						3000 CONTINUE
						
						!PRINT *,'ktr,jyn,ixn, forecast_x9, conusmask(ixn,jyn) = ',&
						!	ktr,jyn,ixn,forecast_x9(ktr,ixa,jya), conusmask(ixn,jyn)
                 	END IF ! jyn>1, jyn < nya, etc.
            	END DO ! ixn
         	END DO ! jyn
      	END IF ! conusmask
   	END DO ! ixa
END DO ! jya

!

!print *,'  forecast_x9(:, 3*nxa/4, nya/2) = ',forecast_x9(:,3*nxa/4, nya/2)

RETURN
END SUBROUTINE cdf_correct_x9_local
