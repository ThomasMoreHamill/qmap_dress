SUBROUTINE grand_ensemble_mean(nxa, nya, nens_cmc, nens_ncep, icm_c_ccpa, inc_c_ccpa,&
           inc_e_ccpa, icm_e_ccpa, ncep_ensemble_ccpa, ncep_control_ccpa, &
           cmc_ensemble_ccpa, cmc_control_ccpa, gem_today)
		   
		   
! purpose:  generate an ensemble mean using both CMC and NCEP ensemble forecasts
!   as well as CMC and NCEP deterministic forecasts, which, since they are
!   slightly more accurate, are given more weight.
!
! coded by Tom Hamill, Apr 2016, tom.hamill@noaa.gov
		   
INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions
INTEGER, INTENT(IN) :: nens_cmc, nens_ncep ! number of ensemble members
INTEGER, INTENT(IN) :: icm_c_ccpa, inc_c_ccpa ! flags for existence of control
INTEGER, INTENT(IN) :: inc_e_ccpa, icm_e_ccpa ! flags for existence of ensembles

REAL, INTENT(IN), DIMENSION(nxa, nya, nens_cmc) :: cmc_ensemble_ccpa 
REAL, INTENT(IN), DIMENSION(nxa, nya, nens_ncep) :: ncep_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa, nya) :: cmc_control_ccpa, ncep_control_ccpa

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: gem_today

INTEGER, DIMENSION(nxa,nya) :: counter

! --- initialize

counter(:,:) = 0
gem_today(:,:) = 0.0

! --- loop over grid points

DO ixa = 1, nxa
	DO jya = 1, nya
	
		! ---- add contribution of cmc ensemble
		
		IF (icm_e_ccpa .gt. 0) THEN
			DO imem = 1, nens_cmc
				IF (cmc_ensemble_ccpa(ixa,jya,imem) .ge. 0.0) THEN
					gem_today(ixa,jya) = gem_today(ixa,jya) + &
						cmc_ensemble_ccpa(ixa,jya,imem)
					counter(ixa,jya) = counter(ixa,jya) + 1
				ENDIF
			END DO
		ENDIF
		
		! ---- add contribution of ncep ensemble
		
		IF (inc_e_ccpa .gt. 0) THEN
			DO imem = 1, nens_ncep
				IF (ncep_ensemble_ccpa(ixa,jya,imem) .ge. 0.0) THEN
					gem_today(ixa,jya) = gem_today(ixa,jya) + &
						ncep_ensemble_ccpa(ixa,jya,imem)
					counter(ixa,jya) = counter(ixa,jya) + 1
				ENDIF
			END DO
		ENDIF
		
		! ---- add contribution of cmc control.  Note 4x greater weight.
		
		IF (icm_c_ccpa .gt. 0 .and. cmc_control_ccpa(ixa,jya) .ge. 0.0) THEN
			gem_today(ixa,jya) = gem_today(ixa,jya) + &
				4.*cmc_control_ccpa(ixa,jya)
			counter(ixa,jya) = counter(ixa,jya) + 4
		ENDIF
		
		! ---- add contribution of ncep control.   Note 4x greater weight.
		
		IF (inc_c_ccpa .gt. 0 .and. ncep_control_ccpa(ixa,jya) .ge. 0.0) THEN
			gem_today(ixa,jya) = gem_today(ixa,jya) + &
				4.*ncep_control_ccpa(ixa,jya)
			counter(ixa,jya) = counter(ixa,jya) + 4
		ENDIF
		
		! ---- finally, calculate the grand mean.
		
		gem_today(ixa,jya) = gem_today(ixa,jya) / REAL(counter(ixa,jya))
		
	END DO ! jya = 1, nya
END DO ! ixa = 1, nxa

RETURN
END SUBROUTINE grand_ensemble_mean

		
