!  f2py -c -m verify_relia_bss_boots verify_relia_bss_boots.f90 ran3.f sort.f

SUBROUTINE verify_relia_bss_boots(cleade, cthresh, nclasses, &
    rthresh, cmodelcombo, date_list_anal, apcp_anal_t, nxa, nya, ndates, &
    relia_raw, relia_raw_05, relia_raw_95, frequse_raw, bss_raw, &
    relia_raw_NCEP, relia_raw_05_NCEP, relia_raw_95_NCEP, frequse_raw_NCEP, bss_raw_NCEP, &
    relia_raw_CMC, relia_raw_05_CMC, relia_raw_95_CMC, frequse_raw_CMC, bss_raw_CMC, &
    relia_raw_ECMWF, relia_raw_05_ECMWF, relia_raw_95_ECMWF, frequse_raw_ECMWF, bss_raw_ECMWF, &
    relia_cdf, relia_cdf_05, relia_cdf_95, frequse_cdf, bss_cdf, &
    relia_cdf_NCEP, relia_cdf_05_NCEP, relia_cdf_95_NCEP, frequse_cdf_NCEP, bss_cdf_NCEP, &
    relia_cdf_CMC, relia_cdf_05_CMC, relia_cdf_95_CMC, frequse_cdf_CMC, bss_cdf_CMC, &
    relia_cdf_ECMWF, relia_cdf_05_ECMWF, relia_cdf_95_ECMWF, frequse_cdf_ECMWF, bss_cdf_ECMWF)

	! purpose:  generate reliability information, frequency of usage, and Brier
	!  skill scores for post-processed and raw forecasts
	
PARAMETER (nresa = 1000)
CHARACTER*3, INTENT(IN) :: cleade, cmodelcombo
CHARACTER*(*), INTENT(IN) :: cthresh
INTEGER, INTENT(IN) :: nclasses, nxa, nya, ndates
REAL, INTENT(IN) :: rthresh ! precip threshold amount
REAL, DIMENSION(nxa,nya,ndates), INTENT(IN) :: apcp_anal_t
CHARACTER*10, DIMENSION(ndates), INTENT(IN) :: date_list_anal

REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_raw, relia_raw_05, relia_raw_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_raw_NCEP, relia_raw_05_NCEP, relia_raw_95_NCEP
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_raw_CMC, relia_raw_05_CMC, relia_raw_95_CMC
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_raw_ECMWF, relia_raw_05_ECMWF, relia_raw_95_ECMWF

REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_cdf, relia_cdf_05, relia_cdf_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_cdf_NCEP, relia_cdf_05_NCEP, relia_cdf_95_NCEP
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_cdf_CMC, relia_cdf_05_CMC, relia_cdf_95_CMC
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_cdf_ECMWF, relia_cdf_05_ECMWF, relia_cdf_95_ECMWF
	
REAL, DIMENSION(nclasses), INTENT(OUT) :: frequse_raw, frequse_raw_NCEP, &
    frequse_raw_CMC, frequse_raw_ECMWF, frequse_cdf, frequse_cdf_NCEP, &
    frequse_cdf_CMC, frequse_cdf_ECMWF
    
REAL, INTENT(OUT) :: bss_raw, bss_raw_NCEP, bss_raw_CMC, bss_raw_ECMWF
REAL, INTENT(OUT) :: bss_cdf, bss_cdf_NCEP, bss_cdf_CMC, bss_cdf_ECMWF

! f2py intent(in) nxa, nya, ndates, cleade, cthresh, nclasses, 
! f2py intent(in) rthresh, cmodelcombo, apcp_anal_t, date_list_anal
! f2py depend(ndates) date_list_anal
! f2py depend(nxa,nya,ndates) apcp_anal_t
! f2py intent(out) relia_raw, relia_raw_05, relia_raw_95
! f2py intent(out) relia_raw_NCEP, relia_raw_05_NCEP, relia_raw_95_NCEP
! f2py intent(out) relia_raw_CMC, relia_raw_05_CMC, relia_raw_95_CMC
! f2py intent(out) relia_raw_ECMWF, relia_raw_05_ECMWF, relia_raw_95_ECMWF
! f2py intent(out) relia_cdf, relia_cdf_05, relia_cdf_95
! f2py intent(out) relia_cdf_NCEP, relia_cdf_05_NCEP, relia_cdf_95_NCEP
! f2py intent(out) relia_cdf_CMC, relia_cdf_05_CMC, relia_cdf_95_CMC
! f2py intent(out) relia_cdf_ECMWF, relia_cdf_05_ECMWF, relia_cdf_95_ECMWF
! f2py intent(out)  frequse_raw, frequse_raw_NCEP, frequse_raw_CMC, frequse_raw_ECMWF
! f2py intent(out)  frequse_cdf, frequse_cdf_NCEP, frequse_cdf_CMC, frequse_cdf_ECMWF
! f2py intent(out)   bss_raw, bss_raw_NCEP, bss_raw_CMC, bss_raw_ECMWF
! f2py intent(out)   bss_cdf, bss_cdf_NCEP, bss_cdf_CMC, bss_cdf_ECMWF
! f2py depend(nclasses) relia_raw, relia_raw_05, relia_raw_95
! f2py depend(nclasses) relia_raw_NCEP, relia_raw_05_NCEP, relia_raw_95_NCEP
! f2py depend(nclasses) relia_raw_CMC, relia_raw_05_CMC, relia_raw_95_CMC
! f2py depend(nclasses) relia_raw_ECMWF, relia_raw_05_ECMWF, relia_raw_95_ECMWF
! f2py depend(nclasses) relia_cdf, relia_cdf_05, relia_cdf_95
! f2py depend(nclasses) relia_cdf_NCEP, relia_cdf_05_NCEP, relia_cdf_95_NCEP
! f2py depend(nclasses) relia_cdf_CMC, relia_cdf_05_CMC, relia_cdf_95_CMC
! f2py depend(nclasses) relia_cdf_ECMWF, relia_cdf_05_ECMWF, relia_cdf_95_ECMWF
! f2py depend(nclasses)  frequse_raw, frequse_raw_NCEP, frequse_raw_CMC, frequse_raw_ECMWF
! f2py depend(nclasses)  frequse_cdf, frequse_cdf_NCEP, frequse_cdf_CMC, frequse_cdf_ECMWF

! --- now local variables

INTEGER*2, DIMENSION(nxa,nya) :: conusmask

REAL, DIMENSION(nxa,nya) :: climo_prob
REAL, DIMENSION(nxa,nya) :: rlonsa
REAL, DIMENSION(nxa,nya) :: rlatsa
REAL, DIMENSION(nxa,nya) :: prob_forecast
REAL, DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, DIMENSION(nxa,nya) :: prob_forecast_raw_CMC
REAL, DIMENSION(nxa,nya) :: prob_forecast_raw_NCEP
REAL, DIMENSION(nxa,nya) :: prob_forecast_raw_ECMWF
REAL, DIMENSION(nxa,nya) :: prob_forecast_cdf
REAL, DIMENSION(nxa,nya) :: prob_forecast_cdf_CMC
REAL, DIMENSION(nxa,nya) :: prob_forecast_cdf_NCEP
REAL, DIMENSION(nxa,nya) :: prob_forecast_cdf_ECMWF




REAL*8, DIMENSION(0:nclasses-1,2) :: contab
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_sum

REAL*8, DIMENSION(0:nclasses-1,2) :: contab_raw
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_raw_NCEP
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_raw_CMC
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_raw_ECMWF

REAL*8, DIMENSION(0:nclasses-1,2) :: contab_cdf
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_cdf_NCEP
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_cdf_CMC
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_cdf_ECMWF

REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_daily

REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_raw_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_raw_daily_NCEP
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_raw_daily_CMC
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_raw_daily_ECMWF

REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_cdf_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_cdf_daily_NCEP
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_cdf_daily_CMC
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_cdf_daily_ECMWF

REAL*8, DIMENSION(0:nclasses-1) :: relia, relia_05, relia_95
REAL*8, DIMENSION(nresa, 0:nclasses-1) :: relia_resa
REAL*8, DIMENSION(0:nclasses-1) :: frequse

REAL*8 :: bs
REAL*8 :: bs_climo

REAL*8 :: bs_raw
REAL*8 :: bs_raw_NCEP
REAL*8 :: bs_raw_CMC
REAL*8 :: bs_raw_ECMWF

REAL*8 :: bs_cdf
REAL*8 :: bs_cdf_NCEP
REAL*8 :: bs_cdf_CMC
REAL*8 :: bs_cdf_ECMWF

REAL, DIMENSION(nresa) :: rsamps

CHARACTER*120 infile
CHARACTER*10 cyyyymmddhh
CHARACTER*3 clead_use

LOGICAL iex


contab_raw = 0.
contab_raw_NCEP = 0.
contab_raw_CMC = 0.
contab_raw_ECMWF = 0.

contab_cdf = 0.
contab_cdf_NCEP = 0.
contab_cdf_CMC = 0.
contab_cdf_ECMWF = 0.

contab_raw_daily = 0.
contab_raw_daily_NCEP = 0.
contab_raw_daily_CMC = 0.
contab_raw_daily_ECMWF = 0.

contab_cdf_daily = 0.
contab_cdf_daily_NCEP = 0.
contab_cdf_daily_CMC = 0.
contab_cdf_daily_ECMWF = 0.

bs_climo = 0.

bs_raw = 0.
bs_raw_NCEP = 0.
bs_raw_CMC = 0.
bs_raw_ECMWF = 0.

bs_cdf = 0.
bs_cdf_NCEP = 0.
bs_cdf_CMC = 0.
bs_cdf_ECMWF = 0.

popthresh = 0.254  ! .4 mm is considered the precipitation threshold for POP
pid180 = 3.1415926/180.
clead_use = cleade
IF (clead_use .eq. '012') clead_use = '12'
IF (clead_use .eq. '024') clead_use = '24'
IF (clead_use .eq. '036') clead_use = '36'
IF (clead_use .eq. '048') clead_use = '48'
IF (clead_use .eq. '060') clead_use = '60'
IF (clead_use .eq. '072') clead_use = '72'
IF (clead_use .eq. '084') clead_use = '84'
IF (clead_use .eq. '096') clead_use = '96'

! ---- loop thru days and verify

!READ (cyyyymmddhh_begin,'(i10)') iyyyymmddhh_begin
DO idate = 1, ndates

    ! ---- read in the forecast for this particular initial date/time and forecast lead and resolution

    !CALL doy(iyyyymmddhh,iyear,imo,iday,ihour,idoy)
    !ioffset = 24*(idate-1)
    !CALL updat(iyyyymmddhh,ioffset,iyyyymmddhh_toread)
    !WRITE (cyyyymmddhh,'(i10)') iyyyymmddhh_toread
    cyyyymmddhh = date_list_anal(idate)
    
    IF (TRIM(cmodelcombo) .eq. 'E') THEN
        infile = '/Users/thamill/precip/ecmwf_data/ECMWF_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ELSEIF (TRIM(cmodelcombo) .eq. 'N') THEN
        infile = '/Users/thamill/precip/ecmwf_data/NCEP_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ELSEIF (TRIM(cmodelcombo) .eq. 'C') THEN
        infile = '/Users/thamill/precip/ecmwf_data/CMC_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ELSEIF (TRIM(cmodelcombo) .eq. 'EC') THEN
        infile = '/Users/thamill/precip/ecmwf_data/ECMWF_CMC_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ELSEIF (TRIM(cmodelcombo) .eq. 'EN') THEN
        infile = '/Users/thamill/precip/ecmwf_data/ECMWF_NCEP_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ELSEIF (TRIM(cmodelcombo) .eq. 'NC') THEN
        infile = '/Users/thamill/precip/ecmwf_data/NCEP_CMC_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ELSE
        infile = '/Users/thamill/precip/ecmwf_data/ECMWF_NCEP_CMC_'//&
            TRIM(clead_use)//'h_IC'//cyyyymmddhh//&
            '_thresh'//TRIM(cthresh)//'.dat'
    ENDIF

    PRINT *, TRIM(infile)
    INQUIRE (file=infile,exist=iex)
    IF (iex) THEN
       
        OPEN (unit=43, file=infile, status='old', form='unformatted')
        READ (43) nxain,nyain
        READ (43) prob_forecast_raw
        READ (43) prob_forecast_raw_CMC
        READ (43) prob_forecast_raw_NCEP
        READ (43) prob_forecast_raw_ECMWF
        READ (43) prob_forecast_cdf
        READ (43) prob_forecast_cdf_CMC
        READ (43) prob_forecast_cdf_NCEP
        READ (43) prob_forecast_cdf_ECMWF
        READ (43) climo_prob
        READ (43) rlonsa
        READ (43) rlatsa
        READ (43) conusmask
        CLOSE(43)    
	    
	    DO itype = 1, 8
		    contab = 0.
		    bs = 0.
		    IF (itype .eq. 1) THEN ! MME
			    prob_forecast(:,:) = prob_forecast_raw
		    ELSE IF (itype .eq. 2) THEN
                prob_forecast(:,:) = prob_forecast_raw_NCEP
		    ELSE IF (itype .eq. 3) THEN
                prob_forecast(:,:) = prob_forecast_raw_CMC
	   	    ELSE IF (itype .eq. 4) THEN
                prob_forecast(:,:) = prob_forecast_raw_ECMWF
		    ELSE IF (itype .eq. 5) THEN
                prob_forecast(:,:) = prob_forecast_cdf
    		ELSE IF (itype .eq. 6) THEN
                prob_forecast(:,:) = prob_forecast_cdf_NCEP                
        	ELSE IF (itype .eq. 7) THEN
                prob_forecast(:,:) = prob_forecast_cdf_CMC     
            ELSE  !IF (itype .eq. 8) THEN
                prob_forecast(:,:) = prob_forecast_cdf_ECMWF     
		    ENDIF  
	   
	   	    ! ---- now let's verify probability forecast

		    DO i = 1, nxa
          	    DO j = 1, nya
             	    IF (conusmask(i,j) .gt. 0 .and. prob_forecast(i,j) .GE. 0.0 .and. &
             	    prob_forecast(i,j) .le. 1.0 .and. apcp_anal_t(i,j,idate) .GE. 0.0) THEN
            	 	    cfac = cos(rlatsa(i,j)*pid180)  ! cos of latitude to acct for unequal grid box size
            	 	    pclimo = climo_prob(i,j)
            	 	    p      = prob_forecast(i,j)
            	 	    ipcat  = nint(p*20)
            	 	    v      = apcp_anal_t(i,j,idate)
				 
            	 	    IF (v .GE. rthresh) THEN
               	    	    contab(ipcat,2) = contab(ipcat,2) + cfac
                    	    bs = bs + cfac*(1.-p)**2
                    	    IF (itype .eq. 1) bs_climo = bs_climo + cfac*(1.-pclimo)**2
            	 	    ELSE
                    	    contab(ipcat,1) = contab(ipcat,1) + cfac
                    	    bs = bs + cfac * p**2
                    	    IF (itype .eq. 1) bs_climo = bs_climo + cfac*pclimo**2
                 	    ENDIF
             	    ENDIF ! conusmask
      	  	    END DO   ! j = 1, nya
   	   	    END DO      ! i = 1, nxa
	   
	   	    IF (itype .eq. 1) THEN
		   	    contab_raw = contab_raw + contab
			    contab_raw_daily(idate,:,:) = contab(:,:)
		   	    bs_raw = bs_raw + bs
	   	    ELSE IF (itype .eq. 2) THEN
		   	    contab_raw_NCEP = contab_raw_NCEP + contab
			    contab_raw_daily_NCEP(idate,:,:) = contab(:,:)
		   	    bs_raw_NCEP = bs_raw_NCEP + bs
	        ELSE IF (itype .eq. 3) THEN
		   	    contab_raw_CMC = contab_raw_CMC + contab
			    contab_raw_daily_CMC(idate,:,:) = contab(:,:)
		   	    bs_raw_CMC = bs_raw_CMC + bs
	        ELSE IF (itype .eq. 4) THEN
		   	    contab_raw_ECMWF = contab_raw_ECMWF + contab
			    contab_raw_daily_ECMWF(idate,:,:) = contab(:,:)
		   	    bs_raw_ECMWF = bs_raw_ECMWF + bs
	        ELSE IF (itype .eq. 5) THEN
		        contab_cdf = contab_cdf + contab
			    contab_cdf_daily(idate,:,:) = contab(:,:)
		        bs_cdf = bs_cdf + bs
    	    ELSE IF (itype .eq. 6) THEN
    		    contab_cdf_NCEP = contab_cdf_NCEP + contab
    			contab_cdf_daily_NCEP(idate,:,:) = contab(:,:)
    		    bs_cdf_NCEP = bs_cdf_NCEP + bs
        	ELSE IF (itype .eq. 7) THEN
        		contab_cdf_CMC = contab_cdf_CMC + contab
        		contab_cdf_daily_CMC(idate,:,:) = contab(:,:)
        		bs_cdf_CMC = bs_cdf_CMC + bs
            ELSE IF (itype .eq. 8) THEN
            	contab_cdf_ECMWF = contab_cdf_ECMWF + contab
            	contab_cdf_daily_ECMWF(idate,:,:) = contab(:,:)
            	bs_cdf_ECMWF = bs_cdf_ECMWF + bs
	        ENDIF 
	   
   	    END DO  ! itype
   ENDIF !(iex)
END DO ! idate
 

DO itype = 1, 8
	
	IF (itype .eq. 1) THEN
		PRINT *,'Raw MME'
		contab = contab_raw
		contab_daily = contab_raw_daily
		bs = bs_raw
	ELSE IF (itype .eq. 2) THEN
		PRINT *, 'raw NCEP'
		contab = contab_raw_NCEP
		contab_daily = contab_raw_daily_NCEP
		bs = bs_raw_NCEP
	ELSE IF (itype .eq. 3) THEN
		print *, 'raw CMC'
		contab = contab_raw_CMC
		contab_daily = contab_raw_daily_CMC
		bs = bs_raw_CMC
	ELSE IF (itype .eq. 4) THEN
		PRINT *, 'raw ECMWF'
		contab = contab_raw_ECMWF
		contab_daily = contab_raw_daily_ECMWF
		bs = bs_raw_ECMWF
	ELSE IF (itype .eq. 5) THEN
		PRINT *, 'calibrated MME'
		contab = contab_cdf
		contab_daily = contab_cdf_daily
	    bs = bs_cdf
	ELSE IF (itype .eq. 6) THEN
		PRINT *, 'calibrated NCEP'
		contab = contab_cdf_NCEP
		contab_daily = contab_cdf_daily_NCEP
		bs = bs_cdf_NCEP
	ELSE IF (itype .eq. 7) THEN
		print *, 'calibrated CMC'
		contab = contab_cdf_CMC
		contab_daily = contab_cdf_daily_CMC
		bs = bs_cdf_CMC
	ELSE IF (itype .eq. 8) THEN
		PRINT *, 'calibrated ECMWF'
		contab = contab_cdf_ECMWF
		contab_daily = contab_cdf_daily_ECMWF
		bs = bs_cdf_ECMWF      
	ENDIF
	
	! ---- with tallied contingency tables, now set mean reliability and frequency of use for mean

	print *,'bs, bs_climo = ',bs, bs_climo
	ctot = SUM(contab)
	bss = 1. - bs / bs_climo
	relia(:) = -99.9999
	PRINT *,'  bss = ',bss
	PRINT *,'  p   reliability   freq of usage'
	DO icat = 0,20
		frequse(icat) = (contab(icat,2)  + contab(icat,1)) / ctot
		IF ((contab(icat,1) + contab(icat,2)) .gt. 0) THEN
			relia(icat) = contab(icat,2) / (contab(icat,1) + contab(icat,2))
			PRINT 203,float(icat*5),relia(icat)*100.,frequse(icat)
		ELSE
			PRINT 203,float(icat*5),relia(icat),frequse(icat)
		ENDIF
  	    203 format(f5.1,3x,2(f8.3,3x))
	END DO  !icat
	
	! ---- perform a resampling to generate the confidence intervals for reliability
	!      sampling the days with replacement
	
	idum = -12345
	relia_resa(:,:) = -99.9999
	DO iresa = 1, nresa
		contab_sum = 0.
		DO idate = 1, ndates
			cran = ran3(idum)
			idate2 = MIN(1+NINT(cran*ndates),ndates)
			contab_sum(:,:) = contab_sum(:,:) + contab_daily(idate2,:,:)
		END DO
		DO icat = 0,20
			IF ((contab_sum(icat,1) + contab_sum(icat,2)) .gt. 0) THEN
				relia_resa(iresa,icat) = contab_sum(icat,2) / (contab_sum(icat,1) + contab_sum(icat,2))
			ENDIF
		END DO ! icat
	END DO  !iresa
	
    ! ---- now find the 5th and 95th percentiles of the resampled distribution, accounting
	!      for the possibility that at high probabilities there may be no samples in some 
	!      cases.

	DO i = 0, nclasses-1
		rsamps(:) = relia_resa(:,i)
		CALL sort(nresa, rsamps)
		ibegin = 1
		DO iresa = 1, nresa
			IF (rsamps(iresa) .gt. -99.) THEN
				ibegin = iresa
				goto 3212
			ENDIF
		END DO ! iresa
3212	isamp = MIN(ibegin + NINT(0.05* (nresa-ibegin+1) ),nresa)
		relia_05(i) = rsamps(isamp)
		isamp = MIN(ibegin + NINT(0.95* (nresa-ibegin+1) ),nresa)
		relia_95(i) = rsamps(isamp)
	END DO ! i 
	
	! ---- copy to output arrays
		
	IF (itype .eq. 1) THEN
		relia_raw = relia
		relia_raw_05 = relia_05
		relia_raw_95 = relia_95
		frequse_raw = frequse
		bss_raw = bss
	ELSE IF (itype .eq. 2) THEN
		relia_raw_NCEP = relia
		relia_raw_05_NCEP = relia_05
		relia_raw_95_NCEP = relia_95
		frequse_raw_NCEP = frequse
		bss_raw_NCEP = bss
	ELSE IF (itype .eq. 3) THEN
		relia_raw_CMC = relia
		relia_raw_05_CMC = relia_05
		relia_raw_95_CMC = relia_95
		frequse_raw_CMC = frequse
		bss_raw_CMC = bss
	ELSE IF (itype .eq. 4) THEN
		relia_raw_ECMWF = relia
		relia_raw_05_ECMWF = relia_05
		relia_raw_95_ECMWF = relia_95
		frequse_raw_ECMWF = frequse
	    bss_raw_ECMWF = bss
	ELSE IF (itype .eq. 5) THEN
		relia_cdf = relia
		relia_cdf_05 = relia_05
		relia_cdf_95 = relia_95
		frequse_cdf = frequse
		bss_cdf = bss
	ELSE IF (itype .eq. 6) THEN
		relia_cdf_NCEP = relia
		relia_cdf_05_NCEP = relia_05
		relia_cdf_95_NCEP = relia_95
		frequse_cdf_NCEP = frequse
		bss_cdf_NCEP = bss
	ELSE IF (itype .eq. 7) THEN
		relia_cdf_CMC = relia
		relia_cdf_05_CMC= relia_05
		relia_cdf_95_CMC = relia_95
		frequse_cdf_CMC = frequse
		bss_cdf_CMC = bss
	ELSE IF (itype .eq. 8) THEN
		relia_cdf_ECMWF = relia
		relia_cdf_05_ECMWF = relia_05
		relia_cdf_95_ECMWF = relia_95
		frequse_cdf_ECMWF = frequse
		bss_cdf_ECMWF = bss
	ENDIF
	
END DO  ! itype
	
RETURN
END SUBROUTINE verify_relia_bss_boots

