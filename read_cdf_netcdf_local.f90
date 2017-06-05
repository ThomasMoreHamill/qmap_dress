SUBROUTINE read_cdf_netcdf_local(nxa, nya, npct, nens_cmc, iyyyymmddhh, &
    cleade, precip_anal_cdf, ecmwf_ensemble_cdf, ncep_ensemble_cdf, &
    cmc_ensemble_cdf, thresh)
    
USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, npct, nens_cmc, iyyyymmddhh
CHARACTER*(*), INTENT(IN) :: cleade
REAL, INTENT(OUT), DIMENSION(nxa,nya,npct) :: precip_anal_cdf
REAL, INTENT(OUT), DIMENSION(nxa,nya,npct) :: ecmwf_ensemble_cdf
REAL, INTENT(OUT), DIMENSION(nxa,nya,npct) :: ncep_ensemble_cdf
REAL, INTENT(OUT), DIMENSION(nxa,nya,nens_cmc, npct) :: cmc_ensemble_cdf
REAL, INTENT(OUT), DIMENSION(npct) :: thresh

CHARACTER*3 clead
CHARACTER*256 infilename_ecmwf, infilename_ncep, infilename_cmc
CHARACTER*10 cyyyymmddhh_start, cyyyymmddhh_end
CHARACTER*20 cfield

!determine which CDF file to use based on the current date

clead = cleade
IF (cleade .eq. '12' .or. cleade .eq. '012') clead = '12'
IF (cleade .eq. '24' .or. cleade .eq. '024') clead = '24'
IF (cleade .eq. '36' .or. cleade .eq. '036') clead = '36'
IF (cleade .eq. '48' .or. cleade .eq. '048') clead = '48'
IF (cleade .eq. '60' .or. cleade .eq. '060') clead = '60'
IF (cleade .eq. '72' .or. cleade .eq. '072') clead = '72'
IF (cleade .eq. '84' .or. cleade .eq. '084') clead = '84'
IF (cleade .eq. '96' .or. cleade .eq. '096') clead = '96'

! ---- cdf files are named by the initial and end dates of the files.  
!      Specify these.

IF (iyyyymmddhh .lt. 2016041000) THEN
    iyyyymmddhh_reference = 2016033000
ELSE IF (iyyyymmddhh .ge. 2016041000 .and. iyyyymmddhh .lt. 2016041000) THEN
    iyyyymmddhh_reference = 2016041000
ELSE IF (iyyyymmddhh .ge. 2016042000 .and. iyyyymmddhh .lt. 2016043000) THEN
    iyyyymmddhh_reference = 2016042000   
ELSE IF (iyyyymmddhh .ge. 2016043000 .and. iyyyymmddhh .lt. 2016051000) THEN
    iyyyymmddhh_reference = 2016043000 
ELSE IF (iyyyymmddhh .ge. 2016051000 .and. iyyyymmddhh .lt. 2016052000) THEN
    iyyyymmddhh_reference = 2016051000 
ELSE IF (iyyyymmddhh .ge. 2016052000 .and. iyyyymmddhh .lt. 2016053000) THEN
    iyyyymmddhh_reference = 2016052000 
ELSE IF (iyyyymmddhh .ge. 2016053000 .and. iyyyymmddhh .lt. 2016061000) THEN
    iyyyymmddhh_reference = 2016053000 
ELSE IF (iyyyymmddhh .ge. 2016061000 .and. iyyyymmddhh .lt. 2016062000) THEN
    iyyyymmddhh_reference = 2016061000 
ELSE 
    iyyyymmddhh_reference = 2016062000 
ENDIF

READ (cleade,'(i3)') ileade
idaysbefore = - 1 -int(ileade/24)
ihoursbefore = idaysbefore*24
call updat(iyyyymmddhh_reference,ihoursbefore,iyyyymmddhh_end)
call updat(iyyyymmddhh_end,-24*61,iyyyymmddhh_start)
!PRINT *,'expecting file dates from ',iyyyymmddhh_start,' to ',iyyyymmddhh_end
WRITE (cyyyymmddhh_start,'(i10)') iyyyymmddhh_start
WRITE (cyyyymmddhh_end,'(i10)') iyyyymmddhh_end

infilename_ecmwf = &
    '/Users/thamill/precip/ecmwf_data/ECMWF_CDF_flead' // TRIM(clead) // &
    '_' // cyyyymmddhh_start // '_to_' // cyyyymmddhh_end // '.nc'
infilename_ncep = &
    '/Users/thamill/precip/ecmwf_data/NCEP_CDF_flead' // TRIM(clead) // &
    '_' // cyyyymmddhh_start // '_to_' // cyyyymmddhh_end // '.nc'   
infilename_cmc = &
    '/Users/thamill/precip/ecmwf_data/CMC_CDF_flead' // TRIM(clead) // &
    '_' // cyyyymmddhh_start // '_to_' // cyyyymmddhh_end // '.nc'    

! --- read ecmwf cdf and precip analysis cdf

PRINT *, TRIM(infilename_ecmwf)
!PRINT *,'nxa, nya, npct = ', nxa, nya, npct
CALL check (nf90_open(infilename_ecmwf,NF90_NOWRITE,netid))
!PRINT *,'netid, reading from ',netid, TRIM(infilename_ecmwf)

cfield='panal_CDF'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,precip_anal_cdf,&
    start=(/1,1,1/),count=(/nxa,nya,npct/)))
    
cfield='pfcst_CDF'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar, ecmwf_ensemble_cdf,&
    start=(/1,1,1/),count=(/nxa,nya,npct/)))

cfield='thrval'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar, thresh,&
    start=(/1/),count=(/npct/)))

CALL check(nf90_close(netid))

! --- read ncep cdf 

!PRINT *,'netid, reading from ',netid, TRIM(infilename_ncep)
PRINT *, TRIM(infilename_ncep)
CALL check (nf90_open(infilename_ncep,NF90_NOWRITE,netid))

cfield='pfcst_CDF'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar, ncep_ensemble_cdf,&
    start=(/1,1,1/),count=(/nxa,nya,npct/)))

CALL check(nf90_close(netid))

! --- read CMC cdf 

!PRINT *,'netid, reading from ',netid, TRIM(infilename_cmc)
PRINT *, TRIM(infilename_cmc)
CALL check (nf90_open(infilename_cmc,NF90_NOWRITE,netid))

cfield='pfcst_CDF'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
!print *,'netid = ', netid
!print *,'nxa, nya, npct = ', nxa, nya, npct
CALL check(nf90_get_var(netid,ivar, cmc_ensemble_cdf,&
    start=(/1,1,1,1/),count=(/nxa,nya,nens_cmc,npct/)))
!print *,'done reading'

CALL check(nf90_close(netid))

!print *,'ecmwf_ensemble_cdf(nxa/2,nya/4,1:30) = ', ecmwf_ensemble_cdf(nxa/2,nya/4,1:30)
!print *,'ecmwf_ensemble_cdf(nxa/2,nya/2,1:30) = ', ecmwf_ensemble_cdf(nxa/2,nya/2,1:30)
!print *,'ecmwf_ensemble_cdf(nxa/2,3*nya/4,1:30) = ', ecmwf_ensemble_cdf(nxa/2,3*nya/4,1:30)


RETURN
END SUBROUTINE  read_cdf_netcdf_local
