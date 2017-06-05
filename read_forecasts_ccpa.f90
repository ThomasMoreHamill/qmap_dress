SUBROUTINE read_forecasts_ccpa (nxf, nyf, nxa, nya, nens_ecmwf, nens_cmc, nens_ncep, iyyyymmddhh, &
     cresolution, pfcst_infile, ecmwf_deterministic_ccpa, cmc_control_ccpa, ncep_control_ccpa, &
     ecmwf_ensemble_ccpa, cmc_ensemble_ccpa, ncep_ensemble_ccpa, &
     iec_d_ccpa, icm_c_ccpa, inc_c_ccpa, &
     iec_e_ccpa, inc_e_ccpa, icm_e_ccpa)

! ---- purpose: read in ensemble/deterministic/control forecasts for the day of interest, but
!      here on the 1/8-degree CCPA grid.  In this way one can compare the interpolated 
!      ensemble predictions from the direct ensemble to those with additional statistical 
!      downscaling.  

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, nens_ecmwf, nens_cmc, nens_ncep, iyyyymmddhh
CHARACTER*2, INTENT(IN) :: cresolution
CHARACTER*(*), INTENT(IN) :: pfcst_infile

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: ecmwf_deterministic_ccpa
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: cmc_control_ccpa
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: ncep_control_ccpa
!REAL, INTENT(OUT), DIMENSION(nxa,nya) :: ukmo_control_ccpa
REAL, INTENT(OUT), DIMENSION(nxa,nya,nens_ecmwf) :: ecmwf_ensemble_ccpa
REAL, INTENT(OUT), DIMENSION(nxa,nya,nens_cmc)  :: cmc_ensemble_ccpa
REAL, INTENT(OUT), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa
!REAL, INTENT(OUT), DIMENSION(nxa,nya,nens_ukmo) :: ukmo_ensemble_ccpa

INTEGER, INTENT(OUT):: iec_d_ccpa, icm_c_ccpa, inc_c_ccpa, &
     iec_e_ccpa, inc_e_ccpa, icm_e_ccpa
     ! these are flags, 0 if data not available, 1 if it is.

REAL, DIMENSION(nxa,nya) :: rlonsa
REAL, DIMENSION(nxa,nya) :: rlatsa

INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_anal
REAL, DIMENSION(nxf,nyf) :: pfcst
REAL, DIMENSION(nxa,nya) :: pfcst_ccpa
CHARACTER*20 cfield
CHARACTER*4 cname

! ---- open the file, use the number of times in the file later to allocate a date array

netid = 0
PRINT *,'netid, reading from ',netid, TRIM(pfcst_infile)
CALL check (nf90_open(pfcst_infile,NF90_NOWRITE,netid))
   
! ---- read in the list of dates/times in yyyymmddhh format associated with 
!      each time index stored in the netcdf file; also lats and lons 

cfield ='time'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_inquire_dimension(netid,ivar,cname,ntimes))
print *,'number of forecast dates to read = ',ntimes
ALLOCATE (iyyyymmddhh_anal(ntimes))

cfield ='yyyymmddhh_date_anal'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,iyyyymmddhh_anal,&
  start=(/1/),count=(/ntimes/)))

cfield = 'lons_fcst_ccpa'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
print *,'ivar = ',ivar
CALL check (nf90_get_var(netid,ivar,rlonsa,&
   start=(/1,1/),count=(/nxa,nya/)))

cfield ='lats_fcst_ccpa'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,rlatsa,&
   start=(/1,1/),count=(/nxa,nya/)))

! ---- loop thru the samples, find the day that matches the forecast

iec_d_ccpa = 0  ! these are flags that, if in your own application, you will not reliably have
icm_c_ccpa = 0  ! all the data sources available, then set the flag to zero
inc_c_ccpa = 0  ! so as to tell the downscaling algorithm to skip this data source this time.
!iuk_c_ccpa = 0
iec_e_ccpa = 0
inc_e_ccpa = 0
icm_e_ccpa = 0
!iuk_e_ccpa = 0

ktrday = 1
DO itime = 1, ntimes
   IF (iyyyymmddhh_anal(itime) .eq. iyyyymmddhh) THEN
      

      pfcst_ccpa = -99.99
      cfield = 'fcst_ecmwf_det_ccpa'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst_ccpa,&
         start=(/1,1,itime/),count=(/nxa,nya,1/)))
      CALL make_tiny_negatives_zero(nxa,nya,pfcst_ccpa)
      ecmwf_deterministic_ccpa(:,:) = pfcst_ccpa(:,:)
      iec_d_ccpa = 1

      pfcst_ccpa = -99.99
      cfield = 'fcst_ncep_cntl_ccpa'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst_ccpa,&
         start=(/1,1,itime/),count=(/nxa,nya,1/)))
      CALL make_tiny_negatives_zero(nxa,nya,pfcst_ccpa)
      ncep_control_ccpa(:,:) = pfcst_ccpa(:,:)
      inc_c_ccpa = 1

      pfcst_ccpa = -99.99
      cfield = 'fcst_cmc_cntl_ccpa'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst_ccpa,&
         start=(/1,1,itime/),count=(/nxa,nya,1/)))
      CALL make_tiny_negatives_zero(nxa,nya,pfcst_ccpa)
      cmc_control_ccpa(:,:) = pfcst_ccpa(:,:)
      icm_c_ccpa = 1

!      pfcst_ccpa = -99.99
!      cfield = 'fcst_ukmo_cntl_ccpa'
!      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
!      CALL check (nf90_get_var(netid,ivar,pfcst_ccpa,&
!         start=(/1,1,itime/),count=(/nxa,nya,1/)))
!      CALL make_tiny_negatives_zero(nxa,nya,pfcst_ccpa)
!      ukmo_control_ccpa(:,:) = pfcst_ccpa(:,:)
!      iuk_c_ccpa = 1

      ecmwf_ensemble_ccpa = -99.99
      cfield = 'fcst_ecmwf_ens_ccpa'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar, ecmwf_ensemble_ccpa,&
         start=(/1,1,1,itime/), count=(/nxa,nya,nens_ecmwf,1/)))
      DO imem = 1, nens_ecmwf
         CALL make_tiny_negatives_zero(nxa,nya,ecmwf_ensemble_ccpa(1,1,imem))
      END DO
      iec_e_ccpa = 1

      ncep_ensemble_ccpa = -99.99
      cfield = 'fcst_ncep_ens_ccpa'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar, ncep_ensemble_ccpa,&
         start=(/1,1,1,itime/), count=(/nxa,nya,nens_ncep,1/)))
      DO imem = 1, nens_ncep
         CALL make_tiny_negatives_zero(nxa,nya,ncep_ensemble_ccpa(1,1,imem))
      END DO
      inc_e_ccpa = 1

      cmc_ensemble_ccpa = -99.99
      cfield = 'fcst_cmc_ens_ccpa'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar, cmc_ensemble_ccpa,&
         start=(/1,1,1,itime/), count=(/nxa,nya,nens_cmc,1/)))
      DO imem = 1, nens_cmc
         CALL make_tiny_negatives_zero(nxa,nya,cmc_ensemble_ccpa(1,1,imem))
      END DO
      icm_e_ccpa = 1

!      ukmo_ensemble_ccpa = -99.99
!      cfield = 'fcst_ukmo_ens_ccpa'
!      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
!      CALL check (nf90_get_var(netid, ivar, ukmo_ensemble_ccpa,&
!         start=(/1,1,1,itime/), count=(/nxa,nya,nens_ukmo,1/)))
!      DO imem = 1, nens_ukmo
!         CALL make_tiny_negatives_zero(nxa,nya,ukmo_ensemble_ccpa(1,1,imem))
!      END DO
!      iuk_e_ccpa = 1

   ENDIF
END DO

! --- close netcdf file.

CALL check(nf90_close(netid))

DEALLOCATE (iyyyymmddhh_anal)

RETURN
END SUBROUTINE read_forecasts_ccpa
