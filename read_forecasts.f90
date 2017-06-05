SUBROUTINE read_forecasts (nxf,nyf, nxa, nya, nens_ecmwf, nens_cmc, nens_ncep, iyyyymmddhh, &
     cresolution, pfcst_infile, ecmwf_deterministic, cmc_control, ncep_control, &
     ecmwf_ensemble, cmc_ensemble, ncep_ensemble, iec_d, icm_c, inc_c, &
     iec_e, inc_e, icm_e)

! ---- purpose: read in ensemble/deterministic/control forecasts for the day of interest.
!      previously these were converted from grib2 files to netcdf files so that they can
!      be loaded up faster and easier.  This was done in the python routine 
!      forecasts_to_polarstereo.py .

USE netcdf

INTEGER, INTENT(IN) :: nxf, nyf, nens_ecmwf, nens_cmc, nens_ncep, iyyyymmddhh
CHARACTER*2, INTENT(IN) :: cresolution
CHARACTER*(*), INTENT(IN) :: pfcst_infile
REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: ecmwf_deterministic
REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: cmc_control
REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: ncep_control
!REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: ukmo_control
REAL, INTENT(OUT), DIMENSION(nxf,nyf,nens_ecmwf) :: ecmwf_ensemble
REAL, INTENT(OUT), DIMENSION(nxf,nyf,nens_cmc)  :: cmc_ensemble
REAL, INTENT(OUT), DIMENSION(nxf,nyf,nens_ncep) :: ncep_ensemble
!REAL, INTENT(OUT), DIMENSION(nxf,nyf,nens_ukmo) :: ukmo_ensemble

INTEGER, INTENT(OUT):: iec_d, icm_c, inc_c, iec_e, inc_e, icm_e
   ! these are flags, 0 if data not available, 1 if it is.

REAL, DIMENSION(nxf,nyf) :: rlonsf
REAL, DIMENSION(nxf,nyf) :: rlatsf

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
!      each time index stored in the netcdf file as well as polar-stereo grid
!      forecast lat/lons 

cfield ='time'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_inquire_dimension(netid,ivar,cname,ntimes))
print *,'number of forecast dates to read = ',ntimes
ALLOCATE (iyyyymmddhh_anal(ntimes))

cfield = 'lons_fcst'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
print *,'ivar = ',ivar
CALL check (nf90_get_var(netid,ivar,rlonsf,&
   start=(/1,1/),count=(/nxf,nyf/)))
print *,'got rlonsf'

cfield ='lats_fcst'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,rlatsf,&
   start=(/1,1/),count=(/nxf,nyf/)))
print *,'check'

cfield ='yyyymmddhh_date_anal'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,iyyyymmddhh_anal,&
  start=(/1/),count=(/ntimes/)))

iec_d = 0  ! These are flags that, if in your own application, you will not reliably have
icm_c = 0  ! all the data sources available, then set the flag to zero
inc_c = 0  ! so as to tell the downscaling algorithm to skip this data source this time.
iec_d = 0  ! Here they get set to 1 automatically when the data is read in.
icm_c = 0
inc_c = 0
iuk_c = 0
iec_e = 0
inc_e = 0
icm_e = 0
iuk_e = 0

! ---- loop thru the samples, find the day that matches the forecast

ktrday = 1
DO itime = 1, ntimes
   PRINT *,'iyyyymmddhh_anal(itime), iyyyymmddhh = ',iyyyymmddhh_anal(itime), iyyyymmddhh
   IF (iyyyymmddhh_anal(itime) .eq. iyyyymmddhh) THEN
      
      pfcst = -99.99
      cfield = 'fcst_ecmwf_det'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst,&
         start=(/1,1,itime/),count=(/nxf,nyf,1/)))
      CALL make_tiny_negatives_zero(nxf,nyf,pfcst)
      ecmwf_deterministic(:,:) = pfcst(:,:)
      iec_d = 1  

      pfcst = -99.99
      cfield = 'fcst_ncep_cntl'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst,&
         start=(/1,1,itime/),count=(/nxf,nyf,1/)))
      CALL make_tiny_negatives_zero(nxf,nyf,pfcst)
      ncep_control(:,:) = pfcst(:,:)
      inc_c = 1

      pfcst = -99.99
      cfield = 'fcst_cmc_cntl'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst,&
         start=(/1,1,itime/),count=(/nxf,nyf,1/)))
      CALL make_tiny_negatives_zero(nxf,nyf,pfcst)
      cmc_control(:,:) = pfcst(:,:)
      icm_c = 1 

!      pfcst = -99.99
!      cfield = 'fcst_ukmo_cntl'
!      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
!      CALL check (nf90_get_var(netid,ivar,pfcst,&
!         start=(/1,1,itime/),count=(/nxf,nyf,1/)))
!      CALL make_tiny_negatives_zero(nxf,nyf,pfcst)
!      ukmo_control(:,:) = pfcst(:,:)
!      iuk_c = 1 

      ecmwf_ensemble = -99.99
      cfield = 'fcst_ecmwf_ens'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar, ecmwf_ensemble,&
         start=(/1,1,1,itime/), count=(/nxf,nyf,nens_ecmwf,1/)))
      DO imem = 1, nens_ecmwf
         CALL make_tiny_negatives_zero(nxf,nyf,ecmwf_ensemble(1,1,imem))
      END DO
      iec_e = 1 

      ncep_ensemble = -99.99
      cfield = 'fcst_ncep_ens'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar, ncep_ensemble,&
         start=(/1,1,1,itime/), count=(/nxf,nyf,nens_ncep,1/)))
      DO imem = 1, nens_ncep
         CALL make_tiny_negatives_zero(nxf,nyf,ncep_ensemble(1,1,imem))
      END DO
      inc_e = 1 

      cmc_ensemble = -99.99
      cfield = 'fcst_cmc_ens'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar, cmc_ensemble,&
         start=(/1,1,1,itime/), count=(/nxf,nyf,nens_cmc,1/)))
      DO imem = 1, nens_cmc
         CALL make_tiny_negatives_zero(nxf,nyf,cmc_ensemble(1,1,imem))
      END DO
      icm_e = 1

!      ukmo_ensemble = -99.99
!      cfield = 'fcst_ukmo_ens'
!      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
!      CALL check (nf90_get_var(netid, ivar, ukmo_ensemble,&
!         start=(/1,1,1,itime/), count=(/nxf,nyf,nens_ukmo,1/)))
!      DO imem = 1, nens_ukmo
!         CALL make_tiny_negatives_zero(nxf,nyf,ukmo_ensemble(1,1,imem))
!      END DO
!      iuk_e = 1

   ENDIF
END DO

! --- close netcdf file.

CALL check(nf90_close(netid))

DEALLOCATE (iyyyymmddhh_anal)

RETURN
END SUBROUTINE read_forecasts
