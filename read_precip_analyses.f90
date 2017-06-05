SUBROUTINE read_precip_analyses(nxa,nya,nxf,nyf,nxs,nys,ndays3mo,ipcpvar,iyyyymmddhh,iendhour,&
                                panal_infile,pmask_infile,precip_anal_fine,precip_anal_coarse,&
                                conusmask,ihavedata)

! --- purpose:  read in the 2002-2013 precipitation analyses, both on the original 1/8-degree
!     grid and coarsened to the resolution of the forecast polar-stereo grid system.  Return
!     only the precipitation data relevant for this month and the two surrounding months.
!     Precip data is in a netcdf file, created with python routine ccpa_to_polarstereo.py.

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf, nxs, nys, ndays3mo, iyyyymmddhh, iendhour,ipcpvar
CHARACTER*(*), INTENT(IN) :: panal_infile,pmask_infile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ndays3mo) :: precip_anal_fine ! precip analysis on the 1/8-deg grid
REAL, INTENT(OUT), DIMENSION(nxf,nyf,ndays3mo) :: precip_anal_coarse ! analysis on the 
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask, ihavedata

CHARACTER*20 cfield
CHARACTER*4 cname

INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_end
INTEGER, DIMENSION(366) :: idxuse     
INTEGER, DIMENSION(12) :: imid ! middle julian day of each month
INTEGER, DIMENSION(12) :: ibegin_noleap ! julian day of beginning of month not in leap year
INTEGER, DIMENSION(12) :: ibegin_leap   ! julian day of beginning of month in leap year
REAL, DIMENSION(nxs,nys,ndays3mo) :: precip_anal_fine_s ! precip analysis on the 1/8-deg grid (smaller)
INTEGER*2, DIMENSION(nxs,nys) :: conusmask_s, ccpamask, ihavedata_s

REAL, DIMENSION(nxf,nyf) :: pverif_f
REAL, DIMENSION(nxs,nys) :: pverif_a
INTEGER :: yoffset, xoffset, g, h, i, j, mask_id


DATA imid /15,46,74,105,135,166,196,227,258,288,319,349/ ! middle julian day of each month
DATA ibegin_noleap /1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/
DATA ibegin_leap   /1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336/

! ---- set an index for which julian days of the year to include in output sample for a given month

iyyyymm = iyyyymmddhh / 10000
imonth = iyyyymm - (iyyyymm/100)*100
PRINT *,'imonth = ',imonth

idxuse(:) = 0
IF (imonth .eq. 1) THEN
   idxuse(1:60) = 1
   idxuse(336:365) = 1
ELSE IF (imonth .gt. 1 .and. imonth .lt. 12) THEN
   idxuse(imid(imonth)-45:imid(imonth)+45) = 1
ELSE
   idxuse(305:365) = 1
   idxuse(1:30) = 1
ENDIF   

! ---- open the file, use the number of times in the file later to allocate a date array

netid = 0
PRINT *,'netid, reading from ',netid, TRIM(panal_infile)
CALL check (nf90_open(panal_infile,NF90_NOWRITE,netid))

! ---- read in the list of dates/times in yyyymmddhh format associated with 
!      each time index stored in the netcdf file 

cfield ='time'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_inquire_dimension(netid,ivar,cname,ntimes))
print *,'ntimes = ',ntimes
ALLOCATE (iyyyymmddhh_end(ntimes))

cfield ='yyyymmddhh_date_end'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,iyyyymmddhh_end,&
  start=(/1/),count=(/ntimes/)))

cfield ='conusmask'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,conusmask_s,&
  start=(/1/),count=(/nxs,nys/)))

cfield ='practical_mask'  ! this is the list of CCPA points where there is valid data
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,ihavedata_s,&
  start=(/1/),count=(/nxs,nys/)))
!print *,'masks read in'

! ---- loop thru the samples, and if this day is flagged for being read in, do so

ktrday = 1
precip_anal_fine = -99.99   ! initialize to missing data value
precip_anal_coarse = -99.99
precip_anal_fine_s = -99.99

DO itime = 1, ntimes
   iyear   = iyyyymmddhh_end(itime) / 1000000
   immddhh = iyyyymmddhh_end(itime) - iyear*1000000
   imo     = immddhh / 10000
   iddhh   = immddhh - imo*10000
   iday    = iddhh / 100
   ihour   = iddhh -  iday*100

   IF (MOD(iyear,4) .eq. 0) THEN
      ijulday = ibegin_leap(imo) + iday - 1
   ELSE
      ijulday = ibegin_noleap(imo) + iday -1
   ENDIF
   IF (idxuse(ijulday) .eq. 1 .and. iyear .ne. 2016 .and. iendhour .eq. ihour) THEN
 
      ! --- read in precipitation analysis information
      PRINT *, 'ktrday,itime = ', ktrday, itime
      pverif_a = -99.99
      if(ipcpvar.eq.0)cfield='apcp_anal_12h'  ! ---- raw 1/8 degree data
      if(ipcpvar.eq.1)cfield='apcp_anal_06h'  ! ---- raw 1/8 degree data
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pverif_a,&
         start=(/1,1,itime/),count=(/nxs,nys,1/)))
      precip_anal_fine_s(:,:,ktrday) = pverif_a(:,:)

      pverif_f = -99.99
      if(ipcpvar.eq.0)cfield='apcp_anal_ps_12h' ! ---- 47 or 95 km upscaled polar stereo data
      if(ipcpvar.eq.1)cfield='apcp_anal_ps_06h' ! ---- 47 or 95 km upscaled polar stereo data
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pverif_f,&
         start=(/1,1,itime/),count=(/nxf,nyf,1/)))
      precip_anal_coarse(:,:,ktrday) = pverif_f(:,:)

      ! ---- increment counter

      ktrday = ktrday+1

   ENDIF ! (idxuse(iday) .eq. 1) 

END DO  ! itime = 1, ntimes
ktrday = ktrday-1

! --- close netcdf file.

CALL check(nf90_close(netid))
DEALLOCATE (iyyyymmddhh_end)

! ---- put on expanded grid & pad with missings


!DO i = 1,nxs
!   DO j = 1,nys
!      IF (ccpamask(i,j) .eq. 0 ) THEN
!          precip_anal_fine_s(i,j,:) = -99.99
!      ELSE
!          precip_anal_fine_s(i,j,:) = precip_anal_fine_s(i,j,:)
!      ENDIF
!   END DO
!END DO

yoffset=22
xoffset=21
conusmask = -99.99
precip_anal_fine = -99.99
ihavedata = 0
DO i=xoffset,nxs+xoffset-1
   DO j=yoffset,nys+yoffset-1
      g=i-xoffset+1
      h=j-yoffset+1
      precip_anal_fine(i,j,:)=precip_anal_fine_s(g,h,:)
      conusmask(i,j)=conusmask_s(g,h)
      ihavedata(i,j)=ihavedata_s(g,h)
   END DO
END DO


RETURN
END subroutine read_precip_analyses
