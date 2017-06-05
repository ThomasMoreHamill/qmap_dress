SUBROUTINE read_precip_analyses_x9(nxa,nya,nxs,nys,ndays3mo,ipcpvar,iyyyymmddhh,iendhour,&
    panal_infile,pmask_infile,precip_anal_fine,conusmask)

! --- purpose:  Read in the 2002-2016 (or possibly later) 
!     precipitation analyses, both on the original 1/8-degree
!     grid and coarsened to the resolution of the forecast polar-stereo grid system.  Return
!     only the precipitation data relevant for this month and the two surrounding months.
!     Precip data is in a netcdf file, created with python routine ccpa_to_polarstereo.py.
!
!     Modified May 2016 so that the mask files are updated to trim US West coast data.

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, nxs, nys, ndays3mo, iyyyymmddhh, iendhour, ipcpvar
CHARACTER*(*), INTENT(IN) :: panal_infile,pmask_infile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ndays3mo) :: precip_anal_fine ! precip analysis on the 1/8-deg grid
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask

CHARACTER*20 cfield
CHARACTER*4 cname

INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_end
INTEGER, DIMENSION(366) :: idxuse     
INTEGER, DIMENSION(12) :: imid ! middle julian day of each month
INTEGER, DIMENSION(12) :: ibegin_noleap ! julian day of beginning of month not in leap year
INTEGER, DIMENSION(12) :: ibegin_leap   ! julian day of beginning of month in leap year
REAL, DIMENSION(nxs,nys,ndays3mo) :: precip_anal_fine_s ! precip analysis on the 1/8-deg grid (smaller)
INTEGER*2, DIMENSION(nxs,nys) :: conusmask_s, ccpamask

REAL, DIMENSION(nxs,nys) :: lons_anal
REAL, DIMENSION(nxs,nys) :: lats_anal


REAL, DIMENSION(nxs,nys) :: pverif_a
INTEGER :: i,j,ii,jj,xoffset,yoffset

DATA imid /15,46,74,105,135,166,196,227,258,288,319,349/ ! middle julian day of each month
DATA ibegin_noleap /1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/
DATA ibegin_leap   /1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336/

! ---- Initialize

conusmask(:,:)=0
conusmask_s(:,:)=0

precip_anal_fine(:,:,:)=-99.99
precip_anal_fine_s(:,:,:)=-99.99

! ---- Set an index for which julian days of the year to include in output sample for a given month
iyyyymm=iyyyymmddhh/10000
imonth=iyyyymm-(iyyyymm/100)*100
write(6,*)'imonth = ',imonth

idxuse(:)=0
if(imonth.eq.1)then
   idxuse(1:60)=1
   idxuse(336:365)=1
elseif(imonth.gt.1.and.imonth.lt.12) then
   idxuse(imid(imonth)-45:imid(imonth)+45)=1
else
   idxuse(305:365)=1
   idxuse(1:30)=1
ENDIF   

! ---- Open the file, use the number of times in the file later to allocate a date array
netid=0
PRINT *,'netid, reading from ',netid, TRIM(panal_infile)
CALL check (nf90_open(panal_infile,NF90_NOWRITE,netid))

! ---- Read in the list of dates/times in yyyymmddhh format associated with 
!      each time index stored in the netcdf file 
cfield ='time'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_inquire_dimension(netid,ivar,cname,ntimes))
print *,'ntimes=',ntimes
ALLOCATE (iyyyymmddhh_end(ntimes))

cfield='yyyymmddhh_date_end'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,iyyyymmddhh_end,&
           start=(/1/),count=(/ntimes/)))

cfield='conusmask'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,conusmask_s,&
           start=(/1/),count=(/nxs,nys/)))

cfield='lons_anal'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,lons_anal,&
           start=(/1/),count=(/nxs,nys/)))

cfield='lats_anal'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,lats_anal,&
           start=(/1/),count=(/nxs,nys/)))

! ---- Loop thru the samples, and if this day is flagged for being read in, do so
ktrday=1
precip_anal_fine=-99.99   ! initialize to missing data value
precip_anal_fine_s=-99.99

DO itime=1, ntimes
   iyear  =iyyyymmddhh_end(itime) / 1000000
   immddhh=iyyyymmddhh_end(itime) - iyear*1000000
   imo    =immddhh / 10000
   iddhh  =immddhh - imo*10000
   iday   =iddhh / 100
   ihour  =iddhh -  iday*100

   IF (MOD(iyear,4) .eq. 0) THEN
      ijulday=ibegin_leap(imo) + iday - 1
   ELSE
      ijulday=ibegin_noleap(imo) + iday -1
   ENDIF
   IF (idxuse(ijulday) .eq. 1 .and. iyear .ne. 2016 .and. iendhour .eq. ihour) THEN
 
      ! ---- Read precipitation analysis data on 1/8 deg. grid
      !PRINT *, 'ktrday,itime=', ktrday, itime
      pverif_a=-99.99
      if(ipcpvar.eq.0)cfield='apcp_anal_12h'  ! ---- 12-hr APCP, raw 1/8 degree data
      if(ipcpvar.eq.1)cfield='apcp_anal_06h'  ! ---- 06-hr APCP, raw 1/8 degree data
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pverif_a,&
         start=(/1,1,itime/),count=(/nxs,nys,1/)))
      precip_anal_fine_s(:,:,ktrday)=pverif_a(:,:)

      ! ---- Increment counter
      ktrday=ktrday+1

   ENDIF ! (idxuse(iday) .eq. 1) 

END DO  ! itime=1, ntimes
ktrday=ktrday-1

! ---- Close netcdf file.
CALL check(nf90_close(netid))
DEALLOCATE (iyyyymmddhh_end)

! ---- Perform regridding of precip analysis and conusmask from the original ccpa grid to
!      the expanded. Here we iterate over the original grid (464x224) and create ii,jj
!      coordinates for the expanded grid (515x262).  Using this method, we should not
!      eclipse array bounds.  Also note that the expanded grid is first initialized
!      to missing, so those gridpoints outside the original CCPA will remain missing.

xoffset=21
yoffset=22
do j=1,nys
   do i=1,nxs
      ii=i+xoffset
      jj=j+yoffset
      precip_anal_fine(ii,jj,:)=precip_anal_fine_s(i,j,:)
      conusmask(ii,jj)=conusmask_s(i,j)
   end do
end do

print *,'done reading read_precip_analyses_x9'
RETURN
END subroutine read_precip_analyses_x9
