subroutine read_cdf_netcdf(nxa,nya,npct,ipcpvar,thresh,precip_anal_cdf,&
                           ecmwf_deterministic_cdf,ncep_control_cdf,&
                           cmc_control_cdf,ecmwf_ensemble_cdf,ncep_ensemble_cdf,&
                           cmc_ensemble_cdf,gem_cdf)

! ---- Purpose:  Read CDF data from NetCDF

use netcdf
implicit none

integer, intent(in) :: nxa,nya,npct,ipcpvar
real, intent(out), dimension(npct) :: thresh
real*8, intent(out), dimension(nxa,nya,npct) :: precip_anal_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: ecmwf_deterministic_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: ncep_control_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: cmc_control_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: ecmwf_ensemble_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: ncep_ensemble_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: cmc_ensemble_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: gem_cdf

!real*8, dimension(nxa,nya,npct) :: ukmo_control_cdf
!real*8, dimension(nxa,nya,npct) :: ukmo_ensemble_cdf

! ---- Local Variables
character(len=:), allocatable :: cvar
character(len=256) :: cdf_infile
integer :: ios,idimid,ivarid,ncid,nxain,nyain,npctin

! ---- Initialize
ios=0

! ---- Get input CDF NetCDF filename from Fortran unit number environment variable
call get_environment_variable('FORT30',cdf_infile)

! ---- Open the file, use the number of times in the file later to allocate a date array
write(6,*)'Opening NetCDF file, ',trim(cdf_infile),' to read CDF data'
call check(nf90_open(trim(cdf_infile),NF90_NOWRITE,ncid))
write(6,*)'NetCDF file opened with ncid = ',ncid

! ---- Check NetCDF dimension values against nxa,nya,npct
cvar="nxa"
call check(nf90_inq_dimid(ncid,cvar,idimid))
call check(nf90_inquire_dimension(ncid,idimid,cvar,nxain))
print *,'nxa = ',nxain
cvar="nya"
call check(nf90_inq_dimid(ncid,cvar,idimid))
call check(nf90_inquire_dimension(ncid,idimid,cvar,nyain))
print *,'nya = ',nyain
cvar="npct"
call check(nf90_inq_dimid(ncid,cvar,idimid))
call check(nf90_inquire_dimension(ncid,idimid,cvar,npctin))
print *,'npctin = ',npctin
if(nxain.ne.nxa.or.nyain.ne.nya.or.npctin.ne.npct)then
   write(6,*)'****ERROR IN READ_CDF_NETCDF. NETCDF DIMENSIONS DO NOT MATCH EXPECTED VALUES.'
   call w3tage('BLEND_PRECIP_DOWNSCALE')
   stop
endif 

! ---- Read CDF data
cvar="thresh"
call check(nf90_inq_varid(ncid,cvar,ivarid))
call check(nf90_get_var(ncid,ivarid,thresh))
PRINT *,'thresh = ',thresh

! ---- Read CDF data according to ipcpvar.
if(ipcpvar.eq.0)then

   cvar="precip_anal_cdf"
   call check(nf90_inq_varid(ncid,cvar,ivarid))
   call check(nf90_get_var(ncid,ivarid,precip_anal_cdf))
   PRINT *,'precip_anal_cdf = ',precip_anal_cdf(nxa/2,nya/2,:)

   cvar="cmc_control_cdf"
   call check(nf90_inq_varid(ncid,cvar,ivarid))
   call check(nf90_get_var(ncid,ivarid,cmc_control_cdf))
   print *,'cmc_control_cdf = ',cmc_control_cdf(nxa/2, nya/2, :)

   cvar="cmc_ensemble_cdf"
   call check(nf90_inq_varid(ncid,cvar,ivarid))
   call check(nf90_get_var(ncid,ivarid,cmc_ensemble_cdf))
   print *,'cmc_ensemble_cdf = ',cmc_ensemble_cdf(nxa/2,nya/2,:)

   cvar="ncep_control_cdf"
   call check(nf90_inq_varid(ncid,cvar,ivarid))
   call check(nf90_get_var(ncid,ivarid,ncep_control_cdf))
   print *,'ncep_control_cdf = ',ncep_control_cdf(nxa/2,nya/2,:)

   cvar="ncep_ensemble_cdf"
   call check(nf90_inq_varid(ncid,cvar,ivarid))
   call check(nf90_get_var(ncid,ivarid,ncep_ensemble_cdf))
   print *,'ncep_ensemble_cdf = ',ncep_ensemble_cdf(nxa/2,nya/2,:)

elseif(ipcpvar.eq.1)then

   cvar="gem_cdf"
   call check(nf90_inq_varid(ncid,cvar,ivarid))
   call check(nf90_get_var(ncid,ivarid,gem_cdf))

endif

! ---- Close NetCDF file.
CALL check(nf90_close(ncid))

! ---- Set ECMWF stuff to missing.
ecmwf_ensemble_cdf=-99.99
ecmwf_deterministic_cdf=-99.99

return
end subroutine read_cdf_netcdf
