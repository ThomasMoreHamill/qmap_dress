
subroutine read_cdf_netcdf(nxa,nya,npct,nens_cmc,ipcpvar,thresh,precip_anal_cdf,&
    ecmwf_ensemble_cdf,ncep_ensemble_cdf,cmc_ensemble_cdf,gem_cdf,rlatsa, rlonsa)

! ---- Purpose:  Read CDF data from NetCDF

use netcdf
implicit none

integer, intent(in) :: nxa,nya,npct,nens_cmc,ipcpvar
real, intent(out), dimension(npct) :: thresh
real*8, intent(out), dimension(nxa,nya,npct) :: precip_anal_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: ecmwf_ensemble_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: ncep_ensemble_cdf
real*8, intent(out), dimension(nxa,nya,nens_cmc,npct) :: cmc_ensemble_cdf
real*8, intent(out), dimension(nxa,nya,npct) :: gem_cdf
real*4, intent(out), dimension(nxa,nya) :: rlatsa
real*4, intent(out), dimension(nxa,nya) :: rlonsa

! ---- Local Variables
character(len=:), allocatable :: cvar
character(len=256) :: cdf_infile
integer :: ios,idimid,ivarid,ncid,nxain,nyain,npctin,nens_cmcin

! ---- Initialize
ios=0
print *,'ipcpvar in read_cdf_netcdf = ',ipcpvar

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
cvar="nya"
call check(nf90_inq_dimid(ncid,cvar,idimid))
call check(nf90_inquire_dimension(ncid,idimid,cvar,nyain))
cvar="npct"
call check(nf90_inq_dimid(ncid,cvar,idimid))
call check(nf90_inquire_dimension(ncid,idimid,cvar,npctin))

if (nxain .ne. nxa .or. nyain .ne. nya .or. npctin .ne. npct) then
   write(6,*)'****ERROR IN READ_CDF_NETCDF. NETCDF DIMENSIONS DO NOT MATCH EXPECTED VALUES.'
   write(6,*)'nxa, nya, npct = ',nxa, nya, npct
   write(6,*)'nxain, nyain, npctin = ',nxain, nyain, npctin
   call w3tage('BLEND_PRECIP_DOWNSCALE')
   stop
endif 

! ---- Check nens_cmc dimension if input CDF NetCDF file is POP12.

if(ipcpvar.eq.0)then
   cvar="nens_cmc"
   call check(nf90_inq_dimid(ncid,cvar,idimid))
   call check(nf90_inquire_dimension(ncid,idimid,cvar,nens_cmcin))
   if(nens_cmcin.ne.nens_cmc)then
      write(6,*)'****ERROR IN READ_CDF_NETCDF. CMC ENSEMBLE COUNT DIM DOES NOT MATCH EXPECTED VALUE.'
      call w3tage('BLEND_PRECIP_DOWNSCALE')
      stop
   endif
endif 

! ---- Read CDF data

cvar="thresh"
call check(nf90_inq_varid(ncid,cvar,ivarid))
call check(nf90_get_var(ncid,ivarid,thresh))
write(6,*)'Read CDF Threshold values.'

! ---- Read CDF data according to ipcpvar.

print *,'in read_cdf_netcdf, ipcpvar = ',ipcpvar

cvar="precip_anal_cdf"
call check(nf90_inq_varid(ncid,cvar,ivarid))
call check(nf90_get_var(ncid,ivarid,precip_anal_cdf))
write(6,*)'Read Precip Analysis CDF.'

if(ipcpvar.eq.0)then

    cvar="ecmwf_ensemble_cdf"
    call check(nf90_inq_varid(ncid,cvar,ivarid))
    call check(nf90_get_var(ncid,ivarid,ncep_ensemble_cdf))
    write(6,*)'Read GEFS CDF.'

    cvar="ncep_ensemble_cdf"
    call check(nf90_inq_varid(ncid,cvar,ivarid))
    call check(nf90_get_var(ncid,ivarid,ncep_ensemble_cdf))
    write(6,*)'Read GEFS CDF.'

    ! ---- CMC Ensemble CDF is a 4D array. The start and count
    !      arguments are needed.  Otherwise, reading fails.

    cvar="cmc_ensemble_cdf"
    call check(nf90_inq_varid(ncid,cvar,ivarid))
    call check(nf90_get_var(ncid,ivarid,cmc_ensemble_cdf,start=(/1,1,1,1/),count=(/nxa,nya,nens_cmc,npct/)))
    write(6,*)'Read CMC Ensemble CDF.'

else if (ipcpvar.eq.1)then

    cvar="gem_cdf"
    call check(nf90_inq_varid(ncid,cvar,ivarid))
    call check(nf90_get_var(ncid,ivarid,gem_cdf))
    write(6,*)'Read GEFS CDF.'

endif

cvar="rlatsa"
call check(nf90_inq_varid(ncid,cvar,ivarid))
call check(nf90_get_var(ncid,ivarid,rlatsa))
write(6,*)'Read output grid latitudes'

cvar="rlonsa"
call check(nf90_inq_varid(ncid,cvar,ivarid))
call check(nf90_get_var(ncid,ivarid,rlonsa))
write(6,*)'Read output grid longitudes'

! ---- Close NetCDF file.

CALL check(nf90_close(ncid))

return
end subroutine read_cdf_netcdf
