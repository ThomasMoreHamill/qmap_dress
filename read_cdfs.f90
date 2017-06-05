subroutine read_cdfs(nxa,nya,npct,ipcpvar,infile,thresh,precip_anal_cdf,&
                     ecmwf_deterministic_cdf,ncep_control_cdf,&
                     cmc_control_cdf,ecmwf_ensemble_cdf,ncep_ensemble_cdf,&
                     cmc_ensemble_cdf,gem_cdf)
! ---- purpose:  read in the cdfs previously calculated in create_fcst_anal_cdfs.f90
implicit none

integer, intent(in) :: nxa,nya,npct,ipcpvar
real, intent(out), dimension(npct) :: thresh
character(len=*), intent(in) :: infile
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
integer :: ios
integer :: nxain,nyain,npctin
integer(kind=2), dimension(nxa,nya) :: conusmask

! ---- Initialize
ios=0

open(unit=30,status='old',form='unformatted',iostat=ios)
if(ipcpvar.eq.0)then ! CDFs for POP12 for all models
   read(30,iostat=ios)nxain,nyain,npctin
   read(30,iostat=ios)thresh
   read(30,iostat=ios)precip_anal_cdf
   read(30,iostat=ios)ncep_control_cdf
   read(30,iostat=ios)cmc_control_cdf
   read(30,iostat=ios)ncep_ensemble_cdf
   read(30,iostat=ios)cmc_ensemble_cdf
   read(30,iostat=ios)conusmask
   !read(30,iostat=ios)rlonsa
   !read(30,iostat=ios)rlatsa
elseif(ipcpvar.eq.1)then ! CDFs for QPF06 for gem
   read(30,iostat=ios)nxain,nyain,npctin
   read(30,iostat=ios)thresh
   read(30,iostat=ios)precip_anal_cdf
   read(30,iostat=ios)gem_cdf
   read(30,iostat=ios)conusmask
   !read(30,iostat=ios)rlonsa
   !read(30,iostat=ios)rlatsa
endif
close(30,iostat=ios)

! ---- Set ECMWF stuff to missing.
ecmwf_ensemble_cdf=-99.99
ecmwf_deterministic_cdf=-99.99

return
end subroutine read_cdfs
