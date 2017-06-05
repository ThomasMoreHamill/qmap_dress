subroutine read_forecasts_tdlp(nxf,nyf,nens_cmc,nens_ncep,ipcpvar,idate,ifhr,cmc_ctl,ncep_ctl,&
                                    cmc_ens,ncep_ens,icm_c,inc_c,icm_e,inc_e,ier)
implicit none
! ---- Input/Output Variables
integer, intent(in) :: nxf,nyf,nens_cmc,nens_ncep,idate,ifhr,ipcpvar
real, intent(out), dimension(nxf,nyf) :: cmc_ctl,ncep_ctl
real, intent(out), dimension(nxf,nyf,nens_cmc) :: cmc_ens
real, intent(out), dimension(nxf,nyf,nens_ncep) :: ncep_ens
integer, intent(out) :: icm_c,inc_c,icm_e,inc_e
integer, intent(out) :: ier

! ---- Local Variables
integer :: dd,id1,j,kunit
integer, save :: ifirst
integer, dimension(4) :: id
real, dimension(nxf,nyf) :: xdata
character(len=200) :: fname, cdet_infile2

data ifirst/0/

! ---- Initialize
ier=0
icm_c=0
icm_e=0
inc_c=0
inc_e=0

! ---- id(2:4) are independent of model data.  id(2) = 0; id(3) = forecast projection; and 
!      id(4) = 300 since these data have been interpolated from 47km grid PS to 1/8 deg. LL.
if(ipcpvar.eq.0) id1=003220000 ! 12-hr Precip Amount
if(ipcpvar.eq.1) id1=003210000 ! 6-Precip Amount
id(2)=0
id(3)=ifhr
id(4)=300

cmc_ctl(:,:)=-99.99
cmc_ens(:,:,:)=-99.99
ncep_ctl(:,:)=-99.99
ncep_ens(:,:,:)=-99.99

write(6,100)idate   
100 format(/'Reading 47km Gridded TDLPACK model data for date = ',I10) 

! ---- Get CMC Deterministic (Unit = 50), Model DD = 63.
kunit=50
if(ifirst.eq.0)open(unit=kunit,form='unformatted',status='old',convert='big_endian')
dd=63
id(1)=id1+dd
call tdlp_get_grid(kunit,id,idate,nxf,nyf,xdata,ier)
if(ier.eq.-1)then
   close(kunit)
elseif(ier.eq.0)then
   write(6,110)id(:),' found for CMC Control for date ',idate
   110 format('ID ',3(I9.9,1X),I10.10,A,I10.10) 
   cmc_ctl(:,:)=xdata(:,:)
   icm_c=1 
else
   write(6,110)id(:),' missing for CMC Control for date ',idate
endif

! ---- Get CMC Ensemble (Unit = 51); 20 members, Model DD = 41-60
kunit=51
if(ifirst.eq.0)open(unit=kunit,form='unformatted',status='old',convert='big_endian')
do j=1,nens_cmc
   dd=j+40
   id(1)=id1+dd
   call tdlp_get_grid(kunit,id,idate,nxf,nyf,xdata,ier)
   if(ier.eq.-1)then
      close(kunit)
   elseif(ier.eq.0)then
      write(6,110)id(:),' found for CMC Ensemble for date ',idate
      cmc_ens(:,:,j)=xdata(:,:)
      icm_e=1
   else 
      write(6,110)id(:),' missing for CMC Ensemble for date ',idate
   endif
end do

! ---- Get ECMWF Deterministic (Unit = 52)
! ---- Get ECMWF Ensemble (Unit = 53)

! ---- Get NCEP Control (i.e. GFS) (Unit = 54), Model DD = 08
kunit=54
if(ifirst.eq.0)open(unit=kunit,form='unformatted',status='old',convert='big_endian')
dd=08
id(1)=id1+dd
call tdlp_get_grid(kunit,id,idate,nxf,nyf,xdata,ier)
if(ier.eq.-1)then
   close(kunit)
elseif(ier.eq.0)then
   write(6,110)id(:),' found for NCEP Control for date ',idate
   ncep_ctl(:,:)=xdata(:,:)
   inc_c=1
else
   write(6,110)id(:),' missing for NCEP Control for date ',idate
endif

! ---- Get NCEP Ensemble (i.e. GEFS) (Unit = 55); 20 members, Model DD = 11-30
kunit=55
if(ifirst.eq.0)open(unit=kunit,form='unformatted',status='old',convert='big_endian')
do j=1,nens_ncep
   dd=j+10
   id(1)=id1+dd
   call tdlp_get_grid(kunit,id,idate,nxf,nyf,xdata,ier)
   if(ier.eq.-1)then
      close(kunit)
   elseif(ier.eq.0)then
      write(6,110)id(:),' found for NCEP Ensemble for date ',idate
      ncep_ens(:,:,j)=xdata(:,:)
      inc_e=1
   else
      write(6,110)id(:),' missing for NCEP Ensemble for date ',idate
   endif
end do

! ---- Get UKMet Deterministic (Unit = 56)
! ---- Get UKMet Ensemble (Unit = 57)

ifirst=ifirst+1
return
end subroutine read_forecasts_tdlp
