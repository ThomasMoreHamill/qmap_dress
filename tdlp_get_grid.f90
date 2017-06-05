subroutine tdlp_get_grid(kunit,id,idate,nx,ny,data_out,ier)
implicit none

integer, parameter :: l3264b=32
integer, parameter :: nbywd=l3264b/8
integer, parameter :: nd5=500000
integer, parameter :: nd7=54

integer, intent(in) :: kunit,nx,ny,idate
character*200 :: infile
integer, intent(in), dimension(4) :: id
real, intent(out), dimension(nx,ny) :: data_out
integer, intent(inout) :: ier

integer :: i,j,ij,ios
logical :: data_found

! ---- Internal Variabled for TDLPACK I/O
integer :: igive,ioctet,misspx,misssx,ntrash
integer, dimension(nd5) :: ipack,iwork
integer, dimension(nd5) :: is0,is1,is2,is4
real, dimension(nd5) :: xdata

! ---- Initialize
data_found=.false.
ier=0
igive=2

do while(.not.data_found)

   if(nbywd.eq.4)then
      read(kunit,iostat=ios)ntrash,ioctet,(ipack(i),i=1,ioctet/nbywd)
   elseif(nbywd.eq.8)then
      read(kunit,iostat=ios)ioctet,(ipack(i),i=1,ioctet/nbywd)
   endif

   ! ---- Check ios.
   !         = -1 - End of file. If data not found, then rewind the file.
   !         >  0 - An error has occurred.
   if(ios.eq.-1)then
      exit
      !if(.not.data_found)rewind(kunit) 
   elseif(ios.gt.0)then
      ier=ios
      return
   endif

   ! ---- Unpack data
   call unpack(6,ipack,iwork,xdata,nd5,is0,is1,is2,is4,nd7,misspx,misssx,igive,l3264b,ier)
   ! ---- Check for matching date and IDs
   if(idate.eq.is1(8).and.&
      id(1).eq.is1(9).and.&
      id(2).eq.is1(10).and.&
      id(3).eq.is1(11).and.&
      id(4).eq.is1(12))then
      do j=1,ny
         do i=1,nx
            ij=(nx*(j-1))+i
            data_out(i,j)=xdata(ij)
         end do
      end do
      data_found=.true.
   endif

end do

call make_tiny_negatives_zero(nx,ny,data_out)

return
end subroutine tdlp_get_grid
