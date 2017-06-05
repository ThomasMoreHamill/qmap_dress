subroutine tdlp_get_vect(kunit,id,idate,nx,ny,data_out,ier)
implicit none

integer, parameter :: l3264b=32
integer, parameter :: nbywd=l3264b/8
integer, parameter :: nd5=500000
integer, parameter :: nd7=54

integer, intent(in) :: kunit,nx,ny,idate
integer, intent(in), dimension(4) :: id
real, intent(out), dimension(nx,ny) :: data_out
integer, intent(inout) :: ier

character(len=4) :: ctemp
character(len=8), dimension(nx*ny) :: ccall
integer :: i,j,ij,ios,ircnt,ix,jy,n
logical :: data_found

! ---- Internal Variabled for TDLPACK I/O
integer :: igive,ioctet,misspx,misssx,ntrash,nsta
integer, dimension(nd5) :: ipack,iwork
integer, dimension(nd5) :: is0,is1,is2,is4
real, dimension(nd5) :: xdata

! ---- Initialize
data_found=.false.
ier=0
igive=2

ircnt=0
do while(.not.data_found)

   ! ---- Read from the input TDLPACK file.
   if(nbywd.eq.4)then
      read(kunit,iostat=ios)ntrash,ioctet,(ipack(i),i=1,ioctet/nbywd)
   elseif(nbywd.eq.8)then
      read(kunit,iostat=ios)ioctet,(ipack(i),i=1,ioctet/nbywd)
   endif
!   PRINT *, 'ioctet = ', ioctet
   ! ---- Check ios.
   !         = -1 - End of file. If data not found, then rewind the file.
   !         >  0 - An error has occurred.
   if(ios.eq.0)then
      ircnt=ircnt+1
   elseif(ios.eq.-1)then
      ier=ios
      exit
      !if(.not.data_found)rewind(kunit) 
   elseif(ios.gt.0)then
      ier=ios
      return
   endif

   ! ---- Check for station call letter record
   CTEMP=TRANSFER(ipack(1),CTEMP)
   IF(CTEMP.EQ.'PLDT')THEN
      ! ---- Record is TDLPACK data, unpack it.
      call unpack(6,ipack,iwork,xdata,nd5,is0,is1,is2,is4,nd7,misspx,misssx,igive,l3264b,ier)
   ELSE
      ! ---- Record is not TDLPACK data, could station call letter record or trailer record.
      if(ipack(1).eq.0)cycle ! trailer record, cycle to the next iteration
      backspace(kunit)
      !PRINT *, 'nbywd = ', nbywd
      if(nbywd.eq.4)then
         read(kunit,iostat=ios)ntrash,ioctet,(ccall(i),i=1,ioctet/8)
      elseif(nbywd.eq.8)then
         read(kunit,iostat=ios)ioctet,(ccall(i),i=1,ioctet/8)
      endif
      nsta=ioctet/8
!      PRINT *, 'nsta = ', nsta
!      write(6,fmt='(14(A8,1X))')(ccall(n),n=1,nsta)
   ENDIF

   ! ---- Check for matching date and IDs
   if(idate.eq.is1(8).and.&
      id(1).eq.is1(9).and.&
      id(2).eq.is1(10).and.&
      id(3).eq.is1(11).and.&
      id(4).eq.is1(12))then
      do n=1,nsta
         read(ccall(n)(1:4),fmt='(I)') ix
         read(ccall(n)(5:8),fmt='(I)') jy
         data_out(ix,jy)=xdata(n)
         !PRINT *, 'data_out = ', data_out(ix, jy), ix, jy
      end do
      data_found=.true.
   endif

end do

call make_tiny_negatives_zero(nx,ny,data_out)

return
end subroutine tdlp_get_vect
