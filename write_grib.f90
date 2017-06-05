SUBROUTINE WRITE_GRIB(iyyyymmddhh,ileadb,ileade,nxa,nya,raw,ipcpvar,KFILGO)
IMPLICIT NONE

INTEGER, PARAMETER :: lcgrib = 6000000
INTEGER, PARAMETER :: numcoord = 0

INTEGER, INTENT(IN) :: ileadb, ileade, iyyyymmddhh, nxa, nya, ipcpvar
INTEGER, DIMENSION(2) :: listsec0
INTEGER, DIMENSION(13) :: listsec1
INTEGER, DIMENSION(5) :: igds
INTEGER, DIMENSION(20) :: igdstmpl
INTEGER, DIMENSION(36) :: ipdstmpl
INTEGER, DIMENSION(18) :: idrstmpl
REAL, DIMENSION(numcoord) :: coordlist
REAL, DIMENSION(nxa*nya) :: raw, fld
INTEGER :: missVal
CHARACTER(LEN=1) :: cgrib(lcgrib)

INTEGER :: ierr, i
INTEGER :: iyear, imonth, iday, ihour, idoy
INTEGER :: jyear, jmonth, jday, jhour, jdoy
INTEGER :: ipdsnum, idrsnum, ibmap
INTEGER :: igdstmplen, ipdstmplen, idrstmplen
INTEGER :: ideflist, idefnum, ngrdpts
INTEGER :: lengrib
INTEGER :: KFILGO
INTEGER :: jyyyymmddhh
INTEGER :: ileadb_s

LOGICAL*1 bmap(nxa*nya)

ierr=0
ngrdpts=nxa*nya

! Set Section 0 and Section 1 Arrays
listsec0(1) = 0           ! Discipline-GRIB Master Table Number
listsec0(2) = 2           ! GRIB Edition Number (Currently 2)

! Parsing of iyyyymmddhh for Reference Time
call doy(iyyyymmddhh,iyear,imonth,iday,ihour,idoy)

! Determine valid date
call updat(iyyyymmddhh,ileade,jyyyymmddhh)
call doy(jyyyymmddhh,jyear,jmonth,jday,jhour,jdoy)

! Now fill listsec1()
listsec1(1) = 7           ! ID of Originating Center
listsec1(2) = 14          ! ID of Originating Sub-Center
listsec1(3) = 1           ! GRIB Master Tables Version Number
listsec1(4) = 0           ! GRIB Local Tables Version Number
listsec1(5) = 1           ! Significance of Reference Time
listsec1(6) = iyear       ! Reference Time (4-Digit Year)
listsec1(7) = imonth      ! Reference Time (Month)
listsec1(8) = iday        ! Reference Time (Day)
listsec1(9) = ihour       ! Reference Time (Hour)
listsec1(10) = 0          ! Reference Time (Minute)
listsec1(11) = 0          ! Reference Time (Second)
listsec1(12) = 0          ! Production Status of Data
listsec1(13) = 1          ! Type of Processed Data

! Set Section 3 Array
igds(1) = 0               ! Source of Grid Definition
igds(2) = nxa*nya         ! Number of Points in the Defined Grid
igds(3) = 0               ! Number of Octets Needed for Additional Grid Def
igds(4) = 0               ! Interpretation of List for Optional Points Def
igds(5) = 0               ! Grid Definition Template Number

! Below values for igdstmpl() assume use of Template 3.0 for Lat/Lon Grids.  May
! want to expand this further in the future.
igdstmpl(1) = 6           ! Shape of the Earth
igdstmpl(2) = 0           ! Scale Factor of Radius of Spherical Earth
igdstmpl(3) = 6371200     ! Scaled Value of Radius of Spherical Earth
igdstmpl(4) = 255         ! Scale Factor of Major Axis of Oblate Spheroid
igdstmpl(5) = 255         ! Scaled Value of Major Axis of Oblate Spheroid
igdstmpl(6) = 255         ! Scale Factor of Minor Axis of Oblate Spheroid
igdstmpl(7) = 255         ! Scaled Value of Minor Axis of Oblate Spheroid
igdstmpl(8) = nxa         ! Ni - Number of Points Along a Parallel
igdstmpl(9) = nya         ! Nj - Number of Points Alogn a Meridian
igdstmpl(10) = 0          ! Basic Angle of Initial Production Domain
igdstmpl(11) = 0          ! Subdivisions of Basic Angle
igdstmpl(12) = 22313000   ! La1 - Latitude of First Grid Point
igdstmpl(13) = 232437000  ! Lo1 - Longitude of First Grid Point
igdstmpl(14) = 48         ! Resolution and Component Flags
igdstmpl(15) = 54938000   ! La2 - Latitude of Last Grid Point
igdstmpl(16) = 296687000  ! Lo2 - Longitude of Last Grid Point
igdstmpl(17) = 125000     ! Di - I Direction Increment
igdstmpl(18) = 125000     ! Dj - J Direction Increment
igdstmpl(19) = 64         ! Scanning Mode
igdstmpl(20) = 255        ! List of Number Pts for Quasi-Regular Grid
igdstmplen = 20           ! Max dimension of igdstmpl()

ileadb_s=ileadb
CALL MKIEEE(FLOAT(ileadb_s),ileadb_s,1)

! Assign ending of time interval for template 4.9
! ipdstmpl(23:28) settings copied from u13xlib/packgrib2.f
!
! Dummy settings from pop12 Section 4 header file, commented out now.
!ipdstmpl(23) = 65535      ! Year of End of Overall Time Interval
!ipdstmpl(24) = 255        ! Month of End of Overall Time Interval
!ipdstmpl(25) = 255        ! Day of End of Overall Time Interval
!ipdstmpl(26) = 255        ! Hour of End of Overall Time Interval
!ipdstmpl(27) = 255        ! Minute of End of Overall Time Interval
!ipdstmpl(28) = 255        ! Second of End of Overall Time Interval
!
! Exact MOS20000 CODEBLOCK:
!            IPDSTMPL(23)=JDATE/1000000
!            IPDSTMPL(24)=JDATE/10000-IPDSTMPL(23)*100
!            IPDSTMPL(25)=JDATE/100-IPDSTMPL(23)*10000-
!     1                             IPDSTMPL(24)*100
!            IPDSTMPL(26)=JDATE-IPDSTMPL(23)*1000000-
!     1                         IPDSTMPL(24)*10000-IPDSTMPL(25)*100
!            IPDSTMPL(27)=0  !MINUTES ARE SET TO ZERO.
!            IPDSTMPL(28)=0  !SECONDS ARE SET TO ZERO.
!

if(ipcpvar.eq.0)then

! Set Section 4 Array
ipdsnum = 9
ipdstmpl(1) = 1           ! Parameter Category
ipdstmpl(2) = 8           ! Parameter Number
ipdstmpl(3) = 2           ! Type of Generating Process
ipdstmpl(4) = 255         ! Background Generating Process Identifier
ipdstmpl(5) = 104         ! Forecast Generating Process Identifier
ipdstmpl(6) = 65535       ! Hours after Reference Time Data Cutoff
ipdstmpl(7) = 255         ! Minutes after Reference Time Data Cutoff
ipdstmpl(8) = 1           ! Indicator of Unit of Time Range
ipdstmpl(9) = ileadb      ! Forecast Time in Units Above
ipdstmpl(10) = 1          ! Type of First Fixed Surface
ipdstmpl(11) = 0          ! Scale Factor of First Fixed Surface
ipdstmpl(12) = 0          ! Scaled Value of First Fixed Surface
ipdstmpl(13) = 255        ! Type of Second Fixed Surface
ipdstmpl(14) = 0          ! Scale Factor of Second Fixed Surface
ipdstmpl(15) = 0          ! Scaled Vaule of Second Fixed Surface
ipdstmpl(16) = 255        ! Forecast Probability Number
ipdstmpl(17) = 255        ! Total Number of Forecast Probabilities
ipdstmpl(18) = 1          ! Probability Type
ipdstmpl(19) = 255        ! Scale Factor of Lower Limit
ipdstmpl(20) = 255        ! Scaled Value of Lower Limit
ipdstmpl(21) = 3          ! Scale Factor of Upper Limit
ipdstmpl(22) = 254        ! Scaled Value of Upper Limit
ipdstmpl(23) = jyear      ! Year - Time of end of overall time interval
ipdstmpl(24) = jmonth     ! Month - Time of end of overall time interval
ipdstmpl(25) = jday       ! Day - Time of end of overall time interval
ipdstmpl(26) = jhour      ! Hour - Time of end of overall time interval
ipdstmpl(27) = 0          ! Minutes (set to zero).
ipdstmpl(28) = 0          ! Seconds (set to zero).
ipdstmpl(29) = 1          ! Number of Time Range Specifications fo Field
ipdstmpl(30) = 0          ! Total Number of Data Values Missing
ipdstmpl(31) = 1          ! Statistical Process Used for Field
ipdstmpl(32) = 255        ! Type of Time Increment Between Fields
ipdstmpl(33) = 1          ! Indicator of Unit of Time Range for Processing
ipdstmpl(34) = 12         ! Length of Time Range for Processing
ipdstmpl(35) = 1          ! Indicator of Unit of Time Between Fields
ipdstmpl(36) = 0          ! Time Increment Between Successive Fields
ipdstmplen = 36           ! Max dimension of ipdstmpl()

elseif(ipcpvar.eq.1)then  ! QPF

! Set Section 4 Array
ipdsnum = 8               ! PDS Template for QPF is 8.
ipdstmpl(1) = 1           ! Parameter Category
ipdstmpl(2) = 8           ! Parameter Number
ipdstmpl(3) = 2           ! Type of Generating Process
ipdstmpl(4) = 255         ! Background Generating Process Identifier
ipdstmpl(5) = 104         ! Forecast Generating Process Identifier
ipdstmpl(6) = 65535       ! Hours after Reference Time Data Cutoff
ipdstmpl(7) = 255         ! Minutes after Reference Time Data Cutoff
ipdstmpl(8) = 1           ! Indicator of Unit of Time Range
ipdstmpl(9) = ileadb      ! Forecast Time in Units Above
ipdstmpl(10) = 1          ! Type of First Fixed Surface
ipdstmpl(11) = 0          ! Scale Factor of First Fixed Surface
ipdstmpl(12) = 0          ! Scaled Value of First Fixed Surface
ipdstmpl(13) = 255        ! Type of Second Fixed Surface
ipdstmpl(14) = 0          ! Scale Factor of Second Fixed Surface
ipdstmpl(15) = 0          ! Scaled Vaule of Second Fixed Surface
ipdstmpl(16) = jyear      ! Year - Time of end of overall time interval
ipdstmpl(17) = jmonth     ! Month - Time of end of overall time interval
ipdstmpl(18) = jday       ! Day - Time of end of overall time interval
ipdstmpl(19) = jhour      ! Hour - Time of end of overall time interval
ipdstmpl(20) = 0          ! Minutes (set to zero).
ipdstmpl(21) = 0          ! Seconds (set to zero).
ipdstmpl(22) = 1          ! Number of Time Range Specifications fo Field
ipdstmpl(23) = 0          ! Total Number of Data Values Missing
ipdstmpl(24) = 1          ! Statistical Process Used for Field
ipdstmpl(25) = 255        ! Type of Time Increment Between Fields
ipdstmpl(26) = 1          ! Indicator of Unit of Time Range for Processing
ipdstmpl(27) = 6          ! Length of Time Range for Processing
ipdstmpl(28) = 1          ! Indicator of Unit of Time Between Fields
ipdstmpl(29) = 0          ! Time Increment Between Successive Fields
ipdstmplen = 36           ! Max dimension of ipdstmpl()

endif

! Use the g2 lib routine MKIEEE to turn the missing value into an
! IEEE-compliant 32-bit float.
CALL MKIEEE(FLOAT(9999), missVal, 1)

! Set Section 5 Array
idrsnum = 3               ! Data Representation Section Template Number
idrstmpl(1) = 0           ! Reference Value
idrstmpl(2) = 0           ! Binary Scale Factor
if(ipcpvar.eq.0)idrstmpl(3) = 0           ! Decimal Scale Factor
if(ipcpvar.eq.1)idrstmpl(3) = 3           ! Decimal Scale Factor
idrstmpl(4) = 255         ! Number of Bits Used For Group Ref. Value
idrstmpl(5) = 0           ! Type of Original Field Values
idrstmpl(6) = 1           ! Group Splitting Method
idrstmpl(7) = 1           ! Missing Value Management Used
!idrstmpl(8) = 11176255488 ! Primary Missing Value Substitute
idrstmpl(8) = missVal     ! Primary Missing Value Substitute
!idrstmpl(8) = 9999       ! Primary Missing Value Substitute
idrstmpl(9) = 255         ! Secondary Missing Vaule Substitute
idrstmpl(10) = 0          ! Number of Groups Field is Split Into
idrstmpl(11) = 255        ! Reference for Group Widths
idrstmpl(12) = 255        ! Number of Bits Used for Group Widths
idrstmpl(13) = 0          ! Reference for Group Lengths
idrstmpl(14) = 255        ! Length Increment for Group Lenghts
idrstmpl(15) = 0          ! True Length of Last Group
idrstmpl(16) = 255        ! Number of Bits Used for Scaled Group Lengths
idrstmpl(17) = 2          ! Order of Spatial Difference
idrstmpl(18) = 255        ! Number of Octets Required in Data Section for Extra
idrstmplen = 18           ! Max dimension of idrstmpl()

! Set Section 6 Array -- BITMAP INFO; NOT CURRENTLY USED
ibmap = 255

if(ipcpvar.eq.0)then
   ! Data are in decimal format -- inflate to whole percent:
   do i=1,nxa*nya
      if(raw(i).gt.5000.0)then
         fld(i)=9999.0
      else
         fld(i)=raw(i)*100.
      endif
   end do
elseif(ipcpvar.eq.1)then
   do i=1,nxa*nya
      if(raw(i).lt.0.0)then
         !fld(i)=9999.0
         fld(i)=0.0
      else
         fld(i)=raw(i)
      endif
   end do
endif

! Create GRIB2 Message and encode Sec0, Sec1
CALL GRIBCREATE(cgrib,lcgrib,listsec0,listsec1,ierr)

! Encode Sec3
CALL ADDGRID(cgrib,lcgrib,igds,igdstmpl,igdstmplen,ideflist,idefnum,ierr)

! Encode Sec4, Sec5, Sec6, Sec7
CALL ADDFIELD(cgrib,lcgrib,ipdsnum,ipdstmpl,ipdstmplen,coordlist,numcoord,&
             idrsnum,idrstmpl,idrstmplen,fld,ngrdpts,ibmap,bmap,ierr)

! Close GRIB2 Message
CALL GRIBEND(cgrib,lcgrib,lengrib,ierr)

! Write GRIB2 Message to file using WRYTE from the BACIO lib.
CALL WRYTE(KFILGO,lengrib,cgrib)

END SUBROUTINE WRITE_GRIB
