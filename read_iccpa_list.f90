SUBROUTINE read_iccpa_list (nxa, nya, nxf, nyf, nxs, nys, ncount, nnearest, &
   ilist, jlist, rlonsf, rlatsf, rlonsa, rlatsa)

! purpose:  previously, a data set of which 1/8-degree grid points were nearest
!    to each polar stereo grid point was calculated.  Read this in.

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf, ncount
INTEGER, INTENT(OUT), DIMENSION(nxf,nyf) :: nnearest ! # grid pts nearest to this polar stereo point
INTEGER, INTENT(OUT), DIMENSION(nxf,nyf,ncount) :: ilist ! list of i coordinates nearest
INTEGER, INTENT(OUT), DIMENSION(nxf,nyf,ncount) :: jlist ! list of j coordinates nearest
REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: rlonsf, rlatsf ! lon/lat of forecast polar-stereo grid
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlonsa, rlatsa ! lon/lat of 1/8-degree analysis
INTEGER :: xoffset, yoffset, i, j, g, h, k
INTEGER, DIMENSION(nxf,nyf,ncount) :: ilist_old ! list of i coordinates nearest
INTEGER, DIMENSION(nxf,nyf,ncount) :: jlist_old ! list of j coordinates nearest
CHARACTER*200 :: JUNK, JUNK2
!REAL, DIMENSION(nxs,nys) :: rlatsa_s,rlonsa_s

! ---- read the data from the file, previously calculated by subroutine find_nearest_polarst
!      (in the file fortran_routine.f90) called from python routine find_nearest_latlons.py

OPEN(unit=22,status='old',form='unformatted')
READ (22) nxfin, nyfin, nxain, nyain, ncountin
IF (nxfin .ne. nxf .or. nyfin .ne. nyf .or. nxain .ne. nxa .or. nyain .ne. nya .or. &
     ncountin .ne. ncount) THEN
   PRINT *,'error in read_iccpa_list: array dimensions not right.'
   PRINT *,'nxfin, nyfin, nxain, nyain, ncountin = ',nxfin, nyfin, nxain, nyain, ncountin 
   PRINT *,'nxf, nyf, nxa, nya, ncount = ',nxf, nyf, nxa, nya, ncount
   PRINT *,'Stopping.'
   STOP
ENDIF 
READ (22) nnearest
READ (22) ilist
READ (22) jlist
READ (22) rlonsf
READ (22) rlatsf
!READ (22) rlonsa
!READ (22) rlatsa
!READ (22) conusmask
CLOSE (22)

!DO i = 1,nxf
!   DO j = 1,nyf
!      DO k = 1,ncount
!         ilist(i,j,k) = ilist_old(i,j,k) + 21
!         jlist(i,j,k) = jlist_old(i,j,k) + 22
!      END DO
!   END DO
!END DO




rlatsa(:,:) = -99.99
rlonsa(:,:) = -99.99

OPEN(unit=23,status='old',form='formatted')
DO i = 1, nxa
   DO j = 1, nya
        READ(23,1010) rlatsa(i,j), rlonsa(i,j)
1010    FORMAT (56X,F7.4,2X,F8.4,73X)
!        PRINT *,'Lat, Lon, I, J = ', rlatsa(i,j), rlonsa(i,j), i, j
   END DO
END DO
CLOSE(23)

!PRINT *, 'rlatsa = ', rlatsa(1,:)
!PRINT *, 'rlonsa = ', rlonsa(:,1)

RETURN
END SUBROUTINE read_iccpa_list
