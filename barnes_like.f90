SUBROUTINE barnes_like(ihavedata, precip_in, nxa, nya, ideterministic, precip_out)

! ---- perform a Barnes-like objective analysis of points over water, i.e., a weighted sum
!      of the values of nearby land points, with the weights being an exponential function of distance.

INTEGER, INTENT(IN) :: nxa, nya, ideterministic
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: ihavedata
REAL, INTENT(INOUT), DIMENSION(nxa,nya) :: precip_in
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: precip_out

!f2py intent(in) nxa, nya, ihavedata, precip_in
!f2py intent(out) precip_out
!f2py depend(nxa,nya) precip_in, precip_out

REAL*8, DIMENSION(nxa,nya) :: numer, denom, work
INTEGER, DIMENSION(3) :: cutradius
REAL*8 weight

pavg = SUM(precip_in*REAL(ihavedata))/SUM(real(ihavedata))

! ---- now do something akin to a Barnes filtering, weighting the data by exp(-dist**2/lengthscale**2)
!      this loop tallies up the numerator and denominator needed for weighted sum calculation

!PRINT *,'Barnes_like: ihavedata(42:82,93) = ',ihavedata(42:82,93)

icut = 10
rlenscale = 3.
rlenscale2 = rlenscale**2
numer = 0.
denom = 0.
DO ixa = 1, nxa
   imin = MAX(1,ixa-icut)
   imax = MIN(nxa,ixa+icut)
   DO jya = 1, nya
      jmin = MAX(1,jya-icut)
      jmax = MIN(nya,jya+icut)
      IF (ihavedata(ixa,jya) .eq. 0) THEN
         DO iloop = imin, imax
            DO jloop = jmin, jmax
               dist = SQRT(REAL((ixa-iloop)**2 + (jya-jloop)**2))
               weight = exp(-dist**2/rlenscale2)


               !IF (ixa .eq. 62 .and. jya .eq. 93) THEN
               !   PRINT 2056,'ixa,jya,iloop,jloop,ihavedata,dist,wt,precip=',&
               !     ixa,jya,iloop,jloop,ihavedata(iloop,jloop),dist,weight, precip_in(iloop,jloop)
               !   2056 FORMAT(a43,5(i3,1x),3(f9.6,1x))
               !ENDIF

               IF (dist .le. icut .and. ihavedata(iloop,jloop) .eq. 1 .and. precip_in(iloop,jloop) .ge. 0.0) THEN
                  numer(ixa,jya) = numer(ixa,jya)  + weight*precip_in(iloop,jloop)
                  denom(ixa,jya) = denom(ixa,jya)  + weight

                  !IF (ixa .eq. 62 .and. jya .eq. 93) THEN
                  !   PRINT *,'iloop,jloop,dist,wt,precip=',&
                  !     ixa,jya,iloop,jloop,dist,weight, precip_in(iloop,jloop)
                  !   PRINT *,'  numer, denom = ', numer(ixa,jya), denom(ixa,jya)   
                  !ENDIF

               ENDIF
            END DO
         END DO
      ENDIF
   END DO ! jya
END DO ! ixa

!PRINT *,'precip_out(62,93) = ',numer(62,93) /denom (62,93)

! --- set the output over water to the weighted sum.  over land use original values.  If
!     at a point that is too far from land to have any points counted in weight, just set the 
!     value to the domain-averaged land probability.

DO ixa = 1, nxa
   DO jya = 1, nya
      IF (denom(ixa,jya) .gt. 0.0) THEN 
         precip_out(ixa,jya) = numer(ixa, jya) /denom (ixa, jya) 
      ELSE
         precip_out(ixa,jya) = precip_in(ixa,jya)
      ENDIF
      IF (ideterministic .eq. 0) THEN !probabilistic data, 0 to 1
         IF (precip_out(ixa,jya) .lt. 0.0 .or. precip_out(ixa,jya) .gt. 1.0) precip_out(ixa,jya) = pavg
      ELSE ! deterministic data >=0
         IF (precip_out(ixa,jya) .lt. 0.0) precip_out(ixa,jya) = pavg
      ENDIF
   END DO
END DO

RETURN
END SUBROUTINE barnes_like
