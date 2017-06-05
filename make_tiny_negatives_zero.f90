SUBROUTINE make_tiny_negatives_zero(nxf,nyf,pfcst)
INTEGER, INTENT(IN) :: nxf, nyf
REAL, INTENT(INOUT), DIMENSION(nxf,nyf) :: pfcst
DO i = 1,nxf
   DO j = 1,nyf
      IF (pfcst(i,j) .lt. 0. .and. pfcst(i,j) .gt. -1.) pfcst(i,j) = 0.0
   END DO
END DO
RETURN
END SUBROUTINE make_tiny_negatives_zero
