SUBROUTINE rejigger_weights(nmembers, nens, closest_histogram, member_weights)
    !
    ! purpose: this subroutine is called in the eventuality that not all member
    !    forecasts that were expected were actually available.  In this 
    !    instance, then set the member_weights output array such that 
    !    it reflects the same general shape as closest_histogram but
    !    has nens rather than nmember elements.  Coded by Tom Hamill,
    !    Dec 2016.
    INTEGER, INTENT(IN) :: nmembers ! dimension of input closest_histogram
    INTEGER, INTENT(IN) :: nens ! dimension of output member_weights
    REAL, INTENT(IN), DIMENSION(nmembers) :: closest_histogram
    REAL, INTENT(OUT), DIMENSION(nens) :: member_weights
    
    REAL :: rstart  ! the starting coordinate of the boundary 
    !                 for summing information from closest_histogram
    REAL :: rend    ! ending coordinate of the boundary 
    !                 for summing information from closest_histogram
    
    !print *,'nmembers, nens = ', nmembers, nens
    !print *, 'closest_histogram = ', closest_histogram
    
    slope = REAL(nmembers) / REAL(nens) 
    !print *,'slope = ', slope
    member_weights(:) = 0.0
    DO iens = 1, nens  
        rstart = REAL(iens-1)*slope ! lower integration bound in closest_histogram [0,nmembers]
        rend = REAL(iens)*slope ! upper integration bound of closest_histogram [0,nmembers]
        iboundary_start = 1 + INT(rstart) ! start, end indices in 
        iboundary_fin = 1 + INT(rend) ! closest_histogram array
        IF (iboundary_fin .gt. nmembers) iboundary_fin = nmembers
        !PRINT *,'iens, rstart, rend, iboundary_start, iboundary_fin = ', &
        !    iens, rstart, rend, iboundary_start, iboundary_fin
        IF (iboundary_start .eq. iboundary_fin) THEN
            member_weights(iens) = (rend-rstart)*closest_histogram(iboundary_start)
        ELSE
            DO ib = iboundary_start, iboundary_fin
                IF (ib .eq. iboundary_start) THEN
                    rbound_up = REAL(ib)
                    rbound_dn = rstart
                ELSE IF (ib .eq. iboundary_fin) THEN
                    rbound_up = rend
                    rbound_dn = REAL(ib-1)
                ELSE
                    rbound_up = REAL(ib)
                    rbound_dn = REAL(ib-1)
                ENDIF   
                !print *,'   ib, rbound_dn, rbound_up = ', ib, rbound_dn, rbound_up                
                member_weights(iens) = member_weights(iens) + &
                    (rbound_up-rbound_dn)*closest_histogram(ib)
            END DO 
        END IF
        !print *,'iens, member_weights(iens) = ',iens,member_weights(iens)
    ENDDO
    !PRINT *,'closest_histogram = ', closest_histogram
    !PRINT *,'member_weights = ', member_weights
    RETURN
END SUBROUTINE rejigger_weights
