SUBROUTINE read_dressing_parameters (n_amounts, infile, ramt, &
    fraczero, gamma_shape, gamma_scale)
    
    USE netcdf
    
    INTEGER, INTENT(IN) :: n_amounts
    CHARACTER*(*), INTENT(IN) :: infile
    REAL, INTENT(OUT), DIMENSION(n_amounts) :: ramt ! amts associated
    ! with indices of fraczero, gamma_shape, gamma_scale
    REAL, INTENT(OUT), DIMENSION(n_amounts) :: fraczero ! fraction
    ! of best-member dressing samples with zero precip
    REAL, INTENT(OUT), DIMENSION(n_amounts) :: gamma_shape ! alpha 
    ! parameter in the fitted Gamma distribution
    REAL, INTENT(OUT), DIMENSION(n_amounts) :: gamma_scale ! beta
    ! parameter in the fitted Gamma distribution                

    PRINT *,TRIM(infile)
    CALL check(nf90_open(TRIM(infile), NF90_NOWRITE, netid))
    CALL check(nf90_inq_varid(netid, "precip_values", ivar))
    CALL check(nf90_get_var(netid, ivar, ramt, &
        start=(/1/), count=(/n_amounts/)))
        
    CALL check(nf90_inq_varid(netid, "fraction_zeros", ivar))
    CALL check(nf90_get_var(netid, ivar, fraczero, &
        start=(/1/), count=(/n_amounts/)))
            
    CALL check(nf90_inq_varid(netid, "gamma_shapes", ivar))
    CALL check(nf90_get_var(netid, ivar, gamma_shape, &
        start=(/1/), count=(/n_amounts/)))
    
    CALL check(nf90_inq_varid(netid, "gamma_scales", ivar))
    CALL check(nf90_get_var(netid, ivar, gamma_scale, &
        start=(/1/), count=(/n_amounts/)))
        
    CALL check(nf90_close(netid))
    
RETURN
END SUBROUTINE read_dressing_parameters
