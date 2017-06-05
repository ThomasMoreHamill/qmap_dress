SUBROUTINE read_dressing_parameters_v2 (n_amounts, n_climocats, infile, ramt, &
    fraczero, fraczero_fclimpop, climo_pop_thresholds, gamma_shape, gamma_scale)

    ! read statistics that are used to perform best-member ensemble dressing

    USE netcdf
    
    INTEGER, INTENT(IN) :: n_amounts, n_climocats
    CHARACTER*(*), INTENT(IN) :: infile
    REAL, INTENT(OUT), DIMENSION(n_amounts) :: ramt ! amts associated
    ! with indices of fraczero, gamma_shape, gamma_scale
    REAL, INTENT(OUT), DIMENSION(n_amounts) :: fraczero ! fraction
    ! of best-member dressing samples with zero precip
    REAL, INTENT(OUT), DIMENSION(n_climocats) :: fraczero_fclimpop ! fraction
    ! of best-member dressing samples with zero precip for situation when
    ! forecast = 0, stratified by bins of climatological probability
    REAL, INTENT(OUT), DIMENSION(n_climocats-1) :: climo_pop_thresholds
    ! these are the climatological POP boundaries that serve as boundaries
    ! between the categories in the fraczero_fclimpop vector

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
            
    CALL check(nf90_inq_varid(netid, "fraction_zeros_fclimpop", ivar))
    CALL check(nf90_get_var(netid, ivar, fraczero_fclimpop, &
        start=(/1/), count=(/n_climocats/)))
            
    CALL check(nf90_inq_varid(netid, "climo_pop_thresholds", ivar))
    CALL check(nf90_get_var(netid, ivar, climo_pop_thresholds, &
         start=(/1/), count=(/n_climocats-1/)))

    CALL check(nf90_inq_varid(netid, "gamma_shapes", ivar))
    CALL check(nf90_get_var(netid, ivar, gamma_shape, &
        start=(/1/), count=(/n_amounts/)))
    
    CALL check(nf90_inq_varid(netid, "gamma_scales", ivar))
    CALL check(nf90_get_var(netid, ivar, gamma_scale, &
        start=(/1/), count=(/n_amounts/)))
        
    CALL check(nf90_close(netid))
    
RETURN
END SUBROUTINE read_dressing_parameters_v2
