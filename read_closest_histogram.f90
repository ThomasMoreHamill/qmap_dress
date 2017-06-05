SUBROUTINE read_closest_histogram(nmembers, histofile, closest_histogram)

    USE netcdf
    INTEGER, INTENT(IN) :: nmembers
    CHARACTER*(*), INTENT(IN) :: histofile
    REAL, INTENT(OUT), DIMENSION(nmembers) :: closest_histogram


    PRINT *,TRIM(histofile)
    CALL check(nf90_open(TRIM(histofile), NF90_NOWRITE, netid))
    CALL check(nf90_inq_varid(netid, "closest_histogram", ivar))
    CALL check(nf90_get_var(netid, ivar, closest_histogram, &
        start=(/1/), count=(/nmembers/)))
    CALL check(nf90_close(netid))

    RETURN

END SUBROUTINE read_closest_histogram
	
