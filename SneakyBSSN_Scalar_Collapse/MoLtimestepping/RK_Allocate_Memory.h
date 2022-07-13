// Code snippet allocating gridfunction memory for "RK4" method:
REAL *restrict y_n_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
REAL *restrict y_nplus1_running_total_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
REAL *restrict k_odd_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
REAL *restrict k_even_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
REAL *restrict diagnostic_output_gfs = y_nplus1_running_total_gfs;
