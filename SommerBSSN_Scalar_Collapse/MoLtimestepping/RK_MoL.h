// C code implementation of RK4 Method of Lines timestepping.
// ***k1 substep:***

Ricci_eval(&rfmstruct, &params, y_n_gfs, auxevol_gfs);
rhs_eval(&rfmstruct, &params, auxevol_gfs, y_n_gfs, k_odd_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = k_odd_gfs[i]*dt*(1.0/6.0);
  k_odd_gfs[i] = y_n_gfs[i] + k_odd_gfs[i]*dt*(1.0/2.0);
}

//apply_bcs_curvilinear(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_odd_gfs);
apply_bcs_sommerfeld(&params, &xx[3], &bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_odd_gfs, &rhs_gfs); 
enforce_detgammahat_constraint(&rfmstruct, &params,                     k_odd_gfs);


// ***k2 substep:***

Ricci_eval(&rfmstruct, &params, k_odd_gfs, auxevol_gfs);
rhs_eval(&rfmstruct, &params, auxevol_gfs, k_odd_gfs, k_even_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = y_nplus1_running_total_gfs[i] + k_even_gfs[i]*dt*(1.0/3.0);
  k_even_gfs[i] = y_n_gfs[i] + k_even_gfs[i]*dt*(1.0/2.0);
}

//apply_bcs_curvilinear(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_even_gfs);
apply_bcs_sommerfeld(&params, &xx[3], &bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_even_gfs, &rhs_gfs); 
enforce_detgammahat_constraint(&rfmstruct, &params,                     k_even_gfs);


// ***k3 substep:***

Ricci_eval(&rfmstruct, &params, k_even_gfs, auxevol_gfs);
rhs_eval(&rfmstruct, &params, auxevol_gfs, k_even_gfs, k_odd_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = y_nplus1_running_total_gfs[i] + k_odd_gfs[i]*dt*(1.0/3.0);
  k_odd_gfs[i] = y_n_gfs[i] + k_odd_gfs[i]*dt;
}

//apply_bcs_curvilinear(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_odd_gfs);
apply_bcs_sommerfeld(&params, &xx[3], &bcstruct, NUM_EVOL_GFS, evol_gf_parity, k_odd_gfs, &rhs_gfs); 
enforce_detgammahat_constraint(&rfmstruct, &params,                     k_odd_gfs);


// ***k4 substep:***

Ricci_eval(&rfmstruct, &params, k_odd_gfs, auxevol_gfs);
rhs_eval(&rfmstruct, &params, auxevol_gfs, k_odd_gfs, k_even_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_n_gfs[i] = y_n_gfs[i] + y_nplus1_running_total_gfs[i] + k_even_gfs[i]*dt*(1.0/6.0);
}

//apply_bcs_curvilinear(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs);
apply_bcs_sommerfeld(&params, &xx[3], &bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs, &rhs_gfs); 
enforce_detgammahat_constraint(&rfmstruct, &params,                     y_n_gfs);


