/*
 * (c) 2021 Leo Werneck
 * This function takes as input either (x,y,z) or (r,th,ph) and outputs all
 * scalar field quantities in the Cartesian or Spherical basis, respectively.
 */
void ID_scalarfield_spherical(const REAL xyz_or_rthph[3],const ID_inputs other_inputs,REAL *restrict sf,REAL *restrict sfM) {


  const REAL r  = xyz_or_rthph[0];
  const REAL th = xyz_or_rthph[1];
  const REAL ph = xyz_or_rthph[2];

  REAL sf_star,psi4_star,alpha_star;

  scalarfield_interpolate_1D(r,
                             other_inputs.interp_stencil_size,
                             other_inputs.numlines_in_file,
                             other_inputs.r_arr,
                             other_inputs.sf_arr,
                             other_inputs.psi4_arr,
                             other_inputs.alpha_arr,
                             &sf_star,&psi4_star,&alpha_star);

  // Update varphi
  *sf  = sf_star;
  // Update Pi
  *sfM = 0;
}
