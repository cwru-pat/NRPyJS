/*
 * (c) 2021 Leo Werneck
 * This function takes as input either (x,y,z) or (r,th,ph) and outputs all
 * scalar field quantities in the Cartesian or Spherical basis, respectively.
 */
void ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(const paramstruct *restrict params,const REAL xx0xx1xx2[3],
ID_inputs other_inputs,
REAL *restrict sf, REAL *restrict sfM) {
#include "./set_Cparameters.h"


  const REAL xx0 = xx0xx1xx2[0];
  const REAL xx1 = xx0xx1xx2[1];
  const REAL xx2 = xx0xx1xx2[2];
  REAL rthph[3];
  rthph[0] = xx0;
  rthph[1] = xx1;
  rthph[2] = xx2;

  ID_scalarfield_spherical(rthph,other_inputs,sf,sfM);
}
