/*
 * (c) 2021 Leo Werneck
 * This is the scalar field initial data driver functiono.
 */
void ID_scalarfield(const paramstruct *restrict params,REAL *restrict xx[3],
                ID_inputs other_inputs,REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        const REAL xx0 = xx[0][i0];
        
          const int idx = IDX3S(i0,i1,i2);
          const REAL xx0xx1xx2[3] = {xx0,xx1,xx2};
          ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(params,xx0xx1xx2,other_inputs,
                                                     &in_gfs[IDX4ptS(SFGF,idx)],
                                                     &in_gfs[IDX4ptS(SFMGF,idx)]);
        
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
