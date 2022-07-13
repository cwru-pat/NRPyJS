/*
 * Driver function for ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(),
 * which writes BSSN variables in terms of ADM variables at a given point xx0,xx1,xx2
 */
void ID_BSSN__ALL_BUT_LAMBDAs(const paramstruct *restrict params,REAL *restrict xx[3],ID_inputs other_inputs,REAL *in_gfs) {
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
        ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(params, xx0xx1xx2,other_inputs,
                            &in_gfs[IDX4ptS(HDD00GF,idx)],&in_gfs[IDX4ptS(HDD01GF,idx)],&in_gfs[IDX4ptS(HDD02GF,idx)],
                            &in_gfs[IDX4ptS(HDD11GF,idx)],&in_gfs[IDX4ptS(HDD12GF,idx)],&in_gfs[IDX4ptS(HDD22GF,idx)],
                            &in_gfs[IDX4ptS(ADD00GF,idx)],&in_gfs[IDX4ptS(ADD01GF,idx)],&in_gfs[IDX4ptS(ADD02GF,idx)],
                            &in_gfs[IDX4ptS(ADD11GF,idx)],&in_gfs[IDX4ptS(ADD12GF,idx)],&in_gfs[IDX4ptS(ADD22GF,idx)],
                            &in_gfs[IDX4ptS(TRKGF,idx)],
                            &in_gfs[IDX4ptS(VETU0GF,idx)],&in_gfs[IDX4ptS(VETU1GF,idx)],&in_gfs[IDX4ptS(VETU2GF,idx)],
                            &in_gfs[IDX4ptS(BETU0GF,idx)],&in_gfs[IDX4ptS(BETU1GF,idx)],&in_gfs[IDX4ptS(BETU2GF,idx)],
                            &in_gfs[IDX4ptS(ALPHAGF,idx)],&in_gfs[IDX4ptS(CFGF,idx)]);
        
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
