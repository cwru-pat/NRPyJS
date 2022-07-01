/*
 * Enforce det(gammabar) = det(gammahat) constraint.
 */
void enforce_detgammahat_constraint(const rfm_struct *restrict rfmstruct, const paramstruct *restrict params, REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    #include "rfm_files/rfm_struct__read2.h"
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      #include "rfm_files/rfm_struct__read1.h"
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        #include "rfm_files/rfm_struct__read0.h"
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          const double hDD00 = in_gfs[IDX4S(HDD00GF, i0,i1,i2)];
          const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
          const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
          const double hDD11 = in_gfs[IDX4S(HDD11GF, i0,i1,i2)];
          const double hDD12 = in_gfs[IDX4S(HDD12GF, i0,i1,i2)];
          const double hDD22 = in_gfs[IDX4S(HDD22GF, i0,i1,i2)];
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          const double FDPart3_0 = hDD00 + 1;
          const double FDPart3_2 = ((f0_of_xx0)*(f0_of_xx0)*(f0_of_xx0)*(f0_of_xx0))*((f1_of_xx1)*(f1_of_xx1));
          const double FDPart3_3 = ((f0_of_xx0)*(f0_of_xx0));
          const double FDPart3_4 = FDPart3_3*hDD11 + FDPart3_3;
          const double FDPart3_5 = FDPart3_3*((f1_of_xx1)*(f1_of_xx1));
          const double FDPart3_6 = FDPart3_5*hDD22 + FDPart3_5;
          const double FDPart3_7 = cbrt(fabs(FDPart3_2)/(-FDPart3_0*FDPart3_2*((hDD12)*(hDD12)) + FDPart3_0*FDPart3_4*FDPart3_6 + 2*FDPart3_2*hDD01*hDD02*hDD12 - FDPart3_3*FDPart3_6*((hDD01)*(hDD01)) - FDPart3_4*FDPart3_5*((hDD02)*(hDD02))));
          in_gfs[IDX4S(HDD00GF, i0, i1, i2)] = FDPart3_0*FDPart3_7 - 1;
          in_gfs[IDX4S(HDD01GF, i0, i1, i2)] = FDPart3_7*hDD01;
          in_gfs[IDX4S(HDD02GF, i0, i1, i2)] = FDPart3_7*hDD02;
          in_gfs[IDX4S(HDD11GF, i0, i1, i2)] = FDPart3_7*(hDD11 + 1) - 1;
          in_gfs[IDX4S(HDD12GF, i0, i1, i2)] = FDPart3_7*hDD12;
          in_gfs[IDX4S(HDD22GF, i0, i1, i2)] = FDPart3_7*(hDD22 + 1) - 1;
        
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
