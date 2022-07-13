/*
 * Output lambdaU[i] for BSSN, built using finite-difference derivatives.
 */
void ID_BSSN_lambdas(const paramstruct *restrict params,REAL *restrict xx[3],REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const REAL xx0 = xx[0][i0];
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          const double hDD00_i0m4_i1_i2 = in_gfs[IDX4S(HDD00GF, i0-4,i1,i2)];
          const double hDD00_i0m3_i1_i2 = in_gfs[IDX4S(HDD00GF, i0-3,i1,i2)];
          const double hDD00_i0m2_i1_i2 = in_gfs[IDX4S(HDD00GF, i0-2,i1,i2)];
          const double hDD00_i0m1_i1_i2 = in_gfs[IDX4S(HDD00GF, i0-1,i1,i2)];
          const double hDD00 = in_gfs[IDX4S(HDD00GF, i0,i1,i2)];
          const double hDD00_i0p1_i1_i2 = in_gfs[IDX4S(HDD00GF, i0+1,i1,i2)];
          const double hDD00_i0p2_i1_i2 = in_gfs[IDX4S(HDD00GF, i0+2,i1,i2)];
          const double hDD00_i0p3_i1_i2 = in_gfs[IDX4S(HDD00GF, i0+3,i1,i2)];
          const double hDD00_i0p4_i1_i2 = in_gfs[IDX4S(HDD00GF, i0+4,i1,i2)];
          const double hDD01_i0m4_i1_i2 = in_gfs[IDX4S(HDD01GF, i0-4,i1,i2)];
          const double hDD01_i0m3_i1_i2 = in_gfs[IDX4S(HDD01GF, i0-3,i1,i2)];
          const double hDD01_i0m2_i1_i2 = in_gfs[IDX4S(HDD01GF, i0-2,i1,i2)];
          const double hDD01_i0m1_i1_i2 = in_gfs[IDX4S(HDD01GF, i0-1,i1,i2)];
          const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
          const double hDD01_i0p1_i1_i2 = in_gfs[IDX4S(HDD01GF, i0+1,i1,i2)];
          const double hDD01_i0p2_i1_i2 = in_gfs[IDX4S(HDD01GF, i0+2,i1,i2)];
          const double hDD01_i0p3_i1_i2 = in_gfs[IDX4S(HDD01GF, i0+3,i1,i2)];
          const double hDD01_i0p4_i1_i2 = in_gfs[IDX4S(HDD01GF, i0+4,i1,i2)];
          const double hDD02_i0m4_i1_i2 = in_gfs[IDX4S(HDD02GF, i0-4,i1,i2)];
          const double hDD02_i0m3_i1_i2 = in_gfs[IDX4S(HDD02GF, i0-3,i1,i2)];
          const double hDD02_i0m2_i1_i2 = in_gfs[IDX4S(HDD02GF, i0-2,i1,i2)];
          const double hDD02_i0m1_i1_i2 = in_gfs[IDX4S(HDD02GF, i0-1,i1,i2)];
          const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
          const double hDD02_i0p1_i1_i2 = in_gfs[IDX4S(HDD02GF, i0+1,i1,i2)];
          const double hDD02_i0p2_i1_i2 = in_gfs[IDX4S(HDD02GF, i0+2,i1,i2)];
          const double hDD02_i0p3_i1_i2 = in_gfs[IDX4S(HDD02GF, i0+3,i1,i2)];
          const double hDD02_i0p4_i1_i2 = in_gfs[IDX4S(HDD02GF, i0+4,i1,i2)];
          const double hDD11_i0m4_i1_i2 = in_gfs[IDX4S(HDD11GF, i0-4,i1,i2)];
          const double hDD11_i0m3_i1_i2 = in_gfs[IDX4S(HDD11GF, i0-3,i1,i2)];
          const double hDD11_i0m2_i1_i2 = in_gfs[IDX4S(HDD11GF, i0-2,i1,i2)];
          const double hDD11_i0m1_i1_i2 = in_gfs[IDX4S(HDD11GF, i0-1,i1,i2)];
          const double hDD11 = in_gfs[IDX4S(HDD11GF, i0,i1,i2)];
          const double hDD11_i0p1_i1_i2 = in_gfs[IDX4S(HDD11GF, i0+1,i1,i2)];
          const double hDD11_i0p2_i1_i2 = in_gfs[IDX4S(HDD11GF, i0+2,i1,i2)];
          const double hDD11_i0p3_i1_i2 = in_gfs[IDX4S(HDD11GF, i0+3,i1,i2)];
          const double hDD11_i0p4_i1_i2 = in_gfs[IDX4S(HDD11GF, i0+4,i1,i2)];
          const double hDD12_i0m4_i1_i2 = in_gfs[IDX4S(HDD12GF, i0-4,i1,i2)];
          const double hDD12_i0m3_i1_i2 = in_gfs[IDX4S(HDD12GF, i0-3,i1,i2)];
          const double hDD12_i0m2_i1_i2 = in_gfs[IDX4S(HDD12GF, i0-2,i1,i2)];
          const double hDD12_i0m1_i1_i2 = in_gfs[IDX4S(HDD12GF, i0-1,i1,i2)];
          const double hDD12 = in_gfs[IDX4S(HDD12GF, i0,i1,i2)];
          const double hDD12_i0p1_i1_i2 = in_gfs[IDX4S(HDD12GF, i0+1,i1,i2)];
          const double hDD12_i0p2_i1_i2 = in_gfs[IDX4S(HDD12GF, i0+2,i1,i2)];
          const double hDD12_i0p3_i1_i2 = in_gfs[IDX4S(HDD12GF, i0+3,i1,i2)];
          const double hDD12_i0p4_i1_i2 = in_gfs[IDX4S(HDD12GF, i0+4,i1,i2)];
          const double hDD22_i0m4_i1_i2 = in_gfs[IDX4S(HDD22GF, i0-4,i1,i2)];
          const double hDD22_i0m3_i1_i2 = in_gfs[IDX4S(HDD22GF, i0-3,i1,i2)];
          const double hDD22_i0m2_i1_i2 = in_gfs[IDX4S(HDD22GF, i0-2,i1,i2)];
          const double hDD22_i0m1_i1_i2 = in_gfs[IDX4S(HDD22GF, i0-1,i1,i2)];
          const double hDD22 = in_gfs[IDX4S(HDD22GF, i0,i1,i2)];
          const double hDD22_i0p1_i1_i2 = in_gfs[IDX4S(HDD22GF, i0+1,i1,i2)];
          const double hDD22_i0p2_i1_i2 = in_gfs[IDX4S(HDD22GF, i0+2,i1,i2)];
          const double hDD22_i0p3_i1_i2 = in_gfs[IDX4S(HDD22GF, i0+3,i1,i2)];
          const double hDD22_i0p4_i1_i2 = in_gfs[IDX4S(HDD22GF, i0+4,i1,i2)];
          const double FDPart1_Rational_4_5 = 4.0/5.0;
          const double FDPart1_Rational_4_105 = 4.0/105.0;
          const double FDPart1_Rational_1_5 = 1.0/5.0;
          const double FDPart1_Rational_1_280 = 1.0/280.0;
          const double hDD_dD000 = invdx0*(FDPart1_Rational_1_280*(hDD00_i0m4_i1_i2 - hDD00_i0p4_i1_i2) + FDPart1_Rational_1_5*(hDD00_i0m2_i1_i2 - hDD00_i0p2_i1_i2) + FDPart1_Rational_4_105*(-hDD00_i0m3_i1_i2 + hDD00_i0p3_i1_i2) + FDPart1_Rational_4_5*(-hDD00_i0m1_i1_i2 + hDD00_i0p1_i1_i2));
          const double hDD_dD010 = invdx0*(FDPart1_Rational_1_280*(hDD01_i0m4_i1_i2 - hDD01_i0p4_i1_i2) + FDPart1_Rational_1_5*(hDD01_i0m2_i1_i2 - hDD01_i0p2_i1_i2) + FDPart1_Rational_4_105*(-hDD01_i0m3_i1_i2 + hDD01_i0p3_i1_i2) + FDPart1_Rational_4_5*(-hDD01_i0m1_i1_i2 + hDD01_i0p1_i1_i2));
          const double hDD_dD020 = invdx0*(FDPart1_Rational_1_280*(hDD02_i0m4_i1_i2 - hDD02_i0p4_i1_i2) + FDPart1_Rational_1_5*(hDD02_i0m2_i1_i2 - hDD02_i0p2_i1_i2) + FDPart1_Rational_4_105*(-hDD02_i0m3_i1_i2 + hDD02_i0p3_i1_i2) + FDPart1_Rational_4_5*(-hDD02_i0m1_i1_i2 + hDD02_i0p1_i1_i2));
          const double hDD_dD110 = invdx0*(FDPart1_Rational_1_280*(hDD11_i0m4_i1_i2 - hDD11_i0p4_i1_i2) + FDPart1_Rational_1_5*(hDD11_i0m2_i1_i2 - hDD11_i0p2_i1_i2) + FDPart1_Rational_4_105*(-hDD11_i0m3_i1_i2 + hDD11_i0p3_i1_i2) + FDPart1_Rational_4_5*(-hDD11_i0m1_i1_i2 + hDD11_i0p1_i1_i2));
          const double hDD_dD120 = invdx0*(FDPart1_Rational_1_280*(hDD12_i0m4_i1_i2 - hDD12_i0p4_i1_i2) + FDPart1_Rational_1_5*(hDD12_i0m2_i1_i2 - hDD12_i0p2_i1_i2) + FDPart1_Rational_4_105*(-hDD12_i0m3_i1_i2 + hDD12_i0p3_i1_i2) + FDPart1_Rational_4_5*(-hDD12_i0m1_i1_i2 + hDD12_i0p1_i1_i2));
          const double hDD_dD220 = invdx0*(FDPart1_Rational_1_280*(hDD22_i0m4_i1_i2 - hDD22_i0p4_i1_i2) + FDPart1_Rational_1_5*(hDD22_i0m2_i1_i2 - hDD22_i0p2_i1_i2) + FDPart1_Rational_4_105*(-hDD22_i0m3_i1_i2 + hDD22_i0p3_i1_i2) + FDPart1_Rational_4_5*(-hDD22_i0m1_i1_i2 + hDD22_i0p1_i1_i2));
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          const double FDPart3_0 = sin(xx1);
          const double FDPart3_1 = hDD12*((xx0)*(xx0)*(xx0));
          const double FDPart3_2 = ((xx0)*(xx0));
          const double FDPart3_3 = FDPart3_2*hDD11 + FDPart3_2;
          const double FDPart3_4 = FDPart3_0*hDD02;
          const double FDPart3_5 = FDPart3_0*FDPart3_1*hDD01 - FDPart3_3*FDPart3_4*xx0;
          const double FDPart3_7 = ((FDPart3_0)*(FDPart3_0));
          const double FDPart3_8 = FDPart3_7*((xx0)*(xx0)*(xx0)*(xx0));
          const double FDPart3_9 = hDD00 + 1;
          const double FDPart3_10 = FDPart3_8*((hDD12)*(hDD12));
          const double FDPart3_11 = FDPart3_2*FDPart3_7;
          const double FDPart3_12 = FDPart3_11*((hDD02)*(hDD02));
          const double FDPart3_13 = FDPart3_11*hDD22 + FDPart3_11;
          const double FDPart3_14 = FDPart3_2*((hDD01)*(hDD01));
          const double FDPart3_16 = (1.0/(-FDPart3_10*FDPart3_9 - FDPart3_12*FDPart3_3 - FDPart3_13*FDPart3_14 + FDPart3_13*FDPart3_3*FDPart3_9 + 2*FDPart3_8*hDD01*hDD02*hDD12));
          const double FDPart3_17 = cos(xx1);
          const double FDPart3_18 = FDPart3_17*FDPart3_2*hDD12;
          const double FDPart3_20 = 2*xx0;
          const double FDPart3_23 = -FDPart3_2*hDD_dD110 - FDPart3_20*hDD11 - FDPart3_20;
          const double FDPart3_24 = FDPart3_16*(-FDPart3_10 + FDPart3_13*FDPart3_3);
          const double FDPart3_25 = (1.0/2.0)*FDPart3_24;
          const double FDPart3_26 = FDPart3_16*(-FDPart3_12 + FDPart3_13*FDPart3_9);
          const double FDPart3_27 = FDPart3_17*hDD02*xx0;
          const double FDPart3_28 = FDPart3_0*FDPart3_2;
          const double FDPart3_31 = FDPart3_0*FDPart3_20*hDD12;
          const double FDPart3_32 = FDPart3_28*hDD_dD120 + FDPart3_31;
          const double FDPart3_33 = FDPart3_27 + FDPart3_32;
          const double FDPart3_34 = (1.0/2.0)*FDPart3_16;
          const double FDPart3_35 = FDPart3_34*FDPart3_5;
          const double FDPart3_36 = FDPart3_2*hDD_dD110 + FDPart3_20*hDD11 + FDPart3_20;
          const double FDPart3_37 = FDPart3_1*FDPart3_7*hDD02 - FDPart3_13*hDD01*xx0;
          const double FDPart3_38 = FDPart3_34*FDPart3_37;
          const double FDPart3_39 = 2*FDPart3_16;
          const double FDPart3_40 = FDPart3_37*FDPart3_39;
          const double FDPart3_41 = sin(2*xx1);
          const double FDPart3_43 = 2*FDPart3_17*FDPart3_28*hDD22;
          const double FDPart3_44 = FDPart3_2*FDPart3_41 + FDPart3_43;
          const double FDPart3_45 = FDPart3_27 - FDPart3_28*hDD_dD120 - FDPart3_31;
          const double FDPart3_46 = FDPart3_2*FDPart3_4*hDD01 - FDPart3_28*FDPart3_9*hDD12;
          const double FDPart3_47 = FDPart3_39*FDPart3_46;
          const double FDPart3_49 = 2*FDPart3_7*xx0;
          const double FDPart3_52 = FDPart3_11*hDD_dD220 + FDPart3_49*hDD22 + FDPart3_49;
          const double FDPart3_53 = -FDPart3_27 + FDPart3_32;
          const double FDPart3_54 = FDPart3_39*FDPart3_5;
          const double FDPart3_55 = -FDPart3_2*FDPart3_41 - FDPart3_43;
          const double FDPart3_56 = -FDPart3_11*hDD_dD220 - FDPart3_49*hDD22 - FDPart3_49;
          const double FDPart3_57 = FDPart3_16*(-FDPart3_14 + FDPart3_3*FDPart3_9);
          const double FDPart3_58 = FDPart3_0*FDPart3_20*hDD_dD020 + 2*FDPart3_4;
          const double FDPart3_59 = FDPart3_20*hDD_dD010 + 2*hDD01;
          const double FDPart3_60 = FDPart3_34*FDPart3_46;
          const double FDPart3_61 = (1.0/2.0)*FDPart3_26;
          const double FDPart3_62 = -1/xx0;
          const double FDPart3_64 = (1.0/2.0)*FDPart3_57;
          in_gfs[IDX4S(LAMBDAU0GF, i0, i1, i2)] = FDPart3_24*(FDPart3_25*hDD_dD000 + FDPart3_35*FDPart3_58 + FDPart3_38*FDPart3_59) + FDPart3_26*(FDPart3_16*FDPart3_18*FDPart3_5 + FDPart3_23*FDPart3_25 + xx0) + FDPart3_40*(FDPart3_33*FDPart3_35 + FDPart3_36*FDPart3_38) + FDPart3_47*(FDPart3_25*FDPart3_45 + FDPart3_35*FDPart3_44) + FDPart3_54*(FDPart3_35*FDPart3_52 + FDPart3_38*FDPart3_53) + FDPart3_57*(FDPart3_25*FDPart3_56 + FDPart3_38*FDPart3_55 + FDPart3_7*xx0);
          in_gfs[IDX4S(LAMBDAU1GF, i0, i1, i2)] = xx0*(FDPart3_24*(FDPart3_38*hDD_dD000 + FDPart3_58*FDPart3_60 + FDPart3_59*FDPart3_61) + FDPart3_26*(FDPart3_16*FDPart3_18*FDPart3_46 + FDPart3_23*FDPart3_38) + FDPart3_40*(FDPart3_33*FDPart3_60 + FDPart3_36*FDPart3_61 + FDPart3_62) + FDPart3_47*(FDPart3_38*FDPart3_45 + FDPart3_44*FDPart3_60) + FDPart3_54*(FDPart3_52*FDPart3_60 + FDPart3_53*FDPart3_61) + FDPart3_57*(FDPart3_38*FDPart3_56 + (1.0/2.0)*FDPart3_41 + FDPart3_55*FDPart3_61));
          in_gfs[IDX4S(LAMBDAU2GF, i0, i1, i2)] = FDPart3_0*xx0*(FDPart3_24*(FDPart3_35*hDD_dD000 + FDPart3_58*FDPart3_64 + FDPart3_59*FDPart3_60) + FDPart3_26*(FDPart3_18*FDPart3_57 + FDPart3_23*FDPart3_35) + FDPart3_40*(FDPart3_33*FDPart3_64 + FDPart3_36*FDPart3_60) + FDPart3_47*(FDPart3_35*FDPart3_45 - 1.0/2.0*FDPart3_41/FDPart3_7 + FDPart3_44*FDPart3_64) + FDPart3_54*(FDPart3_52*FDPart3_64 + FDPart3_53*FDPart3_60 + FDPart3_62) + FDPart3_57*(FDPart3_35*FDPart3_56 + FDPart3_55*FDPart3_60));
        
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
