/*
 * Write BSSN variables in terms of ADM variables at a given point xx0,xx1,xx2
 */
void ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(const paramstruct *restrict params, const REAL xx0xx1xx2[3],ID_inputs other_inputs,
                    REAL *hDD00,REAL *hDD01,REAL *hDD02,REAL *hDD11,REAL *hDD12,REAL *hDD22,
                    REAL *aDD00,REAL *aDD01,REAL *aDD02,REAL *aDD11,REAL *aDD12,REAL *aDD22,
                    REAL *trK,
                    REAL *vetU0,REAL *vetU1,REAL *vetU2,
                    REAL *betU0,REAL *betU1,REAL *betU2,
                    REAL *alpha,  REAL *cf) {
#include "./set_Cparameters.h"


      REAL gammaSphorCartDD00,gammaSphorCartDD01,gammaSphorCartDD02,
           gammaSphorCartDD11,gammaSphorCartDD12,gammaSphorCartDD22;
      REAL KSphorCartDD00,KSphorCartDD01,KSphorCartDD02,
           KSphorCartDD11,KSphorCartDD12,KSphorCartDD22;
      REAL alphaSphorCart,betaSphorCartU0,betaSphorCartU1,betaSphorCartU2;
      REAL BSphorCartU0,BSphorCartU1,BSphorCartU2;
      const REAL xx0 = xx0xx1xx2[0];
      const REAL xx1 = xx0xx1xx2[1];
      const REAL xx2 = xx0xx1xx2[2];
      REAL xyz_or_rthph[3];
  xyz_or_rthph[0] = xx0;
  xyz_or_rthph[1] = xx1;
  xyz_or_rthph[2] = xx2;
      ID_scalarfield_ADM_quantities(xyz_or_rthph, other_inputs,
                       &gammaSphorCartDD00,&gammaSphorCartDD01,&gammaSphorCartDD02,
                       &gammaSphorCartDD11,&gammaSphorCartDD12,&gammaSphorCartDD22,
                       &KSphorCartDD00,&KSphorCartDD01,&KSphorCartDD02,
                       &KSphorCartDD11,&KSphorCartDD12,&KSphorCartDD22,
                       &alphaSphorCart,&betaSphorCartU0,&betaSphorCartU1,&betaSphorCartU2,
                       &BSphorCartU0,&BSphorCartU1,&BSphorCartU2);
      // Next compute all rescaled BSSN curvilinear quantities:
  const double tmp_1 = gammaSphorCartDD00*gammaSphorCartDD11*gammaSphorCartDD22;
  const double tmp_3 = 2*gammaSphorCartDD01*gammaSphorCartDD02*gammaSphorCartDD12;
  const double tmp_5 = gammaSphorCartDD00*((gammaSphorCartDD12)*(gammaSphorCartDD12));
  const double tmp_7 = ((gammaSphorCartDD01)*(gammaSphorCartDD01))*gammaSphorCartDD22;
  const double tmp_9 = ((gammaSphorCartDD02)*(gammaSphorCartDD02))*gammaSphorCartDD11;
  const double tmp_10 = tmp_1 + tmp_3 - tmp_5 - tmp_7 - tmp_9;
  const double tmp_11 = (1.0/(tmp_10));
  const double tmp_12 = sin(xx1);
  const double tmp_13 = cbrt(tmp_11)*pow(fabs(tmp_12), 2.0/3.0)*pow(fabs(xx0), 4.0/3.0);
  const double tmp_14 = tmp_13/xx0;
  const double tmp_15 = (1.0/(tmp_12));
  const double tmp_17 = ((xx0)*(xx0));
  const double tmp_18 = (1.0/(tmp_17));
  const double tmp_20 = tmp_13*tmp_15*tmp_18;
  const double tmp_21 = ((tmp_12)*(tmp_12));
  const double tmp_22 = tmp_18/tmp_21;
  const double tmp_23 = 2*tmp_11;
  const double tmp_24 = KSphorCartDD00*tmp_11*(gammaSphorCartDD11*gammaSphorCartDD22 - ((gammaSphorCartDD12)*(gammaSphorCartDD12))) + KSphorCartDD01*tmp_23*(-gammaSphorCartDD01*gammaSphorCartDD22 + gammaSphorCartDD02*gammaSphorCartDD12) + KSphorCartDD02*tmp_23*(gammaSphorCartDD01*gammaSphorCartDD12 - gammaSphorCartDD02*gammaSphorCartDD11) + KSphorCartDD11*tmp_11*(gammaSphorCartDD00*gammaSphorCartDD22 - ((gammaSphorCartDD02)*(gammaSphorCartDD02))) + KSphorCartDD12*tmp_23*(-gammaSphorCartDD00*gammaSphorCartDD12 + gammaSphorCartDD01*gammaSphorCartDD02) + KSphorCartDD22*tmp_11*(gammaSphorCartDD00*gammaSphorCartDD11 - ((gammaSphorCartDD01)*(gammaSphorCartDD01)));
  const double tmp_25 = (1.0/3.0)*tmp_24;
  const double tmp_27 = tmp_11*tmp_21*((xx0)*(xx0)*(xx0)*(xx0));
  *hDD00 = gammaSphorCartDD00*tmp_13 - 1;
  *hDD01 = gammaSphorCartDD01*tmp_14;
  *hDD02 = gammaSphorCartDD02*tmp_14*tmp_15;
  *hDD11 = tmp_18*(gammaSphorCartDD11*tmp_13 - tmp_17);
  *hDD12 = gammaSphorCartDD12*tmp_20;
  *hDD22 = tmp_22*(gammaSphorCartDD22*tmp_13 - tmp_17*tmp_21);
  *aDD00 = tmp_13*(KSphorCartDD00 - gammaSphorCartDD00*tmp_25);
  *aDD01 = tmp_14*(KSphorCartDD01 - gammaSphorCartDD01*tmp_25);
  *aDD02 = tmp_14*tmp_15*(KSphorCartDD02 - gammaSphorCartDD02*tmp_25);
  *aDD11 = tmp_13*tmp_18*(KSphorCartDD11 - gammaSphorCartDD11*tmp_25);
  *aDD12 = tmp_20*(KSphorCartDD12 - gammaSphorCartDD12*tmp_25);
  *aDD22 = tmp_13*tmp_22*(KSphorCartDD22 - gammaSphorCartDD22*tmp_25);
  *trK = tmp_24;
  *vetU0 = betaSphorCartU0;
  *vetU1 = betaSphorCartU1*xx0;
  *vetU2 = betaSphorCartU2*tmp_12*xx0;
  *betU0 = BSphorCartU0;
  *betU1 = BSphorCartU1*xx0;
  *betU2 = BSphorCartU2*tmp_12*xx0;
  *alpha = alphaSphorCart;
  *cf = pow(tmp_10/(tmp_1*tmp_27 + tmp_27*tmp_3 - tmp_27*tmp_5 - tmp_27*tmp_7 - tmp_27*tmp_9), -1.0/6.0);
}
