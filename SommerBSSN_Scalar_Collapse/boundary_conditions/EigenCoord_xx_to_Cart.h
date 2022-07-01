
static inline void EigenCoord_xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
#include "../set_Cparameters.h"
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = xx0*sin(xx1)*cos(xx2),
   *    xCart[1] = xx0*sin(xx1)*sin(xx2),
   *    xCart[2] = xx0*cos(xx1)]"
   */
  {
    const double tmp_0 = xx0*sin(xx1);
    xCart[0] = tmp_0*cos(xx2);
    xCart[1] = tmp_0*sin(xx2);
    xCart[2] = xx0*cos(xx1);
  }
}
