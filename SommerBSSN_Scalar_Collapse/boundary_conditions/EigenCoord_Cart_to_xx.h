/*
 *  Original SymPy expressions:
 *  "[Cart_to_xx0_inbounds = sqrt(Cartx**2 + Carty**2 + Cartz**2),
 *    Cart_to_xx1_inbounds = acos(Cartz/sqrt(Cartx**2 + Carty**2 + Cartz**2)),
 *    Cart_to_xx2_inbounds = atan2(Carty, Cartx)]"
 */
{
  const double tmp_0 = sqrt(((Cartx)*(Cartx)) + ((Carty)*(Carty)) + ((Cartz)*(Cartz)));
  Cart_to_xx0_inbounds = tmp_0;
  Cart_to_xx1_inbounds = acos(Cartz/tmp_0);
  Cart_to_xx2_inbounds = atan2(Carty, Cartx);
}
