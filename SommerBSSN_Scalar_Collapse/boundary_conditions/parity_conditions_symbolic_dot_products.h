/*
 *  Original SymPy expressions:
 *  "[parity[0] = 1,
 *    parity[1] = sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds),
 *    parity[2] = sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds),
 *    parity[3] = sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds),
 *    parity[4] = (sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds))**2,
 *    parity[5] = (sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds)),
 *    parity[6] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds)),
 *    parity[7] = (sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))**2,
 *    parity[8] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)),
 *    parity[9] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))**2]"
 */
{
  const double tmp_0 = cos(xx1)*cos(xx1_inbounds);
  const double tmp_1 = sin(xx1)*sin(xx1_inbounds);
  const double tmp_2 = sin(xx2)*sin(xx2_inbounds);
  const double tmp_3 = cos(xx2)*cos(xx2_inbounds);
  const double tmp_4 = tmp_0 + tmp_1*tmp_2 + tmp_1*tmp_3;
  const double tmp_5 = tmp_0*tmp_2 + tmp_0*tmp_3 + tmp_1;
  const double tmp_6 = tmp_2 + tmp_3;
  parity[0] = 1;
  parity[1] = tmp_4;
  parity[2] = tmp_5;
  parity[3] = tmp_6;
  parity[4] = ((tmp_4)*(tmp_4));
  parity[5] = tmp_4*tmp_5;
  parity[6] = tmp_4*tmp_6;
  parity[7] = ((tmp_5)*(tmp_5));
  parity[8] = tmp_5*tmp_6;
  parity[9] = ((tmp_6)*(tmp_6));
}
