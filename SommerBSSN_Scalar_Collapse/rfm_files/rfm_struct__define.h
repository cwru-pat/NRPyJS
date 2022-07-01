for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
  const REAL xx0 = xx[0][i0];
  rfmstruct.f0_of_xx0[i0] = xx0;
}

for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
  const REAL xx1 = xx[1][i1];
  rfmstruct.f1_of_xx1[i1] = sin(xx1);
}

for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
  const REAL xx1 = xx[1][i1];
  rfmstruct.f1_of_xx1__D1[i1] = cos(xx1);
}

for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) {
  const REAL xx1 = xx[1][i1];
  rfmstruct.f1_of_xx1__DD11[i1] = -sin(xx1);
}

