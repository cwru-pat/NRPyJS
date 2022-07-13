rfm_struct rfmstruct;
rfmstruct.f0_of_xx0 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
rfmstruct.f1_of_xx1 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
rfmstruct.f1_of_xx1__D1 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
rfmstruct.f1_of_xx1__DD11 = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
