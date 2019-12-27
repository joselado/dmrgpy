if (numprod==1) {
  string op0; 
  int i0; 
  hfile >> cr >> ci >> op0 >> i0; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0;
}

else if (numprod==2) {
  string op0,op1; 
  int i0,i1; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1;
}

else if (numprod==3) {
  string op0,op1,op2; 
  int i0,i1,i2; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2;
}

else if (numprod==4) {
  string op0,op1,op2,op3; 
  int i0,i1,i2,i3; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3;
}

else if (numprod==5) {
  string op0,op1,op2,op3,op4; 
  int i0,i1,i2,i3,i4; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4;
}

else if (numprod==6) {
  string op0,op1,op2,op3,op4,op5; 
  int i0,i1,i2,i3,i4,i5; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5;
}

else if (numprod==7) {
  string op0,op1,op2,op3,op4,op5,op6; 
  int i0,i1,i2,i3,i4,i5,i6; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6;
}

else if (numprod==8) {
  string op0,op1,op2,op3,op4,op5,op6,op7; 
  int i0,i1,i2,i3,i4,i5,i6,i7; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7;
}

else if (numprod==9) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8;
}

else if (numprod==10) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9;
}

else if (numprod==11) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10;
}

else if (numprod==12) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11;
}

else if (numprod==13) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12;
}

else if (numprod==14) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13;
}

else if (numprod==15) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14;
}

else if (numprod==16) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15;
}

else if (numprod==17) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16;
}

else if (numprod==18) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17;
}

else if (numprod==19) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18;
}

else if (numprod==20) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19;
}

else if (numprod==21) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20;
}

else if (numprod==22) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21;
}

else if (numprod==23) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22;
}

else if (numprod==24) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23;
}

else if (numprod==25) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24;
}

else if (numprod==26) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25;
}

else if (numprod==27) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26;
}

else if (numprod==28) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27;
}

else if (numprod==29) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28;
}

else if (numprod==30) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29;
}

else if (numprod==31) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30;
}

else if (numprod==32) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31;
}

else if (numprod==33) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32;
}

else if (numprod==34) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32,op33; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32 >> op33 >> i33; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32,op33,i33;
}

else if (numprod==35) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32,op33,op34; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32 >> op33 >> i33 >> op34 >> i34; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32,op33,i33,op34,i34;
}

else if (numprod==36) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32,op33,op34,op35; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,i35; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32 >> op33 >> i33 >> op34 >> i34 >> op35 >> i35; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32,op33,i33,op34,i34,op35,i35;
}

else if (numprod==37) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32,op33,op34,op35,op36; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,i35,i36; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32 >> op33 >> i33 >> op34 >> i34 >> op35 >> i35 >> op36 >> i36; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32,op33,i33,op34,i34,op35,i35,op36,i36;
}

else if (numprod==38) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32,op33,op34,op35,op36,op37; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,i35,i36,i37; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32 >> op33 >> i33 >> op34 >> i34 >> op35 >> i35 >> op36 >> i36 >> op37 >> i37; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32,op33,i33,op34,i34,op35,i35,op36,i36,op37,i37;
}

else if (numprod==39) {
  string op0,op1,op2,op3,op4,op5,op6,op7,op8,op9,op10,op11,op12,op13,op14,op15,op16,op17,op18,op19,op20,op21,op22,op23,op24,op25,op26,op27,op28,op29,op30,op31,op32,op33,op34,op35,op36,op37,op38; 
  int i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,i35,i36,i37,i38; 
  hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >> i3 >> op4 >> i4 >> op5 >> i5 >> op6 >> i6 >> op7 >> i7 >> op8 >> i8 >> op9 >> i9 >> op10 >> i10 >> op11 >> i11 >> op12 >> i12 >> op13 >> i13 >> op14 >> i14 >> op15 >> i15 >> op16 >> i16 >> op17 >> i17 >> op18 >> i18 >> op19 >> i19 >> op20 >> i20 >> op21 >> i21 >> op22 >> i22 >> op23 >> i23 >> op24 >> i24 >> op25 >> i25 >> op26 >> i26 >> op27 >> i27 >> op28 >> i28 >> op29 >> i29 >> op30 >> i30 >> op31 >> i31 >> op32 >> i32 >> op33 >> i33 >> op34 >> i34 >> op35 >> i35 >> op36 >> i36 >> op37 >> i37 >> op38 >> i38; 
  auto cz = cr + ci*1i;
  ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3,op4,i4,op5,i5,op6,i6,op7,i7,op8,i8,op9,i9,op10,i10,op11,i11,op12,i12,op13,i13,op14,i14,op15,i15,op16,i16,op17,i17,op18,i18,op19,i19,op20,i20,op21,i21,op22,i22,op23,i23,op24,i24,op25,i25,op26,i26,op27,i27,op28,i28,op29,i29,op30,i30,op31,i31,op32,i32,op33,i33,op34,i34,op35,i35,op36,i36,op37,i37,op38,i38;
}

else exit(EXIT_FAILURE) ;
