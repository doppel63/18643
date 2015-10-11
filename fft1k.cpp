//-------------------------------------------
//- CMU 18-643 Fall 2015, Vivado HLS Lab
//-------------------------------------------
#include <math.h>
#include "fft1k.h"
#include "luts.h"

// nn must be a two-power
// 0<=x<nn
// returns bit-reverse of x
unsigned int reverseBits(unsigned int x, int nn) {
  unsigned int y=0;
  nn>>=1;

  while (nn) {
    y<<=1;
    if (x&1) {
      y|=1;
    }
    nn>>=1;
    x>>=1;
  }

  return y;
}

// permuate a nn-entry complex valued array into bit-reverse order
// i'th real at data[2*i], i'th imanginary at data[2*i+1]
void br(float data[], int nn)
{
  float tempr;
  int i, j;

  for (i = 0; i < nn; i += 1) {
    j=reverseBits(i,nn);
    if (j > i) {
      tempr = data[2*j];     data[2*j] = data[2*i];     data[2*i] = tempr;
      tempr = data[2*j+1]; data[2*j+1] = data[2*i+1]; data[2*i+1] = tempr;
    }
  }
}

// derived from four1, index complex pairs from 0
// assumes input already in bit reverse order
// this is only the bufferfly portion
void fft1k_br(float data[2*SIZE])
{
#pragma HLS ARRAY_PARTITION variable=data cyclic factor=2 dim=1
  ap_uint<10> i, j;
  COEFF_TYPE wr, wi;
  COEFF_TYPE d2i, d2j, d2i1, d2j1;
  COEFF_TYPE tempr, tempi;

  int table, entry;

  // constants?
  Stage_Loop: for (table=0; table<10; table++) {
	  Row_Loop: for (entry=0; entry<512; entry++) {
#pragma HLS DEPENDENCE variable=data inter false
#pragma HLS PIPELINE
		  // do all "loading" into "registers"
		  consts_t en = coeff[table][entry];
		  i = en.i;
		  j = en.j;
		  wr = en.wr;
		  wi = en.wi;

		  d2i = data[2*i];
		  d2j = data[2*j];
		  d2i1 = data[2*i+1];
		  d2j1 = data[2*j+1];

		  tempr = wr*d2j  - wi*d2j1;
		  tempi = wr*d2j1 + wi*d2j;
		  data[2*j]   = d2i  - tempr;
		  data[2*j+1] = d2i1 - tempi;
		  data[2*i] = d2i + tempr;
		  data[2*i+1] = d2i1 + tempi;
	  }
  }
}
