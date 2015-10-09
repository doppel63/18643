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
  int i, j;
  float wr, wi;
  float tempr, tempi;
  int table, entry;

  // constants?
  for (table=0; table<10; table++) {
	  for (entry=0; entry<512; entry++) {
		  consts_t en = coeff[table][entry];
		  i = en.i;
		  j = en.j;
		  wr = en.wr;
		  wi = en.wi;
		  tempr = wr*data[2*j]   - wi*data[2*j+1];
		  tempi = wr*data[2*j+1] + wi*data[2*j];
		  data[2*j]   = data[2*i]   - tempr;
		  data[2*j+1] = data[2*i+1] - tempi;
		  data[2*i] += tempr;
		  data[2*i+1] += tempi;
	  }
  }
}
