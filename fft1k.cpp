//-------------------------------------------
//- CMU 18-643 Fall 2015, Vivado HLS Lab
//-------------------------------------------
#include <math.h>
#include "fft1k.h"
#include "luts.h"
#include <hls_stream.h>


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
// this is only the butterfly portion
void fft1k_br(float data[2*SIZE])
{
#pragma HLS ARRAY_RESHAPE variable=data cyclic factor=2 dim=1

  ap_uint<10> i, j;
  float wr, wi;
  float d_2i, d_2j, d_2i_p1, d_2j_p1; //temporary storage for reads
  float tempr, tempi;

  int table, entry;

  Stage_Loop: for (table=0; table<10; table += 2) {
	  float data2[2*SIZE]; // for double buffering
#pragma HLS ARRAY_RESHAPE variable=data2 cyclic factor=2 dim=1

	  // this loop consumes memory and writes to buffer
	  Row_Loop1: for (entry=0; entry<512; entry++) {
#pragma HLS DEPENDENCE variable=data intra false
#pragma HLS DEPENDENCE variable=data inter false
#pragma HLS DEPENDENCE variable=data2 intra false
#pragma HLS DEPENDENCE variable=data2 inter false
#pragma HLS PIPELINE

		  consts_t en = coeff[table][entry];
#pragma HLS data_pack variable=coeff

		  // fetch from lookup table
		  i = en.i;
		  j = en.j;
		  wr = en.wr;
		  wi = en.wi;

		  // do all "loading" into "registers"
		  // avoid double accesses later
		  d_2j = data[2*j];
		  d_2j_p1 = data[2*j+1];
		  d_2i = data[2*i];
		  d_2i_p1 = data[2*i+1];

		  // do the butterfly
		  tempr = wr*d_2j - wi*d_2j_p1;
		  tempi = wr*d_2j_p1 + wi*d_2j;

		  float outj = d_2i - tempr;
		  float outj1 = d_2i_p1 - tempi;
		  float outi = d_2i + tempr;
		  float outi1 = d_2i_p1 + tempi;

		  // write internal (buffer1)
		  data2[2*i] = outi;
		  data2[2*i+1] = outi1;
		  data2[2*j] = outj;
		  data2[2*j+1] = outj1;
	  }

	  // this loop consumes buffer and writes to memory
	  Row_Loop2: for (entry=0; entry<512; entry++) {
#pragma HLS DEPENDENCE variable=data intra false
#pragma HLS DEPENDENCE variable=data inter false
#pragma HLS DEPENDENCE variable=data2 intra false
#pragma HLS DEPENDENCE variable=data2 inter false
#pragma HLS PIPELINE

		  consts_t en = coeff[table+1][entry];
#pragma HLS data_pack variable=coeff

		  // fetch from lookup table
		  i = en.i;
		  j = en.j;
		  wr = en.wr;
		  wi = en.wi;

		  // do all "loading" into "registers"
		  // avoid double accesses later
		  d_2j = data2[2*j];
		  d_2j_p1 = data2[2*j+1];
		  d_2i = data2[2*i];
		  d_2i_p1 = data2[2*i+1];

		  // do the butterfly
		  tempr = wr*d_2j  - wi*d_2j_p1;
		  tempi = wr*d_2j_p1 + wi*d_2j;

		  float outj = d_2i - tempr;
		  float outj1 = d_2i_p1 - tempi;
		  float outi = d_2i + tempr;
		  float outi1 = d_2i_p1 + tempi;

		  // write to "real output" (buffer2)
		  data[2*i] = outi;
		  data[2*i+1] = outi1;
		  data[2*j] = outj;
		  data[2*j+1] = outj1;
	  }
  }
}
