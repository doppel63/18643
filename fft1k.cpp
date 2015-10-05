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
  int nn=SIZE;
  int isign=FORWARD;
  int n, mmax, m, j, istep, i;
  float wtemp, wr, wpr, wpi, wi, theta;
  float tempr, tempi;

  n = nn;
MMAX_LOOP:
  for(mmax = 1; mmax < n; mmax = mmax*2) {
    istep = 2*mmax;
    theta = theta_lookup[own_log2(mmax)];
//    theta = M_PI/(isign*mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
M_LOOP:
	for (m = 0; m < mmax; m += 1) {
I_LOOP:
      for (i = m; i < n; i += istep) {
		j =i + mmax;
		tempr = wr*data[2*j]   - wi*data[2*j+1];
		tempi = wr*data[2*j+1] + wi*data[2*j];
		data[2*j]   = data[2*i]   - tempr;
		data[2*j+1] = data[2*i+1] - tempi;
		data[2*i] += tempr;
		data[2*i+1] += tempi;
      }
      wr = (wtemp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
  }
}

// helper function for log2(mmax) in fft1k_br
// only covers possible values of mmax with given constants
int own_log2(int in) {
  switch (in) {
  case 1:	return 0;
  case 2:	return 1;
  case 4:	return 2;
  case 8:	return 3;
  case 16:	return 4;
  case 32:	return 5;
  case 64:	return 6;
  case 128:	return 7;
  case 256:	return 8;
  case 512:	return 9;
  default:	return 0;
  }
}
