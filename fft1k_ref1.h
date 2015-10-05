//-------------------------------------------
//- CMU 18-643 Fall 2015, Vivado HLS Lab
//-------------------------------------------

#include "fft1k.h"

#define ERROROK (0.0001)

float data[2*SIZE] = {
  1, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0,
  0, 0
};

float result[2*SIZE] = {
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0,
  1, 0
};

