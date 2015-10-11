//-------------------------------------------
//- CMU 18-643 Fall 2015, Vivado HLS Lab
//-------------------------------------------

#define SIZE (1024)
#define FORWARD 1

// make it easier to play with coefficient data types
#define COEFF_TYPE float
//#define COEFF_TYPE ap_fixed<1,1>

void br(float data[], int nn);

void fft1k_br(float data[2*SIZE]);
