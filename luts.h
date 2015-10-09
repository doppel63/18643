typedef struct{
	int i;	// replace with 10 bits!
	int j;
	float wr;
	float wi;
} consts_t;

// LUTs
static const consts_t coeff[10][512] = {
#include "fft_coeff_table.txt"
};