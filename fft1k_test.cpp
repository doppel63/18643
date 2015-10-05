//-------------------------------------------
//- CMU 18-643 Fall 2015, Vivado HLS Lab
//-------------------------------------------

#include <iostream>
#include <math.h>

using namespace std;

#include "fft1k.h"
#include "fft1k_ref0.h"
//#include "fft1k_ref1.h"
//#include "fft1k_ref2.h"
//#include "fft1k_ref3.h"
//#include "fft1k_ref4.h"

int main(int argc, char **argv) {
  int error_count = 0;
  
  // perform bit reverse
  br(&data[0], SIZE);
 
  // perform butterfly
  fft1k_br(&data[0]);
  
  // check against result
  {
    int i;
    for(i=0;i<SIZE;i++) {
      float realEr=((fabsf(result[2*i])<ERROROK)?
		    fabsf(data[2*i]-result[2*i]):
		    fabsf((data[2*i]-result[2*i])/result[2*i]));
      float imagineEr=((fabsf(result[2*i+1])<ERROROK)?
		       fabsf(data[2*i+1]-result[2*i+1]):
		       fabsf((data[2*i+1]-result[2*i+1])/result[2*i+1]));
		       
      if ((realEr>ERROROK)||(imagineEr>ERROROK)) {
	cout << "i=" << i << ": " << data[2*i] << "+" << data[2*i+1] << "i! ="
	     << result[2*i] << "+" << result[2*i+1] << "i" << endl;
	error_count++;
      }
    }
  }
  
  if (error_count)
    cout << "TEST FAIL: " << error_count << " Results do not match!" << endl;
  else
    cout << "Test passed!" << endl;

  return error_count;
}

