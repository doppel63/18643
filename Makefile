all:
	g++ -Wall -lm -g fft1k.cpp fft1k_test.cpp

clean:
	rm a.out *~

