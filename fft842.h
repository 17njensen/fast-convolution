#include <stdio.h>
#include <stdlib.h>

typedef struct {
	float re;
	float im;
}complx;

// void r8tx(int nxtlt, int nthpo, int lengt, int i, complx *c0, complx *c1,
// 		  complx *c2, complx *c3, complx *c4, complx *c5, complx *c6,complx *c7);

// void r4tx(int nthpo, complx *c0, complx *c1, complx *c2, complx *c3);

// void r2tx(int nthpo, complx *c0, complx *c1);

// void gentabfft(int n);

void fft842(int in, int n, complx *x);

