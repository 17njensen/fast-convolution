/*
(C) Copyright Utah State University 1992-2002.  All rights reserved.
No part of this program may be photocopied, reproduced, or translated
to another program language without the prior written consent of
Utah State University. 

(program) fft842
(class) utility
(section) 2
(distribution) yes
(last_modified) Thu Jan 31 16:21:24 2002 by Scott Budge <scott@goga.ece.usu.edu>
(status) under development.
(special) routine
(description) Computes the single-precision FFT of a complex array.
(keywords) fft mixed-radix

(long_d)
This program replaces the complex array x by its finite
discrete, complex fourier transform.  It performs as many base
8 iterations as possible and then finishes with a base 4 iteration
or a base 2 iteration if needed.  (Define KR if K&R C is used.)

fft842(in, n, x)
int in,n;
complx x[];

in
If in=0, the forward fft is computed in place.  If in=1, the inverse
fft is computed.

n
This is the number of points in the fft.  It must be a power of two.

x
This is a complex array of points to be transformed.  The complx data
structure is defined by:

typedef struct {
	float re;
	float im;
	}complx;
(long_d) 

(see_also) fft842d

(bugs)
There are no known bugs at this time
(bugs)

(author) Scott Budge

(modifications)
6/6/90 SEB - eliminated the need for calling gentabfft.
Written 6/6/90 by Scott Budge
(modifications)

 SCCS Info: @(#) LUfft842.c 1.3@(#)   USU
 Version: 1.3
 Version created: 16:27:20, 01/31/02 
 This copy created: 16:27:22, 01/31/02

*/

static char SccsId[] = "@(#)LUfft842.c	1.3\t01/31/02\tUSU";

#define PI 3.14159265358979
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*  The following defines the largest size FFT that can be computed.
    This is 2^19 in length.  If this is changed, the code in fft842()
    *must* be changed.
*/
#define MAX_SIZE_LOG_2 19

typedef struct {
	float re;
	float im;
}complx;

static complx *stable;
static float rsr2d;
static int oldn = 0;


/************************************************************
  
  Function:
  This subroutine generates a table of
  constants for the computation of the FFT
  
  Parameter:
  n: dimension of the FFT
  
  -----------------------------------------------------------*/
#ifdef KR		/* K & R C */
gentabfft(n)
	 int n;
#else
void gentabfft(int n)
#endif
{
	double w,wi;
	register int i,m;
	register complx *pstable;
	
	if(n != oldn){
		m = n;
		if(m < 2){
			fprintf(stderr,
					"fft842: array size must be greater than 1.\n");
			exit(-2);
		}
		do{
			if(((m % 2) != 0) && (m != 1)){
				fprintf(stderr,
						"fft842: array size not a power of two.\n");
				exit(-2);
			}
		} while((m /= 2) != 0);
		rsr2d = 1.0 / sqrt(2.0);
		w = 2.0 * PI / n;
		if(oldn == 0)
			stable = (complx *) malloc(n * sizeof(complx));
		else {
            free(stable);
			stable = (complx *) malloc(n * sizeof(complx));
        }
		if(stable == NULL){
			fprintf(stderr,"fft842: array could not be allocated.\n");
			exit(-3);
		}
		oldn = n;
		for(i=0,pstable=stable ; i<n ; i++,pstable++) {
			wi = w * i;
			pstable->re = cos(wi);
			pstable->im = sin(wi);
		}
	}
}

/*----------------------------------------------------------------------
  C SUBROUTINE:  r2tx
  C RADIX 2 ITERATION SUBROUTINE
  C---------------------------------------------------------------------*/
#ifdef KR
r2tx(nthpo, c0, c1)
	 int nthpo;
	 complx *c0,*c1;
#else
void r2tx(int nthpo, complx *c0, complx *c1)
#endif
{
	register int k;
	register complx *pc0,*pc1;
	float r1,i1;
	
	for(k=0,pc0=c0,pc1=c1 ; k<nthpo ; k+=2,pc0+=2,pc1+=2) {
		r1 = pc0->re + pc1->re;
		pc1->re = pc0->re - pc1->re;
		pc0->re = r1;
		i1 = pc0->im + pc1->im;
		pc1->im = pc0->im - pc1->im;
		pc0->im = i1;
	}
}

/*----------------------------------------------------------------------
  C SUBROUTINE:  r4tx
  C RADIX 4 ITERATION SUBROUTINE
  C---------------------------------------------------------------------*/
#ifdef KR
r4tx(nthpo, c0, c1, c2, c3)
	 int nthpo;
	 complx *c0,*c1,*c2,*c3;
#else
void r4tx(int nthpo, complx *c0, complx *c1, complx *c2, complx *c3)
#endif
{
	register int k;
	float r1,r2,r3,r4,i1,i2,i3,i4;
	register complx *pc0,*pc1,*pc2,*pc3;
	
	for(k=0,pc0=c0,pc1=c1,pc2=c2,pc3=c3 ; k<nthpo ; k+=4,pc0+=4,pc1+=4,pc2+=4,pc3+=4) {
		r1 = pc0->re + pc2->re;
		r2 = pc0->re - pc2->re;
		r3 = pc1->re + pc3->re;
		r4 = pc1->re - pc3->re;
		i1 = pc0->im + pc2->im;
		i2 = pc0->im - pc2->im;
		i3 = pc1->im + pc3->im;
		i4 = pc1->im - pc3->im;
		pc0->re = r1 + r3;
		pc0->im = i1 + i3;
		pc1->re = r1 - r3;
		pc1->im = i1 - i3;
		pc2->re = r2 - i4;
		pc2->im = i2 + r4;
		pc3->re = r2 + i4;
		pc3->im = i2 - r4;
	}
}

/*----------------------------------------------------------------------
  C SUBROUTINE:  r8tx
  C RADIX 8 ITERATION SUBROUTINE
  C---------------------------------------------------------------------*/
#ifdef KR
r8tx(nxtlt, nthpo, lengt, i, c0, c1, c2, c3, c4, c5, c6, c7)
	 int nxtlt,nthpo,lengt,i;
	 complx *c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7;
#else
void r8tx(int nxtlt, int nthpo, int lengt, int i, complx *c0, complx *c1,
		  complx *c2, complx *c3, complx *c4, complx *c5, complx *c6,complx *c7)
#endif
{
	complx *pc0,*pc1,*pc2,*pc3,*pc4,*pc5,*pc6,*pc7;
	register int j,k;
	int arg,index;
	complx  a0,a1,a2,a3,a4,a5,a6,a7;
	complx  b0,b1,b2,b3,b4,b5,b6,b7;
	float co1,co2,co3,co4,co5,co6,co7,si1,si2,si3,si4,si5,si6,si7;
	float ti,tr;
	
	index = 1 << (3 * (i - 1));
	for(j=0 ; j<nxtlt ; j++) {
		arg = j * index;
		co1 = stable[arg].re;
		si1 = stable[arg].im;
		co2 = co1 * co1 - si1 * si1;
		si2 = co1*si1 + co1*si1;
		co3 = co1*co2 - si1*si2;
		si3 = co2*si1 + si2*co1;
		co4 = co2 * co2 - si2 * si2;
		si4 = co2*si2 + co2*si2;
		co5 = co2*co3 - si2*si3;
		si5 = co3*si2 + si3*co2;
		co6 = co3 * co3 - si3 * si3;
		si6 = co3*si3 + co3*si3;
		co7 = co3*co4 - si3*si4;
		si7 = co4*si3 + si4*co3;
		for(k=j,pc0=c0+j,pc1=c1+j,pc2=c2+j,pc3=c3+j,pc4=c4+j,pc5=c5+j,pc6=c6+j,pc7=c7+j ; k<nthpo ; k+=lengt,pc0+=lengt,pc1+=lengt,pc2+=lengt,pc3+=lengt,pc4+=lengt,pc5+=lengt,pc6+=lengt,pc7+=lengt) {
			a0.re = pc0->re + pc4->re;
			a1.re = pc1->re + pc5->re;
			a2.re = pc2->re + pc6->re;
			a3.re = pc3->re + pc7->re;
			a4.re = pc0->re - pc4->re;
			a5.re = pc1->re - pc5->re;
			a6.re = pc2->re - pc6->re;
			a7.re = pc3->re - pc7->re;
			a0.im = pc0->im + pc4->im;
			a1.im = pc1->im + pc5->im;
			a2.im = pc2->im + pc6->im;
			a3.im = pc3->im + pc7->im;
			a4.im = pc0->im - pc4->im;
			a5.im = pc1->im - pc5->im;
			a6.im = pc2->im - pc6->im;
			a7.im = pc3->im - pc7->im;
			b0.re = a0.re + a2.re;
			b1.re = a1.re + a3.re;
			b2.re = a0.re - a2.re;
			b3.re = a1.re - a3.re;
			b4.re = a4.re - a6.im;
			b5.re = a5.re - a7.im;
			b6.re = a4.re + a6.im;
			b7.re = a5.re + a7.im;
			b0.im = a0.im + a2.im;
			b1.im = a1.im + a3.im;
			b2.im = a0.im - a2.im;
			b3.im = a1.im - a3.im;
			b4.im = a4.im + a6.re;
			b5.im = a5.im + a7.re;
			b6.im = a4.im - a6.re;
			b7.im = a5.im - a7.re;
			pc0->re = b0.re + b1.re;
			pc0->im = b0.im + b1.im;
			if(j>0) {
				pc1->re = co4*(b0.re-b1.re) - si4*(b0.im-b1.im);
				pc1->im = co4*(b0.im-b1.im) + si4*(b0.re-b1.re);
				pc2->re = co2*(b2.re-b3.im) - si2*(b2.im+b3.re);
				pc2->im = co2*(b2.im+b3.re) + si2*(b2.re-b3.im);
				pc3->re = co6*(b2.re+b3.im) - si6*(b2.im-b3.re);
				pc3->im = co6*(b2.im-b3.re) + si6*(b2.re+b3.im);
				tr = rsr2d*(b5.re-b5.im);
				ti = rsr2d*(b5.re+b5.im);
				pc4->re = co1*(b4.re+tr) - si1*(b4.im+ti);
				pc4->im = co1*(b4.im+ti) + si1*(b4.re+tr);
				pc5->re = co5*(b4.re-tr) - si5*(b4.im-ti);
				pc5->im = co5*(b4.im-ti) + si5*(b4.re-tr);
				tr = -rsr2d*(b7.re+b7.im);
				ti = rsr2d*(b7.re-b7.im);
				pc6->re = co3*(b6.re+tr) - si3*(b6.im+ti);
				pc6->im = co3*(b6.im+ti) + si3*(b6.re+tr);
				pc7->re = co7*(b6.re-tr) - si7*(b6.im-ti);
				pc7->im = co7*(b6.im-ti) + si7*(b6.re-tr);
  			} else {
  				pc1->re = b0.re - b1.re;
				pc1->im = b0.im - b1.im;
				pc2->re = b2.re - b3.im;
				pc2->im = b2.im + b3.re;
				pc3->re = b2.re + b3.im;
				pc3->im = b2.im - b3.re;
				tr = rsr2d*(b5.re-b5.im);
				ti = rsr2d*(b5.re+b5.im);
				pc4->re = b4.re + tr;
				pc4->im = b4.im + ti;
				pc5->re = b4.re - tr;
				pc5->im = b4.im - ti;
				tr = -rsr2d*(b7.re+b7.im);
				ti = rsr2d*(b7.re-b7.im);
				pc6->re = b6.re + tr;
				pc6->im = b6.im + ti;
				pc7->re = b6.re - tr;
				pc7->im = b6.im - ti;
  			}
		}
	}
}

/*---------------------------------------------------------------------
  
  This routine is based in the IEEE standard routine.
  
  ---------------------------------------------------------------------*/
#ifdef KR
fft842(in, n, x)
	 int in,n;
	 complx x[];
#else
void fft842(int in, int n, complx *x)
#endif
{
	register int i;
	int m,nt,n2pow,nthpo,n8pow,nxtlt,lengt,l[MAX_SIZE_LOG_2],*pl;
	int ij,ji,j,j1a,j2,j3,j4,j5,j6,j7,j8,j9,j1a0,j1a1,j1a2,j1a3,j1a4;
    int j1b,j1c,j1d,j1e;
	float r,fi;
	register complx *px,*px1;
	
	if(n != oldn)
		gentabfft(n);
	nt = 1;
	i=1;
	while(n!=nt) {
		m = i;
		nt <<= 1;
		i++;
	}
	n2pow = m;

    if(m > MAX_SIZE_LOG_2){
        fprintf(stderr,
                "fft842: array size must be smaller than 2^20!\n");
        exit(-2);
    }
        
	nthpo = n;
	if(in != 1)
		for(i=0,px=x ; i<nthpo ; i++,px++)
			px->im = -px->im;
	n8pow = n2pow/3;
	if(n8pow != 0)
		/*
		  C RADIX 8 PASSES,IF ANY.
		  */
		for(i=1 ; i<=n8pow ; i++) {
			nxtlt = 1<<(n2pow-3*i);
			lengt = 8*nxtlt;
			r8tx(nxtlt, nthpo, lengt, i, x, x+nxtlt, x+2*nxtlt,x+3*nxtlt,
				 x+4*nxtlt, x+5*nxtlt, x+6*nxtlt,x+7*nxtlt);
		}
	/*
	  C IS THERE A FOUR FACTOR LEFT
	  */
	j = n2pow - 3 * n8pow - 1;
	if(j==0)
		/*
		  C GO THROUGH THE BASE 2 ITERATION
		  */
		r2tx(nthpo, x, x+1);
	if(j>0)
		/*
		  C GO THROUGH THE BASE 4 ITERATION
		  */
		r4tx(nthpo, x, x+1, x+2, x+3);
	for(i=0,pl=l ; i<MAX_SIZE_LOG_2 ; i++,pl++) {
		*pl = 1;
		if(i<n2pow)
			*pl = 1<<(n2pow-i);
	}
	px1 = x;
	ij = 1;
	for(j1e=1 ; j1e<=l[18] ; j1e++)
	 for(j1d=j1e ; j1d<=l[17] ; j1d+=l[18])
	  for(j1c=j1d ; j1c<=l[16] ; j1c+=l[17])
	   for(j1b=j1c ; j1b<=l[15] ; j1b+=l[16])
	    for(j1a=j1b ; j1a<=l[14] ; j1a+=l[15])
	     for(j2=j1a ; j2<=l[13] ; j2+=l[14])
		  for(j3=j2 ; j3<=l[12] ; j3+=l[13])
		   for(j4=j3 ; j4<=l[11] ; j4+=l[12])
			for(j5=j4 ; j5<=l[10] ; j5+=l[11])
			 for(j6=j5 ; j6<=l[9] ; j6+=l[10])
			   for(j7=j6 ; j7<=l[8] ; j7+=l[9])
				 for(j8=j7 ; j8<=l[7] ; j8+=l[8])
				   for(j9=j8 ; j9<=l[6] ; j9+=l[7])
					 for(j1a0=j9 ; j1a0<=l[5] ; j1a0+=l[6])
					   for(j1a1=j1a0 ; j1a1<=l[4] ; j1a1+=l[5])
						 for(j1a2=j1a1 ; j1a2<=l[3] ; j1a2+=l[4])
						   for(j1a3=j1a2 ; j1a3<=l[2] ; j1a3+=l[3])
							 for(j1a4=j1a3 ; j1a4<=l[1] ; j1a4+=l[2])
							   for(ji=j1a4,px=x+j1a4-1 ; ji<= l[0] ; ji+=l[1],px+=l[1],px1++,ij++)
								  if(ij < ji) {
									  r = px1->re;
									  px1->re = px->re;
									  px->re = r;
									  fi = px1->im;
									  px1->im = px->im;
									  px->im = fi;
								  }
	if(in==1)
		for(i=0,px=x ; i<nthpo ; i++,px++) {
			px->re /= n;
			px->im /= n;
		}
	else
		for(i=0,px=x ; i<nthpo ; i++,px++)
			px->im = -px->im;
}
