#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "filter.h"
#include <time.h>
#include "fft842.h"
/*
    Classic convolution:
        -Number of multiplies and adds:
            Lh-1 adds
            Lh-1 multiplys




*/
typedef struct 
{ 
 int ndim; //number of dimensions 
 int nchan; //number of channels 
 int d0;  //length of the first dimension 
 int d1;  //length of second dimension or sample rate if audio 
 int d2;  //length of third dimension 
} dsp_file_header; 


/**
 * @brief Convolution function
 * @param x is input
 * @param h is what x gets convovled with
 * @param y is the output
 * @param Lx is len of input
 * @param Lh is len of what gets convolved
 * @param ly is the zero padded output length
 */
float* conv(float* x, float* h, float* y, int Lx, int Lh, int Ly){
    printf("Start Convolution\n");
    int i,j;
    for (i = 0; i < Ly; i++) { 
        for (j = 0; j < Lh; j++) { 
            //printf("x[%d+%d] = %f\n",i,j, x[i+j]);
            //printf("x[10] = %f\n", x[10]);
            y[i] += h[j] * x[i + j]; //multiply and accumulate (MAC) 
            //impulse is assumed to be in time reverse order 
        } 
        // printf("y[i] = %f\n", y[i]);
    }
    return y;
}

unsigned int nextPowerOf2(unsigned int n){
    unsigned count = 0;
    
    // First n in the below condition
    // is for the case where n is 0
    if (n && !(n & (n - 1)))
        return n;
    
    while( n != 0)
    {
        n >>= 1;
        count += 1;
    }
    
    return 1 << count;
}
 

void main(int argc, char** argv){
    //read in file
    FILE* fx;
    FILE* fy;
    FILE* ff;
    if (NULL == (fx = fopen(argv[1], "rb"))) { //error check and open file
        printf("error: Cannot open input file.\n");
        return;
    }
    if (NULL == (fy = fopen(argv[2], "wb"))) { //error check and open file
        printf("error: Cannot open output file for writing.\n");
        return;
    }
    if (NULL == (ff = fopen(argv[3], "wb"))) { //error check and open file
        printf("error: Cannot open output file for writing.\n");
        return;
    }
    //grab headers of each file 
    dsp_file_header h0, ho, hf; 
    fread(&h0, sizeof(dsp_file_header), 1, fx); 
    memcpy(&ho, &h0, sizeof(dsp_file_header));
    memcpy(&hf, &h0, sizeof(dsp_file_header));
    fwrite(&ho, sizeof(dsp_file_header), 1, fy); 
    fwrite(&hf, sizeof(dsp_file_header), 1, ff); 
    printf("ndim = %d, nchan = %d, d0 = %d, d1 = %d, d2 = %d\n", h0.ndim, h0.nchan, h0.d0, h0.d1, h0.d2);
    int Lh = sizeof(filter)/sizeof(filter[0]);//length of impulse signal 
    int Lx = h0.d0; //length of input signal 
    Lx = nextPowerOf2(Lx);
    int diff = (Lx - h0.d0)/2;
    int Ly = h0.d0 + (Lh - 1); //len of conv result 
    int Lw = Lx + (Lh-1);
    int Lz = h0.d0 + 2 * (Lh - 1); //len of zero padded input 
    printf("Lh = %d, Lx = %d, Ly = %d, Lz = %d, h0.d0 = %d\n", Lh, Lx, Ly, Lz, h0.d0); 

    //allocate data space
    float* x = calloc(sizeof(float), Lz); //for conv
    float* pad_h = calloc(sizeof(float), h0.d0); 
    float* z = calloc(sizeof(float), Lx); //for fft
    float* y = calloc(sizeof(float), Ly);

    while (!feof(fx)) { 
        fread((x + Lh - 1), sizeof(float), h0.d0, fx); 
        // fread((x), sizeof(float), h0.d0, fx); 
    } 
    for(int i = 0; i < Lx; i++){
        z[i] = x[i+Lh-1]; //z is no zero padding x
    }
    int Lpad = (h0.d0 - Lh)/2;
    printf("Lpad = %d\n", Lpad);
    for(int i = 0; i < Lh; i++){
        pad_h[i+Lpad-1] = filter[i];
    }
    //-------------------filter work - linear convolution----------------------------------
    // float* h_1 = calloc(sizeof(float), Lh);
    // for(int i = 0; i < Lh; i++){
    //     h_1[i] = filter[(Lh-1)-i];
    // }
    // y = conv(x, h_1, y, Lz, Lh, Ly);
    // fwrite(y, sizeof(float), Ly, fy); 
    // printf("Done.\n");

    // //------------------filter work - circular convolution --------------------------------
    printf("Start Circ Conv\n");
    clock_t begin = clock();
    for (int i = 0; i < h0.d0; i++) { 
        for (int j = 0; j < Lh; j++) { 
            y[i] +=  filter[j]*x[(i-j)%h0.d0];
        } 
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;

    printf("Time of circular convolution = %0.20fs\n", time_spent);
    fwrite(y, sizeof(float), Ly, fy);
    free(y);free(pad_h);
    //-------------------fft842.c work--------------------------------
    // complx* cosine = calloc(sizeof(complx), Lx);
    complx* h = calloc(sizeof(complx), Lx);
    // complx* output = calloc(sizeof(complx), Lx); 
    float* fft_out = calloc(sizeof(float), Lx);
    // for(int i = 0; i < Lx; i++){
    //     cosine[i].re = z[i];
    //     cosine[i].im = 0;
    // }
    for(int i = 0; i < Lh; i++){
        h[i].re = filter[i];
        h[i].im = 0;
    }
    //overlap add
    int L = 256; //to find step size
    int k = Lx/Lh; //number of segments divided out for overlap add
    printf("Lx = %d, L = %d, k = %d\n",Lx,L,k);
    // float** x_add = calloc(sizeof(float*), k); //2d allocated array
    float ** x_add = (float **)calloc(k, sizeof(float *));
    // complx** x_add_fft = calloc(sizeof(complx*), k);
    complx ** x_add_fft = (complx **)calloc(k, sizeof(complx*));
    // complx** y_add = calloc(sizeof(complx*), k);
    complx** y_add = calloc(k, sizeof(complx *));
    for(int i = 0; i < k; i++){
        x_add[i] = (float *)calloc(L, sizeof(float));
        x_add_fft[i] = (complx *)calloc(L, sizeof(complx));
        y_add[i] = (complx *)calloc(L, sizeof(complx));
    }
    printf("Start overlap.\n");
    for (int m = 0; m < k; m++){
        for(int i = 0; i < L; i++){
            // printf("i = %d\n",i);
            x_add[m][i] = z[i+m*L]; //separate x into segments of L, no zero padding input
            x_add_fft[m][i].re = x_add[m][i];
            x_add_fft[m][i].im = 0;
        }
        printf("m = %d\n", m);
    }
    printf("Start FFT\n");
    clock_t begin1 = clock();
    for(int m = 0; m < k; m++){ //fft of each segmet
        fft842(0, L, x_add_fft[m]); 
    }
    fft842(0, Lx, h);
    // fft842(0, Lx, cosine);
    for(int m = 0; m < k; m++){
        for(int i = 0; i < Lx; i++){
            y_add[m][i].re = x_add_fft[m][i].re*h[i].re;
            y_add[m][i].im = x_add_fft[m][i].im*h[i].im;
        }
        fft842(1, Lx, y_add[m]);
    }
    clock_t end1 = clock();
    double time_spent1 = (double)(end1 - begin1)/CLOCKS_PER_SEC;
    printf("Time of fft = %0.20fs\n", time_spent1);
    for(int m = 0; m < k; m++){
        for(int i = 0; i < Lx; i++){
            fft_out[i] += y_add[m][i-m*L].re;
        }
    }
    fwrite(fft_out, sizeof(float), Lx, ff);

    printf("Done.\n");
    fclose(fx); 
    fclose(fy); 
    fclose(ff);
}
