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
    int Ly = h0.d0 + (Lh - 1); //len of conv result 
    int Lw = Lx + (Lh-1);
    int Lz = h0.d0 + 2 * (Lh - 1); //len of zero padded input 
    printf("Lh = %d, Lx = %d, Ly = %d, Lz = %d, h0.d0 = %d\n", Lh, Lx, Ly, Lz, h0.d0); 

    //allocate data space
    float* x = calloc(sizeof(float), Lz); //for conv
    float* pad_h = calloc(sizeof(float), h0.d0); 
    float* z = calloc(sizeof(float), h0.d0+Lh-1); //for fft
    float* y = calloc(sizeof(float), Ly);

    while (!feof(fx)) { 
        fread((x + Lh - 1), sizeof(float), h0.d0, fx); 
        // fread((x), sizeof(float), h0.d0, fx); 
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
    // printf("Start Circ Conv\n");
    // clock_t begin = clock();
    // for (int i = 0; i < h0.d0; i++) { 
    //     for (int j = 0; j < Lh; j++) { 
    //         y[i] +=  filter[j]*x[(i-j)%h0.d0];
    //     } 
    // }
    // clock_t end = clock();
    // double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;

    // printf("Time of circular convolution = %0.20fs\n", time_spent);
    // fwrite(y, sizeof(float), Ly, fy);
    // free(y);free(pad_h);
    //-------------------fft842.c work--------------------------------
    //overlap save
    int L = 769;
    int N = L+Lh-1;
    int k = (h0.d0+Lh-1)/(L); //number of segments divided out for overlap add
    printf("Lx = %d, L = %d, k = %d\n",Lx,L,k);
    float ** x_s = (float **)calloc(k, sizeof(float *));
    complx ** x_s_fft = (complx **)calloc(k, sizeof(complx*));
    complx** y_s = calloc(k, sizeof(complx *));
    for(int i = 0; i < k; i++){
        x_s[i] = (float *)calloc(N, sizeof(float)); //make into 2d array
        x_s_fft[i] = (complx *)calloc(N, sizeof(complx)); //make into 2d array
        y_s[i] = (complx *)calloc(N, sizeof(complx)); //make into 2d array
    }
    printf("Start overlap.\n");
    //insert Lh-1 zeros
    for(int i = 0; i < h0.d0+Lh-1; i++){
        z[i] = x[i];
    }
    //segment input x
    for(int r = 0; r < k; r++){
        // printf("r = %d\n",r);
        for(int i = 0; i < N; i++){
            if(r==0){
                x_s[r][i] = z[i];
                // printf("z[%d] = %f\t",i,z[i]);
            }
            else{
                x_s[r][i] = z[i+r*(L+Lh-1)-Lh+1]; //segment x into blocks of length L+Lh-1
                // printf("z[%d] = %f\t",i+r*(L+Lh-1)-Lh+1,z[i+r*(L+Lh-1)-Lh+1]);
            }
            // printf("x_s[%d][%d] = %f\n",r,i,x_s[r][i]);
        }
        // fwrite(x_s[r], sizeof(float), N, ff);
    }
    printf("Start FFT\n");
    complx* h = calloc(sizeof(complx), N);//does not allocate zeros!
    float* h_1 = calloc(sizeof(float), N);
    int j = 0;
    for(int i = 0; i < N; i++){
        // printf("i = %d\n",i);
        if((i < (L-1)/2) || (i > ((L-1)/2)+Lh-1)){
            h[i].re = 0;
            h[i].im = 0;
        }
        else{
            h[i].re = filter[j];
            h[i].im = 0;
            j++;
        }
        // printf("h[%d].re = %f\n",i,h[i].re);
    }
    fft842(0, N, h); //N-point dft
    for(int m = 0; m < k; m++){
        for(int i = 0; i < N; i++){
            x_s_fft[m][i].re = x_s[m][i];
            printf("x_s_fft[%d][%d].re = %f\tx_s[%d][%i] = %f\n",m,i,x_s_fft[m][i].re,m,i,x_s[m][i]);
            x_s_fft[m][i].im = 0;
        }
    }
    clock_t begin1 = clock();
    for(int m = 0; m < k; m++){ //fft of each segment
        fft842(0, N, x_s_fft[m]); 
        // fwrite(x_s_fft[m], sizeof(float), N, ff);
        // printf("m = %d\n", m);
    }
    
    for(int m = 0; m < k; m++){
        for(int i = 0; i < N; i++){
            y_s[m][i].re = x_s_fft[m][i].re*h[i].re;
            y_s[m][i].im = x_s_fft[m][i].im*h[i].im;
            printf("y_s[%d][%d].re = %f\tx_s_fft[%d][%i].re = %f\th[%d] = %f\n",m,i,y_s[m][i].re,m,i,x_s_fft[m][i].re,i,h[i].re);
        }
        fft842(1, N, y_s[m]);
        
    }
    clock_t end1 = clock();
    double time_spent1 = (double)(end1 - begin1)/CLOCKS_PER_SEC;
    printf("Time of fft = %0.20fs\n", time_spent1);
    float* y_out = calloc(sizeof(float), Lx);
    for(int m = 0; m < k; m++){
        for(int i = 0; i < L; i++){
            y_out[i+m*(L-1)] = y_s[m][i+Lh-1].re;
            // printf("m = %d\ti = %d\ty_out[%d] = %f\ty_s[%d][%d].re = %f\n",m,i,i+m*(L-1),y_out[i+m*(L-1)],m,i+Lh-1,y_s[m][i+Lh-1].re);
            // printf("i = %d\n", i+Lh-1);

            //CURRENT PROBLEM: OUTPUT DOES NOT GO TO LENGTH 4096. 


        }
    }
    fwrite(y_out, sizeof(float), Lx, ff);

    printf("Done.\n");
    fclose(fx); 
    fclose(fy); 
    fclose(ff);
}
