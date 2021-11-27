#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include "filter.h"
// #include "filter_scale.h"
#include "filter_3.h"
#include "cos_16.h"
#include "cos_8.h"
#include "cos_4.h"
#include <time.h>

/*
This is for the second computer assignment in DSIP. 
Sample rate of input is 11,025 Hz
Period = T = 1/11,025 Hz = 9.07x10^-5

3 Sections: 
1. Upsample by L
2. Convolve with a lowpass
3. Downsample by M 

L/M = 3/2

Upsampling: 
insert L-1 zeros into the input signal
Period becomes T/L
Signal size becomes 3*Fx.d0

Low Pass Filter info:
->  L > M so our cutoff frequency is pi/L = pi/3
->  cutoff is Fc = 1837.5
->  Has a gain of L = 3
->  Sampling freq is Fs*(3/2) = 11025 *(3/2) = 16,537.5 Hz

https://en.wikipedia.org/wiki/Downsampling_(signal_processing)
https://en.wikipedia.org/wiki/Upsampling
*/

#define L 3
#define M 2

typedef struct 
{ 
 int ndim; //number of dimensions 
 int nchan; //number of channels 
 int d0;  //length of the first dimension 
 int d1;  //length of second dimension or sample rate if audio 
 int d2;  //length of third dimension 
} dsp_file_header; 

/**
 * @brief Upsampler
 * @param l L - our sampling constant
 * @param x is our input
 * @param xl is our output
 * @param Lzl is the length of our output (with padding)
 * @param Lx is our input length
 * @param Lh is our impulse length
 */
float* upsample(int l, float* x, float* xl, int Lzl, int Lx, int Lh){
    printf("Start Upsampling\n");
    int j = 0;
    int i;
    for(i = Lh-1; i < (Lzl-Lh-1); i++){
    // for(int i = Lh-1, j = 0; i < (Lzl-Lh-1), j < Lx; i++, j++){ 
        // printf("i = %d\tj = %d\n",i,j);
        if(i%(L) == 0){
            // printf("i = %d, j = %d\n",i,j);
            xl[i] = x[j];
            j++;
        }
    }
    printf("i = %d, j = %d\n",i,j);
    // printf("%d\n",sizeof(xl)); //i = 144175      j = 47899
    return xl;
}

/**
 * @brief Downsampler
 * @param m M - the downsampler constant
 * @param x is the input
 * @param y is the output
 * @param Lx is the length of the input
 * @param Ly is the length of our output
 */
float* downsample(int m, float* x, float* y, int Lx, int Ly){
    printf("Start Downsampling\n");
    int j = 0;
    for(int i = 0; i < Lx; i++){
    // for(int i = 0, j = 0; i < Lx, j < Ly; i++, j++){//grab every M-1 sample, store in every y[j]
        if(i%(m) == 0){
            // printf("i = %d, j = %d\n",i,j);
            y[j] = x[i];
            j++;
        }
    }
    return y;
}

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
    int j;
    for (int i = 0; i < Ly; i++) { 
        // printf("check %d\n", i);
        for (j = 0; j < Lh; j++) { 
            y[i] += h[j] * x[i + j]; //multiply and accumulate (MAC) 
            // printf("x[%d] = %f\th[%d] = %f\n",i,x[i+j],i,h[i]);
            //impulse is assumed to be in time reverse order 
        } 
        //printf("i = %d\t j = %d\n", i, j); //i = 61703        j = 479
    }
    return y;
}

/**
 * @brief Delay function?
 * @param input is input
 * @param output is output 
 * @param in_length is input length
 * @param out_length is output length
 * @param delay_val is how much the input gets delayed by
 */
float* delay(float* input, float* output, int in_length, int out_length, int delay_val){
    for(int i = 0; i < in_length; i++){
        output[i] = input[i-delay_val]; //for out example, I believe we are seperating the input into 3 different outcomes 
    }
}

void main(int argc, char** argv){
    //read in file
    FILE* fx;
    FILE* fy;
    if (NULL == (fx = fopen(argv[1], "rb"))) { //error check and open file
        printf("error: Cannot open input file.\n");
        return;
    }
    if (NULL == (fy = fopen(argv[2], "wb"))) { //error check and open file
        printf("error: Cannot open output file for writing.\n");
        return;
    }
    //grab headers of each file 
    dsp_file_header h0, h1, ho; 
    fread(&h0, sizeof(dsp_file_header), 1, fx); 
    memcpy(&ho, &h0, sizeof(dsp_file_header));
    int a = L;
    int b = M;
    float fs = (float)h0.d1;
    float fs_out = ((float)a/(float)b)*fs;
    printf("%d, %d\n",a,b);
    printf("%f\n",fs_out);
    // ho.d1 = (int)16537.5;
    ho.d1 = (int)fs_out;
    fwrite(&ho, sizeof(dsp_file_header), 1, fy); 
    printf("ndim = %d, nchan = %d, d0 = %d, d1 = %d, d2 = %d\n", h0.ndim, h0.nchan, h0.d0, h0.d1, h0.d2);

    int Lh = sizeof(h)/sizeof(h[0]); //coefficients from filter.h
    
    int Lx = h0.d0;
    int Lxl = h0.d0 * L;
    int Lz = Lxl + 2*(Lh-1); //length output for upsampled
    int Lv = Lxl + (Lh-1); //length of filtered upsampled signal 
    int Ly = Lv / M; //final output length
    printf("Lh = %d, Lx = %d, Ly = %d, Lv = %d, Lzl = %d\n", Lh, Lx, Ly, Lv, Lz); 
    float* x = calloc(sizeof(float), Lx); //allocate data for file store 
    float* v = calloc(sizeof(float), Lz);
    float* y = calloc(sizeof(float), Ly); 
    float* xl = calloc(sizeof(float), Lz); 
    
    while (!feof(fx)) {
        fread(x, sizeof(float), Lx, fx);
    }
    clock_t begin = clock(); 
    xl = upsample(L, x, xl, Lz, Lx, Lh);
    v = conv(xl, h, v, Lz, Lh, Lv);
    y = downsample(M, v, y, Lv, Lx);
    clock_t end = clock();
    double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
    printf("Time of standard form = %f\n", time_spent);
    fwrite(y, sizeof(float), Ly, fy); 
    printf("Lh = %d, Lx = %d, Ly = %d, Lv = %d, Lz = %d\n", Lh, Lx, Ly, Lv, Lz); 
    printf("ndim = %d, nchan = %d, d0 = %d, d1 = %d, d2 = %d\n", ho.ndim, ho.nchan, ho.d0, ho.d1, ho.d2);
    //------------------------------------------------------------------------------------------------------cos_16
    int Lx_16 = sizeof(x_16)/sizeof(x_16[0]);
    int Lxl_16 = Lx_16 * L;
    int Lz_16 = Lxl_16 + 2*(Lh-1); //length output for upsampled
    int Lv_16 = Lxl_16 + (Lh-1); //length of filtered upsampled signal 
    int Ly_16 = Lv_16 / M; //final output length
    float* v_16 = calloc(sizeof(float), Lz_16);
    float* y_16 = calloc(sizeof(float), Ly_16); 
    float* xl_16 = calloc(sizeof(float), Lz_16);
    printf("Lh_16 = %d, Lx_16 = %d, Ly_16 = %d, Lv_16 = %d, Lzl_16 = %d\n", Lh, Lx_16, Ly_16, Lv_16, Lz_16); 

    xl_16 = upsample(L, x_16, xl_16, Lz_16, Lx_16, Lh);
    // fwrite(xl_16, sizeof(float), Lxl_16, fy);
    v_16 = conv(xl_16, h, v_16, Lz_16, Lh, Lv_16);
    // fwrite(v_16, sizeof(float), Lv_16, fy);
    y_16 = downsample(M, v_16, y_16, Lv_16, Ly_16);
    fwrite(y_16, sizeof(float), Ly_16, fy); 
    //------------------------------------------------------------------------------------------------------cos_8
    int Lx_8 = sizeof(x_8)/sizeof(x_8[0]);
    int Lxl_8 = Lx_8 * L;
    int Lz_8 = Lxl_8 + 2*(Lh-1); //length output for upsampled
    int Lv_8 = Lxl_8 + (Lh-1); //length of filtered upsampled signal 
    int Ly_8 = Lv_8 / M; //final output length
    float* v_8 = calloc(sizeof(float), Lz_8);
    float* y_8 = calloc(sizeof(float), Ly_8); 
    float* xl_8 = calloc(sizeof(float), Lz_8);
    printf("Lh_8 = %d, Lx_8 = %d, Ly_8 = %d, Lv_8 = %d, Lzl_8 = %d\n", Lh, Lx_8, Ly_8, Lv_8, Lz_8); 

    xl_8 = upsample(L, x_8, xl_8, Lz_8, Lx_8, Lh);
    v_8 = conv(xl_8, h, v_8, Lz_8, Lh, Lv_8);
    y_8 = downsample(M, v_8, y_8, Lv_8, Ly_8);
    fwrite(y_8, sizeof(float), Ly_8, fy); 
    //------------------------------------------------------------------------------------------------------cos_4
    int Lx_4 = sizeof(x_4)/sizeof(x_4[0]);
    int Lxl_4 = Lx_4 * L;
    int Lz_4 = Lxl_4 + 2*(Lh-1); //length output for upsampled
    int Lv_4 = Lxl_4 + (Lh-1); //length of filtered upsampled signal 
    int Ly_4 = Lv_4 / M; //final output length
    float* v_4 = calloc(sizeof(float), Lz_4);
    float* y_4 = calloc(sizeof(float), Ly_4); 
    float* xl_4 = calloc(sizeof(float), Lz_4);
    printf("Lh_4 = %d, Lx_4 = %d, Ly_4 = %d, Lv_4 = %d, Lzl_4 = %d\n", Lh, Lx_4, Ly_4, Lv_4, Lz_4); 

    xl_4 = upsample(L, x_4, xl_4, Lz_4, Lx_4, Lh);
    v_4 = conv(xl_4, h, v_4, Lz_4, Lh, Lv_4);
    printf("v_4[10] = %f\n", v_4[10]);
    y_4 = downsample(M, v_4, y_4, Lv_4, Ly_4);
    printf("y_4[10] = %f\n", y_4[10]);
    fwrite(y_4, sizeof(float), Ly_4, fy); 
    //------------------------------------------------------------------------------------------------------
    // //output to file
    fclose(fx); 
    fclose(fy); 
}



/*
Matlab code to generate cosine
Lh = 479;
N = 100*Lh;
n = 0:1:N-1;
fo = 1/16; %fo = 1/8; %fo = 1/4;
x = cos(2*pi*fo*n);
plot(x)
fid = fopen('cos_in.txt', 'w');
fprintf(fid,'%0.10f,\r\n', x);
fclose(fid);

fid = fopen('cos_in.txt','r');
x = textscan(fid, '%f');
fclose(fid);
plot(x{1})


Matlab code to plot outputted cosine
fid = fopen('cos_16_out.bin');
ndim  = fread(fid,1,'int');
nchan = fread(fid,1,'int');
dim0  = fread(fid,1,'int');
dim1  = fread(fid,1,'int');
dim2  = fread(fid,1,'int');
[x,fs] = fread(fid,inf,'float');
fclose(fid);
plot(x);
*/