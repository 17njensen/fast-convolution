#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "filter.h"
#include <time.h>
#include "cos_500.h"
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

// Function to reverse elements of an array
void reverse(float arr[], int n)
{
    float aux[n];
 
    for (int i = 0; i < n; i++) {
        aux[n - 1 - i] = arr[i];
    }
 
    for (int i = 0; i < n; i++) {
        arr[i] = aux[i];
    }
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
    dsp_file_header h0, ho; 
    // fread(&h0, sizeof(dsp_file_header), 1, fx); 
    memcpy(&ho, &h0, sizeof(dsp_file_header));
    //fwrite(&ho, sizeof(dsp_file_header), 1, fy); 
    //printf("ndim = %d, nchan = %d, d0 = %d, d1 = %d, d2 = %d\n", h0.ndim, h0.nchan, h0.d0, h0.d1, h0.d2);
    int Lh = sizeof(filter)/sizeof(filter[0]);//length of impulse signal 
    int Lx = sizeof(cosx)/sizeof(cosx[0]);
    // int Lx = h0.d0; //length of input signal 
    int Ly = Lx + (Lh - 1); //len of conv result 
    int Lz = Lx + 2 * (Lh - 1); //len of zero padded input 
    printf("Lh = %d, Lx = %d, Ly = %d, Lz = %d\n", Lh, Lx, Ly, Lz); 
    //allocate data space
    float* x = calloc(sizeof(float), Lz);
    float* y = calloc(sizeof(float), Ly);
    for(int i = 0; i < Lx; i++){ //add zero padding
        x[i+Lh-1] = cosx[i];
    }
    // while (!feof(fx)) { 
    //     fread((x + Lh - 1), sizeof(float), Lx, fx); //pulls in data when H0 is present, hence x+Lh-1 
    // } 
    //reverse(filter, Lh);
    //double check that the input are placed in time reverse order (flip and truncate for convolution)
    y = conv(x, filter, y, Lz, Lh, Ly);
    printf("Done.\n");
    printf("y[10] = %f\n", y[10]);
    fwrite(y, sizeof(float), Ly, fy); 
    fclose(fx); 
    fclose(fy); 
}
