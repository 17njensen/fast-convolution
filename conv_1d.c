#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
 
typedef struct 
{ 
 int ndim; //number of dimensions 
 int nchan; //number of channels 
 int d0;  //length of the first dimension 
 int d1;  //length of second dimension or sample rate if audio 
 int d2;  //length of third dimension 
} dsp_file_header; 
 
int conv() { 
    FILE* fh, * fx, * fy; 
    if (NULL == (fh = fopen("lpf_260_400_44100_80dB.bin", "rb"))) 
    { 
        perror(fh); 
        printf("error: Can't open lpf_260_400_44100_80dB.bin for input.\n"); 
        return 0; 
    } 
    if (NULL == (fx = fopen("fireflyintro.bin", "rb"))) 
    { 
        printf("ERROR: Can't open fireflyintro.bin for input.\n"); 
        return 0; 
    } 
    if (NULL == (fy = fopen("conv_firefly.bin", "wb"))) 
    { 
        printf("ERROR: cant open fstereo.bin for output.\n"); 
        return 0; 
    } 
    //grab headers of each file 
    dsp_file_header h0, h1, ho; 
    fread(&h0, sizeof(dsp_file_header), 1, fx); 
    fread(&h1, sizeof(dsp_file_header), 1, fh); 
    memcpy(&ho, &h0, sizeof(dsp_file_header)); 
    fwrite(&ho, sizeof(dsp_file_header), 1, fy); 
 
    int Lh = h1.d0;//length of impulse signal 
    int Lx = h0.d0; //length of input signal 
    int Ly = Lx + (Lh - 1); //len of conv result 
    int Lz = Lx + 2 * (Lh - 1); //len of zero padded input 
    printf("Lh = %d, Lx = %d, Ly = %d, Lz = %d\n", Lh, Lx, Ly, Lz); 
 
    float* x = calloc(sizeof(float), Lz); //allocate data for file store 
    float* h = calloc(sizeof(float), Lz);  
    float* y = calloc(sizeof(float), Ly); 
  
    while (!feof(fx)) { 
        fread((x + Lh - 1), sizeof(float), Lx, fx); //pulls in data when H0 is present, hence x+Lh-1 
    } 
    while (!feof(fh)) { 
        fread(h, sizeof(float), Lh, fh); 
    } 
 
    int i, j; 
    for (i = 0; i < Ly; i++) { 
        for (j = 0; j < Lh; j++) { 
            y[i] += h[j] * x[i + j]; //multiply and accumulate (MAC) 
            //impulse is assumed to be in time reverse order 
        } 
    } 
    fwrite(y, sizeof(float), Ly, fy); 
    fclose(fx); 
    fclose(fh); 
    fclose(fy); 
    return 0; 
} 