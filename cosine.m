Lh = 256;
N = 4.095; %to get a power of 2
n = 0:0.001:N;

fid = fopen('cos_500.bin','r');
ndim1  = fread(fid1,1,'int');
nchan1 = fread(fid1,1,'int');
dim01  = fread(fid1,1,'int');
dim11  = fread(fid1,1,'int');
dim21  = fread(fid1,1,'int');
[x,fs] = fread(fid,inf,'float');
fclose(fid);
nfft = 2^9;
nvals = (0:1:nfft-1)/nfft;

fid1 = fopen('cos_500_conv.bin','r');
ndim  = fread(fid1,1,'int');
nchan = fread(fid1,1,'int');
dim0  = fread(fid1,1,'int');
dim1  = fread(fid1,1,'int');
dim2  = fread(fid1,1,'int');
[x_in,fs_in] = fread(fid1,inf,'float');
fclose(fid1);

fid2 = fopen('cos_500_fft.bin','r');
ndim  = fread(fid1,1,'int');
nchan = fread(fid1,1,'int');
dim0  = fread(fid1,1,'int');
dim1  = fread(fid1,1,'int');
dim2  = fread(fid1,1,'int');
[x_fft,fs_fft] = fread(fid2,inf,'float');
fclose(fid2);

plot(nvals,20*log10(abs(fftshift(fft(x_in,nfft)))));
%plot(x_in);
hold on;
plot(nvals,20*log10(abs(fftshift(fft(x,nfft)))));
%plot(x);
hold off;
hold on;
plot(nvals,20*log10(abs(fftshift(fft(x_fft,nfft)))));
%plot(x_fft);
hold off;
xlabel('Frequency(Hz)', 'FontSize', 10);
ylabel('Amplitude(dB)','FontSize', 10);
title('Cosine Frequency Response','FontSize', 10);
legend('Conv', 'Input of 10', 'fft842 result')
whos