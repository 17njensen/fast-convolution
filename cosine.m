Lh = 256;
N = 10*Lh;
n = 0:0.0001:N-1;
%n = 0:1:N-1;
%fo = 1/16;
%fo = 10;
%fo = 40; 
%fo = 150; 
%fo = 350; 
fo = 500;
x = cos(2*pi*fo*n);
%plot(n,x)
fid = fopen('cos_10.txt', 'w');
% ndim  = 1;
% nchan = 1;
% dim0  = N;
% dim1  = 11025;
% dim2  = 0;
%fwrite(fid,[1 size(x) 11025 0], 'int');
%fwrite(fid,x(:),'float');
%fprintf(fid, '%d\n%d\n%d\n%d\n%d', ndim, nchan, dim0, dim1, dim2);
%fprintf(fid,'%0.10f,\r\n', x);
fclose(fid);

nfft = 2^9;
nvals = (0:1:nfft-1)/nfft;
fid1 = fopen('cos_500_test.bin','r');
% ndim1  = fread(fid1,1,'int');
% nchan1 = fread(fid1,1,'int');
% dim01  = fread(fid1,1,'int');
% dim11  = fread(fid1,1,'int');
% dim21  = fread(fid1,1,'int');
[x_in,fs_in] = fread(fid1,inf,'float');
fclose(fid1);
%plot(x_in)
plot(nvals,20*log10(abs(fftshift(fft(x_in,nfft)))));
hold on;
plot(nvals,20*log10(abs(fftshift(fft(x,nfft)))));
%plot(fft(x));
%plot(x_in)
hold off;
xlabel('Frequency(Hz)', 'FontSize', 10);
ylabel('Amplitude(dB)','FontSize', 10);
title('Cosine Frequency Response','FontSize', 10);
legend('Output', 'Input of 500')
whos
