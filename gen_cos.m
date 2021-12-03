Lh = 256;
N = 4.095; %to get a power of 2
%N = 0.4095; %for 500
n = 0:0.001:N;
%n = 0:0.0001:N; %for 500
%fo = 1/16;
fo = 10;
%fo = 40; 
%fo = 150; 
%fo = 350; 
%fo = 500;
x = cos(2*pi*fo*n);
nfft = 2^9;
nvals = (0:1:nfft-1)/nfft;
%nvals = (0:1:nfft-1);
plot(nvals,20*log10(abs(fftshift(fft(x,nfft)))));
%plot(n,x)
s1 = 'cos_';
s2 = string(fo);
s3 = '.bin';
s_1 = strcat(s1,s2);
s = strcat(s_1,s3);
fid = fopen(s, 'w');
fwrite(fid,[1 size(x) 11025 0], 'int');
fwrite(fid,x(:),'float');
fclose(fid);