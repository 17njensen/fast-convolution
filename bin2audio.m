function bin2audio(infile,outfile)
if(nargin == 1)
    outfile = ['NEW_',infile(1:max(strfind(infile, '.'))), 'wav'];
end
fprintf('Input file = %s\n', infile);
%fprintf('Output file = %s\n', outfile);

%[y,fsy] = audioread(infile2);

%read file and gather data for conversion
fid = fopen(infile, 'rb');
ndim  = fread(fid,1,'int');
nchan = fread(fid,1,'int');
dim0  = fread(fid,1,'int');
dim1  = fread(fid,1,'int');
dim2  = fread(fid,1,'int');
[x,fs] = fread(fid,inf,'float');
fclose(fid);
%second file
%fid1 = fopen(infile1, 'rb');
%ndim1  = fread(fid1,1,'int');
%nchan1 = fread(fid1,1,'int');
%dim01  = fread(fid1,1,'int');
%dim11  = fread(fid1,1,'int');
%dim21  = fread(fid1,1,'int');
%[x1,fs1] = fread(fid1,inf,'float');
%fclose(fid1);

%adjust data to fit audio file wave format
fprintf('ndim = %d  nchan = %d  dim0 = %d, dim1 = %d, dim2 = %d\n',ndim,nchan,dim0,dim1,dim2);
x = reshape(x,1,[]).';
fs = dim1;
x = x/max(abs(x(:)));
fprintf('sample rate = %d\n', fs);
%second file
%fprintf('ndim = %d  nchan = %d  dim0 = %d, dim1 = %d, dim2 = %d\n',ndim1,nchan1,dim01,dim11,dim21);
%x1 = reshape(x1,1,[]).';
%fs1 = dim11;
%x1 = x1/max(abs(x1(:)));
%fprintf('sample rate 2 = %d\n', fs1);

%fft
%[x,fs] = audioread(infile);
time = [0:(length(x)-1)]/fs;
t1 = 1;
t2 = 1.01;
i1 = round(t1*fs);
i2 = round(t2*fs);
nfft = 2^12;
freq = ([0:nfft-1]/nfft - 0.5)*fs;
X = fft(x(i1:i2),nfft);
%second fft data
%X1 = fft(x1(i1:i2),nfft);


%spectrogram
%[a,fs1]= audioread(infile);
% nfft2 = 2^8;
% overlap = round(0.8*nfft2);
% window = hamming(nfft2);

%plot and output
soundsc(x,fs)
%subplot(2,2,2);
a = 0:length(x)-1;
%plot(a*(10^-2), x)
%xlabel('Time(seconds)', 'FontSize', 10);
%ylabel('Normalized Amplitude','FontSize', 10);
%title('voice recording filter result','FontSize', 10);
%sound(x,fs);

%subplot(2,2,3);
%zplane(x);
 
%plot(freq,20*log10(abs(fftshift(X1))));
 hold on;
 plot(freq,20*log10(abs(fftshift(X))));
 hold off;
 %hold on;
 %plot(freq,20*log10(abs(fftshift(X2))));
 %hold off;
 %hold on;
 %plot(freq,20*log10(abs(fftshift(X3))));
 %hold off;
 xlabel('frequency [Hz]','FontSize',10);
 ylabel('magnitude [dB]','FontSize',10);
 title('Bandpass filter for fireflyintro.wav','FontSize', 10);
 legend('Original','Bandpass');
 xlim([0 22050]);
 %grid on;

% %subplot(2,2,4);
%spectrogram(x(:),window,overlap,nfft,fs);
spectrogram(x(:), 'yaxis');
ylabel('time','FontSize',10);
xlabel('frequency (Hz)','FontSize',10);
title('Original Spectrogram of fireflies.wav','FontSize', 10);

zplane(X);
title('Zplane of Bandpass filter');

audiowrite(outfile,x,fs);
whos
return;