function audio2bin(infile,outfile)
if(nargin == 1)
    outfile = ['test_',infile(1:max(strfind(infile, '.'))), 'bin'];
end
fprintf('Input file = %s\n', infile);
fprintf('Output file = %s\n', outfile);

%fft
[x,fs] = audioread(infile);
%for our signals in lab 4
%sound(x,fs);
time = [0:(length(x)-1)]/fs;
t1 = 1;
t2 = 1.01;
i1 = round(t1*fs);
i2 = round(t2*fs);
nfft = 2^20;
freq = ([0:nfft-1]/nfft - 0.5)*fs;
X = fft(x(i1:i2),nfft);
fmax = max(X(:))
%X = lowpass(X,500,fs);
%spectrogram
[a,fs1]= audioread(infile);
nfft2 = 2^12;
overlap = round(0.8*nfft2);
window = hamming(nfft2);
fprintf('sample rate = %d\n', fs);

%file writing
x = x.';
sx = size(x);
y = [1 size(x) fs 0];
fid = fopen(outfile, 'wb');
fwrite(fid,[1 size(x) fs 0], 'int');
fwrite(fid,x(:),'float');
fclose(fid);

% %oscilloscope
% plot(time,x);
% title('Oscilloscope of flute22.wav');
% win_sec = 0.01;
% win_sam = round(win_sec*fs);
% step_sec = 0.001;
% step_sam = round(step_sec*fs);
% han = plot(time(1:win_sam),x(1:win_sam)); drawnow;
% ylim([-1 1]);
% for i=win_sam:step_sam:length(x)
%     ind = [i-win_sam+1:i];
%     set(han,'XData',time(ind),'YData',x(ind));
%     xlim(time(ind([1,end])));
%     drawnow;
%     pause(0.01);
% end

%Plotting

%subplot(2,2,[3,4]);
% plot(time,x)
% xlabel('Time(seconds)', 'FontSize', 10);
% ylabel('Amplitude','FontSize', 10);
% title('impulseclap.wav','FontSize', 10);
% %xlim([0 10.5]);
% %ylim([-0.7 0.3]);
% grid on;
%soundsc(x)

%zoomed plot
% subplot(2,2,2);
% plot(time,x)
% xlim([1 1.01])
% xlabel('Time(seconds)', 'FontSize', 10);
% ylabel('Amplitude','FontSize', 10);
% grid on;

% FFT plot
%subplot(2,2,[3 4]);
plot(freq,20*log10(abs(fftshift(X))));
xlabel('frequency [Hz]','FontSize',10);
ylabel('magnitude [dB]','FontSize',10);
title('Frequency plot of note.wav - Problem 2','FontSize', 10);
xlim([0 500]);
grid on;

%spectrogram plot
%subplot(2,2,4);
% spectrogram(a(:),window,overlap,nfft2,fs1);
% spectrogram(a(:), 'yaxis');
% ylabel('time','FontSize',10);
% xlabel('frequency (Hz)','FontSize',10);
% title('voicerecording.wav','FontSize', 10);
% grid on;

whos
return;