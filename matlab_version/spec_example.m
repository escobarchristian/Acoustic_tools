%% Spectrogram Example
% Example to use the spectrogram function

clc
clear
close all

read_audio=true; % true to read in a wav file or false to generate a sample signal

if read_audio
    % You can just read the audiofile as
    [x,fs]=audioread("chirp_example.wav");
    t = (0:length(x)-1) /fs;
else
    % Generate a sample signal
    fs=8000; % in Hz, sampling frequency
    duration= 10; % Duration in seconds
    N = fs * duration; % Number of samples
    t = (0:N) /fs;

    % Example
    % Generate a linear frequency modulated signal from 100 to 2000 Hz with a duration of 50 seconds
    f0=100; % Freq at time 0
    f1=2000; % Freq at time t1
    t1=10; % time in seconds
    x = chirp(t,f0,t1,f1,'linear'); % Linear frequency modulated signal.
    % Add some noise
    mu=0; % Mean
    sigma= 0.2; % std
    noise = mu + sigma*randn(size(x));
    x=x+noise;

    % Add a tone
    tone_freq=3500.0; % in Hz
    x2=sin(2*pi*tone_freq.*t);
    x=x+x2;
    x = x./(max(abs(x))*1.001);
end
% Calculate spectrogram using the time series x and sampling rate fs
nfft=2^12; % number of samples for each FFT (2^12=4096). df=fs/nfft
tWindow = 0.1; % Integration time in seconds, ideally it is nfft/fs
NWindow = fs*tWindow; % Window size, i.e., number of points to be considered in the fft
window = hamming(NWindow); % Apply the window you would like to use.
overlap=0.5; % Overlap (0.0 - 0.99).Typical value: 0.5
NOverlap = floor(NWindow*overlap); % Number of Overlap points

[~,F,T,P] =spectrogram(x,window,NOverlap,nfft,fs,'yaxis');
spec=10*log10(P);

file_root = 'chirp';
folder_name='spectrogram_figs/';
% Create folder if it does not exist
if ~exist(folder_name, 'dir')
    mkdir(folder_name)
end
% Plot time series
figure
plot(t,x)
xlabel('Time (s)')
ylabel('Amplitude')
title('Time series')

% saveas(gcf,[results_folder1,fname_rec(1:end-8),'P',num2str(i_p),'.jpg'],'jpg')
exportgraphics(gca,[folder_name, file_root, '_time_series.png'],'Resolution',150)

imagesc(T,F,10*log10(P))
caxis([-100 -20 ]) %%%%%% Make sure you modify these values of dynamic
% range for the spectrogram plot based on your intensity levels.
h=colorbar;
ylabel(h,'Power Spectral Density (dB)','fontsize',12)
xlabel("Time (sec)")
ylabel("Frequency (Hz)")
title('Spectrogram')

exportgraphics(gca,[folder_name,file_root,'_spec.png'],'Resolution',150)

%% Save audiofile
filename='chirp_example.wav';
x = x./(max(abs(x))*1.001);
audiowrite(filename,x,fs)




