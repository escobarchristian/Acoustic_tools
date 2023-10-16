%% Example2 of Fast fourier transform (with chirp)

clc
clear
close all
Fs = 2000;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
x1 = cos(2*pi*15*t);          % First row : CW
x2 = cos(2*pi*450*t);         % Second row: CW
x3= chirp(t,80,0.5,350);      % Third row : chirp 80 - 350 Hz
titles_{1}='x1 = cos(15*2\pit)';
titles_{2}='x2 = cos(450*2\pit)';
titles_{3}='x3 = Chirp: 80 - 350 Hz';
titles_{4}='x4 = x1 + x2 + x3 ';

X = [x1;x2;x3;x2+x1+x3];
lim_t=500;
figure
for i = 1:4
    subplot(4,1,i)
    plot(t(1:lim_t),X(i,1:lim_t),'linewidth',1.3)
    title(titles_{i})
    set(gca,'Fontsize',13)
    ylabel('Amplitude')
%     title(['Row ',num2str(i),' in the Time Domain'])
xlim([0 t(lim_t)])
end
xlabel('Time(s)')
set(gcf, 'Units', 'Normalized', 'OuterPosition',  [  0.00 0.035 0.5 0.99])

n = 2^nextpow2(L);
dim = 2;

Y = fft(X,n,dim);
P2 = abs(Y/L);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
figure
for i=1:4
    subplot(4,1,i)
    plot(0:(Fs/n):(Fs/2-Fs/n),P1(i,1:n/2),'linewidth',1.3)
    title(['fft of: ','x',num2str(i)])
    ylim([0 1])
    set(gca,'Fontsize',13)
    xticks([0:100:1000])
    grid on
end
xlabel('Frequency')
set(gcf, 'Units', 'Normalized', 'OuterPosition',  [  0.00 0.035 0.5 0.99])

