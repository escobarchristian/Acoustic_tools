%% Example1 of Fast fourier transform with filter

clc
clear
close all
Fs = 2000;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
x1 = cos(2*pi*50*t);          % First row wave
x2 = cos(2*pi*150*t);         % Second row wave
x3 = cos(2*pi*300*t);         % Third row wave
titles_{1}='x1 = cos(50*2\pit)';
titles_{2}='x2 = cos(150*2\pit)';
titles_{3}='x3 = cos(300*2\pit)';
titles_{4}='x4 = x1 + x2 + x3';
X = [x1; x2; x3;x2+x3+x1];
figure
for i = 1:4
    subplot(4,1,i)
    plot(t(1:100),X(i,1:100),'linewidth',1.3)
    title(titles_{i})
    set(gca,'Fontsize',13)
    ylabel('Amplitude')
%     title(['Row ',num2str(i),' in the Time Domain'])
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
%     title(['Row ',num2str(i),' in the Frequency Domain'])
end
xlabel('Frequency')
set(gcf, 'Units', 'Normalized', 'OuterPosition',  [  0.00 0.035 0.5 0.99])

%% Filter signals
fc = 70; % Cut off frequency
n = 3; % order
Wn = fc/(Fs/2); % Get frequency 
[b,a] = butter(n, Wn, 'low'); % Generate Butterworth filter
filteredSignal = filter(b, a, X(4,:)); % Filter signal
%% Plot filtered signal
figure
plot(t,X(4,:),'k-','Linewidth',0.5)
hold on
plot(t,filteredSignal,'b-','Linewidth',1.2)
xlim([0 0.15])
legend('X4=X1+X2+X3','Low pass Filtered signal at 70Hz')
set(gca,'Fontsize',13)
ylabel('Amplitude')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.125, 0.8, 0.7]);
%% Plot frequency response
figure
freqz(b,a,[],Fs)
subplot(2,1,1)
hold on
plot([70 70],[-210 0],'r--')
set(gca,'Fontsize',13)

subplot(2,1,2)
hold on
plot([70 70],[-300 0],'r--')
set(gca,'Fontsize',13)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.125, 0.8, 0.7]);