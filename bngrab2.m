clc
close all
clear all

% plotting parameters
fplot = 60; % in KHz
tplot = 1/fplot;
tmin = -40;
tmax = 40;
t = tmin:tplot:tmax;
f = linspace(-fplot/2, fplot/2, length(t));

% message signal
fm = 1; % in KHz
Am = 5;
mt = Am.*cos(2*pi*fm*t);
mf = (1./length(t)).*(fftshift(fft(mt)));

% time domain sampling
% sampling at Fs = 10 KHz
Fs_10 = 10;
Ts_10 = 1/Fs_10;
nmin10 = ceil(tmin/Ts_10);
nmax10 = floor(tmax/Ts_10);
n10 = nmin10:nmax10;
sampled_10k = Am.*cos(2*pi*fm*n10*Ts_10);

% sampling at Fs = 2 KHz
Fs_2 = 2;
Ts_2 = 1/Fs_2;
nmin2 = ceil(tmin/Ts_2);
nmax2 = floor(tmax/Ts_2);
n2 = nmin2:nmax2;
sampled_2k = Am.*cos(2*pi*fm*n2*Ts_2);

% sampling at Fs = 1 KHz
Fs_1 = 1;
Ts_1 = 1/Fs_1;
nmin1 = ceil(tmin/Ts_1);
nmax1 = floor(tmax/Ts_1);
n1 = nmin1:nmax1;
sampled_1k = Am.*cos(2*pi*fm*n1*Ts_1);

% frequency domain sampling
% sampling at Fs = 10 KHz
y10 = zeros(size(f));
indexes10 = find(mod(f, 10) == 0);
y10(indexes10) = 1;
s10 = conv(y10, mf);

% sampling at Fs = 2 KHz
y2 = zeros(size(f));
indexes2 = find(mod(f, 2) == 0);
y2(indexes2) = 1;
s2 = conv(y2, mf);

% sampling at Fs = 1 KHz
y1 = zeros(size(f));
indexes1 = find(mod(f, 1) == 0);
y1(indexes1) = 1;
s1 = conv(y1, mf);

% plotting
%1. message waveform and spectrum
figure('Name','Message signal','NumberTitle','off')
subplot(2,1,1)
plot(t,mt)
xlim([-5/fm 5/fm])
title("Message waveform")
xlabel("time(msec)")
ylabel("Amplitude")

subplot(2,1,2)
stem(f,abs(mf))
xlim([-2.*fm 2.*fm])
ylim([0 3])
title('Message spectrum')
xlabel('frequency(kHz)')
ylabel('Magintude')

%2. Time Domain representatoin
% Sampled at 10 kHz
figure('Name','Time Domain Representation','NumberTitle','off') 
subplot(3,1,1)
hold on
plot(t,mt)
plot(n10*Ts_10, sampled_10k, '.')
hold off
xlim([-5/fm 5/fm])
title("Sampling at 10 kHz")
xlabel("time(msec)")
ylabel("Amplitude")

% Sampled at 2 kHz
subplot(3,1,2)
hold on
plot(t,mt)
plot(n2*Ts_2, sampled_2k, '.')
hold off
xlim([-5/fm 5/fm])
title("Sampling at 2 kHz")
xlabel("time(msec)")
ylabel("Amplitude")

% Sampled at 2 kHz
subplot(3,1,3)
hold on
plot(t,mt)
plot(n1*Ts_1, sampled_1k, '.')
hold off
xlim([-5/fm 5/fm])
title("Sampling at 1 kHz")
xlabel("time(msec)")
ylabel("Amplitude")

%3. Frequency Domain representation
% Sampled at 10 KHz
figure('Name','Frequency Domain Representation','NumberTitle','off') 
subplot(3,1,1)
stem(linspace(-fplot,fplot,length(s10)),abs(s10))
xlim([-20 20])
title("Sampling at 10 kHz")
xlabel("frequency(KHz)")
ylabel("Magnitude")

subplot(3,1,2)
stem(linspace(-fplot,fplot,length(s2)),abs(s2))
xlim([-20 20])
title("Sampling at 2 kHz")
xlabel("frequency(KHz)")
ylabel("Magnitude")

subplot(3,1,3)
stem(linspace(-fplot,fplot,length(s1)),abs(s1))
xlim([-20 20])
title("Sampling at 1 kHz")
xlabel("frequency(KHz)")
ylabel("Magnitude")