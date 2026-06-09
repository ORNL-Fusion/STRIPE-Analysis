close all
clearvars

% Signal 1
% hm = cd;
% cd('E:/LeCroy')
Path = '/Users/pawpiotr/Box Sync/ORNL/Proto-MPEX/Helicon/Experimental Data/3wave/LeCroy/';
FileName = 'C2_00001';
f = ReadLeCroyBinaryWaveform(sprintf('%s%s.trc',Path,FileName));
t = f.x;
v = f.y;

figure(1)
plot(t,v)

Fs = 1/mean(diff(t));

tic
% [s,freq,time] = spectrogram(v,200,10,10000,Fs);
% [q,nd] = max(10*log10(p));
[omega,signal] = FourierTransform(v,t);
toc

clear v t f
subplot(2,1,1)
plot(omega./2./pi,log10(signal/max(signal)))
axis([0 50e6 -6 0])
set(gca,'fontsize',18,'fontname','times')
xlabel('Frequeny [Hz]')
ylabel('log_1_0(Normalized Amplitude)')
subplot(2,1,2)
plot(omega./2./pi,log10(signal/max(signal)))
axis([5e6 20e6 -6 0])
set(gca,'fontsize',18,'fontname','times')
xlabel('Frequeny [Hz]')
ylabel('log_1_0(Normalized Amplitude)')

