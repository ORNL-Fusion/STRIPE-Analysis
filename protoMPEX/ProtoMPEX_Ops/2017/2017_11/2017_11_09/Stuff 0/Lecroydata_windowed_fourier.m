close all
clearvars

% Signal 1
% hm = cd;
% cd('E:/LeCroy')
Path = '/Users/pawpiotr/Box Sync/ORNL/Proto-MPEX/Helicon/Experimental Data/3wave/LeCroy/';
FileName = 'C2_00004';
f = ReadLeCroyBinaryWaveform(sprintf('%s%s.trc',Path,FileName));
t = f.x;
v = f.y;

figure(1)
plot(t,v)

Fs = 1/mean(diff(t));

tic
[s,freq,time] = spectrogram(v,1000,100,10000,Fs);
% [q,nd] = max(10*log10(p));
% [omega,signal] = FourierTransform(v,t);
toc

clear v t f

tmax = 20e-3;
s    = s(:,time<tmax);
time = time(time<tmax);

figure(2)
hold on
A = max(max(abs(s)));
pcolor(time,freq,log10(abs(s)./A))
axis([min(time) max(time) min(freq) 30e6 ])

shading interp
mean(diff(time))

set(gca,'fontsize',18,'fontname','times')
xlabel('Frequeny [Hz]')
ylabel('log_1_0(Normalized Amplitude)')

print(figure(2),sprintf('%s',FileName),'-dpng')

