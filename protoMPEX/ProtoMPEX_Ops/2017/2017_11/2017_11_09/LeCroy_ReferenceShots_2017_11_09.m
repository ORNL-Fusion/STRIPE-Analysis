close all
clear all

shot = 17483;
shot = 17484;

switch shot
    case 17483 % ==========================================================
% Reference signal, no demodulator
FileName_1 = 'C1_00000.trc';
f1 = ReadLeCroyBinaryWaveform(FileName_1);
ty1 = f1.x;
y1 = f1.y;
% Bdot probe signal
FileName_2 = 'C2_00000.trc';
f2 = ReadLeCroyBinaryWaveform(FileName_2);
ty2 = f2.x;
y2 = f2.y;

figure 
hold on
rng = find(ty1>=20e-3 & ty1<=50e-3+0.3e-6);
rng = find(ty1>=-0.1e-3 & ty1<=40e-3+0.3e-6);
plot(ty1(rng),y1(rng),'k-')
plot(ty2(rng),y2(rng),'r-')
ylim([-1,1])

y1h = y1(rng) + 1i*paddedhilbert(y1(rng));
y2h = y2(rng) + 1i*paddedhilbert(y2(rng));

figure; hold on
% plot(unwrap((1*angle(y1h)-0*angle(y2h))))
% plot(unwrap((0*angle(y1h)+1*angle(y2h))),'r')
plot( (angle(y1h) - angle(y2h))/(2*pi), 'k.' )
grid on

rng = find(ty1>=20e-3 & ty1<=50e-3+0.3e-6);
[freq1,Fy1,~] = myfft(ty1(rng),y1(rng),'0');
[freq2,Fy2,~] = myfft(ty2(rng),y2(rng),'0');

figure; hold on
plot(freq2,abs(Fy2)/max(y2),'r')
plot(freq1,abs(Fy1)/max(y2))
xlim([0,30e6])
set(gca,'Yscale','log')
ylim([1e-6,1e-1])

    case 17484
 % Reference signal, no demodulator
FileName_1 = 'C1_00001.trc';
f1 = ReadLeCroyBinaryWaveform(FileName_1);
ty1 = f1.x;
y1 = f1.y;
% Bdot probe signal
FileName_2 = 'C1_00000.trc';
f2 = ReadLeCroyBinaryWaveform(FileName_2);
ty2 = f2.x;
y2 = f2.y;

figure 
hold on
rng = find(ty1>=20e-3 & ty1<=50e-3+0.3e-6);
rng = find(ty1>=0e-3 & ty1<=40e-3+0.3e-6);
plot(ty1(rng),y1(rng),'k-')
plot(ty2(rng),y2(rng),'r-')
ylim([-1,1])

y1h = y1(rng) + 1i*paddedhilbert(y1(rng));
y2h = y2(rng) + 1i*paddedhilbert(y2(rng));

figure; hold on
% plot(unwrap((1*angle(y1h)-0*angle(y2h))))
% plot(unwrap((0*angle(y1h)+1*angle(y2h))),'r')
plot( unwrap(angle(y1h) - angle(y2h))/(2*pi) ,'k.')
grid on
ylim([-pi,pi])
        
end

