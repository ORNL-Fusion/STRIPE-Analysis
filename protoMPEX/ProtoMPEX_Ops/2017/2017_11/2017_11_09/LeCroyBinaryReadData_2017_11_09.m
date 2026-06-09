close all
clear all

% Signal 1
%hm = cd;
%cd('E:\LeCroy')
FileName = 'C1test00001.trc';
f = ReadLeCroyBinaryWaveform(FileName);
ty = f.x;
y = f.y;
%yh = y + 1i*paddedhilbert(y);

% Signal 2
FileName = 'C2test00001.trc'; % Reference
r = ReadLeCroyBinaryWaveform(FileName);
tr = r.x;
r = r.y;
%rh = r + 1i*paddedhilbert(r);

%cd(hm)

if 0
    figure 
    hold on
    plot(ty,real(yh),'k')
    plot(ty,imag(yh),'r')
    % xlim([40e-3,40e-3 + 40*70e-9])
    ylim([-0.1,0.1])
end
if 1 
    figure 
    hold on
    plot(ty,y,'k.-')
    plot(tr,r,'r.-')
    xlim([0.02,0.02 + 0.3e-6])
    ylim([-1,1])
end
