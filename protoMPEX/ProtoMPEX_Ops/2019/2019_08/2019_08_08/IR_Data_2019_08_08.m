close all 
clear all


% NOTES:
% - Determine the best shot to do determine the emissivity calibration
% - For a shot, need to create the WinArea variable

 
ViewType = 3;

switch ViewType
    case 1 
    shot(1) = 27377; 
    TempThreshold = 34.4;

    case 2 % 164 kW net but tripped
    shot(1) = 27372;
    TempThreshold = 31;

    case 3 
    shot(1) = 27383;
    TempThreshold = 35;

    case 4 % 3nd 500 ms pulse at 102 kW net
    shot(1) = 26635;
    TempThreshold = 38;
    
    case 5 % 330 ms pulse at 124 kW net
    shot(1) = 26637;
    TempThreshold = 38;
    
    case 6 % 330 ms pulse at 138 kW net
    shot(1) = 26639;
    TempThreshold = 40;
    
    case 7 % 330 ms pulse at 138 kW net
    shot(1) = 26641;
    TempThreshold = 39;

end
        
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_08_08';

mov = 1;

for s = 1:length(shot)
PATHNAME = [a{s},'\'];
try
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
catch
end
% FILENAME = ['Shot-0',num2str(shot(s)),'.seq'];    
videoFileName=[PATHNAME FILENAME];


% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file'
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
seq = file.ThermalSequencePlayer();
%Get the pixels
img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
RawData{s}(:,:,1) = im;
fr = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
        RawData{s}(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;
    end
end
end

% =========================================================================
%  Calculate surface temperature
% =========================================================================
emissivity = 0.8;
[TempMeasured,dTMeasured] = IntensityTempConv2(emissivity,0.2*emissivity,RawData{1},seq);
% framerate = 100;
framerate = seq.FrameRate;
dt = 1/framerate;
time = 0:dt:(fr-2)*dt;

% =========================================================================
%  Gather FP window surface temperature
% =========================================================================
b = 0;
% To convert voltage to temperature, multiply by 20
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

% Gather RawData
DA{1} = [RootAddress,'FLUOROPT_1']; 
DA{2} = [RootAddress,'FLUOROPT_2']; 
DA{3} = [RootAddress,'FLUOROPT_3']; 
DA{4} = [RootAddress,'FLUOROPT_4']; 

[f1,t_f1]   = my_mdsvalue_v2(shot(1),DA(1));
[f2,t_f2]   = my_mdsvalue_v2(shot(1),DA(2));
[f3,t_f3]   = my_mdsvalue_v2(shot(1),DA(3));
[f4,t_f4]   = my_mdsvalue_v2(shot(1),DA(4));

% =========================================================================
% Determine the start of the RF pulse:
% =========================================================================
Nx = size(TempMeasured(:,:,:),1);
Ny = size(TempMeasured(:,:,:),2);
Nz = size(TempMeasured(:,:,:),3);

PixelCount = Nx*Ny;

for s = 1:Nz
   MeanTempMeasured(s) = sum(sum(TempMeasured(:,:,s)))/PixelCount; 
end

[~,n1] = max(diff(MeanTempMeasured));
[~,n2] = min(diff(MeanTempMeasured));

rng = [n1-10:n2+60];
D0 = TempMeasured(:,:,n1-10);
D1 = TempMeasured(:,:,n2+100);

% D0 = TempMeasured(:,:,b+20);

figure; 
subplot(1,2,1)
hold on
plot(time,MeanTempMeasured,'k.')
plot(time(rng),MeanTempMeasured(rng),'g')
title('Sum')
ylabel('surface T [C]')
xlabel('t [s]')
xlim([3,7])

subplot(1,2,2)
hold on
plot(time(1:end-1),diff(MeanTempMeasured),'k.')
plot(time(rng(1:end-1)),diff(MeanTempMeasured(rng)),'g')
title('dSum')
ylabel('surface T [C]')
xlabel('t [s]')
xlim([3,7])

%% Image prior to RF

WinArea = double(D0>TempThreshold);
NotWinArea = double(D0<TempThreshold);

figure;
hold on
surf(D0,'LineStyle','none')

figure;
hold on
surf(D0.*WinArea,'LineStyle','none')
% surf(D0,'LineStyle','none')

view([0,90])
caxis([31,40])
colormap('hot')
title(['shot: ',num2str(shot(1))])
set(gca,'PlotBoxAspectRatio',[1 0.5 1])
colormap(flipud('hot'))
colorbar
% zlim([31,33])
% caxis([31,33])

% Window temperature
C = {'k','bl','r','g'};
figure; hold on
b = 0;
plot(t_f1{1}(1:end-1),(f1{1}-b*min(f1{1}))*20,C{1},'LineWidth',2)
plot(t_f2{1}(1:end-1),(f2{1}-b*min(f2{1}))*20,C{2},'LineWidth',2) % Ground Side
plot(t_f3{1}(1:end-1),(f3{1}-b*min(f3{1}))*20,C{3},'LineWidth',2) % High Voltage side
plot(t_f4{1}(1:end-1),(f4{1}-b*min(f4{1}))*20,C{4},'LineWidth',2)

ylim([0,100])
ylabel('${\Delta}T$ $[C]$','Interpreter','latex','FontSize',14)
xlabel('$t$ $[sec]$','Interpreter','latex','FontSize',14)
box on
set(gcf,'color','w')
grid on

return

%% Window heating with RF
figure
D0max = max(max(D0));
% rngD0Mask = find(D0>14650);

for s = rng
     surf((TempMeasured(:,:,rng(1)).*WinArea)-(D0.*WinArea),'LineStyle','none')
      hold on
     surf((TempMeasured(:,:,s).*WinArea)-(D0.*WinArea),'LineStyle','none')

      view([0,90])
      caxis([0,15])
      colormap('hot')
      title(['shot: ',num2str(shot(1)),', frame: ',num2str(s)])
      if s == n1
                title(['RF start', 'frame: ',num2str(s)])
      elseif s == n2
                title(['RF end', 'frame: ',num2str(s)])
      end
      set(gca,'PlotBoxAspectRatio',[1 0.5 1])
      colormap(flipud(hot))
      colorbar
      drawnow
%      pause(0.01)
     hold off
end

dx = 16;
dy = 16;
y = 257;
x = 119;

rngx = [x - dx/2 : 1 : x + dx/2];
rngy = [y - dy/2 : 1 : y + dy/2];

hold on ; line(rngy,rngx,10*ones(size(rngx)))

%%
for s = 1:Nz
    TempHV(s) = mean(mean(TempMeasured(rngx,rngy,s)));
end

figure
plot(time,TempHV,'k.-')

Tmin = min(TempHV);
t0 = time(n1);
figure
T2 = (TempHV-Tmin).^2;
plot(time-t0,T2,'k.-')
xlim([0,1])

dT2dt = diff(T2)./diff(time-t0);
rho = 3300; 
kt   = 180 ;
cp  = 740 ;
a = kt/(rho*cp);

L = 6/1000; % ALN wall thickness in ProtoMPEX
ts = L*L/a;

figure; 
plot(time(1:end-1)-t0,sqrt((pi/4)*kt*rho*cp*dT2dt)*1e-3,'k.-')
xlim([-0.5,1])
ylabel('kWm^{-2}')
