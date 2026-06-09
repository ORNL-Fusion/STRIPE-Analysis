% Determining the emissivity of the helicon window
clear all
close all

% notes:
% The best shots for calibration identified so far are the following:

% - 26610: 06_18, first pulse of the day (RF trip). all FP at 24.4 deg C,
% associated IR data available, uniform IR emission 

% - 26646: 06_19, 1st plasma shot of the day 56 kW net, 100  ms pulse, all
% FP are identical at 23.8 deg C. Associated IR data available with uniform
% IR emission prior to the RF. RF starts on frame 289

% - 26674: 09_19, there was a 17 min gap between this shot and the previous plasma
% thus the temperature fileds on the windwo have spread, the temperature
% range of the helicon window is [28,34] deg C. the RF starts on frame 189
% to 190 on the IR data

shot       = 26000 + [610,646,674];
AddressLoc =         [1  ,2  ,2  ];

a{1} = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2019\2019_06\2019_06_18\IR_RawData';
a{2} = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2019\2019_06\2019_06_19\IR_RawData';

for s = 1:length(shot)
% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file'
PATHNAME = [a{AddressLoc(s)},'\'];
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
videoFileName=[PATHNAME FILENAME];
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);

seq{s} = file.ThermalSequencePlayer();
%Get the pixels
img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);

% Thermal Parameters:
seq{s}.ThermalImage.ThermalParameters.Emissivity = 1;
seq{s}.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.Transmission = 1;
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 0.7;
seq{s}.ThermalImage.ThermalParameters.RelativeHumidity = 0.1;
seq{s}.ThermalImage.ThermalParameters.ReferenceTemperature = 24;

RawData{s}(:,:,1) = im;
fr = 1;
if(seq{s}.Count > 1)
    while(seq{s}.Next())
        img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
        RawData{s}(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;
    end
end
end

%%
% =========================================================================
% Determine the start of the RF pulse:
% =========================================================================
'Data subset'
n1_offset = 10;
n2_offset = 40;
n_Before = 3;
n_After = 60;
FrameRate = 100;
dt = 1/FrameRate;

tic
for s = 1:length(shot)
    [Nx(s),Ny(s),Nz(s)] = size(RawData{s}); 

    for ii = 1:Nz(s)
       MeanRawData{s}(ii) = mean(RawData{s}(125,:,ii)); 
    end
    t_MeanRawData{s} = 0:dt:(Nz(s)-1)*dt;
    
    % Find the start and end of the RF pulse:
    [~,n1(s)] = max(diff(MeanRawData{s}(n1_offset:end),1));
    n1(s) = n1(s) + n1_offset - 1;
    [~,n2(s)] = min(diff(MeanRawData{s}(n2_offset:end),1));
    n2(s) = n2(s) + n2_offset - 1;
    
    t0_MeanRawData(s) = t_MeanRawData{s}(n1(s));
    
    % Define the time window over which the RF is on:
    % Include some points before and after the RF
    rng{s} = [n1(s)-n_Before:n2(s)+n_After];
    % Extract the relevent data from all the raw data:
    Data{s} = RawData{s}(:,:,rng{s});
    t_Data{s} = 0:dt:(length(rng{s})-1)*dt;
end
toc

%% Plot the data
figure; 
subplot(2,1,1)
hold on
for s = 1:length(shot)
    MeanRaw{s} = MeanRawData{s}-min(MeanRawData{s});
    plot(t_MeanRawData{s}-t0_MeanRawData(s),MeanRaw{s}.^1,'LineWidth',0.5);
    hIR(s) = plot(t_MeanRawData{s}(rng{s})-t0_MeanRawData(s),MeanRaw{s}(rng{s}).^1,'LineWidth',2);
end
title('Mean RawData')
ylabel('{\Delta}Intensity')
xlabel('t [s]')
xlim([-0.5,2])
legend(hIR,num2str(shot'))

subplot(2,1,2)
hold on
for s = 1:length(shot)
    plot(t_MeanRawData{s}-t0_MeanRawData(s),MeanRaw{s}.^2,'LineWidth',0.5);
    hIR(s) = plot(t_MeanRawData{s}(rng{s})-t0_MeanRawData(s),MeanRaw{s}(rng{s}).^2,'LineWidth',2);
end
title('Mean RawData^2')
ylabel('{\Delta}Intensity^2')
xlabel('t [s]')
xlim([-0.5,2])
legend(hIR,num2str(shot'))

%%
% =========================================================================
%  Calculate surface temperature
% =========================================================================
tic
emissivity = 1;
IntensityOffset = 0; 
for s = 1:length(shot)
    [TempMeasured{s}] = IntensityTempConv(emissivity,Data{s}-IntensityOffset,seq{s});
    t_TempMeasured{s} = t_Data{s};
    T0{s} = TempMeasured{s}(:,:,1);
end
toc

%% Plot IR inferred temperatures

figure;
hold on
for k = 1:length(shot)
    surf(T0{k},'LineStyle','none')
    view([0,0])
    zlim([20,35])
    title(['shot: ',num2str(shot(k))])
end

%%
DataFromServer = 0;
if DataFromServer
    % To convert voltage to temperature, multiply by 20
    Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

    % Gather data
    DA{1} = [RootAddress,'FLUOROPT_1']; 
    DA{2} = [RootAddress,'FLUOROPT_2']; 
    DA{3} = [RootAddress,'FLUOROPT_3']; 
    DA{4} = [RootAddress,'FLUOROPT_4']; 

    [f1,t_f1]   = my_mdsvalue_v2(shot,DA(1));
    [f2,t_f2]   = my_mdsvalue_v2(shot,DA(2));
    [f3,t_f3]   = my_mdsvalue_v2(shot,DA(3));
    [f4,t_f4]   = my_mdsvalue_v2(shot,DA(4));

    save('EmissivityCalibration_FP_Data.mat','f1','f2','f3','f4','t_f1','t_f2','t_f3','t_f4');
else
    load('EmissivityCalibration_FP_Data.mat')
end

C1 = {'b','b.','b:','b--'};
C2 = {'r','r.','r:','r--'};
C3 = {'g','g.','g:','g--'};
C4 = {'c','c.','c:','c--'};

k = 1;
b = 0;
figure;
for k = 1:length(shot)
    subplot(2,2,k);
    hold on
    plot(t_f1{k}(1:end-1),(f1{k}-b*min(f1{k}))*20,C1{k},'LineWidth',2)
    plot(t_f2{k}(1:end-1),(f2{k}-b*min(f2{k}))*20,C2{k},'LineWidth',2) % Ground Side
    plot(t_f3{k}(1:end-1),(f3{k}-b*min(f3{k}))*20,C3{k},'LineWidth',2) % High Voltage side
    plot(t_f4{k}(1:end-1),(f4{k}-b*min(f4{k}))*20,C4{k},'LineWidth',2)

    ylim([0,100])
    ylabel('${\Delta}T$ $[C]$','Interpreter','latex','FontSize',14)
    xlabel('$t$ $[sec]$','Interpreter','latex','FontSize',14)
    box on
    set(gcf,'color','w')
    grid on
    title(['shot: ',num2str(shot(k))])
end




