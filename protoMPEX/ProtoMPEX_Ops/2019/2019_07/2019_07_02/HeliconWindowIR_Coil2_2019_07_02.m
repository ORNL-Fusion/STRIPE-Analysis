% Helicon window IR imaging as a function of RF power in Proto-MPEX
% Script created 2019_07_23

clear all
% close all

XP = 2; 

switch XP
    case 1
        shot       = 26000 + [765,771,774];
        Coil2      =         [0  ,520,520];
        AddressLoc =         [1  ,1  ,1  ];
        ViewType = 'Top';
     
    case 2
        shot       = 26000 + [766];
        Coil2      =         [0  ];
        AddressLoc =         [1  ];
        ViewType = 'Middle';
        
%         shot       = 26000 + [766,769,770];
%         Coil2      =         [0  ,200,520];
%         AddressLoc =         [1  ,1  ,1  ];
%         ViewType = 'Middle';
        
    case 3
        shot       = 26000 + [767];
        Coil2      =         [0  ];
        AddressLoc =         [1  ];
        ViewType = 'Bottom';   
        
end

%% 
% =========================================================================
% Load seq files and create RawData
% =========================================================================

RootAddress = cd;
a{1} = [RootAddress,'\IR_RawData'];

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
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 1;
seq{s}.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.Distance = 1;
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
% seq{s}.ThermalImage.ThermalParameters.ReferenceTemperature = 0;
seq{s}.ThermalImage.ThermalParameters.Transmission = 1;
seq{s}.ThermalImage.ThermalParameters.RelativeHumidity = 0;

% seq{s}.ThermalImage.ThermalParameters.

%Get the pixels
img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
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
       MeanRawData{s}(ii) = mean(mean(RawData{s}(:,:,ii))); 
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

% =========================================================================
% Calculating threshold gradient
% =========================================================================
'Gradient'
tic
for s = 1:length(shot)
    D0{s} = Data{s}(:,:,1);
    D1{s} = Data{s}(:,:,end);
    [fx,fy] = gradient(D0{s});
    SqNormD0 = fy.^2 + fx.^2;
    
    % Need to remove effect from edges from SqNormD0:
    SqNormD0(1,:) = 0;
    SqNormD0(end,:) = 0;
    SqNormD0(:,1) = 0;
    SqNormD0(:,end) = 0;
    
    [~,nx(s)] = max(max(SqNormD0,[],2));
    [~,ny(s)] = max(max(SqNormD0,[],1));
    ThresholdGradient(s) = D0{s}(nx(s),ny(s));
end
toc

%%
% =========================================================================
%  Calculate surface temperature
% =========================================================================
tic
emissivity = 0.81;
emissivity = 0.88;
for s = 1:length(shot)
    [TempMeasured{s}] = IntensityTempConv(emissivity,Data{s}-0*700,seq{s});
%     t_TempMeasured{s} = (-1*dt):dt:(length(rng{s})-(1+1))*dt;
    t_TempMeasured{s} = t_Data{s};
    T0{s} = TempMeasured{s}(:,:,1);
end
toc

%%
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
legend(hIR,num2str(Coil2'))

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
legend(hIR,num2str(Coil2'))

%% Image prior to RF

for s = 1:length(shot)
    WinArea{s} = double(D0{s}>ThresholdGradient(s));
    NotWinArea{s} = double(D0{s}<ThresholdGradient(s));
end

nRows = ceil(sqrt(length(shot)));
nColu = nRows;

figure;
set(gcf,'name','Raw Image prior to RF')
for s = 1:length(shot)
    subplot(nRows,nColu,s)
    surf(D0{s},'LineStyle','none')
    hold on
    plot3(ny(s),nx(s),D0{s}(nx(s),ny(s)),'ko');
    title(['shot: ',num2str(shot(s))])
    view([0,90])
end

figure;
set(gcf,'name','Masked Image prior to RF')
for s = 1:length(shot)
    subplot(nRows,nColu,s)
    surf(D0{s}.*WinArea{s},'LineStyle','none')
    view([0,90])
    colormap('hot')
    title(['shot: ',num2str(shot(s))])
    set(gca,'PlotBoxAspectRatio',[1 0.5 1])
    colormap(flipud('hot'))
%     colorbar
end

%% Window heating with RF

% Select shot
k = 1;

switch ViewType
    case 'Top'
        % AddressLoc = 1
        % Probe 1:
        y{1}(1) = 240; x{1}(1) = 108; dx{1}(1) = 20; dy{1}(1) = 20;
        % Probe 2:
        y{1}(2) = 375; x{1}(2) = 35 ; dx{1}(2) = 20; dy{1}(2) = 20;
        % AddressLoc = 2
        % Probe 1:
        y{2}(1) = 247; x{2}(1) = 109; dx{2}(1) = 20; dy{2}(1) = 20;
        % Probe 2:
        y{2}(2) = 355; x{2}(2) = 25 ; dx{2}(2) = 20; dy{2}(2) = 20;
        
    case 'Bottom'
        % AddressLoc = 1
        % Probe 1:
        y{1}(1) = 298; x{1}(1) = 140; dx{1}(1) = 20; dy{1}(1) = 20;
        % AddressLoc = 2
        % Probe 1:
        y{2}(1) = 318; x{2}(1) = 107; dx{2}(1) = 20; dy{2}(1) = 20;
   case 'Middle'
        % AddressLoc = 1
        % Probe 1:
        y{1}(1) = 353; x{1}(1) = 44; dx{1}(1) = 20; dy{1}(1) = 20;
        y{1}(2) = 348; x{1}(2) = 13; dx{1}(2) = 20; dy{1}(2) = 20;
        % AddressLoc = 2
        % Probe 1:
         y{2}(1) = 348; x{2}(1) = 13; dx{2}(1) = 20; dy{2}(1) = 20;
         y{2}(2) = 348; x{2}(2) = 13; dx{2}(2) = 20; dy{2}(2) = 20;

end

for j = 1:2;
    for r = 1:length(y{1})
        rngx{j}{r} = [x{j}(r) - dx{j}(r)/2 : 1 : x{j}(r) + dx{j}(r)/2];
        rngy{j}{r} = [y{j}(r) - dy{j}(r)/2 : 1 : y{j}(r) + dy{j}(r)/2];
    end
end

%%
k = 1;
% k = 9;
if 0
figure
% D0max = max(max(D0{k}));
set(gcf,'color','w')
for s = 1:size(TempMeasured{k},3)
    
    surf((TempMeasured{k}(:,:,s).*WinArea{k})-(T0{k}.*WinArea{k}),'LineStyle','none')
    hold on
    j = AddressLoc(k);
    for r = 1:length(y{1})
        line(rngy{j}{r}([1,end])  ,rngx{j}{r}([1,1])    ,10*ones(size(rngx{j}{r}([1,1]))))
        line(rngy{j}{r}([1,end])  ,rngx{j}{r}([end,end]),10*ones(size(rngx{j}{r}([1,1]))))
        line(rngy{j}{r}([1,1])    ,rngx{j}{r}([1,end])  ,10*ones(size(rngx{j}{r}([1,1]))))
        line(rngy{j}{r}([end,end]),rngx{j}{r}([1,end])  ,10*ones(size(rngx{j}{r}([1,1]))))
    end
    
%     if s == 50
%         break
%     end

      view([0,90])
      caxis([0,8])
      colormap('hot')
      title(['shot: ',num2str(shot(k)),', frame: ',num2str(s)])
      if s == n1(k)
                title(['RF start', 'frame: ',num2str(s)])
      elseif s == n2(k)
                title(['RF end', 'frame: ',num2str(s)])
      end
      set(gca,'PlotBoxAspectRatio',[1 0.5 1])
      colormap(flipud(hot))
      colorbar
      hold off
%       drawnow
     pause(0.00001)
end
end

%
for k = 1:length(shot)
    for s = 1:size(TempMeasured{k},3)
        for r = 1:length(y{1})
            j = AddressLoc(k);
            TempProbe{k}{r}(s) = mean(mean(TempMeasured{k}(rngx{j}{r},rngy{j}{r},s)));
            Tmin{k}{r} = min(TempProbe{k}{r});
        end
    end
    
end

%%
r = 2;
figure;
set(gcf,'Position',[  360.3333  197.6667  314.0000  420.0000])
subplot(3,1,1)
hold on
for s = 1:length(shot)
    hIR(s) = plot(t_TempMeasured{s}-t_TempMeasured{s}(1),TempProbe{s}{r}-Tmin{s}{r},'LineWidth',2);
    %     title(['shot: ',num2str(shot(1))])
    xlim([-0.1,1])
    ylabel('\DeltaT')
    ylim([0,10])
end
legend(hIR,num2str(Coil2'))
box on 
set(gcf,'color','w')

subplot(3,1,2);
hold on
for s = 1:length(shot)
    T2{s} = (TempProbe{s}{r}-Tmin{s}{r}).^2;
    hIR(s) = plot(t_TempMeasured{s}-t_TempMeasured{s}(1),T2{s},'LineWidth',2)
    xlim([-0.1,1])
    ylabel('(T-Tmin)^{2}')
%     title(['shot: ',num2str(shot(1))])
    ylim([0,100])
end
legend(hIR,num2str(Coil2'))
box on 
set(gcf,'color','w')

% ALN properties
rho = 3300; 
kt   = 180 ;
cp  = 740 ;
a = kt/(rho*cp);
L = 6/1000; % ALN wall thickness in ProtoMPEX
ts = L*L/a;

subplot(3,1,3)
hold on
for s = 1:length(shot)
    dT2dt{s} = diff(T2{s})./diff(t_TempMeasured{s}-t_TempMeasured{s}(1));
    HeatFlux{s} = sqrt((pi/4)*kt*rho*cp*dT2dt{s});
    hIR(s) = plot(t_TempMeasured{s}(1:end-1)-t_TempMeasured{s}(1),HeatFlux{s}*1e-3,'LineWidth',2);
    xlim([-0.1,1])
    ylabel('kWm^{-2}')
%     title(['shot: ',num2str(shot(1))])
    ylim([0,500])
end
legend(hIR,num2str(Coil2'))
box on 
xlabel('t [s]')
set(gcf,'color','w')

% plot the mean heat flux versus power
tStart{1} = 0.14; tEnd{1} = 0.19;
tStart{2} = 0.3; tEnd{2} = 0.4;
% tStart{2} = 0.14; tEnd{2} = 0.19;

for s = 1:length(shot)
    j = AddressLoc(s);
    rngMean = find((t_TempMeasured{s}-t_TempMeasured{s}(1))>=tStart{j} & (t_TempMeasured{s}-t_TempMeasured{s}(1))<=tEnd{j});
    MeanHeatFlux(s) = mean(HeatFlux{s}(rngMean));
    StdHeatFlux(s) =  std(HeatFlux{s}(rngMean),1);
end

figure; 
L(1) = errorbar(RFpNet,real(MeanHeatFlux)*1e-3,real(StdHeatFlux)*1e-3);
set(L(1),'Marker','sq','LineStyle','none','MarkerSize',5)
xlim([0,200])
ylim([0,500])
title(['Peak Surface Heatflux vs RF power, ',ViewType])
ylabel('kWm^{-2}')
xlabel('RF [kW]')
box on
grid on
set(gcf,'color','w')
set(gcf,'position',[360.3333  356.3333  388.6667  261.3333])