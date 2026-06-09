close all 
clear all

XP = 3; 
% On 2019_06_19 SAIREM + FRT-86 operations
switch XP
    case 1
        shot    = 26000 + [658,659,660,661,663,664];
        RFpNet  =         [104,122,135,145,160,169];
        ViewType = 'Bottom';
        n1_offset = 150;
    case 2
        shot    = 26000 + [670,668,665];
        RFpNet  =         [96 ,135,169];
        shot    = 26000 + [665,666,667,668,669,670,674];
        RFpNet  =         [169,161,154,135,114,96 ,144];
        ViewType = 'Middle';
        n1_offset = 1;
    case 3
        shot    = 26000 + [646,647,648,649,650,653,655,656];
        RFpNet  =         [56 ,89 ,138,153,150,134,165,104];
        ViewType = 'Top';
        n1_offset = 150;
end

%% 
% =========================================================================
% Load seq files and create RawData
% =========================================================================

% a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_06_19';
RootAddress = cd;
a{1} = [RootAddress,'\IR_RawData'];

for s = 1:length(shot)
PATHNAME = [a{1},'\'];
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
seq{s} = file.ThermalSequencePlayer();
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
emissivity = 0.8;
for s = 1:length(shot)
    [TempMeasured{s}] = IntensityTempConv(emissivity,Data{s},seq{s});
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
legend(hIR,num2str(RFpNet'))

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
legend(hIR,num2str(RFpNet'))


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

% Select Probe. This selects the region of space to sample
RegionType = 2;

switch ViewType
    case 'Top'
        % For top view (312 deg)
        % Probe 1:
        y(1) = 257; x(1) = 119; dx(1) = 20; dy(1) = 20;
        % Probe 2:
        y(2) = 375; x(2) = 35 ; dx(2) = 20; dy(2) = 20;

    case 'Bottom'
        % Probe 1:
        y(1) = 318; x(1) = 107; dx(1) = 20; dy(1) = 20;
   case 'Middle'
        % Probe 1:
        y(1) = 375; x(1) = 81; dx(1) = 20; dy(1) = 20;
end

for r = 1:length(y)
    rngx{r} = [x(r) - dx(r)/2 : 1 : x(r) + dx(r)/2];
    rngy{r} = [y(r) - dy(r)/2 : 1 : y(r) + dy(r)/2];
end

if 0
figure
D0max = max(max(D0{k}));
set(gcf,'color','w')
for s = 1:size(TempMeasured{k},3)
    
    surf((TempMeasured{k}(:,:,s).*WinArea{k})-(T0{k}.*WinArea{k}),'LineStyle','none')
    hold on
    for r = 1:length(y)
        line(rngy{r}([1,end]),rngx{r}([1,1]),10*ones(size(rngx{r}([1,1]))))
        line(rngy{r}([1,end]),rngx{r}([end,end]),10*ones(size(rngx{r}([1,1]))))
        line(rngy{r}([1,1]),rngx{r}([1,end]),10*ones(size(rngx{r}([1,1]))))
        line(rngy{r}([end,end]),rngx{r}([1,end]),10*ones(size(rngx{r}([1,1]))))
    end
    
%     if s == 50
%         break
%     end

      view([0,90])
      caxis([0,8])
      colormap('hot')
      title(['shot: ',num2str(shot(k)),', frame: ',num2str(s)])
      if s == n1
                title(['RF start', 'frame: ',num2str(s)])
      elseif s == n2
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

%%
for k = 1:length(shot)
    for s = 1:size(TempMeasured{k},3)
        for r = 1:length(y)
            TempProbe{k}{r}(s) = mean(mean(TempMeasured{k}(rngx{r},rngy{r},s)));
            Tmin{k}{r} = min(TempProbe{k}{r});
        end
    end
    
end

%%
r = 1;
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
legend(hIR,num2str(RFpNet'))
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
legend(hIR,num2str(RFpNet'))
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
legend(hIR,num2str(RFpNet'))
box on 
xlabel('t [s]')
set(gcf,'color','w')

% plot the mean heat flux versus power
tStart = 0.14;
tEnd = 0.19;
for s = 1:length(shot)
    rngMean = find((t_TempMeasured{s}-t_TempMeasured{s}(1))>=tStart & (t_TempMeasured{s}-t_TempMeasured{s}(1))<=tEnd);
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