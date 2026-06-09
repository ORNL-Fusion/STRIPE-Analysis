close all 
clear all

shot    = 26000 + [634 ,635 ,639 ];
RFpNet  =         [88  ,102 ,138 ];

ViewType = 'Top';
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_06_18';
a{2} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_06_18';
a{3} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_06_18';

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

% =========================================================================
%  Calculate surface temperature
% =========================================================================
emissivity = 0.8;
for s = 1:length(shot)
    [TempMeasured{s},dTMeasured{s}] = IntensityTempConv2(emissivity,0.2*emissivity,RawData{s},seq{s});
    % framerate = 100;
    framerate = seq{s}.FrameRate;
    dt = 1/framerate;
    time{s} = 0:dt:(fr-2)*dt;
end

% =========================================================================
% Determine the start of the RF pulse:
% =========================================================================
for s = 1:length(shot)
    Nx(s) = size(TempMeasured{s}(:,:,:),1);
    Ny(s) = size(TempMeasured{s}(:,:,:),2);
    Nz(s) = size(TempMeasured{s}(:,:,:),3);

    PixelCount(s) = Nx(s)*Ny(s);

    for ii = 1:Nz(s)
       MeanTempMeasured{s}(ii) = sum(sum(TempMeasured{s}(:,:,ii)))/PixelCount(s); 
    end

    [~,n1(s)] = max(diff(MeanTempMeasured{s}));
    [~,n2(s)] = min(diff(MeanTempMeasured{s}));
    t0(s) = time{s}(n1(s));
    
    rng{s} = [n1(s)-10:n2(s)+60];
    D0{s} = TempMeasured{s}(:,:,n1(s)-10);
    D1{s} = TempMeasured{s}(:,:,n2(s)+100);
    
    % Calculating threshold Temp
    [fx,fy] = gradient(D0{s});
    [~,nx] = max(max(fy.^2 + fx.^2,[],2));
    [~,ny] = max(max(fy.^2 + fx.^2,[],1));
    TempThresholdCalc(s) = D0{s}(nx,ny);
end

%%
figure; 
subplot(2,1,1)
hold on
for s = 1:length(shot)
    MeanTemp{s} = MeanTempMeasured{s}-min(MeanTempMeasured{s});
    plot(time{s}-t0(s),MeanTemp{s}.^1,'LineWidth',0.5);
    hIR(s) = plot(time{s}(rng{s})-t0(s),MeanTemp{s}(rng{s}).^1,'LineWidth',2);
end
title('Mean Front Surface T')
ylabel('\Delta T [C]')
xlabel('t [s]')
xlim([-0.5,2])
legend(hIR,num2str(RFpNet'))

subplot(2,1,2)
hold on
for s = 1:length(shot)
    MeanTemp{s} = MeanTempMeasured{s}-min(MeanTempMeasured{s});
    plot(time{s}-t0(s),MeanTemp{s}.^2,'LineWidth',0.5);
    hIR(s) = plot(time{s}(rng{s})-t0(s),MeanTemp{s}(rng{s}).^2,'LineWidth',2);
end
title('Mean Front Surface T^2')
ylabel('\Delta T^2')
xlabel('t [s]')
xlim([-0.5,2])
legend(hIR,num2str(RFpNet'))


%% Image prior to RF
% 
% WinArea = double(D0>TempThreshold);
% NotWinArea = double(D0<TempThreshold);

for s = 1:length(shot)
    WinArea{s} = double(D0{s}>TempThresholdCalc(s));
    NotWinArea{s} = double(D0{s}<TempThresholdCalc(s));
end

nRows = ceil(sqrt(length(shot)));
nColu = nRows;

figure;
for s = 1:length(shot)
    subplot(nRows,nColu,s)
    surf(D0{s},'LineStyle','none')
end

figure;
for s = 1:length(shot)
    subplot(nRows,nColu,s)
    surf(D0{s}.*WinArea{s},'LineStyle','none')
    view([0,90])
    caxis([3,50])
    colormap('hot')
    title(['shot: ',num2str(shot(s))])
    set(gca,'PlotBoxAspectRatio',[1 0.5 1])
    colormap(flipud('hot'))
    colorbar
    % zlim([31,33])
    % caxis([31,33])
end

%% Window heating with RF

% Select shot
k = 3;

% Select Probe. This selects the region of space to sample
RegionType = 2;

% For top view (312 deg)
switch RegionType
    case 1
        y = 257;
        x = 119;
        dx = 20;
        dy = 20;
    case 2
        y = 375;
        x = 35;
        dx = 20;
        dy = 20;
end

% Probe 1:
y(1) = 257; x(1) = 119; dx(1) = 20; dy(1) = 20;
% Probe 2:
y(2) = 375; x(2) = 35 ; dx(2) = 20; dy(2) = 20;

for r = 1:length(y)
    rngx{r} = [x(r) - dx(r)/2 : 1 : x(r) + dx(r)/2];
    rngy{r} = [y(r) - dy(r)/2 : 1 : y(r) + dy(r)/2];
end


figure
D0max = max(max(D0{k}));
set(gcf,'color','w')
for s = rng{k}

     surf((TempMeasured{k}(:,:,s).*WinArea{k})-(D0{k}.*WinArea{k}),'LineStyle','none')
    hold on
    for r = 1:length(y)
        line(rngy{r}([1,end]),rngx{r}([1,1]),10*ones(size(rngx{r}([1,1]))))
        line(rngy{r}([1,end]),rngx{r}([end,end]),10*ones(size(rngx{r}([1,1]))))
        line(rngy{r}([1,1]),rngx{r}([1,end]),10*ones(size(rngx{r}([1,1]))))
        line(rngy{r}([end,end]),rngx{r}([1,end]),10*ones(size(rngx{r}([1,1]))))
    end
    
    if s == 50
        break
    end

      view([0,90])
      caxis([0,15])
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
     pause(0.01)
end


%%
for k = 1:length(shot)
    for s = 1:Nz
        for r = 1:length(y)
            TempProbe{k}{r}(s) = mean(mean(TempMeasured{k}(rngx{r},rngy{r},s)));
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
    hIR(s) = plot(time{s}-t0(s),TempProbe{s}{r}-Tmin{s}{r},'LineWidth',2);
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
    hIR(s) = plot(time{s}-t0(s),T2{s},'LineWidth',2)
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
    dT2dt{s} = diff(T2{s})./diff(time{s}-t0(s));
    HeatFlux{s} = sqrt((pi/4)*kt*rho*cp*dT2dt{s});
    hIR(s) = plot(time{s}(1:end-1)-t0(s),HeatFlux{s}*1e-3,'LineWidth',2);
    xlim([-0.1,1])
    ylabel('kWm^{-2}')
%     title(['shot: ',num2str(shot(1))])
    ylim([0,400])
end
legend(hIR,num2str(RFpNet'))
box on 
xlabel('t [s]')
set(gcf,'color','w')

% plot the mean heat flux versus power
for s = 1:length(shot)
    rngMean = find((time{s}-t0(s))>=0.38 & (time{s}-t0(s))<=0.48);
    MeanHeatFlux(s) = mean(HeatFlux{s}(rngMean));
end

figure; 
plot(RFpNet,real(MeanHeatFlux)*1e-3,'ko-')
xlim([50,150])
ylim([50,300])
title('Peak Surface Heatflux vs RF power')
ylabel('kWm^{-2}')
xlabel('RF [kW]')
box on
grid on
set(gcf,'color','w')