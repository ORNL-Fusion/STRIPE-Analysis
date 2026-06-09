% Postprocessor 1
% This scripts performs the following tasks:
% 1 - load data from pre-processor 1
% 2 - Calculate the surface temperature given emissivity and optical
% transmission
% 3 - Animation of surface temperature
% 4 - Selection of temperature probe
% 5 - Plot temperature probe data
% 6 - Heat flux reconstruction with finite slab geometry
%       + Define material properties
%       + Define slab thickness
%       + Normalization and non-dimensionalization of experimental data
%       + Define convolution operator and kernel
%       + Apply the conjugate gradient method
% 7 - Compare experimental temperature and reconstruction
% 8 - Plot time evolution of surface heat flux
% 9 - Plot peak heat flux vs RF power
% 10- Compare finite and infinite slab model calculations 

clear all 
close all

% =========================================================================
% Select dataset
% =========================================================================
SaveFig = 0;


% =========================================================================
% Load data
% =========================================================================
load(['preprocess_Dataset_3.mat'])

% =========================================================================
% Load seq files
% =========================================================================
% seq files cannot saved correctly into the .mat format hence we need to
% reload them in order to the temperature data
a{1} = ['C:\Users\nfc\Documents\ProtoMPEX_Ops\2019\2019_11\2019_11_25\IR_Data'];

for s = 1:length(shot)
% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file'
PATHNAME = [a{addressLoc(s)},'\'];
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
videoFileName=[PATHNAME FILENAME];
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);

% Define the seq file and associated thermal parameters
seq{s} = file.ThermalSequencePlayer();
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 0.7;
seq{s}.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.Distance = 1;
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.ReferenceTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.Transmission = 1;
seq{s}.ThermalImage.ThermalParameters.RelativeHumidity = 0;
seq{s}.ThermalImage.ThermalParameters.ReflectedTemperature = 24;
end

% =========================================================================
%  Calculate surface temperature
% =========================================================================
tic
emissivity = 0.97;
for s = 1:length(shot)
    [temperature{s}] = IntensityTempConv(emissivity,intensity{s},seq{s});
    t_temperature{s} = t_intensity{s};
    temperature0{s} = temperature{s}(:,:,1);
end
toc

%%
% =========================================================================
% Animate data
% =========================================================================
try
    hf = findobj('Tag','temperatureVideo');
    close(hf)
end

xCenter{1} = 113;
yCenter{1} = 376;
xCenter{2} = 113;
yCenter{2} = 376;

for s = 1:length(shot)
    rngy_temperature_probe{s} = [(yCenter{s}-10):(yCenter{s}+10)];
    rngx_temperature_probe{s} = [(xCenter{s}-10):(xCenter{s}+10)];

    [X,Y] = meshgrid(rngy_temperature_probe{s},rngx_temperature_probe{s});
    Z = 16*ones(size(Y));
    for fr = 1:1:size(temperature{s},3)
        temperature_windowArea_noOffset = (temperature{s}(:,:,fr) - temperature0{s}).*windowArea{s};
        temperature_probe_noOffset{s}(fr) = mean(mean(temperature_windowArea_noOffset(rngx_temperature_probe{s},rngy_temperature_probe{s})));
        
        temperature_windowArea = temperature{s}(:,:,fr).*windowArea{s};
        temperature_probe{s}(fr) = mean(mean(temperature_windowArea(rngx_temperature_probe{s},rngy_temperature_probe{s})));
        
        intensity_windowArea = intensity{s}(:,:,fr).*windowArea{s};
        intensity_probe{s}(fr) = mean(mean(intensity_windowArea(rngx_temperature_probe{s},rngy_temperature_probe{s})));
    end
end

s = 1;
figure;
set(gcf,'color','w')
set(gcf,'Tag','temperatureVideo')
for fr = 1:3:size(temperature{s},3)
    s = 1;
    subplot(2,1,s)
    temperature_windowArea_noOffset = (temperature{s}(:,:,fr) - temperature0{s}).*windowArea{s};
    surf(temperature_windowArea_noOffset,'LineStyle','none')
    view([0,90])
    zlim([0,20])
    caxis([0,15])
    colormap(flipud(hot))
    h_title = title(['$${\Delta}T$$ [C]',' , frame: ',num2str(fr)]);
    set(h_title,'interpreter','latex','FontSize',12)
    axis('equal')
    xlim([0,size(temperature{s},2)])
    ylim([0,size(temperature{s},1)])
    hold on
    plot3(X(:,end),Y(:,end),Z(1,:),'k.','MarkerSize',0.5)
    plot3(X(1,:  ),Y(1,:  ),Z(1,:),'k.','MarkerSize',0.5)
    plot3(X(:,1  ),Y(:,1  ),Z(1,:),'k.','MarkerSize',0.5)
    plot3(X(end,:  ),Y(end,:  ),Z(1,:),'k.','MarkerSize',0.5)
    hold off
    h_cb = colorbar;
    h_text = text(50,150,{viewType{s};['shot: ',num2str(shot(s))];num2str(scanVar(s))});
    set(h_text,'interpreter','Latex','FontSize',11)
    set(gca,'YTick',[],'XTick',[])
    box on
    
    s = 2;
    subplot(2,1,s)
    temperature_windowArea_noOffset = (temperature{s}(:,:,fr) - temperature0{s}).*windowArea{s};
    surf(temperature_windowArea_noOffset,'LineStyle','none')
    view([0,90])
    zlim([0,20])
    caxis([0,15])
    colormap(flipud(hot))
    h_title = title(['$${\Delta}T$$ [C]',' , frame: ',num2str(fr)]);
    set(h_title,'interpreter','latex','FontSize',12)
    axis('equal')
    xlim([0,size(temperature{s},2)])
    ylim([0,size(temperature{s},1)])
    hold on
    plot3(X(:,end),Y(:,end),Z(1,:),'k.','MarkerSize',0.5)
    plot3(X(1,:  ),Y(1,:  ),Z(1,:),'k.','MarkerSize',0.5)
    plot3(X(:,1  ),Y(:,1  ),Z(1,:),'k.','MarkerSize',0.5)
    plot3(X(end,:  ),Y(end,:  ),Z(1,:),'k.','MarkerSize',0.5)
    hold off
    h_cb = colorbar;
    h_text = text(50,150,{viewType{s};['shot: ',num2str(shot(s))];num2str(scanVar(s))});
    set(h_text,'interpreter','Latex','FontSize',11)
    set(gca,'YTick',[],'XTick',[])
    box on
    drawnow
end

if SaveFig
    saveas(gcf,['temp surface ',viewType{s}],'tiffn')
end

%%
% =========================================================================
% Time evolution of probe data
% =========================================================================

figure
set(gcf,'Tag','temperatureVideo')
hold on
[~,b] = sort(scanVar);
for s = 1:length(shot)
    h_probe(s) = plot(t_temperature{b(s)},temperature_probe_noOffset{b(s)},'LineWidth',2);
end
% title(viewType{:})
box on
grid on
set(gcf,'color','w')
ylim([0,16])
xlim([0,0.8])
% l_intensity_probe = legend([h_probe],num2str(scanVar(b)'));
l_intensity_probe = legend([h_probe],viewType{scanVar(b)'});
l_intensity_probe.Location = 'northwest';
ylabel('$${\Delta}T$$ [C]','interpreter','Latex','FontSize',12)
xlabel('time [s]','interpreter','Latex','FontSize',12)
set(gca,'FontName','times','FontSize',11)

if SaveFig
    saveas(gcf,['delta T'],'tiffn')
end

%%
% =========================================================================
% Reconstruct heat flux from temperature probe
% =========================================================================

% ALN material properties
rho = 3300; 
kt   = 180 ;
cp  = 740 ;
a = kt/(rho*cp); % thermal diffusivity

% Geometry
Lz = 6/1000; % ALN wall thickness in ProtoMPEX

% Dimensionless and normalized data
% Order of magnitude heat flux:
q0 = 1e6; 
% Characteristic time scale:
t_star = Lz*Lz/a;
% Characteristic temperature:
T_star = q0*Lz/kt;
% Normalized measured temperature and time:
for s = 1:length(shot)
    z{s}    = temperature_probe_noOffset{s}'/T_star; 
    t_z{s}  = t_temperature{s}'/t_star;
    dt_z(s) = t_z{s}(2) - t_z{s}(1);
end
% Select characteritic noise level in Kelvin:
dT = 0.25; % K
% Create dimensionless noise data:
dnoise = (dT/(q0*Lz/kt));

% =========================================================================
% Create the convolution operator
% =========================================================================
xx0 = 0.0005;
% To capture the details of the front surface we need to use 1e4
% partial sum terms and use xx0 = 0.01 to 0.05

% Partial sum
% "n" is the index
Ns = 1e4;

for s = 1:length(shot)
     for ii = 1:length(t_z{s})
           K0{s}(ii) = G_Impulse_1D(xx0,t_z{s}(ii),Ns);
     end
     P{s} = toeplitz(K0{s},zeros(size(K0{s})));
end

%
% =========================================================================
%  Conjugate gradient method:
% =========================================================================

% Input data:
% P: convolution operator
% q{1}: initial guess at heat flux, column vector
% z: experimental data, column vector
% Ni: number of iteration
% Output data:
% u: calculated temperature, column vector
% J: residual, column vector
% q: minimized heat flux, structure, column vector

tic
for s = 1:length(shot)
    qdummy{1} = zeros(size(t_z{s}));
    [a1,a2,a3] = IHCP_ConjugateGradient(P{s},qdummy,z{s},dt_z(s),1900,3,[]);
    q{s} = a1{end};
    u{s} = a2{end};
    J{s} = a3;
end
toc

% -------------------------------------------------------------------------
% Compare xp data with reconstruction
try
    hf = findobj('Tag','xp vs reconstruction');
    close(hf)
end
figure; 
set(gcf,'Tag','xp vs reconstruction')
hold on
switch 1
    case 1
        instances = [1,2];
    case 2
        instances = [1];
end
colorForLine  = {'k','k'};
markerForPlot = {'sq','sq'};
ii = 1;
for s = instances
    h_reconstruction(ii) = plot(t_z{s}*t_star,u{s}*T_star,['r',markerForPlot{ii}],'MarkerSize',5,...
        'LineWidth',2)
    h_xpdata(ii) = plot(t_z{s}*t_star,z{s}*T_star,colorForLine{ii},'LineWidth',3)
    ii = ii + 1;
end
legend([h_xpdata(1),h_reconstruction(1)],'Exp.','Reconstruction')
xlim([0,0.8])
ylim([0,16])
box on
grid on
set(gcf,'color','w')
ylabel('$${\Delta}T$$ [C]','interpreter','Latex','FontSize',12)
xlabel('time [s]','interpreter','Latex','FontSize',12)
set(gca,'FontName','times','FontSize',11)


% -------------------------------------------------------------------------
% Time dependent recontructed heat flux
try
    hf = findobj('Tag','heatflux_timeEvolution');
    close(hf)
end

tStart = 0.14; tEnd = 0.19;
tStart = 0.2; tEnd = 0.5;

for s = 1:length(shot)
    rng_ss = find(t_temperature{s}>=tStart & t_temperature{s}<=tEnd);
    q_ss(s) = mean(q0*q{s}(rng_ss));
    q_ss_std(s) =  std(q0*q{s}(rng_ss),1);
    dq_ss(s) = sqrt((0.125*q_ss(s))^2 + q_ss_std(s)^2);
end

figure; 
hold on
set(gcf,'Tag','heatflux')
[~,b] = sort(scanVar);
for s = 1:length(shot)
    h_q(s) = plot(t_z{b(s)}*t_star,q{b(s)}*q0*1e-6,'LineWidth',2);
end
box on
grid on
set(gcf,'color','w')
ylim([0,0.7])
xlim([0,0.7])
l_heatflux = legend([h_q],viewType{scanVar(b)'});
set(gca,'FontName','times','FontSize',12)

if SaveFig
    saveas(gcf,['heat flux '],'tiffn')
end

%% -------------------------------------------------------------------------
% Comparing linear and non-linear recontruction
try
    hf = findobj('Tag','methodComparison');
    close(hf)
end

for s = 1:length(shot)
    dt = t_temperature{s}(2) - t_temperature{s}(1);
    q_infiniteSlab{s} = real(sqrt(0.25*pi*kt*rho*cp*diff(temperature_probe_noOffset{s}.^2)/dt));
end

switch numel(shot)
    case 2
        nx_sp = 1;
        ny_sp = 2;
    case {3,4}
        nx_sp = 2;
        ny_sp = 2;
     case {5,6}
        nx_sp = 3;
        ny_sp = 2;
end
        
figure; 
set(gcf,'Tag','methodComparison','color','w')
[~,s] = min(scanVar);
subplot(nx_sp,ny_sp,1)
hold on
h_comp(1) = plot(t_z{s}*t_star,q{s}*q0*1e-3,'k','LineWidth',3);
h_comp(2) = plot(t_z{s}(1:end-1)*t_star,q_infiniteSlab{s}*1e-3,'r','LineWidth',3);
ylim([0,700])
xlim([0,0.7])
axis square
box on
grid on
set(gca,'FontName','times','FontSize',11)
xlabel('time [s]','interpreter','Latex','FontSize',12)
ylabel('[kWm$$^{-2}$$]','Interpreter','latex')
l_comp = legend(h_comp,'Finite slab','Infinite slab');
set(l_comp,'Interpreter','latex','FontSize',11)
title(['P$_{RF}$: ',num2str(scanVar(s)),' kW'],'Interpreter','latex','FontSize',11)

[~,s] = max(scanVar);
subplot(nx_sp,ny_sp,2)
hold on
h_comp(1) = plot(t_z{s}*t_star,q{s}*q0*1e-3,'k','LineWidth',3);
h_comp(2) = plot(t_z{s}(1:end-1)*t_star,q_infiniteSlab{s}*1e-3,'r','LineWidth',3);
ylim([0,700])
xlim([0,0.7])
axis square
box on
grid on
set(gca,'FontName','times','FontSize',11)
xlabel('time [s]','interpreter','Latex','FontSize',12)
ylabel('[kWm$$^{-2}$$]','Interpreter','latex')
title(['P$_{RF}$: ',num2str(scanVar(s)),' kW'],'Interpreter','latex','FontSize',11)

if SaveFig
    saveas(gcf,['Heat transfer model ',magConfig{xp}],'tiffn')
end 