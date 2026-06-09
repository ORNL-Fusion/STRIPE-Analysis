% Understanding how the emissivity controls the relationship between pixel
% intensity and temperature

clear all
close all

% notes:
% The best shots for calibration identified so far are the following:

% - 26610: 06_18, first pulse of the day (RF trip). all FP at 24.4 deg C,
% associated IR data available, uniform IR emission 

% - 26664: 06_19, bottom view, 169 kW net.

% Process:
% We use two reference shots, the first one corresponds to conditions at
% the start of the day when no plasma has been produced. the window
% temperature according to the FP is uniform at about 24 deg C. the second
% shot corresponds to a case with the highest net RF power produced in
% Proto-MPEX at 169 kW.
% First, we identify the frame with the largest surface temperature which
% at the enf of the RF pulse. We exclude the intensities at the edges of
% the window where reflections from the stainless steel affect the
% measurement.
% Second, knowing the range of intensities, we calculate a surface for the
% estimated temperature vs emissivity and intensity.
% The aim is to fit linear equations to the surface to approximate the
% partial derivatives.

shot       = 26000 + [610,664];
AddressLoc =         [1  ,2  ];

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
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 0.7;
seq{s}.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.Distance = 1;
seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.ReferenceTemperature = 24;
seq{s}.ThermalImage.ThermalParameters.Transmission = 1;
seq{s}.ThermalImage.ThermalParameters.RelativeHumidity = 0;
seq{s}.ThermalImage.ThermalParameters.ReflectedTemperature = 24;

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
end
toc


if 1
    figure; 
    subplot(2,1,1)
    hold on
    for s = 1:length(shot)
        MeanRaw{s} = MeanRawData{s}-min(MeanRawData{s});
        hIR(s) = plot(t_MeanRawData{s}-t0_MeanRawData(s),MeanRaw{s}.^1,'LineWidth',2);
    end
    title('Mean RawData')
    ylabel('{\Delta}Intensity')
    xlabel('t [s]')
    xlim([-0.5,2])
    legend(hIR,num2str(shot'))
    box on

    subplot(2,1,2)
    hold on
    for s = 1:length(shot)
        hIR(s) = plot(t_MeanRawData{s}-t0_MeanRawData(s),MeanRaw{s}.^2,'LineWidth',2);
    end
    title('Mean RawData^2')
    ylabel('{\Delta}Intensity^2')
    xlabel('t [s]')
    xlim([-0.5,2])
    legend(hIR,num2str(shot'))
    box on
end

%% Determine min and max intensities
% Min temperature
k = 1;
MySurf1(RawData{k}(:,:,1));

% The mininum intensities are obtained from shot 26610 which has no RF and
% is the first shot of the day and is given by: 14394

% max temperature
k = 2;
MySurf1(RawData{k}(:,:,n2(k)))

% The maximum intensity is obtained from shot 26664. To do this we need to
% neglect the large intensities that occur at the edges of the window which
% are likely due to reflection from the stainless steel flanges. the
% largest intensities observed in the window are can be found by
% identifying the end of the RF pulse and are in the vecinity of 1.7e4

% Thus the range of intensities that we are concerned with is [1.4,1.7]*1e4

%%
% close all
MinInten = 14e3;
MaxInten = 17e3;
RefInten = 14500;
N_Inten = 1e2;

MaxEmiss = 1;
MinEmiss = 0.8;
N_Emiss = 0.8e2;

Inten = linspace(MinInten,MaxInten,N_Inten);
Emiss = linspace(MinEmiss,MaxEmiss,N_Emiss);

% Inten2 = linspace(0,64e3);
Inten2 = logspace(2,log10(64000),100);

% Extract the measured temperature given the emissivity and intensity
for ii=1:N_Inten 
    for jj = 1:length(Emiss)
        T(jj,ii) = seq{1}.ThermalImage.GetValueFromEmissivity(Emiss(jj),Inten(ii));
    end
end

for ii = 1:length(Inten2)
    T2(ii) = seq{1}.ThermalImage.GetValueFromEmissivity(1,Inten2(ii));
end

% =========================================================================
% Plot the temperature vs Intensity over the whole range
% =========================================================================
ScaleType = 'log';
figure;
hold on
[a,n] = max(diff(T2));

h(1) = plot(Inten2-Inten2(n),T2-T2(n),'r.-','LineWidth',2);
h(2) = plot((0.71e-6)*(T2-T2(n)).^4,T2-T2(n),'g.-','LineWidth',2);

box on
grid on
set(gcf,'color','w')
set(gca,'YScale',ScaleType);
set(gca,'XScale',ScaleType);

% =========================================================================
% Plot the temperature vs Intensity for different emissivities
% =========================================================================
figure; hold on
h(1) = plot(Inten,T(1,:),'k','LineWidth',2);
h(2) = plot(Inten,T(round(N_Emiss/2),:),'r','LineWidth',2);
h(3) = plot(Inten,T(N_Emiss,:),'g','LineWidth',2);
box on
set(gcf,'color','w')
legend(h,'\epsilon = 0.8','\epsilon = 0.9','\epsilon = 1')
xlabel('Intensity')
ylabel('T [C]')

% Fit linear functions to these plots:
p{1} = polyfit(Inten,T(1,:),1); 
p{2} = polyfit(Inten,T(round(N_Emiss/2),:),1);
p{3} = polyfit(Inten,T(N_Emiss,:),1);

% Emissivity values:
e = [Emiss(1),Emiss(round(N_Emiss/2)),Emiss(N_Emiss)];
% Intercept:
a = [p{1}(2) ,p{2}(2)                ,p{3}(2)       ];
% Slope:
b = [p{1}(1) ,p{2}(1)                ,p{3}(1)       ];

c1c2 = polyfit(e,a,1);
c1 = c1c2(2)
c2 = c1c2(1)

c3c4 = polyfit(e,b,1);
c3 = c3c4(2)
c4 = c3c4(1)


% =========================================================================
% Plot the temperature vs Emissivity for different Intensities
% =========================================================================
figure; hold on
h(1) = plot(Emiss,T(:,1),'k','LineWidth',2);
h(2) = plot(Emiss,T(:,round(N_Inten/2)),'r','LineWidth',2);
h(3) = plot(Emiss,T(:,N_Inten),'g','LineWidth',2);
box on
set(gcf,'color','w')
legend(h,'I = 14e3','I = 15.4e3','I = 17e3')
xlabel('Emissivity')
ylabel('T [C]')

% Fit linear functions to these plots:
q{1} = polyfit(Emiss,T(:,1)',1);
q{2} = polyfit(Emiss,T(:,round(N_Inten/2))',1);
q{3} = polyfit(Emiss,T(:,N_Inten)',1);

% For each emissivity value, we can approximate the value of T as a linear
% function of Intensity where the intercept and slope are functions of
% emissivity. The slope appears to be linear with emissivity.
% so we can approximate the data T as follows:

% T(e,i) = a(e) + b(e)*i

% where: a(e) = c1 + c2*e ; b(e) = c3 + c4*e

% Thus we get the following:

% T(e,i) = c1 + c2*e + c3*i + c4*e*i

% The above expression is very similar to a plane if one neglects the term
% c4*e*i. Provided this term is small we can use a plane.

% Method A: Fit a plane using 3 points
% =============================================================
n_pe = 1;
n_pi = 1;
P = [Emiss(n_pe),Inten(n_pi),T(n_pe,n_pi)];
n_qe = N_Emiss;
n_qi = N_Inten;
Q = [Emiss(n_qe),Inten(n_qi),T(n_qe,n_qi)];
n_re = N_Emiss;
n_ri = 1;
R = [Emiss(n_re),Inten(n_ri),T(n_re,n_ri)];

normal = cross(R,P) + cross(P,Q) + cross(Q,R);

Ta = @(E,I) ((P*cross(Q,R)'/normal(3)) -(normal(1)/normal(3)).*E - (normal(2)/normal(3)).*I)';
[E,I] = meshgrid(Emiss,Inten);


% Method B: Fit a surface with nonlinear term e*i
% =============================================================
C0(1) = c1;
C0(2) = c2;
C0(3) = c3;
C0(4) = c4;

Tb = @(C) (C(1) + C(2)*E + C(3)*I + C(4)*E.*I)';
dTbdi = @(C) C(3) + C(4)*E;
res = @(C) sum(sum(Tb(C) - T));
[C1,~,~] = LMFnlsq(res,C0,'MaxIter',30);

% =========================================================================
% Calculate and plot both methods:
% =========================================================================
figure
hold on
surf(Emiss,Inten,T','LineStyle','none')
xlabel('Emissivity')
ylabel('Inten')
title('T interpolated')
plot3(P(1),P(2),P(3),'ko')
plot3(Q(1),Q(2),Q(3),'ro')
plot3(R(1),R(2),R(3),'go')
% plot3(E,I,Ta(E,I)','m.')
surf(Emiss,Inten,Tb(C1)')
view([45,45])

% =========================================================================
% Calculate and plot the error for both methods:
% =========================================================================
res_a = Ta(E,I) - T;
res_b = Tb(C1) - T;

figure; 
subplot(1,2,1); hold on
surf(E,I,100*abs(res_a')./T','LineStyle','none')
colorbar
title('% Error in T_a')
xlabel('Emissivity')
ylabel('Intensity')
box on
axis(gca,'square')
caxis([0,3])

subplot(1,2,2); hold on
surf(E,I,100*abs(res_b')./T','LineStyle','none')
colorbar
title('% Error in T_b')
xlabel('Emissivity')
ylabel('Intensity')
box on
set(gcf,'color','w')
axis(gca,'square')
caxis([0,3])

% The calculation shows that method b has a much lower error. the largest
% error is about 2.5 %

% =========================================================================
% Calculate the gradient of T and Tb
% =========================================================================
de = Emiss(2)-Emiss(1);
di = Inten(2)-Inten(1);
[Te ,Ti ] = gradient(T'     ,de,di);
[Tbe,Tbi] = gradient(Tb(C1)',de,di);

figure; hold on
surf(Emiss,Inten,Ti ,'LineStyle','none')
surf(Emiss,Inten,Tbi,'LineStyle','none')

xlabel('Emissivity')
ylabel('Inten')
title('$\partial T/\partial i$','Interpreter','latex')
zlim([0,0.01])
view([45,45])

% =========================================================================
% Calculate and plot the error in the gradient of Tb
% =========================================================================
res_dTbdi = Tbi - Ti;
figure
surf(E,I,100*abs(res_dTbdi./Ti),'LineStyle','none')
xlabel('Emissivity')
ylabel('Intensity')
colorbar
caxis([0,15])
view([0,90])


figure
dE = 0.2;
plot(Emiss,abs(dE*C1(4)./(C1(3) + C1(4)*Emiss)))
ylim([0,0.4])
grid on
xlabel('Emissivity')
ylabel('Error due to uncertainty in \epsilon')


figure
dE = 0.2;
plot(Emiss,C1(3) + C1(4)*Emiss)
ylim([0,1e-2])
grid on
xlabel('Emissivity')
ylabel('b(\epsilon)')

% Within the small interval of Intensities that we are concerned with, the
% temperature depends linearly on Intensity hence, the change in
% temperature is linearly related to the change in Intensity and is not
% dependent on having absolute values of temperatures.
% Under these circumstances, the main uncertainty comes from the
% unceratinty on the slope of T vs Inten due to uncertainty in the value of
% emissivity.
% This slope is a linear function of the emissivity, hence for an interval
% of emissivities, we will end up with range of slopes and thus an interval
% of heat fluxes for a given change in intensity
