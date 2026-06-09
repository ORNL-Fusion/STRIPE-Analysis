% Step 0:
% =========================================================================
% Create diagrams of device:
% a: Whole device diagram
% b: Magnetic field profile:
% c: ECH region
% d: Target-IR camera region

clearvars
clc
close all

saveFig = 1;

%% Calculate magnetic field and fluxes:

% Magnetic configuration of interest:
% =========================================================================
confType = 'conf_G';

% Read "CoilSetup" spreadsheet:
% =========================================================================
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet',confType);
disp('Reading complete!')
toc(dum1);

% Display "coilSetup" on CLI:
% =========================================================================
coilSetup

% Define computational domain:
% =========================================================================
z_Dump = 0.5;
z_Target = 4.2;
r1D = linspace(1e-3,0.13,40 );
z1D = linspace(z_Dump,z_Target,300);

% Assignment of currents per power supply:
% =========================================================================
% Peak heat flux shot:
ii = 1;
shot{ii} = 29655;
label{ii} = 'shot: ';
coilCurrents{ii}.TR1 = 530;
coilCurrents{ii}.TR2 = 2200;
coilCurrents{ii}.PS1 = 3500;
coilCurrents{ii}.PS2 = 1630;
coilCurrents{ii}.PS3 = 550;

% Calculate the magnetic field and magnetic vector potential:
% =========================================================================
dum1 = tic;
disp(['Calculating magnetic field for ',num2str(numel(coilCurrents)),' cases ...'])
for ii = 1:numel(coilCurrents)
    % Create "coil" structure":
    [coil] = CreateCoilStructure(coilSetup,coilCurrents{ii});
    % Calculate magnetic field and vector potential:
    [Br2D,Bz2D,~,phi2D{ii},z2D,r2D] = CalculateMagField(coil,z1D,r1D,'grid');
    % Magnetic field magnitude:
    B2D{ii} = sqrt(Br2D.*Br2D  + Bz2D.*Bz2D);
    disp(['Case ',num2str(ii),' complete!'])
end
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

% Calculate the magnetic field on-axis all the way to pump:
% =========================================================================
z_1 = 0;
z_2 = 6;
r1Db = linspace(1e-3,0.05,5 );
z1Db = linspace(z_1,z_2,300);
for ii = 1:numel(coilCurrents)
    % Calculate magnetic field:
    [Br2Db,Bz2Db,~,~,~,~] = CalculateMagField(coil,z1Db,r1Db,'grid');
    % Magnetic field magnitude:
    B2Db{ii} = sqrt(Br2Db.*Br2Db  + Bz2Db.*Bz2Db);
end
B_onAxis = B2Db{ii}(:,1);

% Draw Proto-MPEX vacuum vessel:
% =========================================================================
DrawVacuumVessel_Conf_G_2021_04_15

% Calculate reference flux:
% =========================================================================
% Define reference flux based on phi at vacuum vessel boundary
for ii = 1:numel(coilCurrents)
    % Select region of interest:
    rng_jj = find(vessel_1_U.z > 0.5 & vessel_1_U.z < 3.7);
    % Initialize variable:
    phiBoundary = ones(size(vessel_1_U.z));
    % Interpolate phi along vaccum vessel contour
    zq = vessel_1_U.z(rng_jj);
    rq = vessel_1_U.r(rng_jj);
    a = interp2(z1D,r1D,phi2D{ii}',zq,rq);
    phiBoundary(rng_jj) = a;
    % Find location of minimim phi along contour
    [~,jj] = min(phiBoundary);
    % Physical location of limit:
    rlimit(ii) = vessel_1_U.r(jj);
    zlimit(ii) = vessel_1_U.z(jj);
    % Extract reference magnetic flux at the limiting location:
    nr = find(r1D > rlimit(ii),1,'first');
    nz = find(z1D > zlimit(ii),1);
    phi0 = interp2(z1D,r1D,phi2D{ii}',zlimit(ii),rlimit(ii));
    % Flux coordinate
    xi{ii} = phi2D{ii}/phi0;
end

% Calculate magnetic flux lines:
% =========================================================================
% Magnetic field field lines up to the plasma edge
for ii = 1:numel(coilCurrents)
    % Define the number of flux lines to plot:
    xi_lines = linspace(1e-2,1,1);
    % Calculate the flux line trajectory r(z) in physical space:
    for jj = 1:numel(xi_lines)
        C = contour(z2D,r2D,xi{ii},[1,1]*xi_lines(jj));
        z_fluxline{ii}{jj} = C(1,2:end);
        r_fluxline{ii}{jj} = C(2,2:end);
    end
end
close(gcf)

% Calculate cyclotron regions:
% =========================================================================
% Flag to draw resonance layers:
drawResLayer = 1;

% RF frequency:
f_RF = [28]*1e9;

% Particle mass:
m_e = 9.1094e-31;

% Harmonics:
n_harmonic = [2,3];

% nth harmonic cyclotron resonance: 
e_c = 1.6020e-19;
B_res = (2*pi*f_RF*m_e/e_c)./n_harmonic;

% How to draw doppler shifted resonance layers?

% Define contour lines:
rng_z = find(z1D > 3.1 & z1D < 3.35);
rng_r = find(r1D > 0 & r1D < 0.15);
for ii = 1:numel(coilCurrents)
    for nn = 1:numel(n_harmonic)
        dum1 = contour(z2D(rng_z,rng_r),r2D(rng_z,rng_r),B2D{ii}(rng_z,rng_r),[1,1]*B_res(nn));
        z_resLayer{nn}{ii} = dum1(1,2:end);
        r_resLayer{nn}{ii} = dum1(2,2:end);
    end
end
close all

%% Figure_0a: Whole device

% Magnetic flux mapping:
% =========================================================================
figure('color','w')
hfig(1) = gcf;
hold on
% Magnetic coils:
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,+coil{ii}.rfil,'r.');
    plot(coil{ii}.zfil,-coil{ii}.rfil,'r.');
end
% Flux lines:
lineColor = {'k','r','g','bl','m','c'};
for ii = 1:numel(coilCurrents)
    for jj = 1:1:numel(xi_lines)
        dum1 = plot((z_fluxline{ii}{jj}),+r_fluxline{ii}{jj},lineColor{ii});
        dum2 = plot((z_fluxline{ii}{jj}),-r_fluxline{ii}{jj},lineColor{ii});
        if jj == numel(xi_lines)
            set(dum1,'LineStyle','-','LineWidth',2)
            set(dum2,'LineStyle','-','LineWidth',2)
            hFlux(ii) = [dum1];
            legendText{ii} = [label{ii},', ',num2str(shot{ii})];
        end
    end
end

% Cyclotron resonances:
color = {'r','g','m','c','k'};
if drawResLayer
    for ii = 1:numel(coilCurrents)
        for nn = 1:numel(n_harmonic)
            plot(z_resLayer{nn}{ii},+r_resLayer{nn}{ii},color{nn},'Marker','.','LineStyle','none')
            plot(z_resLayer{nn}{ii},-r_resLayer{nn}{ii},color{nn},'Marker','.','LineStyle','none')
        end
    end
end

% Draw vacuum vessel:
% =========================================================================
% Vacuuum vessel:
plot(vessel_1_U.z,+vessel_1_U.r,'r','LineWidth',1)
plot(vessel_1_L.z,-vessel_1_L.r,'r','LineWidth',1)
plot(vessel_0_U.z,+vessel_0_U.r,'k-','LineWidth',1)
plot(vessel_0_L.z,-vessel_0_L.r,'k-','LineWidth',1)

% Helicon window:
w = diff(heliconWindow.z(2:3));
h = diff(heliconWindow.r(3:end));
x = heliconWindow.z(1);
y = heliconWindow.r;
pos = [x +y(2) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'g';
hrect.EdgeColor = 'g';
pos = [x -y(1) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'g';
hrect.EdgeColor = 'g';

% ICRF window:
w = diff(ichSleeve.z(2:3));
h = diff(heliconWindow.r(3:end));
x = ichSleeve.z(1);
y = ichSleeve.r;
pos = [x +y(2) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'c';
hrect.EdgeColor = 'c';
pos = [x -y(2)-h w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'c';
hrect.EdgeColor = 'c';

% Limiter:
w = diff(limiter.z(2:3));
h = diff(limiter.r(3:end));
x = limiter.z(1);
y = limiter.r;
pos = [x +y(2) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'r';
hrect.EdgeColor = 'r';
pos = [x -y(1) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'r';
hrect.EdgeColor = 'r';

% Limiting location:
for ii = 1:numel(coilCurrents)
    hdum1 = plot(zlimit(ii),+rlimit(ii));
    hdum2 = plot(zlimit(ii),-rlimit(ii));
    set(hdum1,'Marker','sq','MarkerFaceColor',lineColor{ii},'Color',lineColor{ii},'MarkerSize',5)
    set(hdum2,'Marker','sq','MarkerFaceColor',lineColor{ii},'Color',lineColor{ii},'MarkerSize',5)
end

% Target:
hT = line(z_Target*[1,1],0.045*[-1,+1]);
set(hT,'color','k','LineWidth',4)
    
% Dump:
hT = line(z_Dump*[1,1],[-0.2,0.2]);
set(hT,'color','k','LineWidth',4)

% Formatting:
set(gca,'FontName','times')
zoomType = 1;
switch zoomType
    case 1
        set(gcf,'position',[203.0000  203.6667  943.3333  328.6667])
        xlim([0.4,5.7])
        ylim(0.4*[-1,+1])
end
box on
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('y [m]','Interpreter','Latex','FontSize',13)
grid off

legendType = 2;
switch legendType
    case 1
        hLeg = legend(hFlux,legendText);
        set(hLeg,'FontSize',11,'Location','NorthWest')
    case 2
        hLeg = title(['shot: ',num2str(shot{ii})]);
        set(hLeg,'FontSize',11,'Interpreter','Latex')
end

% Arrows:
% =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 11;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% Dump plate:
fields.String = 'A';
x =      [0.65 ,0.53];
y =      [0.25,0.17];
hta = myTextArrow(gca,x,y,fields);

% Vacuum boundary:
fields.String = 'B';
x =      [0.9,0.9   ];
y =      [0.25,0.075];
hta = myTextArrow(gca,x,y,fields);

% LUFS:
fields.String = 'C';
x =      [1.16,1.16];
y =      [0.25 ,0.035];
hta = myTextArrow(gca,x,y,fields);

% Magnetic coil:
fields.String = 'D';
x =  0.8+[0.6 ,0.53];
y =      [0.25,0.14  ];
hta = myTextArrow(gca,x,y,fields);

% Helicon Window:
fields.String = 'E';
x = 1.75*[1  ,1   ];
y =      [0.25,0.07];
hta = myTextArrow(gca,x,y,fields);

% Limiter:
fields.String = 'F';
x = 2.05*[1  ,1   ];
y =      [0.25,0.07];
hta = myTextArrow(gca,x,y,fields);

% Skimmer 1:
fields.String = 'G1';
x = 2.265*[1  ,1   ];
y =       [0.25,0.07];
hta = myTextArrow(gca,x,y,fields);

% Skimmer 2:
fields.String = 'G2';
x =       [2.7,2.88];
y =       [0.2,0.07];
hta = myTextArrow(gca,x,y,fields);

% Strike point:
fields.String = 'H';
x =       [1.9,zlimit(end)-0.02];
y =       [0.01,rlimit(end)-0.01];
hta = myTextArrow(gca,x,y,fields);

% Gas injection:
fields.String = 'I';
x = 1.50*[1    ,1   ];
y =      [-0.2,-0.07];
hta = myTextArrow(gca,x,y,fields);
x = 1.16*[1    ,1   ];
y =      [-0.2,-0.07];
hta = myTextArrow(gca,x,y,fields);

% Turbo molecular pump 1 (TMPs):
fields.String = 'J1';
x = 2.67*[1    ,1   ];
y =      [-0.25,-0.4];
hta = myTextArrow(gca,x,y,fields);

% Turbo molecular pump 2 (TMPs):
fields.String = 'J2';
x =      [5.23  ,5.23   ];
y =      [-0.15,-0.4];
hta = myTextArrow(gca,x,y,fields);

% ICRF window:
fields.String = 'K';
x =   3.5*[1  ,1   ];
y =       [0.25,0.04];
hta = myTextArrow(gca,x,y,fields);

% Target plate:
fields.String = 'L';
x = [4.1 ,4.2];
y = [0.25,0.06];
hta = myTextArrow(gca,x,y,fields);

% 28 GHz injection port:
fields.String = 'M';
x = [3.25 ,3.25];
y = [-0.3,-0.2];
hta = myTextArrow(gca,x,y,fields);

% 2nd harmonic resonance:
fields.String = 'N1';
x = [3.0 ,3.13];
y = [-0.3 ,-0.07];
hta = myTextArrow(gca,x,y,fields);

% 3nd harmonic resonance:
fields.String = 'N2';
x = [3.24,3.207];
y = [0.26,0.08];
hta = myTextArrow(gca,x,y,fields);

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0a_WholeDevice';
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Save data: LUFS
% =========================================================================
varList = {'z_fluxline','r_fluxline','coil','vessel_1_U','z_resLayer','r_resLayer','ichSleeve','heliconWindow','z1Db','B_onAxis'};
fileName = 'step_0a_Profile_A_LUFS.mat';
save(fileName,varList{:});

%% Figure_0b: Magnetic field profile

% Magnetic field profiles:
% =========================================================================
% fontSize

figure('color','w')
hfig(2) = gcf;
hold on
for ii = 1:numel(coilCurrents)
    hBz(ii) = plot(z1Db,B_onAxis,lineColor{ii},'LineWidth',2);
end
box on 
grid on
set(gca,'PlotBoxAspectRatio',[2 1 1])
set(gca,'FontName','times')
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('B$_0$ [T]','Interpreter','Latex','FontSize',13)

set(gcf,'position',[203.0000  203.6667  943.3333  328.6667])
set(gca,'PlotBoxAspectRatio',[1.0000    0.3573/2    0.3573])
xlim([0.4,5.7])
ylim([0,0.6])

% Target:
hT = line(z_Target*[1,1],[0,0.3]);
set(hT,'color','k','LineWidth',4)

% Dump:
hT = line(z_Dump*[1,1],[0,0.3]);
set(hT,'color','k','LineWidth',4)

% Coils:
for cc = 1:numel(coil)
    x = coil{cc}.z - 0.5*coil{cc}.dz;
    y = 0;
    w = coil{cc}.dz;
    h = 0.08;
    hc(cc) = rectangle('position',[x y w h]);
    
    if cc > 9
        dx = 0.008;
    else
        dx = 0.025;
    end
    ht(cc) = text(x + dx,0.04,num2str(cc));
end

set(hc,'LineWidth',2,'EdgeColor','r')
set(ht,'color','r','FontWeight','bold','BackgroundColor','none','fontSize',8)


% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0a_MagneticfieldProfile';
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Figure 0c: ECH heating region

try
    clear ray
end

% Copy figure:
% =========================================================================
figure(hfig(1));
hold on
a1 = [gca];
f2 = figure('color','w');
a2 = copyobj(a1,f2);

% Formatting:
% =========================================================================
axis image
ylim([-0.3,0.3])
xlim([3,3.5])
title([])
set(gca,'FontSize',12)

% Vacuum vessel:
% =========================================================================
plot(vessel_0_U.z,+vessel_0_U.r,'k-','LineWidth',2)
plot(vessel_0_L.z,-vessel_0_L.r,'k-','LineWidth',2)

% ICRF window:
% =========================================================================
w = diff(ichSleeve.z(2:3));
h = diff(heliconWindow.r(3:end));
x = ichSleeve.z(1);
y = ichSleeve.r;
pos = [x +y(2) w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'c';
hrect.EdgeColor = 'c';
pos = [x -y(2)-h w h];
hrect = rectangle('Position',pos);
hrect.FaceColor = 'c';
hrect.EdgeColor = 'c';

% Draw waveguide:
% =========================================================================
p0.x = 3.25;
p0.y = -6.16*in2m;
rad = 1.3*0.5*63.5/1000;
waveguide = DrawCircle(p0,rad,100);

% Rays:
% =========================================================================
theta = 30 + 90;
nRays = 11;
dx     = linspace(-rad,rad,nRays)*cosd(30)*0.9;
dy     = linspace(-rad,rad,nRays)*sind(30)*0.9;
dtheta = linspace(-9,9,nRays);

rng_z = find(z_fluxline{1}{1} > 3.1 & z_fluxline{1}{1} < 3.4);
absorptionSurface.x = +z_fluxline{1}{1}(fliplr(rng_z));
absorptionSurface.y = -r_fluxline{1}{1}(fliplr(rng_z));

for jj = 1:numel(dy)
    ray{jj}.x(1)  = p0.x + dx(jj);
    ray{jj}.y(1)  = p0.y + dy(jj);
    ray{jj}.theta(1) = theta + dtheta(jj);  
    ray{jj} = intersectPoint(ray{jj},absorptionSurface);
    
    plot(ray{jj}.x,ray{jj}.y,'g.-','lineWidth',1)
end

for jj = 1:numel(dy)
end

% Waveguide formatting:
% =========================================================================
fill(waveguide.x,waveguide.y,'w');
plot(waveguide.x,waveguide.y,'k','LineWidth',2);
hfig(3) = gcf;

% Sectional view:
% =========================================================================
hl(1) = line([p0.x,p0.x],[-0.29,0.29],'color','k','LineWidth',1,'LineStyle','--');

% Arrows:
% =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 11;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% sectional plane:
fields.String = 'A';
x =      [p0.x + 0.03 ,p0.x];
y =      [+0.275 ,+0.275];
hta = myTextArrow(gca,x,y,fields);
x =      [p0.x + 0.03 ,p0.x];
y =      [-0.275 ,-0.275];
hta = myTextArrow(gca,x,y,fields);

% plasma LUFS:
fields.String = 'B';
x =      [3.12,3.12];
y =      [0.06,0.025];
hta = myTextArrow(gca,x,y,fields);

% 2nd harmonic resonance:
fields.String = 'C1';
x =      [3.15,3.175];
y =      [-0.06,-0.06];
hta = myTextArrow(gca,x,y,fields);

% 3rd harmonic resonance:
fields.String = 'C2';
x =      [3.208,3.208]+0.014;
y =      [0.15,0.07708];
hta = myTextArrow(gca,x,y,fields)

fields.String = [];
x =      [3.208,3.28]+0.014;
y =      [0.15,0.07708];
hta = myTextArrow(gca,x,y,fields)

% Wavguide:
fields.String = 'D';
x =      [3.28,3.26];
y =      [-0.1,-0.14];
hta = myTextArrow(gca,x,y,fields);

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0a_ECH_Region';
    % PDF figure:
    exportgraphics(hfig(3),[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(hfig(3),[figureName,'.tiff'],'Resolution',600) 
end

%% Figure 0d: Target-IR camera setup

% Copy figure:
% =========================================================================
figure(hfig(1));
hold on
a1 = [gca];
f2 = figure('color','w');
a2 = copyobj(a1,f2);

% Formatting:
axis image
ylim([-0.4,0.4])
xlim([3.8,5.85])
title([])
set(gca,'FontSize',11,'YTick',[-0.3:0.1:0.3])
title(['Side-view'],'Interpreter','Latex','FontSize',13)

% Vacuum vessel:
plot(vessel_0_U.z,+vessel_0_U.r,'k-','LineWidth',2)
plot(vessel_0_L.z,-vessel_0_L.r,'k-','LineWidth',2)

% ZnSe window:
hr(1) = DrawRectangle(gca,4,[5.546           ,0],0.5*in2m,2.75*in2m);
hr(2) = DrawRectangle(gca,4,[5.546 + 0.5*in2m,0],0.5*in2m,2.75*in2m);
set(hr,'LineWidth',1)
set(hr(2),'LineWidth',1,'Edgecolor','k')

% IR camera:
hc(1) = DrawRectangle(gca,4,[5.546 + 1.5*in2m,0],0.5*in2m,2.75*in2m);
hc(2) = DrawRectangle(gca,4,[5.546 + 2.0*in2m,0],0.2,3.0*in2m);
set(hc,'LineWidth',1)

% Lines of sight:
hls(1) = line([4.022,5.546],[+0.1085,0]);
hls(2) = line([4.022,5.546],[-0.1085,0]);
set(hls(1:2),'LineWidth',1,'color','bl')
hls(3) = line([4.2  ,5.546],[0,0]);
set(hls(3),'LineWidth',1,'color','bl','LineStyle','--')

% Arrows:
% =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 11;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% Plasma:
fields.String = 'A';
x =      [4.07,4.07];
y =      [0.14,0.02];
hta = myTextArrow(gca,x,y,fields);

% Target:
fields.String = 'B';
x =      [4.20,4.20];
y =      [0.14,0.02];
hta = myTextArrow(gca,x,y,fields);

% Vacuum boundary:
fields.String = 'C';
x =      [4.6 ,4.6 ];
y =      [0.14,0.065];
hta = myTextArrow(gca,x,y,fields);

% IR FOV:
fields.String = 'D';
x =      [4.8 ,4.8 ];
y =      [0.14,0.025];
hta = myTextArrow(gca,x,y,fields);

% ZnSe window:
fields.String = 'E';
x =      [5.45,5.52];
y =      [0.14,0.02];
hta = myTextArrow(gca,x,y,fields);

% IR camera:
fields.String = 'F';
x =      [5.7 ,5.7];
y =      [0.14,0.02];
hta = myTextArrow(gca,x,y,fields);

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0a_IR_section';
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Functions:
function circle = DrawCircle(p0,rad,N)
% p0: vector of center of circle
% r : radius of circle
% N: number of points

t = linspace(0,2*pi,N);
r.x = rad*cos(t);
r.y = rad*sin(t);

circle.x = p0.x + r.x;
circle.y = p0.y + r.y;
end

function ray = intersectPoint(ray,surface)
% Ray conditions:
n  = numel(ray.y);
m0 = tand(ray.theta(n));

% Starting location of input ray:
y0 = ray.y(n);
x0 = ray.x(n);

% Components of the unit vector of input ray:
cos0 = 1/sqrt(1 + m0^2);
sin0 = m0/sqrt(1 + m0^2);

% Input unit vector:
r0 = [cos0;sin0];

% Surface conditions:
mi = diff(surface.y)./diff(surface.x);
xi = surface.x;
yi = surface.y;

% Calculate intersection points:
for ii = 1:(numel(mi) - 1)
    MI = [-m0 1; -mi(ii) 1];
    r  = inv(MI)*[y0 - m0*x0; yi(ii) - mi(ii)*xi(ii)];
    xI(ii) = r(1);
    yI(ii) = r(2);    
end

for ii = 1:(numel(mi) - 1)
    if ( xI(ii) >= xi(ii) && xI(ii) < xi(ii+1) )
       
        ray.x(n+1) = xI(ii);
        ray.y(n+1) = yI(ii);
        
        % Slope of intersecting segment:
        m1 = mi(ii);
        
        % Components of the unit vector along intersecting segment:
        cos1 = 1/sqrt(1 + m1^2);
        sin1 = m1/sqrt(1 + m1^2);
      
        % Rotation matrix:
        R1 = [+cos1, -sin1 ;...
              +sin1, +cos1];
        
        % Express ray unit vector in terms of segment's coordinate system:
        r1 = transpose(R1)*r0;
        
        % Reflect "y'" component of unit vector:
        r1(2) = -r1(2);
        
        % Transform back to standard coordinate system:
        r2 = R1*r1;
        
        % Calculate slope of reflecting ray:
        ray.theta(n+1) = atan2d(r2(2),r2(1));
        break
    end
end
end

function hij = DrawLine(ax,pi,pj,lineColor,LineWidth,LineStyle)
% Draws a line between point pi and pj
X = [pi(1),pj(1)];
Y = [pi(2),pj(2)];
hij = line(ax,X,Y);
set(hij,'Color',lineColor,'LineWidth',LineWidth,'LineStyle',LineStyle);
end

function hr = DrawRectangle(ax,datum,pdat,w,h)
% Select datum:
% p1 --- p2 --- p3
% |      |      |
% p4 --- p5 --- p6
% |      |      |
% p7 --- p8 --- p9

switch datum
    case {1,2,3}
        dy = - 1.0*h;
    case {4,5,6}
        dy = - 0.5*h;
    case {7,8,9}
        dy = - 0.0*h;
end

switch datum
    case {3,6,9}
        dx = - 1.0*w;
    case {2,5,8}
        dx = - 0.5*w;
    case {1,4,7}
        dx = - 0.0*w;
end

x = pdat(1) + dx;
y = pdat(2) + dy;

hr = rectangle(ax,'Position',[x,y,w,h]);
% set(hr,'Color',lineColor,'LineWidth',LineWidth,'LineStyle',LineStyle);
end