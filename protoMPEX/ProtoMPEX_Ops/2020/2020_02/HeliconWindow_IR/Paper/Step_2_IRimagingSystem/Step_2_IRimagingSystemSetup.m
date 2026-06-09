% Describes the IR imaging system setup in Proto-MPEX using the IR camera
% and the stainless steel mirror:

% Written by J.F. Caneses Marin
% Created on 2020-09-03

%% SECTION 1: Read "CoilSetup" spreadsheet
clearvars
clc
close all
 
% =========================================================================
% Magnetic configuration of interest:
confType = 'conf_G';

% =========================================================================
% Read "CoilSetup" spreadsheet:
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet',confType);
disp('Reading complete!')
toc(dum1);

% =========================================================================
% Display "coilSetup" on CLI:
coilSetup

% #########################################################################
% IMPORTANT QUANTITIES FROM THIS SECTION:
% #########################################################################
% - coilSetup: structure that contains all the geometrical info of the
% coils
% - coilCurrents: strucuture that contains the power supply currents
% - coil: object that describes the specific coil setup and is the input
% for the magnetic field calculator function

%% SECTION 2: Calculate magnetic field
% =========================================================================
% Define the area to evaluate the fields at:
z_Dump = 0.5;
z_Target = 4.2;
r1D = linspace(1e-3,0.13,40 );
z1D = linspace(z_Dump-1,z_Target+1,300);

% #########################################################################
%                       INPUT FROM USER:
% #########################################################################

% Assignment of currents per power supply:
% =========================================================================

% Window-limiter:
ii = 1;
shot{ii} = 29097;
label{ii} = 'Window-limiter';
coilCurrents{ii}.TR1 = 530;
coilCurrents{ii}.TR2 = 2300;
coilCurrents{ii}.PS1 = 6800;
coilCurrents{ii}.PS2 = 4000;
coilCurrents{ii}.PS3 = 160;

% MPEX-limiter:
ii = 1;
shot{ii} = 29128;
label{ii} = 'MPEX-limiter';
coilCurrents{ii}.TR1 = 530;
coilCurrents{ii}.TR2 = 2300;
coilCurrents{ii}.PS1 = 3500;
coilCurrents{ii}.PS2 = 4000;
coilCurrents{ii}.PS3 = 430;

% #########################################################################

% #########################################################################

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

%% Draw Proto-MPEX vacuum vessel:
if ~strcmpi(confType,'conf_G')
%     error('Change "confType" to config_G')
end
DrawVacuumVessel_Conf_G

%% SECTION 3: Calculate reference flux:
% =========================================================================
% Define reference flux based on phi at vacuum vessel boundary
for ii = 1:numel(coilCurrents)
    % Select region of interest:
    rng_jj = find(vessel_1.z > 0.5 & vessel_1.z < 3.7);
    % Initialize variable:
    phiBoundary = ones(size(vessel_1.z));
    % Interpolate phi along vaccum vessel contour
    zq = vessel_1.z(rng_jj);
    rq = vessel_1.r(rng_jj);
    a = interp2(z1D,r1D,phi2D{ii}',zq,rq);
    phiBoundary(rng_jj) = a;
    % Find location of minimim phi along contour
    [~,jj] = min(phiBoundary);
    % Physical location of limit:
    rlimit(ii) = vessel_1.r(jj);
    zlimit(ii) = vessel_1.z(jj);
    % Extract reference magnetic flux at the limiting location:
    nr = find(r1D > rlimit(ii),1,'first');
    nz = find(z1D > zlimit(ii),1);
    phi0 = interp2(z1D,r1D,phi2D{ii}',zlimit(ii),rlimit(ii));
    % Flux coordinate
    xi{ii} = phi2D{ii}/phi0;
end
%% SECTION 4: MAGNETIC FIELD LINES AND PLASMA EDGE
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
clearvars ii

%% SECTION 6: CYCLOTRON RESONANCE REGIONS
% Flag to draw resonance layers:
drawResLayer = 0;

% RF frequency:
f_RF = [28]*1e9;

% Particle mass:
m_e = 9.1094e-31;

% Harmonics:
n_harmonic = [1,2,3,4,5];

% nth harmonic cyclotron resonance: 
e_c = 1.6020e-19;
B_res = (2*pi*f_RF*m_e/e_c)./n_harmonic;

% Define contour lines:
for ii = 1:numel(coilCurrents)
    for nn = 1:numel(n_harmonic)
        dum1 = contour(z2D,r2D,B2D{ii},[1,1]*B_res(nn));
        z_resLayer{nn}{ii} = dum1(1,2:end);
        r_resLayer{nn}{ii} = dum1(2,2:end);
    end
end
close all

%% SECTION 7: PLOT DATA

% Figure 1:
% =========================================================================
% Magnetic flux mapping:
figureName = 'Step_2_IRimagingSystem_OpticalSetup';
figure('color','w','Tag',figureName)
hold on
ax = gca;
% Magnetic coils:
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,+coil{ii}.rfil,'r.');
    plot(coil{ii}.zfil,-coil{ii}.rfil,'r.');
end

% Flux lines:
if 0
    lineColor = {'r','k','g','bl','m','c'};
else
    lineColor = {'k','g','bl','m','c'};
end
for ii = 1:numel(coilCurrents)
    for jj = 1:1:numel(xi_lines)
        dum1 = plot((z_fluxline{ii}{jj}),+r_fluxline{ii}{jj},lineColor{ii});
        dum2 = plot((z_fluxline{ii}{jj}),-r_fluxline{ii}{jj},lineColor{ii});
        if jj == numel(xi_lines)
            set(dum1,'LineStyle','-','LineWidth',1)
            set(dum2,'LineStyle','-','LineWidth',1)
            hFlux(ii) = [dum1];
            legendText{ii} = [label{ii},', shot: ',num2str(shot{ii})];
        end
    end
end

% Fill in the plasma volume:
ZZ = [z_fluxline{1}{1},fliplr(z_fluxline{1}{1})];
RR = [r_fluxline{1}{1},fliplr(-r_fluxline{1}{1})];
patch(ZZ,RR,[0.94,0.94,0.94])

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

% Draw hardware:
% =========================================================================
% Draw vacuuum vessel:
% -------------------------------------------------------------------------
plot(vessel_1.z,+vessel_1.r,'r','LineWidth',2)
plot(vessel_1.z,-vessel_1.r,'r','LineWidth',2)
plot(vessel_0.z,+vessel_0.r,'k-','LineWidth',2)
plot(vessel_0.z,-vessel_0.r,'k-','LineWidth',2)

% Define Helicon window:
% -------------------------------------------------------------------------
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

% Define Limiter:
% -------------------------------------------------------------------------
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

% Draw Limiting location:
% -------------------------------------------------------------------------
for ii = 1:numel(coilCurrents)
    hdum1 = plot(zlimit(ii),+rlimit(ii));
    hdum2 = plot(zlimit(ii),-rlimit(ii));
    set(hdum1,'Marker','o','MarkerFaceColor',lineColor{ii},'Color',lineColor{ii},'MarkerEdgeColor','k')
    set(hdum2,'Marker','o','MarkerFaceColor',lineColor{ii},'Color',lineColor{ii},'MarkerEdgeColor','k')
end

% Draw Target:
% -------------------------------------------------------------------------
hT = line(z_Target*[1,1],0.045*[-1,+1]);
set(hT,'color','k','LineWidth',4)
    
% Draw Dump:
% -------------------------------------------------------------------------
hT = line(z_Dump*[1,1],[0,0.2]);
set(hT,'color','k','LineWidth',4)

% Draw Spool 2.5 flanges and nipples:
% -------------------------------------------------------------------------
spool_2_u.pdat = [1.2,+rVac1]';
spool_2_d.pdat = [1.2,-rVac1]';
flange.w  = 06.98/100;
flange.h  = 01.27/100;
nipple.w  = 04.15/100;
nipple.h  = 10.16/100;
hfl(1) = DrawRectangle(ax,2,spool_2_d.pdat                        ,flange.w,flange.h);
hfl(2) = DrawRectangle(ax,2,[spool_2_d.pdat(1),hfl(1).Position(2)],flange.w,flange.h);
hfl(3) = DrawRectangle(ax,2,[spool_2_d.pdat(1),hfl(2).Position(2)],nipple.w,nipple.h);
hfl(4) = DrawRectangle(ax,2,[spool_2_d.pdat(1),hfl(3).Position(2)],flange.w,flange.h);
hfl(5) = DrawRectangle(ax,2,[spool_2_d.pdat(1),hfl(4).Position(2)],flange.w,flange.h);
hfl(6) = DrawRectangle(ax,8,spool_2_u.pdat                        ,flange.w,flange.h);
hfl(7) = DrawRectangle(ax,8,[spool_2_u.pdat(1),sum(hfl(6).Position([2,4]))],flange.w,flange.h);
set(hfl(:),'LineWidth',1,'EdgeColor','k')
set(hfl(5),'LineWidth',2,'EdgeColor','r')

% Draw camera:
% -------------------------------------------------------------------------
cam.pdat = [spool_2_u.pdat(1),hfl(5).Position(2) - 1/100]';
cam.body.w = 75/1000;
cam.body.h = 202/1000;
cam.lense.w = 67/1000;
cam.lense.h = 14/1000;
hcam(1) = DrawRectangle(ax,2,[cam.pdat]                       ,cam.lense.w,cam.lense.h);
hcam(2) = DrawRectangle(ax,2,[cam.pdat(1),hcam(1).Position(2)],cam.body.w ,cam.body.h );
set(hcam,'LineWidth',2,'EdgeColor','k')

% Define mirror position:
% -------------------------------------------------------------------------
mirror.pdat = [cam.pdat(1),0.06]';

% Define principle position:
% -------------------------------------------------------------------------
% Define camera IR chip location:
p1 = [cam.pdat(1),cam.pdat(2)-(24.5/1000) ]';
% Define center of mirror:
p2 = mirror.pdat;
% Define center of helicon window:
heliconWindow.th = heliconWindow.r(1) - heliconWindow.r(2);
p3 = [mean(heliconWindow.z),-heliconWindow.r(2)]';

% Define principle rays:
% =========================================================================
% p12 ray: Principal ray
% -------------------------------------------------------------------------
p12 = p1-p2;
a12 = mag(p12);
x12 = p12(1);
y12 = p12(2);
% p32 ray: Principal ray
% -------------------------------------------------------------------------
p32 = p3-p2;
a32 = mag(p32);
x32 = p32(1);
y32 = p32(2);

% Define mirror normal:
% -------------------------------------------------------------------------
% p24 normal:
m42 = -(a32*x12 - a12*x32)/(a32*y12 - a12*y32);
t42 = atan(m42);
p4 = p2 + (15/100)*e_hat(t42);
p42 = p4 - p2;
x42 = p42(1);
y42 = p42(2);

% Define edges of mirror:
% -------------------------------------------------------------------------
% Length along principle axis of mirror:
mirror.diam_max = 1.8*2.54/100;
% Half width of the mirror
a52 = mirror.diam_max/2;
% Angle of the mirror surface relative to horizontal
t52 = t42+pi/2;
% End points of the mirror
p5 = p2 + a52*e_hat(t52);
p6 = p2 - a52*e_hat(t52);

% Define envelope rays:
% -------------------------------------------------------------------------
% p16 and p86 ray:
p16 = p1 - p6;
aa = dot(p16,p42)/mag(p16);
p8(2) = -rVac1;
y86 = p8(2) - p6(2);
bb = y42*y86;
A = x42^2 - aa^2;
B = 2*x42*bb;
C = bb^2 - (aa*y86)^2;
x86 = (-B - sqrt(B^2 - 4*A*C))/(2*A);
p86 = [x86,y86]';
p8 = p6 + p86;

% p15 and p95 ray:
p15 = p1 - p5;
aa = dot(p15,p42)/mag(p15);
p9(2) = p3(2);
y95 = p9(2) - p5(2);
bb = y42*y95;
A = x42^2 - aa^2;
B = 2*x42*bb;
C = bb^2 - (aa*y95)^2;
x95 = (-B - sqrt(B^2 - 4*A*C))/(2*A);
p95 = [x95,y95]';
p9 = p5 + p95;

% Draw rays:
% =========================================================================
% Formating:
principalRay.LineStyle = '--';
principalRay.Color = 'bl';
SecondaryRay.LineStyle = '-';
SecondaryRay.Color = 'bl';

% Principal rays:
h12 = DrawLine(ax,p1,p2,principalRay.Color,1,principalRay.LineStyle);
h23 = DrawLine(ax,p2,p3,principalRay.Color,1,principalRay.LineStyle);

% Mirror:
% h24 = DrawLine(ax,p2,p4,'r',1,'--');
h25 = DrawLine(ax,p2,p5,'r',3,'-');
h26 = DrawLine(ax,p2,p6,'r',3,'-');

% Envelope rays:
h16 = DrawLine(ax,p1,p6,SecondaryRay.Color,2,SecondaryRay.LineStyle);
h68 = DrawLine(ax,p6,p8,SecondaryRay.Color,2,SecondaryRay.LineStyle);
h15 = DrawLine(ax,p1,p5,SecondaryRay.Color,2,SecondaryRay.LineStyle);
h95 = DrawLine(ax,p5,p9,SecondaryRay.Color,2,SecondaryRay.LineStyle);

% Figure zoom formatting:
% =========================================================================
set(gca,'FontName','times')
zoomType = 4;
switch zoomType
    case 1
        set(gca,'PlotBoxAspectRatio',[2 1.5 1])
        xlim([0.25,4.5])
        ylim(0.18*[0,+1])
    case 2
        set(gca,'PlotBoxAspectRatio',[2.3 1.5 1])
        xlim([0.8,2.7])
        ylim(0.11*[-1,+1])
        axisFontsize = 14;
        legendFontSize = 11;
        textArrowFontSize = 14;
    case 3
        set(gca,'PlotBoxAspectRatio',[2.3 1.5 1])
        xlim([0.25,4.5])
        ylim(0.35*[-1,+1])
    case 4
        axis image
        xlim([0.9,2.3])
        ylim([-0.35,0.2])
        axisFontsize = 12;
        legendFontSize = 11;
        textArrowFontSize = 12;
        set(gcf,'position',[203.0000  225.0000  900.6667  328.6667])
end
box on
xlabel('z [m]','Interpreter','Latex','FontSize',axisFontsize)
ylabel('r [m]','Interpreter','Latex','FontSize',axisFontsize)
set(gca,'FontSize',axisFontsize)
grid on

% Arrow annotations:
% =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = textArrowFontSize;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% Helicon Window:
fields.String = {'Helicon';'Window'};
x = 1.75*[1  ,1   ];
y =      [-0.2,-0.085];
hta = myTextArrow(gca,x,y,fields);

% Limiter:
fields.String = 'Limiter';
x = 1.98*[1  ,1   ];
y =      [-0.2,-0.085];
hta = myTextArrow(gca,x,y,fields);

% strike point:
fields.String = 'Strike point';
x = 0.987*zlimit*[1   ,1   ];
y =        [0.12,1.15*rlimit];
hta = myTextArrow(gca,x,y,fields);

% Camera:
fields.String = 'IR camera';
x = [1.32,1.27];
y = [-0.3,-0.3];
hta = myTextArrow(gca,x,y,fields);

% Mirror:
fields.String = 'Mirror';
x = [1.20,1.22];
y = [0.13,0.065];
hta = myTextArrow(gca,x,y,fields);

% Last Uninterrupted Flux Surface (LUFS):
fields.String = '$\Phi_0$';
fields.Interpreter = 'Latex';
x = [1.11,1.11];
y = [0.13,0.02536];
hta = myTextArrow(gca,x,y,fields);

% ZnSe :
fields.String = 'ZnSe window';
x = [+1.32,+1.26];
y = [-0.22,-0.22];
hta = myTextArrow(gca,x,y,fields);

% Electromagnets:
fields.String = {'Electro-';'magnet coil'};
x = [+1.05,+1.05];
y = [-0.25,-0.18];
hta = myTextArrow(gca,x,y,fields);

% Legend formatting:
% =========================================================================
hLeg = legend(hFlux,legendText);
set(hLeg,'FontSize',legendFontSize,'Location','SouthEast')
hLeg.Box = 'off';

% #########################################################################
% #########################################################################
% #########################################################################

% Figure 2:
% =========================================================================
% Magnetic field profiles:
figureName = 'Step_2_IRimagingSystem_FieldProfile';
figure('color','w','Tag',figureName)
hold on
for ii = 1:numel(coilCurrents)
    hBz(ii) = plot(z1D,B2D{ii}(:,1),lineColor{ii},'LineWidth',2);
end
box on 
grid on
set(gca,'PlotBoxAspectRatio',[2 1 1])
set(gca,'FontName','times')
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('B$_0$ [T]','Interpreter','Latex','FontSize',13)
xlim([0,5])

% Cyclotron resonance layer:
if drawResLayer
    for nn = 1:numel(n_harmonic)
        hR = line([0,5],[1,1]*B_res(nn));
        set(hR,'color',color{nn},'LineWidth',1)
    end
end

% Target:
hT = line(z_Target*[1,1],[0,0.3]);
set(hT,'color','k','LineWidth',4)

% Dump:
hT = line(z_Dump*[1,1],[0,0.3]);
set(hT,'color','k','LineWidth',4)

hLeg = legend(hBz,legendText);
set(hLeg,'FontSize',11,'Location','best')

%% SECTION 8: SAVE PICTURE
% =========================================================================
% Saving figure:

InputStructure.prompt = {['Would you like to save figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    figureName = 'Step_2_IRimagingSystem_OpticalSetup';
    hdum1 = findobj('Tag',figureName);
    saveas(hdum1,figureName,'tiffn')
    
%     figureName = 'Step_2_IRimagingSystem_FieldProfile';
%     hdum1 = findobj('Tag',figureName);
%     saveas(hdum1,figureName,'tiffn')
end

% =========================================================================
disp('End of script')

%% Example: drawing optical setup and rays

runExample = 0;

if runExample 
close all

% Define principle position:
% =========================================================================
% Define Camera lens location:
p1 = [+1,-5.25]';
% Define center of mirror:
p2 = [+1,+2]';
% Define center of helicon window:
p3 = [+9,-2]';

% Define principle rays:
% =========================================================================
% p12 ray: Principal ray
p12 = p1-p2;
a12 = mag(p12);
x12 = p12(1);
y12 = p12(2);
% p32 ray: Principal ray
p32 = p3-p2;
a32 = mag(p32);
x32 = p32(1);
y32 = p32(2);

% Define mirror normal:
% =========================================================================
% p24 normal:
m42 = -(a32*x12 - a12*x32)/(a32*y12 - a12*y32);
t42 = atan(m42);
p4 = p2 + 2.0*e_hat(t42);
p42 = p4 - p2;
x42 = p42(1);
y42 = p42(2);

% Define edges of mirror:
% =========================================================================
% Half width of the mirror
a52 = 0.5;
% Angle of the mirror surface relative to horizontal
t52 = t42+pi/2;
% End points of the mirror
p5 = p2 + a52*e_hat(t52);
p6 = p2 - a52*e_hat(t52);

% Define envelope rays:
% =========================================================================
% p16 and p86 ray:
p16 = p1 - p6;
aa = dot(p16,p42)/mag(p16);
p8(2) = p3(2);
y86 = p8(2) - p6(2);
bb = y42*y86;
A = x42^2 - aa^2;
B = 2*x42*bb;
C = bb^2 - (aa*y86)^2;
x86 = (-B - sqrt(B^2 - 4*A*C))/(2*A);
p86 = [x86,y86]';
p8 = p6 + p86;

% p15 and p95 ray:
p15 = p1 - p5;
aa = dot(p15,p42)/mag(p15);
p9(2) = p3(2);
y95 = p9(2) - p5(2);
bb = y42*y95;
A = x42^2 - aa^2;
B = 2*x42*bb;
C = bb^2 - (aa*y95)^2;
x95 = (-B - sqrt(B^2 - 4*A*C))/(2*A);
p95 = [x95,y95]';
p9 = p5 + p95;

% Plot data:
% =========================================================================
figure('color','w')
ax = gca;
hold on
grid on
principalRay.LineStyle = '--';
principalRay.Color = 'bl';
SecondaryRay.LineStyle = '-';
SecondaryRay.Color = 'bl';

% Draw vacuum vessel:
hV  = DrawRectangle(ax,4,[-1,0],15,4)

% Draw rays:
h12 = DrawLine(ax,p1,p2,principalRay.Color,1,principalRay.LineStyle);
h23 = DrawLine(ax,p2,p3,principalRay.Color,1,principalRay.LineStyle);
h24 = DrawLine(ax,p2,p4,'r',1,'--');
h25 = DrawLine(ax,p2,p5,'r',3,'-');
h26 = DrawLine(ax,p2,p6,'r',3,'-');
h16 = DrawLine(ax,p1,p6,SecondaryRay.Color,1,SecondaryRay.LineStyle);
h68 = DrawLine(ax,p6,p8,SecondaryRay.Color,1,SecondaryRay.LineStyle);
h15 = DrawLine(ax,p1,p5,SecondaryRay.Color,1,SecondaryRay.LineStyle);
h95 = DrawLine(ax,p5,p9,SecondaryRay.Color,1,SecondaryRay.LineStyle);

% Draw camera:
hcam(1) = DrawRectangle(ax,2,p1,1,1.25);
hcam(2) = DrawRectangle(ax,8,p1,0.75,0.25);

% Draw flange:
th_fl = 0.25;
hfl(1) = DrawRectangle(ax,2,[p1(1),hV.Position(2)    ],1.5,th_fl)
hfl(2) = DrawRectangle(ax,2,[p1(1),hfl(1).Position(2)],1.5,th_fl)
hfl(3) = DrawRectangle(ax,2,[p1(1),hfl(2).Position(2)],1.0,1.75)
hfl(4) = DrawRectangle(ax,2,[p1(1),hfl(3).Position(2)],1.5,th_fl)
hfl(5) = DrawRectangle(ax,2,[p1(1),hfl(4).Position(2)],1.5,th_fl)
set(hfl,'LineWidth',2)

% Formatting:
axis image
ylim([-7,4])
xlim([-2,13])
end

%% Functions:
% All these functions are for 2D vector algebra:

function r = mag(v)
% r: magnitude of vector v
r = sqrt(v'*v);
end

function t = arg(v)
% t: angle of vector v in radians
t = atan(v(2)/v(1));
end

function s = e_hat(t)
% t is angle in radians
s = [cos(t),sin(t)]';
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
