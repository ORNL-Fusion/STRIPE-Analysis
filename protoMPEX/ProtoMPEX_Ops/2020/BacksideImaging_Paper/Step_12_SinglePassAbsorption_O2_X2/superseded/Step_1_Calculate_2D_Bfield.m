% Step 1: Produce the magnetic field profiles used during the magnetic
% field optimization experiments

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
r1D = linspace(1E-3,0.2,200);
z1D = linspace(3,3.5,300)';

% Assignment of currents per power supply:
% =========================================================================
shotType = 1;
switch shotType
    case 1
        % 29778, PS2: 1630 A
        coilCurrents{1}.TR1 = 530;
        coilCurrents{1}.TR2 = 2200;
        coilCurrents{1}.PS1 = 3500;
        coilCurrents{1}.PS2 = 1630;
        coilCurrents{1}.PS3 = 650;
    case 2
        % 29775, PS2: 2360 A
        coilCurrents{2}.TR1 = 530;
        coilCurrents{2}.TR2 = 2200;
        coilCurrents{2}.PS1 = 3500;
        coilCurrents{2}.PS2 = 2360;
        coilCurrents{2}.PS3 = 650;
    case 3
        % 29771, PS2: 3370 A
        coilCurrents{3}.TR1 = 530;
        coilCurrents{3}.TR2 = 2500;
        coilCurrents{3}.PS1 = 3500;
        coilCurrents{3}.PS2 = 3370;
        coilCurrents{3}.PS3 = 650;
end

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
    
    % Assemble column:
    B{ii}   = [fliplr(B2D{ii}),B2D{ii}];   
    phi{ii} = [fliplr(phi2D{ii}),phi2D{ii}];    

end
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

coord.z = z1D;
coord.r = linspace(-max(r1D),max(r1D),2*numel(r1D));


%% 2nd harmonic resoance:
% =========================================================================
% RF frequency:
f_RF = [28]*1e9;

% Particle mass:
m_e = 9.1094e-31;

% Harmonics:
n_harmonic = [2];

% nth harmonic cyclotron resonance: 
e_c = 1.6020e-19;
B_res = (2*pi*f_RF*m_e/e_c)./n_harmonic;

% How to draw doppler shifted resonance layers?

% Define contour lines:
rng_z = find(coord.z > 3.1 & coord.z < 3.3);
rng_r = find(coord.r > -0.15 & coord.r < 0.15);
for ii = 1:numel(coilCurrents)
    for nn = 1:numel(n_harmonic)
%         dum1 = contour(z2D(rng_z,rng_r),r2D(rng_z,rng_r),B{ii}(rng_z,rng_r),[1,1]*B_res(nn));
        dum1 = contour(coord.z(rng_z),coord.r(rng_r),B{ii}(rng_z,rng_r)',[1,1]*B_res(nn));
        z_resLayer{nn}{ii} = dum1(1,2:end);        
        r_resLayer{nn}{ii} = dum1(2,2:end);
    end
end

for ii = 1:numel(coilCurrents)
    phi0 = 7E-4;
    % Flux coordinate
    xi{ii} = phi{ii}/phi0;
end

% Calculate magnetic flux lines:
% =========================================================================
% Magnetic field field lines up to the plasma edge
for ii = 1:numel(coilCurrents)
    % Define the number of flux lines to plot:
    xi_lines = linspace(1e-2,1,1);
    % Calculate the flux line trajectory r(z) in physical space:
    for jj = 1:numel(xi_lines)
        C = contour(coord.z,coord.r,xi{ii}',[1,1]*xi_lines(jj));
        z_fluxline{ii}{jj} = C(1,2:end);
        r_fluxline{ii}{jj} = C(2,2:end);
    end
end
close(gcf)

% Calculate density profile:
neMax  = 1e18;
neMin = neMax/100;
xip = 0.6;
dxi = 0.19;
R = (tanh((xi{1} + xip)/dxi) - tanh((xi{1} - xip)/dxi) )/2;
ne = (neMax - neMin)*R + neMin;


%% Plot data:
close all

figure;
mesh(coord.z,coord.r,ne')

figure('color','w')
plot(z1D,Bz2D(:,1))
hold on
% ni = find(Bz2D(:,1) >

figure('color','w')
hold on
hM = surf(coord.z,coord.r,B{1}','LineStyle','none');
zlim([0,0.8])
caxis([0,0.8])
view([0 90])

figure('color','w')
drawResLayer = 1;

hold on
contourf(z1D,+r1D,B2D{1}',400,'LineStyle','--','LineColor','k');
contourf(z1D,-r1D,B2D{1}',400,'LineStyle','--','LineColor','k');
caxis([0,0.8])
set(gca,'fontName','times','fontSize',12)

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

% Magnetic coils:
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,+coil{ii}.rfil,'ro','MarkerSize',4,'MarkerFaceColor','r');
    plot(coil{ii}.zfil,-coil{ii}.rfil,'ro','MarkerSize',4,'MarkerFaceColor','r');
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
%             legendText{ii} = [label{ii},', ',num2str(shot{ii})];
        end
    end
end

xlim([min(z1D),max(z1D)]);
ylim([min(-r1D),max(r1D)]);

% Path of Xmode beam:
beam.Npoints = 100;
beam.z1 = +3.259;
beam.z0 = +3.179;
beam.r1 = -0.16;
beam.r0 = 0.11;

% From experiment:
beam.z0 = +3.184;
beam.r0 = -0.052;
beam.z1 = +3.135;
beam.r1 = 0.096;

% Through 2nd harmonic:
beam.z0 = +3.164;
beam.r0 = -0.061;
beam.z1 = +3.08;
beam.r1 = 0.079;

% Through 2nd harmonic and launch from experiment:
beam.zS = +3.25;
beam.rS = -0.1861;

beam.z0 = +3.165;
beam.r0 = -0.053;
beam.z1 = +3.08;
beam.r1 = 0.079;

% beam.z1 = 3.121;
% beam.z0 = 3.121;
% beam.r1 = +0*0.05;
% beam.r0 = -0.1;

beam.z = linspace(beam.z0,beam.z1,beam.Npoints);
beam.r = linspace(beam.r0,beam.r1,beam.Npoints);

plot(beam.z,beam.r,'k','LineWidth',3)
line([beam.zS,beam.z1],[beam.rS,beam.r1])

% Calculate the ne along path:
beam.path.ne = interp2(coord.z,coord.r,ne',beam.z,beam.r);
beam.path.B = interp2(coord.z,coord.r,B{1}',beam.z,beam.r);
beam.path.s  = sqrt((beam.z-beam.z0).^2 + (beam.r - beam.r0).^2);

figure; 
subplot(2,1,1)
plot(beam.path.s,beam.path.ne);

subplot(2,1,2)
plot(beam.path.s,beam.path.B);

%% Calculate the absorption:

w   = 2*pi*f_RF;
wp  = w_pe(beam.path.ne);
wc  = w_ce(beam.path.B);
wuh = sqrt(wp.^2 + wc.^2);
Te  = 4;
n   = 2;
d   = sqrt( (n*wc/w).^2 - 1 );
G   = Te*e_c/(m_e*c_light^2);

t1 = sqrt(2*pi)/(30*c_light);
t2 = (wp.^2).*sqrt( ((wp.^2) - (w^2)).^2 - (wc*w).^2 )./((w^3)*wc);
t3 = ((w^2 - (wc.^2)).^2).*((w^2 + w*wc - wp.^2).^2).*(G^(-5/2))./((w+wc).^2).*((w^2 - (wuh).^2).^(5/2));
t4 = (d.^5).*exp( -(d.^2)./(2*G) ).*theta(d.^2);

kIa = t1.*t2.*t3.*t4; 
eta = 2*cumtrapz(beam.path.s,kIa);
P   = exp(-eta);

figure;
plot(beam.path.s,P)
ylim([0,1])

figure;
mesh(coord.z,coord.r,2*e_c*B{1}'/(w*m_e)  - 1);
zlim([0,1])
view([0,90])

% return

%% Scale lenghts:
close all

x = linspace(0,0.2);
L = 0.55; 
B0 = beam.path.B(1);

y = B0*(1 + (x/L));

figure; 
hold on
plot(x,y)
plot(beam.path.s,beam.path.B)
yyaxis right
plot(beam.path.s,beam.path.ne)

w   = 2*pi*f_RF;
wp  = w_pe(1e18);
wc  = w_ce(beam.path.B);
wc0 = w_ce(0.5);
Te  = logspace(0,4);

% Second harmonic:
eta_O2 = pi*((wp.^2)/wc0).*(L/c_light).*((1 - 0.25*(wp/wc0).^2).^(3/2)).*(e_c*Te/(m_e*c_light^2)).^(2);
P =  exp(-eta_O2);

figure; 
hold on
plot(Te,eta_O2)
title('optical length [m], Chu 1982')
set(gca,'Xscale','log')
set(gca,'Yscale','log')

%% Bornatici 1983: Optical thickness and absorption

close all

t0  = 60;
n = 2;
f_RF = 28e9;
ne0 = 1e18;
B0  = 0.5;
Te0 = 4; 
Lb = 0.55;
L0 = c_light/f_RF;

vT = @(Te) sqrt(e_c*Te/m_e);
At = @(n) (pi*pi*n^(2*(n-1)))/((2^(n-1))*factorial(n-1));
Bt = @(ne,B) (w_pe(ne)/w_ce(B))^2;
Ct = @(n,Te) (vT(Te)/c_light).^(2*(n-1));
Dt = @(n,t) sind(t)^(2*(n-1)); 
Et = @(t) 1 + (cosd(t)^2);

A1  = @(n,t) (sind(t)^4)/(4*n) + cosd(t)^2;
A2  = @(n,t) sqrt( ((sind(t)^4)/(4*n*n)) + cosd(t)^2  )*Et(t);
mu_OX = @(n,t,OX) 0.5 - OX*A1(n,t)/A2(n,t); 

tau_OX = @(n,t,Te,ne,B,LB,L,OX) At(n)*Bt(ne,B)*Ct(n,Te)*Dt(n,t)*Et(t)*mu_OX(n,t,OX)*LB/L;

Te = logspace(0,log10(1e3));

tau_O2 = tau_OX(2,t0,Te,ne0,B0,Lb,L0,+1);
tau_X2 = tau_OX(2,t0,Te,ne0,B0,Lb,L0,-1);
tau_O3 = tau_OX(3,t0,Te,ne0,B0,Lb,L0,-1);

figure; 
hold on
htau(1) = plot(Te,tau_O2,'k');
htau(2) = plot(Te,tau_X2,'r');
htau(3) = plot(Te,tau_O3,'g');
title('optical length [m]')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
ylim([1e-8,1])
legend(htau,'O2','X2','O3')

Pabs_O2 = 1 - exp(-tau_O2);
Pabs_X2 = 1 - exp(-tau_X2);
Pabs_O3 = 1 - exp(-tau_O3);

figure; 
hold on
plot(Te,Pabs_O2)
plot(Te,Pabs_X2,'r')
plot(Te,Pabs_O3,'g')
title('P_{abs}/P_{inj} [m]')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
ylim([1e-8,1])


%%

X = (wp/w).^2;
Y = wc/w;

t1 = (X./Y).*sqrt( (X - 1).^2 - Y.^2);
t2 = ((1 - Y.^2).^2).*(1 + Y - X).^2;
t3 = ((1 + Y).^2).*(1 - Y.^2 - X).^(5/2);
t4 = (d.^5).*exp( -(d.^2)./(2*G) ).*theta(d.^2);
kIb = (w/c_light)*(sqrt(2*pi)/30).*t1.*(t2./t3)*(G^(-5/2)).*t4


%% functions:
function y = theta(x)
   y = zeros(size(x));
   rng = find(x >= 0);
   y(rng) = 1;
end