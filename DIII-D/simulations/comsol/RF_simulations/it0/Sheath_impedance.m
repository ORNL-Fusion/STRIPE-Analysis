function maxmeanmed=Sheath_impedance()
%% Sheath Impedance Fits

%The code is based upon work described in

% H. Kohno, J.R. Myra, A finite element procedure for radio-frequency 
%sheath-plasma interactions based on a sheath impedance model, 
%Comput. Phys. Commun. 220 (2017) 129–142. 
%doi:https://doi.org/10.1016/j.cpc.2017.06.025.
% 

% filename is input electric field file from sheath layer
% init is to set initial value of epsilon and sigma to 1 and 0 (1 is True)
% unmagnetized is to set unmagnetized (1) or magnetized (0) 
%% Code Start

%Read in table

%EdgeData=xlsread('C:\Users\cxe\Documents\Proto-MPEX\COMSOLmodel\Helicon\Inputs\Beers_Inputs\HigherTe\Beers_Helicon_3D_HighDensityHighTe_SheathEdge.xlsx');

close all;

init = 0;
unmagnetized = 0;
prefix='test2';
iteration=0;
data=readmatrix([prefix '_it' num2str(iteration)],'NumHeaderLines',9);

Inputs.theta= data(:,4);
Inputs.d= 0.002; %[m] assuming data is at middle of layer
Inputs.phi = data(:,5);
Inputs.E_norm = data(:,9); % normal to sheath layer
Inputs.Vlayer = abs(Inputs.E_norm) * Inputs.d; %2 is peak-to-peak voltage
Inputs.Density=data(:,6);
Inputs.Br = data(:,8);
Inputs.B0 = data(:,7);

figure
dtri=delaunay(Inputs.phi,Inputs.theta);
boundary=[-0.356770278297659,0.351334138834533;
0.379149237370599,0.351334138834532;
0.352729894569156,0.610560916039146;
-0.384173821346535,0.610560916039146];
okIdx=[];
for i=1:size(dtri,1)
    if ~inpolygon(mean(Inputs.phi(dtri(i,:))),mean(Inputs.theta(dtri(i,:))),boundary(:,1),boundary(:,2))
        okIdx=[okIdx i];
    end
end
triplot(dtri(okIdx,:),Inputs.phi,Inputs.theta)
dtri=dtri(okIdx,:);

figure
trisurf(dtri,Inputs.phi,Inputs.theta,Inputs.Br./Inputs.B0)
shading interp
view(0,90)
xlabel('phi')
ylabel('theta')
title('B direction dot surface normal')

%% Plots
figure('Position', [10 10 1500 500])
tiledlayout(1,2)
nexttile
trisurf(dtri,Inputs.phi,Inputs.theta,Inputs.Vlayer)
shading interp
view(0,90)
xlabel('phi')
ylabel('theta')
title(['case ' prefix 'V_{RF}'], 'Interpreter', 'none')
colorbar
%set(gca,'ColorScale','log')
%saveas(gcf,'Enorm.png')

%figure
%trisurf(delaunay(Inputs.Z,Inputs.phi),Inputs.Z,Inputs.phi,Inputs.Density./mean(Inputs.Density) - 1)
%shading interp
%view(0,90)
%xlabel('z coordinate')
%ylabel('azimuthal angle')
%title(['Density/mean(Density) - 1. mean(Density)=',num2str(mean(Inputs.Density),3)])
%colorbar


%% Variable Creation

mu = 24.17; % for Deuterium; upgrade to arbitrary single ion species is in progress
e=1.602E-19; %J to eV, also elementary charge
epsilon_0org=8.8419E-12;  %F/m or C^2/J or s^2*C^2/m^3*kg %Permittivity of free space
epsilon_0=epsilon_0org*e; %C^2/eV %Permittivity of free space
m_D_kg=3.34E-27; %kg %mass of D ion
w=(476e6)*2*pi; %rad/s
%omega=0; %0 for unmagnetized case; omega_ci for magnetized case
%bx=1*sind(90);
%xi=14;
%bn=1; %1 for unmagnetized limit
bn=abs(Inputs.Br ./ Inputs.B0);
LayerThickness=Inputs.d; %m

Te=5; %eV
Ti=0; %0; %eV
Z=1; %D
A=2; %D

thick_debye=sqrt((epsilon_0*Te)./(Inputs.Density*e^2)); %m %Chen2015

omegapi=1.32E3*Z*sqrt(Inputs.Density*1e-6./A); %cgs units %rad/s
omega_hat=w./omegapi;
w=omega_hat;

%omegaci=9.58e3.*Z.*Inputs.B.*10e4./A; %rad/s
%omegaci=Z*e*Inputs.B./(m_D_kg); %rad/s
omegaci=7.63e6.*Inputs.B0.*6.28; %rad/s for D, Dolan2013 eq.5.33
omega_ci_hat=omegaci./omegapi;
omega=omega_ci_hat;

bx=bn;
xi=1*Inputs.Vlayer./Te;
if unmagnetized == 1
    omega = 0; 
    bx = 1; 
end

% w=0.2;
% omega=0.3;
% bx=0.4;
% xi=13;

j=0;
upar0=1.1;

%% phi0avg

% w=0.4;
% xi=6;

a1 = 3.70285;
a2 = 3.81991;
b1 = 1.13352;
b2 = 1.24171;
a3 = 2.*b2./pi;

c0=0.966463;
c1=0.141639;

gg=c0+c1.*tanh(w);
xi1=gg.*xi;
ff=((log(mu)+xi1.*a1+xi1.^2.*a2+xi1.^3.*a3)./(1+xi1.*b1+xi1.^2.*b2))-log(1-(j./upar0))+log(mu./24.17);

phi0avg=ff;


nexttile
trisurf(dtri,Inputs.phi,Inputs.theta,phi0avg*Te)
shading interp
view(0,90)
xlabel('phi')
ylabel('theta')
title(['case ' prefix 'V_{DC}'], 'Interpreter', 'none')
colorbar
%set(gca,'ColorScale','log')

writematrix([data(:,1),data(:,2),data(:,3),Inputs.Density,phi0avg*Te],[prefix '_VDC.csv']);

maxmeanmed = [max(phi0avg*Te),mean(phi0avg*Te),median(phi0avg*Te)]

%% ni1avg
% niw = ion density at the wall for a static sheath
% vv is the total dc potential drop 
%niw(omega,bx, xi)
%niww(w,omega,bx,xi)

% w=0.2;
% omega=0.3;
% bx=0.4;
% xi=13;

k0=3.7616962640756197;
k1=0.2220204461728174;

phiavg=phi0avg;
philowomega=k0+k1.*(xi-k0)-log(1-(j./upar0));
phimod=philowomega+(phiavg-philowomega).*tanh(w);


d0=0.7944430930529499;
d1=0.803531266389172;
d2=0.18237897510951012;
d3=0.9957212047604492;
nu1=1.4555923231100891;

arg=sqrt(((mu.^2.*bx.^2+1)./(mu.^2+1)));

fff=-log(arg)./(1+d3.*omega.^2);

Phip=phimod-fff;
Phip1=zeros(1);

for ii=1:length(Phip)
    Phip1(ii,1)=Phip(ii,1);
if Phip(ii,1) < 0
    Phip1(ii,1)=0;
end
end

omegaPhi=omega.*Phip1.^(1/4);
d4=(d2.^2)./((mu.^2.*d0.^2)-d2.^2);

niw=(d0./(d2+sqrt(Phip1))).*sqrt((bx.^2+d4+d1.^2.*omegaPhi.^(2.*nu1))./(1+d4+d1.^2.*omegaPhi.^(2.*nu1)));

% niwOmega is the ion density at the wall for an rf sheath
% Xi is the 0-peak rf voltage

niwomega=real(niw);

%% yd

s0=1.1241547327789232;
phi0a=phi0avg;
niwomegaa=niwomega;
Delta=sqrt(phi0a./niwomegaa);

yd=-sqrt(-1).*s0.*w./Delta;

%% ye

% bx=0.4;
% xi=3.6;

h1=0.607405123251634;
h2=0.3254965671158986;
g1=0.6243920388599393;
g2=0.5005946718280853;
g3=(pi./4).*h2;


he=(1+xi.*h1+xi.^2*h2)./(1+xi.*g1+xi.^2.*g2+xi.^3.*g3);

h0=1.05704235;

ye=h0.*abs(bx).*he.*(1-(j./upar0));

%% yi

parp0=1.0555369617763768;
parp1=0.7976591020008023;
parp2=1.47404874815277;
parp3=0.8096145628336325;

wcup=parp3.*w./sqrt(niwomegaa);
ycup=abs(bx)./(niwomegaa.*sqrt(phi0a));
epsilon=0.0001;

gsmall=(w.^2-bx.^2.*omega.^2+sqrt(-1).*epsilon)./(w.^2-omega.^2+sqrt(-1).*epsilon);
yi0=niwomegaa./sqrt(phi0a);

yi=parp0.*yi0.*((sqrt(-1).*wcup))./((wcup.^2./gsmall)-parp1+sqrt(-1).*parp2.*ycup.*wcup);

%% ytot and ztot

ytot=yi+ye+yd;

ztot=1./ytot;

%% Epsilon and Sigma Calculations

DataOut(:,1)=-(imag(ytot)./omega_hat).*(LayerThickness./thick_debye); %epsilon

DataOut(:,2)=epsilon_0org.*omegapi.*(LayerThickness./thick_debye).*real(ytot); %sigma

Epsilon(:,1)=Inputs.phi;
Epsilon(:,2)=Inputs.theta;
Epsilon(:,3)=DataOut(:,1);
Sigma(:,1)=Inputs.phi;
Sigma(:,2)=Inputs.theta;
Sigma(:,3)=DataOut(:,2);
if init == 1
    Epsilon(:,3)=DataOut(:,1) * 0 + 1;
    Sigma(:,3)=DataOut(:,2) * 0;
end

%% print txt files for COMOSL

csvwrite('sigma.csv', Sigma);
csvwrite('epsilon.csv', Epsilon);





end