function densplot()

filename = ["Eout_sheath1.txt" "Eout_sheath2.txt" "Eout_sheath3.txt" "Eout_sheath4.txt"];
fileNormalsB = ["new_sheath1_coords.txt" "new_sheath2_coords.txt" "new_sheath3_coords.txt" "new_sheath4_coords.txt"];
fileNe = 'plasma_profiles_lowdens_equilBaseline2010.txt';
init = 0;
unmagnetized = 0;

NeTeData = readmatrix(fileNe,'numHeaderLines',1);
L=size(NeTeData,1)
p=polyfit(NeTeData(L-20:L,1),log(NeTeData(L-20:L,3)),1);
newx=linspace(max(NeTeData(:,1)),8.5,10);
extrapolatedDensityX=[NeTeData(1:L-1,1);newx'];
extrapolatedDensityN=[NeTeData(1:L-1,3);(exp(p(2)+p(1)*newx))'];


LayerData=[];
normalsBData=[];
coord2DTor=[];
coord2DPol=[];
dataRanges=[];
simpleCoord1=[];
simpleCoord2=[];

for i=1:length(filename)
    ldi=readmatrix(filename(i),'numHeaderLines',9);
    ndi=readmatrix(fileNormalsB(i),'numHeaderLines',9);
    dataRanges=[dataRanges;[size(LayerData,1)+1,size(LayerData,1)+size(ldi,1)]];

    if i==1 || i==3
        coordInterp = readmatrix('C:\Users\wti\Documents\ITERCOMSOL\s1CoordinatesDef.csv');
        vecv=[0.119286, 0., -0.99286];
        vecv1=[-0.99286, 0., -0.119286];
        vecv2=[0,-1,0];
        if i==3
            vecv2=-vecv2;
        end
    else
        coordInterp = readmatrix('C:\Users\wti\Documents\ITERCOMSOL\s2CoordinatesDef.csv');
        vecv=[0.0847387, 0., 0.996403];
        vecv1=[0.996403, 0., -0.0847387];
        vecv2=[0,-1,0];
        if i==4
            vecv2=-vecv2;
        end
    end
    proj2dX = ldi(:,1:3)*vecv1';
    proj2dY = ldi(:,1:3)*vecv2';

    coord2DTor = [coord2DTor;griddata(coordInterp(:,1),coordInterp(:,2),coordInterp(:,3),proj2dX,proj2dY)];
    coord2DPol = [coord2DPol;ldi(:,1:3)*vecv'];
    simpleCoord1=[simpleCoord1;proj2dX];
    simpleCoord2=[simpleCoord2;ldi(:,1:3)*vecv'];

    LayerData=[LayerData; ldi];
    normalsBData=[normalsBData; ndi];
end



Inputs.X= LayerData(:,1);
Inputs.d= 0.002; %[m] assuming data is at middle of layer
Inputs.Y = LayerData(:,2);
%Inputs.phi = atan2(Inputs.Y, Inputs.X);
Inputs.Z = LayerData(:,3);
%Inputs.phi_uniq = unique(Inputs.phi);
%Inputs.Z_uniq = unique(Inputs.Z);
Inputs.E_norm = (LayerData(:,4) + 1i*LayerData(:,5)) .* normalsBData(:,4) ...
              + (LayerData(:,6) + 1i*LayerData(:,7)) .* normalsBData(:,5) ...
              + (LayerData(:,8) + 1i*LayerData(:,9)) .* normalsBData(:,6); % normal to sheath layer
Inputs.Vlayer = abs(Inputs.E_norm) * Inputs.d; %2 is peak-to-peak voltage
Inputs.Density=interp1(extrapolatedDensityX,extrapolatedDensityN,normalsBData(:,10));
Inputs.Br =     normalsBData(:,7) .* normalsBData(:,4) ...
              + normalsBData(:,8) .* normalsBData(:,5) ...
              + normalsBData(:,9) .* normalsBData(:,6); % normal to sheath layer;
Inputs.B0 =  vecnorm(normalsBData(:,7:9),2,2);


%% Plot normal E
writematrix([LayerData(:,1:3),abs(Inputs.E_norm)],'Enormal4.txt')
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,abs(Inputs.E_norm),'abs normal component of RF E')
saveas(gcf,'Enormal4.png')

%% Plot VRF
writematrix([LayerData(:,1:3),abs(Inputs.E_norm)*Inputs.d],'VRF4.txt')
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,abs(Inputs.E_norm)*Inputs.d,'VRF')
saveas(gcf,'VRF4.png')

%% Plot B0
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,abs(Inputs.B0),'Norm of the confining magnetic field')
saveas(gcf,'absB4.png')

%% Plot B_normal
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,Inputs.Br,'Normal component of the confining magnetic field')
saveas(gcf,'Br4.png')

%% Plot density
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,log10(Inputs.Density),'Log10(Density)')
saveas(gcf,'ne4.png')
%return

%% Variable Creation

mu = 24.17; % for Deuterium; upgrade to arbitrary single ion species is in progress
e=1.602E-19; %J to eV, also elementary charge
epsilon_0org=8.8419E-12;  %F/m or C^2/J or s^2*C^2/m^3*kg %Permittivity of free space
epsilon_0=epsilon_0org*e; %C^2/eV %Permittivity of free space
m_D_kg=3.34E-27; %kg %mass of D ion
w=13.56e6*2*pi; %rad/s
%omega=0; %0 for unmagnetized case; omega_ci for magnetized case
%bx=1*sind(90);
%xi=14;
%bn=1; %1 for unmagnetized limit
bn=abs(Inputs.Br ./ Inputs.B0);
LayerThickness=Inputs.d; %m

Te=7.09857; %eV
Ti=18.58276; %0; %eV
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

%% phi0avg -> Myra says this is DC voltage

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

%% plot VDC
writematrix([LayerData(:,1:3),phi0avg],'VDC4.txt')
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,phi0avg,'VDC')
saveas(gcf,'VDC4.png')

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

Epsilon(:,1)=simpleCoord1;
Epsilon(:,2)=simpleCoord2;
Epsilon(:,3)=DataOut(:,1);
Sigma(:,1)=simpleCoord1;
Sigma(:,2)=simpleCoord2;
Sigma(:,3)=DataOut(:,2);
if init == 1
    Epsilon(:,3)=DataOut(:,1) * 0 + 1;
    Sigma(:,3)=DataOut(:,2) * 0;
end

%% Plot sigma
writematrix([LayerData(:,1:3),Sigma(:,3)],'sigma4.txt')
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,Sigma(:,3),'Computed sheath-equivalent \sigma')
saveas(gcf,'sigma4.png')

%% Plot epsilon
writematrix([LayerData(:,1:3),Epsilon(:,3)],'epsilon4.txt')
plot4(dataRanges,coord2DTor,coord2DPol,LayerData,Epsilon(:,3),'Computed sheath-equivalent \epsilon')
saveas(gcf,'epsilon4.png')

%% print txt files for COMOSL

for i=1:size(dataRanges,1)
    csvwrite(['sigma' num2str(i) '.csv'], Sigma(dataRanges(i,1):dataRanges(i,2),:));
    csvwrite(['epsilon' num2str(i) '.csv'], Epsilon(dataRanges(i,1):dataRanges(i,2),:));
end

end