clear all
close all

filename = 'ADAS_Rates_O.nc';

IonizationData.Temp = ncread(filename,'gridTemperature_Ionization');
IonizationData.Temp = 10.^IonizationData.Temp;
IonizationData.Density = ncread(filename,'gridDensity_Ionization');
IonizationData.Density = 10.^IonizationData.Density;
IonizationData.RateCoeff = ncread(filename,'IonizationRateCoeff');
IonizationData.RateCoeff = 10.^IonizationData.RateCoeff;

RecombinationData.Temp = ncread(filename,'gridTemperature_Recombination');

RecombinationData.Density = ncread(filename,'gridDensity_Recombination');
RecombinationData.Temp = 10.^RecombinationData.Temp;

RecombinationData.Density = 10.^RecombinationData.Density;
RecombinationData.RateCoeff = ncread(filename,'RecombinationRateCoeff');
RecombinationData.RateCoeff = 10.^RecombinationData.RateCoeff;
nT = 10;
T = 15;%linspace(0,100,nT);
n = 6e18;
tions = zeros(1,10);
trecs = zeros(1,10);

Zmax = 8;
S = zeros(nT,Zmax+1);
a = zeros(nT,Zmax+1);
b = zeros(Zmax+1,1);
b(end) = 1;
for i=1:Zmax+1
    Z = i-1;
    if Z< (Zmax-1)
        S(:,i) = 1e5*interpn(IonizationData.Density,IonizationData.Temp,IonizationData.RateCoeff(:,:,Z+1),n,T,'linear',0);
    end
    if Z>0
        a(:,i) = 1e5*interpn(RecombinationData.Density,RecombinationData.Temp,RecombinationData.RateCoeff(:,:,Z),n,T,'linear',0);
    end
end

A = zeros(Zmax+1,Zmax+1);
conc = zeros(Zmax+1,nT);

for i=1:nT
    
    A(1,1) = 1;%-S(i,1);
    A(1,2) = 0;%a(i,2);
    for j=2:Zmax
        A(j,j-1) = S(i,j-1);
        A(j,j) = -(S(i,j)+a(i,j));
        A(j,j+1) = a(i,j+1);
    end
    A(Zmax+1,end-1) = S(i,Zmax-1);
    A(Zmax+1,end) = -(S(i,Zmax)+a(i,Zmax));
%     A(Zmax+1,:) = ones(1,Zmax);
    
    b(1) = 1 ;%-a(i,Zmax);
    conc1 = A\b;
end
conc1 = abs(conc1);
conc1 = conc1/sum(conc1);

plot(0:Zmax,conc1,'--o','LineWidth',2)
title({'Equilibrium Charge State Distribution','of O in 15eV 6e18m^{-3} electron Plasma','dt=1e-6s nP=1e5 nT=1e4'})
xlabel('Charge State [#]') % x-axis label
ylabel('Distribution Fraction') % y-axis label
set(gca,'fontsize',16)
axis([0 Zmax 0 0.8])
% legend('neutral','1','2')
hold on
m = 184;
T = 20;
ti0 = 4;
vTh = sqrt(2*ti0*1.602e-19/m/1.66e-27);
n = 1e19;
tion = 1/(n*interpn(IonizationData.Density,IonizationData.Temp,IonizationData.RateCoeff(:,:,i),n,T,'linear',0));

mfp = vTh*tion
% file = '/Users/tyounkin/Docs/elder/d3d/tungsten/tests/values_1e5p_1e4t_1en6';
% fileID = fopen(file,'r');
% formatSpec = '%f';
% A = fscanf(fileID,formatSpec)
% hold on
% plot(0:1:19,A,':*','LineWidth',1)
% legend('Equilibrium Values','GITR')
% % x = ncread(file,'x');
% % y = ncread(file,'y');
% % z = ncread(file,'z');
% % vx = ncread(file,'vx');
% % vy = ncread(file,'vy');
% % vz = ncread(file,'vz');
% % charge = ncread(file,'charge');
% % weight = ncread(file,'weight');
% % sizeArray = size(x);
% % nP = sizeArray(2);
% % figure(2)
% % 
% % h1=histogram(charge)
% % vals=h1.Values;
% % vals = vals./sum(vals);
% % figure(1)
% % hold on
% % plot(h1.BinEdges(1:end-1)+0.5,vals,'-o')
% % legend('Equilibrium Values','GITR')
% 
% nP = 1e4;
% nT = 1e4;
% dt = 1e-6;
% charge = zeros(1,nP);
% Ss = zeros(1,nP);
% aa = zeros(1,nP);
% S = 1e-5*S;
% a = 1e-5*a;
% tic
% for i=1:nT
%     Ss = 0*Ss;
%     aa = 0*aa;
%     
% Ss = interp1(linspace(0,Zmax-1,Zmax),S,charge,'linear',0);
% tion = 1./(n*Ss);
% Pion = exp(-dt./tion);
% randIon = rand(1,nP);
% whereIonize = find(randIon <= (1-Pion));
% charge(whereIonize) = charge(whereIonize) + 1;
% 
% charge1 = find(charge > 0);
% aa(charge1) = interp1(linspace(1,Zmax,Zmax),a,charge(charge1)+1,'linear',0);
% 
% trec = 1./(n*aa);
% Prec = exp(-dt./trec);
% randRec = rand(1,nP);
% whereRec = find((randRec <= (1-Prec)) & (charge > 0));
% charge(whereRec) = charge(whereRec) - 1;
% 
% end
% toc
% figure(100)
% h1=histogram(charge)
% vals=h1.Values;
% vals = vals./sum(vals);
% 
% figure(1)
% hold on
% plot(h1.BinEdges(1:end-1)+0.5,vals,'-o')