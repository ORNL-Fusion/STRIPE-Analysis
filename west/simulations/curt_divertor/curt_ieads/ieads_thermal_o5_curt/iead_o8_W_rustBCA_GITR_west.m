close all
clear all
ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;
Tarray = [5,6,7,8,9,10]
narray = [1e16, 1e17, 5e17, 1e18, 5e18, 1e19]
Varray = [30, 160, 290, 420, 550, 680, 810, 940, 1070, 1200]


file = 'RustBCA_OonW.nc';
ncid = netcdf.open(file,'NC_NOWRITE');
[dimname, nE] = netcdf.inqDim(ncid,0);
[dimname, nA] = netcdf.inqDim(ncid,1);
if strcmp(file,'ftridynBackground.nc')
    [dimname, nS] = netcdf.inqDim(ncid,2);
else
    nS = 1;
end
energy = ncread(file,'E');
angle = ncread(file,'A');
spyld = ncread(file,'spyld');
rfyld = ncread(file,'rfyld');
% cosxDistRef = ncread(file,'cosXDistRef');
% cosyDist = ncread(file,'cosYDist');
% eDist = ncread(file,'energyDist');
% % eDistRef = ncread(file,'energyDistRef');
% eDistEgrid = ncread(file,'eDistEgrid');
% phiGrid = ncread(file,'phiGrid');
% thetaGrid = ncread(file,'thetaGrid');
% cosxDist = ncread(file,'cosXDist');

if( nS > 1)
    for i=1:nS
        this_spyld = spyld(:,:,i);
        figure(i)
        h = pcolor(energy,angle,log10(this_spyld))
        h.EdgeColor = 'none';
        colorbar
        set(gca, 'XScale', 'log')
    end
else
    this_spyld = spyld;
    figure(1)
    h = pcolor(energy,angle,this_spyld)
    h.EdgeColor = 'none';
    colorbar
    set(gca, 'XScale', 'log')
    set(gca, 'ColorScale', 'log')
    xlabel('Energy [eV]')
    ylabel('Angle [degrees]')
    
    title({'Sputtering Yield','O on W'})
 
end
energy_dist = linspace(0,1000,1000);
angle_dist = linspace(0,90,90);

[ee aa] = meshgrid(energy_dist,angle_dist);
Y0=interpn(energy,angle,spyld',ee,aa);

figure; plot(energy, this_spyld)


A = readmatrix('Targets.txt','NumHeaderLines',1);
% r_mean = readmatrix('Centroid.csv');
         ne =A(:,2);
         te = A(:,3);
         ti = A(:,4);
         vp = A(:,5);
         btot = A(:,6);
         br = A(:,7);
         bt = A(:,8);
         bz = A(:,9);

angle_imp = A(:,10);

nLoc = 1387;
nCharge = 5;
yields = zeros(nLoc,nCharge);
meanE = zeros(nLoc,nCharge);
meanA = zeros(nLoc,nCharge);

for i=1:nLoc
for j=nCharge
filename = strcat('surface_C',string(j),'_loc_',string(i-1),'.nc');
surfEDist = ncread(filename,'surfEDist');

eff_yield = surfEDist.*Y0;
eff_yield(isnan(eff_yield)) = 0;

% figure
% h = pcolor(linspace(0,90,90),linspace(0,1000,1000),surfEDist')
% h.EdgeColor = 'none';

eff_yield = sum(eff_yield(:))./sum(surfEDist(:));
%             yield = mean(Y0);
            yields(i,j) = eff_yield;
            
%             meanE(i,j) = mean(E);
%             meanA(i,j) = mean(impact_angle);
        end
end


yields(isnan(yields))=0;


writematrix(yields,"yields.csv");
figure
plot(yields(:,5),'lineWidth',2)
% title({'ne'})
xlabel('location #') % x-axis label
ylabel('Yield') % y-axis label
set(gca,'fontsize',16)
return

Dflux = vp.*ne.*tand(angle_imp);
Dflux(isnan(Dflux))=0;
writematrix(Dflux,'Dflux.csv')
figure
semilogy(vp.*ne,'lineWidth',2)
hold on
semilogy(Dflux,'lineWidth',2)
C4_flux= 0.05*Dflux;
semilogy(C4_flux,'lineWidth',2)
% title({'ne'})
xlabel('location #') % x-axis label
ylabel('Flux [m-2s-1]') % y-axis label
set(gca,'fontsize',16)
legend('Parallel Flux','W','O')
eroded_flux=yields(:,5).*C4_flux(1:end);
writematrix(eroded_flux,'eroded_flux.csv')
figure
plot(yields(:,5).*C4_flux(1:end),'lineWidth',2)
% title({'ne'})
xlabel('location #') % x-axis label
ylabel('Eroded W Flux [m-2s-1]') % y-axis label
set(gca,'fontsize',16)

figure(5)
h = pcolor(narray,Tarray,yields(:,:,8)')
h.EdgeColor = 'none';
colorbar
set(gca, 'XScale', 'log')
%     set(gca, 'ColorScale', 'log')
xlabel('Density [m^{-3}]')
ylabel('Temperature [eV]')

title({'Sputtering Yield','O on Al 940 Volts'})
set(gca,'TickDir','out');
set(gca,'fontsize',14)
ax = gca;

k = 0.02
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 10*k;

figure(6)
hold on
% for i=1:length(narray)
%     for j=1:length(Tarray)
        
        plot(Varray,reshape(meanE(3,3,:),1,[]),'Linewidth',2)
%     end
%     
% end
title({'Mean Impact Energy As Function of Surface Potential','n=10^{18} m^{-3} T = 7 eV'})
xlabel('Surface Potential [V]')
ylabel('Energy [eV]')
set(gca,'FontSize',14)
figure(7)
hold on
% for i=1:length(narray)
%     for j=1:length(Tarray)
        
        plot(Varray,reshape(meanA(3,3,:),1,[]),'Linewidth',2)
%     end
%     
% end
title({'Mean Impact Angle As Function of Surface Potential','n=10^{18} m^{-3} T = 7 eV'})
xlabel('Surface Potential [V]')
ylabel('Angle [degrees]')
set(gca,'FontSize',14)
% impact_angle = acosd(abs(vz)./sqrt(vx.^2 + vy.^2 + vz.^2));
% E = 0.5*amu*1.66e-27*(v(hit,1).^2+v(hit,2).^2+v(hit,3).^2)/Q;
% figure(120)
% hist3([90-impact_angle,E] ,'Edges',{0:90/180:90',0:1:2000'},'EdgeColor','none')
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% colorbar
% title({'Energy Angle Distribution of All Surfaces [Counts]','5 degrees'})
% ylabel('Energy [eV]')
% xlabel('Angle [degrees]')
% set(gca,'FontSize',10)
% axis([0 90 0 1000])
% view(2)
save('yields_ntv_OAl.mat','narray','Tarray','Varray','yields')