close all
clearvars;
ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;
Tarray = [5,6,7,8,9,10]
narray = [1e16, 1e17, 5e17, 1e18, 5e18, 1e19]
Varray = [30, 160, 290, 420, 550, 680, 810, 940, 1070, 1200]


file = 'ftridyn_NeonW_80keV.nc';
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


A = readmatrix('Targets_Ne3+.txt','NumHeaderLines',1);
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

nLoc = size(A,1);
nCharge = 3;
yields = zeros(nLoc,nCharge);
surfCounts = zeros(nLoc,nCharge);
meanE = zeros(nLoc,nCharge);
meanA = zeros(nLoc,nCharge);

for i=1:nLoc
for j=3
filename = strcat('surface_C',string(j),'_loc_',string(i-1),'.nc');
if ~isfile(filename)
    continue
end
surfEDist = ncread(filename,'surfEDist');
surfCounts(i,j) = sum(surfEDist(:));

eff_yield = surfEDist.*Y0;
eff_yield(isnan(eff_yield)) = 0;

% figure
% h = pcolor(linspace(0,90,90),linspace(0,1000,1000),surfEDist')
% h.EdgeColor = 'none';

%             yield = mean(Y0);
            denom = sum(surfEDist(:));
            if denom > 0
                yields(i,j) = sum(eff_yield(:))./denom;
            else
                yields(i,j) = 0;
            end
            
%             meanE(i,j) = mean(E);
%             meanA(i,j) = mean(impact_angle);
        end
end


yields(isnan(yields))=0;


writematrix(yields,"yields_Ne3+.csv");
figure
plot(yields(:,3),'lineWidth',2)
% title({'ne'})
xlabel('location #') % x-axis label
ylabel('Yield') % y-axis label
title('IEADS Effective Yield (Ne3+, C3)')
set(gca,'fontsize',16)
grid on

if exist('exportgraphics','file') == 2
    exportgraphics(gcf,'ieads_yield_Ne3.png','Resolution',300);
else
    saveas(gcf,'ieads_yield_Ne3.png');
end

% Plot integrated surfEDist on the full geometry surface.
if (exist('x1','var') == 0)
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    fclose(fid);
end

subset = 1:length(x1);
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];

if numel(x1) == nLoc + 1
    surfGeomVals = [0; surfCounts(:,3)];
elseif numel(x1) == nLoc
    surfGeomVals = surfCounts(:,3);
else
    nGeom = numel(x1);
    surfGeomVals = zeros(nGeom,1);
    nCopy = min(nGeom,nLoc);
    surfGeomVals(1:nCopy) = surfCounts(1:nCopy,3);
end

figure
patch(transpose(X),transpose(Y),transpose(Z),surfGeomVals,'FaceAlpha',1,'EdgeAlpha',0.2)
title('Integrated surfEDist on GITR Geometry (Ne3+, C3)')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')
set(gca,'fontsize',14)
if exist('exportgraphics','file') == 2
    exportgraphics(gcf,'surfEDist_geometry_Ne3.png','Resolution',300);
else
    saveas(gcf,'surfEDist_geometry_Ne3.png');
end

save('ieads_ITER_Ne3.mat','yields','surfCounts','nLoc','nCharge','ne','te','ti','vp','btot','br','bt','bz','angle_imp')
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
h = pcolor(narray,Tarray,yields(:,:,2)')
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
