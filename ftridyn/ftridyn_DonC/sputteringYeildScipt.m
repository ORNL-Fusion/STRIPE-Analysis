file = 'ftridynSelf.nc';
ncid = netcdf.open(file,'NC_NOWRITE');
[dimname, nE] = netcdf.inqDim(ncid,0);
[dimname, nA] = netcdf.inqDim(ncid,1);
if strcmp(file,'assets/ftridynBackground.nc')
[dimname, nS] = netcdf.inqDim(ncid,2);
else
    nS = 1;
end
energy = ncread(file,'E');
angle = ncread(file,'A');
spyld = ncread(file,'spyld');
rfyld = ncread(file,'rfyld');
cosxDist = ncread(file,'cosXDist');
cosxDistRef = ncread(file,'cosXDistRef');
cosyDist = ncread(file,'cosYDist');
% coszDist = ncread(file,'cosZDist');
eDist = ncread(file,'energyDist');
eDistRef = ncread(file,'energyDistRef');
eDistEgrid = ncread(file,'eDistEgrid');
eDistEgridRef = ncread(file,'eDistEgridRef');
phiGrid = ncread(file,'phiGrid');
thetaGrid = ncread(file,'thetaGrid');
thisEdistRef = reshape(eDistRef(:,1,:),length(eDistEgridRef),[]);
% figure(100)
% plot(eDistEgridRef,thisEdistRef)

spyld=reshape(spyld(:,:,1),length(angle),length(energy));

figure(113)
h = pcolor(energy,angle,spyld)
h.EdgeColor = 'none';
colorbar
set(gca,'ColorScale','log')
% set(gca, 'YDir', 'normal')
 set(gca, 'XScale', 'log')
title({'Sputtering Yield D on Al','As a Function of Energy and Angle'})
xlabel('E [eV]') % x-axis label
ylabel('Angle [degrees]') % y-axis label
set(gca,'fontsize',16)