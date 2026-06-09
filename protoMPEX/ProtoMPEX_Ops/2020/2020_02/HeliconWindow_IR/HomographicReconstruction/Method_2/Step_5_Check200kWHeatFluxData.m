% Script to check content of heat flux data structure:

clear all
close all

%
% Load file and view comments
%
load HeatFlux2D_ExtrapolatedTo_200kW_created_26-May-2020.mat 
load HeatFlux2D_ExtrapolatedTo_200kW_created_24-Oct-2020.mat 
d.comment
d.La
d.Ra
d.RFpower

%
% Get heat flux
%
phi = d.phi_q2D;
flux = d.q2D;
z = d.z_q2D;

%
% Plot for confirmation
%
figure(1)
image(phi(1,:)*180/pi,z(:,1)*100,flux/1000,'CDataMapping','scaled')
colorbar
colormap(hot)
caxis([0 1000]);
set(gca,'XDir','reverse','YDir','normal')
xlabel('Angle [deg]')
ylabel('z [cm]')
title('Proto-MPEX 200kW Scaled Heat flux [kW/m^2]')

% 
% Save figure
%
figureName = ['CheckContentsOfHeatFluxData200kW'];
saveas(gcf,figureName,'tiffn')