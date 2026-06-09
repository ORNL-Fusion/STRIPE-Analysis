% Planck's blackbody radiation law:
clear all
close all

% Physical constants:
% =========================================================================
c = c_light;
h = h_planck;

% Spectral range:
% =========================================================================
lm = linspace(0,30e-6,1e3);

% Black-body spectral radiance:
% =========================================================================
L_lm = @(T) ((2*pi*h*c^2)./(lm.^5))./( (exp(h*c./(lm.*k_B*T))) - 1 );

% FLIR IR camera spectral range:
% =========================================================================
lmIR = linspace(7,14)*1e-6;

% Plot data:
% =========================================================================
figureName = 'Step_4_SpectralRadiance';

figure('color','w')
set(gca,'FontName','times','Fontsize',12)
hold on
area(lmIR*1e6,ones(size(lmIR))*5e8,'FaceColor','g','LineStyle','none')
hL(1) = plot(lm*1e6,L_lm(300),'k','LineWidth',3);
hL(2) = plot(lm*1e6,L_lm(400),'r','LineWidth',3);
hL(3) = plot(lm*1e6,L_lm(500),'bl','LineWidth',3);
box on
hLeg = legend(hL,'300 [K]','400 [K]','600 [K]')
grid on
ylabel('$L_\lambda^{bb}$ [$\mathrm{W m^{-2} sr^{-1} m^{-1}}$]'...
    ,'interpreter','Latex','Fontsize',13)
title('Spectral radiance $L_\lambda^{bb}$'...
    ,'interpreter','Latex','Fontsize',13)
xlabel('$\lambda$ [m$^{-1}$]','interpreter','Latex','Fontsize',13)
set(gcf,'Position',[423.0000  263.6667  497.0000  354.3333]);

% Save figure:
% =========================================================================
saveas(gcf,figureName,'tiffn')