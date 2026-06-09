% Understanding the O-X mode conversion process:
clear all
close all

saveFig = 0;

% Define refractive index
N = 1;

% Define dimensionless parameters:
X = linspace(0,2,1e3);
Y = w_ce(0.40)/(2*pi*28E9);
% Y = w_ce(0.52)/(2*pi*28E9);

% Define incidence angle:
theta = 0*pi/180;
N_par(1) = N*cos(pi/2 - theta);
N_par(2) = 0.98*sqrt(Y/(1+Y));

% Derived quantities:
X_UHR = 1 - Y^2;

% Electron bernstein wave:
Te = 6; 
vTe = sqrt(2*e_c*Te/m_e);
b = vTe/c_light;

% FLR term:
L = 0.5*(linspace(1E-4,100,500).^2)*(b/Y)^2;

% Solve bessel function sum for every value of L:
N = 50;
n = -N:1:N;
f = zeros(size(L));

for ii = 1:numel(L)
    for nn = 1:numel(n)
        I_n = besseli(n(nn),L(ii));
        Num(nn) = (n(nn)^2)*I_n/L(ii);
        Den(nn) = 1 - (n(nn)*Y);
        S_n(ii,nn) = Num(nn)/Den(nn);
    end
    f(ii) = sum(S_n(ii,:)).*exp(-L(ii));
end

X_Bmode = 1./f;
n_perp_Bmode_a = sqrt(2*L*(Y/b)^2);

figure; 
plot(X_Bmode,n_perp_Bmode_a)
ylim([0,max(n_perp_Bmode_a)])
xlim([0,max(X)])


% Solve for L instead:
% clear n_perp_Bmode

L = 0.5*(linspace(1E-4,600,500).^2)*(b/Y)^2;

% Function:
N = 150;
n = -N:1:N;
S = @(LL,YY) (n.^2).*(besseli(n,LL)/LL)./(1 - (n*YY));
g = @(LL,XX,YY) 1 - XX*exp(-LL).*sum(S(LL,YY));

for ss = 1:numel(X)
    L_Bmode(ss) = fzero(@(LL) g(LL,X(ss),Y),1-Y^2);
    n_perp_Bmode(ss) = sqrt(2*L_Bmode(ss)*(Y/b)^2);
end

% Cold plasma waves:
for ii = 1:numel(N_par)
    
    % X-mode right hand cutoff based on Ram 2000:
    X_R = (1 - N_par(ii)^2)*(1-Y);
    X_L = (1 - N_par(ii)^2)*(1+Y);

    % Assemble Stix terms:
    S = 1 - X*(1./(1 - Y^2));
    P = 1 - X;
    R = 1 - X*(1./(1 - Y));
    L = 1 - X*(1./(1 + Y));

    % Assemble terms of biquadratic dispersion relation:
    Aperp = S;
    Bperp = -((R.*L) + (P.*S) - (N_par(ii).^2).*(P + S)); 
    Cperp = P.*( (N_par(ii).^2) - R).*( (N_par(ii).^2) - L);

    % Solve squared refractive index:
    n_perp_Omode = ( sqrt((-Bperp + sqrt(Bperp.^2 - 4*Aperp.*Cperp))./(2*Aperp)) );
    n_perp_Xmode = ( sqrt((-Bperp - sqrt(Bperp.^2 - 4*Aperp.*Cperp))./(2*Aperp)) );

    % Plot data:
    figure('color','w')
    box on
    hold on

    % O-mode:
    rng = find(imag(n_perp_Omode) < -0.1 & X > 1);
    n_perp_Omode(rng) = -1;
    hN(1) = plot(X,n_perp_Omode,'r','lineWidth',3);

    % X-mode:
    rng = find(imag(n_perp_Xmode) > 0.1 & X > 1);
    n_perp_Xmode(rng) = -1;
    hN(2) = plot(X,n_perp_Xmode,'k','lineWidth',3);

    % B-mode:
    hN(7) = plot(X,n_perp_Bmode,'bl','lineWidth',3);
    
    % O-mode cutoff:
    hN(3) = line([1,1],[1E-2,100]);
    set(hN(3),'color','r','lineStyle','--','lineWidth',2)

    % Right hand cutoff
    hN(4) = line([1,1]*X_R,[1E-2,100]);
    set(hN(4),'color','k','lineStyle','--','lineWidth',2)

    % Left hand cutoff:
    hN(5) = line([1,1]*X_L,[1E-2,100]);
    set(hN(5),'color','k','lineStyle','--','lineWidth',2)

    % UH resonance layer:
    hN(6) = line([1,1]*X_UHR,[1E-2,100]);
    set(hN(6),'color','g','lineStyle','--','lineWidth',2)

    set(gca,'YScale','log')
    set(gca,'fontName','times','fontSize',12)
    xlabel('$ \omega_{pe}/\omega$','interpreter','latex','fontSize',14)
    ylabel('$ n_{\perp} $','interpreter','latex','fontSize',16)
    grid on 
    %hL = legend(hN,'O-mode','X-mode','O-mode cutoff','X-mode R cutoff','X-mode L cutoff','UH resonance','B-mode');
    %hL.Location = 'northwest';
    ylim([0.01,100])
    xlim([0,max(X)])
end

%%
% Labels:
% fontSize = 16;
fontSize.label = 16;
fontSize.title = 15;

fontWeight = 'normal';
figure(2)
title(['Incidence angle = ',num2str(0),'$^{\circ}$, ','$ \omega_{ce}/\omega$ = ',num2str(Y,2)],'fontSize',fontSize.title,'interpreter','latex');
hT{1} = text(0.15,1.2,'O-mode','fontSize',fontSize.label,'interpreter','latex','color','r','FontWeight',fontWeight);
hT{2} = text(0.1,0.3,'FX-mode','fontSize',fontSize.label,'interpreter','latex','color','k','FontWeight',fontWeight);
hT{3} = text(0.42,0.012,'FX-cutoff','fontSize',fontSize.label,'interpreter','latex','color','k','rotation',90,'FontWeight',fontWeight);
hT{4} = text(0.67,0.012,'UHR','fontSize',fontSize.label,'interpreter','latex','color','g','rotation',90,'FontWeight',fontWeight);
hT{5} = text(1.078,0.012,'O-cutoff','fontSize',fontSize.label,'interpreter','latex','color','r','rotation',90,'FontWeight',fontWeight);
hT{6} = text(1.597,0.012,'SX-cutoff','fontSize',fontSize.label,'interpreter','latex','color','k','rotation',90,'FontWeight',fontWeight);
hT{7} = text(1.08,1.09,'SX-mode','fontSize',fontSize.label,'interpreter','latex','color','k','rotation',0,'FontWeight',fontWeight);

if Y < 0.5
    hT{8} = text(0.19,40,'B-mode','fontSize',fontSize.label,'interpreter','latex','color','bl','rotation',0,'FontWeight',fontWeight);
else
    hT{8} = text(1.08,19,'B-mode','fontSize',fontSize.label,'interpreter','latex','color','bl','rotation',0,'FontWeight',fontWeight);
end

% Save figure:
% =========================================================================
if saveFig
    Y_str = num2str(Y,2);
    Y_str(strfind(Y_str,'.')) = 'p'
    figureName = ['Step_11a_OX_dispersion_0deg_Y_',Y_str];
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

figure(3)
title(['Incidence angle = ',num2str( asind(N_par(2)),2),'$^{\circ}$, ','$ \omega_{ce}/\omega$ = ',num2str(Y,2)],'fontSize',fontSize.title,'interpreter','latex');
hT{1} = text(0.4,0.8,'O-mode','fontSize',fontSize.label,'interpreter','latex','color','r','FontWeight',fontWeight,'rotation',-8);
hT{2} = text(0.1,0.3,'FX','fontSize',fontSize.label,'interpreter','latex','color','k','FontWeight',fontWeight);
hT{3} = text(0.26,0.012,'FX-cutoff','fontSize',fontSize.label,'interpreter','latex','color','k','rotation',90,'FontWeight',fontWeight);
hT{4} = text(0.67,0.012,'UHR','fontSize',fontSize.label,'interpreter','latex','color','g','rotation',90,'FontWeight',fontWeight);
hT{5} = text(0.95,1.7,'O-cutoff','fontSize',fontSize.label,'interpreter','latex','color','r','rotation',90,'FontWeight',fontWeight);
hT{6} = text(1.08,1.5,'SX-cutoff','fontSize',fontSize.label,'interpreter','latex','color','k','rotation',90,'FontWeight',fontWeight);
hT{7} = text(1.15,0.4,'SX-mode','fontSize',fontSize.label,'interpreter','latex','color','k','rotation',0,'FontWeight',fontWeight);

if Y < 0.5
    hT{8} = text(0.19,40,'B-mode','fontSize',fontSize.label,'interpreter','latex','color','bl','rotation',0,'FontWeight',fontWeight);
else
    hT{8} = text(1.2,22,'B-mode','fontSize',fontSize.label,'interpreter','latex','color','bl','rotation',0,'FontWeight',fontWeight);
end

% Save figure:
% =========================================================================
if saveFig
    Y_str = num2str(Y,2);
    Y_str(strfind(Y_str,'.')) = 'p'
    figureName = ['Step_11a_OX_dispersion_30deg_Y_',Y_str];
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

