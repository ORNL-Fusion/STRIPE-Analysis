% Develop density profile for OXB process:
% Expressed in normalized coordinates:

clear all
close all

% Save flag:
saveFig = 1;

% Normalized "x" coordinate:
x = linspace(0,3,1e3);

% Radial gradient length scale:
L_star = 0.6;

% Over dense factor:
a = 3;

% Inflection point:
x_star = 1 - L_star*log(sqrt( a-1 ));

% Normalized density profile:
phi   = @(X) ( (X - x_star) / L_star );
n_hat = @(X,A) A*(1 - tanh(phi(X)) )/2;

% Ratio wce/wRF:
b = 0.75;

% Density at the UHR:
n_hat_UHR = 1 - b^2;

% Find radius at UHR:
i_UHR = find(n_hat(x,a) < n_hat_UHR,1);
r_n_hat_UHR = x(i_UHR);
r_hat_max   = 2.75;

% Plot figure:
% =========================================================================
figure('color','w')
set(gcf,'position'     ,[450 288 505 300])
% set(gcf,'position'     ,[647 328 505 233])


% plasma radial profile:
% -------------------------------------------------------------------------
subplot(1,6,[1,2])
hold on

% Overdense plasma region:
rng = find(x < 1);
hA1 = area(n_hat(x(rng),a),x(rng));
hA2 = area([0 1],[1 1]);
set(hA1,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none')
set(hA2,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none')

% Density profile:
f = n_hat(x,a);
plot(f,x,'k','LineWidth',2)

% O-mode cutoff:
plot(1       ,1          ,'ko','MarkerFaceColor','k')

% UH resonance:
plot(n_hat_UHR,r_n_hat_UHR,'ro','MarkerFaceColor','r')

% Axes position:
ax_pos = get(gca,'position');
ax_pos(2) = ax_pos(2)*1.3;
ax_pos = set(gca,'position',ax_pos);

% Formatting:
set(gca,'fontSize',10,'fontName','Times')
set(gca,'YTick',[0 1 2])
box on
grid off
xlim([0,a])
ylim([0,1/L_star])
ylim([0,r_hat_max])
xlabel('$\hat{n}(\hat{r})$','interpreter','latex','fontsize',14)
ylabel('$\hat{r}$','interpreter','latex','fontsize',15,'Rotation',0,'HorizontalAlignment','right')

% Plasma r-z profile:
% -------------------------------------------------------------------------
subplot(1,6,[3,6])
hold on

% Share over-dense plasma:
hA = area([0,1],[1,1]);
hA.FaceColor = [0.9 0.9 0.9];
hA.EdgeColor = 'none';

% Cut-off and resonance layers:
hLayers(1) = line([0,1]    ,[1,1]*r_n_hat_UHR); % UHR
hLayers(2) = line([0,1]    ,[1,1]            ); % O-mode
hLayers(3) = line([1,1]*0.2,[0,1]*r_hat_max  ); % w_ce

% Axes position:
ax_pos = get(gca,'position');
ax_pos(2) = ax_pos(2)*1.3;
ax_pos = set(gca,'position',ax_pos);

% Formatting:
set(hLayers(1),'color','r','LineWidth',3)
set(hLayers(2),'color','k','LineWidth',3)
set(hLayers(3),'color','k','LineWidth',2,'LineStyle','--')
% set(gca,'PlotBoxAspectRatio',[1 0.82 1])
ylim([0,r_hat_max])
set(gca,'YTickLabel',[])
set(gca,'fontSize',10,'fontName','Times')
box on
grid off
xlabel('$\hat{z}$','interpreter','latex','fontsize',15,'Rotation',0,'HorizontalAlignment','right')
hLeg = legend([hLayers,hA],'UH resonance','O-mode cutoff','Cyclotron resonance','Overdense plasma');
set(hLeg,'Location','NorthWest','interpreter','latex','fontSize',10)

% O-mode ray:
fields.Color = 'bl';
fields.LineWidth = 3;
hta = myTextArrow(gca,[1 0.6],[r_hat_max 1],fields);
hText = text(0.9,2,'O','FontWeight','bold','fontSize',14,'Color','bl');

% X-mode ray:
fields.Color = 'g';
fields.LineWidth = 3;
hta = myTextArrow(gca,[0.6 0.45],[1 r_n_hat_UHR],fields);
hText = text(0.46,r_n_hat_UHR*1.1,'X','FontWeight','bold','fontSize',14,'Color','g');

% B-mode ray:
fields.Color = 'r';
fields.LineWidth = 3;
hta = myTextArrow(gca,[0.45 0.2],[r_n_hat_UHR 0.3],fields);
hText = text(0.25,0.3,'B','FontWeight','bold','fontSize',14,'Color','r');

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_1_OXB_diagram'];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end
