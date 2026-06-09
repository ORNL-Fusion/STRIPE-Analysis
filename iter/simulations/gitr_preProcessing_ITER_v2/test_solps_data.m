%% Example 5 – 2-D profiles for electrons, ions, and total Neon
run_path = '/Users/78k/ORNL Dropbox/Atul Kumar/123361/run';
Case = load_solps_case(run_path);

%% === Load geometry centroids ============================================
% Load ITER geometry (must contain 'centroid' and 'r_centroid')
load('final_ITERGeom.mat','centroid');
r_centroid = sqrt(centroid(:,1).^2 + centroid(:,2).^2);

% --- Find midplane (Z≈0) centroid at minimum R ---
zTol = 0.01;  % ±1 cm band around Z=0
mid_idx = abs(centroid(:,3)) <= zTol;

if any(mid_idx)
    mid_centroids = centroid(mid_idx,:);
    r_mid = r_centroid(mid_idx);
    z_mid = centroid(mid_idx,3);
    [R_min_mid, iMin] = min(r_mid);
    centroid_minR = mid_centroids(iMin,:);
    fprintf('Innermost midplane centroid: R = %.3f m, Z = %.3f m\n', ...
            R_min_mid, centroid_minR(3));
else
    warning('No centroids found near Z=0 within ±%.3f m', zTol);
    r_mid = []; z_mid = [];
    R_min_mid = NaN; centroid_minR = [NaN NaN NaN];
end

%% === Extract and compute SOLPS data =====================================
if Case.Geo.isUnstructured
    ne = Case.State.ne;                        % electron density [m^-3]
    Te = Case.State.te./1.602e-19;             % electron temp [eV]
    Ti = Case.State.ti./1.602e-19;
    Ua_e = Case.State.ua(:,2);                 % electron parallel flow
else
    ne = Case.State.ne;
    Te = Case.State.te./1.602e-19;
    Ti = Case.State.ti./1.602e-19;
    Ua_e = Case.State.ua(:,:,2);
end
flux_e = ne'.*abs(Ua_e);

%% --- D+ ion --------------------------------------------------------------
if Case.Geo.isUnstructured
    nD = Case.State.na(:,2);
else
    nD = Case.State.na(:,:,2);
end

%% --- total Neon (Ne0–Ne10+) ---------------------------------------------
idxNe = 6:16;
if Case.Geo.isUnstructured
    nNe_total = sum(Case.State.na(:,idxNe),2);
    flux_Ne_total = zeros(size(nNe_total));
    for k = idxNe
        flux_Ne_total = flux_Ne_total + Case.State.na(:,k).*abs(Case.State.ua(:,k));
    end
else
    nNe_total = sum(Case.State.na(:,:,idxNe),3);
    flux_Ne_total = zeros(size(nNe_total));
    for k = idxNe
        flux_Ne_total = flux_Ne_total + Case.State.na(:,:,k).*abs(Case.State.ua(:,:,k));
    end
end

%% === Individual 2-D plots ===============================================
plot_solps_2d_profile(Case,log10(ne),'log_{10}(n_e) [m^{-3}]'); 
colormap(jet); f_ne = gcf;

plot_solps_2d_profile(Case,log10(Te),'log_{10}(T_e) [eV]'); 
colormap(jet); f_Te = gcf;

plot_solps_2d_profile(Case,flux_e,'Electron flux n_e×|U_{||,e}| [m^{-2}s^{-1}]'); 
colormap(jet); f_fluxE = gcf;

plot_solps_2d_profile(Case,log10(nNe_total),'log_{10}(∑ n_{Ne^q}) [m^{-3}]'); 
colormap(jet); f_NeDen = gcf;

plot_solps_2d_profile(Case,log10(Ti),'log_{10}(T_i) [eV]'); 
colormap(jet); f_Ti = gcf;

plot_solps_2d_profile(Case,flux_Ne_total,'Total Ne flux ∑ n_{Ne^q}×|U_{||,Ne^q}| [m^{-2}s^{-1}]'); 
colormap(jet); f_fluxNe = gcf;

%% === Combine exactly those six plots (in the right order) ================
orderedFigs = [f_ne, f_Te, f_fluxE, f_NeDen, f_Ti, f_fluxNe];

figure('Color','w','Position',[100 100 1600 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for i = 1:6
    nexttile;
    axOld = findobj(orderedFigs(i),'Type','axes');
    if isempty(axOld), continue; end
    copyobj(allchild(axOld), gca);
    axis equal tight;
    box on;
    set(gca,'FontSize',12,'LineWidth',1);
    xlabel('R [m]','FontSize',11);
    ylabel('Z [m]','FontSize',11);
    colormap(jet);
    title(axOld.Title.String,'FontSize',12);
    
    % Add colorbar to each subplot
    cb = colorbar;
    cb.Box = 'on';
    cb.LineWidth = 0.8;
    cb.FontSize = 10;

    % --- Overlay centroids ---
    hold on;
    plot(r_centroid, centroid(:,3), '.', 'Color',[0.5 0.5 0.5], 'MarkerSize', 5);  % all centroids (gray)
    if ~isempty(r_mid)
        plot(r_mid, z_mid, 'ro', 'MarkerSize', 4, 'DisplayName','Midplane Centroids');
    end
    if ~isnan(R_min_mid)
        plot(R_min_mid, centroid_minR(3), 'yp', 'MarkerFaceColor','y', ...
             'MarkerSize', 10, 'DisplayName','Min-R Midplane Centroid');
    end
end

sgtitle('SOLPS 2-D Profiles: Electrons (Top Row) and Total Neon (Bottom Row)','FontSize',14);

% Optional: close the individual figure windows
% close(orderedFigs);