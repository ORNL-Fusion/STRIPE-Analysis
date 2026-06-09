%% revised_script_physical_eV_angle_with_true_combined.m
% Energy grid: 0 -> 5000 bins corresponds to 0 -> 500 eV
% Angle grid: 0 -> 90 degrees

clc; clear; close all;

%% ----------------------------
% CONFIG
%% ----------------------------
tilt = 85;
file = sprintf('../tilted_targets/test/%d_degrees/MPEX_runs_v1/gitrm-surface.nc', tilt);

SKIP_FIRST_NSURFACES = 6;
DROP_LAST_NSURFACES  = 3;

TARGET_N_E_BINS = 5000;
COARSEN_ANGLE_FACTOR  = 1;
COARSEN_ENERGY_FACTOR = 5;

TARGET_ANGLE_DEG = 85;

scale = 1;

% Ta bins to combine explicitly
BIN_TA = 3:6;
ta_labels = {'Ta1^+','Ta2^+','Ta3^+','Ta4^+'};

% % Including neutrals
% BIN_TA = 2:6;
% ta_labels = {'Ta0+','Ta1^+','Ta2^+','Ta3^+','Ta4^+'};


CLIM_2D_INC  = [1e-5 1e1];
CLIM_2D_REFL = [1e-5 1e-2];

E_ZOOM_MIN = 0;
E_ZOOM_MAX = 500;

COLOR_LOG_2D_INC  = false;
COLOR_LOG_2D_REFL = true;

PLOT_REFLECTED = true;
NORMALIZE_REFLECTED = true;

%% ----------------------------
% READ DATA
%% ----------------------------
surfEDist_all    = ncread(file,'surfEDist');
surfReflDist_all = ncread(file,'surfReflDist');

[nA0, nE0, nSurf0, nSpecies0, nBins0] = size(surfEDist_all);
fprintf('surfEDist dims: A=%d, E=%d, S=%d, species=%d, bins=%d\n', ...
    nA0,nE0,nSurf0,nSpecies0,nBins0);

%% ----------------------------
% TRUE COMBINED INCIDENT IEAD
%% ----------------------------
tmp_comb = sum(surfEDist_all(:,:,:,:,BIN_TA),5);  % sum selected Ta bins
AES_comb = squeeze(sum(tmp_comb,4));              % sum over species dimension

[AE_comb, angle_deg, energy_eV] = prep_AES_to_AE( ...
    AES_comb, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
    TARGET_N_E_BINS, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale);

%% ----------------------------
% INDIVIDUAL Ta IEADs
%% ----------------------------
nTa = numel(BIN_TA);
inc_AE_list = cell(nTa,1);

for j = 1:nTa
    b = BIN_TA(j);

    tmp = surfEDist_all(:,:,:,:,b);
    AES = squeeze(sum(tmp,4));

    [AEj, angle_j, energy_j] = prep_AES_to_AE( ...
        AES, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
        TARGET_N_E_BINS, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale);

    inc_AE_list{j} = AEj;
end

nA = numel(angle_deg);
nE = numel(energy_eV);

AE_comb = AE_comb(1:nA,1:nE);

for j = 1:nTa
    inc_AE_list{j} = inc_AE_list{j}(1:nA,1:nE);
end

[Eidx, ~] = get_zoom_indices(energy_eV,E_ZOOM_MIN,E_ZOOM_MAX);
[~,angle_idx] = min(abs(angle_deg - TARGET_ANGLE_DEG));

%% ----------------------------
% REFLECTED PREP
%% ----------------------------
refl_AES = squeeze(sum(surfReflDist_all,4));

[refl_AE, ~, ~] = prep_AES_to_AE( ...
    refl_AES, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, ...
    TARGET_N_E_BINS, COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale);

refl_AE = refl_AE(1:nA,1:nE);

%% ----------------------------
% INCIDENT 2x3 PLOT
%% ----------------------------
figure('Color','w','Position',[80 80 1400 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
colormap(jet);

sgtitle(sprintf('Incident IEADs at MPEX target, tilt = %d^\\circ', tilt), ...
    'FontSize',16,'FontWeight','bold');

% Combined
nexttile(1)
pcolor(energy_eV(Eidx),angle_deg,AE_comb(:,Eidx)); shading flat
set(gca,'YDir','normal')
if COLOR_LOG_2D_INC
    set(gca,'ColorScale','log');
end
colorbar
xlabel('Energy [eV]')
ylabel('Incidence angle [deg]')
title('True Combined Ta IEAD')
xlim([0 500])
ylim([0 90])
hold on
plot([energy_eV(Eidx(1)),energy_eV(Eidx(end))], ...
     [angle_deg(angle_idx),angle_deg(angle_idx)],'w--','LineWidth',1.5)

% Ta1-Ta2
for j = 1:2
    nexttile(j+1)
    AE = inc_AE_list{j};
    pcolor(energy_eV(Eidx),angle_deg,AE(:,Eidx)); shading flat
    set(gca,'YDir','normal')
    if COLOR_LOG_2D_INC
        set(gca,'ColorScale','log');
    end
    caxis(CLIM_2D_INC)
    colorbar
    xlabel('Energy [eV]')
    ylabel('Incidence angle [deg]')
    title(ta_labels{j})
    xlim([0 500])
    ylim([0 90])
end

% 1D spectra
nexttile(4)
hold on; box on; grid on

spec_comb = nanmean(AE_comb(:,Eidx),1);
plot(energy_eV(Eidx),spec_comb,'k','LineWidth',2.5)

for j = 1:4
    AEj = inc_AE_list{j};
    specj = nanmean(AEj(:,Eidx),1);
    plot(energy_eV(Eidx),specj,'LineWidth',1.6)
end

xlabel('Energy [eV]')
ylabel('\langle IEAD \rangle_\theta')
title('Angle-averaged incident IEADs')
legend(['True Combined',ta_labels],'Location','northeast')
xlim([0 500])

% Ta3-Ta4
for j = 3:4
    nexttile(j+2)
    AE = inc_AE_list{j};
    pcolor(energy_eV(Eidx),angle_deg,AE(:,Eidx)); shading flat
    set(gca,'YDir','normal')
    if COLOR_LOG_2D_INC
        set(gca,'ColorScale','log');
    end
    caxis(CLIM_2D_INC)
    colorbar
    xlabel('Energy [eV]')
    ylabel('Incidence angle [deg]')
    title(ta_labels{j})
    xlim([0 500])
    ylim([0 90])
end

%% ----------------------------
% SHARPER SMOOTH TRUE COMBINED IEAD PLOT
%% ----------------------------
Z = AE_comb(:,Eidx);
E = energy_eV(Eidx);
A = angle_deg;

Zsmooth = Z;
Zsmooth(~isfinite(Zsmooth)) = 0;

smoothSigma = 0.5;
kernelSize  = 2;

[xk,yk] = meshgrid(-floor(kernelSize/2):floor(kernelSize/2), ...
                   -floor(kernelSize/2):floor(kernelSize/2));

G = exp(-(xk.^2 + yk.^2)/(2*smoothSigma^2));
G = G/sum(G(:));

Zsmooth = conv2(Zsmooth,G,'same');

Zplot = Zsmooth;
Zplot(Zplot <= 0) = NaN;

figure('Color','w','Position',[200 200 850 650])

imagesc(E,A,Zplot)
set(gca,'YDir','normal')

colormap(turbo)
colorbar

xlabel('Incident Energy [eV]','FontSize',16,'FontWeight','bold')
ylabel('Incidence Angle [deg from normal]','FontSize',16,'FontWeight','bold')
title(sprintf('True Combined Incident IEAD, tilt = %d^\\circ',tilt), ...
    'FontSize',18,'FontWeight','bold')

xlim([0 100])
ylim([0 90])
caxis([0 prctile(Zplot(isfinite(Zplot)),99.5)])

set(gca,'FontSize',14,'LineWidth',1.2)
box on

%% ----------------------------
% REFLECTED EAD
%% ----------------------------
if PLOT_REFLECTED

    refl2 = refl_AE(:,Eidx);

    if NORMALIZE_REFLECTED
        tot_ref = nansum(refl2(:));
        if tot_ref > 0
            refl2 = refl2 ./ tot_ref;
        end
    end

    figure('Color','w','Position',[200 200 900 650])

    pcolor(energy_eV(Eidx),angle_deg,refl2); shading flat
    set(gca,'YDir','normal')
    set(gca,'ColorScale','log')
    caxis(CLIM_2D_REFL)
    colorbar
    colormap(jet)

    xlabel('Energy [eV]')
    ylabel('Reflection angle [deg]')
    title(sprintf('Reflected EAD, species summed, tilt = %d^\\circ',tilt))

    xlim([0 500])
    ylim([0 90])
end

fprintf('Plotting complete.\n');

%% ========================================================================
% Local functions
%% ========================================================================

function [AEmat, angle_deg, energy_eV] = prep_AES_to_AE( ...
    AESmat, SKIP_FIRST_NSURFACES, DROP_LAST_NSURFACES, TARGET_N_E_BINS, ...
    COARSEN_ANGLE_FACTOR, COARSEN_ENERGY_FACTOR, scale)

    nSurf = size(AESmat,3);

    i0 = SKIP_FIRST_NSURFACES + 1;
    i1 = nSurf - DROP_LAST_NSURFACES;

    if i0 > i1
        error('Invalid surface selection.')
    end

    AESsel = AESmat(:,:,i0:i1);
    AE = squeeze(sum(AESsel,3));

    [nA,nE_orig] = size(AE);

    angle_grid0  = linspace(0,90,nA).';
    energy_grid0 = linspace(0,500,nE_orig);

    if ~isempty(TARGET_N_E_BINS) && TARGET_N_E_BINS ~= nE_orig

        energy_new = linspace(0,500,TARGET_N_E_BINS);
        AEi = nan(nA,numel(energy_new));

        for ia = 1:nA
            y = AE(ia,:);
            mask = ~isnan(y);

            if any(mask)
                y2 = y;
                y2(~mask) = 0;

                yi = interp1(energy_grid0,y2,energy_new,'linear',0);
                wi = interp1(energy_grid0,double(mask),energy_new,'linear',0);

                yi(wi < 0.5) = NaN;
                AEi(ia,:) = yi;
            end
        end

        AE = AEi;
        energy_grid0 = energy_new;
    end

    fa = max(1,round(COARSEN_ANGLE_FACTOR));
    fe = max(1,round(COARSEN_ENERGY_FACTOR));

    nA2 = floor(size(AE,1)/fa)*fa;
    nE2 = floor(size(AE,2)/fe)*fe;

    AE = AE(1:nA2,1:nE2);
    angle_grid0  = angle_grid0(1:nA2);
    energy_grid0 = energy_grid0(1:nE2);

    AE = reshape(AE,fa,nA2/fa,fe,nE2/fe);
    AEc = squeeze(mean(mean(AE,1,'omitnan'),3,'omitnan'));
    AEc = reshape(AEc,nA2/fa,nE2/fe);

    angle_deg  = block_mean_1d(angle_grid0,fa).';
    energy_eV  = block_mean_1d(energy_grid0,fe);

    angle_deg = angle_deg(:);

    AEc = scale .* AEc;

    % Keep weak values for diagnosis. Do not threshold aggressively here.
    AEc(AEc <= 0) = NaN;

    AEmat = AEc;
end

function y = block_mean_1d(x,f)
    x = x(:);
    n = numel(x);
    n2 = floor(n/f)*f;
    x = x(1:n2);
    X = reshape(x,f,n2/f);
    y = mean(X,1,'omitnan');
end

function [idx,xvec] = get_zoom_indices(xgrid,xmin,xmax)
    [~,i1] = min(abs(xgrid - xmin));
    [~,i2] = min(abs(xgrid - xmax));

    if i2 < i1
        tmp = i1;
        i1 = i2;
        i2 = tmp;
    end

    idx = i1:i2;
    xvec = xgrid(idx);
end