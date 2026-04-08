%% chargeScaled_areaIntegrated_flux_full.m
% Full script: charge-scale yields, compute area averages and
% area-integrated erosion flux per charge state, normalize to target total.
clc; clear; close all; format long;

%% ========== USER SETTINGS ==========
Gamma_int_target = 5.89e19;     % target total integrated erosion flux [particles/s]
save_csv = true;
save_fig = true;

%% ========== FILE LISTS ==========
labels = {'D+','Ne1+','Ne2+','Ne3+','Ne4+','Ne5+','Ne6+','Ne7+','Ne8+','Ne9+','Ne10+'};

yield_files = { ...
    "../ieads_D+/yields_D+.csv", ...
    "../ieads_Ne1+/yields_Ne1+.csv", ...
    "../ieads_Ne2+/yields_Ne2+.csv", ...
    "../ieads_Ne3+/yields_Ne3+.csv", ...
    "../ieads_Ne4+/yields_Ne4+.csv", ...
    "../ieads_Ne5+/yields_Ne5+.csv", ...
    "../ieads_Ne6+/yields_Ne6+.csv", ...
    "../ieads_Ne7+/yields_Ne7+.csv", ...
    "../ieads_Ne8+/yields_Ne8+.csv", ...
    "../ieads_Ne9+/yields_Ne9+.csv", ...
    "../ieads_Ne10+/yields_Ne10+.csv" };

target_files = { ...
    "../ieads_D+/Targets_D+.txt", ...
    "../ieads_Ne1+/Targets_Ne1+.txt", ...
    "../ieads_Ne2+/Targets_Ne2+.txt", ...
    "../ieads_Ne3+/Targets_Ne3+.txt", ...
    "../ieads_Ne4+/Targets_Ne4+.txt", ...
    "../ieads_Ne5+/Targets_Ne5+.txt", ...
    "../ieads_Ne6+/Targets_Ne6+.txt", ...
    "../ieads_Ne7+/Targets_Ne7+.txt", ...
    "../ieads_Ne8+/Targets_Ne8+.txt", ...
    "../ieads_Ne9+/Targets_Ne9+.txt", ...
    "../ieads_Ne10+/Targets_Ne10+.txt" };

nspec = numel(labels);

% charge numbers: D+ -> 1, Ne1+..Ne10+ -> 1..10
charge = zeros(nspec,1);
charge(1) = 1;
for s = 2:nspec, charge(s) = s-1; end

%% ========== GEOMETRY ==========
if ~exist('x1','var')
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    if fid == -1, error('Failed to open gitrGeometryPointPlane3d.cfg'); end
    % skip two header lines
    fgetl(fid); fgetl(fid);
    for ii = 1:18
        tline = fgetl(fid);
        if ischar(tline), evalc(tline); end
    end
    fclose(fid);
end

subset = 1:length(x1);
Xtri = [x1(subset).', x2(subset).', x3(subset).'];
Ytri = [y1(subset).', y2(subset).', y3(subset).'];
Ztri = [z1(subset).', z2(subset).', z3(subset).'];
Nfaces = size(Xtri,1);

% face areas
p1 = [Xtri(:,1), Ytri(:,1), Ztri(:,1)];
p2 = [Xtri(:,2), Ytri(:,2), Ztri(:,2)];
p3 = [Xtri(:,3), Ytri(:,3), Ztri(:,3)];
cp = cross(p2-p1, p3-p1, 2);
faceArea = 0.5 * sqrt(sum(cp.^2,2));
faceArea(faceArea==0) = eps;

Atot = sum(faceArea,'omitnan');
if ~isfinite(Atot) || Atot <= 0, error('Total area non-positive. Check geometry.'); end

fprintf('Nfaces = %d, Atot = %.6e m^2\n', Nfaces, Atot);

%% ========== HELPER: robust triangular yield loader ==========
function y = loadYieldTriangular(fname, Nfaces, col)
    M = readmatrix(fname);
    if isempty(M), error('Empty yield file: %s', fname); end
    [nr,nc] = size(M);
    if col > nc
        error('Yield file %s has %d columns; requested column %d', fname, nc, col);
    end
    yraw = M(:,col);
    yraw(~isfinite(yraw)) = 0;
    if numel(yraw) == Nfaces
        y = yraw;
    elseif numel(yraw) == Nfaces-1
        % common case: files missing leading zero -> prepend zero
        y = [0; yraw];
    elseif numel(yraw) > Nfaces
        y = yraw(1:Nfaces);
    else
        % shorter than Nfaces: pad zeros at end
        y = zeros(Nfaces,1);
        y(1:numel(yraw)) = yraw;
        warning('Yield %s shorter than Nfaces (%d). Padded with zeros.', fname, numel(yraw));
    end
    y = y(:);
end

%% ========== LOAD YIELDS and TARGET (n, v) ==========
Y = zeros(Nfaces, nspec);   % per-face yield
n = zeros(Nfaces, nspec);   % per-face density
v = zeros(Nfaces, nspec);   % per-face velocity (signed)

for s = 1:nspec
    % choose yield column: D -> col1, Neq+ -> col = q
    if s == 1
        col = 1;
    else
        col = charge(s);
    end
    Y(:,s) = loadYieldTriangular(yield_files{s}, Nfaces, col);

    % load target file
    D = readmatrix(target_files{s});
    if size(D,1) < Nfaces
        error('Target file %s has %d rows < Nfaces=%d', target_files{s}, size(D,1), Nfaces);
    end
    if size(D,2) < 5, error('Target file %s has <5 cols (need col5 for v).', target_files{s}); end

    v(:,s) = D(1:Nfaces,5);
    if s == 1
        if size(D,2) < 2, error('D target %s missing density col2.', target_files{s}); end
        n(:,s) = D(1:Nfaces,2);
    else
        if size(D,2) < 11, error('Ne target %s missing density col11.', target_files{s}); end
        n(:,s) = D(1:Nfaces,11);
    end

    % clean NaNs -> zeros
    Y(~isfinite(Y(:,s)),s) = 0;
    n(~isfinite(n(:,s)),s) = 0;
    v(~isfinite(v(:,s)),s) = 0;

    % debug summary
    fprintf('%-5s: yield col=%2d  min=%g  max=%g  mean=%g\n', labels{s}, col, min(Y(:,s)), max(Y(:,s)), mean(Y(:,s)));
end

%% ========== COMPUTE ION FLUX & SCALED EROSION FLUX ==========
Gamma_ion = n .* abs(v);               % per-face ion flux n|v| [#/m^2/s]

% scale yield by charge (use repmat to be robust)
Y_scaled = Y .* repmat(charge.', Nfaces, 1);   % q * Y_{q,f}
Gamma_ero_scaled = Y_scaled .* Gamma_ion;      % per-face erosion flux (qY n|v|) [#/m^2/s]

%% ========== AREA-AVERAGED QUANTITIES (total area) ==========
Yscaled_Aavg    = zeros(nspec,1);    % <qY>_A (dimensionless)
Geroscaled_Aavg = zeros(nspec,1);    % <qY n|v|>_A [#/m^2/s]

for s = 1:nspec
    valid = isfinite(faceArea); % faceArea should be finite
    Yscaled_Aavg(s)    = sum( Y_scaled(valid,s) .* faceArea(valid), 'omitnan') / Atot;
    Geroscaled_Aavg(s)= sum( Gamma_ero_scaled(valid,s) .* faceArea(valid), 'omitnan') / Atot;
end

%% ========== AREA-INTEGRATED FLUX PER CHARGE (raw) ==========
Phi_raw = zeros(nspec,1);   % [particles/s]
for s = 1:nspec
    Phi_raw(s) = sum( Gamma_ero_scaled(:,s) .* faceArea , 'omitnan' );
end
Phi_total_raw = sum(Phi_raw,'omitnan');

fprintf('\nRaw total integrated erosion flux (sum Phi_raw) = %.6e s^-1\n', Phi_total_raw);

%% ========== NORMALIZE TO TARGET TOTAL ==========
if Phi_total_raw <= 0
    error('Phi_total_raw <= 0, cannot normalize to Gamma_int_target.');
end

scaleFactor = Gamma_int_target / Phi_total_raw;
Phi_norm = Phi_raw * scaleFactor;   % normalized per-charge integrated flux [particles/s]
fprintf('Normalization factor = %.6e\n', scaleFactor);
fprintf('Normalized total (sum Phi_norm) = %.6e s^-1 (target %.6e)\n', sum(Phi_norm), Gamma_int_target);

%% ========== PRINT & SAVE TABLE ==========
% Build table with useful columns
T = table(labels(:), charge(:), Yscaled_Aavg, Geroscaled_Aavg, Phi_raw, Phi_norm, ...
    'VariableNames', {'Species','q','Yscaled_areaAvg_totalA','EroFlux_areaAvg_totalA_per_m2_s','Phi_raw_per_s','Phi_norm_per_s'});

fprintf('\n=== Results (first rows) ===\n');
disp(T);

if save_csv
    writetable(T, 'chargeScaled_areaIntegrated_results.csv');
    fprintf('Wrote chargeScaled_areaIntegrated_results.csv\n');
end

%% ========== PLOTS ==========
% 1) Area-averaged yields and area-averaged erosion flux (linear)
figure('Color','w','Name','Area-averaged yield & erosion flux','Position',[100 100 1100 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
bar(Yscaled_Aavg,'FaceColor',[0.2 0.6 0.9]);
set(gca,'XTick',1:nspec,'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('<qY>_A'); title('Area-averaged yield (qY), total-area averaging'); grid on;

nexttile;
bar(Geroscaled_Aavg,'FaceColor',[0.9 0.4 0.2]);
set(gca,'XTick',1:nspec,'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('<qY n|v|>_A [#/m^2/s]'); title('Area-averaged erosion flux (total-area averaging)'); grid on;

% 2) Area-integrated (normalized) erosion flux per charge [particles/s]
figure('Color','w','Name','Normalized integrated erosion flux','Position',[120 120 900 420]);
bar(Phi_norm,'FaceColor',[0.3 0.5 0.8]);
set(gca,'XTick',1:nspec,'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('Area-integrated erosion flux [particles/s]');
title(sprintf('Phi_q (normalized)  — total = %.3e s^{-1}', Gamma_int_target));
grid on;

% annotate percentages over bars
for i = 1:nspec
    pct = 100 * Phi_norm(i) / Gamma_int_target;
    text(i, Phi_norm(i), sprintf('%.2f%%', pct), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end

if save_fig
    saveas(gcf,'Phi_norm_bar.png');
    fprintf('Saved figure Phi_norm_bar.png\n');
end

%% ========== OPTIONAL: per-charge face histograms (if desired) ==========
% If you want per-charge face-distributions (histograms of per-face qY or per-face qY n|v|),
% uncomment the block below. For large Nfaces this can produce many figures, so it's off by default.
do_face_histograms = false;
if do_face_histograms
    for s = 1:nspec
        figure; 
        subplot(1,2,1);
        hist(Y_scaled(:,s), 100);
        title([labels{s} ' face histogram: qY']);
        xlabel('qY'); ylabel('count');

        subplot(1,2,2);
        hist(Gamma_ero_scaled(:,s), 100);
        title([labels{s} ' face histogram: qY n|v|']);
        xlabel('qY n|v| [#/m^2/s]'); ylabel('count');
    end
end

fprintf('\nDone.\n');