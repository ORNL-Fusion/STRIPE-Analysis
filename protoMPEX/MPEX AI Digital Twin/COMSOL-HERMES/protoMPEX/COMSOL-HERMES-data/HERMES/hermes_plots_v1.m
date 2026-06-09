%% hermes_plots_v1
% Plot HERMES plasma profiles for all cases and iterations.
%
% Outputs (written to plots/<case>/<iteration>/):
%   overview_rz.png           — Ne and Te in the theta=0 R-Z plane
%   profiles_1d.png           — radial and axial 1D profiles at theta=0
%   azimuthal_overview.png    — Ne / Te / q_Bohm in azimuthal planes across device
%   helicon_window_thetaz.png — q_Bohm / Ne / Te as theta-z map at the wall
%   <case>_<iter>_xyz_ne_te.csv

clear; clc; close all;

%% ---------- user settings ----------

zWindowMin       = 1.6;    % helicon window start [m]
zWindowMax       = 1.9;    % helicon window end   [m]
nAzimuthalSlices = 6;      % evenly-spaced z-slices across device
m_ion            = 3.344e-27;  % ion mass [kg]: 1.672e-27 H, 3.344e-27 D

windowLabels = {'it0'; 'it1'; 'it2'};
windowFiles  = {
    'window-limited/260507/time_average.nc'
    'window-limited/260528-it1/time_average.nc'
    'window-limited/260607-it2/time_average.nc'
};

limiterLabels = {'it0'; 'it1'; 'it2'; 'it3'};
limiterFiles  = {
    'limiter-limited/260507/time_average.nc'
    'limiter-limited/260515-it1/time_average.nc'
    'limiter-limited/260528-it2/time_average.nc'
    'limiter-limited/260607-it3/time_average.nc'
};

caseName   = {'window_limited';  'limiter_limited'};
caseTitle  = {'Window-limited';  'Limiter-limited'};
iterFiles  = {windowFiles;  limiterFiles};
iterLabels = {windowLabels; limiterLabels};

% ------------------------------------

for iCase = 1:numel(caseName)
    for iIter = 1:numel(iterLabels{iCase})
        ncFile    = iterFiles{iCase}{iIter};
        iterLabel = iterLabels{iCase}{iIter};
        outDir    = fullfile('plots', caseName{iCase}, iterLabel);

        if exist(ncFile, 'file') ~= 2
            fprintf('Skipping (file not found): %s\n', ncFile);
            continue;
        end

        plotCase(ncFile, caseName{iCase}, caseTitle{iCase}, iterLabel, ...
                 outDir, zWindowMin, zWindowMax, nAzimuthalSlices, m_ion);
    end
end

%% ======================= local functions =======================

function plotCase(ncFile, caseName, caseTitle, iterLabel, outDir, ...
                  zWindowMin, zWindowMax, nAzimuthalSlices, m_ion)

    makeFolder(outDir);
    tag = sprintf('%s  %s', caseTitle, iterLabel);   % e.g. "Window-limited it1"

    %% --- read and orient ---
    r     = ncread(ncFile, 'r');
    z     = ncread(ncFile, 'z');
    theta = ncread(ncFile, 'angle');

    fileInfo   = ncinfo(ncFile);
    getVarDims = @(n) {fileInfo.Variables(strcmp({fileInfo.Variables.Name}, n)).Dimensions.Name};
    Ne = timeAverage(ncread(ncFile, 'electron_density'),     getVarDims('electron_density'));
    Te = timeAverage(ncread(ncFile, 'electron_temperature'), getVarDims('electron_temperature'));

    [r, z, theta, reverseDims] = orientCoordinates(r, z, theta);
    if reverseDims
        Ne = reverseDimensions(Ne);
        Te = reverseDimensions(Te);
    end

    %% --- diagnostics ---
    fprintf('\n=== %s ===\n', tag);
    fprintf('Grid (nr x nz x ntheta): %d x %d x %d\n', size(r,1), size(r,2), size(r,3));
    fprintf('r: [%.4g, %.4g] m   z: [%.4g, %.4g] m   theta: [%.4g, %.4g] rad\n', ...
            min(r(:)), max(r(:)), min(z(:)), max(z(:)), min(theta(:)), max(theta(:)));
    fprintf('Ne: [%.4g, %.4g] m^-3   Te: [%.4g, %.4g] eV\n', ...
            min(Ne(:)), max(Ne(:)), min(Te(:)), max(Te(:)));
    if zWindowMin < min(z(:)) || zWindowMax > max(z(:))
        fprintf('WARNING: helicon window [%.2f, %.2f] m outside z range [%.4g, %.4g] m\n', ...
                zWindowMin, zWindowMax, min(z(:)), max(z(:)));
    end

    %% --- CSV ---
    X = r .* cos(theta);
    Y = r .* sin(theta);
    csvFile = fullfile(outDir, sprintf('%s_%s_xyz_ne_te.csv', caseName, iterLabel));
    writeCSV(csvFile, X(:), Y(:), z(:), Ne(:), Te(:));
    fprintf('Wrote CSV.\n');

    %% Figure 1 — R-Z overview
    thetaInd = thetaZeroIndex(theta);
    R0  = r(:, :, thetaInd);
    Z0  = z(:, :, thetaInd);
    Ne0 = Ne(:, :, thetaInd);
    Te0 = Te(:, :, thetaInd);

    fig = figure('Color', 'w', 'Position', [100 100 1100 450]);

    subplot(1, 2, 1);
    contourf(R0, Z0, Ne0, 50, 'LineStyle', 'none');
    cb = colorbar; cb.Label.String = 'N_e [m^{-3}]';
    grid on; xlabel('r [m]'); ylabel('z [m]');
    title([tag '  |  N_e at \theta = 0'], 'Interpreter', 'tex');

    subplot(1, 2, 2);
    contourf(R0, Z0, Te0, 50, 'LineStyle', 'none');
    cb = colorbar; cb.Label.String = 'T_e [eV]';
    grid on; xlabel('r [m]'); ylabel('z [m]');
    title([tag '  |  T_e at \theta = 0'], 'Interpreter', 'tex');

    saveFigure(fig, fullfile(outDir, 'overview_rz'));

    %% Figure 2 — 1D profiles
    zMid = round(size(R0, 2) / 2);
    rMid = round(size(R0, 1) / 2);

    fig = figure('Color', 'w', 'Position', [100 100 1000 750]);

    subplot(2, 2, 1);
    plot(R0(:, zMid), Ne0(:, zMid), 'LineWidth', 2);
    grid on; xlabel('r [m]'); ylabel('N_e [m^{-3}]');
    title(sprintf('Radial N_e  at z = %.3g m', Z0(rMid, zMid)));

    subplot(2, 2, 2);
    plot(R0(:, zMid), Te0(:, zMid), 'LineWidth', 2);
    grid on; xlabel('r [m]'); ylabel('T_e [eV]');
    title(sprintf('Radial T_e  at z = %.3g m', Z0(rMid, zMid)));

    subplot(2, 2, 3);
    plot(Z0(rMid, :)', Ne0(rMid, :)', 'LineWidth', 2);
    grid on; xlabel('z [m]'); ylabel('N_e [m^{-3}]');
    title(sprintf('Axial N_e  at r = %.3g m', R0(rMid, zMid)));

    subplot(2, 2, 4);
    plot(Z0(rMid, :)', Te0(rMid, :)', 'LineWidth', 2);
    grid on; xlabel('z [m]'); ylabel('T_e [eV]');
    title(sprintf('Axial T_e  at r = %.3g m', R0(rMid, zMid)));

    sgtitle([tag '  |  1D profiles at \theta = 0'], 'Interpreter', 'tex');
    saveFigure(fig, fullfile(outDir, 'profiles_1d'));

    %% Figure 3 — azimuthal overview across device
    plotAzimuthalOverview(r, z, theta, Ne, Te, m_ion, nAzimuthalSlices, tag, outDir);

    %% Figure 4 — helicon window theta-z map
    plotHeliconThetaZ(r, z, theta, Ne, Te, zWindowMin, zWindowMax, m_ion, tag, outDir);

    fprintf('Done: %s\n', outDir);
end


% ---------------------------------------------------------------
function plotAzimuthalOverview(r, z, theta, Ne, Te, m_ion, nSlices, tag, outDir)
% 3-row x nSlices-column azimuthal (X-Y) contour maps.
% Row 1: Ne,  Row 2: Te,  Row 3: q_Bohm.
% Each row shares one global color scale. Azimuthal ring is closed at 360 deg.

    e_charge = 1.602e-19;
    gamma    = 3.0;
    kTe = Te * e_charge;
    q   = gamma .* Ne .* kTe .* sqrt(kTe ./ m_ion);

    X = r .* cos(theta);
    Y = r .* sin(theta);

    thetaInd = thetaZeroIndex(theta);
    zLine    = squeeze(z(1, :, thetaInd));
    nz       = numel(zLine);

    zIndices = unique(round(linspace(1, nz, nSlices)));
    nSlices  = numel(zIndices);

    % Extract slices; close the ring by appending the first theta column.
    Xs = cell(1, nSlices);  Ys = cell(1, nSlices);
    Ns = cell(1, nSlices);  Ts = cell(1, nSlices);  Qs = cell(1, nSlices);
    for k = 1:nSlices
        iz = zIndices(k);
        Xsl = squeeze(X(:, iz, :));  Ysl = squeeze(Y(:, iz, :));
        Xs{k} = [Xsl, Xsl(:,1)];    Ys{k} = [Ysl, Ysl(:,1)];
        d = squeeze(Ne(:, iz, :));   Ns{k} = [d, d(:,1)];
        d = squeeze(Te(:, iz, :));   Ts{k} = [d, d(:,1)];
        d = squeeze(q(:,  iz, :));   Qs{k} = [d, d(:,1)];
    end

    % Global color limits per variable (finite values only).
    finLim = @(cells) [min(cellfun(@(c) min(c(isfinite(c(:)))), cells)), ...
                        max(cellfun(@(c) max(c(isfinite(c(:)))), cells))];
    limNe = finLim(Ns);  limTe = finLim(Ts);  limQ = finLim(Qs);

    varData  = {Ns,      Ts,      Qs};
    varLims  = {limNe,   limTe,   limQ};
    varLabel = {'N_e [m^{-3}]', 'T_e [eV]', 'q_{Bohm} [W m^{-2}]'};
    varCmap  = {'parula', 'parula', 'hot'};

    fig = figure('Color', 'w', 'Position', [50 50 max(1200, 240*nSlices+100) 870]);

    for iVar = 1:3
        for k = 1:nSlices
            iz   = zIndices(k);
            zVal = zLine(iz);

            ax = subplot(3, nSlices, (iVar-1)*nSlices + k);
            contourf(Xs{k}, Ys{k}, varData{iVar}{k}, 30, 'LineStyle', 'none');
            colormap(ax, varCmap{iVar});
            set(ax, 'CLim', varLims{iVar});
            axis equal tight; grid on;

            if iVar == 1,  title(sprintf('z = %.3g m', zVal));  end
            if k == 1,     ylabel(varLabel{iVar});
            else,          set(ax, 'YTickLabel', []);  end
            xlabel('X [m]');

            if k == nSlices
                cb = colorbar(ax, 'eastoutside');
                cb.Label.String = varLabel{iVar};
            end
        end
    end

    sgtitle([tag '  |  Azimuthal planes across device'], 'Interpreter', 'none');
    saveFigure(fig, fullfile(outDir, 'azimuthal_overview'));
    fprintf('Wrote azimuthal_overview.png\n');
end


% ---------------------------------------------------------------
function plotHeliconThetaZ(r, z, theta, Ne, Te, zMin, zMax, m_ion, tag, outDir)
% Theta-z map of q_Bohm / Ne / Te at r = r_max (outer wall).
% x-axis: theta [0, 360 deg],  y-axis: z [m] in the helicon window.

    e_charge = 1.602e-19;
    gamma    = 3.0;
    kTe = Te * e_charge;
    q   = gamma .* Ne .* kTe .* sqrt(kTe ./ m_ion);

    thetaRaw = squeeze(theta(1, 1, :));
    thetaDeg = mod(thetaRaw * (180/pi), 360);
    [thetaDeg, sortIdx] = sort(thetaDeg);

    thetaInd = thetaZeroIndex(theta);
    zLine    = squeeze(z(1, :, thetaInd));

    winIdx = find(zLine >= zMin & zLine <= zMax);
    if isempty(winIdx)
        fprintf('Helicon window [%.2f, %.2f] m: no z-points — skipping theta-z plot.\n', zMin, zMax);
        return;
    end
    zWin = zLine(winIdx);
    fprintf('Helicon window: %d z-points (%.4g – %.4g m).\n', numel(winIdx), zWin(1), zWin(end));

    q_wall  = squeeze(q(end,  winIdx, sortIdx));
    Ne_wall = squeeze(Ne(end, winIdx, sortIdx));
    Te_wall = squeeze(Te(end, winIdx, sortIdx));

    [THETA, ZMAP] = meshgrid(thetaDeg, zWin);

    fig = figure('Color', 'w', 'Position', [80 80 1350 430]);

    subplot(1, 3, 1);
    contourf(THETA, ZMAP, q_wall, 50, 'LineStyle', 'none');
    colormap(gca, hot);
    cb = colorbar; cb.Label.String = 'q_{Bohm} [W m^{-2}]';
    xlabel('\theta [deg]'); ylabel('z [m]');
    set(gca, 'XTick', 0:45:360);
    title('Heat flux q_{Bohm}');

    subplot(1, 3, 2);
    contourf(THETA, ZMAP, Ne_wall, 50, 'LineStyle', 'none');
    cb = colorbar; cb.Label.String = 'N_e [m^{-3}]';
    xlabel('\theta [deg]'); ylabel('z [m]');
    set(gca, 'XTick', 0:45:360);
    title('N_e');

    subplot(1, 3, 3);
    contourf(THETA, ZMAP, Te_wall, 50, 'LineStyle', 'none');
    cb = colorbar; cb.Label.String = 'T_e [eV]';
    xlabel('\theta [deg]'); ylabel('z [m]');
    set(gca, 'XTick', 0:45:360);
    title('T_e');

    sgtitle(sprintf('%s  |  Helicon window wall (z = %.1f–%.1f m,  r = r_{max})', ...
            tag, zMin, zMax), 'Interpreter', 'none');
    saveFigure(fig, fullfile(outDir, 'helicon_window_thetaz'));
    fprintf('Wrote helicon_window_thetaz.png\n');
end


% ---------------------------------------------------------------
% Orientation detection (same as compare_profile_convergence.m)
% ---------------------------------------------------------------
function [r, z, theta, flip] = orientCoordinates(r, z, theta)
    rawScore  = orientationScore(r, z, theta);
    rF = reverseDimensions(r);  zF = reverseDimensions(z);  tF = reverseDimensions(theta);
    flip = orientationScore(rF, zF, tF) > rawScore;
    if flip,  r = rF;  z = zF;  theta = tF;  end
end

function score = orientationScore(r, z, theta)
    if ndims(r) ~= 3 || ndims(z) ~= 3 || ndims(theta) ~= 3
        score = -Inf; return;
    end
    i = round(size(r,1)/2);  j = round(size(r,2)/2);  k = round(size(r,3)/2);
    score = dataRange(theta(i,j,:)) + dataRange(r(:,j,k)) + dataRange(z(i,:,k)) ...
          - dataRange(theta(:,j,k)) - dataRange(theta(i,:,k)) ...
          - dataRange(r(i,j,:))     - dataRange(z(i,j,:));
end

function A = reverseDimensions(A)
    if ndims(A) > 1,  A = permute(A, ndims(A):-1:1);  end
end

function v = dataRange(x)
    x = x(isfinite(x(:)));
    if isempty(x), v = -Inf; else, v = max(x) - min(x); end
end


% ---------------------------------------------------------------
function field = timeAverage(field, dimNames)
    if ~isempty(dimNames) && strcmp(dimNames{1}, 't')
        field = squeeze(mean(field, 1));
    end
end

function thetaInd = thetaZeroIndex(theta)
    [~, thetaInd] = min(abs(squeeze(theta(1, 1, :))));
end

function writeCSV(csvFile, x, y, z, Ne, Te)
    fid = fopen(csvFile, 'w');
    fprintf(fid, 'x[m],y[m],z[m],Ne[m^-3],Te[eV]\n');
    fclose(fid);
    writematrix([x, y, z, Ne, Te], csvFile, 'WriteMode', 'append', 'Delimiter', ',');
end

function saveFigure(fig, baseName)
    persistent hasExportGraphics hasSavefig
    if isempty(hasExportGraphics)
        hasExportGraphics = exist('exportgraphics', 'file') == 2;
        hasSavefig        = exist('savefig',        'file') == 2;
    end
    pngFile = [baseName '.png'];
    if hasExportGraphics
        exportgraphics(fig, pngFile, 'Resolution', 300);
    else
        print(fig, pngFile, '-dpng', '-r300');
    end
    if hasSavefig,  savefig(fig, [baseName '.fig']);  end
    close(fig);
end

function makeFolder(folderName)
    if ~exist(folderName, 'dir'),  mkdir(folderName);  end
end
