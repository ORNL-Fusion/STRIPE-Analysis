%% compare_profile_convergence
% Compare radial and axial HERMES plasma profiles between iterations.
%
% Outputs:
%   profile_convergence/<case>/radial_profiles_by_iteration.png
%   profile_convergence/<case>/axial_profiles_by_iteration.png
%   profile_convergence/<case>/convergence_by_iteration.png
%   profile_convergence/convergence_metrics.csv

clear; clc; close all;

%% User Settings

outDir = 'profile_convergence';

% Fractions through the theta = 0 R-Z plane.
% Radial profiles are plotted at these axial locations.
% Axial profiles are plotted at these radial locations.
zFractions = [0.25 0.50 0.75];
rFractions = [0.25 0.50 0.75];

windowLabels = {'it0'; 'it1'; 'it2'};
windowFiles = {
    'window-limited/260507/time_average.nc'
    'window-limited/260528-it1/time_average.nc'
    'window-limited/260607-it2/time_average.nc'
};

limiterLabels = {'it0'; 'it1'; 'it2'; 'it3'};
limiterFiles = {
    'limiter-limited/260507/time_average.nc'
    'limiter-limited/260515-it1/time_average.nc'
    'limiter-limited/260528-it2/time_average.nc'
    'limiter-limited/260607-it3/time_average.nc'
};

caseName = {'window_limited'; 'limiter_limited'};
caseTitle = {'Window-limited'; 'Limiter-limited'};
iterationLabel = {windowLabels; limiterLabels};
netcdfFile = {windowFiles; limiterFiles};

%% Run Comparison

makeFolder(outDir);
allMetrics = {};

for iCase = 1:numel(caseName)
    fprintf('Loading %s...\n', caseTitle{iCase});

    thisOutDir = fullfile(outDir, caseName{iCase});
    makeFolder(thisOutDir);

    nIter = numel(iterationLabel{iCase});
    profile = repmat(emptyProfile(), 1, nIter);

    for iIter = 1:nIter
        profile(iIter) = readHermesProfile( ...
            netcdfFile{iCase}{iIter}, iterationLabel{iCase}{iIter});
    end

    checkMatchingGrid(profile);

    zIndex = indexFromFraction(zFractions, size(profile(1).R, 2));
    rIndex = indexFromFraction(rFractions, size(profile(1).R, 1));

    plotRadialProfiles(profile, zIndex, caseTitle{iCase}, thisOutDir);
    plotAxialProfiles(profile, rIndex, caseTitle{iCase}, thisOutDir);

    allMetrics = [allMetrics; convergenceMetrics(caseName{iCase}, profile, zIndex, rIndex)]; %#ok<AGROW>
end

metrics = cell2table(allMetrics, 'VariableNames', ...
    {'case', 'comparison', 'variable', 'profile', 'location', ...
     'rel_L2', 'rel_Linf', 'n_points'});

writetable(metrics, fullfile(outDir, 'convergence_metrics.csv'));

fprintf('\nConvergence summary for theta = 0 R-Z profiles\n');
for iCase = 1:numel(caseName)
    printConvergenceSummary(metrics, caseName{iCase}, caseTitle{iCase});
end

for iCase = 1:numel(caseName)
    thisOutDir = fullfile(outDir, caseName{iCase});
    plotConvergenceHistory(metrics, caseName{iCase}, caseTitle{iCase}, thisOutDir);
end

fprintf('\nWrote profile plots to: %s\n', outDir);
fprintf('Wrote convergence metrics to: %s\n', fullfile(outDir, 'convergence_metrics.csv'));

%% Local Functions

function profile = emptyProfile()
    profile = struct( ...
        'label', '', ...
        'file', '', ...
        'thetaIndex', NaN, ...
        'R', [], ...
        'Z', [], ...
        'Ne', [], ...
        'Te', []);
end


function profile = readHermesProfile(ncFile, label)
    if exist(ncFile, 'file') ~= 2
        error('File not found: %s', ncFile);
    end

    r = ncread(ncFile, 'r');
    z = ncread(ncFile, 'z');
    theta = ncread(ncFile, 'angle');

    [r, z, theta, reverseDims] = orientCoordinates(r, z, theta);

    Ne = ncread(ncFile, 'electron_density');
    Te = ncread(ncFile, 'electron_temperature');

    if reverseDims
        Ne = reverseDimensions(Ne);
        Te = reverseDimensions(Te);
    end

    thetaIndex = thetaZeroIndex(theta);

    % Read file metadata once instead of once per variable inside timeAverageField.
    fileInfo = ncinfo(ncFile);
    getVarDims = @(name) {fileInfo.Variables(strcmp({fileInfo.Variables.Name}, name)).Dimensions.Name};

    profile = emptyProfile();
    profile.label = label;
    profile.file = ncFile;
    profile.thetaIndex = thetaIndex;
    profile.R = r(:, :, thetaIndex);
    profile.Z = z(:, :, thetaIndex);
    profile.Ne = timeAverageField(Ne, getVarDims('electron_density'));
    profile.Te = timeAverageField(Te, getVarDims('electron_temperature'));

    if ~isequal(size(profile.R), size(profile.Z), size(profile.Ne), size(profile.Te))
        error('Profile array sizes do not match in %s.', ncFile);
    end
end


function field2D = timeAverageField(field, dimNames)
    if numel(dimNames) == 3 && strcmp(dimNames{1}, 't')
        field2D = squeeze(mean(field, 1));
    elseif ndims(field) == 3
        field2D = field(:, :, 1);
    elseif ismatrix(field)
        field2D = field;
    else
        error('Unsupported field dimensions.');
    end
end


function [r, z, theta, reverseDims] = orientCoordinates(r, z, theta)
    rawScore = orientationScore(r, z, theta);

    rReverse = reverseDimensions(r);
    zReverse = reverseDimensions(z);
    thetaReverse = reverseDimensions(theta);
    reverseScore = orientationScore(rReverse, zReverse, thetaReverse);

    reverseDims = reverseScore > rawScore;

    if reverseDims
        r = rReverse;
        z = zReverse;
        theta = thetaReverse;
    end
end


function score = orientationScore(r, z, theta)
    if ndims(r) ~= 3 || ndims(z) ~= 3 || ndims(theta) ~= 3
        score = -Inf;
        return;
    end

    i = round(size(r, 1) / 2);
    j = round(size(r, 2) / 2);
    k = round(size(r, 3) / 2);

    thetaAlongK = dataRange(theta(i, j, :));
    rAlongI = dataRange(r(:, j, k));
    zAlongJ = dataRange(z(i, :, k));

    thetaLeak = dataRange(theta(:, j, k)) + dataRange(theta(i, :, k));
    rLeak = dataRange(r(i, j, :));
    zLeak = dataRange(z(i, j, :));

    score = thetaAlongK + rAlongI + zAlongJ - thetaLeak - rLeak - zLeak;
end


function A = reverseDimensions(A)
    if ndims(A) > 1
        A = permute(A, ndims(A):-1:1);
    end
end


function value = dataRange(x)
    x = x(:);
    x = x(isfinite(x));

    if isempty(x)
        value = -Inf;
    else
        value = max(x) - min(x);
    end
end


function thetaIndex = thetaZeroIndex(theta)
    thetaLine = squeeze(theta(1, 1, :));
    [~, thetaIndex] = min(abs(thetaLine));
end


function index = indexFromFraction(fraction, n)
    index = round(1 + fraction(:)' .* (n - 1));
    index = max(1, min(n, index));
    index = unique(index, 'stable');
end


function checkMatchingGrid(profile)
    rRef = profile(1).R;
    zRef = profile(1).Z;
    tol = 1e-10;

    for k = 2:numel(profile)
        if ~isequal(size(rRef), size(profile(k).R))
            error('Grid size mismatch between %s and %s.', profile(1).label, profile(k).label);
        end

        maxRDiff = max(abs(rRef(:) - profile(k).R(:)));
        maxZDiff = max(abs(zRef(:) - profile(k).Z(:)));

        if maxRDiff > tol || maxZDiff > tol
            error('Grid mismatch between %s and %s: max dR = %.3e, max dZ = %.3e.', ...
                profile(1).label, profile(k).label, maxRDiff, maxZDiff);
        end
    end
end


function plotRadialProfiles(profile, zIndex, caseTitle, outDir)
    nRows = numel(zIndex);
    color = lines(numel(profile));
    midR = round(size(profile(1).Z, 1) / 2);

    fig = figure('Color', 'w', 'Position', [100 100 1100 320*nRows]);

    for iLoc = 1:nRows
        j = zIndex(iLoc);
        zValue = profile(1).Z(midR, j);

        subplot(nRows, 2, 2*iLoc - 1);
        hold on;
        for k = 1:numel(profile)
            [x, y] = cleanProfile(profile(k).R(:, j), profile(k).Ne(:, j));
            plot(x, y, 'LineWidth', 1.8, 'Color', color(k, :), ...
                'DisplayName', profile(k).label);
        end
        grid on;
        xlabel('r [m]');
        ylabel('N_e [m^{-3}]');
        title(sprintf('%s radial N_e at z = %.4g m', caseTitle, zValue));
        legend('Location', 'best');

        subplot(nRows, 2, 2*iLoc);
        hold on;
        for k = 1:numel(profile)
            [x, y] = cleanProfile(profile(k).R(:, j), profile(k).Te(:, j));
            plot(x, y, 'LineWidth', 1.8, 'Color', color(k, :), ...
                'DisplayName', profile(k).label);
        end
        grid on;
        xlabel('r [m]');
        ylabel('T_e [eV]');
        title(sprintf('%s radial T_e at z = %.4g m', caseTitle, zValue));
        legend('Location', 'best');
    end

    saveFigure(fig, fullfile(outDir, 'radial_profiles_by_iteration'));
end


function plotAxialProfiles(profile, rIndex, caseTitle, outDir)
    nRows = numel(rIndex);
    color = lines(numel(profile));
    midZ = round(size(profile(1).R, 2) / 2);

    fig = figure('Color', 'w', 'Position', [150 150 1100 320*nRows]);

    for iLoc = 1:nRows
        i = rIndex(iLoc);
        rValue = profile(1).R(i, midZ);

        subplot(nRows, 2, 2*iLoc - 1);
        hold on;
        for k = 1:numel(profile)
            [x, y] = cleanProfile(profile(k).Z(i, :), profile(k).Ne(i, :));
            plot(x, y, 'LineWidth', 1.8, 'Color', color(k, :), ...
                'DisplayName', profile(k).label);
        end
        grid on;
        xlabel('z [m]');
        ylabel('N_e [m^{-3}]');
        title(sprintf('%s axial N_e at r = %.4g m', caseTitle, rValue));
        legend('Location', 'best');

        subplot(nRows, 2, 2*iLoc);
        hold on;
        for k = 1:numel(profile)
            [x, y] = cleanProfile(profile(k).Z(i, :), profile(k).Te(i, :));
            plot(x, y, 'LineWidth', 1.8, 'Color', color(k, :), ...
                'DisplayName', profile(k).label);
        end
        grid on;
        xlabel('z [m]');
        ylabel('T_e [eV]');
        title(sprintf('%s axial T_e at r = %.4g m', caseTitle, rValue));
        legend('Location', 'best');
    end

    saveFigure(fig, fullfile(outDir, 'axial_profiles_by_iteration'));
end


function [x, y] = cleanProfile(x, y)
    x = x(:);
    y = y(:);

    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    [x, order] = sort(x);
    y = y(order);
end


function rows = convergenceMetrics(caseName, profile, zIndex, rIndex)
    nIter = numel(profile) - 1;
    nRowsPerIter = 2 * (1 + numel(zIndex) + numel(rIndex));
    rows = cell(nIter * nRowsPerIter, 8);
    iRow = 0;

    % Grid is verified identical across iterations, so midpoints are constant.
    midR = round(size(profile(1).Z, 1) / 2);
    midZ = round(size(profile(1).R, 2) / 2);

    for k = 2:numel(profile)
        old = profile(k - 1);
        new = profile(k);
        comparison = [old.label '_to_' new.label];

        [rows, iRow] = addMetric(rows, iRow, caseName, comparison, 'Ne', 'theta0_rz_slice', 'all', old.Ne, new.Ne);
        [rows, iRow] = addMetric(rows, iRow, caseName, comparison, 'Te', 'theta0_rz_slice', 'all', old.Te, new.Te);

        for j = zIndex
            location = sprintf('z=%.6g m', new.Z(midR, j));
            [rows, iRow] = addMetric(rows, iRow, caseName, comparison, 'Ne', 'radial', location, old.Ne(:, j), new.Ne(:, j));
            [rows, iRow] = addMetric(rows, iRow, caseName, comparison, 'Te', 'radial', location, old.Te(:, j), new.Te(:, j));
        end

        for i = rIndex
            location = sprintf('r=%.6g m', new.R(i, midZ));
            [rows, iRow] = addMetric(rows, iRow, caseName, comparison, 'Ne', 'axial', location, old.Ne(i, :), new.Ne(i, :));
            [rows, iRow] = addMetric(rows, iRow, caseName, comparison, 'Te', 'axial', location, old.Te(i, :), new.Te(i, :));
        end
    end

    rows = rows(1:iRow, :);
end


function [rows, iRow] = addMetric(rows, iRow, caseName, comparison, variable, profileType, location, oldValue, newValue)
    oldValue = oldValue(:);
    newValue = newValue(:);

    keep = isfinite(oldValue) & isfinite(newValue);
    oldValue = oldValue(keep);
    newValue = newValue(keep);

    nPoints = numel(oldValue);

    if nPoints == 0
        relL2 = NaN;
        relLinf = NaN;
    else
        difference = newValue - oldValue;
        relL2 = norm(difference, 2) / max([norm(oldValue, 2), norm(newValue, 2), eps]);
        relLinf = norm(difference, Inf) / max([norm(oldValue, Inf), norm(newValue, Inf), eps]);
    end

    iRow = iRow + 1;
    rows(iRow, :) = {caseName, comparison, variable, profileType, location, relL2, relLinf, nPoints};
end


function plotConvergenceHistory(metrics, caseName, caseTitle, outDir)
    caseRows = metrics(strcmp(metrics.case, caseName), :);

    if isempty(caseRows)
        return;
    end

    comparison = unique(caseRows.comparison, 'stable');
    comparisonLabel = prettyComparisonLabels(comparison);
    x = 1:numel(comparison);

    % Compute unique locations once — they are the same for Ne and Te.
    radialLocations = unique(caseRows.location(strcmp(caseRows.profile, 'radial')), 'stable');
    axialLocations  = unique(caseRows.location(strcmp(caseRows.profile, 'axial')),  'stable');

    fig = figure('Color', 'w', 'Position', [200 200 1100 450]);
    variableName = {'Ne', 'Te'};
    variableLabel = {'N_e', 'T_e'};

    for iVar = 1:numel(variableName)
        subplot(1, 2, iVar);
        hold on;

        variableRows = caseRows(strcmp(caseRows.variable, variableName{iVar}), :);

        plotMetricGroup(variableRows, comparison, x, ...
            'theta0_rz_slice', 'all', 'o-', 2.4, 'R-Z slice');

        for iLoc = 1:numel(radialLocations)
            plotMetricGroup(variableRows, comparison, x, ...
                'radial', radialLocations{iLoc}, '.-', 1.2, ...
                ['radial ' prettyLocation(radialLocations{iLoc})]);
        end

        for iLoc = 1:numel(axialLocations)
            plotMetricGroup(variableRows, comparison, x, ...
                'axial', axialLocations{iLoc}, 'x--', 1.2, ...
                ['axial ' prettyLocation(axialLocations{iLoc})]);
        end

        grid on;
        set(gca, 'XTick', x, 'XTickLabel', comparisonLabel);
        xtickangle(25);
        xlim([0.8 numel(comparison) + 0.2]);
        ylabel('relative L2 change');
        title(sprintf('%s %s convergence', caseTitle, variableLabel{iVar}));
        legend('Location', 'eastoutside');
    end

    saveFigure(fig, fullfile(outDir, 'convergence_by_iteration'));
end


function printConvergenceSummary(metrics, caseName, caseTitle)
    caseRows = metrics(strcmp(metrics.case, caseName), :);
    comparison = unique(caseRows.comparison, 'stable');

    % Filter to theta=0 R-Z slice rows once; then split by variable inside the loop.
    sliceRows = caseRows(strcmp(caseRows.profile, 'theta0_rz_slice') & ...
                         strcmp(caseRows.location, 'all'), :);

    fprintf('\n%s\n', caseTitle);
    fprintf('%-12s %12s %12s %12s %12s\n', ...
        'step', 'Ne L2', 'Ne Linf', 'Te L2', 'Te Linf');

    for k = 1:numel(comparison)
        compRows = sliceRows(strcmp(sliceRows.comparison, comparison{k}), :);
        neRow = compRows(strcmp(compRows.variable, 'Ne'), :);
        teRow = compRows(strcmp(compRows.variable, 'Te'), :);

        fprintf('%-12s %12.4g %12.4g %12.4g %12.4g\n', ...
            prettyComparisonLabels(comparison{k}), ...
            neRow.rel_L2(1), neRow.rel_Linf(1), ...
            teRow.rel_L2(1), teRow.rel_Linf(1));
    end
end


function label = prettyComparisonLabels(comparison)
    label = strrep(comparison, '_to_', ' -> ');
end


function label = prettyLocation(location)
    token = regexp(location, '^([rz])=([0-9eE+\-.]+) m$', 'tokens', 'once');

    if isempty(token)
        label = location;
        return;
    end

    coordinate = token{1};
    value = str2double(token{2});

    if strcmp(coordinate, 'r')
        label = sprintf('r = %.3f m', value);
    else
        label = sprintf('z = %.2f m', value);
    end
end


function plotMetricGroup(rows, comparison, x, profileType, location, lineSpec, lineWidth, label)
    % Filter to the target profile/location once, then align to comparison order.
    subset = rows(strcmp(rows.profile, profileType) & strcmp(rows.location, location), :);

    y = NaN(size(x));
    if ~isempty(subset)
        [~, ia, ib] = intersect(comparison, subset.comparison, 'stable');
        y(ia) = subset.rel_L2(ib);
    end

    if any(isfinite(y))
        plot(x, y, lineSpec, 'LineWidth', lineWidth, 'DisplayName', label);
    end
end


function saveFigure(fig, baseName)
    persistent hasExportGraphics hasSavefig
    if isempty(hasExportGraphics)
        hasExportGraphics = exist('exportgraphics', 'file') == 2;
        hasSavefig        = exist('savefig',        'file') == 2;
    end

    pngFile = [baseName '.png'];
    figFile = [baseName '.fig'];

    if hasExportGraphics
        exportgraphics(fig, pngFile, 'Resolution', 300);
    else
        print(fig, pngFile, '-dpng', '-r300');
    end

    if hasSavefig
        savefig(fig, figFile);
    end

    close(fig);
end


function makeFolder(folderName)
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
end
