%% Neon coronal charge-state distribution (Ne0..Ne10+)
% Computes coronal fractions f_q(Te,ne) from ionization/recombination rates:
%   f_{q+1}/f_q = S_q / alpha_eff_{q+1}
% where S_q is ionization q->q+1 and alpha_{q+1} is recombination q+1->q.
%
% Rate inputs can be:
%   1) ADAS ADF11 files (recommended): scdXX_ne.dat and acdXX_ne.dat
%      Optional: ccdXX_ne.dat for CX recombination contribution
%   2) Manual rate arrays supplied in this script.
%
% Outputs:
%   - coronal fractions vs Te plot
%   - table print at Te_ref
%
% Notes:
%   - ADAS ADF11 SCD/ACD are in cm^3/s; converted internally to m^3/s.
%   - If CCD is enabled, effective recombination is:
%       alpha_eff = ACD + ne * CCD
%     where CCD has units cm^6/s and ACD cm^3/s (after conversion to SI:
%       ACD [m^3/s], CCD [m^6/s], ne [m^-3], so ne*CCD [m^3/s]).
%   - For ADF11 mapping here:
%       S_q uses SCD block Z1=q+1  (q=0..9 for Neon)
%       alpha_{q+1} uses ACD block Z1=q+1 (q+1=1..10)
%   - If your ACD file uses a different block convention, adjust map in get_ne_rates().

clear; clc;

%% -------------------- User controls --------------------
mode = 'adas';                      % 'adas' or 'manual'
fixed_Te_mode = true;               % true: evaluate only at Te_fixed_eV
Te_fixed_eV = 4;                    % fixed Te [eV] when fixed_Te_mode=true
Te_eV = logspace(log10(1), log10(1000), 200).';   % temperature grid for sweep mode
ne_m3 = 1e17;                       % electron density [m^-3]
Te_ref = 4;                         % report fractions at this Te [eV]
adas_set = 89;                      % 85 or 89
use_ccd = false;                    % pure coronal default; set true to include CCD term

% ADAS file search locations (first existing file is used)
adas_dirs = { ...
    '/Users/78k/ORNL Dropbox/Atul Kumar/work/STRIPE-Analysis/iter/simulations/solps_iter/ADAS', ...
    pwd ...
    };

% Manual fallback rates if mode='manual'
% Dimensions: [nTe x 10], columns are q=0..9 for S and q=1..10 for alpha.
S_manual = [];
A_manual = [];

if fixed_Te_mode
    Te_eV = Te_fixed_eV;
    Te_ref = Te_fixed_eV;
end

%% -------------------- Get rates --------------------
switch lower(mode)
    case 'adas'
        scd_file = resolve_adas_file(sprintf('scd%d_ne.dat', adas_set), adas_dirs);
        acd_file = resolve_adas_file(sprintf('acd%d_ne.dat', adas_set), adas_dirs);
        if use_ccd
            ccd_file = resolve_adas_file(sprintf('ccd%d_ne.dat', adas_set), adas_dirs, true);
        else
            ccd_file = '';
        end
        [S_m3s, A_m3s] = get_ne_rates_from_adf11(scd_file, acd_file, ccd_file, Te_eV, ne_m3);
    case 'manual'
        if isempty(S_manual) || isempty(A_manual)
            error('Set S_manual and A_manual for mode=''manual''.');
        end
        S_m3s = S_manual;
        A_m3s = A_manual;
    otherwise
        error('Unknown mode: %s', mode);
end

if size(S_m3s,2) ~= 10 || size(A_m3s,2) ~= 10
    error('Expected 10 Neon transitions for S and A.');
end

%% -------------------- Coronal fractions --------------------
f = zeros(numel(Te_eV), 11);   % columns q=0..10
for i = 1:numel(Te_eV)
    f(i,:) = coronal_fractions_row(S_m3s(i,:), A_m3s(i,:));
end

%% -------------------- Plot --------------------
cmap = turbo(11);
if numel(Te_eV) == 1
    fq = f(1,:);
    fq = fq ./ max(sum(fq), 1e-60);  % ensure normalized abundance

    figure('Color','w');
    bar(0:10, fq, 'FaceColor', [0.2 0.45 0.8]);
    grid on; box on;
    xlabel('Ne charge state q');
    ylabel('Fractional abundance f_q');
    title(sprintf('Neon fractional abundance at T_e = %.3f eV (n_e = %.2e m^{-3})', Te_eV, ne_m3));

    figure('Color','w');
    bar(0:10, 100*fq, 'FaceColor', [0.85 0.4 0.2]);
    grid on; box on;
    xlabel('Ne charge state q');
    ylabel('Fractional abundance [%]');
    title(sprintf('Neon fractional abundance [%%] at T_e = %.3f eV (n_e = %.2e m^{-3})', Te_eV, ne_m3));

    figure('Color','w');
    semilogy(0:10, max(fq,1e-40), 'o-', 'LineWidth', 2);
    grid on; box on;
    xlabel('Ne charge state q');
    ylabel('Fractional abundance f_q (log)');
    title(sprintf('Neon fractional abundance (log) at fixed T_e = %.3f eV', Te_eV));
    ylim([1e-40 1]);
else
    figure('Color','w');
    hold on; box on; grid on;
    for q = 0:10
        semilogx(Te_eV, f(:,q+1), 'LineWidth', 2, 'Color', cmap(q+1,:), ...
            'DisplayName', sprintf('Ne^{%d+}', q));
    end
    xlabel('T_e [eV]');
    ylabel('Coronal fraction f_q');
    title(sprintf('Neon Coronal Charge-State Distribution (linear scale, n_e = %.2e m^{-3})', ne_m3));
    legend('Location','eastoutside');

    % Log-scale view to expose low-abundance charge states
    figure('Color','w');
    hold on; box on; grid on;
    for q = 0:10
        semilogx(Te_eV, max(f(:,q+1), 1e-40), 'LineWidth', 2, 'Color', cmap(q+1,:), ...
            'DisplayName', sprintf('Ne^{%d+}', q));
    end
    set(gca, 'YScale', 'log');
    xlabel('T_e [eV]');
    ylabel('Coronal fraction f_q (log scale)');
    title(sprintf('Neon Coronal Charge-State Distribution (log scale, n_e = %.2e m^{-3})', ne_m3));
    ylim([1e-40 1]);
    legend('Location','eastoutside');
end

%% -------------------- Table at reference Te --------------------
[~, iref] = min(abs(Te_eV - Te_ref));
fprintf('\nNeon coronal fractions at Te = %.3f eV, ne = %.3e m^-3\n', Te_eV(iref), ne_m3);
fprintf(' q    f_q\n');
for q = 0:10
    fprintf('%2d   %.6e\n', q, f(iref,q+1));
end
fprintf('sum f_q = %.6f\n', sum(f(iref,:)));
fprintf('\nCharge-state ranges over Te grid:\n');
for q = 0:10
    fprintf('Ne^{%d+}: min=%.3e  max=%.3e\n', q, min(f(:,q+1)), max(f(:,q+1)));
end

%% -------------------- Local functions --------------------
function fq = coronal_fractions_row(Sq, Aq1)
    % Sq: ionization rates q->q+1 for q=0..9  [1x10]
    % Aq1: recombination rates q+1->q for q+1=1..10 [1x10]
    epsv = 1e-60;
    g = zeros(1,11);
    g(1) = 1;  % Ne0 baseline
    for q = 0:9
        g(q+2) = g(q+1) * Sq(q+1) / max(Aq1(q+1), epsv);
    end
    s = sum(g);
    if ~isfinite(s) || s <= 0
        fq = zeros(1,11);
        fq(1) = 1;
    else
        fq = g / s;
    end
end

function [S_m3s, Aeff_m3s] = get_ne_rates_from_adf11(scd_file, acd_file, ccd_file, Te_eV, ne_m3)
    if ~isfile(scd_file)
        error('SCD file not found: %s', scd_file);
    end
    if ~isfile(acd_file)
        error('ACD file not found: %s', acd_file);
    end
    if ~isempty(ccd_file) && ~isfile(ccd_file)
        error('CCD file not found: %s', ccd_file);
    end

    scd = read_adf11_blocks(scd_file);
    acd = read_adf11_blocks(acd_file);
    if ~isempty(ccd_file)
        ccd = read_adf11_blocks(ccd_file);
    else
        ccd = [];
    end

    ne_cm3 = ne_m3 * 1e-6;
    logTe = log10(Te_eV(:));
    logNe = log10(ne_cm3);

    nT = numel(Te_eV);
    S_cm3s = zeros(nT,10);
    A_cm3s = zeros(nT,10);
    C_cm6s = zeros(nT,10);
    for q = 0:9
        z1 = q + 1;
        if z1 > numel(scd.blocks) || isempty(scd.blocks{z1})
            avail = find(~cellfun(@isempty, scd.blocks));
            error('SCD block Z1=%d missing in %s (available Z1: %s)', ...
                z1, scd_file, mat2str(avail));
        end
        if z1 > numel(acd.blocks) || isempty(acd.blocks{z1})
            avail = find(~cellfun(@isempty, acd.blocks));
            error('ACD block Z1=%d missing in %s (available Z1: %s)', ...
                z1, acd_file, mat2str(avail));
        end
        S_cm3s(:,q+1) = interp_block(scd, z1, logTe, logNe);
        A_cm3s(:,q+1) = interp_block(acd, z1, logTe, logNe);
        if ~isempty(ccd)
            if z1 > numel(ccd.blocks) || isempty(ccd.blocks{z1})
                avail = find(~cellfun(@isempty, ccd.blocks));
                error('CCD block Z1=%d missing in %s (available Z1: %s)', ...
                    z1, ccd_file, mat2str(avail));
            end
            C_cm6s(:,q+1) = interp_block(ccd, z1, logTe, logNe);
        end
    end

    % cm^3/s -> m^3/s
    S_m3s = S_cm3s * 1e-6;
    A_m3s = A_cm3s * 1e-6;
    C_m6s = C_cm6s * 1e-12;

    % effective recombination
    Aeff_m3s = A_m3s + ne_m3 .* C_m6s;
end

function fpath = resolve_adas_file(fname, dirs, allow_missing)
    if nargin < 3
        allow_missing = false;
    end
    fpath = '';
    for i = 1:numel(dirs)
        cand = fullfile(dirs{i}, fname);
        if isfile(cand)
            fpath = cand;
            return;
        end
    end
    if allow_missing
        return;
    end
    error('Required ADAS file not found: %s', fname);
end

function out = interp_block(adf, z1, logTe, logNe)
    M = adf.blocks{z1};  % [nTe x nNe], values are log10(rate[cm3/s])
    out = zeros(size(logTe));
    for i = 1:numel(logTe)
        v = bilinear(adf.logTe, adf.logNe, M, logTe(i), logNe);
        out(i) = 10.^v;
    end
end

function v = bilinear(xg, yg, M, x, y)
    [ix0, ix1, tx] = bracket(xg, x);
    [iy0, iy1, ty] = bracket(yg, y);
    v00 = M(ix0, iy0); v01 = M(ix0, iy1);
    v10 = M(ix1, iy0); v11 = M(ix1, iy1);
    v0 = v00*(1-ty) + v01*ty;
    v1 = v10*(1-ty) + v11*ty;
    v  = v0*(1-tx) + v1*tx;
end

function [i0, i1, t] = bracket(g, x)
    if x <= g(1), i0 = 1; i1 = 1; t = 0; return; end
    if x >= g(end), i0 = numel(g); i1 = numel(g); t = 0; return; end
    i1 = find(g >= x, 1, 'first');
    i0 = i1 - 1;
    t = (x - g(i0)) / (g(i1) - g(i0));
end

function adf = read_adf11_blocks(filename)
    txt = readlines(filename);
    txt = cellstr(txt);
    txt = txt(:);

    % header dims from first line
    nums = sscanf(txt{1}, '%f');
    if numel(nums) < 5
        error('Cannot parse ADF11 header in %s', filename);
    end
    nBlocks = round(nums(1));
    nNe = nums(2);
    nTe = nums(3);

    % Start after first header line and skip separators/blank lines.
    i = 2;
    while i <= numel(txt) && (isempty(strtrim(txt{i})) || contains(txt{i}, '-'))
        i = i + 1;
    end

    % read log10(ne)
    logNe = [];
    while numel(logNe) < nNe
        if i > numel(txt)
            error('Unexpected EOF while reading logNe in %s', filename);
        end
        vals = sscanf(txt{i}, '%f').';
        logNe = [logNe vals]; %#ok<AGROW>
        i = i + 1;
    end
    logNe = logNe(1:nNe);

    % read log10(Te)
    logTe = [];
    while numel(logTe) < nTe
        if i > numel(txt)
            error('Unexpected EOF while reading logTe in %s', filename);
        end
        vals = sscanf(txt{i}, '%f').';
        logTe = [logTe vals]; %#ok<AGROW>
        i = i + 1;
    end
    logTe = logTe(1:nTe);

    % Use indexed cell storage by Z1.
    blocks = cell(1, max(16, nBlocks));
    nLoaded = 0;

    while i <= numel(txt) && nLoaded < nBlocks
        ln = txt{i};
        if ~contains(ln, 'Z1=')
            i = i + 1;
            continue;
        end

        % Robust Z1 parser for headers like:
        % "----/ Z1= 1   / DATE= ..."
        z = sscanf(extractAfter(ln, 'Z1='), '%f', 1);
        if isempty(z) || ~isfinite(z)
            tok = regexp(ln, 'Z1=\s*([0-9]+)', 'tokens', 'once');
            if isempty(tok)
                i = i + 1;
                continue;
            end
            z = str2double(tok{1});
        end
        z = round(z);

        i = i + 1; % move to first data line in this block
        M = zeros(nTe, nNe);
        for it = 1:nTe
            row = [];
            while numel(row) < nNe
                if i > numel(txt)
                    error('Unexpected EOF in block Z1=%d of %s', z, filename);
                end
                if contains(txt{i}, 'Z1=')
                    error('Malformed ADF11 block Z1=%d in %s (early next header).', z, filename);
                end
                vals = sscanf(txt{i}, '%f').';
                if ~isempty(vals)
                    row = [row vals]; %#ok<AGROW>
                end
                i = i + 1;
            end
            M(it,:) = row(1:nNe);
        end

        if z > numel(blocks)
            blocks{z} = [];
        end
        blocks{z} = M;
        nLoaded = nLoaded + 1;
    end

    adf.logNe = logNe;
    adf.logTe = logTe;
    adf.blocks = blocks;
end
