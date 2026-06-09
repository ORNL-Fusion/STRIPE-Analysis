%% -------------------- Neon fractional abundance vs R at Z=0.63 --------------------
z_target = 0.0;   % [m]

[~, iz0] = min(abs(Zg - z_target));

fprintf('\nNeon fractional abundance vs R\n');
fprintf('Requested Z = %.6g m, nearest grid Z = %.6g m\n', z_target, Zg(iz0));
fprintf('Using R range: %.6g to %.6g m\n', R(1), R(end));
fprintf('niA size = [%d %d %d]\n', size(niA,1), size(niA,2), size(niA,3));

%% ----- Find Neon species -----
idx_neon = find(Z_all == 10);
q_neon   = q_all(idx_neon);

% If no Neon metadata found, use fallback block
if isempty(idx_neon)
    warning('No Z_all == 10 entries found. Using fallback indices 6:16.');
    idx_neon = 6:16;
    idx_neon = idx_neon(idx_neon <= ns);
    q_neon = 0:(numel(idx_neon)-1);
else
    % If all Neon charges are identical, charge metadata is bad
    if numel(unique(q_neon)) == 1
        warning('Neon charge_number metadata is not usable. Relabeling Neon block as Ne0+...Ne10+.');
        q_neon = 0:(numel(idx_neon)-1);
    end
end

% Ensure column vectors
idx_neon = idx_neon(:);
q_neon   = q_neon(:);

% Sort by plotted charge state
[q_neon, isort] = sort(q_neon);
idx_neon = idx_neon(isort);

%% ----- Print species mapping -----
fprintf('\nUsing Neon species:\n');
fprintf('Curve\tSpeciesIndex\tPlotLabel\tMetadata_Z\tMetadata_q\n');

for jj = 1:numel(idx_neon)
    fprintf('%d\t%d\t\tNe%d+\t\t%d\t\t%d\n', ...
        jj, idx_neon(jj), q_neon(jj), ...
        Z_all(idx_neon(jj)), q_all(idx_neon(jj)));
end

%% ----- Extract Neon density vs R -----
nNe_R = zeros(nR, numel(idx_neon));

for jj = 1:numel(idx_neon)
    nNe_R(:,jj) = niA(:, iz0, idx_neon(jj));
end

nNe_R(~isfinite(nNe_R)) = 0;
nNe_R = max(nNe_R, 0);

%% ----- Fractional abundance -----
nNe_total_R = sum(nNe_R, 2);

fracNe_R = zeros(size(nNe_R));
valid = nNe_total_R > 0;

fracNe_R(valid,:) = nNe_R(valid,:) ./ nNe_total_R(valid);

%% ----- Dominant charge diagnostic -----
[~, imax] = max(fracNe_R, [], 2);
q_dominant = q_neon(imax);

fprintf('\nAt outer R = %.6g m, dominant Neon charge = Ne%d+\n', ...
        R(end), q_dominant(end));

fprintf('\nMATLAB outer-R raw Neon densities:\n');
for jj = 1:numel(idx_neon)
    fprintf('idx=%d, label=Ne%d+, n=%.6e, frac=%.4f\n', ...
        idx_neon(jj), q_neon(jj), ...
        niA(end, iz0, idx_neon(jj)), fracNe_R(end,jj));
end

%% ----- Optional direct checksum-style diagnostic -----
fprintf('\nMATLAB direct species check at outer R, Z index %d:\n', iz0);
for idx = 6:min(16, ns)
    fprintf('idx=%d, Z=%d, q=%d, n=%.6e\n', ...
        idx, Z_all(idx), q_all(idx), niA(end, iz0, idx));
end

%% ----- Plot -----
figNeR = figure('Color','w','Name','Neon Fractional Abundance vs R at Z=0.63');
axNeR = axes(figNeR);

plot(axNeR, R, fracNe_R, 'LineWidth', 2);

xlabel(axNeR,'R [m]');
ylabel(axNeR,'Fractional abundance n_{Ne^{q+}} / \Sigma_q n_{Ne^{q+}}');

title(axNeR, sprintf('Neon charge-state fractional abundance vs R at Z = %.3f m', ...
      Zg(iz0)), ...
      'FontSize',14,'FontWeight','bold');

grid(axNeR,'on');
xlim(axNeR,[min(R), max(R)]);
ylim(axNeR,[0, 1]);

legend_labels = arrayfun(@(q) sprintf('Ne^{%d+}', q), q_neon, ...
                         'UniformOutput', false);

legend(axNeR, legend_labels, 'Location','bestoutside');

if save_png
    exportgraphics(figNeR, 'neon_fractional_abundance_vs_R_z063.png', ...
                   'Resolution', 300);
end

%% ----- Bar chart at outer R -----
figNeBar = figure('Color','w','Name','Neon Fractional Abundance Bar Chart at outer R');
axNeBar = axes(figNeBar);

bar(axNeBar, q_neon, fracNe_R(end,:), 0.75);

xlabel(axNeBar,'Neon charge state q');
ylabel(axNeBar,'Fractional abundance');

title(axNeBar, sprintf('Neon ion fractional abundance at R = %.3f m, Z = %.3f m', ...
      R(end), Zg(iz0)), ...
      'FontSize',14,'FontWeight','bold');

grid(axNeBar,'on');
ylim(axNeBar,[0, 1]);
xlim(axNeBar,[min(q_neon)-0.5, max(q_neon)+0.5]);

xticks(axNeBar, q_neon);
xticklabels(axNeBar, arrayfun(@(q) sprintf('Ne^{%d+}', q), q_neon, ...
                              'UniformOutput', false));

if save_png
    exportgraphics(figNeBar, 'neon_fractional_abundance_bar_outerR_z063.png', ...
                   'Resolution', 300);
end