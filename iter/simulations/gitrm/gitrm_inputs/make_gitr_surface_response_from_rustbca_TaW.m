function make_gitr_surface_response_from_rustbca_TaW(outfile, rustbca_dir)
% make_gitr_surface_response_from_rustbca_TaW
% (Same header & behavior as before, but with W->W comparison plots appended)
%
% Usage:
%   make_gitr_surface_response_from_rustbca_TaW;
%   make_gitr_surface_response_from_rustbca_TaW('TaW_surface_response.nc','/path/to/rustBCA_dir');

    if nargin < 1 || isempty(outfile)
        outfile = 'TaW_surface_response.nc';
    end
    if nargin < 2 || isempty(rustbca_dir)
        rustbca_dir = '.'; % current folder
    end

    %% ---- Output grids (edit if you want different resolution) ----------
    angle_dim      = 361;       % e.g. 0–89 deg
    energy_dim     = 200;       % e.g. 10–1e4 eV
    projectile_dim = 2;         % Ta, W
    target_dim     = 2;         % Ta, W

    Angles   = linspace(0.0, 89.0, angle_dim);                 % degrees
    Energies = logspace(log10(10), log10(1e4), energy_dim);    % eV

    % Species labels (atomic numbers)
    Z_Ta = 73;
    Z_W  = 74;
    Projectiles = single([Z_Ta, Z_W]);   % index 1=Ta, 2=W
    Targets     = single([Z_Ta, Z_W]);   % index 1=Ta, 2=W

    projNames = {'Ta','W'};
    targNames = {'Ta','W'};

    %% ---- Allocate arrays (MATLAB logical order) ------------------------
    Physical_Sputtering      = zeros(projectile_dim, target_dim, angle_dim, energy_dim);
    Reflection_Yield         = zeros(projectile_dim, target_dim, angle_dim, energy_dim);
    Energy_Reflection_Yield  = zeros(projectile_dim, target_dim, angle_dim, energy_dim);

    %% ---- Fill arrays from RustBCA tables -------------------------------
    for ip = 1:projectile_dim
        for it = 1:target_dim

            pName = projNames{ip};
            tName = targNames{it};

            rustfile = fullfile(rustbca_dir, sprintf('RustBCA_%son%s.nc', pName, tName));

            if exist(rustfile,'file') ~= 2
                fprintf('[WARN] Missing RustBCA file for %s->%s: %s (leaving zeros)\n', pName, tName, rustfile);
                continue
            end

            fprintf('Loading RustBCA: %s\n', rustfile);

            % --- Read RustBCA variables (native grids) ---
            E_r = ncread(rustfile,'E'); E_r = E_r(:);    % [nE x 1]
            A_r = ncread(rustfile,'A'); A_r = A_r(:);    % [nA x 1]

            spyld_r = ncread(rustfile,'spyld');          % [nA x nE]
            rfyld_r = ncread(rustfile,'rfyld');          % [nA x nE]

            % Energy reflection coefficient from reflected energy distribution
            % energyDistRef dims: [nEdistBinsRef x nA x nE]
            has_ref = true;
            try
                egrid_ref  = ncread(rustfile,'eDistEgrid_ref'); egrid_ref = egrid_ref(:);
                eDistRef   = ncread(rustfile,'energyDistRef');
            catch
                has_ref = false;
                fprintf('[WARN] %s: missing energyDistRef/eDistEgrid_ref, setting Energy_Reflection_Yield = 0 for this pair.\n', rustfile);
            end

            % --- Build RE_r (nA x nE) if possible ---
            if has_ref
                RE_r = compute_energy_reflection(E_r, A_r, rfyld_r, egrid_ref, eDistRef); % [nA x nE]
            else
                RE_r = zeros(size(rfyld_r));
            end

            % --- Interpolate (A,E) -> (Angles,Energies) on output grid ---
            % Note: RustBCA arrays are [nA x nE] => griddedInterpolant({A,E}, V)
            F_sp = griddedInterpolant({A_r, E_r}, spyld_r, 'linear', 'none');
            F_rf = griddedInterpolant({A_r, E_r}, rfyld_r, 'linear', 'none');
            F_re = griddedInterpolant({A_r, E_r}, RE_r,    'linear', 'none');

            [EE, AA] = meshgrid(Energies, Angles);  % AA: angle grid, EE: energy grid

            Yout  = F_sp(AA, EE);
            RNout = F_rf(AA, EE);
            REout = F_re(AA, EE);

            % replace NaNs (outside convex hull) with 0
            Yout(isnan(Yout))  = 0;
            RNout(isnan(RNout))= 0;
            REout(isnan(REout))= 0;

            % clip to [0,1] to be safe (especially reflection)
            RNout = max(min(RNout,1),0);
            REout = max(min(REout,1),0);

            % --- Store into arrays: [proj x targ x angle x energy] ---
            Physical_Sputtering(ip,it,:,:)     = Yout;
            Reflection_Yield(ip,it,:,:)        = RNout;
            Energy_Reflection_Yield(ip,it,:,:) = REout;

        end
    end

    %% ---- Quick range checks --------------------------------------------
    fprintf('Sputtering range: [%.4e, %.4e]\n', min(Physical_Sputtering(:)), max(Physical_Sputtering(:)));
    fprintf('RN range:         [%.4e, %.4e]\n', min(Reflection_Yield(:)),    max(Reflection_Yield(:)));
    fprintf('RE range:         [%.4e, %.4e]\n', min(Energy_Reflection_Yield(:)), max(Energy_Reflection_Yield(:)));

    %% ---- Write NetCDF (reversed dims for C++ reader) -------------------
    if exist(outfile,'file'), delete(outfile); end

    % coordinate variables
    nccreate(outfile,'Angles',      'Dimensions',{'angle_dim', angle_dim}, 'Datatype','double');
    nccreate(outfile,'Energies',    'Dimensions',{'energy_dim', energy_dim}, 'Datatype','double');
    nccreate(outfile,'Projectiles', 'Dimensions',{'projectile_dim', projectile_dim}, 'Datatype','single');
    nccreate(outfile,'Targets',     'Dimensions',{'target_dim', target_dim}, 'Datatype','single');

    % response variables stored as (energy, angle, target, projectile)
    nccreate(outfile,'Physical_Sputtering', ...
        'Dimensions',{'energy_dim', energy_dim, 'angle_dim', angle_dim, 'target_dim', target_dim, 'projectile_dim', projectile_dim}, ...
        'Datatype','double');

    nccreate(outfile,'Reflection_Yield', ...
        'Dimensions',{'energy_dim', energy_dim, 'angle_dim', angle_dim, 'target_dim', target_dim, 'projectile_dim', projectile_dim}, ...
        'Datatype','double');

    nccreate(outfile,'Energy_Reflection_Yield', ...
        'Dimensions',{'energy_dim', energy_dim, 'angle_dim', angle_dim, 'target_dim', target_dim, 'projectile_dim', projectile_dim}, ...
        'Datatype','double');

    % write coordinate vars
    ncwrite(outfile,'Angles',Angles);
    ncwrite(outfile,'Energies',Energies);
    ncwrite(outfile,'Projectiles',Projectiles);
    ncwrite(outfile,'Targets',Targets);

    % permute MATLAB logical order (P,T,A,E) -> (E,A,T,P)
    ncwrite(outfile,'Physical_Sputtering',     permute(Physical_Sputtering,     [4 3 2 1]));
    ncwrite(outfile,'Reflection_Yield',        permute(Reflection_Yield,        [4 3 2 1]));
    ncwrite(outfile,'Energy_Reflection_Yield', permute(Energy_Reflection_Yield, [4 3 2 1]));

    fprintf('Wrote RustBCA-based surface model NetCDF: %s\n', outfile);

    %% ---- Simple visualization (optional quick check) -------------------
    % Show RN and RE vs angle for Ta->W at a representative energy index
    ip_example = 1; it_example = 2;  % Ta -> W example
    [~, iErep] = min(abs(Energies - 1000)); % nearest 1 keV
    RN_angle = squeeze(Reflection_Yield(ip_example,it_example,:,iErep));
    RE_angle = squeeze(Energy_Reflection_Yield(ip_example,it_example,:,iErep));

    figure('Name','Quick check: Ta->W RN/RE vs angle (E~1keV)','Color','w');
    plot(Angles, RN_angle, 'LineWidth', 1.8); hold on;
    plot(Angles, RE_angle, '--', 'LineWidth', 1.8);
    xlabel('\theta (deg)'); ylabel('Yield');
    legend({'RN','RE'},'Location','best'); grid on;
    title(sprintf('Ta \\rightarrow W (E \\approx %.0f eV)', Energies(iErep)));

    %% ==== W -> W comparison & plotting (only) ===========================
    rust_wonw = fullfile(rustbca_dir, 'RustBCA_WonW.nc');
    if exist(rust_wonw,'file') ~= 2
        fprintf('[WARN] RustBCA_WonW.nc not found at %s — skipping W->W comparison plots.\n', rustbca_dir);
        return
    end

    fprintf('Performing W->W comparison using: %s\n', rust_wonw);

    % Read W->W RustBCA
    [A_r, E_r, spyld_r, rfyld_r] = read_rustbca_pair(rust_wonw);

    % Interpolate RustBCA -> surface model grids
    [Spy_r_onSM, Rfy_r_onSM] = interp_rust_to_sm(A_r, E_r, spyld_r, rfyld_r, Angles, Energies);

    % Extract W->W from generated arrays (indices: ip=2, it=2)
    ipW = 2; itW = 2;
    Spy_sm = squeeze(Physical_Sputtering(ipW,itW,:,:));   % [angle x energy]
    Rfy_sm = squeeze(Reflection_Yield(ipW,itW,:,:));
    Efy_sm = squeeze(Energy_Reflection_Yield(ipW,itW,:,:));

    % Plot maps: Sputtering (SM, Rust, Diff)
    doLogMap_sp = true;
    doLogX = true;

    figure('Name','W->W Sputtering: SM | Rust | Diff','Color','w','Units','normalized','Position',[0.05 0.1 0.9 0.7]);
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    nexttile; plot_map_safe(Energies, Angles, Spy_sm, doLogMap_sp, doLogX); title('SurfaceModel: Physical\_Sputtering'); ylabel('Angle (deg)');
    nexttile; plot_map_safe(Energies, Angles, Spy_r_onSM, doLogMap_sp, doLogX); title('RustBCA (interp): spyld'); ylabel('Angle (deg)');
    nexttile; plot_map_diff(Energies, Angles, Spy_sm - Spy_r_onSM, doLogX); title('DIFF: SM - RustBCA'); ylabel('Angle (deg)');
    sgtitle('W \rightarrow W: Sputtering (maps)','FontWeight','bold');

    % Plot maps: Reflection (SM, Rust, Diff)
    figure('Name','W->W Reflection: SM | Rust | Diff','Color','w','Units','normalized','Position',[0.05 0.1 0.9 0.7]);
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    nexttile; plot_map_safe(Energies, Angles, Rfy_sm, false, doLogX); title('SurfaceModel: Reflection\_Yield'); ylabel('Angle (deg)');
    nexttile; plot_map_safe(Energies, Angles, Rfy_r_onSM, false, doLogX); title('RustBCA (interp): rfyld'); ylabel('Angle (deg)');
    nexttile; plot_map_diff(Energies, Angles, Rfy_sm - Rfy_r_onSM, doLogX); title('DIFF: SM - RustBCA'); ylabel('Angle (deg)');
    sgtitle('W \rightarrow W: Reflection (maps)','FontWeight','bold');

    % Linecuts vs Energy at a few angles
    angles_to_plot_deg  = [0 15 30 45 60 75];
    aIdx = nearest_indices(Angles, angles_to_plot_deg);

    figure('Name','W->W Linecuts vs Energy','Color','w','Units','normalized','Position',[0.08 0.1 0.85 0.75]);
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    for k=1:numel(aIdx)
        plot(Energies, Spy_sm(aIdx(k),:), 'LineWidth', 1.6);
    end
    set(gca,'XScale','log'); set(gca,'YScale','log'); grid on; box on;
    xlabel('Energy (eV)'); ylabel('Sputtering yield'); title('SM: Y vs Energy at selected angles');
    legend(compose('A=%.1f°', Angles(aIdx)), 'Location','best');

    nexttile; hold on;
    for k=1:numel(aIdx)
        plot(Energies, Spy_r_onSM(aIdx(k),:), 'LineWidth', 1.6);
    end
    set(gca,'XScale','log'); set(gca,'YScale','log'); grid on; box on;
    xlabel('Energy (eV)'); ylabel('Sputtering yield'); title('RustBCA (interp): Y vs Energy at selected angles');
    legend(compose('A=%.1f°', Angles(aIdx)), 'Location','best');

    % Linecuts vs Angle at selected energies
    energies_to_plot_eV = [20 100 500 1000 5000];
    eIdx = nearest_indices(Energies, energies_to_plot_eV);

    figure('Name','W->W Linecuts vs Angle','Color','w','Units','normalized','Position',[0.08 0.1 0.85 0.75]);
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    for k=1:numel(eIdx)
        plot(Angles, Spy_sm(:,eIdx(k)), 'LineWidth', 1.6);
    end
    set(gca,'YScale','log'); grid on; box on;
    xlabel('Angle (deg)'); ylabel('Sputtering yield'); title('SM: Y vs Angle at selected energies');
    legend(compose('E=%.0f eV', Energies(eIdx)), 'Location','best');

    nexttile; hold on;
    for k=1:numel(eIdx)
        plot(Angles, Spy_r_onSM(:,eIdx(k)), 'LineWidth', 1.6);
    end
    set(gca,'YScale','log'); grid on; box on;
    xlabel('Angle (deg)'); ylabel('Sputtering yield'); title('RustBCA (interp): Y vs Angle');
    legend(compose('E=%.0f eV', Energies(eIdx)), 'Location','best');

    % Print error stats for sputtering & reflection
    fprintf('\n=== Error stats for W->W ===\n');
    print_stats('Sputtering', Spy_sm, Spy_r_onSM);
    print_stats('Reflection', Rfy_sm, Rfy_r_onSM);

    fprintf('W->W comparison plots complete.\n');

end

% ======================================================================
% Compute energy reflection yield from RustBCA reflected energy distribution
% ======================================================================
function RE = compute_energy_reflection(E_r, A_r, RN, egrid_ref, eDistRef)
% RE(angle,energy) = RN(angle,energy) * <E_reflected>/E_incident
%
% Inputs:
%   E_r        [nE x 1] incident energies grid
%   A_r        [nA x 1] angles grid
%   RN         [nA x nE] particle reflection yield
%   egrid_ref  [nBins x 1] reflected-energy bin centers/grid (eV)
%   eDistRef   [nBins x nA x nE] distribution of reflected energies
%
% Output:
%   RE         [nA x nE] energy reflection yield

    nE = numel(E_r);
    nA = numel(A_r);

    if ~isequal(size(RN), [nA nE])
        error('RN must be [nA x nE], got %s', mat2str(size(RN)));
    end
    if size(eDistRef,2) ~= nA || size(eDistRef,3) ~= nE
        error('energyDistRef must be [nBins x nA x nE], got %s', mat2str(size(eDistRef)));
    end

    RE = zeros(nA, nE);

    % ensure column
    egrid_ref = egrid_ref(:);

    for ia = 1:nA
        for ie = 1:nE
            p = eDistRef(:, ia, ie);
            s = sum(p);

            if s > 0
                p = p ./ s; % normalize if not already
                meanEref = sum(p .* egrid_ref);
            else
                meanEref = 0;
            end

            Einc = E_r(ie);
            if Einc > 0
                RE(ia, ie) = RN(ia, ie) * (meanEref / Einc);
            else
                RE(ia, ie) = 0;
            end
        end
    end

    % clip to [0,1]
    RE = max(min(RE,1),0);

end

% ======================================================================
% Helper: read a RustBCA pair file (assumes spyld/rfyld shapes)
% ======================================================================
function [A_r, E_r, spyld_r, rfyld_r] = read_rustbca_pair(rustfile)
    E_r = ncread(rustfile,'E'); E_r = E_r(:);
    A_r = ncread(rustfile,'A'); A_r = A_r(:);
    spyld_r = ncread(rustfile,'spyld'); % [nA x nE]
    rfyld_r = ncread(rustfile,'rfyld'); % [nA x nE]

    if ~isequal(size(spyld_r), [numel(A_r) numel(E_r)])
        error('RustBCA spyld shape mismatch in %s: got %s expected [%d %d]', ...
            rustfile, mat2str(size(spyld_r)), numel(A_r), numel(E_r));
    end
end

% ======================================================================
% Helper: interpolate RustBCA (A_r,E_r) -> (A_sm,E_sm)
% ======================================================================
function [Spy_onSM, Rfy_onSM] = interp_rust_to_sm(A_r, E_r, spyld_r, rfyld_r, A_sm, E_sm)
    F_sp = griddedInterpolant({A_r, E_r}, spyld_r, 'linear', 'none');
    F_rf = griddedInterpolant({A_r, E_r}, rfyld_r, 'linear', 'none');

    [EE, AA] = meshgrid(E_sm, A_sm); % AA: angle rows, EE: energy cols

    Spy_onSM = F_sp(AA, EE);  Spy_onSM(isnan(Spy_onSM)) = 0;
    Rfy_onSM = F_rf(AA, EE);  Rfy_onSM(isnan(Rfy_onSM)) = 0;

    Rfy_onSM = max(min(Rfy_onSM,1),0);
    Spy_onSM = max(Spy_onSM,0);
end

% ======================================================================
% Safe plotting helpers & diagnostics
% ======================================================================
function plot_map_safe(E, A, Z, doLogColor, doLogX)
    h = pcolor(E, A, Z); h.EdgeColor = 'none';
    set(gca,'YDir','normal');
    if doLogX, set(gca,'XScale','log'); end

    if doLogColor
        zpos = Z(isfinite(Z) & Z>0);
        if isempty(zpos)
            set(gca,'ColorScale','linear');
            caxis([0 1]);
        else
            set(gca,'ColorScale','log');
            lo = max(min(zpos), 1e-30);
            hi = max(zpos);
            if hi <= lo, hi = lo*10; end
            caxis([lo hi]);
        end
    end
    colorbar; axis tight; box on;
end

function plot_map_diff(E, A, D, doLogX)
    h = pcolor(E, A, D); h.EdgeColor = 'none';
    set(gca,'YDir','normal');
    if doLogX, set(gca,'XScale','log'); end
    colorbar; axis tight; box on;
end

function idx = nearest_indices(grid, values)
    idx = zeros(size(values));
    for k=1:numel(values)
        [~, idx(k)] = min(abs(grid - values(k)));
    end
    idx = unique(idx);
end

function print_stats(name, SM, RB)
    M = SM(:); R = RB(:);
    ok = isfinite(M) & isfinite(R);
    M = M(ok); R = R(ok);
    D = M - R;

    Linf = max(abs(D));
    L2   = sqrt(mean(D.^2));
    denom = sqrt(mean(R.^2));
    relL2 = L2 / max(denom, eps);

    fprintf('%s: Linf=%.3e, L2=%.3e, relL2=%.3e\n', name, Linf, L2, relL2);
end