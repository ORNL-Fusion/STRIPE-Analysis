function compare_surfaceModel_vs_rustBCA(surfaceModelFile, rustbca_dir)
% compare_surfaceModel_vs_rustBCA
%
% Plots ALL sputtering + reflection yields from:
%   (1) GITR surface-model NetCDF: surface_model_*.nc
%   (2) RustBCA tables: RustBCA_<Proj>on<Targ>.nc
%
% and does "check" plots:
%   - side-by-side maps (SurfaceModel vs RustBCA) for each pair
%   - difference maps (SurfaceModel - RustBCA)
%   - linecuts vs Energy at selected Angles
%   - linecuts vs Angle at selected Energies
%   - basic error stats (Linf, L2, relL2)
%
% Assumptions:
%   Surface model variables per your ncdisp:
%     Angles [nA], Energies [nE], Projectiles [nP], Targets [nT]
%     Physical_Sputtering, Reflection_Yield, Energy_Reflection_Yield
%     with dims (energy, angle, target, projectile) in-file
%   RustBCA files provide:
%     E [nE_r], A [nA_r], spyld [nA_r x nE_r], rfyld [nA_r x nE_r]
%
% Usage:
%   compare_surfaceModel_vs_rustBCA( ...
%     'surface_model_GITRm_rustbca_C_W_d.nc', ...
%     '/Users/78k/ORNL Dropbox/Atul Kumar/work/STRIPE-Analysis/rustBCA_data');
%
% Tip:
%   If your RustBCA filenames differ, edit rustFileForPair() near bottom.

    if nargin < 1 || isempty(surfaceModelFile)
        error('Provide surfaceModelFile');
    end
    if nargin < 2 || isempty(rustbca_dir)
        rustbca_dir = '.';
    end

    %% ---------------- Read surface model ----------------
    Angles_sm    = ncread(surfaceModelFile,'Angles');    Angles_sm = Angles_sm(:);
    Energies_sm  = ncread(surfaceModelFile,'Energies');  Energies_sm = Energies_sm(:);
    Proj_sm      = ncread(surfaceModelFile,'Projectiles'); Proj_sm = double(Proj_sm(:));
    Targ_sm      = ncread(surfaceModelFile,'Targets');     Targ_sm = double(Targ_sm(:));

    PS_sm_in  = ncread(surfaceModelFile,'Physical_Sputtering');
    RN_sm_in  = ncread(surfaceModelFile,'Reflection_Yield');
    RE_sm_in  = ncread(surfaceModelFile,'Energy_Reflection_Yield');

    % In-file dims are (energy, angle, target, projectile)
    % Convert to MATLAB logical order: (projectile, target, angle, energy)
    PS_sm = permute(PS_sm_in, [4 3 2 1]);
    RN_sm = permute(RN_sm_in, [4 3 2 1]);
    RE_sm = permute(RE_sm_in, [4 3 2 1]);

    nP = numel(Proj_sm);
    nT = numel(Targ_sm);
    nA = numel(Angles_sm);
    nE = numel(Energies_sm);

    fprintf('Surface model: nP=%d nT=%d nA=%d nE=%d\n', nP,nT,nA,nE);

    %% ---------------- Plot settings ----------------
    doLogMap = true;          % log colors for sputtering maps (often helpful)
    doLogX   = true;          % log x-axis for energy
    angles_to_plot_deg  = [0 15 30 45 60 75];    % linecuts vs energy
    energies_to_plot_eV = [20 100 500 1000 5000];% linecuts vs angle

    %% ---------------- Loop all pairs ----------------
    for ip = 1:nP
        for it = 1:nT

            Zp = Proj_sm(ip);
            Zt = Targ_sm(it);

            % --- Load RustBCA for this pair (if exists) ---
            rustfile = rustFileForPair(rustbca_dir, Zp, Zt);
            hasRust = exist(rustfile,'file') == 2;

            if hasRust
                [A_r, E_r, spyld_r, rfyld_r] = read_rustbca_pair(rustfile);
                % Interpolate RustBCA -> surface-model grids (Angles_sm,Energies_sm)
                [Spy_r_onSM, Rfy_r_onSM] = interp_rust_to_sm(A_r, E_r, spyld_r, rfyld_r, Angles_sm, Energies_sm);
            else
                fprintf('[WARN] No RustBCA file for Z=%g -> Z=%g (skipping Rust compare)\n', Zp, Zt);
                Spy_r_onSM = nan(nA,nE);
                Rfy_r_onSM = nan(nA,nE);
            end

            % --- Extract SM matrices for this pair on (angle x energy) ---
            % SM stored as (proj,target,angle,energy)
            Spy_sm = squeeze(PS_sm(ip,it,:,:));  % [nA x nE]
            Rfy_sm = squeeze(RN_sm(ip,it,:,:));  % [nA x nE]
            Efy_sm = squeeze(RE_sm(ip,it,:,:));  % [nA x nE]

            % --- CHECK: plot maps (SM vs Rust) + diff for sputtering + reflection ---
            pairTitle = sprintf('Zp=%g -> Zt=%g (ip=%d,it=%d)', Zp, Zt, ip, it);

            % ===== SPUTTERING =====
            figure('Name',['Sputtering map ' pairTitle],'Color','w','Units','normalized','Position',[0.05 0.08 0.9 0.75]);
            tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

            nexttile; plot_map(Energies_sm, Angles_sm, Spy_sm, doLogMap, doLogX);
            title('SurfaceModel: Physical\_Sputtering'); ylabel('Angle (deg)');

            nexttile; 
            if hasRust
                plot_map(Energies_sm, Angles_sm, Spy_r_onSM, doLogMap, doLogX);
                title('RustBCA (interp): spyld'); 
            else
                axis off; text(0.1,0.5,'No RustBCA file','FontSize',14);
            end

            nexttile;
            if hasRust
                D = Spy_sm - Spy_r_onSM;
                plot_map_diff(Energies_sm, Angles_sm, D, doLogX);
                title('DIFF: SM - RustBCA'); 
            else
                axis off;
            end

            % ===== REFLECTION =====
            nexttile; plot_map(Energies_sm, Angles_sm, Rfy_sm, false, doLogX);
            title('SurfaceModel: Reflection\_Yield'); xlabel('Energy (eV)'); ylabel('Angle (deg)');

            nexttile;
            if hasRust
                plot_map(Energies_sm, Angles_sm, Rfy_r_onSM, false, doLogX);
                title('RustBCA (interp): rfyld'); xlabel('Energy (eV)');
            else
                axis off;
            end

            nexttile;
            if hasRust
                D = Rfy_sm - Rfy_r_onSM;
                plot_map_diff(Energies_sm, Angles_sm, D, doLogX);
                title('DIFF: SM - RustBCA'); xlabel('Energy (eV)');
            else
                axis off;
            end

            sgtitle(['MAP CHECK: ' pairTitle],'FontWeight','bold');

            % --- Error stats printed to console (only if Rust exists) ---
            if hasRust
                fprintf('\n=== CHECK STATS for %s ===\n', pairTitle);
                print_stats('Sputtering', Spy_sm, Spy_r_onSM);
                print_stats('Reflection', Rfy_sm, Rfy_r_onSM);
            end

            % --- Linecuts vs Energy at selected angles (SP + RF) ---
            figure('Name',['Linecuts vs Energy ' pairTitle],'Color','w','Units','normalized','Position',[0.08 0.1 0.85 0.75]);
            tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

            % Sputtering: vs E at fixed angles
            nexttile; hold on;
            aIdx = nearest_indices(Angles_sm, angles_to_plot_deg);
            for k=1:numel(aIdx)
                plot(Energies_sm, Spy_sm(aIdx(k),:), 'LineWidth', 1.6);
            end
            if doLogX, set(gca,'XScale','log'); end
            if doLogMap, set(gca,'YScale','log'); end
            grid on; box on;
            xlabel('Energy (eV)'); ylabel('Yield');
            title('Sputtering (SM): Y vs Energy');
            legend(compose('A=%.1f°', Angles_sm(aIdx)), 'Location','best');

            nexttile; hold on;
            if hasRust
                for k=1:numel(aIdx)
                    plot(Energies_sm, Spy_r_onSM(aIdx(k),:), 'LineWidth', 1.6);
                end
                if doLogX, set(gca,'XScale','log'); end
                if doLogMap, set(gca,'YScale','log'); end
                grid on; box on;
                xlabel('Energy (eV)'); ylabel('Yield');
                title('Sputtering (Rust interp): Y vs Energy');
                legend(compose('A=%.1f°', Angles_sm(aIdx)), 'Location','best');
            else
                axis off; text(0.1,0.5,'No RustBCA file','FontSize',14);
            end

            % Reflection: vs E at fixed angles
            nexttile; hold on;
            for k=1:numel(aIdx)
                plot(Energies_sm, Rfy_sm(aIdx(k),:), 'LineWidth', 1.6);
            end
            if doLogX, set(gca,'XScale','log'); end
            grid on; box on;
            xlabel('Energy (eV)'); ylabel('Yield');
            title('Reflection (SM): RN vs Energy');
            legend(compose('A=%.1f°', Angles_sm(aIdx)), 'Location','best');

            nexttile; hold on;
            if hasRust
                for k=1:numel(aIdx)
                    plot(Energies_sm, Rfy_r_onSM(aIdx(k),:), 'LineWidth', 1.6);
                end
                if doLogX, set(gca,'XScale','log'); end
                grid on; box on;
                xlabel('Energy (eV)'); ylabel('Yield');
                title('Reflection (Rust interp): rfyld vs Energy');
                legend(compose('A=%.1f°', Angles_sm(aIdx)), 'Location','best');
            else
                axis off;
            end

            sgtitle(['LINECUTS vs ENERGY: ' pairTitle],'FontWeight','bold');

            % --- Linecuts vs Angle at selected energies (SP + RF) ---
            figure('Name',['Linecuts vs Angle ' pairTitle],'Color','w','Units','normalized','Position',[0.08 0.1 0.85 0.75]);
            tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

            eIdx = nearest_indices(Energies_sm, energies_to_plot_eV);

            nexttile; hold on;
            for k=1:numel(eIdx)
                plot(Angles_sm, Spy_sm(:,eIdx(k)), 'LineWidth', 1.6);
            end
            if doLogMap, set(gca,'YScale','log'); end
            grid on; box on;
            xlabel('Angle (deg)'); ylabel('Yield');
            title('Sputtering (SM): Y vs Angle');
            legend(compose('E=%.0f eV', Energies_sm(eIdx)), 'Location','best');

            nexttile; hold on;
            if hasRust
                for k=1:numel(eIdx)
                    plot(Angles_sm, Spy_r_onSM(:,eIdx(k)), 'LineWidth', 1.6);
                end
                if doLogMap, set(gca,'YScale','log'); end
                grid on; box on;
                xlabel('Angle (deg)'); ylabel('Yield');
                title('Sputtering (Rust interp): Y vs Angle');
                legend(compose('E=%.0f eV', Energies_sm(eIdx)), 'Location','best');
            else
                axis off;
            end

            nexttile; hold on;
            for k=1:numel(eIdx)
                plot(Angles_sm, Rfy_sm(:,eIdx(k)), 'LineWidth', 1.6);
            end
            grid on; box on;
            xlabel('Angle (deg)'); ylabel('Yield');
            title('Reflection (SM): RN vs Angle');
            legend(compose('E=%.0f eV', Energies_sm(eIdx)), 'Location','best');

            nexttile; hold on;
            if hasRust
                for k=1:numel(eIdx)
                    plot(Angles_sm, Rfy_r_onSM(:,eIdx(k)), 'LineWidth', 1.6);
                end
                grid on; box on;
                xlabel('Angle (deg)'); ylabel('Yield');
                title('Reflection (Rust interp): rfyld vs Angle');
                legend(compose('E=%.0f eV', Energies_sm(eIdx)), 'Location','best');
            else
                axis off;
            end

            sgtitle(['LINECUTS vs ANGLE: ' pairTitle],'FontWeight','bold');

            % --- Optional: Energy reflection yield map (SM only) ---
            figure('Name',['Energy Reflection (SM) ' pairTitle],'Color','w');
            plot_map(Energies_sm, Angles_sm, Efy_sm, false, doLogX);
            title(['SurfaceModel: Energy\_Reflection\_Yield — ' pairTitle]);
            xlabel('Energy (eV)'); ylabel('Angle (deg)');

        end
    end

    fprintf('\nDone. Inspect figures for all pairs; console prints stats where RustBCA exists.\n');
end

%% ===================== Helpers =====================

function plot_map(E, A, Z, doLogColor, doLogX)
    h = pcolor(E, A, Z); 
    h.EdgeColor = 'none';
    set(gca,'YDir','normal');
    if doLogX, set(gca,'XScale','log'); end

    if doLogColor
        % Only use log color if there is at least one positive value
        zpos = Z(isfinite(Z) & Z>0);
        if isempty(zpos)
            % All zeros / nonpositive -> use linear scale, fixed limits
            set(gca,'ColorScale','linear');
            caxis([0 1]);
        else
            set(gca,'ColorScale','log');
            lo = max(min(zpos), 1e-30);
            hi = max(zpos);
            if hi <= lo
                hi = lo*10; % ensure increasing
            end
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

function [Spy_onSM, Rfy_onSM] = interp_rust_to_sm(A_r, E_r, spyld_r, rfyld_r, A_sm, E_sm)
    % Interpolate RustBCA (A_r,E_r) -> (A_sm,E_sm)
    F_sp = griddedInterpolant({A_r, E_r}, spyld_r, 'linear', 'none');
    F_rf = griddedInterpolant({A_r, E_r}, rfyld_r, 'linear', 'none');

    [EE, AA] = meshgrid(E_sm, A_sm); % AA: angle, EE: energy

    Spy_onSM = F_sp(AA, EE);  Spy_onSM(isnan(Spy_onSM)) = 0;
    Rfy_onSM = F_rf(AA, EE);  Rfy_onSM(isnan(Rfy_onSM)) = 0;

    Rfy_onSM = max(min(Rfy_onSM,1),0);
    Spy_onSM = max(Spy_onSM,0);
end

function rustfile = rustFileForPair(rustbca_dir, Zp, Zt)
% Edit this if your RustBCA filenames are different.
% Default naming expected: RustBCA_<ProjName>on<TargName>.nc
% where names are: Ta, W, C, etc.

    nameP = Z_to_name(Zp);
    nameT = Z_to_name(Zt);

    rustfile = fullfile(rustbca_dir, sprintf('RustBCA_%son%s.nc', nameP, nameT));
end

function name = Z_to_name(Z)
% minimal mapping; extend as needed
    switch round(Z)
        case 6,  name = 'C';
        case 73, name = 'Ta';
        case 74, name = 'W';
        case 1,  name = 'H';
        case 2,  name = 'He';
        case 8,  name = 'O';
        otherwise
            name = sprintf('Z%d', round(Z)); % fallback
    end
end