function make_gitr_surface_response_TaW(outfile)
% make_gitr_surface_response_TaW
%
% High-resolution GITR surface-response file for:
%   Projectiles: Ta (Z=73), W (Z=74)
%   Targets:     Ta (Z=73), W (Z=74)
%
% Surface model (from challenge "Surface model" section):
%
%   Y(D -> W) = Y(Ta -> W) = 0
%   Y(D -> Fe) = Y(Ta -> Fe) = 0
%   Y(D -> Ta) = 0.001
%   Y(Ta -> Ta) = 0
%
% For reflection (RN) and energy reflection (RE) for Ta→Ta and Ta→W,
% angle-only dependence (independent of energy), with theta in degrees:
%
%   RN(Ta -> W) = RN(Ta -> Ta)
%              = 0.71 + 0.7 * tanh( 2.09 * (theta*pi/180) - 2.63 )
%
%   RE(Ta -> W) = RE(Ta -> Ta)
%              = 0.93 + 0.93 * tanh( 2.7 * (theta*pi/180) - 4.12 )
%
% RN(Ta -> Fe) = RE(Ta -> Fe) = 0
%
% Here we:
%   - Implement these formulas for projectile Ta on targets Ta and W.
%   - Set RN = RE = 0 for any other projectile–target pair (e.g. W→Ta,W),
%     so that W is present as a projectile species but has no reflection
%     yields specified by the model.
%
% Usage:
%   make_gitr_surface_response_TaW;
%   make_gitr_surface_response_TaW('TaW_surface_response.nc');

    if nargin < 1
        outfile = 'TaW_surface_response.nc';
    end

    %% --- High-resolution grids ------------------------------------------
    angle_dim      = 361;       % 0–89 deg with ~0.25° spacing
    energy_dim     = 200;       % smooth log-scale energy resolution
    projectile_dim = 2;         % Ta, W
    target_dim     = 2;         % Ta, W

    Angles   = linspace(0.0, 89.0, angle_dim);                 % degrees
    Energies = logspace(log10(10), log10(1e4), energy_dim);    % eV

    % Species labels (atomic numbers)
    Z_Ta = 73;
    Z_W  = 74;
    Projectiles = single([Z_Ta, Z_W]);           % index 1=Ta, 2=W
    Targets     = single([Z_Ta, Z_W]);           % index 1=Ta, 2=W

    %% --- Allocate response arrays ---------------------------------------
    Physical_Sputtering      = zeros(projectile_dim, target_dim, ...
                                     angle_dim, energy_dim);
    Reflection_Yield         = zeros(projectile_dim, target_dim, ...
                                     angle_dim, energy_dim);
    Energy_Reflection_Yield  = zeros(projectile_dim, target_dim, ...
                                     angle_dim, energy_dim);
    % NOTE: all sputtering is zero for these heavy-on-heavy pairs,
    % consistent with the simplified surface model.

    %% --- Fill reflection yields (angle-dependent only) -------------------
    for ia = 1:angle_dim
        theta = Angles(ia);  % degrees

        for ip = 1:projectile_dim
            for it = 1:target_dim
                Zp = Projectiles(ip);
                Zt = Targets(it);

                [RN, RE] = reflection_surface_model(Zp, Zt, theta);

                % Replicate along energy dimension (independent of E):
                Reflection_Yield(ip, it, ia, :)        = RN;
                Energy_Reflection_Yield(ip, it, ia, :) = RE;
            end
        end
    end

    %% --- Quick consistency checks ----------------------------------------
    RN_min = min(Reflection_Yield(:));
    RN_max = max(Reflection_Yield(:));
    RE_min = min(Energy_Reflection_Yield(:));
    RE_max = max(Energy_Reflection_Yield(:));

    fprintf('RN range: [%.4f, %.4f]\n', RN_min, RN_max);
    fprintf('RE range: [%.4f, %.4f]\n', RE_min, RE_max);

    if RN_min < -1e-6 || RN_max > 1.1 || RE_min < -1e-6 || RE_max > 1.1
        warning('Reflection yields outside expected [0,1] range. Check formulas.');
    end

   %% --- Create NetCDF file with reversed dimension ordering for C++ ---------
if exist(outfile,'file'), delete(outfile); end

% Coordinate variables
nccreate(outfile,'Angles', ...
         'Dimensions',{'angle_dim', angle_dim}, 'Datatype','double');
nccreate(outfile,'Energies', ...
         'Dimensions',{'energy_dim', energy_dim}, 'Datatype','double');
nccreate(outfile,'Projectiles', ...
         'Dimensions',{'projectile_dim', projectile_dim}, 'Datatype','single');
nccreate(outfile,'Targets', ...
         'Dimensions',{'target_dim', target_dim}, 'Datatype','single');

% *** IMPORTANT ***
% Because the C++ reader REVERSES dimension order,
% we must write the variables in REVERSE order:
%
%     (energy_dim, angle_dim, target_dim, projectile_dim)
%
% so that C++ interprets them as:
%
%     (projectile_dim, target_dim, angle_dim, energy_dim)

nccreate(outfile,'Physical_Sputtering', ...
    'Dimensions',{'energy_dim',     energy_dim, ...
                  'angle_dim',      angle_dim, ...
                  'target_dim',     target_dim, ...
                  'projectile_dim', projectile_dim}, ...
    'Datatype','double');

nccreate(outfile,'Reflection_Yield', ...
    'Dimensions',{'energy_dim',     energy_dim, ...
                  'angle_dim',      angle_dim, ...
                  'target_dim',     target_dim, ...
                  'projectile_dim', projectile_dim}, ...
    'Datatype','double');

nccreate(outfile,'Energy_Reflection_Yield', ...
    'Dimensions',{'energy_dim',     energy_dim, ...
                  'angle_dim',      angle_dim, ...
                  'target_dim',     target_dim, ...
                  'projectile_dim', projectile_dim}, ...
    'Datatype','double');

% Write data with reversed dimension order expected by C++
ncwrite(outfile,'Angles',Angles);
ncwrite(outfile,'Energies',Energies);
ncwrite(outfile,'Projectiles',Projectiles);
ncwrite(outfile,'Targets',Targets);

% Permute data BEFORE writing
ncwrite(outfile,'Physical_Sputtering', ...
    permute(Physical_Sputtering, [4 3 2 1]));    % (E,A,T,P)

ncwrite(outfile,'Reflection_Yield', ...
    permute(Reflection_Yield, [4 3 2 1]));

ncwrite(outfile,'Energy_Reflection_Yield', ...
    permute(Energy_Reflection_Yield, [4 3 2 1]));

fprintf('Wrote reversed-dimension NetCDF file for C++ reader: %s\n', outfile);

    %% --- Visualization: 1D RN, RE vs angle for all pairs -----------------
    projNames = {'Ta', 'W'};
    targNames = {'Ta', 'W'};

    figure('Name','Reflection vs Angle (RN & RE, all pairs)','Color','w');
    tiledlayout(projectile_dim, target_dim, 'Padding','compact', 'TileSpacing','compact');

    for ip = 1:projectile_dim
        for it = 1:target_dim
            nexttile;
            RN_angle = squeeze(Reflection_Yield(ip, it, :, 1));
            RE_angle = squeeze(Energy_Reflection_Yield(ip, it, :, 1));

            plot(Angles, RN_angle, 'LineWidth', 1.8); hold on;
            plot(Angles, RE_angle, '--', 'LineWidth', 1.8);
            xlabel('\theta (deg)');
            ylabel('Yield');
            title(sprintf('%s \\rightarrow %s', projNames{ip}, targNames{it}));
            legend({'RN', 'RE'}, 'Location', 'best');
            grid on;
        end
    end

    %% --- Additional cross-check plot similar to Fig. 3 -------------------
    ip_Ta = 1; it_Ta = 1;  % Ta -> Ta
    RN_TaTa = squeeze(Reflection_Yield(ip_Ta, it_Ta, :, 1));
    RE_TaTa = squeeze(Energy_Reflection_Yield(ip_Ta, it_Ta, :, 1));

    figure('Name','Cross-check: Ta \rightarrow Ta reflection (Fig.3 style)','Color','w');
    hold on;
    mask_RN = RN_TaTa > 0;
    mask_RE = RE_TaTa > 0;

    semilogy(Angles(mask_RN), RN_TaTa(mask_RN), 'LineWidth', 2.0);
    semilogy(Angles(mask_RE), RE_TaTa(mask_RE), '--', 'LineWidth', 2.0);

    xlabel('Angle from normal (deg)');
    ylabel('Reflection yield');
    legend({'R_N (Ta \rightarrow Ta)', 'R_E (Ta \rightarrow Ta)'}, 'Location', 'best');
    grid on;
    title('Angular dependence of Ta reflection yields (model)');
    ylim([1e-3 2]);

    %% --- 2D maps vs angle & energy for each pair -------------------------
    for ip = 1:projectile_dim
        for it = 1:target_dim
            Y2  = squeeze(Physical_Sputtering(ip, it, :, :));
            RN2 = squeeze(Reflection_Yield(ip, it, :, :));
            RE2 = squeeze(Energy_Reflection_Yield(ip, it, :, :));

            % Sputtering (all zeros here, but plotted for completeness)
            figure('Name',sprintf('Physical Sputtering %s->%s', ...
                                  projNames{ip}, targNames{it}), ...
                   'Color','w');
            imagesc(Energies, Angles, Y2);
            set(gca, 'YDir','normal', 'XScale','log');
            xlabel('Energy (eV)');
            ylabel('\theta (deg)');
            title(sprintf('Physical Sputtering Y: %s \\rightarrow %s', ...
                          projNames{ip}, targNames{it}));
            colorbar;

            % RN
            figure('Name',sprintf('Reflection Yield RN %s->%s', ...
                                  projNames{ip}, targNames{it}), ...
                   'Color','w');
            imagesc(Energies, Angles, RN2);
            set(gca, 'YDir','normal', 'XScale','log');
            xlabel('Energy (eV)');
            ylabel('\theta (deg)');
            title(sprintf('Reflection Yield R_N: %s \\rightarrow %s', ...
                          projNames{ip}, targNames{it}));
            colorbar;

            % RE
            figure('Name',sprintf('Energy Reflection Yield RE %s->%s', ...
                                  projNames{ip}, targNames{it}), ...
                   'Color','w');
            imagesc(Energies, Angles, RE2);
            set(gca, 'YDir','normal', 'XScale','log');
            xlabel('Energy (eV)');
            ylabel('\theta (deg)');
            title(sprintf('Energy Reflection Yield R_E: %s \\rightarrow %s', ...
                          projNames{ip}, targNames{it}));
            colorbar;
        end
    end
end

% =====================================================================
% === Local helper implementing the surface model =====================
% =====================================================================

function [RN, RE] = reflection_surface_model(Zp, Zt, theta_deg)
% reflection_surface_model
%
% Implements the prescribed surface model for Ta reflection:
%
%   RN(Ta->W) = RN(Ta->Ta)
%             = 0.71 + 0.7 * tanh( 2.09 * (theta*pi/180) - 2.63 )
%
%   RE(Ta->W) = RE(Ta->Ta)
%             = 0.93 + 0.93 * tanh( 2.7 * (theta*pi/180) - 4.12 )
%
% All other projectile–target combinations in this script are set to 0.
%
% INPUT:
%   Zp        projectile atomic number
%   Zt        target atomic number
%   theta_deg incidence angle from normal [deg]
%
% OUTPUT:
%   RN, RE    reflection and energy reflection yields

    Z_Ta = 73;
    Z_W  = 74;

    theta_rad = theta_deg * pi/180;

    if Zp == Z_Ta && (Zt == Z_Ta || Zt == Z_W)
        % Ta projectile on Ta or W target: use given formulas
        RN = 0.71 + 0.7  * tanh( 2.09 * theta_rad - 2.63 );
        RE = 0.93 + 0.93 * tanh( 2.7  * theta_rad - 4.12 );
    else
        % All other combinations (e.g. W projectile) -> zero reflection
        RN = 0.0;
        RE = 0.0;
    end

    % Clip to [0,1] to be safe numerically
    RN = max(min(RN, 1.0), 0.0);
    RE = max(min(RE, 1.0), 0.0);
end