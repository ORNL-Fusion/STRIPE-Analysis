function neon_coronal_balance()
    % Neon coronal equilibrium + transport correction (0D loss-time)
    clearvars -except varargin
    close all
    clc

    Te_target = 10;         % eV
    use_transport = true;
    tau_s = 1e-4;          % s (Inf -> coronal)

    [Te_s_log, ne_s_log, scd_log10_cgs] = ADF11s('scd89_ne.dat');
    [Te_a_log, ne_a_log, acd_log10_cgs] = ADF11a('acd89_ne.dat');

    Te_s = 10.^Te_s_log(:);
    Te_a = 10.^Te_a_log(:);

    iNe_s = 1; iNe_a = 1;
    [~, iTe_s] = min(abs(Te_s - Te_target));
    [~, iTe_a] = min(abs(Te_a - Te_target));

    Te_used = Te_s(iTe_s);
    ne_used_m3 = (10.^ne_s_log(iNe_s)) * 1e6; % cm^-3 -> m^-3
    fprintf('Using Te = %.6g eV, ne = %.6e m^-3 (lowest ADAS density)\n', Te_used, ne_used_m3);

    % Extract and convert to SI (m^3/s)
    S_te     = (10.^squeeze(scd_log10_cgs(iNe_s, iTe_s, :))) / 1e6;
    alpha_te = (10.^squeeze(acd_log10_cgs(iNe_a, iTe_a, :))) / 1e6;

    fprintf('size(scd_log10_cgs) = [%d %d %d]\n', size(scd_log10_cgs));
    fprintf('size(acd_log10_cgs) = [%d %d %d]\n', size(acd_log10_cgs));
    fprintf('numel(S_te) = %d, numel(alpha_te) = %d\n', numel(S_te), numel(alpha_te));

    % --- AUTO-DETECT Z from ADAS transitions ---
    Z = min(numel(S_te), numel(alpha_te));   % number of transitions (0->1 ... (Z-1)->Z)
    nStates = Z + 1;                         % populations (0..Z)

    fprintf('Auto-using Z = %d (states 0..%d)\n', Z, Z);

    % tau vector
    if isscalar(tau_s)
        tau_vec = tau_s * ones(nStates,1);
    else
        tau_vec = tau_s(:);
        if numel(tau_vec) ~= nStates
            error('tau_s must be scalar or vector of length %d.', nStates);
        end
    end

    % --------------------------
    % Coronal ladder
    % --------------------------
    Nz_cor = zeros(nStates,1);
    Nz_cor(1) = 1;
    for z = 1:Z
        Nz_cor(z+1) = Nz_cor(z) * S_te(z) / max(alpha_te(z), realmin);
    end
    Nz_cor = Nz_cor / sum(Nz_cor);

    % --------------------------
    % Transport-corrected solve
    % --------------------------
    if ~use_transport
        Nz = Nz_cor;
    else
        A = zeros(nStates, nStates);

        % Interior z = 1..Z-1  (indices 2..Z)
        for j = 2:Z
            z = j-1;
            A(j, j-1) = ne_used_m3 * S_te(z);          % (z-1)->z
            A(j, j+1) = ne_used_m3 * alpha_te(z+1);    % (z+1)->z
            A(j, j)   = -ne_used_m3*(S_te(z+1) + alpha_te(z)) ...
                        - 1/max(tau_vec(j), realmin);
        end

        % Boundary z=0
        A(1,1) = -ne_used_m3*S_te(1) - 1/max(tau_vec(1), realmin);
        A(1,2) =  ne_used_m3*alpha_te(1);   % 1->0 uses alpha(1)

        % Boundary z=Z
        A(nStates, nStates-1) = ne_used_m3*S_te(Z);    % (Z-1)->Z
        A(nStates, nStates)   = -ne_used_m3*alpha_te(Z) ...
                                - 1/max(tau_vec(nStates), realmin);  % Z->Z-1 uses alpha(Z)

        % Normalization constraint
        Aeq = A;
        beq = zeros(nStates,1);
        Aeq(1,:) = 1;
        beq(1) = 1;

        Nz = Aeq \ beq;
        Nz = abs(Nz);
        Nz = Nz / sum(Nz);
    end

    % --------------------------
    % Plot
    % --------------------------
    figure('Color','w'); hold on
    bar(0:Z, Nz_cor, 0.45);
    bar((0:Z)+0.45, Nz, 0.45);
    xlabel('Charge state (Ne^{z+})');
    ylabel('Fractional abundance');
    title(sprintf('Neon charge balance at T_e = %.4g eV (ne=%.2e m^{-3})', Te_used, ne_used_m3));
    grid on
    legend({'Coronal (ladder)','With transport loss'}, 'Location','best')
    xlim([-0.5, Z+1.0])

    if use_transport
        fprintf('Transport correction: tau = %.3e s\n', tau_s);
    end
end