function oxygen_coronal_balance()
    % Oxygen coronal equilibrium from ADAS acd89/ccd89 files
    % Ensure 'acd89_o.dat' and 'ccd89_o.dat' are in your MATLAB path

    % Electron temperature [eV]
    Te_target = 8;

    % Load electron temperature grid (log10 scale)
    Te_log = linspace(log10(0.1), log10(1e4), 101);  % log10 Te [eV]
    Te_eV = 10.^Te_log;

    % Read ADAS ionization and recombination rates
    S = read_adas_rates('acd89_o.dat');     % Ionization from O^z to O^{z+1}
    alpha = read_adas_rates('ccd89_o.dat'); % Recombination from O^{z+1} to O^z

    % Find nearest Te index
    [~, idx] = min(abs(Te_eV - Te_target));
    disp(['Using Te = ', num2str(Te_eV(idx)), ' eV']);

    % Extract ionization and recombination rates at Te
    S_te = S(:, idx);        % Length 8 (O0 to O7)
    alpha_te = alpha(:, idx); % Length 8 (O1+ to O8+)

    % Compute coronal equilibrium
    Z = 8;  % Max charge state for oxygen
    Nz = zeros(1, Z+1);  % Populations from O0 to O8+
    Nz(1) = 1;  % Arbitrary starting population

    for z = 1:Z
        Nz(z+1) = Nz(z) * S_te(z) / alpha_te(z);
    end

    % Normalize populations
    Nz = Nz / sum(Nz);

    % Plot
    figure;
    bar(0:Z, Nz, 'FaceColor', [0.2 0.6 0.8]);
    xlabel('Charge state (O^{z+})');
    ylabel('Fractional abundance');
    title(['Coronal Equilibrium for Oxygen at T_e = ', num2str(Te_eV(idx)), ' eV']);
    grid on;
end

function rates = read_adas_rates(filename)
    % Read ADAS acd/ccd file format for oxygen (Z=8, 101 points)
    fid = fopen(filename, 'r');
    if fid == -1
        error(['Failed to open ', filename]);
    end

    Z = 8;       % Oxygen charge states: O0 to O8+
    Nt = 101;    % Number of temperature points
    rates = zeros(Z, Nt);

    % Skip 3 header lines
    for i = 1:3
        fgetl(fid);
    end

    % Read each charge state block
    for z = 1:Z
        fgetl(fid);  % Skip label line

        count = 0;
        while count < Nt
            line = fgetl(fid);
            if ~ischar(line), break; end
            values = sscanf(line, '%e');
            n_values = length(values);
            rates(z, count+1:count+n_values) = values';
            count = count + n_values;
        end
    end

    fclose(fid);
end