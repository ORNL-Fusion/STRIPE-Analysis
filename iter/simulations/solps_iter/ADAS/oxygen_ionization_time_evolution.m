function oxygen_ionization_time_evolution()
    % -------------------------------
    % Parameters
    % -------------------------------
    Te_target = 50;        % Electron temperature in eV
    ne = 1e19;             % Electron density [m^-3]
    tspan = [0, 1e-3];     % Time range in seconds

    % -------------------------------
    % Step 1: Detect actual Nt and temperature grid
    % -------------------------------
    S_temp = read_adas_rates('acd89_o.dat', 200);  % Try reading up to 200 columns
    Nt = size(S_temp, 2);                          % Actual number of Te points
    Te_log = linspace(-1, 4, Nt);                  % Assumes log10(Te/eV) grid
    Te_eV = 10.^Te_log;

    % -------------------------------
    % Step 2: Read rates using detected Nt
    % -------------------------------
    S_all = read_adas_rates('acd89_o.dat', Nt);      % Ionization rates
    alpha_all = read_adas_rates('ccd89_o.dat', Nt);  % Recombination rates

    % -------------------------------
    % Step 3: Interpolate rates at target Te
    % -------------------------------
    Z = 8;
    S = zeros(Z, 1);
    alpha = zeros(Z, 1);
    for i = 1:Z
        S(i) = interp1(Te_eV, S_all(i,:), Te_target, 'pchip');
        alpha(i) = interp1(Te_eV, alpha_all(i,:), Te_target, 'pchip');
    end

    % -------------------------------
    % Step 4: Initial condition — fully stripped (O8+)
    % -------------------------------
    N0 = zeros(Z+1, 1);
    N0(end) = 1;  % All population in O8+

    % -------------------------------
    % Step 5: Solve time-dependent ionization balance
    % -------------------------------
    [t, Nsol] = ode15s(@(t, N) ion_balance_rhs(t, N, ne, S, alpha), tspan, N0);

    % -------------------------------
    % Step 6: Plot results
    % -------------------------------
    figure;
    plot(t * 1e3, Nsol, 'LineWidth', 1.5);  % Time in ms
    xlabel('Time [ms]');
    ylabel('Fractional abundance');
    legend(arrayfun(@(z) sprintf('O^{%d+}', z), 0:Z, 'UniformOutput', false));
    title(['Oxygen Ionization Evolution at T_e = ', num2str(Te_target), ...
           ' eV, n_e = ', sprintf('%.1e', ne), ' m^{-3}']);
    grid on;
end

function dNdt = ion_balance_rhs(~, N, ne, S, alpha)
    Z = 8;
    dNdt = zeros(Z+1, 1);

    for z = 0:Z
        gain = 0;
        loss = 0;

        if z > 0
            gain = gain + ne * S(z) * N(z);         % From z-1
        end
        if z < Z
            gain = gain + ne * alpha(z+1) * N(z+2); % From z+1
        end
        if z < Z
            loss = loss + ne * S(z+1);              % To z+1
        end
        if z > 0
            loss = loss + ne * alpha(z);            % To z-1
        end
        dNdt(z+1) = gain - loss * N(z+1);
    end
end

function rates = read_adas_rates(filename, Nt_max)
    % Reads ADAS ACD/CCD file, returns [Z × Nt] rate matrix
    fid = fopen(filename, 'r');
    if fid == -1
        error(['Could not open ', filename]);
    end

    Z = 8;                     % Oxygen charge states (O0 to O7+)
    rates = zeros(Z, Nt_max);  % Over-allocate, trim later

    for i = 1:3
        fgetl(fid);            % Skip headers
    end

    for z = 1:Z
        fgetl(fid);            % Skip label
        count = 0;
        while count < Nt_max
            line = fgetl(fid);
            if ~ischar(line), break; end
            values = sscanf(line, '%e');
            rates(z, count + (1:length(values))) = values';
            count = count + length(values);
        end
    end

    % Trim to actual size by removing trailing zeros
    nonzero_cols = find(any(rates, 1));
    rates = rates(:, 1:max(nonzero_cols));

    fclose(fid);
end