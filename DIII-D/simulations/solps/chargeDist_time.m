%% Initialization
close all; clear; clc;

% Constants
Zmax = 10;          % Maximum charge state for Neon (Ne0 to Ne10+)
ne = 1e19;          % Electron density (m^-3)
Te_target = 100;    % Electron temperature in eV
tspan = [0, 1e-3];  % Time range (0 to 1 ms)

% Load ADAS Ionization (SCD) and Recombination (ACD) Coefficients
[Te_acd, ACD] = read_ADAS_acd('acd85_ne.dat');  % Recombination rates
[Te_scd, SCD] = read_ADAS_acd('acd89_ne.dat');  % Ionization rates

% Interpolate rates for the target temperature
R = interp1(Te_acd, ACD, Te_target, 'linear', 'extrap'); % Recombination
S = interp1(Te_scd, SCD, Te_target, 'linear', 'extrap'); % Ionization

% Initial condition: all neutral Neon atoms
nZ0 = zeros(Zmax+1, 1);
nZ0(1) = 1;  % Assume all starts as Ne0

% Solve the ODE system for charge state evolution
[t, nZ] = ode45(@(t, nZ) ionization_ode(t, nZ, S, R, ne, Zmax), tspan, nZ0);

%% Plot results
figure;
plot(t * 1e3, nZ, 'LineWidth', 2);
legend(arrayfun(@(z) sprintf('Ne^{%d+}', z), 0:Zmax, 'UniformOutput', false));
xlabel('Time (ms)');
ylabel('Fractional abundance');
title('Time-Dependent Ionization Balance of Neon');
grid on;

%% ODE Function: Ionization Balance
function dndt = ionization_ode(~, nZ, S, R, ne, Zmax)
    dndt = zeros(Zmax+1, 1);
    
    for Z = 0:Zmax
        % Ionization to next state
        ionization = 0;
        if Z < Zmax
            ionization = ne * S(Z+1) * nZ(Z+1);
        end
        
        % Recombination from higher state
        recombination = 0;
        if Z > 0
            recombination = ne * R(Z) * nZ(Z+1);
        end
        
        % Evolution equation
        if Z == 0
            dndt(Z+1) = -ionization + recombination;
        elseif Z == Zmax
            dndt(Z+1) = -recombination + ne * S(Z) * nZ(Z);
        else
            dndt(Z+1) = -ionization + recombination + ne * S(Z) * nZ(Z) - ne * R(Z+1) * nZ(Z+1);
        end
    end
end

%% Function to Read ADAS ACD/SCD Data
function [Te_vals, coeff_vals] = read_ADAS_acd(filename)
    % Open the file
    fid = fopen(filename, 'r');
    
    % Read header lines (adjust as needed)
    for i = 1:5
        header_line = fgetl(fid);
        disp(header_line);  % Display header info
    end

    % Read numeric data (assuming Te in first column, coefficients in second)
    data = fscanf(fid, '%f', [inf]); % Read all data
    fclose(fid);

    % Split into temperature and coefficient arrays
    num_Te = length(data) / 2; % Assuming two columns: Te and coefficients
    Te_vals = data(1:num_Te);
    coeff_vals = data(num_Te+1:end);
end