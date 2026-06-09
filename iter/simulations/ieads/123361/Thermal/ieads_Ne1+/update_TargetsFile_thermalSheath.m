%% Replace first column with thermal sheath voltage for D+ ions

clear; clc;

infile  = 'Targets_Ne1+.txt';
outfile = 'Targets_Ne1+_thermal_sheath.txt';

data = readmatrix(infile);

% Columns
% 1 = old sheath voltage
% 3 = Te [eV]
% 4 = Ti [eV]

Te = data(:,3);
Ti = data(:,4);

% D+ ion mass in electron-mass-normalized units
me = 1/2000;
mD = 2;   % D+ mass number

% Avoid bad values
Te(~isfinite(Te) | Te <= 0) = 0;
Ti(~isfinite(Ti) | Ti <  0) = 0;

thermal_sheath_Dplus = zeros(size(Te));

good = Te > 0;

thermal_sheath_Dplus(good) = Te(good) .* abs( ...
    0.5 .* log( ...
    (2*pi*me/mD) .* (1 + Ti(good)./Te(good)) ...
    ) ...
    );

% Replace first column
data(:,1) = thermal_sheath_Dplus;

% Save revised target file
writematrix(data, outfile, 'Delimiter', ' ');

fprintf('Wrote revised file: %s\n', outfile);
fprintf('Min thermal sheath voltage = %.6e V\n', min(thermal_sheath_Dplus));
fprintf('Max thermal sheath voltage = %.6e V\n', max(thermal_sheath_Dplus));
fprintf('Mean thermal sheath voltage = %.6e V\n', mean(thermal_sheath_Dplus));