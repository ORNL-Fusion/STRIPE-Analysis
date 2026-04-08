function out = calculate_W_ionization_mfp(Te_eV, ne_m3, neutralEnergy_eV, chargeState, filename)
% calculate_W_ionization_mfp
% Compute W ionization lifetime and mean free path from ADAS_Rates_W.nc.
%
% out = calculate_W_ionization_mfp(Te_eV, ne_m3, neutralEnergy_eV, chargeState, filename)
%
% Inputs:
%   Te_eV            Electron temperature [eV]
%   ne_m3            Electron density [m^-3]
%   neutralEnergy_eV Neutral tungsten kinetic energy [eV] (used for speed)
%   chargeState      Initial W charge state z for z -> z+1 (default 0)
%   filename         NetCDF file (default 'ADAS_Rates_W.nc')
%
% Output struct fields:
%   .rateCoeff_m3s   Interpolated ionization rate coefficient <sigma v> [m^3/s]
%   .tauIon_s        Ionization lifetime 1/(ne*<sigma v>) [s]
%   .vNeutral_ms     Neutral speed from kinetic energy [m/s]
%   .lambdaIon_m     Ionization mean free path v*tau [m]
%   .transition      [z z+1] charge-state transition used

if nargin < 5 || isempty(filename)
    filename = 'ADAS_Rates_W.nc';
end
if nargin < 4 || isempty(chargeState)
    chargeState = 0;
end
if nargin < 3 || isempty(neutralEnergy_eV)
    neutralEnergy_eV = 3.0;
end

% ADAS grids are stored as log10 values in this project.
Te_grid = 10.^ncread(filename, 'gridTemperature_Ionization');
ne_grid = 10.^ncread(filename, 'gridDensity_Ionization');

% ncread returns IonizationRateCoeff with dimensions [density, temperature, charge]
S_all = 10.^ncread(filename, 'IonizationRateCoeff');
chargePairs = ncread(filename, 'gridChargeState_Ionization');

idx = find(chargePairs(1,:) == chargeState & chargePairs(2,:) == (chargeState + 1), 1, 'first');
if isempty(idx)
    error('Transition W^{%d+} -> W^{%d+} not found in %s.', chargeState, chargeState + 1, filename);
end

S_slice = S_all(:,:,idx); % [density, temperature]
S_interp = interpn(ne_grid, Te_grid, S_slice, ne_m3, Te_eV, 'linear', 0);

if S_interp <= 0
    error(['Interpolated ionization rate is zero or outside the ADAS grid. ', ...
           'Check Te/ne range or use nearest-neighbor extrapolation.']);
end

tau_ion = 1.0 / (ne_m3 * S_interp);

amu = 1.66053906660e-27;      % [kg]
eV_to_J = 1.602176634e-19;    % [J/eV]
mW = 183.84 * amu;            % Tungsten atomic mass [kg]
vW = sqrt(2.0 * neutralEnergy_eV * eV_to_J / mW);
lambda_ion = vW * tau_ion;

out.rateCoeff_m3s = S_interp;
out.tauIon_s = tau_ion;
out.vNeutral_ms = vW;
out.lambdaIon_m = lambda_ion;
out.transition = [chargeState, chargeState + 1];
out.Te_eV = Te_eV;
out.ne_m3 = ne_m3;
out.neutralEnergy_eV = neutralEnergy_eV;

fprintf('W^{%d+} -> W^{%d+}\\n', chargeState, chargeState + 1);
fprintf('Te = %.6g eV, ne = %.6g m^-3, E_W = %.6g eV\\n', Te_eV, ne_m3, neutralEnergy_eV);
fprintf('<sigma v> = %.6e m^3/s\\n', S_interp);
fprintf('tau_ion   = %.6e s\\n', tau_ion);
fprintf('v_W       = %.6e m/s\\n', vW);
fprintf('lambda_ion= %.6e m\\n', lambda_ion);

end
