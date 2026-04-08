%% Carbon ADAS rates -> NetCDF + ionization mean free path example
% Self-contained reader for ADF11-style files (SCD/ACD/CCD).
% Stores log10(SI) rate coefficients on [density, temperature, charge] grids.
%
% IMPORTANT:
%   SCD, ACD, and CCD files can have different density/temperature grids.
%   This script preserves separate grids for each dataset in the NetCDF file.

close all;
clear;
clc;

%% -------------------- constants --------------------
ME   = 9.10938356e-31; %#ok<NASGU>
MI   = 1.6737236e-27;
Q    = 1.60217662e-19;
EPS0 = 8.854187e-12; %#ok<NASGU>
Z    = 6; % Carbon atomic number

%% -------------------- user controls --------------------
adas_dir = pwd;

% ADAS files
file_inz  = fullfile(adas_dir, 'scd89_c.dat');
file_rcmb = fullfile(adas_dir, 'acd85_c.dat');
file_cx   = fullfile(adas_dir, 'ccd89_c.dat');

out_nc = 'ADAS_Rates_C.nc';

% CCD order in ADF11 is typically cm^6/s for 3-body/CX-like tables
ccd_power_cm = 6; % cm^6/s -> m^6/s conversion => log10 shift -6

% Mean free path example inputs
impurity_mass_amu          = 12;   % Carbon mass [amu]
impurity_charge            = 0;    % C^0 -> C^1+
impurity_kinetic_energy_eV = 7.41;
electron_temperature_eV    = 10;
electron_density_m3        = 1e17;

% Gyroradius example
B_field_T                  = 1.4;    % [T]
impurity_Tperp_eV          = 3;    % [eV]
impurity_charge_state_gyro = 1;    % must be > 0

%% -------------------- read ADAS blocks --------------------
assert(isfile(file_inz),  'Missing file: %s', file_inz);
assert(isfile(file_rcmb), 'Missing file: %s', file_rcmb);

scd = read_adf11_blocks(file_inz);
acd = read_adf11_blocks(file_rcmb);

has_cx = isfile(file_cx);
if has_cx
    ccd = read_adf11_blocks(file_cx);
else
    warning('CCD file not found: %s (continuing without CX table)', file_cx);
    ccd = [];
end

[IonizationChargeState,   IonizationRateCoeff_cgs_log]    = blocks_to_cube(scd);
[RecombinationChargeState, RecombinationRateCoeff_cgs_log] = blocks_to_cube(acd);

if has_cx
    [CXChargeState, CXRateCoeff_cgs_log] = blocks_to_cube(ccd);
else
    CXChargeState = zeros(0,2);
    CXRateCoeff_cgs_log = zeros(0,0,0);
end

%% -------------------- grids for each dataset --------------------
% ADF11 grids are usually log10(ne[cm^-3]) and log10(Te[eV])

DensityIon_cgs_log = scd.logNe(:).';
TempIon_log        = scd.logTe(:).';

DensityRec_cgs_log = acd.logNe(:).';
TempRec_log        = acd.logTe(:).';

DensityIon_si_log = DensityIon_cgs_log + 6; % cm^-3 -> m^-3 in log10
DensityRec_si_log = DensityRec_cgs_log + 6;

if has_cx
    DensityCX_cgs_log = ccd.logNe(:).';
    TempCX_log        = ccd.logTe(:).';
    DensityCX_si_log  = DensityCX_cgs_log + 6;
end

%% -------------------- convert rate coefficients to SI --------------------
% SCD/ACD: cm^3/s -> m^3/s => log10 shift -6
IonizationRateCoeff_si_log    = IonizationRateCoeff_cgs_log - 6;
RecombinationRateCoeff_si_log = RecombinationRateCoeff_cgs_log - 6;

% CCD: cm^N/s -> m^N/s => log10 shift -N
if has_cx
    CXRateCoeff_si_log = CXRateCoeff_cgs_log - ccd_power_cm;
else
    CXRateCoeff_si_log = zeros(0,0,0);
end

% Reorder from [temperature, density, charge] -> [density, temperature, charge]
IonizationRateCoeff_si_log    = permute(IonizationRateCoeff_si_log,    [2 1 3]);
RecombinationRateCoeff_si_log = permute(RecombinationRateCoeff_si_log, [2 1 3]);

if has_cx
    CXRateCoeff_si_log = permute(CXRateCoeff_si_log, [2 1 3]);
end

%% -------------------- quick consistency checks --------------------
assert(size(IonizationRateCoeff_si_log,1) == numel(DensityIon_si_log), ...
    'Ionization density dimension mismatch.');
assert(size(IonizationRateCoeff_si_log,2) == numel(TempIon_log), ...
    'Ionization temperature dimension mismatch.');

assert(size(RecombinationRateCoeff_si_log,1) == numel(DensityRec_si_log), ...
    'Recombination density dimension mismatch.');
assert(size(RecombinationRateCoeff_si_log,2) == numel(TempRec_log), ...
    'Recombination temperature dimension mismatch.');

if has_cx
    assert(size(CXRateCoeff_si_log,1) == numel(DensityCX_si_log), ...
        'CX density dimension mismatch.');
    assert(size(CXRateCoeff_si_log,2) == numel(TempCX_log), ...
        'CX temperature dimension mismatch.');
end

fprintf('Loaded SCD grid: ne=%d, Te=%d, charge states=%d\n', ...
    numel(DensityIon_si_log), numel(TempIon_log), size(IonizationChargeState,1));
fprintf('Loaded ACD grid: ne=%d, Te=%d, charge states=%d\n', ...
    numel(DensityRec_si_log), numel(TempRec_log), size(RecombinationChargeState,1));
if has_cx
    fprintf('Loaded CCD grid: ne=%d, Te=%d, charge states=%d\n', ...
        numel(DensityCX_si_log), numel(TempCX_log), size(CXChargeState,1));
end

%% -------------------- write NetCDF --------------------
if exist(out_nc, 'file')
    delete(out_nc);
end

ncid = netcdf.create(out_nc, 'NC_CLOBBER');

% Dimensions
dimScalar = netcdf.defDim(ncid, 'scalar', 1);
dimPair   = netcdf.defDim(ncid, 'pair', 2);

dimTemp_Ionize       = netcdf.defDim(ncid, 'n_Temperatures_Ionize',    numel(TempIon_log));
dimDensity_Ionize    = netcdf.defDim(ncid, 'n_Densities_Ionize',       numel(DensityIon_si_log));
dimChargeState_Ionize = netcdf.defDim(ncid, 'n_ChargeStates_Ionize',   size(IonizationChargeState,1));

dimTemp_Recombine       = netcdf.defDim(ncid, 'n_Temperatures_Recombine',  numel(TempRec_log));
dimDensity_Recombine    = netcdf.defDim(ncid, 'n_Densities_Recombine',     numel(DensityRec_si_log));
dimChargeState_Recombine = netcdf.defDim(ncid, 'n_ChargeStates_Recombine', size(RecombinationChargeState,1));

if has_cx
    dimTemp_CX        = netcdf.defDim(ncid, 'n_Temperatures_CX',   numel(TempCX_log));
    dimDensity_CX     = netcdf.defDim(ncid, 'n_Densities_CX',      numel(DensityCX_si_log));
    dimChargeState_CX = netcdf.defDim(ncid, 'n_ChargeStates_CX',   size(CXChargeState,1));
end

% Scalars and grids
Z_ID = netcdf.defVar(ncid, 'Atomic_Number', 'int', dimScalar);

TempGridIonization       = netcdf.defVar(ncid, 'gridTemperature_Ionization',    'double', dimTemp_Ionize);
DensityGridIonization    = netcdf.defVar(ncid, 'gridDensity_Ionization',        'double', dimDensity_Ionize);
ChargeStateGridIonization = netcdf.defVar(ncid, 'gridChargeState_Ionization',   'double', [dimChargeState_Ionize dimPair]);

TempGridRecombination       = netcdf.defVar(ncid, 'gridTemperature_Recombination',    'double', dimTemp_Recombine);
DensityGridRecombination    = netcdf.defVar(ncid, 'gridDensity_Recombination',        'double', dimDensity_Recombine);
ChargeStateGridRecombination = netcdf.defVar(ncid, 'gridChargeState_Recombination',   'double', [dimChargeState_Recombine dimPair]);

IonizeCoeff = netcdf.defVar(ncid, 'IonizationRateCoeff_log10_SI', 'double', ...
    [dimDensity_Ionize dimTemp_Ionize dimChargeState_Ionize]);

RecombineCoeff = netcdf.defVar(ncid, 'RecombinationRateCoeff_log10_SI', 'double', ...
    [dimDensity_Recombine dimTemp_Recombine dimChargeState_Recombine]);

if has_cx
    TempGridCX        = netcdf.defVar(ncid, 'gridTemperature_CX',  'double', dimTemp_CX);
    DensityGridCX     = netcdf.defVar(ncid, 'gridDensity_CX',      'double', dimDensity_CX);
    ChargeStateGridCX = netcdf.defVar(ncid, 'gridChargeState_CX',  'double', [dimChargeState_CX dimPair]);

    CXCoeff = netcdf.defVar(ncid, 'ChargeExchangeRateCoeff_log10_SI', 'double', ...
        [dimDensity_CX dimTemp_CX dimChargeState_CX]);
end

netcdf.endDef(ncid);

% Write metadata and grids
netcdf.putVar(ncid, Z_ID, Z);

netcdf.putVar(ncid, TempGridIonization, TempIon_log);
netcdf.putVar(ncid, DensityGridIonization, DensityIon_si_log);
netcdf.putVar(ncid, ChargeStateGridIonization, IonizationChargeState);

netcdf.putVar(ncid, TempGridRecombination, TempRec_log);
netcdf.putVar(ncid, DensityGridRecombination, DensityRec_si_log);
netcdf.putVar(ncid, ChargeStateGridRecombination, RecombinationChargeState);

if has_cx
    netcdf.putVar(ncid, TempGridCX, TempCX_log);
    netcdf.putVar(ncid, DensityGridCX, DensityCX_si_log);
    netcdf.putVar(ncid, ChargeStateGridCX, CXChargeState);
end

% Write coefficients
netcdf.putVar(ncid, IonizeCoeff, IonizationRateCoeff_si_log);
netcdf.putVar(ncid, RecombineCoeff, RecombinationRateCoeff_si_log);

if has_cx
    netcdf.putVar(ncid, CXCoeff, CXRateCoeff_si_log);
end

netcdf.close(ncid);
fprintf('Wrote %s\n', out_nc);

%% -------------------- interpolation example --------------------
q_idx = impurity_charge + 1;
assert(q_idx >= 1 && q_idx <= size(IonizationRateCoeff_si_log,3), ...
    'impurity_charge=%d is out of bounds for available SCD states (%d).', ...
    impurity_charge, size(IonizationRateCoeff_si_log,3));

Coeff_log10_SI = interpn( ...
    DensityIon_si_log, ...
    TempIon_log, ...
    IonizationRateCoeff_si_log(:,:,q_idx), ...
    log10(electron_density_m3), ...
    log10(electron_temperature_eV), ...
    'linear', NaN);

if ~isfinite(Coeff_log10_SI)
    error(['Requested (ne,Te) point is outside the ionization grid range. ' ...
           'Try values inside the SCD table bounds.']);
end

ionization_time_s = 1 / (10^Coeff_log10_SI * electron_density_m3);
impurity_velocity = sqrt(2 * impurity_kinetic_energy_eV * Q / (impurity_mass_amu * MI));
mean_free_path_m  = ionization_time_s * impurity_velocity;

assert(impurity_charge_state_gyro > 0, ...
    'impurity_charge_state_gyro must be > 0 for gyroradius calculation.');

v_perp_imp = sqrt(2 * impurity_Tperp_eV * Q / (impurity_mass_amu * MI));
gyro_radius_imp_m = (impurity_mass_amu * MI * v_perp_imp) / ...
                    (impurity_charge_state_gyro * Q * B_field_T);

fprintf('\nC^%d at Te=%.3g eV, ne=%.3e m^-3\n', ...
    impurity_charge, electron_temperature_eV, electron_density_m3);
fprintf('log10(k_ion [m^3/s]) = %.6f\n', Coeff_log10_SI);
fprintf('ionization_time      = %.6e s\n', ionization_time_s);
fprintf('mean_free_path       = %.6e m\n', mean_free_path_m);
fprintf('C gyroradius: T_perp=%.3g eV, Z=%d, B=%.3g T -> rho=%.6e m\n', ...
    impurity_Tperp_eV, impurity_charge_state_gyro, B_field_T, gyro_radius_imp_m);

%% -------------------- quick plots --------------------
% Temperature cut at closest density index
[~, iNe] = min(abs(DensityIon_si_log - log10(electron_density_m3)));
figure('Color', 'w');
plot(10.^TempIon_log, squeeze(IonizationRateCoeff_si_log(iNe,:,q_idx)), 'LineWidth', 1.8);
set(gca, 'XScale', 'log');
grid on;
box on;
xlabel('T_e [eV]');
ylabel('log_{10}(k_{ion} [m^3/s])');
title(sprintf('C^%d ionization rate vs T_e at n_e=%.2e m^{-3}', ...
    impurity_charge, electron_density_m3));

% Density cut at closest temperature index
[~, iTe] = min(abs(TempIon_log - log10(electron_temperature_eV)));
figure('Color', 'w');
plot(10.^DensityIon_si_log, squeeze(IonizationRateCoeff_si_log(:,iTe,q_idx)), 'LineWidth', 1.8);
set(gca, 'XScale', 'log');
grid on;
box on;
xlabel('n_e [m^{-3}]');
ylabel('log_{10}(k_{ion} [m^3/s])');
title(sprintf('C^%d ionization rate vs n_e at T_e=%.2f eV', ...
    impurity_charge, electron_temperature_eV));

%% -------------------- local helpers --------------------
function [chargePairs, cube] = blocks_to_cube(adf)
    idx = find(~cellfun(@isempty, adf.blocks));
    ncs = numel(idx);
    if ncs == 0
        error('No Z1 blocks found in ADAS file.');
    end

    nTe = numel(adf.logTe);
    nNe = numel(adf.logNe);

    cube = zeros(nTe, nNe, ncs);
    chargePairs = zeros(ncs, 2);

    for k = 1:ncs
        z1 = idx(k);
        cube(:,:,k) = adf.blocks{z1};
        chargePairs(k,:) = [z1-1, z1];
    end
end

function adf = read_adf11_blocks(filename)
    txt = readlines(filename);
    txt = cellstr(txt);
    txt = txt(:);

    nums = sscanf(txt{1}, '%f');
    if numel(nums) < 5
        error('Cannot parse ADF11 header in %s', filename);
    end

    nBlocks = round(nums(1));
    nNe = nums(2);
    nTe = nums(3);

    i = 2;
    while i <= numel(txt) && (isempty(strtrim(txt{i})) || contains(txt{i}, '-'))
        i = i + 1;
    end

    logNe = [];
    while numel(logNe) < nNe
        if i > numel(txt)
            error('Unexpected EOF while reading logNe in %s', filename);
        end
        vals = sscanf(txt{i}, '%f').';
        logNe = [logNe vals]; %#ok<AGROW>
        i = i + 1;
    end
    logNe = logNe(1:nNe);

    logTe = [];
    while numel(logTe) < nTe
        if i > numel(txt)
            error('Unexpected EOF while reading logTe in %s', filename);
        end
        vals = sscanf(txt{i}, '%f').';
        logTe = [logTe vals]; %#ok<AGROW>
        i = i + 1;
    end
    logTe = logTe(1:nTe);

    blocks = cell(1, max(80, nBlocks));
    nLoaded = 0;

    while i <= numel(txt) && nLoaded < nBlocks
        ln = txt{i};

        if ~contains(ln, 'Z1=')
            i = i + 1;
            continue;
        end

        z1 = sscanf(extractAfter(ln, 'Z1='), '%f', 1);
        if isempty(z1) || ~isfinite(z1)
            tok = regexp(ln, 'Z1=\s*([0-9]+)', 'tokens', 'once');
            if isempty(tok)
                i = i + 1;
                continue;
            end
            z1 = str2double(tok{1});
        end
        z1 = round(z1);

        i = i + 1;
        M = zeros(nTe, nNe);

        for it = 1:nTe
            row = [];
            while numel(row) < nNe
                if i > numel(txt)
                    error('Unexpected EOF in block Z1=%d of %s', z1, filename);
                end
                if contains(txt{i}, 'Z1=')
                    error('Malformed ADF11 block Z1=%d in %s (early next header).', ...
                        z1, filename);
                end
                vals = sscanf(txt{i}, '%f').';
                if ~isempty(vals)
                    row = [row vals]; %#ok<AGROW>
                end
                i = i + 1;
            end
            M(it,:) = row(1:nNe);
        end

        if z1 > numel(blocks)
            blocks{z1} = [];
        end
        blocks{z1} = M;
        nLoaded = nLoaded + 1;
    end

    adf.logNe = logNe;
    adf.logTe = logTe;
    adf.blocks = blocks;
end