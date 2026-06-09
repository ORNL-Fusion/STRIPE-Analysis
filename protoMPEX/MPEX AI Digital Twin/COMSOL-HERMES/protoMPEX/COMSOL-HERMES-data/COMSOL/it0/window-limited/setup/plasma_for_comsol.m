%% ProtoMPEX: Compute psi, ne, Te and write COMSOL-like TXT
% Output format:
% % Description:        Poloidal magnetic flux, Magnetic flux density, r-component, Magnetic flux density, phi-component, Magnetic flux density, z-component, Electron density, Electron temperature
% % Length unit:        m
% % r,z,mf.Psi (Wb),mf.Br (T),mf.Bphi (T),mf.Bz (T),ne (m^-3),te (eV)
%
% Column indices:
% 1 -> r
% 2 -> z
% 3 -> psi
% 4 -> Br
% 5 -> Bt (= Bphi)
% 6 -> Bz
% 7 -> ne
% 8 -> te

close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING B-FIELD FILE (assumed [r x z])
% -------------------------------------------------------------------------
fileB = 'bfield_protoMPEX_shotSeries_5.nc';

r  = ncread(fileB,'r');      % [m]
z  = ncread(fileB,'z');      % [m]
Br = ncread(fileB,'br');     % [r x z]
Bt = ncread(fileB,'bt');     % [r x z]
Bz = ncread(fileB,'bz');     % [r x z]

% Force coordinate vectors to be columns
r = r(:);
z = z(:);

% If orientation is swapped, fix once
if ~isequal(size(Bz), [numel(r) numel(z)])
    Br = permute(Br,[2 1]);
    Bt = permute(Bt,[2 1]);
    Bz = permute(Bz,[2 1]);
end

nR = numel(r);
nZ = numel(z);

% -------------------------------------------------------------------------
% PLOT INPUT Bz FIELD
% -------------------------------------------------------------------------
figure;
imagesc(z, r, Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Input $B_z$ [T]','Interpreter','latex');
colorbar;
axis image tight;
xlim([0.5 4.2]);

% -------------------------------------------------------------------------
% CALCULATE PSI USING ORIGINAL METHOD
% -------------------------------------------------------------------------
disp('Calculating psi using original summation method...');

dr = r(2) - r(1);
psi = zeros(nR, nZ);

for ii = 1:nR
    for jj = 1:nZ
        psi(ii,jj) = 2*pi * sum(Bz(1:ii,jj) .* r(1:ii)) * dr;
    end
end

% Normalize psi using helicon / stagnation location
z0 = 1.745;         % [m]
r0 = 0.0625;        % [m]
psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');
psiN = psi ./ psi_norm_val;

figure;
imagesc(z, r, psi);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$\psi$ [Wb]','Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
disp('Constructing Te and ne profiles...');

qe  = 1.602176634e-19;      % [C]

% Te profile
Te0       = 2;      % eV at core (psiN=0)
Te_edge   = 8;      % eV at LCFS (psiN=1)
Te_floor  = 0.2;    % eV far outside
alpha     = 1.75;   % core shaping
lambdaPsi = 0.15;   % SOL falloff width

Te = zeros(size(psiN));

% Core (psiN <= 1)
Te(psiN <= 1) = Te_edge + (Te0 - Te_edge) .* (1 - psiN(psiN <= 1)).^alpha;

% SOL (psiN > 1)
Te(psiN > 1)  = Te_floor + (Te_edge - Te_floor) .* exp(-(psiN(psiN > 1) - 1) ./ lambdaPsi);

% Density model
gasPuff_Gps = 1.0e20;   % [particles/s]
PRF_kW      = 180.0;    % [kW]
Eion_eV     = 1.0e3;    % [eV/pair]
ne_cap      = 1.0e20;   % [m^-3]

Preq_kW = gasPuff_Gps * Eion_eV * qe / 1e3;
powerFactor = min(PRF_kW / max(Preq_kW,1e-12), 1.0);

ne_min = 1.0e17;        % [m^-3]
ne_max = ne_cap;        % [m^-3]
ne_max = min(max(ne_max, 5.0e18), 5.0e20);

ne_ref_for_power = 1.0e20;   % [m^-3]
Preq_ne_kW = Preq_kW * (ne_max / ne_ref_for_power);
powerMargin = PRF_kW / max(Preq_ne_kW,1e-12);

Ne = (psiN < 1) .* (ne_max - ne_min) .* (1 - psiN).^1.75 + ne_min;

% -------------------------------------------------------------------------
% OPTIONAL CHECK PLOTS
% -------------------------------------------------------------------------
figure;
pcolor(r, z, psi');
shading flat;
colorbar;
xlabel('r');
ylabel('z');
title('Poloidal flux \psi');

figure;
pcolor(r, z, Ne');
shading flat;
colorbar;
xlabel('r');
ylabel('z');
title('Electron density n_e');

figure;
pcolor(r, z, Te');
shading flat;
colorbar;
xlabel('r');
ylabel('z');
title('Electron temperature T_e');

% -------------------------------------------------------------------------
% BUILD COORDINATE GRIDS WITH CORRECT LINEAR INDEXING
% Matches MATLAB linearization of [nR x nZ] arrays:
% (r1,z1), (r2,z1), ..., (rN,z1), (r1,z2), ...
% -------------------------------------------------------------------------
[Rgrid, Zgrid] = ndgrid(r, z);

% Flatten variables
r_flat   = Rgrid(:);
z_flat   = Zgrid(:);
psi_flat = psi(:);
br_flat  = Br(:);
bt_flat  = Bt(:);   % COMSOL names this Bphi
bz_flat  = Bz(:);
ne_flat  = Ne(:);
te_flat  = Te(:);

% Combine into one matrix in exact column order
output = [r_flat, z_flat, psi_flat, br_flat, bt_flat, bz_flat, ne_flat, te_flat];

% -------------------------------------------------------------------------
% WRITE COMSOL-LIKE TXT FILE
% -------------------------------------------------------------------------
disp('Writing COMSOL-like TXT file...');

outTxt = 'protoMPEX_profiles.txt';
if exist(outTxt,'file')
    delete(outTxt);
end

fid = fopen(outTxt, 'w');

fprintf(fid, '%% Model:              ai2_2Dcase.mph\n');
fprintf(fid, '%% Version:            COMSOL 6.3.0.420\n');
fprintf(fid, '%% Date:               Feb 25 2026, 16:38\n');
fprintf(fid, '%% Dimension:          2\n');
fprintf(fid, '%% Nodes:              %d\n', nR*nZ);
fprintf(fid, '%% Expressions:        6\n');
fprintf(fid, '%% Description:        Poloidal magnetic flux, Magnetic flux density, r-component, Magnetic flux density, phi-component, Magnetic flux density, z-component, Electron density, Electron temperature\n');
fprintf(fid, '%% Length unit:        m\n');
fprintf(fid, '%% r,z,mf.Psi (Wb),mf.Br (T),mf.Bphi (T),mf.Bz (T),ne (m^-3),te (eV)\n');

fclose(fid);

% Append comma-separated numeric data
dlmwrite(outTxt, output, '-append', 'delimiter', ',', 'precision', '%.15E');

disp(['✅ Wrote ', outTxt, ' successfully.']);

% -------------------------------------------------------------------------
% COLUMN INDEX GUIDE
% -------------------------------------------------------------------------
disp('Column indices in protoMPEX_profiles.txt:');
disp('1 = r');
disp('2 = z');
disp('3 = psi');
disp('4 = Br');
disp('5 = Bphi (= Bt)');
disp('6 = Bz');
disp('7 = ne');
disp('8 = te');

% -------------------------------------------------------------------------
% QUICK PREVIEW
% -------------------------------------------------------------------------
disp('First 5 rows of exported data:');
disp(output(1:min(5,size(output,1)), :));