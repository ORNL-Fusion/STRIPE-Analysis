%% ProtoMPEX: Compute psi, Br, ne, Te and write COMSOL-like TXT
% Bz is forced flat from the beginning.
% psi, Br, ne, Te, and exported fields all use flat Bz = 0.05 T.

close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING B-FIELD FILE FOR GRID AND Bt
% -------------------------------------------------------------------------
fileB = 'bfield_protoMPEX_shotSeries_1.nc';

r  = ncread(fileB,'r');      % [m]
z  = ncread(fileB,'z');      % [m]
Bt = ncread(fileB,'bt');     % [r x z]
Bz_in = ncread(fileB,'bz');  % only for size/orientation check

r = r(:);
z = z(:);

if ~isequal(size(Bz_in), [numel(r) numel(z)])
    Bt    = permute(Bt,[2 1]);
    Bz_in = permute(Bz_in,[2 1]);
end

nR = numel(r);
nZ = numel(z);

% -------------------------------------------------------------------------
% FORCE FLAT Bz FROM THE BEGINNING
% -------------------------------------------------------------------------
Bz0 = 0.05;                  % [T]
Bz  = Bz0 * ones(nR,nZ);     % [T]

disp('Using flat Bz = 0.05 T from the beginning.');

% -------------------------------------------------------------------------
% CALCULATE PSI USING FLAT Bz
% psi(r,z) = 2*pi * integral_0^r Bz(r,z) * r dr
% -------------------------------------------------------------------------
dr = r(2) - r(1);
psi = zeros(nR, nZ);

for ii = 1:nR
    for jj = 1:nZ
        psi(ii,jj) = 2*pi * sum(Bz(1:ii,jj) .* r(1:ii)) * dr;
    end
end

% -------------------------------------------------------------------------
% CALCULATE Br CONSISTENTLY FROM PSI
% Br = -1/(2*pi*r) * dpsi/dz
%
% Since flat Bz is independent of z, psi is independent of z,
% therefore Br should be zero everywhere.
% -------------------------------------------------------------------------
Br = zeros(nR,nZ);

% Optional numerical form, useful if Bz later varies with z:
% dpsidz = zeros(nR,nZ);
% for ii = 1:nR
%     dpsidz(ii,:) = gradient(psi(ii,:), z);
% end
% Br = zeros(nR,nZ);
% for ii = 1:nR
%     if r(ii) > 0
%         Br(ii,:) = -dpsidz(ii,:) ./ (2*pi*r(ii));
%     else
%         Br(ii,:) = 0;
%     end
% end

% -------------------------------------------------------------------------
% NORMALIZE PSI
% -------------------------------------------------------------------------
z0 = 1.745;
r0 = 0.0625;

psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');
psiN = psi ./ psi_norm_val;

% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
figure;
imagesc(z, r, Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Flat $B_z = 0.05$ T','Interpreter','latex');
colorbar;
axis image tight;
xlim([0.5 4.2]);

figure;
imagesc(z, r, Br);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Calculated $B_r$ [T]','Interpreter','latex');
colorbar;
axis image tight;
xlim([0.5 4.2]);

figure;
imagesc(z, r, psi);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$\psi$ from flat $B_z$ [Wb]','Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
disp('Constructing Te and ne profiles...');

qe = 1.602176634e-19;

Te0       = 2;
Te_edge   = 8;
Te_floor  = 0.2;
alpha     = 1.75;
lambdaPsi = 0.15;

Te = zeros(size(psiN));

Te(psiN <= 1) = Te_edge + ...
    (Te0 - Te_edge) .* (1 - psiN(psiN <= 1)).^alpha;

Te(psiN > 1) = Te_floor + ...
    (Te_edge - Te_floor) .* exp(-(psiN(psiN > 1) - 1) ./ lambdaPsi);

gasPuff_Gps = 1.0e20;
PRF_kW      = 180.0;
Eion_eV     = 1.0e3;
ne_cap      = 1.0e20;

Preq_kW = gasPuff_Gps * Eion_eV * qe / 1e3;
powerFactor = min(PRF_kW / max(Preq_kW,1e-12), 1.0);

ne_min = 1.0e17;
ne_max = ne_cap;
ne_max = min(max(ne_max, 5.0e18), 5.0e20);

ne_ref_for_power = 1.0e20;
Preq_ne_kW = Preq_kW * (ne_max / ne_ref_for_power);
powerMargin = PRF_kW / max(Preq_ne_kW,1e-12);

Ne = (psiN < 1) .* (ne_max - ne_min) .* ...
    (1 - psiN).^1.75 + ne_min;

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
% BUILD COORDINATE GRIDS
% -------------------------------------------------------------------------
[Rgrid, Zgrid] = ndgrid(r, z);

% -------------------------------------------------------------------------
% FLATTEN VARIABLES FOR EXPORT
% -------------------------------------------------------------------------
r_flat   = Rgrid(:);
z_flat   = Zgrid(:);
psi_flat = psi(:);

br_flat  = Br(:);
bt_flat  = Bt(:);
bz_flat  = Bz(:);

ne_flat  = Ne(:);
te_flat  = Te(:);

output = [
    r_flat, ...
    z_flat, ...
    psi_flat, ...
    br_flat, ...
    bt_flat, ...
    bz_flat, ...
    ne_flat, ...
    te_flat
];

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

dlmwrite(outTxt, output, '-append', 'delimiter', ',', 'precision', '%.15E');

disp(['Wrote ', outTxt, ' successfully.']);
disp('NOTE: Bz was forced flat to 0.05 T before psi, Br, ne, Te, and export.');
disp('For flat z-independent Bz, calculated Br = 0 everywhere.');

% -------------------------------------------------------------------------
% COLUMN INDEX GUIDE
% -------------------------------------------------------------------------
disp('Column indices in protoMPEX_profiles.txt:');
disp('1 = r');
disp('2 = z');
disp('3 = psi');
disp('4 = Br calculated consistently from psi');
disp('5 = Bphi (= Bt)');
disp('6 = Bz = 0.05 T');
disp('7 = ne');
disp('8 = te');

% -------------------------------------------------------------------------
% QUICK PREVIEW
% -------------------------------------------------------------------------
disp('First 5 rows of exported data:');
disp(output(1:min(5,size(output,1)), :));