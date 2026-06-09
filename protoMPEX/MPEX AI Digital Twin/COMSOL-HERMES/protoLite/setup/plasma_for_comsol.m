%% ProtoLite: Compute psi, ne, Te and write COMSOL-like TXT

close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING PROTOLITE B-FIELD FILE
% -------------------------------------------------------------------------
fileB = 'bfield_protoLite.nc';

r  = ncread(fileB,'r');
z  = ncread(fileB,'z');
Br = ncread(fileB,'br');
Bt = ncread(fileB,'bt');
Bz = ncread(fileB,'bz');

r = r(:);
z = z(:);

if ~isequal(size(Bz), [numel(r), numel(z)])
    Br = permute(Br,[2 1]);
    Bt = permute(Bt,[2 1]);
    Bz = permute(Bz,[2 1]);
end

nR = numel(r);
nZ = numel(z);

zL = z(1);
zR = z(end);

% -------------------------------------------------------------------------
% READ OR COMPUTE PSI
% -------------------------------------------------------------------------
infoB = ncinfo(fileB);
varNames = {infoB.Variables.Name};

if any(strcmp(varNames,'psi'))
    disp('Reading psi directly from B-field file...');
    psi = ncread(fileB,'psi');

    if ~isequal(size(psi), [nR nZ])
        psi = permute(psi,[2 1]);
    end
else
    disp('psi not found in B-field file. Computing psi from Bz...');

    dr = r(2) - r(1);
    psi = zeros(nR,nZ);

    for ii = 1:nR
        for jj = 1:nZ
            psi(ii,jj) = 2*pi * sum(Bz(1:ii,jj) .* r(1:ii)) * dr;
        end
    end
end

% -------------------------------------------------------------------------
% NORMALIZE PSI USING PROTOLITE HELICON / STAGNATION LOCATION
% -------------------------------------------------------------------------
z0 = min(max(1.0, zL), zR);       % helicon center / stagnation point
r0 = min(max(0.06, r(1)), r(end));

psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');

if ~isfinite(psi_norm_val) || abs(psi_norm_val) < 1e-30
    error('Bad psi normalization value. Check z0/r0 and B-field file.');
end

psiN = psi ./ psi_norm_val;

% -------------------------------------------------------------------------
% PLOT INPUT FIELDS
% -------------------------------------------------------------------------
figure;
imagesc(z,r,Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $B_z$ [T]','Interpreter','latex');
colorbar;
xlim([zL zR]);

figure;
imagesc(z,r,Br);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $B_r$ [T]','Interpreter','latex');
colorbar;
xlim([zL zR]);

figure;
imagesc(z,r,psiN);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite normalized $\psi_N$','Interpreter','latex');
colorbar;
xlim([zL zR]);

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES
% -------------------------------------------------------------------------
disp('Constructing ProtoLite Te and ne profiles...');

qe = 1.602176634e-19;

Te0       = 2.0;
Te_edge   = 8.0;
Te_floor  = 0.2;
alpha     = 1.75;
lambdaPsi = 0.15;

Te = zeros(size(psiN));

inside  = psiN <= 1;
outside = psiN > 1;

Te(inside) = Te_edge + ...
    (Te0 - Te_edge) .* (1 - psiN(inside)).^alpha;

Te(outside) = Te_floor + ...
    (Te_edge - Te_floor) .* exp(-(psiN(outside) - 1) ./ lambdaPsi);

gasPuff_Gps = 1.0e20;
PRF_kW      = 180.0;
Eion_eV     = 1.0e3;

ne_min = 1.0e17;
ne_cap = 4.5e19;
ne_max = min(max(ne_cap, 5.0e18), 4.5e19);

Preq_kW = gasPuff_Gps * Eion_eV * qe / 1e3;
powerFactor = min(PRF_kW / max(Preq_kW,1e-12), 1.0);

ne_ref_for_power = 1.0e20;
Preq_ne_kW = Preq_kW * (ne_max / ne_ref_for_power);
powerMargin = PRF_kW / max(Preq_ne_kW,1e-12);

Ne = ne_min .* ones(size(psiN));
Ne(inside) = (ne_max - ne_min) .* (1 - psiN(inside)).^1.75 + ne_min;

disp(sprintf(['ProtoLite ne_max = %.3e m^-3\n', ...
    'Required power for ne_max: %.2f kW, available PRF = %.1f kW\n', ...
    'Power margin = %.3f'], ...
    ne_max, Preq_ne_kW, PRF_kW, powerMargin));

% -------------------------------------------------------------------------
% CHECK PLOTS
% -------------------------------------------------------------------------
figure;
imagesc(z,r,Ne);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $n_e$ [m$^{-3}$]','Interpreter','latex');
colorbar;
xlim([zL zR]);

figure;
imagesc(z,r,Te);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $T_e$ [eV]','Interpreter','latex');
colorbar;
xlim([zL zR]);

% -------------------------------------------------------------------------
% BUILD COMSOL-LIKE OUTPUT MATRIX
% -------------------------------------------------------------------------
[Rgrid, Zgrid] = ndgrid(r,z);

output = [ ...
    Rgrid(:), ...
    Zgrid(:), ...
    psi(:), ...
    Br(:), ...
    Bt(:), ...
    Bz(:), ...
    Ne(:), ...
    Te(:)];

% -------------------------------------------------------------------------
% WRITE COMSOL-LIKE TXT FILE
% -------------------------------------------------------------------------
disp('Writing ProtoLite COMSOL-like TXT file...');

outTxt = 'protoLite_profiles.txt';

if exist(outTxt,'file')
    delete(outTxt);
end

fid = fopen(outTxt,'w');

fprintf(fid,'%% Model:              protoLite_generated_profiles\n');
fprintf(fid,'%% Version:            COMSOL-like export from MATLAB\n');
fprintf(fid,'%% Date:               %s\n', datestr(now));
fprintf(fid,'%% Dimension:          2\n');
fprintf(fid,'%% Nodes:              %d\n', nR*nZ);
fprintf(fid,'%% Expressions:        6\n');
fprintf(fid,'%% Description:        Poloidal magnetic flux, Magnetic flux density r-component, Magnetic flux density phi-component, Magnetic flux density z-component, Electron density, Electron temperature\n');
fprintf(fid,'%% Length unit:        m\n');
fprintf(fid,'%% r,z,mf.Psi (Wb),mf.Br (T),mf.Bphi (T),mf.Bz (T),ne (m^-3),te (eV)\n');

fclose(fid);

dlmwrite(outTxt, output, '-append', 'delimiter', ',', 'precision', '%.15E');

disp(['Wrote ', outTxt, ' successfully.']);

disp('Column indices:');
disp('1 = r');
disp('2 = z');
disp('3 = psi');
disp('4 = Br');
disp('5 = Bphi = Bt');
disp('6 = Bz');
disp('7 = ne');
disp('8 = te');

disp('First 5 rows:');
disp(output(1:min(5,size(output,1)),:));

disp('End of ProtoLite COMSOL-like TXT export script.');