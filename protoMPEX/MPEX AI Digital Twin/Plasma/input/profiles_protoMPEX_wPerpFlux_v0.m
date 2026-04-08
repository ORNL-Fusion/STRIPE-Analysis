%% ProtoMPEX: Compute ne, te using original psi formulation and write to one NetCDF
% Perp transport CHECK with CONSTANT coefficients:
%   D_perp   = 0.5  m^2/s
%   chi_perp = 1.0  m^2/s
%
% Perp flux definition (radial, positive):
%   Gamma_perp = D_perp * |dn/dr|
%   q_perp     = chi_perp * ne * e * |dTe/dr|
%
% Diagnostics:
%   plots of dn/dr and dTe/dr to locate stripe source

close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING B-FIELD FILE (assumed [r × z])
% -------------------------------------------------------------------------
fileB = 'bfield_protoMPEX.nc';
r  = ncread(fileB,'r');                % [m]
z  = ncread(fileB,'z');                % [m]
Br = ncread(fileB,'br');               % [r × z]
Bt = ncread(fileB,'bt');               % [r × z]
Bz = ncread(fileB,'bz');               % [r × z]

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
figure; imagesc(z, r, Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Input $B_z$ [T]','Interpreter','latex');
colorbar; axis image tight; xlim([0.5 4.2]);

% -------------------------------------------------------------------------
% CALCULATE PSI USING ORIGINAL METHOD (YOUR LOOP)
% -------------------------------------------------------------------------
disp('Calculating psi using original summation method...');
dx = r(2) - r(1);
psi = zeros(nR, nZ);
for ii = 1:nR
    for jj = 1:nZ
        psi(ii,jj) = 2*pi * sum(Bz(1:ii,jj) .* r(1:ii)) * dx;
    end
end

% Normalize psi using helicon / stagnation location (z=1.745, r=0.0625)
z0 = 1.745;         % [m]
r0 = 0.0625;        % [m]
psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');
psiN = psi ./ psi_norm_val;

figure; imagesc(z, r, psiN);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Normalized $\psi_N$','Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
disp('Constructing Te and ne profiles...');
te_min = 8;   % eV
te_max = 2;   % eV
Te = (psiN < 1) .* (te_max - te_min) .* (1 - psiN).^1.75 + te_min;
Te(psiN > 1) = te_max;   % flatten beyond LCFS

ne_min = 1.0E18;   % m^-3
ne_max = 1.0E20;   % m^-3
Ne = (psiN < 1) .* (ne_max - ne_min) .* (1 - psiN).^1.75 + ne_min;

% -------------------------------------------------------------------------
% BUILD PARALLEL VELOCITY PROFILE FROM Te (Mach ~ 1 with stagnation at z0)
% -------------------------------------------------------------------------
disp('Constructing parallel velocity profile from Te...');

qe  = 1.602176634e-19;      % [C]
mp  = 1.67262192369e-27;    % [kg] proton mass (assume H+)

Cs = sqrt(qe .* Te ./ mp);  % [m/s]
[Zmat, ~] = meshgrid(z, r); % [nR × nZ]
L_M = 0.10;
M = tanh((Zmat - z0) ./ L_M);
Vpar = M .* Cs;

% -------------------------------------------------------------------------
% ANALYTICAL PARTICLE AND HEAT FLUXES (PARALLEL)
% -------------------------------------------------------------------------
disp('Constructing analytical particle and heat fluxes...');
GammaPar = Ne .* Vpar;                        % [m^-2 s^-1]
Qpar     = (5/2) * qe .* Te .* GammaPar;      % [W m^-2]

% -------------------------------------------------------------------------
% PERPENDICULAR FLUXES (CONSTANT COEFFICIENT CHECK)
% -------------------------------------------------------------------------

% --- Smoothing ONLY for gradient evaluation (reduces z-striping artifacts)
% window sizes in index space (tune if needed)
wz = 9;     % z-window (odd)   try 7–21
wr = 1;     % keep radial untouched (set 1) unless you want slight r smoothing

Ne_g = Ne;
Te_g = Te;

% Smooth along z (dimension 2) using Savitzky-Golay (good for derivatives)
Ne_g = smoothdata(Ne_g, 2, 'sgolay', wz);
Te_g = smoothdata(Te_g, 2, 'sgolay', wz);

% (optional) if still noisy, add a light moving mean along z
Ne_g = smoothdata(Ne_g, 2, 'movmean', wz);
Te_g = smoothdata(Te_g, 2, 'movmean', wz);
disp('Computing perpendicular particle + heat fluxes (constant coefficients)...');

Dperp   = 0.5;   % [m^2/s] constant check
chiPerp = 1.0;   % [m^2/s] constant check

dr = r(2) - r(1);
dz = z(2) - z(1);

% IMPORTANT: use gradient with correct spacings
[dndr, ~] = gradient(Ne_g, dr, dz);
[dTdr, ~] = gradient(Te_g, dr, dz);

% Radial-only perpendicular fluxes, forced positive
GammaPerp_r = Dperp .* abs(dndr);                 % [m^-2 s^-1]
GammaPerp_z = zeros(size(GammaPerp_r));
GammaPerp   = GammaPerp_r;

Qperp_r = chiPerp .* Ne .* (qe .* abs(dTdr));     % [W/m^2]
Qperp_z = zeros(size(Qperp_r));
Qperp   = Qperp_r;



% -------------------------------------------------------------------------
% DIAGNOSTICS: IF STRIPES EXIST HERE, THEY COME FROM dn/dr or dTe/dr
% -------------------------------------------------------------------------
figure; imagesc(z, r, abs(dndr));
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('$|dn_e/dr|$ (diagnostic)','Interpreter','latex'); colorbar;

figure; imagesc(z, r, abs(dTdr));
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('$|dT_e/dr|$ (diagnostic)','Interpreter','latex'); colorbar;

% -------------------------------------------------------------------------
% PLOT PROFILES
% -------------------------------------------------------------------------
figure; imagesc(z, r, Te);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('$T_e$ [eV]','Interpreter','latex'); colorbar;

figure; imagesc(z, r, Ne);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('$n_e$ [m$^{-3}$]','Interpreter','latex'); colorbar;

figure; imagesc(z, r, GammaPar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('$\Gamma_\parallel$ [m$^{-2}$ s$^{-1}$]','Interpreter','latex'); colorbar;

figure; imagesc(z, r, Qpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('$q_\parallel$ [W m$^{-2}$]','Interpreter','latex'); colorbar;

figure; imagesc(z, r, GammaPerp);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title(['$\Gamma_\perp$ (radial, const $D_\perp$) [m$^{-2}$ s$^{-1}$], $D_\perp$=', num2str(Dperp)],'Interpreter','latex');
colorbar;

figure; imagesc(z, r, Qperp);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title(['$q_\perp$ (radial, const $\chi_\perp$) [W m$^{-2}$], $\chi_\perp$=', num2str(chiPerp)],'Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% WRITE TO ONE NetCDF FILE
% -------------------------------------------------------------------------
disp('Writing combined ProtoMPEX profiles to NetCDF...');

outFile = 'protoMPEX_profiles_constPerp.nc';
if exist(outFile,'file'); delete(outFile); end
cmpr = {'Format','netcdf4','DeflateLevel',4,'Shuffle',true};

% Dimensions
nccreate(outFile,'r','Dimensions',{'r',nR},'Datatype','double',cmpr{:});
nccreate(outFile,'z','Dimensions',{'z',nZ},'Datatype','double',cmpr{:});
ncwrite(outFile,'r',r);
ncwrite(outFile,'z',z);

% Variables (all [r × z])
nccreate(outFile,'br','Dimensions',{'r',nR,'z',nZ},'Datatype','single',cmpr{:});
nccreate(outFile,'bt','Dimensions',{'r',nR,'z',nZ},'Datatype','single',cmpr{:});
nccreate(outFile,'bz','Dimensions',{'r',nR,'z',nZ},'Datatype','single',cmpr{:});
nccreate(outFile,'ne','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'te','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'vpar','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_par','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qpar','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Constant coefficients stored as scalars
nccreate(outFile,'d_perp','Datatype','double',cmpr{:});
nccreate(outFile,'chi_perp','Datatype','double',cmpr{:});

% Perpendicular particle fluxes (radial, positive)
nccreate(outFile,'gamma_perp','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_perp_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_perp_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Perpendicular heat fluxes (radial, positive)
nccreate(outFile,'qperp','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qperp_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qperp_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Write data
ncwrite(outFile,'br',Br);
ncwrite(outFile,'bt',Bt);
ncwrite(outFile,'bz',Bz);
ncwrite(outFile,'ne',Ne);
ncwrite(outFile,'te',Te);
ncwrite(outFile,'vpar',Vpar);
ncwrite(outFile,'gamma_par',GammaPar);
ncwrite(outFile,'qpar',Qpar);

ncwrite(outFile,'d_perp',Dperp);
ncwrite(outFile,'chi_perp',chiPerp);

ncwrite(outFile,'gamma_perp',GammaPerp);
ncwrite(outFile,'gamma_perp_r',GammaPerp_r);
ncwrite(outFile,'gamma_perp_z',GammaPerp_z);

ncwrite(outFile,'qperp',Qperp);
ncwrite(outFile,'qperp_r',Qperp_r);
ncwrite(outFile,'qperp_z',Qperp_z);

% Metadata
ncwriteatt(outFile,'/','title','ProtoMPEX profiles with constant perpendicular transport coefficients (radial)');
ncwriteatt(outFile,'/','layout','All 2D variables are [r × z]');
ncwriteatt(outFile,'/','perp_definition','Gamma_perp = D_perp*|dn/dr|, q_perp = chi_perp*ne*e*|dTe/dr| (radial only)');
ncwriteatt(outFile,'/','perp_note','Perp fluxes positive by abs(). d_perp and chi_perp stored as scalars.');

disp(['✅ Wrote ', outFile, ' successfully.']);