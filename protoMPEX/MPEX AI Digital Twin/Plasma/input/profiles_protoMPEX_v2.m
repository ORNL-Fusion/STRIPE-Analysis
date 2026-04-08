%% ProtoMPEX: Compute ne, te using original psi formulation and write to one NetCDF
close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING B-FIELD FILE (assumed [r × z])
% -------------------------------------------------------------------------
fileB = 'bfield_protoMPEX.nc';
r = ncread(fileB,'r');                % [m]
z = ncread(fileB,'z');                % [m]
Br = ncread(fileB,'br');              % [r × z]
Bt = ncread(fileB,'bt');              % [r × z]
Bz = ncread(fileB,'bz');              % [r × z]

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
z0 = 1.745;         % stagnation / helicon location [m]
r0 = 0.0625;        % [m]

psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');
psiN = psi ./ psi_norm_val;

figure; imagesc(z, r, psiN);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Normalized $\psi_N$','Interpreter','latex');
colorbar; 
% axis image tight;

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
disp('Constructing Te and ne profiles...');
% te_min = 8;   % eV
% te_max = 2;   % eV
% Te = (psiN < 1) .* (te_max - te_min) .* (1 - psiN).^1.75 + te_min;
% Te(psiN > 1) = te_max;   % flatten beyond LCFS

ne_min = 1.0E18;   % m^-3
ne_max = 1.0E20;   % m^-3
Ne = (psiN < 1) .* (ne_max - ne_min) .* (1 - psiN).^1.75 + ne_min;

% --------- Te: max in core, lower at LCFS, then exponential falloff ----------
Te0       = 2;      % eV at core (psiN=0)  <-- MAX
Te_edge   = 8;      % eV at LCFS (psiN=1)
Te_floor  = 0.2;    % eV far outside
alpha     = 1.75;   % core shaping
lambdaPsi = 0.15;   % SOL falloff width (smaller = steeper)

Te = zeros(size(psiN));

% Core (psiN <= 1): Te(0)=Te0, Te(1)=Te_edge
Te(psiN <= 1) = Te_edge + (Te0 - Te_edge) .* (1 - psiN(psiN <= 1)).^alpha;

% SOL (psiN > 1): exponential starting from Te_edge at psiN=1
Te(psiN > 1)  = Te_floor + (Te_edge - Te_floor) .* exp(-(psiN(psiN > 1) - 1) ./ lambdaPsi);

% -------------------------------------------------------------------------
% BUILD PARALLEL VELOCITY PROFILE FROM Te (Mach ~ 1 with stagnation at z0)
% -------------------------------------------------------------------------
disp('Constructing parallel velocity profile from Te...');

% Physical constants
qe  = 1.602176634e-19;      % [C] electron charge
mp  = 1.67262192369e-27;    % [kg] proton mass (assume H+ ions)

% Ion sound speed [m/s], using Te in eV and isothermal e- fluid:
% c_s = sqrt( e * Te / m_i )
Cs = sqrt(qe .* Te ./ mp);  % same size as Te  [nR × nZ]

% Build a Z-mesh to define Mach as a function of z only
[Zmat, Rmat] = meshgrid(z, r);  %#ok<ASGLU>  % Zmat, Rmat both [nR × nZ]

L_M = 0.10;                     % characteristic length [m] for stagnation region

% Mach number profile (only depends on z): M(z) = tanh((z - z0)/L_M)
M = tanh((Zmat - z0) ./ L_M);   % same size as Cs: [nR × nZ]

% Parallel velocity (sign set by z-z0, |M| -> 1 far from stagnation)
Vpar = M .* Cs;                 % [m/s], [nR × nZ]

% -------------------------------------------------------------------------
% ANALYTICAL PARTICLE AND HEAT FLUXES
% -------------------------------------------------------------------------
disp('Constructing analytical particle and heat fluxes...');

% Parallel particle flux: Gamma_par = n_e * v_par
GammaPar = Ne .* Vpar;          % [m^-2 s^-1]

% Parallel convective heat flux:
% q_par = (5/2) * e * Te(eV) * Gamma_par
Qpar = (5/2) * qe .* Te .* GammaPar;   % [W m^-2]

% -------------------------------------------------------------------------
% PLOT PROFILES
% -------------------------------------------------------------------------
figure; imagesc(z, r, Te);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$T_e$ [eV]','Interpreter','latex');
colorbar; 
% axis image tight;

figure; imagesc(z, r, Ne);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$n_e$ [m$^{-3}$]','Interpreter','latex');
colorbar; 
% axis image tight;

figure; imagesc(z, r, Vpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$v_\parallel$ [m/s]','Interpreter','latex');
colorbar;
% axis image tight;

figure; imagesc(z, r, GammaPar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$\Gamma_\parallel$ [m$^{-2}$ s$^{-1}$]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Qpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$q_\parallel$ [W m$^{-2}$]','Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% WRITE TO ONE NetCDF FILE
% -------------------------------------------------------------------------
disp('Writing combined ProtoMPEX profiles to NetCDF...');

outFile = 'protoMPEX_profiles.nc';
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

% Write data
ncwrite(outFile,'br',Br);
ncwrite(outFile,'bt',Bt);
ncwrite(outFile,'bz',Bz);
ncwrite(outFile,'ne',Ne);
ncwrite(outFile,'te',Te);
ncwrite(outFile,'vpar',Vpar);
ncwrite(outFile,'gamma_par',GammaPar);
ncwrite(outFile,'qpar',Qpar);

% Add metadata
ncwriteatt(outFile,'/','title','ProtoMPEX profiles with ne, Te, B-field, v_parallel, and analytical fluxes');
ncwriteatt(outFile,'/','layout','All 2D variables are [r × z]');
ncwriteatt(outFile,'r','units','m');
ncwriteatt(outFile,'z','units','m');
ncwriteatt(outFile,'br','units','tesla');
ncwriteatt(outFile,'bt','units','tesla');
ncwriteatt(outFile,'bz','units','tesla');
ncwriteatt(outFile,'ne','units','m^-3');
ncwriteatt(outFile,'te','units','eV');
ncwriteatt(outFile,'vpar','units','m s^-1');
ncwriteatt(outFile,'vpar','description','Parallel flow speed, ~Mach 1 away from z0 with stagnation at z0');
ncwriteatt(outFile,'gamma_par','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_par','description','Parallel particle flux: ne * v_parallel');
ncwriteatt(outFile,'qpar','units','W m^-2');
ncwriteatt(outFile,'qpar','description','Convective parallel heat flux: (5/2)*e*Te*Gamma_par');

disp(['✅ Wrote ', outFile, ' successfully.']);

% -------------------------------------------------------------------------
% OPTIONAL: QUICK CHECK
% -------------------------------------------------------------------------
bz_check = ncread(outFile,'bz');
figure; imagesc(z, r, bz_check);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex'); ylabel('$r$ [m]','Interpreter','latex');
title('Read-back $B_z$','Interpreter','latex');
colorbar; 
% axis image tight;
