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

figure; plot(z, Bz(1,:))

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

% Normalize psi using helicon location (z=1.745, r=0.0625)
[m, n] = size(psi);
psi_norm_val = interp2(z, r, psi, 1.745, 0.0625, 'linear');
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
te_min = 8;   % eV
te_max = 2;   % eV
Te = (psiN < 1) .* (te_max - te_min) .* (1 - psiN).^1.75 + te_min;
Te(psiN > 1) = te_max;   % flatten beyond LCFS

ne_min = 1.0E18;   % m^-3
ne_max = 1.0E20;   % m^-3
Ne = (psiN < 1) .* (ne_max - ne_min) .* (1 - psiN).^1.75 + ne_min;

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

% Write data
ncwrite(outFile,'br',Br);
ncwrite(outFile,'bt',Bt);
ncwrite(outFile,'bz',Bz);
ncwrite(outFile,'ne',Ne);
ncwrite(outFile,'te',Te);

% Add metadata
ncwriteatt(outFile,'/','title','ProtoMPEX profiles with ne, Te, and B-field');
ncwriteatt(outFile,'/','layout','All 2D variables are [r × z]');
ncwriteatt(outFile,'r','units','m');
ncwriteatt(outFile,'z','units','m');
ncwriteatt(outFile,'br','units','tesla');
ncwriteatt(outFile,'bt','units','tesla');
ncwriteatt(outFile,'bz','units','tesla');
ncwriteatt(outFile,'ne','units','m^-3');
ncwriteatt(outFile,'te','units','eV');

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