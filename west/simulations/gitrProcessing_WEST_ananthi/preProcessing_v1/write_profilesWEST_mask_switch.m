% close all
% clear all

%% Write SOLEDGE plasma data for GITR input
% Includes in profilesWEST.nc:
%   Coordinates:
%       x, z
%   Electrons:
%       ne, te
%   D+:
%       ni, ti, vr, vt, vz, vp
%   B-field:
%       br, bt, bz
%   Oxygen O1+--O8+:
%       no1--no8, to1--to8,
%       vopar_o1--vopar_o8,
%       vro1--vro8, vto1--vto8, vzo1--vzo8, vpo1--vpo8
%   Geometry:
%       wall_r, wall_z
%       sep1_r, sep1_z
%       sep2_r, sep2_z
%       antenna_centroid_r, antenna_centroid_z
%
% Important:
%   viz_q and vzo_q{q} are treated as PARALLEL PARTICLE FLUXES.
%   They are converted to velocity by dividing by density:
%       v_parallel_D  = viz_q    ./ ni_q
%       v_parallel_Oq = vzo_q{q} ./ no_q{q}

disp('>>>> loading SOLEDGE data')
load('psep1p0mw.mat')

%% ---------------------------
%  User settings
%  ---------------------------
outfile = './profilesWEST.nc';

% Wall masking switch
apply_wall_mask = false;   % true  = zero data outside wall before writing
                           % false = write profiles without wall masking

% Geometry files
wall_file = 'wall_from_mesh.txt';
sep1_file = 'separatrix1_from_mesh.txt';
sep2_file = 'separatrix2_from_mesh.txt';
centroid_file = 'surf_fields_rot.txt';

% Oxygen charge states
nO = 8;

% Experimental OMP ne comparison
compare_exp_omp_ne = false;
exp_ne_file = 'ICRH_57877_8s.mat';
exp_t_min = 5;
exp_t_max = 10;
omp_z_target = 0.0;

% Plot settings
save_plots = false;
plot_dir = 'profilesWEST_plots';

if save_plots && ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% ---------------------------
%  Load geometry from text files
%  ---------------------------
geomIn = loadGeometryFromFiles(wall_file, sep1_file, sep2_file, centroid_file);

%% ---------------------------
%  Coordinates and magnetic field
%  ---------------------------
x0  = double(rgrid(:));
z0  = double(zgrid(:));

nR = length(x0);
nZ = length(z0);

br0 = double(br_q);
bt0 = double(bphi_q);
bz0 = double(bz_q);

x0(~isfinite(x0))   = 0;
z0(~isfinite(z0))   = 0;
br0(~isfinite(br0)) = 0;
bt0(~isfinite(bt0)) = 0;
bz0(~isfinite(bz0)) = 0;

%% ---------------------------
%  Magnetic field unit vectors
%  ---------------------------
b_mag = sqrt(br0.^2 + bt0.^2 + bz0.^2);

ubr = zeros(size(br0));
ubt = zeros(size(bt0));
ubz = zeros(size(bz0));

idxB = b_mag > 0 & isfinite(b_mag);
ubr(idxB) = br0(idxB) ./ b_mag(idxB);
ubt(idxB) = bt0(idxB) ./ b_mag(idxB);
ubz(idxB) = bz0(idxB) ./ b_mag(idxB);

ubr(~isfinite(ubr)) = 0;
ubt(~isfinite(ubt)) = 0;
ubz(~isfinite(ubz)) = 0;

%% ---------------------------
%  Electron and D+ plasma parameters
%  ---------------------------
% Electrons
te0 = 0.5 .* double(te_q);
ne0 = 0.5 .* double(ne_q);

% D+
ti0 = double(ti_q);
ni0 = double(ni_q);

% viz_q is D+ parallel particle flux. Convert flux to velocity.
flux_D = double(viz_q);
vpar_D = flux_D ./ max(ni0, eps);

te0(~isfinite(te0)) = 0;
ne0(~isfinite(ne0)) = 0;
ti0(~isfinite(ti0)) = 0;
ni0(~isfinite(ni0)) = 0;
vpar_D(~isfinite(vpar_D)) = 0;

% Optional old SOLEDGE wall/limiter mask for velocity only
if exist('rq','var') && exist('zq','var') && exist('rW','var') && exist('zW','var')
    try
        maskWallOld = (rq >= rW) & (zq >= zW);
        vpar_D(maskWallOld) = 0;
    catch
        warning('Could not apply rq/zq/rW/zW velocity mask to D+.');
    end
end

% D+ velocity components
vr0 = vpar_D .* ubr;
vt0 = vpar_D .* ubt;
vz0 = vpar_D .* ubz;
vp  = sqrt(vr0.^2 + vt0.^2 + vz0.^2);

vr0(~isfinite(vr0)) = 0;
vt0(~isfinite(vt0)) = 0;
vz0(~isfinite(vz0)) = 0;
vp(~isfinite(vp))   = 0;

%% ---------------------------
%  Oxygen O1+ to O8+
%  ---------------------------
no0     = cell(nO,1);
to0     = cell(nO,1);
vopar0  = cell(nO,1);
vro0    = cell(nO,1);
vto0    = cell(nO,1);
vzo0    = cell(nO,1);
vpo0    = cell(nO,1);

have_oxygen = exist('no_q','var') && exist('to_q','var') && exist('vzo_q','var');

if have_oxygen
    fprintf('>>>> Building oxygen O1+ to O%d+ profiles\n', nO);

    for q = 1:nO
        if q <= numel(no_q) && q <= numel(to_q) && q <= numel(vzo_q) && ...
                ~isempty(no_q{q}) && ~isempty(to_q{q}) && ~isempty(vzo_q{q})

            no0{q} = 0.5 .* double(no_q{q});
            to0{q} = 0.5 .* double(to_q{q});

            % vzo_q{q} is oxygen parallel particle flux. Convert flux to velocity.
            flux_O = double(vzo_q{q});
            vopar0{q} = flux_O ./ max(no0{q}, eps);

            no0{q}(~isfinite(no0{q})) = 0;
            to0{q}(~isfinite(to0{q})) = 0;
            vopar0{q}(~isfinite(vopar0{q})) = 0;

            if exist('maskWallOld','var')
                try
                    vopar0{q}(maskWallOld) = 0;
                catch
                    warning('Could not apply old wall mask to O%d+.', q);
                end
            end

            vro0{q} = vopar0{q} .* ubr;
            vto0{q} = vopar0{q} .* ubt;
            vzo0{q} = vopar0{q} .* ubz;
            vpo0{q} = sqrt(vro0{q}.^2 + vto0{q}.^2 + vzo0{q}.^2);

            vro0{q}(~isfinite(vro0{q})) = 0;
            vto0{q}(~isfinite(vto0{q})) = 0;
            vzo0{q}(~isfinite(vzo0{q})) = 0;
            vpo0{q}(~isfinite(vpo0{q})) = 0;
        else
            warning('Missing oxygen data for O%d+. Writing zeros.', q);
            no0{q}    = zeros(size(ni0));
            to0{q}    = zeros(size(ti0));
            vopar0{q} = zeros(size(ni0));
            vro0{q}   = zeros(size(ni0));
            vto0{q}   = zeros(size(ni0));
            vzo0{q}   = zeros(size(ni0));
            vpo0{q}   = zeros(size(ni0));
        end
    end
else
    warning('Oxygen variables no_q/to_q/vzo_q not found. Writing oxygen fields as zeros.');
    for q = 1:nO
        no0{q}    = zeros(size(ni0));
        to0{q}    = zeros(size(ti0));
        vopar0{q} = zeros(size(ni0));
        vro0{q}   = zeros(size(ni0));
        vto0{q}   = zeros(size(ni0));
        vzo0{q}   = zeros(size(ni0));
        vpo0{q}   = zeros(size(ni0));
    end
end

%% ---------------------------
%  Optional wall polygon mask before writing GITR file
%  ---------------------------
if apply_wall_mask
    assert(geomIn.have_wall, ...
        'apply_wall_mask = true, but wall file is missing or invalid.');

    [RR, ZZ] = meshgrid(x0, z0);
    mask_inside_wall_ZR = inpolygon(RR, ZZ, geomIn.wall_r, geomIn.wall_z);  % [nZ x nR]

    br0 = applyWallMaskKeepOrientation(br0, mask_inside_wall_ZR, 'br0');
    bt0 = applyWallMaskKeepOrientation(bt0, mask_inside_wall_ZR, 'bt0');
    bz0 = applyWallMaskKeepOrientation(bz0, mask_inside_wall_ZR, 'bz0');

    te0 = applyWallMaskKeepOrientation(te0, mask_inside_wall_ZR, 'te0');
    ne0 = applyWallMaskKeepOrientation(ne0, mask_inside_wall_ZR, 'ne0');

    ti0 = applyWallMaskKeepOrientation(ti0, mask_inside_wall_ZR, 'ti0');
    ni0 = applyWallMaskKeepOrientation(ni0, mask_inside_wall_ZR, 'ni0');

    vr0 = applyWallMaskKeepOrientation(vr0, mask_inside_wall_ZR, 'vr0');
    vt0 = applyWallMaskKeepOrientation(vt0, mask_inside_wall_ZR, 'vt0');
    vz0 = applyWallMaskKeepOrientation(vz0, mask_inside_wall_ZR, 'vz0');
    vp  = applyWallMaskKeepOrientation(vp,  mask_inside_wall_ZR, 'vp');

    for q = 1:nO
        no0{q}    = applyWallMaskKeepOrientation(no0{q},    mask_inside_wall_ZR, sprintf('no%d',q));
        to0{q}    = applyWallMaskKeepOrientation(to0{q},    mask_inside_wall_ZR, sprintf('to%d',q));
        vopar0{q} = applyWallMaskKeepOrientation(vopar0{q}, mask_inside_wall_ZR, sprintf('vopar_o%d',q));
        vro0{q}   = applyWallMaskKeepOrientation(vro0{q},   mask_inside_wall_ZR, sprintf('vro%d',q));
        vto0{q}   = applyWallMaskKeepOrientation(vto0{q},   mask_inside_wall_ZR, sprintf('vto%d',q));
        vzo0{q}   = applyWallMaskKeepOrientation(vzo0{q},   mask_inside_wall_ZR, sprintf('vzo%d',q));
        vpo0{q}   = applyWallMaskKeepOrientation(vpo0{q},   mask_inside_wall_ZR, sprintf('vpo%d',q));
    end

    fprintf('Applied wall mask: all data outside wall set to zero before writing %s.\n', outfile);
    fprintf('Wall-mask cells inside wall: %d / %d\n', nnz(mask_inside_wall_ZR), numel(mask_inside_wall_ZR));
else
    fprintf('Wall mask disabled: writing profiles without zeroing outside wall.\n');
end

%% ---------------------------
%  Sanity check before write
%  ---------------------------
vars  = {'x0','z0','br0','bt0','bz0','te0','ne0','ti0','ni0','vr0','vt0','vz0','vp'};
vals  = { x0 , z0 , br0 , bt0 , bz0 , te0 , ne0 , ti0 , ni0 , vr0 , vt0 , vz0 , vp };

for k = 1:length(vars)
    nbad = sum(~isfinite(vals{k}(:)));
    fprintf('%s : bad values = %d\n', vars{k}, nbad);
end

for q = 1:nO
    fprintf('O%d+ bad values: no=%d, to=%d, vopar=%d\n', ...
        q, sum(~isfinite(no0{q}(:))), sum(~isfinite(to0{q}(:))), sum(~isfinite(vopar0{q}(:))));
end

%% ---------------------------
%  Write GITR NetCDF profiles
%  ---------------------------
if isfile(outfile)
    delete(outfile)
end

ncid = netcdf.create(outfile,'NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

% Geometry dimensions
nWall = max(numel(geomIn.wall_r), 1);
nSep1 = max(numel(geomIn.sep1_r), 1);
nSep2 = max(numel(geomIn.sep2_r), 1);
nAnt  = max(numel(geomIn.antenna_centroid_r), 1);

dimWall = netcdf.defDim(ncid,'nWall',nWall);
dimSep1 = netcdf.defDim(ncid,'nSep1',nSep1);
dimSep2 = netcdf.defDim(ncid,'nSep2',nSep2);
dimAnt  = netcdf.defDim(ncid,'nAntennaCentroid',nAnt);

gridRnc = netcdf.defVar(ncid,'x' ,'float',dimR);
gridZnc = netcdf.defVar(ncid,'z' ,'float',dimZ);

% Electron and D+ variables
Ne2Dnc  = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Te2Dnc  = netcdf.defVar(ncid,'te','float',[dimR dimZ]);

Ni2Dnc  = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Ti2Dnc  = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);

vrnc    = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
vtnc    = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
vznc    = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
vpnc    = netcdf.defVar(ncid,'vp','float',[dimR dimZ]);

brnc    = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc    = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc    = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

% Oxygen variables O1+ to O8+
no_nc     = zeros(nO,1);
to_nc     = zeros(nO,1);
vopar_nc  = zeros(nO,1);
vro_nc    = zeros(nO,1);
vto_nc    = zeros(nO,1);
vzo_nc    = zeros(nO,1);
vpo_nc    = zeros(nO,1);

for q = 1:nO
    no_nc(q)    = netcdf.defVar(ncid, sprintf('no%d',q),      'float', [dimR dimZ]);
    to_nc(q)    = netcdf.defVar(ncid, sprintf('to%d',q),      'float', [dimR dimZ]);
    vopar_nc(q) = netcdf.defVar(ncid, sprintf('vopar_o%d',q), 'float', [dimR dimZ]);
    vro_nc(q)   = netcdf.defVar(ncid, sprintf('vro%d',q),     'float', [dimR dimZ]);
    vto_nc(q)   = netcdf.defVar(ncid, sprintf('vto%d',q),     'float', [dimR dimZ]);
    vzo_nc(q)   = netcdf.defVar(ncid, sprintf('vzo%d',q),     'float', [dimR dimZ]);
    vpo_nc(q)   = netcdf.defVar(ncid, sprintf('vpo%d',q),     'float', [dimR dimZ]);
end

% Geometry variables
wallRnc = netcdf.defVar(ncid,'wall_r','float',dimWall);
wallZnc = netcdf.defVar(ncid,'wall_z','float',dimWall);

sep1Rnc = netcdf.defVar(ncid,'sep1_r','float',dimSep1);
sep1Znc = netcdf.defVar(ncid,'sep1_z','float',dimSep1);

sep2Rnc = netcdf.defVar(ncid,'sep2_r','float',dimSep2);
sep2Znc = netcdf.defVar(ncid,'sep2_z','float',dimSep2);

antRnc  = netcdf.defVar(ncid,'antenna_centroid_r','float',dimAnt);
antZnc  = netcdf.defVar(ncid,'antenna_centroid_z','float',dimAnt);

netcdf.endDef(ncid);

% Coordinates
netcdf.putVar(ncid, gridRnc, x0);
netcdf.putVar(ncid, gridZnc, z0);

% Electron and D+
netcdf.putVar(ncid, Ne2Dnc, orientForNetcdfRZ(ne0, nR, nZ, 'ne0'));
netcdf.putVar(ncid, Te2Dnc, orientForNetcdfRZ(te0, nR, nZ, 'te0'));

netcdf.putVar(ncid, Ni2Dnc, orientForNetcdfRZ(ni0, nR, nZ, 'ni0'));
netcdf.putVar(ncid, Ti2Dnc, orientForNetcdfRZ(ti0, nR, nZ, 'ti0'));

netcdf.putVar(ncid, vrnc, orientForNetcdfRZ(vr0, nR, nZ, 'vr0'));
netcdf.putVar(ncid, vtnc, orientForNetcdfRZ(vt0, nR, nZ, 'vt0'));
netcdf.putVar(ncid, vznc, orientForNetcdfRZ(vz0, nR, nZ, 'vz0'));
netcdf.putVar(ncid, vpnc, orientForNetcdfRZ(vp,  nR, nZ, 'vp'));

netcdf.putVar(ncid, brnc, orientForNetcdfRZ(br0, nR, nZ, 'br0'));
netcdf.putVar(ncid, btnc, orientForNetcdfRZ(bt0, nR, nZ, 'bt0'));
netcdf.putVar(ncid, bznc, orientForNetcdfRZ(bz0, nR, nZ, 'bz0'));

% Oxygen O1+ to O8+
for q = 1:nO
    netcdf.putVar(ncid, no_nc(q),    orientForNetcdfRZ(no0{q},     nR, nZ, sprintf('no%d',q)));
    netcdf.putVar(ncid, to_nc(q),    orientForNetcdfRZ(to0{q},     nR, nZ, sprintf('to%d',q)));
    netcdf.putVar(ncid, vopar_nc(q), orientForNetcdfRZ(vopar0{q},  nR, nZ, sprintf('vopar_o%d',q)));
    netcdf.putVar(ncid, vro_nc(q),   orientForNetcdfRZ(vro0{q},    nR, nZ, sprintf('vro%d',q)));
    netcdf.putVar(ncid, vto_nc(q),   orientForNetcdfRZ(vto0{q},    nR, nZ, sprintf('vto%d',q)));
    netcdf.putVar(ncid, vzo_nc(q),   orientForNetcdfRZ(vzo0{q},    nR, nZ, sprintf('vzo%d',q)));
    netcdf.putVar(ncid, vpo_nc(q),   orientForNetcdfRZ(vpo0{q},    nR, nZ, sprintf('vpo%d',q)));
end

% Geometry
netcdf.putVar(ncid, wallRnc, padIfEmpty(geomIn.wall_r));
netcdf.putVar(ncid, wallZnc, padIfEmpty(geomIn.wall_z));

netcdf.putVar(ncid, sep1Rnc, padIfEmpty(geomIn.sep1_r));
netcdf.putVar(ncid, sep1Znc, padIfEmpty(geomIn.sep1_z));

netcdf.putVar(ncid, sep2Rnc, padIfEmpty(geomIn.sep2_r));
netcdf.putVar(ncid, sep2Znc, padIfEmpty(geomIn.sep2_z));

netcdf.putVar(ncid, antRnc, padIfEmpty(geomIn.antenna_centroid_r));
netcdf.putVar(ncid, antZnc, padIfEmpty(geomIn.antenna_centroid_z));

netcdf.close(ncid);

disp('>>>> profilesWEST.nc written successfully')

%% ---------------------------
%  Read back profiles and geometry from GITR NetCDF
%  ---------------------------
disp('>>>> reading profilesWEST.nc for verification and plotting')

R  = double(ncread(outfile,'x'));
z  = double(ncread(outfile,'z'));

bz = double(ncread(outfile,'bz'));
br = double(ncread(outfile,'br'));
bt = double(ncread(outfile,'bt'));

ne = double(ncread(outfile,'ne'));
te = double(ncread(outfile,'te'));
ni = double(ncread(outfile,'ni'));
ti = double(ncread(outfile,'ti'));

vt = double(ncread(outfile,'vt'));
vz = double(ncread(outfile,'vz'));
vr = double(ncread(outfile,'vr'));
vp = double(ncread(outfile,'vp'));

geom = readGeometryFromNc(outfile);

ncvars = {'bz','br','bt','ne','te','ni','ti','vt','vz','vr','vp'};
ncvals = { bz , br , bt , ne , te , ni , ti , vt , vz , vr , vp };

for k = 1:length(ncvars)
    nbad = sum(~isfinite(ncvals{k}(:)));
    fprintf('Readback %s : bad values = %d\n', ncvars{k}, nbad);
end

for q = 1:nO
    tmpNo = double(ncread(outfile, sprintf('no%d',q)));
    tmpTo = double(ncread(outfile, sprintf('to%d',q)));
    fprintf('Readback O%d+: max(no)=%.4e, max(to)=%.4e\n', ...
        q, max(tmpNo(:), [], 'omitnan'), max(tmpTo(:), [], 'omitnan'));
end

%% ---------------------------
%  Plotting: profiles read from .nc file with geometry read from .nc file
%  ---------------------------
plot2DProfile(R, z, bz, 'Input B_z', 'B_z [T]', false, geom, save_plots, plot_dir, 'Bz');
plot2DProfile(R, z, br, 'Input B_R', 'B_R [T]', false, geom, save_plots, plot_dir, 'Br');
plot2DProfile(R, z, bt, 'Input B_t', 'B_t [T]', false, geom, save_plots, plot_dir, 'Bt');

plot2DProfile(R, z, ne, 'Input n_e', 'n_e [m^{-3}]', true, geom, save_plots, plot_dir, 'ne');
plot2DProfile(R, z, ni, 'Input n_i D^+', 'n_i [m^{-3}]', true, geom, save_plots, plot_dir, 'ni_D');
plot2DProfile(R, z, te, 'Input T_e', 'T_e [eV]', false, geom, save_plots, plot_dir, 'Te');
plot2DProfile(R, z, ti, 'Input T_i D^+', 'T_i [eV]', false, geom, save_plots, plot_dir, 'Ti_D');

plot2DProfile(R, z, vt, 'Input V_t D^+', 'V_t [m/s]', false, geom, save_plots, plot_dir, 'Vt_D');
plot2DProfile(R, z, vz, 'Input V_z D^+', 'V_z [m/s]', false, geom, save_plots, plot_dir, 'Vz_D');
plot2DProfile(R, z, vr, 'Input V_R D^+', 'V_R [m/s]', false, geom, save_plots, plot_dir, 'Vr_D');
plot2DProfile(R, z, vp, 'Input |V| D^+', '|V| [m/s]', false, geom, save_plots, plot_dir, 'Vp_D');

% Oxygen density plots
for q = 1:nO
    no_q_nc = double(ncread(outfile, sprintf('no%d',q)));
    to_q_nc = double(ncread(outfile, sprintf('to%d',q)));

    plot2DProfile(R, z, no_q_nc, sprintf('Input n_{O%d+}',q), ...
        sprintf('n_{O%d+} [m^{-3}]',q), true, geom, save_plots, plot_dir, sprintf('no%d',q));

    plot2DProfile(R, z, to_q_nc, sprintf('Input T_{O%d+}',q), ...
        sprintf('T_{O%d+} [eV]',q), false, geom, save_plots, plot_dir, sprintf('to%d',q));
end

writematrix(ne, 'ne.csv')
writematrix(te, 'te.csv')

%% ---------------------------
%  Midplane / mid-poloidal ne and Te lineouts
%  ---------------------------
z_mid_value = 0.5 * (min(z) + max(z));
R_mid_value = 0.5 * (min(R) + max(R));

[~, iz_mid] = min(abs(z - z_mid_value));
[~, iR_mid] = min(abs(R - R_mid_value));

fprintf('Lineout ne/Te vs R at Z = %.6g m, index %d\n', z(iz_mid), iz_mid);
fprintf('Lineout ne/Te vs Z at R = %.6g m, index %d\n', R(iR_mid), iR_mid);

figure('Color','w', 'Name','Mid-Z lineout ne and Te vs R');
yyaxis left
plot(R, ne(:,iz_mid), 'LineWidth', 2.0);
ylabel('$n_e$ [m$^{-3}$]', 'Interpreter','latex');
yyaxis right
plot(R, te(:,iz_mid), 'LineWidth', 2.0);
ylabel('$T_e$ [eV]', 'Interpreter','latex');
xlabel('$R$ [m]', 'Interpreter','latex');
title(sprintf('Mid-Z lineout at Z = %.4f m', z(iz_mid)), 'Interpreter','latex');
set(gca,'FontName','times','FontSize',18);
grid on;

if save_plots
    saveas(gcf, fullfile(plot_dir, 'lineout_midZ_ne_Te_vs_R.png'));
end

figure('Color','w', 'Name','Mid-R lineout ne and Te vs Z');
yyaxis left
plot(z, ne(iR_mid,:), 'LineWidth', 2.0);
ylabel('$n_e$ [m$^{-3}$]', 'Interpreter','latex');
yyaxis right
plot(z, te(iR_mid,:), 'LineWidth', 2.0);
ylabel('$T_e$ [eV]', 'Interpreter','latex');
xlabel('$Z$ [m]', 'Interpreter','latex');
title(sprintf('Mid-R lineout at R = %.4f m', R(iR_mid)), 'Interpreter','latex');
set(gca,'FontName','times','FontSize',18);
grid on;

if save_plots
    saveas(gcf, fullfile(plot_dir, 'lineout_midR_ne_Te_vs_Z.png'));
end

%% ---------------------------
%  Experimental OMP reflectometry comparison: ne(R) at OMP
%  ---------------------------
if compare_exp_omp_ne
    plotExperimentVsSimOMPne(exp_ne_file, exp_t_min, exp_t_max, ...
        R, z, ne, omp_z_target, geom, save_plots, plot_dir);
end

%% ========================================================================
%  Local functions
% ========================================================================

function geom = loadGeometryFromFiles(wall_file, sep1_file, sep2_file, centroid_file)
    geom = struct();

    geom.have_wall = false;
    geom.wall_r = [];
    geom.wall_z = [];

    if isfile(wall_file)
        wall_data = readmatrix(wall_file);
        if size(wall_data,2) >= 2
            wall_r = wall_data(:,1);
            wall_z = wall_data(:,2);
            good = isfinite(wall_r) & isfinite(wall_z);
            wall_r = wall_r(good);
            wall_z = wall_z(good);

            if ~isempty(wall_r)
                if wall_r(1) ~= wall_r(end) || wall_z(1) ~= wall_z(end)
                    wall_r(end+1) = wall_r(1);
                    wall_z(end+1) = wall_z(1);
                end

                geom.have_wall = true;
                geom.wall_r = wall_r;
                geom.wall_z = wall_z;
            end
        end
    end

    geom.have_sep1 = false;
    geom.sep1_r = [];
    geom.sep1_z = [];

    if isfile(sep1_file)
        sep1 = readmatrix(sep1_file);
        if size(sep1,2) >= 2
            r = sep1(:,1);
            z = sep1(:,2);
            good = isfinite(r) & isfinite(z);
            geom.sep1_r = r(good);
            geom.sep1_z = z(good);
            geom.have_sep1 = ~isempty(geom.sep1_r);
        end
    end

    geom.have_sep2 = false;
    geom.sep2_r = [];
    geom.sep2_z = [];

    if isfile(sep2_file)
        sep2 = readmatrix(sep2_file);
        if size(sep2,2) >= 2
            r = sep2(:,1);
            z = sep2(:,2);
            good = isfinite(r) & isfinite(z);
            geom.sep2_r = r(good);
            geom.sep2_z = z(good);
            geom.have_sep2 = ~isempty(geom.sep2_r);
        end
    end

    geom.have_antenna_centroid = false;
    geom.antenna_centroid_r = [];
    geom.antenna_centroid_z = [];

    if isfile(centroid_file)
        centroid = readmatrix(centroid_file);
        if size(centroid,2) >= 3
            zc = centroid(:,3);
            rc = sqrt(centroid(:,1).^2 + centroid(:,2).^2);
            good = isfinite(rc) & isfinite(zc);

            geom.antenna_centroid_r = rc(good);
            geom.antenna_centroid_z = zc(good);
            geom.have_antenna_centroid = ~isempty(geom.antenna_centroid_r);
        end
    end
end

function geom = readGeometryFromNc(ncfile)
    geom = struct();

    geom.wall_r = cleanNcVector(double(ncread(ncfile,'wall_r')));
    geom.wall_z = cleanNcVector(double(ncread(ncfile,'wall_z')));
    geom.have_wall = ~isempty(geom.wall_r) && ~isempty(geom.wall_z);

    geom.sep1_r = cleanNcVector(double(ncread(ncfile,'sep1_r')));
    geom.sep1_z = cleanNcVector(double(ncread(ncfile,'sep1_z')));
    geom.have_sep1 = ~isempty(geom.sep1_r) && ~isempty(geom.sep1_z);

    geom.sep2_r = cleanNcVector(double(ncread(ncfile,'sep2_r')));
    geom.sep2_z = cleanNcVector(double(ncread(ncfile,'sep2_z')));
    geom.have_sep2 = ~isempty(geom.sep2_r) && ~isempty(geom.sep2_z);

    geom.antenna_centroid_r = cleanNcVector(double(ncread(ncfile,'antenna_centroid_r')));
    geom.antenna_centroid_z = cleanNcVector(double(ncread(ncfile,'antenna_centroid_z')));
    geom.have_antenna_centroid = ~isempty(geom.antenna_centroid_r) && ~isempty(geom.antenna_centroid_z);
end

function v = cleanNcVector(v)
    v = v(:);
    v = v(isfinite(v));
    % Remove dummy placeholder when empty geometry was stored as NaN
    v = v(~isnan(v));
end

function v = padIfEmpty(v)
    v = v(:);
    if isempty(v)
        v = NaN;
    end
end

function Aout = applyWallMaskKeepOrientation(Ain, mask_inside_wall_ZR, varname)
    Aout = Ain;

    if isequal(size(Aout), size(mask_inside_wall_ZR))
        Aout(~mask_inside_wall_ZR) = 0;

    elseif isequal(size(Aout), size(mask_inside_wall_ZR.'))
        tmp = Aout.';
        tmp(~mask_inside_wall_ZR) = 0;
        Aout = tmp.';

    else
        error('Size mismatch for %s. Field size %s, wall mask size %s or %s expected.', ...
            varname, mat2str(size(Aout)), mat2str(size(mask_inside_wall_ZR)), mat2str(size(mask_inside_wall_ZR.')));
    end

    Aout(~isfinite(Aout)) = 0;
end

function Awrite = orientForNetcdfRZ(Ain, nR, nZ, varname)
    if isequal(size(Ain), [nR nZ])
        Awrite = Ain;
    elseif isequal(size(Ain), [nZ nR])
        Awrite = Ain.';
    else
        error('Cannot orient %s for NetCDF. Size is %s, expected [%d %d] or [%d %d].', ...
            varname, mat2str(size(Ain)), nR, nZ, nZ, nR);
    end
    Awrite(~isfinite(Awrite)) = 0;
end

function plot2DProfile(R, z, A, titleStr, cbarLabel, useLog, geom, save_plots, plot_dir, fname)
    figure('Color','w', 'Name', titleStr);

    Aplot = A.';
    if useLog
        Aplot(Aplot <= 0) = NaN;
        h = pcolor(R, z, Aplot);
        set(h, 'EdgeColor','none');
        set(gca, 'ColorScale','log');
    else
        h = pcolor(R, z, Aplot);
        set(h, 'EdgeColor','none');
    end

    set(gca,'YDir','normal')
    set(gca,'FontName','times','FontSize',18);
    xlabel('$R$ [m]','Interpreter','latex','FontSize',18);
    ylabel('$Z$ [m]','Interpreter','latex','FontSize',18);
    title(titleStr, 'Interpreter','latex');
    axis equal tight;
    c = colorbar;
    ylabel(c, cbarLabel, 'Interpreter','latex');
    hold on;

    if isfield(geom,'have_sep1') && geom.have_sep1
        plot(geom.sep1_r, geom.sep1_z, 'k--', 'LineWidth', 1.5);
    end

    if isfield(geom,'have_sep2') && geom.have_sep2
        plot(geom.sep2_r, geom.sep2_z, 'k--', 'LineWidth', 1.5);
    end

    if isfield(geom,'have_wall') && geom.have_wall
        plot(geom.wall_r, geom.wall_z, 'k', 'LineWidth', 1.5);
    end

    if isfield(geom,'have_antenna_centroid') && geom.have_antenna_centroid
        plot(geom.antenna_centroid_r, geom.antenna_centroid_z, '.', ...
            'Color', [1 0.65 0.65], 'MarkerSize', 4);
    end

    xlim([min(R) max(R)]);
    ylim([min(z) max(z)]);

    if save_plots
        saveas(gcf, fullfile(plot_dir, [fname '.png']));
    end
end

function plotExperimentVsSimOMPne(exp_ne_file, exp_t_min, exp_t_max, R, z, ne, omp_z_target, geom, save_plots, plot_dir)
    disp('>>>> Plotting experimental vs simulation OMP ne comparison')

    if ~isfile(exp_ne_file)
        warning('Experimental file not found: %s', exp_ne_file);
        return;
    end

    S = load(exp_ne_file);

    try
        t_exp = S.data.WDP.S57877.reflec.t;
        idx_t = find(t_exp >= exp_t_min & t_exp <= exp_t_max);

        if isempty(idx_t)
            warning('No experimental points found in selected time window %.3f -- %.3f s.', ...
                exp_t_min, exp_t_max);
            return;
        end

        r_exp  = S.data.WDP.S57877.reflec.position.r(:,idx_t);
        ne_exp = S.data.WDP.S57877.reflec.ne(:,idx_t);

        r_exp(~isfinite(r_exp)) = NaN;
        ne_exp(~isfinite(ne_exp) | ne_exp <= 0) = NaN;

        r_mean  = mean(r_exp , 2, 'omitnan');
        ne_mean = mean(ne_exp, 2, 'omitnan');
        ne_std  = std(ne_exp , 0, 2, 'omitnan');

        good_exp = isfinite(r_mean) & isfinite(ne_mean) & ne_mean > 0;

        [~, iz_omp] = min(abs(z - omp_z_target));
        ne_sim = ne(:,iz_omp);

        good_sim = isfinite(R) & isfinite(ne_sim) & ne_sim > 0;

        fprintf('Using simulation OMP lineout at Z = %.6f m, index %d\n', z(iz_omp), iz_omp);

        figure('Color','w', 'Name','Experimental vs Simulation OMP ne');
        hold on;

        semilogy(r_exp, ne_exp, '-', ...
            'Color', [0.80 0.80 0.80], ...
            'LineWidth', 0.8);

        errorbar(r_mean(good_exp), ne_mean(good_exp), ne_std(good_exp), ...
            'ko', ...
            'MarkerFaceColor','k', ...
            'LineWidth',1.5, ...
            'DisplayName','Reflectometry average');

        semilogy(R(good_sim), ne_sim(good_sim), ...
            'r-', ...
            'LineWidth',3.0, ...
            'DisplayName', sprintf('SOLEDGE/GITR input, Z = %.4f m', z(iz_omp)));

        set(gca,'YScale','log');
        xlim([2.9 3.1]);
        ylim([1e16 1e20]);

        xlabel('$R$ [m]', 'Interpreter','latex');
        ylabel('$n_e$ [m$^{-3}$]', 'Interpreter','latex');
        title(sprintf('OMP density comparison, %.1f -- %.1f s', exp_t_min, exp_t_max), ...
            'Interpreter','latex');

        set(gca, 'FontName','times', 'FontSize',18);
        grid on;
        box on;
        legend('Location','best');

        if isfield(geom,'have_wall') && geom.have_wall
            xline(min(geom.wall_r), 'k:', 'LineWidth',1.2, ...
                'Label','Wall min R', 'LabelOrientation','horizontal');
            xline(max(geom.wall_r), 'k:', 'LineWidth',1.2, ...
                'Label','Wall max R', 'LabelOrientation','horizontal');
        end

        if save_plots
            saveas(gcf, fullfile(plot_dir, 'comparison_OMP_ne_expt_vs_sim_zoom_log.png'));
        end

    catch ME
        warning('Could not process experimental reflectometry density data: %s', ME.message);
    end
end
