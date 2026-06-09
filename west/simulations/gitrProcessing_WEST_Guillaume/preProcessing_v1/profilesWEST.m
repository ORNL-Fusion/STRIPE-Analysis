% close all
% clear all

%% Write data for GITR input
disp('>>>> loading SOLEDGE data')
load('psep1p0mw.mat')

%% ---------------------------
%  User settings
%  ---------------------------
outfile = './profilesWEST.nc';

% Wall coordinate file used to zero data outside wall before writing GITR file
wall_file = 'wall_from_mesh.txt';

% Optional plotting overlays
sep1_file = 'separatrix1_from_mesh.txt';
sep2_file = 'separatrix2_from_mesh.txt';
centroid_file = 'surf_fields_rot.txt';

save_plots = false;
plot_dir = 'profilesWEST_plots';

if save_plots && ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% ---------------------------
%  Initialize / assign arrays
%  ---------------------------
x0  = double(rgrid(:));
z0  = double(zgrid(:));

br0 = double(br_q);
bt0 = double(bphi_q);
bz0 = double(bz_q);

% --- D+ plasma parameters ---
ti0 = double(ti_q);          % Ion temperature, D+
te0 = 0.5 .* double(te_q);   % Electron temperature
ne0 = 0.5 .* double(ne_q);   % Electron density
ni0 = double(ni_q);          % Ion density, D+

% --- Parallel/background flow, D+ ---
v0  = double(viz_q);

%% ---------------------------
%  Clean base fields first
%  ---------------------------
x0(~isfinite(x0))   = 0;
z0(~isfinite(z0))   = 0;

br0(~isfinite(br0)) = 0;
bt0(~isfinite(bt0)) = 0;
bz0(~isfinite(bz0)) = 0;

ti0(~isfinite(ti0)) = 0;
te0(~isfinite(te0)) = 0;
ne0(~isfinite(ne0)) = 0;
ni0(~isfinite(ni0)) = 0;
v0(~isfinite(v0))   = 0;

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
%  Optional old SOLEDGE wall/limiter mask for velocity only
%  ---------------------------
if exist('rq','var') && exist('zq','var') && exist('rW','var') && exist('zW','var')
    try
        maskWallOld = (rq >= rW) & (zq >= zW);
        v0(maskWallOld) = 0;
    catch
        warning('Could not apply rq/zq/rW/zW velocity mask.');
    end
end

v0(~isfinite(v0)) = 0;

%% ---------------------------
%  Velocity components
%  ---------------------------
vz0 = v0 .* ubz;
vr0 = v0 .* ubr;
vt0 = v0 .* ubt;

vp = sqrt(vz0.^2 + vr0.^2 + vt0.^2);

vz0(~isfinite(vz0)) = 0;
vr0(~isfinite(vr0)) = 0;
vt0(~isfinite(vt0)) = 0;
vp(~isfinite(vp))   = 0;

%% ---------------------------
%  Final cleanup before wall mask
%  ---------------------------
br0(~isfinite(br0)) = 0;
bt0(~isfinite(bt0)) = 0;
bz0(~isfinite(bz0)) = 0;

ti0(~isfinite(ti0)) = 0;
te0(~isfinite(te0)) = 0;
ne0(~isfinite(ne0)) = 0;
ni0(~isfinite(ni0)) = 0;

%% ---------------------------
%  Apply wall polygon mask to all 2D profiles before writing GITR file
%  ---------------------------
nR = length(x0);
nZ = length(z0);

assert(isfile(wall_file), 'Wall file not found: %s', wall_file);

wall_data = readmatrix(wall_file);
assert(size(wall_data,2) >= 2, 'Wall file must contain at least two columns: [R Z].');

wall_r = wall_data(:,1);
wall_z = wall_data(:,2);

good_wall = isfinite(wall_r) & isfinite(wall_z);
wall_r = wall_r(good_wall);
wall_z = wall_z(good_wall);

assert(~isempty(wall_r), 'Wall polygon is empty after removing non-finite points.');

if wall_r(1) ~= wall_r(end) || wall_z(1) ~= wall_z(end)
    wall_r(end+1) = wall_r(1);
    wall_z(end+1) = wall_z(1);
end

[RR, ZZ] = meshgrid(x0, z0);
mask_inside_wall_ZR = inpolygon(RR, ZZ, wall_r, wall_z);  % [nZ x nR]

br0 = applyWallMaskKeepOrientation(br0, mask_inside_wall_ZR, 'br0');
bt0 = applyWallMaskKeepOrientation(bt0, mask_inside_wall_ZR, 'bt0');
bz0 = applyWallMaskKeepOrientation(bz0, mask_inside_wall_ZR, 'bz0');

ti0 = applyWallMaskKeepOrientation(ti0, mask_inside_wall_ZR, 'ti0');
te0 = applyWallMaskKeepOrientation(te0, mask_inside_wall_ZR, 'te0');
ne0 = applyWallMaskKeepOrientation(ne0, mask_inside_wall_ZR, 'ne0');
ni0 = applyWallMaskKeepOrientation(ni0, mask_inside_wall_ZR, 'ni0');

v0  = applyWallMaskKeepOrientation(v0,  mask_inside_wall_ZR, 'v0');
vr0 = applyWallMaskKeepOrientation(vr0, mask_inside_wall_ZR, 'vr0');
vt0 = applyWallMaskKeepOrientation(vt0, mask_inside_wall_ZR, 'vt0');
vz0 = applyWallMaskKeepOrientation(vz0, mask_inside_wall_ZR, 'vz0');
vp  = applyWallMaskKeepOrientation(vp,  mask_inside_wall_ZR, 'vp');

fprintf('Applied wall mask: all data outside wall set to zero before writing %s.\n', outfile);
fprintf('Wall-mask cells inside wall: %d / %d\n', nnz(mask_inside_wall_ZR), numel(mask_inside_wall_ZR));

%% ---------------------------
%  Sanity check before write
%  ---------------------------
vars  = {'x0','z0','br0','bt0','bz0','ti0','te0','ne0','ni0','vr0','vt0','vz0','vp'};
vals  = { x0 , z0 , br0 , bt0 , bz0 , ti0 , te0 , ne0 , ni0 , vr0 , vt0 , vz0 , vp };

for k = 1:length(vars)
    nbad = sum(~isfinite(vals{k}(:)));
    fprintf('%s : bad values = %d\n', vars{k}, nbad);
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

gridRnc = netcdf.defVar(ncid,'x' ,'float',dimR);
gridZnc = netcdf.defVar(ncid,'z' ,'float',dimZ);

Ne2Dnc  = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc  = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc  = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc  = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);

vrnc    = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
vtnc    = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
vznc    = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
vpnc    = netcdf.defVar(ncid,'vp','float',[dimR dimZ]);

brnc    = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc    = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc    = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);

netcdf.putVar(ncid, gridRnc, x0);
netcdf.putVar(ncid, gridZnc, z0);

netcdf.putVar(ncid, Ne2Dnc, orientForNetcdfRZ(ne0, nR, nZ, 'ne0'));
netcdf.putVar(ncid, Ni2Dnc, orientForNetcdfRZ(ni0, nR, nZ, 'ni0'));

netcdf.putVar(ncid, Te2Dnc, orientForNetcdfRZ(te0, nR, nZ, 'te0'));
netcdf.putVar(ncid, Ti2Dnc, orientForNetcdfRZ(ti0, nR, nZ, 'ti0'));

netcdf.putVar(ncid, vtnc, orientForNetcdfRZ(vt0, nR, nZ, 'vt0'));
netcdf.putVar(ncid, vrnc, orientForNetcdfRZ(vr0, nR, nZ, 'vr0'));
netcdf.putVar(ncid, vznc, orientForNetcdfRZ(vz0, nR, nZ, 'vz0'));
netcdf.putVar(ncid, vpnc, orientForNetcdfRZ(vp,  nR, nZ, 'vp'));

netcdf.putVar(ncid, brnc, orientForNetcdfRZ(br0, nR, nZ, 'br0'));
netcdf.putVar(ncid, btnc, orientForNetcdfRZ(bt0, nR, nZ, 'bt0'));
netcdf.putVar(ncid, bznc, orientForNetcdfRZ(bz0, nR, nZ, 'bz0'));

netcdf.close(ncid);

disp('>>>> profilesWEST.nc written successfully')

%% ---------------------------
%  Read back and verify
%  ---------------------------
R  = double(ncread(outfile,'x'));
z  = double(ncread(outfile,'z'));

bz = double(ncread(outfile,'bz'));
br = double(ncread(outfile,'br'));
bt = double(ncread(outfile,'bt'));

ne = double(ncread(outfile,'ne'));
ni = double(ncread(outfile,'ni'));
te = double(ncread(outfile,'te'));
ti = double(ncread(outfile,'ti'));

vt = double(ncread(outfile,'vt'));
vz = double(ncread(outfile,'vz'));
vr = double(ncread(outfile,'vr'));
vp = double(ncread(outfile,'vp'));

ncvars = {'bz','br','bt','ne','ni','te','ti','vt','vz','vr','vp'};
ncvals = { bz , br , bt , ne , ni , te , ti , vt , vz , vr , vp };

for k = 1:length(ncvars)
    nbad = sum(~isfinite(ncvals{k}(:)));
    fprintf('Readback %s : bad values = %d\n', ncvars{k}, nbad);
end

[RR2, ZZ2] = meshgrid(R, z);
mask_inside_wall_ZR_read = inpolygon(RR2, ZZ2, wall_r, wall_z);
mask_outside_RZ = ~mask_inside_wall_ZR_read.';  % [nR x nZ]

fprintf('Readback outside-wall max values:\n');
fprintf('  ne = %.4e\n', max(abs(ne(mask_outside_RZ)), [], 'omitnan'));
fprintf('  te = %.4e\n', max(abs(te(mask_outside_RZ)), [], 'omitnan'));
fprintf('  ni = %.4e\n', max(abs(ni(mask_outside_RZ)), [], 'omitnan'));
fprintf('  ti = %.4e\n', max(abs(ti(mask_outside_RZ)), [], 'omitnan'));

%% ---------------------------
%  Load optional plotting overlays
%  ---------------------------
geom = struct();
geom.have_wall = true;
geom.wall_r = wall_r;
geom.wall_z = wall_z;

geom.have_sep = false;
geom.sep1_r = [];
geom.sep1_z = [];
geom.sep2_r = [];
geom.sep2_z = [];

if isfile(sep1_file)
    sep1 = readmatrix(sep1_file);
    if size(sep1,2) >= 2
        good = isfinite(sep1(:,1)) & isfinite(sep1(:,2));
        geom.sep1_r = sep1(good,1);
        geom.sep1_z = sep1(good,2);
        geom.have_sep = true;
    end
end

if isfile(sep2_file)
    sep2 = readmatrix(sep2_file);
    if size(sep2,2) >= 2
        good = isfinite(sep2(:,1)) & isfinite(sep2(:,2));
        geom.sep2_r = sep2(good,1);
        geom.sep2_z = sep2(good,2);
        geom.have_sep = true;
    end
end

geom.have_antenna = false;
geom.antenna_r = [];
geom.antenna_z = [];
geom.centroid_r = [];
geom.centroid_z = [];

if isfile(centroid_file)
    centroid = readmatrix(centroid_file);
    if size(centroid,2) >= 3
        zc = centroid(:,3);
        rc = sqrt(centroid(:,1).^2 + centroid(:,2).^2);

        good = isfinite(rc) & isfinite(zc);
        rc = rc(good);
        zc = zc(good);

        geom.centroid_r = rc;
        geom.centroid_z = zc;

        if numel(rc) >= 3
            pts = unique([rc zc], 'rows');

            if size(pts,1) >= 3
                try
                    k = boundary(pts(:,1), pts(:,2), 0.85);
                catch
                    k = convhull(pts(:,1), pts(:,2));
                end

                geom.antenna_r = pts(k,1);
                geom.antenna_z = pts(k,2);

                if geom.antenna_r(1) ~= geom.antenna_r(end) || geom.antenna_z(1) ~= geom.antenna_z(end)
                    geom.antenna_r(end+1) = geom.antenna_r(1);
                    geom.antenna_z(end+1) = geom.antenna_z(1);
                end

                geom.have_antenna = true;
            end
        end
    end
end

%% ---------------------------
%  Plotting: 2D maps with wall coordinates
%  ---------------------------
plot2DProfile(R, z, bz, 'Input B_z', 'B_z [T]', false, geom, save_plots, plot_dir, 'Bz');
plot2DProfile(R, z, br, 'Input B_R', 'B_R [T]', false, geom, save_plots, plot_dir, 'Br');
plot2DProfile(R, z, bt, 'Input B_t', 'B_t [T]', false, geom, save_plots, plot_dir, 'Bt');

plot2DProfile(R, z, ne, 'Input n_e', 'n_e [m^{-3}]', true, geom, save_plots, plot_dir, 'ne');
plot2DProfile(R, z, ni, 'Input n_i', 'n_i [m^{-3}]', true, geom, save_plots, plot_dir, 'ni');
plot2DProfile(R, z, te, 'Input T_e', 'T_e [eV]', false, geom, save_plots, plot_dir, 'Te');
plot2DProfile(R, z, ti, 'Input T_i', 'T_i [eV]', false, geom, save_plots, plot_dir, 'Ti');

plot2DProfile(R, z, vt, 'Input V_t', 'V_t [m/s]', false, geom, save_plots, plot_dir, 'Vt');
plot2DProfile(R, z, vz, 'Input V_z', 'V_z [m/s]', false, geom, save_plots, plot_dir, 'Vz');
plot2DProfile(R, z, vr, 'Input V_R', 'V_R [m/s]', false, geom, save_plots, plot_dir, 'Vr');
plot2DProfile(R, z, vp, 'Input |V|', '|V| [m/s]', false, geom, save_plots, plot_dir, 'Vp');

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

%% ========================================================================
%  Local functions
% ========================================================================

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

    if isfield(geom,'have_sep') && geom.have_sep
        if ~isempty(geom.sep1_r)
            plot(geom.sep1_r, geom.sep1_z, 'k--', 'LineWidth', 1.5);
        end
        if ~isempty(geom.sep2_r)
            plot(geom.sep2_r, geom.sep2_z, 'k--', 'LineWidth', 1.5);
        end
    end

    if isfield(geom,'have_wall') && geom.have_wall
        plot(geom.wall_r, geom.wall_z, 'k', 'LineWidth', 1.5);
    end

    if isfield(geom,'centroid_r') && ~isempty(geom.centroid_r)
        plot(geom.centroid_r, geom.centroid_z, '.', 'Color', [1 0.65 0.65], 'MarkerSize', 4);
    end

    if isfield(geom,'have_antenna') && geom.have_antenna
        patch(geom.antenna_r, geom.antenna_z, 'r', ...
            'EdgeColor','k', 'LineWidth', 1.5, 'FaceAlpha', 1.0);
    end

    xlim([min(R) max(R)]);
    ylim([min(z) max(z)]);

    if save_plots
        saveas(gcf, fullfile(plot_dir, [fname '.png']));
    end
end
