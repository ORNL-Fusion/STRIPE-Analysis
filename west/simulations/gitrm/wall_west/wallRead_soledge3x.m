clear; clc; close all

% ==========================================================
% Read Soledge wall + LCFS / separatrix from mesh.h5
% Save compact geometry MAT file
% ==========================================================

fname  = 'mesh.h5';
outmat = 'west_57877_wall.mat';

%% ----------------------------------------------------------
% 1) Read equilibrium configuration
% -----------------------------------------------------------
psi     = double(h5read(fname, '/config/psi'));
R2D     = double(h5read(fname, '/config/r'));
Z2D     = double(h5read(fname, '/config/z'));
psicore = double(h5read(fname, '/config/psicore'));
psisep1 = double(h5read(fname, '/config/psisep1'));
psisep2 = double(h5read(fname, '/config/psisep2'));
nsep    = double(h5read(fname, '/config/nsep'));

psicore = psicore(1);
psisep1 = psisep1(1);
psisep2 = psisep2(1);
nsep    = nsep(1);

fprintf('=== mesh.h5 equilibrium info ===\n');
fprintf('psicore = %.8e\n', psicore);
fprintf('psisep1 = %.8e\n', psisep1);
fprintf('psisep2 = %.8e\n', psisep2);
fprintf('nsep    = %d\n', nsep);

%% ----------------------------------------------------------
% 2) Build usable coordinate vectors for contourc
% -----------------------------------------------------------
rvec = R2D(1,:);
zvec = Z2D(:,1);

% fallback if mesh orientation is different
if numel(unique(rvec)) < 2 || numel(unique(zvec)) < 2
    rvec = R2D(:,1).';
    zvec = Z2D(1,:).';
end

rvec = double(rvec(:)).';
zvec = double(zvec(:));

% remove duplicates if needed
if any(diff(rvec) == 0)
    [rvec, ir] = unique(rvec, 'stable');
    psi = psi(:, ir);
end

if any(diff(zvec) == 0)
    [zvec, iz] = unique(zvec, 'stable');
    psi = psi(iz, :);
end

% enforce increasing coordinates
if numel(rvec) > 1 && rvec(2) < rvec(1)
    rvec = fliplr(rvec);
    psi  = psi(:, end:-1:1);
end

if numel(zvec) > 1 && zvec(2) < zvec(1)
    zvec = flipud(zvec(:));
    psi  = psi(end:-1:1, :);
end
zvec = zvec(:);

fprintf('size(psi)   = [%d %d]\n', size(psi,1), size(psi,2));
fprintf('numel(rvec) = %d\n', numel(rvec));
fprintf('numel(zvec) = %d\n', numel(zvec));

%% ----------------------------------------------------------
% 3) Read single wall from /wall
% -----------------------------------------------------------
wall_r = double(h5read(fname, '/wall/R'));
wall_z = double(h5read(fname, '/wall/Z'));

wall_r = wall_r(:);
wall_z = wall_z(:);

good = isfinite(wall_r) & isfinite(wall_z);
wall_r = wall_r(good);
wall_z = wall_z(good);

if isempty(wall_r)
    error('Wall data in /wall/R and /wall/Z is empty after cleaning.');
end

% close wall polygon if needed
if wall_r(1) ~= wall_r(end) || wall_z(1) ~= wall_z(end)
    wall_r(end+1) = wall_r(1);
    wall_z(end+1) = wall_z(1);
end

fprintf('\n=== Wall info ===\n');
fprintf('Wall points = %d\n', numel(wall_r));
fprintf('R range     = [%.6f, %.6f]\n', min(wall_r), max(wall_r));
fprintf('Z range     = [%.6f, %.6f]\n', min(wall_z), max(wall_z));

%% ----------------------------------------------------------
% 4) Extract core / separatrix contours
% -----------------------------------------------------------
Ccore = contourc(rvec, zvec, psi, [psicore psicore]);
Csep1 = contourc(rvec, zvec, psi, [psisep1 psisep1]);
Csep2 = contourc(rvec, zvec, psi, [psisep2 psisep2]);

core_branches = parse_contours(Ccore);
sep1_branches = parse_contours(Csep1);
sep2_branches = parse_contours(Csep2);

[Rcore, Zcore] = longest_branch(core_branches);
[Rsep1, Zsep1] = longest_branch(sep1_branches);
[Rsep2, Zsep2] = longest_branch(sep2_branches);

% choose LCFS as longer separatrix branch
len1 = numel(Rsep1);
len2 = numel(Rsep2);

if len1 >= len2
    Rlcfs = Rsep1;
    Zlcfs = Zsep1;
    psi_lcfs = psisep1;
    lcfs_name = 'psisep1';
else
    Rlcfs = Rsep2;
    Zlcfs = Zsep2;
    psi_lcfs = psisep2;
    lcfs_name = 'psisep2';
end

fprintf('\n=== Contour extraction ===\n');
fprintf('Core contour points   : %d\n', numel(Rcore));
fprintf('Separatrix 1 points   : %d\n', numel(Rsep1));
fprintf('Separatrix 2 points   : %d\n', numel(Rsep2));
fprintf('Chosen LCFS           : %s\n', lcfs_name);
fprintf('Chosen LCFS psi value : %.8e\n', psi_lcfs);

%% ----------------------------------------------------------
% 5) Save output txt files
% -----------------------------------------------------------
if ~isempty(Rcore)
    writematrix([Rcore Zcore], 'core_from_mesh.txt');
end

if ~isempty(Rsep1)
    writematrix([Rsep1 Zsep1], 'separatrix1_from_mesh.txt');
end

if ~isempty(Rsep2)
    writematrix([Rsep2 Zsep2], 'separatrix2_from_mesh.txt');
end

if ~isempty(Rlcfs)
    writematrix([Rlcfs Zlcfs], 'LCFS_from_mesh.txt');
end

writematrix([wall_r wall_z], 'wall_from_mesh.txt');

%% ----------------------------------------------------------
% 6) Save MAT file cleanly
% -----------------------------------------------------------
geom = struct();

geom.source_file = fname;

geom.rvec = rvec;
geom.zvec = zvec;
geom.R2D = R2D;
geom.Z2D = Z2D;
geom.psi = psi;

geom.psicore = psicore;
geom.psisep1 = psisep1;
geom.psisep2 = psisep2;
geom.psi_lcfs = psi_lcfs;
geom.nsep = nsep;

geom.wall_r = wall_r;
geom.wall_z = wall_z;

geom.core_r = Rcore;
geom.core_z = Zcore;

geom.sep1_r = Rsep1;
geom.sep1_z = Zsep1;

geom.sep2_r = Rsep2;
geom.sep2_z = Zsep2;

geom.lcfs_r = Rlcfs;
geom.lcfs_z = Zlcfs;

save(outmat, 'geom');
fprintf('\nSaved MAT file: %s\n', outmat);

%% ----------------------------------------------------------
% 7) Plot equilibrium with wall and contours
% -----------------------------------------------------------
figure('Name', 'mesh.h5 equilibrium: wall and separatrix');
contour(rvec, zvec, psi, 40, 'b');
hold on

plot(wall_r, wall_z, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.5)

if ~isempty(Rcore)
    plot(Rcore, Zcore, 'm-', 'LineWidth', 2.0)
end

if ~isempty(Rsep1)
    plot(Rsep1, Zsep1, 'r-', 'LineWidth', 2.0)
end

if ~isempty(Rsep2)
    plot(Rsep2, Zsep2, 'k--', 'LineWidth', 2.0)
end

if ~isempty(Rlcfs)
    plot(Rlcfs, Zlcfs, 'c-', 'LineWidth', 1.5)
end

xlabel('R [m]')
ylabel('Z [m]')
title('mesh.h5 equilibrium: wall and separatrix contours')
axis equal
grid on
set(gca, 'FontSize', 16)

legend_entries = {'Psi contours', 'Wall'};
if ~isempty(Rcore), legend_entries{end+1} = '\psi_{core}'; end
if ~isempty(Rsep1), legend_entries{end+1} = '\psi_{sep1}'; end
if ~isempty(Rsep2), legend_entries{end+1} = '\psi_{sep2}'; end
if ~isempty(Rlcfs), legend_entries{end+1} = 'Chosen LCFS'; end
legend(legend_entries, 'Location', 'best');

%% ----------------------------------------------------------
% 8) Plot LCFS with wall only
% -----------------------------------------------------------
figure('Name', 'LCFS with wall only');
hold on

plot(wall_r, wall_z, 'k-', 'LineWidth', 1.5)

if ~isempty(Rlcfs)
    plot(Rlcfs, Zlcfs, 'r-', 'LineWidth', 2.0)
end

if ~isempty(Rsep1)
    plot(Rsep1, Zsep1, 'b--', 'LineWidth', 1.5)
end

if ~isempty(Rsep2)
    plot(Rsep2, Zsep2, 'g--', 'LineWidth', 1.5)
end

xlabel('R [m]')
ylabel('Z [m]')
title('Wall with LCFS / separatrix')
axis equal
grid on
set(gca, 'FontSize', 16)

legend_entries2 = {'Wall'};
if ~isempty(Rlcfs), legend_entries2{end+1} = 'Chosen LCFS'; end
if ~isempty(Rsep1), legend_entries2{end+1} = '\psi_{sep1}'; end
if ~isempty(Rsep2), legend_entries2{end+1} = '\psi_{sep2}'; end
legend(legend_entries2, 'Location', 'best');

%% ==========================================================
% Local functions
% ==========================================================
function branches = parse_contours(C)
    branches = {};
    if isempty(C)
        return
    end

    k = 1;
    while k < size(C,2)
        npts = C(2,k);
        seg = C(:, k+1:k+npts);

        rr = seg(1,:)';
        zz = seg(2,:)';

        good = isfinite(rr) & isfinite(zz);
        rr = rr(good);
        zz = zz(good);

        if ~isempty(rr)
            if rr(1) ~= rr(end) || zz(1) ~= zz(end)
                rr(end+1) = rr(1);
                zz(end+1) = zz(1);
            end
            branches{end+1} = [rr'; zz']; %#ok<AGROW>
        end

        k = k + npts + 1;
    end
end

function [Rout, Zout] = longest_branch(branches)
    Rout = [];
    Zout = [];

    if isempty(branches)
        return
    end

    imax = 1;
    nmax = 0;
    for i = 1:numel(branches)
        if size(branches{i},2) > nmax
            nmax = size(branches{i},2);
            imax = i;
        end
    end

    seg = branches{imax};
    Rout = seg(1,:)';
    Zout = seg(2,:)';
end