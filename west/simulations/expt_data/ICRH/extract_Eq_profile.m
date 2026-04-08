clear; clc; close all

shot = 57877;
load(['ICRH_' num2str(shot) '_8s.mat'])

S = data.WDP.(['S' num2str(shot)]);
equi = S.equi;

% -----------------------------
% Choose time index
% -----------------------------
it = 150;
t  = equi.time(it);

% -----------------------------
% 2D equilibrium fields
% -----------------------------
R   = equi.RInterp;
Z   = equi.ZInterp;
Psi = squeeze(equi.Psi_interp(it,:,:));

% -----------------------------
% Axis and x-point
% -----------------------------
Rax = equi.mag_ax_R(it);
Zax = equi.mag_ax_Z(it);
Rx  = equi.xPoint(it,1);
Zx  = equi.xPoint(it,2);

% -----------------------------
% LCFS from boundPlasma
% -----------------------------
R_lcfs = squeeze(equi.boundPlasma(it,1,:));
Z_lcfs = squeeze(equi.boundPlasma(it,2,:));

good = isfinite(R_lcfs) & isfinite(Z_lcfs);
R_lcfs = R_lcfs(good);
Z_lcfs = Z_lcfs(good);

if ~isempty(R_lcfs)
    if R_lcfs(1) ~= R_lcfs(end) || Z_lcfs(1) ~= Z_lcfs(end)
        R_lcfs(end+1) = R_lcfs(1);
        Z_lcfs(end+1) = Z_lcfs(1);
    end
end

% -----------------------------
% Limiter point
% -----------------------------
Rlim = equi.limPoint(it,1);
Zlim = equi.limPoint(it,2);

% -----------------------------
% Get separatrix from Psi at X-point
% Use nearest grid point for psi_x
% -----------------------------
dist2 = (R - Rx).^2 + (Z - Zx).^2;
[~, indmin] = min(dist2(:));
psi_x = Psi(indmin);

% Build contour at psi = psi_x using hidden figure
fig_tmp = figure('Visible','off');
[Csep, hsep] = contour(R, Z, Psi, [psi_x psi_x], 'k--');
close(fig_tmp)

% Parse contour matrix and keep longest segment
segments = {};
k = 1;
while k < size(Csep,2)
    level = Csep(1,k); %#ok<NASGU>
    npts  = Csep(2,k);
    seg   = Csep(:, k+1:k+npts);
    segments{end+1} = seg; %#ok<AGROW>
    k = k + npts + 1;
end

R_sep = [];
Z_sep = [];

if ~isempty(segments)
    imax = 1;
    nmax = 0;
    for i = 1:numel(segments)
        if size(segments{i},2) > nmax
            nmax = size(segments{i},2);
            imax = i;
        end
    end

    sep_seg = segments{imax};
    R_sep = sep_seg(1,:)';
    Z_sep = sep_seg(2,:)';

    good = isfinite(R_sep) & isfinite(Z_sep);
    R_sep = R_sep(good);
    Z_sep = Z_sep(good);

    if ~isempty(R_sep)
        if R_sep(1) ~= R_sep(end) || Z_sep(1) ~= Z_sep(end)
            R_sep(end+1) = R_sep(1);
            Z_sep(end+1) = Z_sep(1);
        end
    end
end

% -----------------------------
% Optional wall
% -----------------------------
wall_found = false;
if isfile('wall_coordinates.txt')
    W = readmatrix('wall_coordinates.txt');
    if size(W,2) >= 2
        Rw = W(:,1);
        Zw = W(:,2);
        goodw = isfinite(Rw) & isfinite(Zw);
        Rw = Rw(goodw);
        Zw = Zw(goodw);

        if ~isempty(Rw)
            if Rw(1) ~= Rw(end) || Zw(1) ~= Zw(end)
                Rw(end+1) = Rw(1);
                Zw(end+1) = Zw(1);
            end
            wall_found = true;
        end
    end
end

% -----------------------------
% Save outputs
% -----------------------------
writematrix([R_lcfs Z_lcfs], sprintf('shot_%d_it_%03d_LCFS.txt', shot, it));

if ~isempty(R_sep)
    writematrix([R_sep Z_sep], sprintf('shot_%d_it_%03d_separatrix.txt', shot, it));
else
    % fallback: use LCFS as separatrix
    writematrix([R_lcfs Z_lcfs], sprintf('shot_%d_it_%03d_separatrix.txt', shot, it));
end

% -----------------------------
% Plot
% -----------------------------
figure;
contour(R, Z, Psi, 30, 'b');
hold on

plot(R_lcfs, Z_lcfs, 'r-', 'LineWidth', 2)

if ~isempty(R_sep)
    plot(R_sep, Z_sep, 'k--', 'LineWidth', 1.8)
end

if wall_found
    plot(Rw, Zw, 'g-', 'LineWidth', 2)
end

plot(Rax, Zax, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8)
plot(Rx, Zx, 'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 8)
plot(Rlim, Zlim, 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 8)

xlabel('R [m]')
ylabel('Z [m]')
title(sprintf('Shot %d equilibrium at t = %.4f s', shot, t))
axis equal
grid on
set(gca,'FontSize',14)

if wall_found && ~isempty(R_sep)
    legend('Psi contours','LCFS (boundPlasma)','Separatrix (\psi_X)','Wall', ...
           'Magnetic axis','X-point','Limiter point', 'Location','best')
elseif wall_found
    legend('Psi contours','LCFS (boundPlasma)','Wall', ...
           'Magnetic axis','X-point','Limiter point', 'Location','best')
elseif ~isempty(R_sep)
    legend('Psi contours','LCFS (boundPlasma)','Separatrix (\psi_X)', ...
           'Magnetic axis','X-point','Limiter point', 'Location','best')
else
    legend('Psi contours','LCFS (boundPlasma)', ...
           'Magnetic axis','X-point','Limiter point', 'Location','best')
end

% Plot limits
allR = [R_lcfs(:); Rx; Rax; Rlim];
allZ = [Z_lcfs(:); Zx; Zax; Zlim];

if ~isempty(R_sep)
    allR = [allR; R_sep(:)];
    allZ = [allZ; Z_sep(:)];
end

if wall_found
    allR = [allR; Rw(:)];
    allZ = [allZ; Zw(:)];
end

xlim([min(allR)-0.05, max(allR)+0.05])
ylim([min(allZ)-0.05, max(allZ)+0.05])

fprintf('Shot %d, it = %d, t = %.6f s\n', shot, it, t);
fprintf('X-point : (%.6f, %.6f)\n', Rx, Zx);
fprintf('psi_x   : %.8e\n', psi_x);
fprintf('LCFS points       : %d\n', numel(R_lcfs));
fprintf('Separatrix points : %d\n', numel(R_sep));