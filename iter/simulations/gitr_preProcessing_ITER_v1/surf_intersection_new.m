clc; close all; clear all;

% Load geometry
fid = fopen('gitrGeometryPointPlane3d.cfg');
centroid = readmatrix('centroid.csv');

tline = fgetl(fid); % Skip lines
tline = fgetl(fid);
for i = 1:18
    tline = fgetl(fid);
    evalc(tline);
end

subset = 1:length(x1);
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

% Initial geometry plot
figure;
patch(transpose(X), transpose(Y), transpose(Z), 'g', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.3);
title('ProtoMPEX Geometry');
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
set(gca, 'fontsize', 13);

% --- Fixed Z value --- %
z_fixed = 1.0;  % ← Set this to desired Z manually

% --- Sweep over Y space --- %
y_space = linspace(-1.14, 1.14, 100);  % More points for better resolution

% Store results: surface index and hit location
intersected_indices = NaN(1, length(y_space));
hit_points = NaN(length(y_space), 3);

p0 = [0 0 0];  % ray origin

for yi = 1:length(y_space)
    p1 = [8.25, y_space(yi), z_fixed];  % sweeping in Y at fixed Z

    intersected_surfaces = [];
    ps = [];

    for i = 1:length(x1)
        [did_hit, p] = surf_intersect(p0, p1, ...
            x1(i), y1(i), z1(i), x2(i), y2(i), z2(i), ...
            x3(i), y3(i), z3(i), a(i), b(i), c(i), d(i), plane_norm(i));

        if did_hit
            intersected_surfaces = [intersected_surfaces; i]
            ps = [ps; p];
        end
    end

    if ~isempty(intersected_surfaces)
        distance = vecnorm(ps, 2, 2);
        [~, index] = min(distance);
        intersected_indices(yi) = intersected_surfaces(index);
        hit_points(yi, :) = ps(index, :);
    end
end

% Plot intersection indices vs Y
figure;
plot(y_space, intersected_indices, 'o-');
xlabel('Y [m]');
ylabel('Intersected Surface Index');
title(['Intersected Surface at Z = ', num2str(z_fixed)]);
grid on;

% Optionally save
% csvwrite('intersected_indices_vs_y.csv', [y_space(:), intersected_indices(:)]);

% ========== Subfunction ========== %
function [did_hit, p] = surf_intersect(p0, p1, ...
    x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d, plane_norm)

    did_hit = 0;
    p = [NaN NaN NaN];

    d0 = (a * p0(1) + b * p0(2) + c * p0(3) + d) / plane_norm;
    d1 = (a * p1(1) + b * p1(2) + c * p1(3) + d) / plane_norm;

    if sign(d0) ~= sign(d1)
        t = -d0 / (a * (p1(1) - p0(1)) + b * (p1(2) - p0(2)) + c * (p1(3) - p0(3)));
        p = p0 + t * (p1 - p0);

        A = [x1, y1, z1];
        B = [x2, y2, z2];
        C = [x3, y3, z3];

        v0 = B - A;
        v1 = C - A;
        v2 = p - A;

        dot00 = dot(v0, v0);
        dot01 = dot(v0, v1);
        dot02 = dot(v0, v2);
        dot11 = dot(v1, v1);
        dot12 = dot(v1, v2);

        denom = dot00 * dot11 - dot01 * dot01;
        if denom == 0
            return;
        end

        u = (dot11 * dot02 - dot01 * dot12) / denom;
        v = (dot00 * dot12 - dot01 * dot02) / denom;

        if u >= 0 && v >= 0 && (u + v) <= 1
            did_hit = 1;
        end
    end
end