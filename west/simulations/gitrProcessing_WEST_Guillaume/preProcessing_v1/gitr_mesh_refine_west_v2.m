%% WEST limiter geometry + surface fields + erosion-based particle source
close all;
clear all;
clc;

%% ======================= GEOMETRY ============================
tic
planes = planes_from_stl('limiters_west_comsol.stl',0.05);

planes_copy = planes;

% Adjust coordinates to geometry convention
planes(:,2) = planes_copy(:,3);
planes(:,5) = planes_copy(:,6);
planes(:,8) = planes_copy(:,9);
planes(:,3) = planes_copy(:,2);
planes(:,6) = planes_copy(:,5);
planes(:,9) = planes_copy(:,8);

refinementnum = 2;
toc

%% ======================= INITIAL GEOMETRY ARRAYS =============
nFaces = length(planes);

abcd       = zeros(nFaces,4);
area       = zeros(nFaces,1);
centroid   = zeros(nFaces,3);
plane_norm = zeros(nFaces,1);
BCxBA      = zeros(nFaces,1);
CAxCB      = zeros(nFaces,1);

for i = 1:nFaces
    A = planes(i,1:3);
    B = planes(i,4:6);
    C = planes(i,7:9);

    AB = B - A;
    AC = C - A;
    BC = C - B;
    BA = A - B;
    CA = -AC;
    CB = -BC;

    area(i) = 0.5 * norm(cross(AB,AC));

    normalVec = cross(AB,AC);
    normalVec = normalVec ./ norm(normalVec);

    d = -dot(normalVec,A);

    abcd(i,:)    = [normalVec,d];
    plane_norm(i)= norm(normalVec);
    BCxBA(i)     = sign(dot(cross(BC,BA),normalVec));
    CAxCB(i)     = sign(dot(cross(CA,CB),normalVec));

    centroid(i,:) = mean([A;B;C],1);
end

%% ======================= REFINE HELICON REGION ===============
y_faces = find(abs(centroid(:,3)-0.16) <= 2.94);

ret = plot_panes(planes,y_faces);

planes_refined = refine_planes(planes(y_faces,:),refinementnum);

ret = plot_panes(planes_refined,1:length(planes_refined));

% Remove old planes and add refined planes
planes(y_faces,:) = [];
planes = [planes; planes_refined];

ret = plot_panes(planes,1:length(planes));

%% ======================= RECALCULATE GEOMETRY ================
nFaces = length(planes);

abcd       = zeros(nFaces,4);
area       = zeros(nFaces,1);
centroid   = zeros(nFaces,3);
plane_norm = zeros(nFaces,1);
BCxBA      = zeros(nFaces,1);
CAxCB      = zeros(nFaces,1);

for i = 1:nFaces
    A = planes(i,1:3);
    B = planes(i,4:6);
    C = planes(i,7:9);

    AB = B - A;
    AC = C - A;
    BC = C - B;
    BA = A - B;
    CA = -AC;
    CB = -BC;

    area(i) = 0.5 * norm(cross(AB,AC));

    normalVec = cross(AB,AC);
    normalVec = normalVec ./ norm(normalVec);

    d = -dot(normalVec,A);

    abcd(i,:)     = [normalVec,d];
    plane_norm(i) = norm(normalVec);
    BCxBA(i)      = sign(dot(cross(BC,BA),normalVec));
    CAxCB(i)      = sign(dot(cross(CA,CB),normalVec));

    centroid(i,:) = mean([A;B;C],1);
end

area(~isfinite(area)) = 0;

%% ======================= SURFACE FLAGS =======================
materialZ = zeros(nFaces,1);
surfs     = zeros(nFaces,1);

HeliconRefined = find(abs(centroid(:,3)-0.16) <= 2.94 & ...
                      centroid(:,1) <= 0.3583 & ...
                      centroid(:,1) >= -0.3583);

materialZ(HeliconRefined) = 13;
surfs(HeliconRefined)     = 1;

surf_inds = find(surfs);

%% ======================= PLOT GEOMETRY =======================
plotSet = 1:nFaces;

X = [planes(plotSet,1), planes(plotSet,4), planes(plotSet,7)];
Y = [planes(plotSet,2), planes(plotSet,5), planes(plotSet,8)];
Z = [planes(plotSet,3), planes(plotSet,6), planes(plotSet,9)];

figure;
patch(X.',Y.',Z.',[0.7 0.7 0.7], ...
    'FaceAlpha',0.3,'EdgeColor','none');
title('GITR Geometry');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
axis equal tight;
hold on;

%% ======================= NORMAL DIRECTIONS ===================
inDir = ones(1,nFaces);
norm_vec = zeros(nFaces,3);

for i = 1:nFaces
    normal = -abcd(i,1:3) ./ plane_norm(i);
    normal(~isfinite(normal)) = 0;

    l_normal = 0.01;
    normal = l_normal .* normal;

    dot_product = dot(normal,centroid(i,:));

    if sign(dot_product) == 1
        inDir(i) = -1;
    end

    normal = inDir(i) .* normal;
    norm_vec(i,:) = normal;
end

%% ======================= READ PLASMA SURFACE PROFILES =========
read_profiles_west_GU;
antenna_data_WEST;

%% ======================= SHEATH CALCULATION ==================
disp('>>>> Calculating surface potentials')

sheathType = 2;

switch sheathType

    case 1
        disp('>>>>>> Using thermal sheath model')

        me = 1/2000;
        background_amu = 2;

        sheath_factor = abs(0.5 .* log((2*pi*me/background_amu) .* ...
                         (1 + ti_surf ./ te_surf)));

        potential_surf = sheath_factor .* te_surf;

    case 2
        disp('>>>>>> Using COMSOL potential map')

        comsol_data = readmatrix('../Guillaume_cases/ne_lim_7p3e18/VDC_out.txt');
        % comsol_data = readmatrix('../Guillaume_cases/ne_lim_1p1e17/VDC_out.txt');

        ne_g = ne_surf;
        te_g = te_surf;

        k = 1.38e-23 * 11604;
        c_bar = sqrt(8 .* k .* te_g ./ pi ./ 4 ./ (12 .* 1.66e-27));

        xx = comsol_data(:,1);
        zz = comsol_data(:,2);
        yy = comsol_data(:,3);

        potential = comsol_data(:,4);

        figure;
        scatter3(xx,yy,zz,[],potential);
        title('COMSOL potential map');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        colorbar;

        potential_surf = zeros(nFaces,1);

        for i = 1:nFaces
            distance = sqrt((centroid(i,1)-xx).^2 + ...
                            (centroid(i,2)-yy).^2 + ...
                            (centroid(i,3)-zz).^2);

            [~,I] = min(distance);
            potential_surf(i) = potential(I);
        end

    case 3
        disp('>>>>>> Using Petra-M potential map')

        petraM_data = readmatrix('sheath.dat');

        xx = petraM_data(:,3);
        yy = -petraM_data(:,4);
        zz = petraM_data(:,5);

        potential = petraM_data(:,6);

        potential_surf = zeros(nFaces,1);

        for i = 1:nFaces
            distance = sqrt((centroid(i,1)-xx).^2 + ...
                            (centroid(i,2)-yy).^2 + ...
                            (centroid(i,3)-zz).^2);

            [~,I] = min(distance);
            potential_surf(i) = potential(I);
        end
end

potential_surf = potential_surf(:);
potential_surf(~isfinite(potential_surf)) = 0;

figure;
patch(X.',Y.',Z.',potential_surf, ...
    'FaceAlpha',1,'EdgeColor','none');
title('Surface potential');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal tight;
colorbar;

%% ======================= LOAD PRECOMPUTED YIELDS / FLUX ======
centroid_x = centroid(:,1);
centroid_y = centroid(:,2);
centroid_z = centroid(:,3);

yields_total = readmatrix('../ieads_guillaume/ne_lim_7p3e18/yields/west_total_yield_sum.txt');
erosion_flux_total = readmatrix('../ieads_guillaume/ne_lim_7p3e18/yields/west_total_erosion_flux.txt');

% yields_total = readmatrix('../ieads_guillaume/ne_lim_1p1e17/yields/west_total_yield_sum.txt');
% erosion_flux_total = readmatrix('../ieads_guillaume/ne_lim_1p1e17/yields/west_total_erosion_flux.txt');


yields_total = pad_or_trim_local(yields_total, nFaces);
erosion_flux_total = pad_or_trim_local(erosion_flux_total, nFaces);

yields_total(~isfinite(yields_total)) = 0;
erosion_flux_total(~isfinite(erosion_flux_total)) = 0;

%% ======================= WRITE SURFACE FIELD FILES ============
surf_fields = table( ...
    centroid_x, ...
    centroid_y, ...
    centroid_z, ...
    potential_surf, ...
    ne_surf(:), ...
    v_mag(:), ...
    yields_total(:), ...
    erosion_flux_total(:));

writetable(surf_fields, 'surf_fields.txt', ...
    'WriteVariableNames', false, ...
    'Delimiter', ' ');

surface_variables = table( ...
    potential_surf(:), ...
    ne_surf(:), ...
    te_surf(:), ...
    te_surf(:), ...
    v_mag(:), ...
    b_mag(:), ...
    bz_surf(:), ...
    br_surf(:), ...
    bt_surf(:), ...
    theta(:), ...
    ni_surf(:));

writetable(surface_variables,'Targets.txt');

%% ======================= WRITE GITR GEOMETRY =================
disp('>>>> Generating GITR geometry')

fileID = fopen('gitrGeometryPointPlane3d.cfg','w');

fprintf(fileID,'geom = \n{ \n   x1 = [');
write_array_cfg(fileID,planes(:,1),'%5e');

fprintf(fileID,' ] \n   y1 = [');
write_array_cfg(fileID,planes(:,2),'%5e');

fprintf(fileID,' ] \n   z1 = [');
write_array_cfg(fileID,planes(:,3),'%5e');

fprintf(fileID,' ] \n   x2 = [');
write_array_cfg(fileID,planes(:,4),'%5e');

fprintf(fileID,' ] \n   y2 = [');
write_array_cfg(fileID,planes(:,5),'%5e');

fprintf(fileID,' ] \n   z2 = [');
write_array_cfg(fileID,planes(:,6),'%5e');

fprintf(fileID,' ] \n   x3 = [');
write_array_cfg(fileID,planes(:,7),'%5e');

fprintf(fileID,' ] \n   y3 = [');
write_array_cfg(fileID,planes(:,8),'%5e');

fprintf(fileID,' ] \n   z3 = [');
write_array_cfg(fileID,planes(:,9),'%5e');

fprintf(fileID,' ] \n   a = [');
write_array_cfg(fileID,abcd(:,1),'%5e');

fprintf(fileID,' ] \n   b = [');
write_array_cfg(fileID,abcd(:,2),'%5e');

fprintf(fileID,' ] \n   c = [');
write_array_cfg(fileID,abcd(:,3),'%5e');

fprintf(fileID,' ] \n   d = [');
write_array_cfg(fileID,abcd(:,4),'%5e');

fprintf(fileID,' ] \n   plane_norm = [');
write_array_cfg(fileID,plane_norm,'%5e');

fprintf(fileID,' ] \n   BCxBA = [');
write_array_cfg(fileID,BCxBA,'%5e');

fprintf(fileID,' ] \n   CAxCB = [');
write_array_cfg(fileID,CAxCB,'%5e');

fprintf(fileID,' ] \n   area = [');
write_array_cfg(fileID,area,'%5e');

fprintf(fileID,' ] \n   Z = [');
write_array_cfg(fileID,materialZ,'%f');

fprintf(fileID,' ] \n   surface = [');
write_array_cfg(fileID,surfs,'%i');

fprintf(fileID,' ] \n   inDir = [');
write_array_cfg(fileID,inDir(:),'%i');

fprintf(fileID,' ] \n   potential = [');
write_array_cfg(fileID,potential_surf,'%5e');

fprintf(fileID,' ] \n');
fprintf(fileID,'periodic = 0;\n');
fprintf(fileID,'theta0 = 0.0;\n');
fprintf(fileID,'theta1 = 0.0\n');
fprintf(fileID,'periodic_bc_x0 = 0.0;\n');
fprintf(fileID,'periodic_bc_x1 = 0.0;\n');
fprintf(fileID,'periodic_bc_x = 0;}\n');

fclose(fileID);

%% ======================= EROSION RATE FROM PRECOMPUTED FLUX ===
disp('>>>> Initializing particles from precomputed erosion flux')

nP = 1E6;

amu = 1.66053906660e-27;   % kg
amu_eroded = 183.84;       % W target mass; change if eroded target is not W
mass_per_particle = amu_eroded * amu;

% erosion_flux_total units assumed: particles/m^2/s
eroded_flux = erosion_flux_total(:);
eroded_flux(~isfinite(eroded_flux)) = 0;

area = pad_or_trim_local(area, nFaces);
area(~isfinite(area)) = 0;

% Area-integrated erosion per face
erosion_all_faces = abs(eroded_flux) .* area;  % particles/s
erosion_all_faces(~isfinite(erosion_all_faces)) = 0;

% Restrict to active source surfaces
erosion_surface = erosion_all_faces(surf_inds);

erosion_rate = sum(erosion_surface,'omitnan');        % particles/s
mass_erosion_rate = erosion_rate * mass_per_particle; % kg/s

fprintf('\n=== Area-integrated erosion from precomputed flux ===\n');
fprintf('Total geometry area       = %.6e m^2\n', sum(area,'omitnan'));
fprintf('Active source area        = %.6e m^2\n', sum(area(surf_inds),'omitnan'));
fprintf('Total erosion rate        = %.6e particles/s\n', erosion_rate);
fprintf('Total mass erosion rate   = %.6e kg/s\n', mass_erosion_rate);
fprintf('Total mass erosion rate   = %.6e g/s\n', mass_erosion_rate*1e3);

%% ======================= CDF FOR PARTICLE SOURCE ==============
erosion_inds = find(erosion_surface > 0);
erosion_sub  = erosion_surface(erosion_inds);

if isempty(erosion_sub)
    error('erosion_flux_total is zero everywhere on active surfaces. Cannot sample particle source.');
end

erosion_sub_cdf = cumsum(erosion_sub);
erosion_sub_cdf = erosion_sub_cdf ./ erosion_sub_cdf(end);

writematrix(eroded_flux,        'source_eroded_flux_particles_m2_s.txt');
writematrix(erosion_all_faces,  'source_erosion_particles_s_all_faces.txt');
writematrix(erosion_surface,    'source_erosion_particles_s_active_surface.txt');
writematrix(erosion_inds,       'source_erosion_indices.txt');
writematrix(erosion_sub,        'source_erosion_cdf_values.txt');
writematrix(erosion_sub_cdf,    'source_erosion_cdf.txt');
writematrix(erosion_rate,       'source_total_erosion_rate_particles_s.txt');
writematrix(mass_erosion_rate,  'source_total_mass_erosion_rate_kg_s.txt');

figure;
plot(erosion_sub_cdf,'LineWidth',2);
grid on;
xlabel('Active eroding surface index');
ylabel('Cumulative erosion fraction');
title('CDF for erosion-based particle source');

%% ======================= SAMPLE PARTICLES ====================
rand1 = rand(nP,1);

element_ceil = discretize(rand1,[0; erosion_sub_cdf(:)]);
element_ceil(isnan(element_ceil)) = length(erosion_sub_cdf);

xP  = zeros(1,nP);
yP  = zeros(1,nP);
zP  = zeros(1,nP);
vxP = zeros(1,nP);
vyP = zeros(1,nP);
vzP = zeros(1,nP);

offset = 1e-5;
v_source = 5000; % m/s

for j = 1:nP

    iFace = surf_inds(erosion_inds(element_ceil(j)));

    normal = -abcd(iFace,1:3) ./ plane_norm(iFace);
    normal = inDir(iFace) .* normal;
    normal(~isfinite(normal)) = 0;

    x_tri = X(iFace,:) + offset .* normal(1);
    y_tri = Y(iFace,:) + offset .* normal(2);
    z_tri = Z(iFace,:) + offset .* normal(3);

    samples = sample_triangle(x_tri,y_tri,z_tri,1);

    xP(j) = samples(1,1);
    yP(j) = samples(1,2);
    zP(j) = samples(1,3);

    vxP(j) = v_source .* normal(1);
    vyP(j) = v_source .* normal(2);
    vzP(j) = v_source .* normal(3);
end

xP(~isfinite(xP))   = 0;
yP(~isfinite(yP))   = 0;
zP(~isfinite(zP))   = 0;
vxP(~isfinite(vxP)) = 0;
vyP(~isfinite(vyP)) = 0;
vzP(~isfinite(vzP)) = 0;

figure;
patch(X.',Y.',Z.','b','FaceAlpha',0.3,'EdgeColor','none');
hold on;
scatter3(xP,yP,zP,2,'r','filled');
title('Particle source sampled from precomputed erosion flux');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
axis equal tight;

%% ======================= WRITE PARTICLE SOURCE ===============
outfile = './particle_source_west.nc';

if isfile(outfile)
    delete(outfile);
end

ncid = netcdf.create(outfile,'NC_WRITE');

dimP  = netcdf.defDim(ncid,'nP',nP);
xVar  = netcdf.defVar(ncid,'x', 'double',dimP);
yVar  = netcdf.defVar(ncid,'y', 'double',dimP);
zVar  = netcdf.defVar(ncid,'z', 'double',dimP);
vxVar = netcdf.defVar(ncid,'vx','double',dimP);
vyVar = netcdf.defVar(ncid,'vy','double',dimP);
vzVar = netcdf.defVar(ncid,'vz','double',dimP);

netcdf.endDef(ncid);

netcdf.putVar(ncid,xVar, xP);
netcdf.putVar(ncid,yVar, yP);
netcdf.putVar(ncid,zVar, zP);
netcdf.putVar(ncid,vxVar,vxP);
netcdf.putVar(ncid,vyVar,vyP);
netcdf.putVar(ncid,vzVar,vzP);

netcdf.close(ncid);

fprintf('\nParticle source written to %s\n', outfile);

%% ======================= FUNCTIONS ===========================
function ret = plot_panes(planes0,plotSet)

X = [planes0(plotSet,1), planes0(plotSet,4), planes0(plotSet,7)];
Y = [planes0(plotSet,2), planes0(plotSet,5), planes0(plotSet,8)];
Z = [planes0(plotSet,3), planes0(plotSet,6), planes0(plotSet,9)];

figure;
patch(X.',Y.',Z.',[0.7 0.7 0.7], ...
    'FaceAlpha',0.3,'EdgeColor','k');

title('GITR Geometry');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
axis equal tight;

ret = 1;
end

function planes = planes_from_stl(filename,hmax)

model = createpde;
importGeometry(model,filename);

figure;
pdegplot(model,'FaceLabels','on');

mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);

figure;
pdeplot3D(model,'FaceAlpha',0.5);

[p,~,t] = meshToPet(mesh);

tess = transpose(t(1:4,:));

faces = [ ...
    tess(:,[1 2 3]); ...
    tess(:,[1 2 4]); ...
    tess(:,[1 3 4]); ...
    tess(:,[2 3 4])];

faces = sort(faces,2);
faces = sortrows(faces);

Y = diff(faces);
zeroRow = [0,0,0];

k = ismember(Y,zeroRow,'rows');
k2 = find(k ~= 0);

faces([k2; k2+1],:) = [];

C = faces;

planes = zeros(length(C),9);
planes(:,:) = [ ...
    transpose(p(1:3,C(:,1))), ...
    transpose(p(1:3,C(:,2))), ...
    transpose(p(1:3,C(:,3)))];

end

function planes_refined = refine_planes(planes,n)

X = [planes(:,1), planes(:,4), planes(:,7)];
Y = [planes(:,2), planes(:,5), planes(:,8)];
Z = [planes(:,3), planes(:,6), planes(:,9)];

[Xrefined,Yrefined,Zrefined] = refineXYZ(X,Y,Z,n);

planes_refined = [ ...
    Xrefined(:,1), Yrefined(:,1), Zrefined(:,1), ...
    Xrefined(:,2), Yrefined(:,2), Zrefined(:,2), ...
    Xrefined(:,3), Yrefined(:,3), Zrefined(:,3)];

end

function [Xrefined,Yrefined,Zrefined] = refineXYZ(X,Y,Z,n)

for j = 1:n

    Xrefined = zeros(2*length(X),3);
    Yrefined = zeros(2*length(X),3);
    Zrefined = zeros(2*length(X),3);

    for i = 1:length(X)

        A = [X(i,1), Y(i,1), Z(i,1)];
        B = [X(i,2), Y(i,2), Z(i,2)];
        C = [X(i,3), Y(i,3), Z(i,3)];

        AB = B - A;
        AC = C - A;
        BC = C - B;

        norms = [norm(AB), norm(BC), norm(AC)];
        [~,maxInd] = max(norms);

        if maxInd == 1
            M = A + 0.5*AB;

            tri1 = [A;M;C];
            tri2 = [M;B;C];

        elseif maxInd == 2
            M = B + 0.5*BC;

            tri1 = [A;B;M];
            tri2 = [A;M;C];

        else
            M = A + 0.5*AC;

            tri1 = [A;B;M];
            tri2 = [M;B;C];
        end

        Xrefined(2*i-1,:) = tri1(:,1).';
        Yrefined(2*i-1,:) = tri1(:,2).';
        Zrefined(2*i-1,:) = tri1(:,3).';

        Xrefined(2*i,:) = tri2(:,1).';
        Yrefined(2*i,:) = tri2(:,2).';
        Zrefined(2*i,:) = tri2(:,3).';
    end

    X = Xrefined;
    Y = Yrefined;
    Z = Zrefined;
end

end

function samples = sample_triangle(x,y,z,nP)

x_transform = x - x(1);
y_transform = y - y(1);
z_transform = z - z(1);

v1 = [x_transform(2), y_transform(2), z_transform(2)];
v2 = [x_transform(3), y_transform(3), z_transform(3)];

a1 = rand(nP,1);
a2 = rand(nP,1);

flipInd = a1 + a2 > 1;
a1(flipInd) = 1 - a1(flipInd);
a2(flipInd) = 1 - a2(flipInd);

samples = a1.*v1 + a2.*v2;

samples(:,1) = samples(:,1) + x(1);
samples(:,2) = samples(:,2) + y(1);
samples(:,3) = samples(:,3) + z(1);

end

function v = pad_or_trim_local(v,N)

v = double(v(:));
v(~isfinite(v)) = 0;

if numel(v) > N
    v = v(1:N);
elseif numel(v) < N
    v = [v; zeros(N-numel(v),1)];
end

end

function write_array_cfg(fileID,arr,fmt)

arr = arr(:);

fprintf(fileID,fmt,arr(1));

for i = 2:length(arr)
    fprintf(fileID,',');
    fprintf(fileID,fmt,arr(i));
end

end