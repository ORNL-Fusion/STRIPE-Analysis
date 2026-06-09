close all;
clear all;
clc;

data = readmatrix('Targets_D+.txt');
potential_data = data(:,1);
ne_data = data(:,11);
te_data = data(:,4);
v_data = data(:,5);

% Prefer IEAD yields from MAT output; fall back to CSV if needed.
matFile = 'ieads_ITER_Ne+.mat';
if ~isfile(matFile)
    candidates = dir('ieads_ITER_*.mat');
    if ~isempty(candidates)
        [~, idx] = max([candidates.datenum]);
        matFile = candidates(idx).name;
    end
end

ieads_data = [];
ieads_label = 'IEADS';
charge_state = 3;
if isfile(matFile)
    charge_match = regexp(matFile, 'Ne(\d+)\+', 'tokens', 'once');
    if ~isempty(charge_match)
        charge_state = str2double(charge_match{1});
    end
end
if isfile(matFile)
    S = load(matFile);
    if ~isfield(S, 'yields')
        error('Variable "yields" was not found in %s.', matFile);
    end
    yields_data = S.yields;
    if ~isvector(yields_data)
        yields_data = yields_data(:,1);
    end

    if isfield(S, 'ieads')
        ieads_data = reduce_ieads_to_scalar(S.ieads);
    elseif isfield(S, 'IEADS')
        ieads_data = reduce_ieads_to_scalar(S.IEADS);
    end
else
    yields_csv = readmatrix('yields_D+.csv');
    yields_data = yields_csv(:,1);
end

yields_data = yields_data(:);
yields_data(isnan(yields_data)) = 0;
yields_data = [0; yields_data];

if isempty(ieads_data)
    ieads_data = build_ieads_from_surface_ncs(numel(yields_data)-1, charge_state);
    if ~isempty(ieads_data)
        ieads_label = sprintf('Integrated surfEDist from surface_C%d_loc_*.nc', charge_state);
    end
end

if isempty(ieads_data)
    error(['IEADS not found. Need either (1) ieads/IEADS variable in MAT, ' ...
           'or (2) readable surface_C*_loc_*.nc files with surfEDist.']);
end
ieads_data = ieads_data(:);
ieads_data(isnan(ieads_data)) = 0;
ieads_data = [0; ieads_data];

n = min([numel(yields_data), numel(ieads_data), numel(ne_data), numel(v_data)]);
yields_data = yields_data(1:n);
ieads_data = ieads_data(1:n);
potential_data = potential_data(1:n);
ne_data = ne_data(1:n);
te_data = te_data(1:n);
v_data = v_data(1:n);
ero_data = yields_data .* ne_data .* v_data;

if (exist('x1') == 0)
    fid = fopen(strcat('gitrGeometryPointPlane3d.cfg'));
    if fid < 0
        error('Could not open gitrGeometryPointPlane3d.cfg');
    end
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    fclose(fid);
    Zsurface = Z;
end

subset = 1:length(x1); % find(r<0.07 & z1> 0.001 & z1 < .20);
% subset = find(r<0.049 & z1 > -0.001 & z1<0.001)

% Match geometry face count to available data.
subset = subset(1:min(numel(subset), n));

figure(1)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),yields_data(1:numel(subset)),'FaceAlpha',1,'EdgeAlpha', 0.3)
title('Yields')
colorbar('eastoutside')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

figure(2)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),ne_data(1:numel(subset)),'FaceAlpha',1,'EdgeAlpha', 0.3)
title('Density')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(3)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),te_data(1:numel(subset)),'FaceAlpha',1,'EdgeAlpha', 0.3)
title('Temperature')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(4)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),ero_data(1:numel(subset)),'FaceAlpha',1,'EdgeAlpha', 0.3)
title('Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(5)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),ieads_data(1:numel(subset)),'FaceAlpha',1,'EdgeAlpha', 0.1)
title({'IEADs on Geometry', ieads_label}, 'Interpreter', 'none')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')
set(gca,'FontSize',14)

if exist('area','var') == 1
    erosion_rate = sum(sum(ero_data(1:numel(subset)).*area(1:numel(subset))));
else
    erosion_rate = sum(ero_data(1:numel(subset)));
end
disp(['erosion_rate = ', num2str(erosion_rate)]);

function out = reduce_ieads_to_scalar(in)
sz = size(in);
if isvector(in)
    out = in(:);
elseif numel(sz) == 2
    if sz(1) >= sz(2)
        out = sum(in, 2);
    else
        out = sum(in, 1).';
    end
else
    nLoc = sz(1);
    tmp = reshape(in, nLoc, []);
    out = sum(tmp, 2);
end
end

function ieads = build_ieads_from_surface_ncs(nLoc, charge_state)
ieads = nan(nLoc,1);
for i = 1:nLoc
    filename = sprintf('surface_C%d_loc_%d.nc', charge_state, i-1);
    filename_out = fullfile('output', filename);
    if isfile(filename)
        read_file = filename;
    elseif isfile(filename_out)
        read_file = filename_out;
    else
        continue
    end
    try
        surfEDist = ncread(read_file, 'surfEDist');
    catch
        continue
    end
    surfEDist(isnan(surfEDist)) = 0;
    ieads(i) = sum(surfEDist(:));
end

valid = ~isnan(ieads);
if ~any(valid)
    ieads = [];
else
    ieads(~valid) = 0;
end
end
