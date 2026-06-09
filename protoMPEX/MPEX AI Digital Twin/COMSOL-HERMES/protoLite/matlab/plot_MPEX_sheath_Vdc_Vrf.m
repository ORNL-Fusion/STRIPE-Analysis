%% plot_MPEX_phi0avg_Vdc_Vrf_with_sheath_width.m
% Computes V_RF, phi0avg, and physical V_DC using local Te from
% window_parameters.csv.
%
% Efield.csv columns:
%   x, y, z, Ex, Ey
%
% window_parameters.csv columns:
%   x, y, z, ne, Br, B0, Te[J]
%
% Important:
%   V_RF = |E_normal| * d_sheath
%   phi0avg is NOT multiplied by Te
%   V_DC_volts = Te_eV * phi0avg

clc; clear; close all;

%% ----------------------------
% USER INPUT
%% ----------------------------
efieldFile = 'Efield.csv';
windowFile = 'window_parameters.csv';

% Use earlier/effective sheath width here
% 0.002 m gives ~60 V for your current Efield file
% 0.010 m gives ~300 V
d_sheath = 0.002;     % [m]

% If you want to force constant Te instead of local Te, set true
USE_CONSTANT_TE = true;
Te_constant_eV  = 8.0;

% Use median density for weak omega_RF/omega_pi correction
USE_MEDIAN_DENSITY_FOR_PHI0AVG = true;

SMOOTH_MAPS  = true;
SMOOTH_SIGMA = 2.0;

e_charge = 1.602176634e-19;

%% ----------------------------
% READ DATA
%% ----------------------------
LayerData = read_complex_csv(efieldFile);
EdgeData  = readmatrix(windowFile,'CommentStyle','%');

x  = real(LayerData(:,1));
y  = real(LayerData(:,2));
z  = real(LayerData(:,3));

Ex = LayerData(:,4);
Ey = LayerData(:,5);

Density = EdgeData(:,4);   % [m^-3]
Br      = EdgeData(:,5);
B0      = EdgeData(:,6);

% Te is column 7 in joules in your file
Te_local_eV = EdgeData(:,7) ./ e_charge;

N = min([numel(x),numel(Density),numel(Te_local_eV)]);

x = x(1:N);
y = y(1:N);
z = z(1:N);
Ex = Ex(1:N);
Ey = Ey(1:N);

Density = Density(1:N);
Br = Br(1:N);
B0 = B0(1:N);
Te_local_eV = Te_local_eV(1:N);

if USE_CONSTANT_TE
    Te_eV = Te_constant_eV .* ones(size(Te_local_eV));
else
    Te_eV = Te_local_eV;
end

%% ----------------------------
% CLEAN DATA
%% ----------------------------
r = sqrt(x.^2 + y.^2);

good = r > 0 & isfinite(x) & isfinite(y) & isfinite(z) & ...
       isfinite(real(Ex)) & isfinite(real(Ey)) & ...
       isfinite(Density) & isfinite(Br) & isfinite(B0) & ...
       isfinite(Te_eV) & Te_eV > 0;

x = x(good);
y = y(good);
z = z(good);
Ex = Ex(good);
Ey = Ey(good);
Density = Density(good);
Br = Br(good);
B0 = B0(good);
Te_eV = Te_eV(good);
r = r(good);

phi_deg = atan2(y,x) * 180/pi;

%% ----------------------------
% RF SHEATH VOLTAGE
%% ----------------------------
E_norm = (Ex.*x + Ey.*y)./r;

V_RF = abs(E_norm) .* d_sheath;     % [V]
xi   = V_RF ./ Te_eV;               % dimensionless

%% ----------------------------
% PHI0AVG AND V_DC
%% ----------------------------
if USE_MEDIAN_DENSITY_FOR_PHI0AVG
    Density_used = median(Density,'omitnan') .* ones(size(Density));
else
    Density_used = Density;
end

phi0avg = compute_phi0avg_getstyle(V_RF,Density_used,Te_eV);

% Physical DC voltage in volts
V_DC_volts = Te_eV .* phi0avg;

ratio_RF_to_VDC = V_RF ./ V_DC_volts;

%% ----------------------------
% SAVE OUTPUT
%% ----------------------------
out = table(phi_deg(:),z(:),Te_eV(:),V_RF(:),xi(:), ...
    phi0avg(:),V_DC_volts(:),ratio_RF_to_VDC(:), ...
    Density(:),Density_used(:),Br(:),B0(:), ...
    'VariableNames',{'phi_deg','z','Te_eV','V_RF_volts', ...
    'xi_VRF_over_Te','phi0avg_not_times_Te', ...
    'V_DC_volts_Te_times_phi0avg','V_RF_over_V_DC', ...
    'Density_original','Density_used','Br','B0'});

writetable(out,'sheath_phi0avg_Vdc_Vrf_with_sheath_width.csv');

fprintf('Saved: sheath_phi0avg_Vdc_Vrf_with_sheath_width.csv\n');
fprintf('d_sheath = %.4g m\n',d_sheath);
fprintf('Te_eV range = %.3g to %.3g eV\n',min(Te_eV),max(Te_eV));
fprintf('V_RF range [V] = %.3g to %.3g\n',min(V_RF),max(V_RF));
fprintf('xi = V_RF/Te range = %.3g to %.3g\n',min(xi),max(xi));
fprintf('phi0avg range [not * Te] = %.3g to %.3g\n',min(phi0avg),max(phi0avg));
fprintf('V_DC range [V] = %.3g to %.3g\n',min(V_DC_volts),max(V_DC_volts));

%% ----------------------------
% SANITY PLOTS
%% ----------------------------
figure('Color','w','Position',[100 100 760 640])
scatter(V_RF,phi0avg,12,'filled','MarkerFaceAlpha',0.35)
xlabel('V_{RF} [V]','FontSize',15,'FontWeight','bold')
ylabel('\phi_{0,avg} [not multiplied by T_e]','FontSize',15,'FontWeight','bold')
title('\phi_{0,avg} versus V_{RF}','FontSize',17,'FontWeight','bold')
grid on; box on
set(gca,'FontSize',13,'LineWidth',1.2)

figure('Color','w','Position',[150 150 760 640])
scatter(V_RF,V_DC_volts,12,'filled','MarkerFaceAlpha',0.35)
xlabel('V_{RF} [V]','FontSize',15,'FontWeight','bold')
ylabel('V_{DC}=T_e\phi_{0,avg} [V]','FontSize',15,'FontWeight','bold')
title('Physical V_{DC} versus V_{RF}','FontSize',17,'FontWeight','bold')
grid on; box on
set(gca,'FontSize',13,'LineWidth',1.2)

%% ----------------------------
% MAPS
%% ----------------------------
make_interpolated_map(phi_deg,z,V_RF, ...
    'Azimuthal coordinate \phi [deg]','Axial coordinate z [m]', ...
    'RF Sheath Voltage V_{RF}(\phi,z)', ...
    'V_{RF} [V]',SMOOTH_MAPS,SMOOTH_SIGMA);

make_interpolated_map(phi_deg,z,phi0avg, ...
    'Azimuthal coordinate \phi [deg]','Axial coordinate z [m]', ...
    '\phi_{0,avg}(\phi,z), not multiplied by T_e', ...
    '\phi_{0,avg}',SMOOTH_MAPS,SMOOTH_SIGMA);

make_interpolated_map(phi_deg,z,V_DC_volts, ...
    'Azimuthal coordinate \phi [deg]','Axial coordinate z [m]', ...
    'Physical DC Sheath Voltage V_{DC}(\phi,z)', ...
    'V_{DC} [V]',SMOOTH_MAPS,SMOOTH_SIGMA);

make_interpolated_map(phi_deg,z,ratio_RF_to_VDC, ...
    'Azimuthal coordinate \phi [deg]','Axial coordinate z [m]', ...
    'V_{RF}/V_{DC}(\phi,z)', ...
    'V_{RF}/V_{DC}',SMOOTH_MAPS,SMOOTH_SIGMA);

%% ----------------------------
% HISTOGRAMS
%% ----------------------------
figure('Color','w','Position',[250 250 850 600])
histogram(V_RF,60)
xlabel('V_{RF} [V]','FontSize',14,'FontWeight','bold')
ylabel('Counts','FontSize',14,'FontWeight','bold')
title('Distribution of V_{RF}','FontSize',16,'FontWeight','bold')
grid on; box on

figure('Color','w','Position',[300 300 850 600])
histogram(V_DC_volts,60)
xlabel('V_{DC} [V]','FontSize',14,'FontWeight','bold')
ylabel('Counts','FontSize',14,'FontWeight','bold')
title('Distribution of Physical V_{DC}','FontSize',16,'FontWeight','bold')
grid on; box on

fprintf('Plotting complete.\n');

%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

function phi0avg = compute_phi0avg_getstyle(Vlayer,Density,Te_eV)

    mu = 24.17;
    Z  = 1;
    A  = 2;

    w_rf = 13.56e6 * 2*pi;

    omegapi = 1.32e3 * Z .* sqrt(Density*1e-6 ./ A);
    omega_hat = w_rf ./ omegapi;

    xi = Vlayer ./ Te_eV;

    j = 0;
    upar0 = 1.1;

    a1 = 3.70285;
    a2 = 3.81991;
    b1 = 1.13352;
    b2 = 1.24171;
    a3 = 2.0*b2/pi;

    c0 = 0.966463;
    c1 = 0.141639;

    gg = c0 + c1*tanh(omega_hat);
    xi1 = gg .* xi;

    ff = ((log(mu) + xi1*a1 + xi1.^2*a2 + xi1.^3*a3) ./ ...
         (1 + xi1*b1 + xi1.^2*b2)) ...
         - log(1 - (j/upar0)) + log(mu/24.17);

    phi0avg = real(ff);
end

function data = read_complex_csv(filename)

    fid = fopen(filename,'r');
    if fid < 0
        error('Could not open file: %s',filename);
    end

    rows = {};

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        if isempty(line) || startsWith(line,'%')
            continue
        end

        parts = split(line,',');

        vals = zeros(1,numel(parts));

        for k = 1:numel(parts)
            s = strtrim(parts{k});
            s = strrep(s,'i','j');

            if contains(s,'j')
                vals(k) = str2num(s); %#ok<ST2NM>
            else
                vals(k) = str2double(s);
            end
        end

        rows{end+1,1} = vals; %#ok<AGROW>
    end

    fclose(fid);

    nrow = numel(rows);
    ncol = max(cellfun(@numel,rows));

    data = nan(nrow,ncol);

    for i = 1:nrow
        vals = rows{i};
        data(i,1:numel(vals)) = vals;
    end
end

function make_interpolated_map(x,y,val,xlab,ylab,titleStr,cbarLabel,doSmooth,sigma)

    good = isfinite(x) & isfinite(y) & isfinite(val);

    x = x(good);
    y = y(good);
    val = val(good);

    nx = 500;
    ny = 350;

    xq = linspace(-180,180,nx);
    yq = linspace(min(y),max(y),ny);

    [Xq,Yq] = meshgrid(xq,yq);

    Vq = griddata(x,y,val,Xq,Yq,'linear');
    Vn = griddata(x,y,val,Xq,Yq,'nearest');

    missing = isnan(Vq);
    Vq(missing) = Vn(missing);

    if doSmooth
        Vq = smooth2a_simple(Vq,sigma);
    end

    figure('Color','w','Position',[300 150 950 650])

    imagesc(xq,yq,Vq)
    set(gca,'YDir','normal')
    colormap(turbo)

    h = colorbar;
    ylabel(h,cbarLabel,'FontSize',13,'FontWeight','bold')

    xlabel(xlab,'FontSize',14,'FontWeight','bold')
    ylabel(ylab,'FontSize',14,'FontWeight','bold')
    title(titleStr,'FontSize',16,'FontWeight','bold')

    xlim([-180 180])
    set(gca,'FontSize',13,'LineWidth',1.2)
    box on
end

function B = smooth2a_simple(A,sigma)

    if sigma <= 0
        B = A;
        return
    end

    ksize = max(3,2*ceil(3*sigma)+1);
    half = floor(ksize/2);

    [xg,yg] = meshgrid(-half:half,-half:half);
    G = exp(-(xg.^2 + yg.^2)/(2*sigma^2));
    G = G ./ sum(G(:));

    mask = isfinite(A);
    A0 = A;
    A0(~mask) = 0;

    W = conv2(double(mask),G,'same');
    B = conv2(A0,G,'same') ./ W;
    B(W <= 0) = NaN;
end