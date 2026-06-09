ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;


%% Read profiles from SOLEDGE for selected density case

densityCase = 'high';
% Options:
%   'low'   -> profilesWEST_ne_lim=1.1e+17.nc
%   'high'  -> profilesWEST_ne_lim=7.3e+18.nc
%   'both'  -> both cases
%   numeric -> e.g. densityCase = 1.1e17;

switch lower(string(densityCase))

    case "low"
        profileFiles = {'profilesWEST_ne_lim=1.1e+17.nc'};
        caseNames    = {'ne_lim_1p1e17'};

    case "high"
        profileFiles = {'profilesWEST_ne_lim=7.3e+18.nc'};
        caseNames    = {'ne_lim_7p3e18'};

    case "both"
        profileFiles = { ...
            'profilesWEST_ne_lim=1.1e+17.nc', ...
            'profilesWEST_ne_lim=7.3e+18.nc'};
        caseNames = { ...
            'ne_lim_1p1e17', ...
            'ne_lim_7p3e18'};

    otherwise
        if isnumeric(densityCase)
            ne_lim = densityCase;
        else
            ne_lim = str2double(densityCase);
        end

        if ~isfinite(ne_lim)
            error('Invalid densityCase. Use "low", "high", "both", or numeric density.');
        end

        profileFiles = {sprintf('profilesWEST_ne_lim=%.1e.nc', ne_lim)};
        caseNames    = {sprintf('ne_lim_%s', density_tag(ne_lim))};
end

prof = struct([]);

for icase = 1:numel(profileFiles)

    ncfile   = profileFiles{icase};
    caseName = caseNames{icase};

    fprintf('\n>>>> Reading %s\n', ncfile)

    %% Grid
    R = ncread(ncfile,'x');
    z = ncread(ncfile,'z');

    %% B-field
    br = ncread(ncfile,'br');
    bt = ncread(ncfile,'bt');
    bz = ncread(ncfile,'bz');

    %% Densities
    ne = ncread(ncfile,'ne');
    ni = ncread(ncfile,'ni');

    %% Temperatures
    te = ncread(ncfile,'te');
    ti = ncread(ncfile,'ti');

    %% Velocity
    vt = ncread(ncfile,'vt');
    vr = ncread(ncfile,'vr');
    vz = ncread(ncfile,'vz');

    %% Clean NaN / Inf
    br(~isfinite(br)) = 0;
    bt(~isfinite(bt)) = 0;
    bz(~isfinite(bz)) = 0;

    ne(~isfinite(ne)) = 0;
    ni(~isfinite(ni)) = 0;

    te(~isfinite(te)) = 0;
    ti(~isfinite(ti)) = 0;

    vr(~isfinite(vr)) = 0;
    vt(~isfinite(vt)) = 0;
    vz(~isfinite(vz)) = 0;

    %% Centroids
    r_centroid = sqrt(centroid(:,1).^2 + centroid(:,2).^2);
    z_centroid = centroid(:,3);

    %% Interpolate profiles onto geometry
    ne_surf = interpn(R,z,ne,r_centroid,z_centroid,'linear',0);
    ni_surf = interpn(R,z,ni,r_centroid,z_centroid,'linear',0);

    te_surf = interpn(R,z,te,r_centroid,z_centroid,'linear',0);
    ti_surf = interpn(R,z,ti,r_centroid,z_centroid,'linear',0);

    br_surf = interpn(R,z,br,r_centroid,z_centroid,'linear',0);
    bt_surf = interpn(R,z,bt,r_centroid,z_centroid,'linear',0);
    bz_surf = interpn(R,z,bz,r_centroid,z_centroid,'linear',0);

    vr_surf = interpn(R,z,vr,r_centroid,z_centroid,'linear',0);
    vt_surf = interpn(R,z,vt,r_centroid,z_centroid,'linear',0);
    vz_surf = interpn(R,z,vz,r_centroid,z_centroid,'linear',0);

    %% Clean interpolated fields
    ne_surf(~isfinite(ne_surf)) = 0;
    ni_surf(~isfinite(ni_surf)) = 0;

    te_surf(~isfinite(te_surf)) = 0;
    ti_surf(~isfinite(ti_surf)) = 0;

    br_surf(~isfinite(br_surf)) = 0;
    bt_surf(~isfinite(bt_surf)) = 0;
    bz_surf(~isfinite(bz_surf)) = 0;

    vr_surf(~isfinite(vr_surf)) = 0;
    vt_surf(~isfinite(vt_surf)) = 0;
    vz_surf(~isfinite(vz_surf)) = 0;

    %% Convert cylindrical B to Cartesian B
    phi_centroid = atan2(centroid(:,2),centroid(:,1));

    bx = double(br_surf .* cos(phi_centroid) - bt_surf .* sin(phi_centroid));
    by = double(br_surf .* sin(phi_centroid) + bt_surf .* cos(phi_centroid));
    bz_cart = double(bz_surf);

    bx(~isfinite(bx)) = 0;
    by(~isfinite(by)) = 0;
    bz_cart(~isfinite(bz_cart)) = 0;

    b_mag = sqrt(bx.^2 + by.^2 + bz_cart.^2);
    b_mag(~isfinite(b_mag)) = 0;

    ubx = zeros(size(b_mag));
    uby = zeros(size(b_mag));
    ubz = zeros(size(b_mag));

    idxB = b_mag > 0;
    ubx(idxB) = bx(idxB) ./ b_mag(idxB);
    uby(idxB) = by(idxB) ./ b_mag(idxB);
    ubz(idxB) = bz_cart(idxB) ./ b_mag(idxB);

    ubx(~isfinite(ubx)) = 0;
    uby(~isfinite(uby)) = 0;
    ubz(~isfinite(ubz)) = 0;

    %% Convert cylindrical velocity to Cartesian velocity
    vx = double(vr_surf .* cos(phi_centroid) - vt_surf .* sin(phi_centroid));
    vy = double(vr_surf .* sin(phi_centroid) + vt_surf .* cos(phi_centroid));
    vz_cart = double(vz_surf);

    vx(~isfinite(vx)) = 0;
    vy(~isfinite(vy)) = 0;
    vz_cart(~isfinite(vz_cart)) = 0;

    v_mag = sqrt(vx.^2 + vy.^2 + vz_cart.^2);
    v_mag(~isfinite(v_mag)) = 0;

    %% Flux on surface
    o8plus_flux_surf = ni_surf .* v_mag;
    o8plus_flux_surf(~isfinite(o8plus_flux_surf)) = 0;

    %% Surface-normal angle to B
    norm_vec_mag = sqrt(norm_vec(:,1).^2 + norm_vec(:,2).^2 + norm_vec(:,3).^2);
    norm_vec_mag(~isfinite(norm_vec_mag)) = 0;

    unorm_vec = zeros(size(norm_vec));
    idxN = norm_vec_mag > 0;

    unorm_vec(idxN,1) = norm_vec(idxN,1) ./ norm_vec_mag(idxN);
    unorm_vec(idxN,2) = norm_vec(idxN,2) ./ norm_vec_mag(idxN);
    unorm_vec(idxN,3) = norm_vec(idxN,3) ./ norm_vec_mag(idxN);

    unorm_vec(~isfinite(unorm_vec)) = 0;

    cosTheta = unorm_vec(:,1).*ubx + unorm_vec(:,2).*uby + unorm_vec(:,3).*ubz;
    cosTheta(~isfinite(cosTheta)) = 0;
    cosTheta = max(min(cosTheta,1),-1);

    theta = abs(cosTheta);
    theta(~isfinite(theta)) = 0;

    % %% Save case-specific outputs
    % writematrix(ne_surf, sprintf('ne_surf_%s.csv', caseName));
    % writematrix(ni_surf, sprintf('ni_surf_%s.csv', caseName));
    % writematrix(te_surf, sprintf('te_surf_%s.csv', caseName));
    % writematrix(ti_surf, sprintf('ti_surf_%s.csv', caseName));
    % 
    % writematrix(br_surf, sprintf('br_surf_%s.csv', caseName));
    % writematrix(bt_surf, sprintf('bt_surf_%s.csv', caseName));
    % writematrix(bz_surf, sprintf('bz_surf_%s.csv', caseName));
    % 
    % writematrix(vr_surf, sprintf('vr_surf_%s.csv', caseName));
    % writematrix(vt_surf, sprintf('vt_surf_%s.csv', caseName));
    % writematrix(vz_surf, sprintf('vz_surf_%s.csv', caseName));
    % 
    % writematrix(b_mag, sprintf('b_mag_%s.csv', caseName));
    % writematrix(theta, sprintf('theta_%s.csv', caseName));
    % writematrix(o8plus_flux_surf, sprintf('o8plus_flux_surf_%s.csv', caseName));

    %% Store outputs
    prof(icase).caseName = caseName;
    prof(icase).ncfile = ncfile;

    prof(icase).ne_surf = ne_surf;
    prof(icase).ni_surf = ni_surf;
    prof(icase).te_surf = te_surf;
    prof(icase).ti_surf = ti_surf;

    prof(icase).br_surf = br_surf;
    prof(icase).bt_surf = bt_surf;
    prof(icase).bz_surf = bz_surf;

    prof(icase).vr_surf = vr_surf;
    prof(icase).vt_surf = vt_surf;
    prof(icase).vz_surf = vz_surf;

    prof(icase).bx = bx;
    prof(icase).by = by;
    prof(icase).bz_cart = bz_cart;
    prof(icase).b_mag = b_mag;

    prof(icase).vx = vx;
    prof(icase).vy = vy;
    prof(icase).vz_cart = vz_cart;
    prof(icase).v_mag = v_mag;

    prof(icase).theta = theta;
    prof(icase).o8plus_flux_surf = o8plus_flux_surf;

    %% Quick plots
    figure;
    histogram(theta);
    title(['theta = |cos(angle normal,B)|: ', caseName], 'Interpreter','none');
    xlabel('|cos(theta)|');
    ylabel('Counts');

    figure;
    quiver3(centroid(:,1),centroid(:,2),centroid(:,3), ...
            norm_vec(:,1),norm_vec(:,2),norm_vec(:,3));
    hold on;
    quiver3(centroid(:,1),centroid(:,2),centroid(:,3), ...
            bx,by,bz_cart);
    title(['Surface normals and B field: ', caseName], 'Interpreter','none');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;

end

%% Select default output variables for downstream scripts
% If one case is selected, expose variables directly.
% If both cases are selected, expose the high-density case by default.

if numel(prof) == 1
    useCase = 1;
else
    useCase = find(strcmp({prof.caseName}, 'ne_lim_7p3e18'), 1);
    if isempty(useCase)
        useCase = numel(prof);
    end
end

caseName = prof(useCase).caseName;

ne_surf = prof(useCase).ne_surf;
ni_surf = prof(useCase).ni_surf;
te_surf = prof(useCase).te_surf;
ti_surf = prof(useCase).ti_surf;

br_surf = prof(useCase).br_surf;
bt_surf = prof(useCase).bt_surf;
bz_surf = prof(useCase).bz_surf;

vr_surf = prof(useCase).vr_surf;
vt_surf = prof(useCase).vt_surf;
vz_surf = prof(useCase).vz_surf;

b_mag = prof(useCase).b_mag;
theta = prof(useCase).theta;
o8plus_flux_surf = prof(useCase).o8plus_flux_surf;

fprintf('\n>>>> Downstream variables set from case: %s\n', caseName);

%% Local helper
function tag = density_tag(ne_lim)
    tag = sprintf('%.1e', ne_lim);
    tag = strrep(tag, '.', 'p');
    tag = strrep(tag, '+', '');
    tag = strrep(tag, '-', 'm');
end