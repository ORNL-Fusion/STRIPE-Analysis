clc
clear
close all
format long

%% ===================== USER INPUTS =====================

run_path = "/Users/78k/ORNL Dropbox/Atul Kumar/work/STRIPE-Analysis/iter/simulations/IMAS-123361/123361/run";

species_name = "NE";        % Neutral neon in Case.Tri.data.atm
plot_log10   = true;

marker_size  = 24;
color_limits = [10 20];     % [] for automatic color scale

axis_limits  = [];
% Example:
% axis_limits = [4 9 -5 5];

make_grid_plot = false;

fprintf("Working on run_path = %s\n", run_path)

if ~isfolder(run_path)
    error("run_path is not a directory: %s", run_path)
end

%% ===================== LOAD SOLPS CASE =====================

tic
Case = load_solps_case(run_path);
toc

%% ===================== CHECK REQUIRED DATA =====================

if ~isfield(Case, "Tri") || isempty(Case.Tri)
    error("Case.Tri is missing or empty.")
end

if ~isfield(Case.Tri, "data") || isempty(Case.Tri.data)
    error("Case.Tri.data is missing or empty.")
end

if ~isfield(Case.Tri.data, "pdena")
    error("Case.Tri.data.pdena not found.")
end

if ~isfield(Case.Tri.data, "atm")
    error("Case.Tri.data.atm not found. Cannot identify neutral atomic species.")
end

if ~isfield(Case, "Geo") || isempty(Case.Geo)
    error("Case.Geo is missing or empty.")
end

if ~isfield(Case.Geo, "cvX") || ~isfield(Case.Geo, "cvY")
    error("Case.Geo.cvX/cvY not found. Cannot use B2 cell-center coordinates.")
end

%% ===================== IDENTIFY NEUTRAL SPECIES INDEX =====================

atm_list = string(Case.Tri.data.atm);

fprintf("\nEIRENE atomic neutral species:\n")
for i = 1:numel(atm_list)
    fprintf("  index %d : %s\n", i, atm_list(i))
end

iAtom = find(strcmpi(atm_list, species_name), 1);

if isempty(iAtom)
    error("Could not find neutral species '%s' in Case.Tri.data.atm.", species_name)
end

fprintf("\nUsing neutral species %s with atomic index iAtom = %d\n", ...
    species_name, iAtom)

if isfield(Case, "Ion") && isfield(Case.Ion, "label") && ...
        isfield(Case.Ion, "indEireneAtom")

    fprintf("\nIon-to-EIRENE atom mapping:\n")
    for i = 1:numel(Case.Ion.label)
        if Case.Ion.indEireneAtom(i) > 0
            fprintf("  %-8s -> atom index %d\n", ...
                string(Case.Ion.label{i}), Case.Ion.indEireneAtom(i))
        end
    end

end

%% ===================== EXTRACT NEUTRAL DENSITY =====================

pdena = Case.Tri.data.pdena;

if iAtom > size(pdena,2)
    error("Atomic index iAtom=%d exceeds pdena species count=%d.", ...
        iAtom, size(pdena,2))
end

n0_all = pdena(:, iAtom);

fprintf("\nsize(Case.Tri.data.pdena) = [%s]\n", num2str(size(pdena)))
fprintf("length(Case.Geo.cvX)      = %d\n", numel(Case.Geo.cvX))
fprintf("length(Case.Geo.cvY)      = %d\n", numel(Case.Geo.cvY))

%% ===================== MAP TO B2 CELL CENTERS =====================

R_all = Case.Geo.cvX(:);
Z_all = Case.Geo.cvY(:);

nUse = min([numel(R_all), numel(Z_all), numel(n0_all)]);

if nUse < numel(n0_all)
    warning("Using first %d entries to match available cell-center coordinates.", nUse)
end

R = R_all(1:nUse);
Z = Z_all(1:nUse);
n0 = n0_all(1:nUse);

valid = isfinite(R) & ...
        isfinite(Z) & ...
        isfinite(n0) & ...
        n0 > 0;

R = R(valid);
Z = Z(valid);
n0 = n0(valid);

if isempty(n0)
    error("No positive finite neutral density values found for species %s.", species_name)
end

if plot_log10
    plot_data = log10(n0);
    cbar_label = sprintf("log_{10}(n_{%s^0}) [m^{-3}]", species_name);
else
    plot_data = n0;
    cbar_label = sprintf("n_{%s^0} [m^{-3}]", species_name);
end

%% ===================== 2D CELL-CENTER PLOT =====================

figure("Color","w", ...
       "Position",[100 100 850 1000])

scatter(R, Z, marker_size, plot_data, "filled")

axis equal tight
box on
grid on

if ~isempty(axis_limits)
    axis(axis_limits)
end

colormap(turbo)

if ~isempty(color_limits)
    clim(color_limits)
end

cb = colorbar;
ylabel(cb, cbar_label)

xlabel("R [m]")
ylabel("Z [m]")

title(sprintf("EIRENE Neutral %s Density", species_name))

set(gca, ...
    "FontSize", 16, ...
    "FontWeight", "bold")

%% ===================== PRINT DATA RANGE =====================

fprintf("\nNeutral %s density range:\n", species_name)
fprintf("min n0 = %.4e m^-3\n", min(n0, [], "omitnan"))
fprintf("max n0 = %.4e m^-3\n", max(n0, [], "omitnan"))

if plot_log10
    fprintf("min log10(n0) = %.4f\n", min(plot_data, [], "omitnan"))
    fprintf("max log10(n0) = %.4f\n", max(plot_data, [], "omitnan"))
end

%% ===================== OPTIONAL OMP-LIKE PROFILE =====================

if isfield(Case.Geo, "imapCv") && ~isempty(Case.Geo.imapCv)

    nx_tot = size(Case.Geo.imapCv, 1);
    ny_tot = size(Case.Geo.imapCv, 2);

    ix = 2:nx_tot-1;
    jomp = round(ny_tot/2);

    fprintf("\nUsing approximate OMP index jomp = %d of %d\n", jomp, ny_tot)

    cv_ids = Case.Geo.imapCv(ix, jomp);

    valid_cv = cv_ids > 0 & cv_ids <= numel(Case.Geo.cvX);

    cv_ids = cv_ids(valid_cv);

    R_omp = Case.Geo.cvX(cv_ids);
    Z_omp = Case.Geo.cvY(cv_ids);

    F_n0 = scatteredInterpolant(R, Z, n0, "natural", "none");

    n0_omp = F_n0(R_omp, Z_omp);

    valid_omp = isfinite(R_omp) & ...
                isfinite(n0_omp) & ...
                n0_omp > 0;

    R_omp_plot = R_omp(valid_omp);
    n0_omp_plot = n0_omp(valid_omp);

    [R_omp_plot, sort_idx] = sort(R_omp_plot);
    n0_omp_plot = n0_omp_plot(sort_idx);

    figure("Color","w")

    plot(R_omp_plot, n0_omp_plot, "o-", "LineWidth", 2)

    grid on
    box on

    xlabel("R_{OMP} [m]")
    ylabel(sprintf("n_{%s^0}^{OMP} [m^{-3}]", species_name))
    title(sprintf("Neutral %s Density at Approximate OMP", species_name))

    set(gca, "FontSize", 14, "FontWeight", "bold")

    figure("Color","w")

    plot(R_omp_plot, log10(n0_omp_plot), "o-", "LineWidth", 2)

    grid on
    box on

    xlabel("R_{OMP} [m]")
    ylabel(sprintf("log_{10}(n_{%s^0}^{OMP}) [m^{-3}]", species_name))
    title(sprintf("Neutral %s Density at Approximate OMP", species_name))

    set(gca, "FontSize", 14, "FontWeight", "bold")

    fprintf("\nApproximate OMP neutral %s density:\n", species_name)
    fprintf("min n0_OMP = %.4e m^-3\n", min(n0_omp_plot, [], "omitnan"))
    fprintf("max n0_OMP = %.4e m^-3\n", max(n0_omp_plot, [], "omitnan"))
    fprintf("Number of valid OMP points = %d\n", numel(n0_omp_plot))

else

    warning("Case.Geo.imapCv not available. Skipping OMP profile.")

end

%% ===================== OPTIONAL GRID PLOT =====================

if make_grid_plot
    try
        figure("Color","w")
        plot_SOLPS_grid_nice(run_path)
        title("SOLPS Grid")
    catch
        warning("Could not plot SOLPS grid.")
    end
end

%% ===================== FLUID PLOT ON ACTUAL B2 CELLS =====================

% Use actual B2 unstructured cell polygons
Rpoly = Case.Geo.pr(:,1:nUse);
Zpoly = Case.Geo.pz(:,1:nUse);

plot_data_cell = nan(nUse,1);
plot_data_cell(valid) = plot_data;

figure("Color","w", ...
    "Position",[100 100 850 1000])

patch( ...
    "XData", Rpoly, ...
    "YData", Zpoly, ...
    "CData", plot_data_cell.', ...
    "FaceColor", "flat", ...
    "EdgeColor", "none");

axis equal tight
box on
grid on

if ~isempty(axis_limits)
    axis(axis_limits)
end

colormap(turbo)

if ~isempty(color_limits)
    clim(color_limits)
end

cb = colorbar;
ylabel(cb, cbar_label)

xlabel("R [m]")
ylabel("Z [m]")

title(sprintf("EIRENE Neutral %s Density", species_name))

set(gca, ...
    "FontSize", 16, ...
    "FontWeight", "bold")