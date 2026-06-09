%% Loop through all variables in the NetCDF file and plot

file = 'profilesWEST.nc';

% Read coordinate arrays
R = ncread(file, 'x');
Z = ncread(file, 'z');

% Get list of variables
info = ncinfo(file);
variables = {info.Variables.Name};

% Skip coordinates
skipVars = {'x','z'};

for i = 1:length(variables)
    varName = variables{i};
    if any(strcmp(varName, skipVars))
        continue
    end

    % Read data
    data = ncread(file, varName);
    data(isnan(data)) = 0;

    % Plot 2D image
    figure;
    imagesc(R, Z, data');
    set(gca, 'YDir', 'normal', 'FontName', 'Times', 'FontSize', 16);
    xlabel('$R$ [m]', 'Interpreter','latex','FontSize',16);
    ylabel('$Z$ [m]', 'Interpreter','latex','FontSize',16);
    title(['Input ', strrep(varName,'_','\_')], 'Interpreter','latex');
    colorbar;

    % Also plot a 1D slice for quick inspection
    figure;
    plot(R, data(:, round(end/2)));
    xlabel('$R$ [m]', 'Interpreter','latex','FontSize',16);
    ylabel(strrep(varName,'_','\_'), 'Interpreter','latex','FontSize',16);
    title(['Radial Profile of ', strrep(varName,'_','\_')], 'Interpreter','latex');
    grid on;
end