% close all;
% clear all;
%% Loop through all variables in the NetCDF file and plot
file = 'profilesWEST.nc';
% Read coordinate arrays
R = ncread(file, 'x'); R_ant_index = find(R>3.01);
z = ncread(file, 'z'); z_mid_index = find(z>0);

% Normalize density to match ne_lim of the RF simulation
ne = ncread(file, 'ne') ;
ne_lim_SOLEDGE = ne(min(R_ant_index),min(z_mid_index)) ;
ne_lim_COMSOL = 7.3e18;    %1.1E17;  % 7.3e18;
ne_GITR = ne * ne_lim_COMSOL / ne_lim_SOLEDGE ;

% Plot 2D density normalized for GITR simu
figure;
imagesc(R, z, ne_GITR');
set(gca, 'YDir', 'normal', 'FontName', 'Times', 'FontSize', 16);
xlabel('$R$ [m]', 'Interpreter','latex','FontSize',16);
ylabel('$Z$ [m]', 'Interpreter','latex','FontSize',16);
% Fixed title with proper LaTeX subscript formatting
title(['Normalized $n_e$ for GITR : $n_{e,\mathrm{lim}}$=', num2str(ne_lim_COMSOL, '%.1e')], 'Interpreter','latex');
colorbar;

%% Create new NetCDF file with modified ne data
output_file = sprintf('profilesWEST_ne_lim=%.1e.nc', ne_lim_COMSOL);

% Get all variable information from original file
info = ncinfo(file);

% Create new file (delete if exists)
if exist(output_file, 'file')
    delete(output_file);
end

% First pass: Create all dimensions
fprintf('Creating dimensions...\n');
for i = 1:length(info.Dimensions)
    dim = info.Dimensions(i);
    fprintf('  Dimension: %s (length: %d)\n', dim.Name, dim.Length);
end

% Second pass: Create all variables with their dimensions
fprintf('Creating variables...\n');
for i = 1:length(info.Variables)
    var = info.Variables(i);
    varName = var.Name;
    fprintf('  Variable: %s', varName);
    
    % Prepare dimension information
    if ~isempty(var.Dimensions)
        dimNames = {var.Dimensions.Name};
        dimLengths = [var.Dimensions.Length];
        fprintf(' (dimensions: %s)\n', strjoin(dimNames, ' x '));
        
        % Create variable with proper dimensions - fix the syntax
        dimSpec = cell(1, 2*length(dimNames));
        for j = 1:length(dimNames)
            dimSpec{2*j-1} = dimNames{j};
            dimSpec{2*j} = dimLengths(j);
        end
        
        nccreate(output_file, varName, ...
            'Dimensions', dimSpec, ...
            'Datatype', var.Datatype, ...
            'Format', 'netcdf4');
    else
        % Create scalar variable
        fprintf(' (scalar)\n');
        nccreate(output_file, varName, ...
            'Datatype', var.Datatype, ...
            'Format', 'netcdf4');
    end
end

% Third pass: Write data and attributes
fprintf('Writing data and attributes...\n');
for i = 1:length(info.Variables)
    var = info.Variables(i);
    varName = var.Name;
    
    % Read original data or use modified data
    if strcmp(varName, 'ne')
        % Use the normalized ne_GITR data instead
        original_ne = ncread(file, varName);
        data = ne_GITR;
        
        % Ensure dimensions match exactly
        if ~isequal(size(data), size(original_ne))
            fprintf('  Warning: Dimension mismatch detected for ne variable\n');
            fprintf('    Original ne size: [%s]\n', num2str(size(original_ne)));
            fprintf('    ne_GITR size: [%s]\n', num2str(size(data)));
            fprintf('    Reshaping ne_GITR to match original dimensions\n');
            data = reshape(data, size(original_ne));
        end
        
        fprintf('  Using normalized ne_GITR data for variable: %s (size: [%s])\n', ...
            varName, num2str(size(data)));
    else
        data = ncread(file, varName);
    end
    
    % Write data
    try
        ncwrite(output_file, varName, data);
    catch ME
        fprintf('  Error writing variable %s: %s\n', varName, ME.message);
        fprintf('  Variable info - Size: [%s], Class: %s\n', ...
            num2str(size(data)), class(data));
        rethrow(ME);
    end
    
    % Copy variable attributes
    for j = 1:length(var.Attributes)
        attr = var.Attributes(j);
        try
            ncwriteatt(output_file, varName, attr.Name, attr.Value);
        catch ME
            fprintf('  Warning: Could not write attribute %s for variable %s: %s\n', ...
                attr.Name, varName, ME.message);
        end
    end
end

% Copy global attributes
fprintf('Copying global attributes...\n');
for i = 1:length(info.Attributes)
    attr = info.Attributes(i);
    try
        ncwriteatt(output_file, '/', attr.Name, attr.Value);
    catch ME
        fprintf('  Warning: Could not write global attribute %s: %s\n', ...
            attr.Name, ME.message);
    end
end

fprintf('Created new NetCDF file: %s\n', output_file);
fprintf('The ''ne'' variable has been replaced with normalized ne_GITR data\n');

%% Plot all variables from original file
% Get list of variables
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
    imagesc(R, z, data');
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