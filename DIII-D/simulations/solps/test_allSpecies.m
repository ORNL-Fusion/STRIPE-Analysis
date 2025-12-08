%% Quick script: plot density for species 4 from interpolated_values.nc
clear; clc;

fname = 'interpolated_values.nc';

% Load grid
r = ncread(fname,'gridr');   % [nr]
z = ncread(fname,'gridz');   % [nz]

% Load species densities
ni = ncread(fname,'ni');  % (nz,nr,ns)

% Choose species index (e.g. 4)
% ispec = 4;
% ni_spec = squeeze(ni_all(:,:,ispec));   % (nz,nr)

% Quick plot
figure;
pcolor(r, z, log10(ni')); shading flat; colorbar;
xlabel('R [m]'); ylabel('Z [m]');
% title(['log_{10} n_i for species ' num2str(ispec)]);