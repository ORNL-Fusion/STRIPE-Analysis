clc;
clear;
close all;
format longE;

filename = "profilesWEST_ne_lim_1p1e17.nc";
x = ncread(filename, 'x');
z = ncread(filename, 'z');
nX = length(x);
nZ = length(z);

zero_data = zeros(nZ, nX);  % (z, x) order
newVars = {'Er', 'Et', 'Ez', 'gradTir', 'gradTit', 'gradTiz', ...
           'gradTer', 'gradTet', 'gradTez'};

for i = 1:length(newVars)
    varName = newVars{i};
    
    % Create variable with correct dimension order
    nccreate(filename, varName, 'Dimensions', {'nZ', nZ, 'nX', nX});
    
    % Write zero data
    ncwrite(filename, varName, zero_data);
end