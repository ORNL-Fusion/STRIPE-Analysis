clc; clear; close all;

%% Files
finalfile = '../code_tests_t500_noIon/finalPos0.txt';
initfile  = '../code_tests_t500_noIon/file0.txt';

% finalfile = '../ITER_antenna_runs/finalPos0.txt';
% initfile  = '../ITER_antenna_runs/file0.txt';

%% Read data
Mf = readmatrix(finalfile,'NumHeaderLines',1);
Mi = readmatrix(initfile ,'NumHeaderLines',1);

% Final positions
x = Mf(:,1);
y = Mf(:,2);
z = Mf(:,3);
charge = Mf(:,8);
max(charge)

% Initial positions
x0 = Mi(:,1);
y0 = Mi(:,2);
z0 = Mi(:,3);

% Ensure same length
N = min(length(x),length(x0));
x = x(1:N); y = y(1:N); z = z(1:N);
x0 = x0(1:N); y0 = y0(1:N); z0 = z0(1:N);
charge = charge(1:N);

%% Remove invalid rows
valid = isfinite(x) & isfinite(y) & isfinite(z) & ...
        isfinite(x0) & isfinite(y0) & isfinite(z0) & ...
        isfinite(charge);

x = x(valid); y = y(valid); z = z(valid);
x0 = x0(valid); y0 = y0(valid); z0 = z0(valid);
charge = charge(valid);

%% Charge-state histogram (unchanged)
states = unique(charge);
edges_hist = [states; max(states)+1];
counts_hist = histcounts(charge, edges_hist);

figure
bar(states, counts_hist)
xlabel('Charge State')
ylabel('Counts')
title('Charge State Distribution')
grid on

%% Radial displacement from initial position
r = sqrt((x-x0).^2 + (y-y0).^2)+ (z-z0).^2;

%% Average charge vs displacement radius
edges = linspace(min(r), max(r), 40);
bin = discretize(r, edges);
rc = 0.5*(edges(1:end-1)+edges(2:end));

good = ~isnan(bin);

counts_r = accumarray(bin(good),1,[numel(rc),1],@sum,0);
charge_sum = accumarray(bin(good),charge(good),[numel(rc),1],@sum,0);

qavg = charge_sum ./ counts_r;
qavg(counts_r==0) = NaN;

figure
plot(rc,qavg,'-o')
xlabel('Radial Displacement  sqrt((x-x_0)^2+(y-y_0)^2)')
ylabel('Average Charge State')
title('Average Charge vs Radial Displacement')
grid on

%% Charge vs displacement histogram
figure
histogram2(r,charge,40,'DisplayStyle','tile','ShowEmptyBins','on')
xlabel('Radial Displacement')
ylabel('Charge State')
title('Charge Distribution vs Radial Displacement')
colorbar
grid onis 