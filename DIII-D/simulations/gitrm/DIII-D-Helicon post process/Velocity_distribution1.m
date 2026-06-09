clc
clear
close all
format long
file_gitrm="../diiid-helicon/DIII-D_helicon_runs_200882/file0.txt";
% file_gitrm="../ITER_antenna_runs/file0.txt";
% Load the data
data = load(file_gitrm);

% Extract velocity components
vx = data(:,4);
vy = data(:,5);
vz = data(:,6);

% Plot histograms
figure;

subplot(3,1,1)
histogram(vx,50)
xlabel('v_x')
ylabel('Counts')
title('Velocity Distribution v_x')
xlim([-4e4 4e4])

subplot(3,1,2)
histogram(vy,50)
xlabel('v_y')
ylabel('Counts')
title('Velocity Distribution v_y')
xlim([-4e4 4e4])

subplot(3,1,3)
histogram(vz,50)
xlabel('v_z')
ylabel('Counts')
title('Velocity Distribution v_z')
xlim([-4e4 4e4])

bins = 100;

figure
histogram(vx,bins,'Normalization','pdf')
hold on
histogram(vy,bins,'Normalization','pdf')
histogram(vz,bins,'Normalization','pdf')

legend('v_x','v_y','v_z')
xlabel('Velocity')
ylabel('PDF')
title('Velocity Distribution')