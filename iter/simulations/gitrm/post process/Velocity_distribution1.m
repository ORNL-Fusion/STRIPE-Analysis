clc
clear
close all
format long
% file_gitrm="../code_tests_t20000/file0.txt";
file_gitrm="../ITER_antenna_runs/file0.txt";

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
xlim([-1e4 1e4])

subplot(3,1,2)
histogram(vy,50)
xlabel('v_y')
ylabel('Counts')
title('Velocity Distribution v_y')
xlim([-1e4 1e4])

subplot(3,1,3)
histogram(vz,50)
xlabel('v_z')
ylabel('Counts')
title('Velocity Distribution v_z')
xlim([-1e4 1e4])

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