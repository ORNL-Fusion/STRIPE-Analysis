clc
clear
close all
format long

% ---- File path ----
file_gitrm = "../tilted_targets/test/0_degrees/MPEX_runs/file0.txt";

% ---- Load data ----
data = load(file_gitrm);

% ---- Extract quantities ----
vx = data(:,4);
vy = data(:,5);
vz = data(:,6);

% Last column: charge state distribution
charge_state = data(:,end);

% ---- Expected MPEX impurity source parameters ----
v_parallel_expected = 1000;   % m/s
Tz_eV = 5;                    % eV

% ---- Plot velocity histograms separately ----
figure('Color','w');

subplot(3,1,1)
histogram(vx,50)
xlabel('v_x [m/s]')
ylabel('Counts')
title('Velocity Distribution: v_x')
xlim([-4e4 4e4])
grid on

subplot(3,1,2)
histogram(vy,50)
xlabel('v_y [m/s]')
ylabel('Counts')
title('Velocity Distribution: v_y')
xlim([-4e4 4e4])
grid on

subplot(3,1,3)
histogram(vz,50)
hold on
xline(v_parallel_expected,'r--','Expected v_{\parallel,Z}=1000 m/s', ...
    'LineWidth',1.5,'LabelOrientation','horizontal')
xlabel('v_z = v_{\parallel} [m/s]')
ylabel('Counts')
title('Velocity Distribution: v_z')
xlim([-4e4 4e4])
grid on

sgtitle('GITRm MPEX Ta Impurity Velocity Distributions')

% ---- Combined normalized velocity PDF ----
bins = 100;

figure('Color','w')
histogram(vx,bins,'Normalization','pdf')
hold on
histogram(vy,bins,'Normalization','pdf')
histogram(vz,bins,'Normalization','pdf')
xline(v_parallel_expected,'r--','Expected v_{\parallel,Z}=1000 m/s', ...
    'LineWidth',1.5,'LabelOrientation','horizontal')

legend('v_x','v_y','v_z','Expected v_{\parallel,Z}', ...
    'Location','best')
xlabel('Velocity [m/s]')
ylabel('PDF')
title('Normalized Velocity Distribution')
xlim([-4e4 4e4])
grid on

% ---- Charge state distribution ----
figure('Color','w')
histogram(charge_state,'BinMethod','integers','Normalization','probability')
xlabel('Charge state')
ylabel('Fraction')
title('Ta Charge State Distribution from Last Column')
grid on

% ---- Print charge-state fractions ----
fprintf('\nCharge state distribution from last column:\n')
unique_Z = unique(charge_state);

for i = 1:length(unique_Z)
    frac = sum(charge_state == unique_Z(i)) / length(charge_state);
    fprintf('Ta%d+ : %.2f %%\n', unique_Z(i), 100*frac)
end

fprintf('\nExpected model charge state distribution:\n')
fprintf('Ta2+ : 8 %%\n')
fprintf('Ta3+ : 62 %%\n')
fprintf('Ta4+ : 30 %%\n')
fprintf('\nExpected v_parallel,Z = %.1f m/s\n', v_parallel_expected)