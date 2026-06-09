clc
clear
close all
format long
file_path="../diiid-helicon/DIII-D_helicon_runs/"
file_gitrm="file0.txt"
aa=dlmread(file_path+file_gitrm);
v_ptcl_sq=aa(:,4).^2 + aa(:,5).^2 + aa(:,6).^2;

%% 
% Constants and convert to energy
m_i = 12*1.66e-27; % Ion mass (e.g., proton)
e = 1.602e-19;  % Elementary charge
%E_eV = v_ptcl_sq;
E_eV = (0.5 * m_i ).* v_ptcl_sq./e;

% Input
xb = 7.41;     % Carbon
xc = 100;      % Cut-off energy (eV)
x = linspace(0, 300, 1000); % Energy axis
% Compute Thompson distribution
f = x ./ (x + xb).^3 .* (1 - sqrt((x + xb) / (xc + xb)));
% Set f(x > xc) = 0
f(x > xc) = 0;
f(f < 0) = 0;  % Clean up negative values due to rounding
% Normalize to get PDF
f_norm = f / (sum(f) * [x(2)-x(1)]);

% Plotting
nbins = 100;
[counts, centers] = hist(E_eV, nbins);
bin_width = centers(2) - centers(1);
pdf_values = counts / (sum(counts));

figure
subplot(1,2,1)
plot(x, f_norm, 'b-', 'LineWidth', 2)
xlabel('Energy (eV)')
ylabel('PDF')
title('Thompson Distribution Carbon-Analytical/Input')
grid on
xlim([0, 50])

subplot(1,2,2)
plot(centers, pdf_values, 'b-', 'LineWidth', 2)
xlabel('Energy (eV)')
ylabel('PDF')
title('Thompson Distribution of Carbon-Output')
grid on
xlim([0, 50])

