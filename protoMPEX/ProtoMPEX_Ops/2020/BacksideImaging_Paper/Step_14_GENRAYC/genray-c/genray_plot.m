%% get genray.nc file
clear
close all
clc

% ne = 2e19:
% filename = 'genray_2300.0_4500.0_1.0_1.0_-30.0_.nc';
% filename = 'genray_2300.0_4500.0_1.0_1.0_-30.0_nocol_.nc';

% ne = 4e19:
% filename = 'genray_2300.0_4500.0_2.0_1.0_-30.0_nocol_.nc';

% Relevant to paper: but wrong coil 7 and 8
% filename = 'genray_3500.0_1630.0_2.0_1.0_-30.0.nc';
% filename = 'genray_3500.0_1630.0_2.0_1.0_-30.0_nocol.nc';

% filename = 'genray_3500.0_2360.0_2.0_1.0_-30.0.nc';
% filename = 'genray_3500.0_2360.0_2.0_1.0_-30.0_nocol.nc';

% filename = 'genray_3500.0_3370.0_2.0_1.0_-30.0.nc';
% filename = 'genray_3500.0_3370.0_2.0_1.0_-30.0_nocol.nc';

% Relevant to paper:
% filename = 'genray_2200.0_1630.0_2.0_1.0_-30.0_.nc';
% filename = 'genray_2200.0_1630.0_2.0_1.0_-30.0_nocol.nc';

% filename = 'genray_2200.0_2360.0_2.0_1.0_-30.0_.nc';
% filename = 'genray_2200.0_2360.0_2.0_1.0_-30.0_nocol.nc';

% filename = 'genray_2200.0_3370.0_2.0_1.0_-30.0_.nc';
% filename = 'genray_2200.0_3370.0_2.0_1.0_-30.0_nocol.nc';

center_distance = 0; %3.22 for Proto-MPEX, 3.4 for ICH
%0-d data
freq =  ncread(filename, 'freqcy');
B_cyl = freq/6.5e9/2; %second harmonic
%1-D data
B_x = ncread(filename, 'eqdsk_x');
B_z = ncread(filename, 'eqdsk_z') + center_distance;
B = ncread(filename, 'bmodprofxz');
ne = ncread(filename, 'densprofxz')*1e6;
ray_x = ncread(filename, 'wx') /100.; %m
ray_y = ncread(filename, 'wy') /100.; %m
ray_z = ncread(filename, 'wz') /100. + center_distance; %m
ray_dist = ncread(filename, 'ws')/100; %m 
ray_pwr = ncread(filename, 'delpwr')*1e-7; %W
ray_Ex = ncread(filename, 'cwexde'); %normalized polarization
ray_Ex = ray_Ex(:,:,1).^2 + ray_Ex(:,:,2).^2;
ray_Ey = ncread(filename, 'cweyde');
ray_Ey = ray_Ey(:,:,1).^2 + ray_Ey(:,:,2).^2;
ray_Ez = ncread(filename, 'cwezde');
ray_Ez = ray_Ez(:,:,1).^2 + ray_Ez(:,:,2).^2;
ray_Bx = ncread(filename, 'sb_x')/1e4;
ray_By = ncread(filename, 'sb_y')/1e4;
ray_Bz = ncread(filename, 'sb_z')/1e4;
ray_B = ncread(filename, 'sbtot')/1e4;
ray_ne = ncread(filename, 'sene')*1e6;
ray_Te = ncread(filename, 'ste')*1e3;
ray_npar = ncread(filename, 'wnpar');
ray_nperp = ncread(filename, 'wnper');
ray_n = sqrt(ray_npar.^2 + ray_nperp.^2);
ray_num = numel(ray_x(1,:));
rho_power = ncread(filename, 'rho_bin_center');
power = ncread(filename, 'powden')*1e-11;
power_e = ncread(filename, 'powden_e')*1e-11;
power_col = ncread(filename, 'powden_cl')*1e-11;
powertot_e = ncread(filename, 'powtot_e')*1e-7; %W
powertol_col = ncread(filename, 'powtot_cl')*1e-7; %W
den_x = ncread(filename,'w_x_densprof_nc');
den_x = den_x(ceil(numel(den_x)/2):end);
den = ncread(filename, 'w_dens_vs_x_nc')*1e6/1e19;
den = den(ceil(numel(den)/2):end);
temp = ncread(filename, 'w_temp_vs_x_nc')*1e3;
temp = temp(ceil(numel(temp)/2):end);
xscan = ncread(filename, 'xscan');
rhoscan = ncread(filename, 'rhoscan');
nscan = numel(rhoscan);
power_x = interp1(rhoscan(1:nscan/2), xscan(1:nscan/2), rho_power);
%2-D data
rgrid = ncread(filename, 'Rgrid');
zgrid = ncread(filename, 'Zgrid') + center_distance;
power_2D = ncread(filename, 'spwr_rz_e');
B_2D = ncread(filename, 'bmodprofxz');
den_1D_l = zeros(530,1);
den_1D_h = den_1D_l;
UH_1D_l = den_1D_h;
UH_1D_h = den_1D_h;
for ii = 1:numel(den_1D_l)
    den_1D_l(ii) = interp1(ne(1:530,ii), B_x(1:530), .98e19);
    den_1D_h(ii) = interp1(ne(531:1060,ii), B_x(531:1060), .98e19);
    UH_1D_l(ii) = interp1(sqrt(81.* ne(1:530,ii) + 784e18* B(1:530,ii).^2), ...
        B_x(1:530), 28e9);
    UH_1D_h(ii) = interp1(sqrt(81.*ne(531:1060,ii) + 784e18* B(531:1060,ii).^2), ...
        B_x(531:1060), 28e9);
end
%write to file
csvwrite(strcat(erase(filename, '.nc'), '.csv'), [ray_x, ray_y, ray_z, ...
    ray_dist, ray_ne, ray_Te ...
    ,ray_pwr, ray_Bx,ray_By,ray_Bz, ray_npar, ray_nperp, ray_Ex, ...
    ray_Ey, ray_Ez]);
% plot genray data
figure
B_arr = B(530,:);
index = find(B_arr > B_cyl);
ax1 = subplot(2,1,1);
title('Magnetic Field (T)');
plot(ax1, B_z, B_arr, 'k', 'LineWidth', 1.5); xlabel(ax1, 'z (m)'); ylabel(ax1, 'B (T)');
ylim([0, max(B_arr) * 1.2]); xlim([center_distance - .4, center_distance + .4]);
%xlim([3.18,3.21]);
test = B_z(min(index));
grid on
%line([test,test], [0, 3.5], 'Color', 'black', 'LineStyle', '--');
%test = B_z(144);
%line([test,test], [0, 2.5], 'Color', 'black', 'LineStyle', '--');
%test = B_z(289);
%line([test,test], [0, 2.5], 'Color', 'black', 'LineStyle', '--');
%test = B_z(329);
%line([test,test], [0, 2.5], 'Color', 'black', 'LineStyle', '--');
%test = B_z(377);
%line([test,test], [0, 2.5], 'Color', 'black', 'LineStyle', '--');
test = B_z(max(index));
%line([test,test], [0, 3.5], 'Color', 'black', 'LineStyle', '--');
%xlim([3.175,3.21]);
ax2 = subplot(2,1,2);
hold on;
plot(B_z, den_1D_l, 'g--');
plot(B_z, den_1D_h, 'g--');
plot(B_z, UH_1D_l, 'b--');
plot(B_z, UH_1D_h, 'b--');
for ii = 1:ray_num
    rayindex = max(find(ray_x(:,ii) ~= 0));
    surface(ax2, [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'],  [-ray_x(1:rayindex,ii)'; -ray_x(1:rayindex,ii)'], ...
    [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'], [ray_pwr(1:rayindex,ii)'./1e3; ray_pwr(1:rayindex,ii)'./1e3], ...
    'facecol', 'no', 'edgecol', 'interp', 'linew', 3);
end
caxis(ax2, [0,1e2]); colormap(flipud(hot));
ylabel(colorbar, 'Power (%) left in Ray');
ylabel(ax2, 'x (m)'); xlabel(ax2, 'z (m)'); 
ylim([-.11, .05]);
xlim([center_distance - .4, center_distance + .4]);
%pbaspect([1 1 1])
%ylim([-.03, -.01]); xlim([3.183,3.203]);
test = B_z(min(index));
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(144);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(289);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(329);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(377);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
test = B_z(max(index));
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');

%ax3 = subplot(3,1,3);
%imagesc(B_z, B_x, ne); ylabel(colorbar, 'Density (m^-^3)');
%ylabel('x (m)'); xlabel('z (m)'); ylim([-.11, .11]);
%xlim([center_distance - .4, center_distance + .4]);
%test = B_z(min(index));
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(144);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(289);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(329);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(377);
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
%test = B_z(max(index));
%line([test,test], [-.12, .12], 'Color', 'black', 'LineStyle', '--');
figure;
subplot(2,1,1);
plot(den_x,den, 'k', 'LineWidth', 1.5);
%O_x = interp1(den, den_x, .98);
%UH_B = interp2(B_z, B_x, B, 3.17, O_x);
%UH_den = (28e9^2 - 28e9^2*UH_B^2) / 8.98^2;
%UH_x = interp1(den, den_x, UH_den/ 1e19);
hold on;
plot(den_x, temp, 'r', 'LineWidth', 1.5);
ylim([0, 1.2 * max([den;temp])]); xlim([0, .04]);
%plot([UH_x, UH_x], [0,1.2 * max([den;temp])], '--b');
%plot([O_x, O_x], [0, 1.2 * max([den;temp])], '--g');
%plot(power_x,power);

lgnd = legend('Density (10^1^9 m^-^3)', 'Temperature (eV)', 'Location', 'NorthWest');
set(lgnd,'color','none');
subplot(2,1,2);
%index = 1+ find(power_col+power_e > .005 * max(power_col+power_e));
%power_x = power_x - double(power_x(max(index))) + UH_x;
plot(power_x,power_e, 'k');
hold on
plot(power_x,power_col, 'r');
%line([0,.12],[freq^2/81./1e19,freq^2/81./1e19],'Color', 'black', 'LineStyle', '--' );
xlim([0, .04]); xlabel ('r (m)');
%ylim([0, 1.2 * max(power_e + power_col)]); ylabel('Power absorbed (arb)');
%plot([UH_x, UH_x], [0, 1.2 * max(power_e + power_col)], '--b');
%plot([O_x, O_x], [0, 1.2 * max(power_e + power_col)], '--g');
lgnd = legend('Power absorbed resonance (arb)', 'Power absorbed collisions (arb)', 'Location', 'NorthWest');
set(lgnd,'color','none');

figure; hold on;
for ii = 1:ray_num
    rayindex = max(find(ray_x(:,ii) ~= 0));
    surface([ray_x(1:rayindex,ii)'; ray_x(1:rayindex,ii)'],  [ray_y(1:rayindex,ii)'; ray_y(1:rayindex,ii)'], ...
    [ray_z(1:rayindex,ii)'; ray_z(1:rayindex,ii)'], [ray_pwr(1:rayindex,ii)'./1e3; ray_pwr(1:rayindex,ii)'./1e3], ...
    'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
end
th = 0:pi/50:2*pi;
rad = .03;
plot(rad*cos(th), rad*sin(th), 'k');
rad = mean(den_1D_h);
plot(rad*cos(th), rad*sin(th), 'g--');
rad = mean(UH_1D_h);
plot(rad*cos(th), rad*sin(th), 'b--');
xlim([-.03,.03]); ylim([-.03, .03]); set(gca, 'FontSize', 12); 
xlabel('X (m)'); ylabel('Y (m)');
xticks([-.03, -.02, -.01, 0, .01, .02, .03]); yticks([-.03, -.02, -.01, 0, .01, .02, .03]); 
title('Power left in ray'); colormap(jet);
%% other  plots
figure;
x = linspace(0,40*pi,400);
y = sin(x).*cos(1.2*x);
h = plot(x,y); % capture the line handle when you plot it
cd = colormap('parula'); % take your pick (doc colormap)
cd = interp1(linspace(min(y),max(y),length(cd)),cd,y); % map color to y values
cd = uint8(cd'*255); % need a 4xN uint8 array
cd(4,:) = 255; % last column is transparency
drawnow
set(h.Edge, 'ColorBinding','interpolated', 'ColorData',cd)



%figure('Name', 'Magnetic Field')
%imagesc(B_z, B_x, B); ylabel(colorbar, 'Magnetic Field (T)');
%ylabel('x (m)'); xlabel('z (m)'); ylim([-.08, .08]);

%for ii = 1:ray_num
%    surface(ax3, [ray_z(:,ii)'; ray_z(:,ii)'],  [ray_x(:,ii)'; ray_x(:,ii)'], ...
%    [ray_z(:,ii)'; ray_z(:,ii)'], [ray_ne(:,ii)'; ray_ne(:,ii)'], ...
 %   'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
%end
%caxis(ax3, [0,6e19]); colormap(flipud(hot));
%title('Density at ray location');
%ylabel(ax3, 'x (m)'); xlabel(ax3, 'z (m)'); ylim([-.11, .11]);
%xlim([center_distance - .5, center_distance + .5]);
%% westerhof 2008 optical depth for 2nd harmonic x-mode
n = 2.;
den = [5e19, 5e19, 5e19, 5e19, 5e19, 5e19, 5e19, 5e19, 5e19];
temp=[1., 2., 5., 10., 20., 50., 100., 1000., 10000.];
ratio = (9 * den.^(1/2) / 28e9 / 1.89).^2;
ratio1 = sqrt(temp *1.6e-19 / 9.1e-31)/2.99e8;
A_n = .95 * 2^(2*n-3) * (1 + ratio / n / (n^2-1-ratio));
tau = pi^2 * n^(2*n-2) / 2^(n-1) / factorial(n-1) * A_n * ratio .* ...
    ratio1.^(2*n-2) * .042 / .003; 
single_pass = 1-exp(-tau);

x=[1., 2., 5., 10., 20., 50., 100., 1000., 10000.];
y1 = [7.01e-5,.00014, .00035, .00070, .00139, .0035, .0069, .0619, .2434]; %X-mode 105, 1e18
y2 = [6.5665e-4,.0013, .0033, .0065, .0130, .0320, .0630, .4397, .9204]; % X-mode 105, 1e19
y3 = [.0026, .0046, .0112, .0223, .0440, .1063, .2013, .8984, .9999]; % X-mode, 5-19
y3_coll = [.403,.183, .0643, .0425, .0516, .1087, .2029,.9018, .9999 ];
y4 = [1.5859e-4, 3.1709e-4, 7.9249e-4, .0016, .0032, .0079, .0157, .1022, .3378]; % O-mode 53 GHz, 1e19
figure
loglog(x,single_pass, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hold on
loglog(x,y3, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
%loglog(x,y3_coll, 'bv', 'MarkerSize', 8);
clear x y1 y2 y3
xlabel('T_e (eV)'); ylabel('Single pass absorption'); 
xlim([.5, 15000]); ylim([.001, 2]); set(gca,'FontSize', 24);
xticks([1, 100, 10000]); yticks([.01, 1]);
legend('GENRAY', 'Theory', 'Location', 'Southeast');
%legend('GENRAY', 'Theory', 'GENRAY + collisions', 'Location', 'Southeast');


%% Westerhof 2010
% O2 optical length:

vT = @(Te) sqrt(e_c*Te/m_e);
tau_On = @(n,N,Te,R,lambda,ne,B) ( ((pi^2)*(n^(2*(n-1))))/((2^(n-1))*factorial(n-1)) )*(N^(2*n - 1))*(w_pe(ne)/w_ce(B))*(R/lambda)*(vT(Te)/c_light).^(2*n)


