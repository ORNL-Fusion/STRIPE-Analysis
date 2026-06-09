% Exploring the gas-dynamic effect in Proto-MPEX
clear all 
close all

ne_a = 7.47e19; 
ne_b_xp = [7.47, 7.73, 8.07, 8.52, 8.72, 9.27]*1e19; 
Te_b_xp = [2.64, 2.69, 2.88, 2.83, 2.87, 2.94]; % eV

B0 = 0.055; 
B2 = 0.6;

B1type = 1;
switch B1type
    case 1
        B1 = [0.5, 0.72, 0.94, 1.17, 1.39, 1.61];
end

R2a = B2/B0;
R1a = B1(1)/B0;

R2b = B2/B0;
R1b = B1/B0;

Ra = (R1a*R2a)./(R1a + R2a);
Rb = (R1b*R2b)./(R1b + R2b);

% Ra = R2a.*(R1a+R2a)./(R1a+(2*R2a));
% Rb = R2b.*(R1b+R2b)./(R1b+(2*R2b));

Csa = C_s(Te_b_xp(1),2);
Csb = C_s(Te_b_xp,2);

ne_b_Te = ne_a*(Rb./Ra).*(Csa./Csb);
ne_b = ne_a*(Rb./Ra);

figure;
hold on
h(1) = plot(R1b,ne_b_xp,'ko-');
h(2) = plot(R1b,ne_b_Te,'ro-');
h(3) = plot(R1b,ne_b,'gsq-');

xlim([5,35])
ylim([0,1.2e20])
xlabel('R, mirror ratio')
ylabel('n_e [m^{-3}]')
box on
grid on

