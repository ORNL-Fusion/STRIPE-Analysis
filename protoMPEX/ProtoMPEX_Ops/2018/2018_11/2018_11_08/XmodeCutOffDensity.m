% Calculate the x-mode upper cutoff frequency
% Based on:
% [1]R. B. Morales, S. Hacquin, S. Heuraux, and R. Sabot, “Density pro?le reconstruction methods for X-mode re?ectometry,” p. 10.

clear all
close all

% The cutoff frequency is given by:
% Right handed X-mode cutoff
% wR = (  wce + sqrt(4*wpe^2 + wce^2) )/2
% Left handed X-mode cuttoff
% wL = ( -wce + sqrt(4*wpe^2 + wce^2) )/2

% We isolate the plasma frenquency to determine what is the cutoff density.
% This cutoff density will in general depend on the radiation frequency and
% magnetic field
% This can be expressed as:

% wpe^2 = wR^2 - wR*wce
% wpe^2 = wL^2 + wL*wce

% Since, wpe^2 = (ne*e_c^2)/(m_e*e_0)

% neR = (wR^2 - wR*wce)*(m_e*e_0)/(e_c^2)
% neL = (wL^2 + wL*wce)*(m_e*e_0)/(e_c^2)

neR    = @(a,b) (a^2 - a*w_ce(b))*(m_e*e_0)/(e_c^2);
neL    = @(a,b) (a^2 + a*w_ce(b))*(m_e*e_0)/(e_c^2);
ne_UHR = @(a,b) (a^2 - w_ce(b).^2)*(m_e*e_0)/(e_c^2);
ne_OC = @(a) (a^2)*(m_e*e_0)/(e_c^2);
ne_OverDense = @(b) (b.^2)*e_0/m_e;

B = linspace(0,4);
w = 2*pi*105*1e9;

figure; hold on
h(1) = plot(neR(w,B),B,'k--');
h(2) = plot(neL(w,B),B,'r--');
h(3) = plot(ne_UHR(w,B),B,'g-');
h(4) = plot(ne_OverDense(B),B,'m:');
h(5) = line(ne_OC(w)*[1,1],[0,4]);
h(6) = line([0,2e20],(0.5*m_e*w/e_c)*[1,1]);

xlim([0,2e20])
ylim([0,4])
legend([h(1),h(2),h(3),h(4)],'X-mode R-Cutoff n_e','X-mode L-Cutoff n_e','UH Resonance n_e','wpe=wce')
box on
set(gcf,'color','w')
ylabel('B_0 [T]')
xlabel('n_e [m^{-3}]')