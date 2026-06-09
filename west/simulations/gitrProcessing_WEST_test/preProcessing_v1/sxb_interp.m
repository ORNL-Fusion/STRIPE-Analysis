close all;
clear all;
clc;
load('w_4009_finer_grid.mat');
ne1=readmatrix('ne.csv');
te1=readmatrix('te.csv');

[NX, NY]=size(ne1);
ne=1E6*ne;
Te=te;
figure;
sxb_new=interpn(te,ne,wi_4009_sxb,te1,ne1);
% imagesc(sxb_new)
% set(gca,'YDir','normal')


% for i=1:NX
% for j=1:NY
% sxb_new=interp2(te,ne,wi_4009_sxb',te1,ne1);
% end
% end
% writematrix(sxb_new,'sxb_interp.csv');
% z0=readmatrix('z.csv');
% x0=readmatrix('r.csv');

% [m,n]=size(sxb_new);
% % x=1:m;y=1:n;
% % psi_new=(interp2(psi,0.0575/(max(x0)-min(x0))*m,1.745/(max(z0)-min(z0))*n));
% psi_new=(interpn(sxb_new,0.1/(max(z0)-min(z0))*n,2.86/(max(x0)-min(x0))*m));


figure; imagesc(r,z,sxb_new)
set(gca,'YDir','normal','FontSize', 18);
xlabel('r [m]');
ylabel('z [m]')
title('S/XB factor for W I 400.9nm emission')

figure; plot(te,wi_4009_sxb(:,1))
xlim([0 20])
ylim([1E-1 21E2]);
ylim([1E-1 2E2]);
set(gca, 'YScale', 'log', 'FontSize', 18);
xlabel('T_e [eV]');
ylabel('S/XB, Inverse photon efficiency')
title('W I (400.9nm) emission')

figure; plot(ne,wi_4009_sxb(1,:))
xlim([0 20])
ylim([1E-1 21E2]);
ylim([1E-1 2E2]);
set(gca, 'YScale', 'log', 'FontSize', 18);
xlabel('T_e [eV]');
ylabel('S/XB, Inverse photon efficiency')
title('W I (400.9nm) emission')
return;

%% Noise Scan

fontSize.axes   = 18;
fontSize.labels = 20;

z=[-0.4; -0.3; -0.2; -0.1; 0; 0.1; 0.2; 0.3; 0.4];
indexL1= [247;243;239;235;225;227;223;210;215];
indexR1= [120;116;106;108;102;98;92;90;85];
erosionFluxL1= [1.0844E16;1.6790E16;1.8553E16;2.7383E16;1.5031E16;1.9433E16;1.9437E16;2.1208E16;1.7168E16];
erosionFluxR1= [1.7042;2.7799;2.8317;2.3794;4.8990;3.0200;3.2742;2.1762;1.4202].*1E16;

erosionFluxR20k= [1.6209;1.7106;2.8040;2.9334;3.7851;3.1951;3.3011;1.9262;1.7107].*1E18;
erosionFluxR50k= [1.5668;1.6658;2.9690;2.9921;3.5459;3.3522;3.1068;1.8700;1.6830].*1E18;
erosionFluxR100k= [1.5066;1.7044;2.9067;3.0347;3.5384;3.3164;3.0682;1.9309;1.7082].*1E18;





sxbValueL=[9.58779;10.8139;12.1153;12.8246;13.4241;13.242;12.4221;11.1808;9.70449];
sxbValueR=[9.58544;10.2572;11.8979;12.9097;13.436;13.1903;12.3775;10.6288;9.92874];

photonFluxL=erosionFluxL1./sxbValueL/(4*pi);
photonFluxR=erosionFluxR1./sxbValueR/(4*pi);

photonFluxR20k=erosionFluxR20k./(sxbValueR.*(4*pi));
photonFluxR50k=erosionFluxR50k./(sxbValueR.*(4*pi));
photonFluxR100k=erosionFluxR100k./(sxbValueR.*(4*pi));

figure; plot(z,photonFluxR20k, 'bo--', 'LineWidth',2);hold on;
plot(z,photonFluxR50k, 'ro--', 'LineWidth',2);
plot(z,photonFluxR100k, 'ko--', 'LineWidth',2);hold on;
% figure; plot(z,erosionFluxR, 'bo--', 'LineWidth',2);
xlabel('z [m]','Interpreter','latex','fontSize',fontSize.labels)
ylabel('$\Gamma_\phi [/m^{2}sSt]$ ','Interpreter','latex','fontSize',fontSize.labels)
set(gca,'fontName','Times','fontSize',fontSize.axes);



%% Potential Scan

fontSize.axes   = 18;
fontSize.labels = 20;

z=[-0.4; -0.3; -0.2; -0.1; 0; 0.1; 0.2; 0.3; 0.4];
indexL1= [247;243;239;235;225;227;223;210;215];
indexR1= [120;116;106;108;102;98;92;90;85];
erosionFluxL1= [1.0844E16;1.6790E16;2.4553E16;2.7383E16;2.6031E16;2.833E16;3.8437E16;2.31208E16;1.8168E16];
erosionFluxR1= [1.7042;2.7799;2.8317;2.3794;4.8990;3.0200;3.2742;2.1762;1.4202].*1E16;

erosionFluxR100k_p1te= [1.4773;1.6232;2.5565;3.1908;4.1388;3.6225;3.2729;1.9029;1.8046].*1E18;
erosionFluxR100k_1te= [1.5066;1.7044;2.9067;4.0347;5.1384;4.3164;4.0682;1.9309;1.7082].*1E18;
erosionFluxR100k_50te= [1.4130;2.1624;3.846;4.7408;5.9162;5.3;4.8;2.1338;1.8338].*1E18;





sxbValueL=[9.58779;10.8139;12.1153;12.8246;13.4241;13.242;12.4221;11.1808;9.70449];
sxbValueR=[9.58544;10.2572;11.8979;12.9097;13.436;13.1903;12.3775;10.6288;9.92874];



photonFluxR100k_p1te=erosionFluxR100k_p1te./sxbValueR/(4*pi);
photonFluxR100k_1te=erosionFluxR100k_1te./sxbValueR/(4*pi);
photonFluxR100k_50te=erosionFluxR100k_50te./sxbValueR/(4*pi);


figure; plot(z,0.04.*photonFluxR100k_p1te, 'ro--', 'LineWidth',2);hold on;
plot(z,0.04.*photonFluxR100k_1te, 'bo--', 'LineWidth',2);
plot(z,0.04.*photonFluxR100k_50te, 'go--', 'LineWidth',2);hold on;
% figure; plot(z,erosionFluxR, 'bo--', 'LineWidth',2);
xlabel('z [m]','Interpreter','latex','fontSize',fontSize.labels)
ylabel('$\Gamma_\phi [/m^{2}sSt]$ ','Interpreter','latex','fontSize',fontSize.labels)
set(gca,'fontName','Times','fontSize',fontSize.axes);
hold on;

plot(z,0.35.*photonFluxR100k_p1te, 'b.-', 'LineWidth',2);hold on;
plot(z,0.035.*photonFluxR100k_1te, 'r.-', 'LineWidth',2);
plot(z,0.035.*photonFluxR100k_50te, 'k.-', 'LineWidth',2);hold on;
% figure; plot(z,erosionFluxR, 'bo--', 'LineWidth',2);
xlabel('z [m]','Interpreter','latex','fontSize',fontSize.labels)
ylabel('$\Gamma_\phi [/m^{2}sSt]$ ','Interpreter','latex','fontSize',fontSize.labels)
set(gca,'fontName','Times','fontSize',fontSize.axes);

hold on;

plot(z,0.005.*photonFluxR100k_p1te, 'b.-', 'LineWidth',2);hold on;
plot(z,0.005.*photonFluxR100k_1te, 'r.-', 'LineWidth',2);
plot(z,0.005.*photonFluxR100k_50te, 'k.-', 'LineWidth',2);hold on;
% figure; plot(z,erosionFluxR, 'bo--', 'LineWidth',2);
xlabel('z [m]','Interpreter','latex','fontSize',fontSize.labels)
ylabel('$\Gamma_\phi [/m^{2}sSt]$ ','Interpreter','latex','fontSize',fontSize.labels)
set(gca,'fontName','Times','fontSize',fontSize.axes);


%% Ohmic and LH cases

fontSize.axes   = 18;
fontSize.labels = 20;

z=[-0.4; -0.3; -0.2; -0.1; 0; 0.1; 0.2; 0.3; 0.4];


erosionFluxR_Ohmic= [1.7351;1.9381;1.9303;1.9228;1.9683;1.8896;1.9168;1.8343;1.8809].*1E18;
erosionFluxR_LH1p5= [5.3378;5.4220;5.5793;5.6826;5.9116;5.8133;6.2118;5.5443;5.8520].*1E18;
erosionFluxR_LH2p0= [9.2946;9.8295;9.9761;10.101;10.05;10.039;10.948;9.5554;10.597].*1E18;






sxbValueR_Ohmic=[9.58544;10.2572;11.8979;12.9097;13.436;13.1903;12.3775;10.6288;9.92874];



photonFluxR_Ohmic=erosionFluxR_Ohmic./sxbValueR/(4*pi);
photonFluxR_LH1p5=erosionFluxR_LH1p5./sxbValueR/(4*pi);
photonFluxR_LH2p0=erosionFluxR_LH2p0./sxbValueR/(4*pi);


figure; plot(z,photonFluxR_Ohmic, 'bo--', 'LineWidth',2);hold on;
plot(z,photonFluxR_LH1p5, 'ro--', 'LineWidth',2);
plot(z,photonFluxR_LH2p0, 'ko--', 'LineWidth',2);hold on;
% figure; plot(z,erosionFluxR, 'bo--', 'LineWidth',2);
xlabel('z [m]','Interpreter','latex','fontSize',fontSize.labels)
ylabel('$\Gamma_\phi [/m^{2}sSt]$ ','Interpreter','latex','fontSize',fontSize.labels)
set(gca,'fontName','Times','fontSize',fontSize.axes);



