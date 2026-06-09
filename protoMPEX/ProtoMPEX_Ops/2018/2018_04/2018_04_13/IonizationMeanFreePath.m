% In this code we calculate the ionization mean free path for D2 molecules
% emitted from a plasma-surface interface.

%% Molecular neutral deuterium

close all
clear all

m_p = 1.6726e-27;
m_e = 9.1094e-31;
k_B = 1.3806e-23;
e_c = 1.6020e-19;

Tg = 300; 
A = 2; 
v_th_D2 = sqrt(k_B*Tg/(2*m_p*A)); % D2 has two protons and two neutrons

ne = [1.5,4,8,16]*1e19; 
Te = linspace(1.5,10);
Rp = 3/100;
D = 2*Rp;
% D = 1;

% Ionization rate of D2 (Hjatarson 2010)
Kiz_D2 = @(t) (1.1e-14)*(t.^0.42).*exp(-16.05./t);
C = {'k','r','bl','g'};

figure; hold on
for n = 1:length(ne)
viz_D2{n} = ne(n).*Kiz_D2(Te);
Liz_D2{n} = v_th_D2./viz_D2{n};

h(n) = plot(Te,viz_D2{n}*1e-6,C{n},'LineWidth',2)
set(gca,'YScale','log','Xscale','lin','FontName','Times')
title('$\nu_{iz}$ [MHz]','Interpreter','Latex','FontSize',12)
xlabel('$T_e$ [eV]','Interpreter','Latex','FontSize',11)
grid on
box on

leg{n} = ['n_e = ',num2str(ne(n)'*1e-19),'x10^{19} [m^{-3}]'];
end
set(gcf,'color','w','Position',[537.0000   84.3333  383.3333  300])
L = legend(h,leg); L.Box = 'off'; L.Location = 'southeast';
% ylim([1e0,1e4])
% set(gca,'YTick'     ,[1e-2,1e-1,1,1e1,1e2,1e3])
% set(gca,'YTickLabel',[1e-2,0.1 ,1,10 ,1e2,1e3])


figure; hold on

for n = 1:length(ne)
viz_D2{n} = ne(n).*Kiz_D2(Te);
Liz_D2{n} = v_th_D2./viz_D2{n};

h(n) = plot(Te,Liz_D2{n}/D,C{n},'LineWidth',2);
set(gca,'YScale','log','Xscale','lin','FontName','Times')
title(['${\lambda_+^{D_2}}/\phi,$ $ T_0^{D_2}=$ ',num2str(Tg),' K'],'Interpreter','Latex','FontSize',12)
xlabel('$T_e$ [eV]','Interpreter','Latex','FontSize',11)
grid on
box on

leg{n} = ['n_e = ',num2str(ne(n)'*1e-19),'x10^{19} [m^{-3}]'];
end
set(gcf,'color','w','Position',[537.0000   84.3333  383.3333  300])
L = legend(h,leg); L.Box = 'off';
ylim([1e-2,1e3])
set(gca,'YTick'     ,[1e-2,1e-1,1,1e1,1e2,1e3])
set(gca,'YTickLabel',[1e-2,0.1 ,1,10 ,1e2,1e3])


%% Atomic deuterium neutral

% In order to calculate the ionization mean free path of an atomic neutral
% coming out of the target, we need to know at what energy it comes out of
% the target. This energy determines the velocity of the atom towards the
% plasma. We then need to know the ionization frequency.

% In this particular case, we are interested in a Stainless steel target
% After talking to Juergen and Josh, the energy of the outgoing neutral
% atom is given by the energy reflection coefficient and the energy of the
% incoming ion accross the sheath.
% Using Reference: "Calculated Sputtering, Reflection and Range Values" by
% Eckstein 2002, The energy reflection coefficient for an ion falling down a
% sheath potential of 7 eV and colliding with a Fe surface leads to a
% energy reflection coefficient of ~ 40% and particle reflection
% coefficient of 60%. Based on the trend in the data, it appears that the
% energy reflection coefficient increases up to 50% for Te < 7eV. In the
% range we are interested, betweem 2 and 10 eV, we can take the energy
% reflection coefficient to be ~ 40 %

R_E = [4.02,3.91,3.81,3.73,3.50,3.27,2.75]*1e-1;
kT  = [7   ,8   ,9   ,10  ,14  ,20  ,50  ];

figure; plot(kT,R_E,'ko-')
xlim([0,50])
ylim([0,1])
set(gcf,'color','w'); box on
xlabel('kT [eV]')
ylabel('R_E')
ylim([0,0.5])
xlim([0,50])

% Note from 2018_11_11, the energy if the ion should be determined based on
% the floating potential + presheath energy

E_D0 = 0.4*3*Te;
uD0 = sqrt(2*e_c*E_D0/(2*m_p));
Kiz_D0 = @(t) (7.89e-15)*(t.^0.41).*exp(-14.23./t);

figure; hold on
for n = 1:length(ne)
viz_D0{n} = ne(n).*Kiz_D0(Te);
Liz_D0{n} = uD0./viz_D0{n};

h(n) = plot(Te,Liz_D0{n}/D,C{n},'LineWidth',2);
set(gca,'YScale','log','Xscale','lin','FontName','Times')
title(['${\lambda_+^{D_0}}/\phi,$ $E_{D_0}=R_E3k_BT_e,$ $R_E = 0.4$ '],'Interpreter','Latex','FontSize',12)
xlabel('$T_e$ [eV]','Interpreter','Latex','FontSize',11)
grid on
box on

leg{n} = ['n_e = ',num2str(ne(n)'*1e-19),'x10^{19} [m^{-3}]'];
end
set(gcf,'color','w','Position',[537.0000   84.3333  383.3333  300])
L = legend(h,leg); L.Box = 'off'; L.Location = 'southwest';
ylim([1e-2,1e3])
set(gca,'YTick'     ,[1e-2,1e-1,1,1e1,1e2,1e3])
set(gca,'YTickLabel',[1e-2,0.1 ,1,10 ,1e2,1e3])

% Charge exchange D --> D+ From ADAS in the range of 0.1 to 1 eV ions
Kcx_D0 = @(t) 0.5e-14;

figure; hold on
for n = 1:length(ne)
vcx_D0{n} = ne(n).*Kcx_D0(Te);
Lcx_D0{n} = uD0./vcx_D0{n};

h(n) = plot(Te,Lcx_D0{n}/D,C{n},'LineWidth',2);
set(gca,'YScale','log','Xscale','lin','FontName','Times')
title(['${\lambda_{cx}^{D_0}}/\phi,$ $E_{D_0}=R_E3k_BT_e,$ $R_E = 0.4$ '],'Interpreter','Latex','FontSize',12)
xlabel('$T_e$ [eV]','Interpreter','Latex','FontSize',11)
grid on
box on

leg{n} = ['n_e = ',num2str(ne(n)'*1e-19),'x10^{19} [m^{-3}]'];
end
set(gcf,'color','w','Position',[537.0000   84.3333  383.3333  300])
L = legend(h,leg); L.Box = 'off'; L.Location = 'NorthWest';
ylim([1e-2,1e3])
set(gca,'YTick'     ,[1e-2,1e-1,1,1e1,1e2,1e3])
set(gca,'YTickLabel',[1e-2,0.1 ,1,10 ,1e2,1e3])



