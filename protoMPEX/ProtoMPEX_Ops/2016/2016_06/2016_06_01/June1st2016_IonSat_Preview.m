close all
clear all

mdsconnect('mpexserver')

Spool = 6.5;
Shots = [8953:8963];
R = [6:-1:-4];

switch Spool 
    case 6.5
        L_tip = 1.3/1000;  % [m]
        D_tip = 0.25/1000; % [m]
    case 9.5
        L_tip = 1.2/1000;  % [m]
        D_tip = 0.25/1000; % [m]
    case 4.5
        L_tip = 1.2/1000;  % [m]
        D_tip = 0.508/1000; % [m] % December 5th 2016, 
end
rp = D_tip/2; % [m]

AreaType = 1;
switch AreaType
    case 1
        Area = L_tip*pi*D_tip + 0.25*pi*D_tip^2;
    case 2
        Area = 2*L_tip*D_tip;
end
AMU = 2; % D2 gas

Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RA = [Stem,Branch];
Data{1} = [RA,'HELICON_LP'];
Data{2} = [Stem,'SHOT_NOTE'];

[Isat,t_Isat] = my_mdsvalue_v2(Shots,Data(1));
[NOTE,~] = my_mdsvalue_v2(Shots,Data(2));

% Filter data:
for s = 1:length(Shots)
    Isat_f{s} = sgolay_t(-2*Isat{s}/50.2,3,31);
end

% High and low density regions:
t_HD_Gather_Start = 4.3;
t_HD_Gather_End   = 4.31;

t_LD_Gather_Start = 4.16;
t_LD_Gather_End   = 4.19;

for s = 1:length(Shots)

    rng2{s} = find(t_Isat{s}>=t_HD_Gather_Start & t_Isat{s}<=t_HD_Gather_End);
    Isat_HD(s) = mean(Isat_f{s}(rng2{s}));
    
    rng3{s} = find(t_Isat{s}>=t_LD_Gather_Start & t_Isat{s}<=t_LD_Gather_End);
    Isat_LD(s) = mean(Isat_f{s}(rng3{s}));
end

%% Plot data: 
close all

figure; hold on
Te = 5; 
for s = 1:length(Shots)
    plot(t_Isat{s},Isat_f{s}./(0.61*e_c*Area*C_s(Te,2)));
end
% ylim([0,0.1]);
ylim([0,12e19]);
xlim([4.13,4.35])

if 1
    figure; hold on
    plot(R,Isat_HD./(0.61*e_c*Area*C_s(Te,2)),'ko-');
end
ylim([0,6e19]);
xlim([-6,6])
set(gcf,'color','w');
set(gca,'Box','on')
xlabel('Vertical radial location [cm]')
ylabel('n_i based on T_e = 8 eV,  [m^{-3}]')
