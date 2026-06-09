% Compare the data from DLP boxes
close all 
clear all

Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % Vs
DataAddress{2} = [RootAddress,'TARGET_LP']; % Is
DataAddress{3} = [RootAddress,'INT_4MM_1']; % Vsx
DataAddress{4} = [RootAddress,'INT_4MM_2']; % Isx

shotlist = 1e4 + 2e3 + [125];
shotlist = 12600 + [43];

[Vs,t_vs]   = my_mdsvalue_v2(shotlist,DataAddress(1))
[Vsx,t_vsx] = my_mdsvalue_v2(shotlist,DataAddress(3))

[Is,t_Is]   = my_mdsvalue_v2(shotlist,DataAddress(2))
[Isx,t_Isx] = my_mdsvalue_v2(shotlist,DataAddress(4))


%%
figure; hold on
V_Att = 2; 
V_cal = 12.05;
plot(t_vs{1}, Vs{1}*V_Att*V_cal,'k.')
plot(t_vsx{1},Vsx{1}/(0.46e-3),'r-')
ylim([-100,100])
xlim([4.12,4.32])
ylabel('voltage [V]')
set(gcf,'color','w'); box on


figure; hold on
V_Att = 5; 
V_cal = 12.05;
plot(t_vs{1}, 1000*Is{1}*V_Att/142.5,'k')
plot(t_vsx{1},-1000*Isx{1},'r')
ylim([-100,100])
xlim([4.12,4.32])
ylabel('current [A]')
set(gcf,'color','w'); box on

figure
subplot(2,1,1)
plot(t_vs{1}, Is{1},'k')
title('I_raw')
ylim([-10,10])
subplot(2,1,2)
plot(t_vs{1}, Vs{1},'k')
title('V_raw')
ylim([-10,10])

figure; 
rng = find(t_vsx{1} >=4.28 & t_vsx{1} <=4.285);
h(1) = plot(Vsx{1}(rng)/(0.46e-3),-1000*Isx{1}(rng),'k.-')
ylim([-120,120])
xlim([-120,120])
ylabel('Current [mA]')
xlabel('voltage [V]')
set(gcf,'color','w'); box on
legend(h,'shot: 11345, time: 4.28 to 4.285')
