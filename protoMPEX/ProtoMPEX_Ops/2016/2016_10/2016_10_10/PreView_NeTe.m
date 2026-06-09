% Preview data for Oct 20th 2016

clear all
close all

option = 2; 

switch option
    case 1
    shotlist = 1e4 + 1000 + [6  , 17 , 20, 23 ,  29];
    pwr      =              [50 , 66 , 75, 90 , 107]; % [kW]
    case 2
    shotlist = 1e4 + 1000 + [6  , 17 , 20, 23, 24 , 25,  29];
    pwr      =              [50 , 66 , 75, 90, 90 , 90, 107]; % [kW]
end

% shotlist = 11000 + [95,96];

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{2} = [RootAddress,'TARGET_LP']; % I
DataAddress{3} = [RootAddress,'FLUOROPT_2']; 
DataAddress{4} = [RootAddress,'FLUOROPT_3']; 
DataAddress{5} = [RootAddress,'RF_FWD_PWR']; 

%DataAddress{3} = [RootAddress,'PWR_28GHz'];
% Data = [RootAddress,'PG2'];

[V,t_V] = my_mdsvalue_v3(shotlist,DataAddress{1});
[I,t_I] = my_mdsvalue_v3(shotlist,DataAddress{2});
[f2,t_f2] = my_mdsvalue_v3(shotlist,DataAddress{3});
[f3,t_f3] = my_mdsvalue_v3(shotlist,DataAddress{4});
[Fwd,t_Fwd] = my_mdsvalue_v3(shotlist,DataAddress{5});


%%
close all

if 0
figure; 
for s = 1:length(shotlist)
subplot(3,3,s); hold on
plot(t_V{s},V{s})
plot(t_I{s},I{s});
ylim([-10,10])
end
end

figure; 
for s = 1:length(shotlist)
    subplot(3,3,s); hold on;
   
    % Probe 2
    h(1) = plot(t_f2{s},f2{s}*20,'k');
    % Probe 3
    h(2) = plot(t_f3{s},f3{s}*20,'r');
    
    % Change in Temperature
    dt1(s) = max(f2{s}*20)-min(f2{s}*20);
    dt2(s) = max(f3{s}*20)-min(f3{s}*20);
    
    % RF power pulse
    plot(t_Fwd{s},20+10*Fwd{s})
    
    ylabel('[deg C]','Fontsize',8)
    xlabel('[sec]','Fontsize',8)
    title([num2str(pwr(s)),' kW'],'Fontsize',8)
    box on
    L = legend([h],num2str(dt1(s)),num2str(dt2(s)));
    set(L,'Fontsize',6)
    ylim([20,80]);
    xlim([0,50])
end
set(gcf,'color','w')

figure; hold on
plot(pwr,dt1,'ko')
plot(pwr,dt2,'ro')

if option == 1
    % save data
    % Probe 2 is at the ground side
    % Probe 3 is at the high voltage side
    
    % Probe 2
    D = [t_f2{1},f2{1}*20,f2{2}*20,f2{3}*20,f2{4}*20,f2{5}*20];
    F = {'t[sec]',[num2str(pwr(1)),' kW'],[num2str(pwr(2)),' kW'],[num2str(pwr(3)),' kW'],...
        [num2str(pwr(4)),' kW'],[num2str(pwr(5)),' kW']};
    FileName = 'FluoropticProbe_2_Temp_vs_time_PowerScan.xlsx';
    xlswrite(FileName,[F;num2cell(D)]);
    
    % Probe 3 
    D = [t_f3{1},f3{1}*20,f3{2}*20,f3{3}*20,f3{4}*20,f3{5}*20];
    F = {'t[sec]',[num2str(pwr(1)),' kW'],[num2str(pwr(2)),' kW'],[num2str(pwr(3)),' kW'],...
        [num2str(pwr(4)),' kW'],[num2str(pwr(5)),' kW']};
    FileName = 'FluoropticProbe_3_Temp_vs_time_PowerScan.xlsx';
    xlswrite(FileName,[F;num2cell(D)]);
end
    
% figure; 
% plot(t_Fwd{1},Fwd{1})
% ylim([0,10])