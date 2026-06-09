% Sept 29th, combined Helicon + ECH + ICH into a tungsten target XP
close all
clear all

shotlist = 1e4 + 600 + [26];
%shotlist = 1e4 + 600 + [34];

% Acquiring Ne and Te data
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

% 28 GHz delivered
Data = [RootAddress,'PWR_28GHz'];
[f,q] = my_mdsvalue_v3(shotlist,Data);
s{1} = f{1}*10; % kW
t{1} = q{1};

% 18 GHz delivered, 0.0302 V/kW : 33.1 kW/V
Data = [RootAddress,'ICH_LP'];
[f,q] = my_mdsvalue_v3(shotlist,Data);
s{2} = abs(f{1})*33.1;
t{2} = q{1};

% ICH fwd
Data = [RootAddress,'GAS_FLOW_1'];
[f,q] = my_mdsvalue_v3(shotlist,Data);
s{3} = f{1}*-1*320;
t{3} = q{1};
% 26 at 12.97

% Helicon fwd:
Data = [RootAddress,'RF_FWD_PWR'];
[f,q] = my_mdsvalue_v3(shotlist,Data);
s{4} = f{1}*210;
t{4} = q{1};
%%
close all

figure;
%subplot(1,2,1);
%title(['shot: ',num2str(shotlist)])
hold on
LineW = 2;
for k = 1:4
    h(k) = plot(t{k},s{k},'LineWidth',LineW);
end
ylim([0,120])
xlim([4.14,4.33])
set(gca,'PlotBoxAspectRatio',[1 1 1])
legend(h,'28 GHz','18 GHz','ICH','Helicon','Location','NorthWest')
ylabel('[kW]','Interpreter','Latex')
xlabel('t [sec]','Interpreter','Latex')
grid on
box on
set(gcf,'color','w')

figure;
subplot(2,1,1); hold on
 hold on
LineW = 2;
for k = 1:4
    h(k) = plot(t{k},s{k},'LineWidth',LineW);
end
ylim([0,120])
xlim([4.245,4.275])
set(gca,'PlotBoxAspectRatio',[1 0.6 1])
legend(h,'28 GHz','18 GHz','ICH','Helicon','Location','NorthWest')
ylabel('[kW]','Interpreter','Latex')
xlabel('t [sec]','Interpreter','Latex')
grid on
box on
set(gcf,'color','w')

% figure; 
subplot(2,1,2); hold on
LineW = 2;

hN = plot(t{1},s{1} + s{2} + s{3} + s{4},'LineWidth',LineW);
line([0,5],[150,150])
ylim([100,160])
xlim([4.245,4.28])
set(gca,'PlotBoxAspectRatio',[1 0.6 1],'YTick',[100,120,140,150,160])
legend(hN,'Net power','Location','SouthWest')
ylabel('Net Power [kW]','Interpreter','Latex')
xlabel('t [sec]','Interpreter','Latex')
grid on
box on
set(gcf,'color','w')
