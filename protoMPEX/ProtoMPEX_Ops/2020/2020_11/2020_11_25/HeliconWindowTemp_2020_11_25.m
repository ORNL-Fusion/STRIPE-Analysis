% Monitor the temperature of the helicon window
clear all
close all

shotlist = [30000 + [925,927,928,929,930,931]]; 
shotlist = [30000 + [944]]; 

deltaT = 0;

% To convert voltage to temperature, multiply by 20
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

% Gather data
DA{1} = [RootAddress,'FLUOROPT_1']; 
DA{2} = [RootAddress,'FLUOROPT_2']; 
DA{3} = [RootAddress,'FLUOROPT_3']; 
DA{4} = [RootAddress,'FLUOROPT_4']; 

[f1,t_f1]   = my_mdsvalue_v2(shotlist,DA(1));
[f2,t_f2]   = my_mdsvalue_v2(shotlist,DA(2));
[f3,t_f3]   = my_mdsvalue_v2(shotlist,DA(3));
[f4,t_f4]   = my_mdsvalue_v2(shotlist,DA(4));

C1 = {'b','b.','b:','b--','b','b.','b:','b--'}
C2 = {'r','r.','r:','r--','b','b.','b:','b--'}
C3 = {'g','g.','g:','g--','b','b.','b:','b--'}
C4 = {'c','c.','c:','c--','b','b.','b:','b--'}

figure
for s = 1:length(shotlist)
    subplot(5,1,s)
    hold on
    h(1) = plot(t_f1{s}(1:end-1),(f1{s}-deltaT*min(f1{s}))*20,C1{s},'LineWidth',1)
    h(2) = plot(t_f2{s}(1:end-1),(f2{s}-deltaT*min(f2{s}))*20,C2{s},'LineWidth',1) % Ground Side
    h(3) = plot(t_f3{s}(1:end-1),(f3{s}-deltaT*min(f3{s}))*20,C3{s},'LineWidth',1) % High Voltage side
    h(4) = plot(t_f4{s}(1:end-1),(f4{s}-deltaT*min(f4{s}))*20,C4{s},'LineWidth',1)
    legend(h,'1','2','3','4')
%     ylim([0,30])
    ylabel('${\Delta}T$ $[C]$','Interpreter','latex','FontSize',14)
    xlabel('$t$ $[sec]$','Interpreter','latex','FontSize',14)
    box on
    set(gcf,'color','w')
    grid on
    title(num2str(shotlist(s)))
end