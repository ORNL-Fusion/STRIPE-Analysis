% Monitor the temperature of the helicon window
clear all
close all

shotlist = 19800 + [61,62];
shotlist = 19800 + [62,63];
shotlist = 19800 + [63,64];
shotlist = 19800 + [64,65];
shotlist = 19800 + [64,65,66];
shotlist = 19800 + [65,66,67];
shotlist = 19800 + [67,68,69];
shotlist = 19800 + [87:90];
shotlist = 19800 + [98,99,100,102];
shotlist = 19800 + [103,104,105];
shotlist = 19800 + [107];
shotlist = 19900 + [59];

shotlist = [20255,19768];

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

C1 = {'b','b.','b:','b--'}
C2 = {'r','r.','r:','r--'}
C3 = {'g','g.','g:','g--'}
C4 = {'c','c.','c:','c--'}

figure; hold on
for s = 1:length(shotlist)
    plot(t_f1{s}(1:end-1),(f1{s}-1*min(f1{s}))*20,C1{s},'LineWidth',2)
    plot(t_f2{s}(1:end-1),(f2{s}-1*min(f2{s}))*20,C2{s},'LineWidth',2) % Ground Side
    plot(t_f3{s}(1:end-1),(f3{s}-1*min(f3{s}))*20,C3{s},'LineWidth',2) % High Voltage side
    plot(t_f4{s}(1:end-1),(f4{s}-1*min(f4{s}))*20,C4{s},'LineWidth',2)
end

ylim([0,100])
ylabel('${\Delta}T$ $[C]$','Interpreter','latex','FontSize',14)
xlabel('$t$ $[sec]$','Interpreter','latex','FontSize',14)
box on
set(gcf,'color','w')
grid on