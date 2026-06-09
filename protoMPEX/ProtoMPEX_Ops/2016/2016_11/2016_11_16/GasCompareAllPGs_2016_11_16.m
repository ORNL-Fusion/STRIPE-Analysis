% Check gas 

close all
clear all

mdsconnect('mpexserver')

Shots = 11576; % ONLY ONE SHOT

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % 9.5, use x2
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % 2.5, use x2
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % 6.5, use x2
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % 4.5, use x10

[PG,t] = my_mdsvalue_v2(Shots,address);
[PS1,t_ps1] = my_mdsvalue_v2(Shots,{'\MPEX::TOP.MACHOPS1:PS1_I'});

MF = [2,2,2,10]; % convert to mTorr;

%%
figure; hold on
C = {'k','r','bl','g'};
for s = 1 :4
h(s) = plot(t{s},(PG{s}-mean(PG{s}(1:100)) )*MF(s),C{s},'LineWidth',2);
end
% Plot RF pulse region
line([4.24,4.24],[0,50],'LineWidth',0.5)
% Plot RF pulse region
line([4.30,4.30],[0,50],'LineWidth',0.5)

% set(h(3),'LineWidth',2)
xlim([4,5]); ylim([0,20])
legend(h,'S9.5','S2.5','S6.5','S4.5')
ylabel('[mTorr]')
xlabel('t [sec]')
box on
set(gcf,'color','w')