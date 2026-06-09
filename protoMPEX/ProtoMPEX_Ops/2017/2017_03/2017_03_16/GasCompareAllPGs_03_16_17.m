% Check gas 

close all
clear all

mdsconnect('mpexserver')

Shots = 13303; % ONLY ONE SHOT

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % 9.5, use x2
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % 2.5, use x2
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % 6.5, use x2
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % 4.5, use x10

[PG,t] = my_mdsvalue_v2(Shots,address);
[PS1,t_ps1] = my_mdsvalue_v2(Shots,{'\MPEX::TOP.MACHOPS1:PS1_I'});

MF = [2,2,2,10]; % convert to mTorr;

%%
figure; hold on
for s = 1 :4
h(s) = plot(t{s},PG{s}*MF(s));
end
set(h(3),'LineWidth',2)
plot(t_ps1{1},PS1{1},'k:','LineWidth',1)
xlim([3,6])
legend(h,'S9.5','S2.5','S6.5','S4.5')
ylabel('[mTorr]')
xlabel('t [sec]')
box on
set(gcf,'color','w')