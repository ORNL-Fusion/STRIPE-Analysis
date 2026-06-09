% close all
% clear all

shotlist = 13358; % March 21st 2017, 4.8 cm diam insulating skimmer
shotlist = 12867; % Jan 12th 2017  , 5.8 cm diam conducting skimmer

mdsconnect('mpexserver')

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % 9.6
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % 2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % 6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % 4.5

[a,t9] = my_mdsvalue_v2(shotlist,address(1));
% [b,~] = my_mdsvalue_v2(13356,address(1));
PG9 = (a{1}-mean(a{1}(1:100)))*2;
[a,t2] = my_mdsvalue_v2(shotlist,address(2));
% [b,~] = my_mdsvalue_v2(13356,address(2));
PG2 = (a{1}-mean(a{1}(1:100)))*2;
[a,t6] = my_mdsvalue_v2(shotlist,address(3));
% [b,~] = my_mdsvalue_v2(13356,address(3));
PG6 = (a{1}-mean(a{1}(1:100)))*2;
[a,t4] = my_mdsvalue_v2(shotlist,address(4));
% [b,~] = my_mdsvalue_v2(13356,address(4));
PG4 = (a{1}-mean(a{1}(1:100)))*10;


%%
close all
figure; hold on
h(1) = plot(t2{1}(1:length(PG2)),PG2,'k')
h(2) = plot(t4{1}(1:length(PG4)),PG4,'bl')
h(3) = plot(t6{1}(1:length(PG6)),PG6,'r')
h(4) = plot(t9{1}(1:length(PG9)),PG9,'g')
set(h,'LineWidth',1)
legend(h,'2.5','4.5','6.5','9.5')
xlim([4,5])

legend(h,{num2str(shotlist')},'location','NorthEast')
ylabel('mTorr')
ylim([-0.5,22])
xlim([4,5])
grid on
set(gcf,'color','w')
set(gcf,'Position',[360 210  340  400])
box on
