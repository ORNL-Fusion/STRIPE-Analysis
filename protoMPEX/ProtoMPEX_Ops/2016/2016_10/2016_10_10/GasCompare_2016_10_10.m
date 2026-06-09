% close all
% clear all

mdsconnect('mpexserver')

shotlist = 10800 + [64,66,70]; 
shotlist = 10800 + [71:80]; 
shotlist = 10800 + [70,79]; 
shotlist = 10800 + [70,74,79]; 

address{1} = '\MPEX::TOP.MACHOPS1:PG4'; 
nPG =  address{1}(end);

% if nPG >= 4
%     f = 10;
% else
%     f = 2;
% end
f = 10; % use 10 for PG4 and 2 for all other

[PG,t] = my_mdsvalue_v2(shotlist,address(1));


%%
figure; hold on

for s = 1 :length(shotlist)
h(s) = plot(t{s}(1:length(PG{s})),PG{s}*f);
end
set(h(1),'LineWidth',2)
set(h(end),'LineWidth',2)

legend(h,{num2str(shotlist')},'location','NorthWest')
ylabel('mTorr')
ylim([-0.5,20])
xlim([3.5,5.5])
grid on
