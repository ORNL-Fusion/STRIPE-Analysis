% close all
% clear all

mdsconnect('mpexserver')

% shotlist = 13700 + [94:97]; 
% shotlist = 13800+ [17,18,19];

address{1} = '\MPEX::TOP.MACHOPS1:PG3'; 
nPG =  address{1}(end);

if strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG3')
     f = 2;
     T = 'PG6.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG4')
     f = 10;
     T = 'PG4.5';
 end

[PG,t] = my_mdsvalue_v2(shotlist,address(1));
[PG0,t0] = my_mdsvalue_v2(13894,address(1));


%%
figure; hold on

for s = 1 :length(shotlist)
h(s) = plot(t{s}(1:length(PG{s})),(PG{s}-PG0{1})*f);
end
set(h(1),'LineWidth',2)
set(h(end),'LineWidth',2)
title(T)

legend(h,{num2str(shotlist')},'location','NorthWest')
ylabel('mTorr')
if f == 10
ylim([-0.5,20])
elseif f == 2
ylim([-0.5,5])
end
xlim([4.0,4.7])
grid on
set(gcf,'color','w')
box on
