close all
clear all

mdsconnect('mpexserver')

shotlist = 14000 + [16,17,18,19];

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; 
nPG =  address{1}(end);

if strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG3')
     f = 2;
     T = 'PG6.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG1')
     f = 2;
     T = 'PG9.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG4')
     f = 10;
     T = 'PG4.5';
 elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG2')
    f = 2;
    T = 'PG2.5';
 end
% f = 2; % use 10 for PG4 and 2 for all other

[PG,t] = my_mdsvalue_v2(shotlist,address(1));
% [PG0,t0] = my_mdsvalue_v2(13894,address(1));
[PG0,t0] = my_mdsvalue_v2(13968,address(1));


%%
figure; hold on
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

for s = 1 :length(shotlist)
h(s) = plot(t{s}(1:length(PG{s})),(PG{s}-PG0{1})*f,C{s});
set(h(s),'LineWidth',2)
end
set(h(1),'LineWidth',2)
set(h(end),'LineWidth',2)
title(T)

legend(h,{num2str(shotlist')},'location','NorthWest')
ylabel('mTorr')
if f == 10
ylim([-0.5,20])
elseif f == 2
ylim([-0.5,20])
end
xlim([0,8])
grid on
set(gcf,'color','w')
box on
