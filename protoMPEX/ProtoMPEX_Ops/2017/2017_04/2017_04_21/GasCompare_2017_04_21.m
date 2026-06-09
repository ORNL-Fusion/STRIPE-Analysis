% close all
% clear all

mdsconnect('mpexserver')

% shotlist = 13700 + [94:97]; 
% shotlist = 13800+ [17,18,19];
% shotlist = 13894;
% shotlist = 13900 + [90,92];

% shotlist = [8129,14029,14040,14041]

address{1} = '\MPEX::TOP.MACHOPS1:PG3'; 
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
[PG0,t0] = my_mdsvalue_v2(14120,address(1));


%%
figure; hold on
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

for s = 1 :length(shotlist)
h(s) = plot(t{s}(1:length(PG{s})),(PG{s}-PG0{1})*f,C{s});
% h(s) = plot(t{s}(1:length(PG{s})),(PG{s})*f,C{s});
set(h(s),'LineWidth',2)
end
%set(h(1),'LineWidth',2)
% set(h(end),'LineWidth',2)
title(T)

legend(h,{num2str(shotlist')},'location','NorthWest')
ylabel('mTorr')
if f == 10
ylim([-0.5,20])
elseif f == 2
ylim([-0.5,4])
end
xlim([4,4.7])
grid on
set(gcf,'color','w')
box on
