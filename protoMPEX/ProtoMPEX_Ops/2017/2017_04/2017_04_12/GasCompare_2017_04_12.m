% close all
% clear all

mdsconnect('mpexserver')

% shotlist = 13700 + [94:97]; 
% shotlist = 13800+ [17,18,19];
% shotlist = 13894;

address{1} = '\MPEX::TOP.MACHOPS1:PG3'; 
nPG =  address{1}(end);

% if nPG >= 4
%     f = 10;
% else
%     f = 2;
% end
f = 10; % use 10 for PG4 and 2 for all other

[PG,t] = my_mdsvalue_v2(shotlist,address(1));
[PG0,t0] = my_mdsvalue_v2(13894,address(1));


%%
figure; hold on

for s = 1 :length(shotlist)
h(s) = plot(t{s}(1:length(PG{s})),(PG{s}-PG0{1})*f);
end
set(h(1),'LineWidth',2)
set(h(end),'LineWidth',2)

legend(h,{num2str(shotlist')},'location','NorthWest')
ylabel('mTorr')
ylim([-0.5,20])
xlim([3.5,5.5])
grid on
