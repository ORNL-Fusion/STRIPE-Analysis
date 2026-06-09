% close all
% clear all

mdsconnect('mpexserver')
 
% shotlist = 13500 + [18:23,24]; % B0 calibration shots at 4000 A

address{1} = '\MPEX::TOP.MACHOPS1:PG3'; 
nPG =  address{1}(end);

% if nPG >= 4
%     f = 10;
% else
%     f = 2;
% end
f = 2; % use 10 for PG4 and 2 for all other

[PG,t] = my_mdsvalue_v2(shotlist,address(1));

% B0 only calibration shots
[PG0a,t0] = my_mdsvalue_v2([13469],address(1)); % 4000 A
[PG0b,t0] = my_mdsvalue_v2([13470],address(1)); % 4400 A

% Gas only calibration shots
[PG1,t1] = my_mdsvalue_v2([13471],address(1));

% Reference 250 ms pulse
[PG0L,t0L] = my_mdsvalue_v2([12975],address(1));


%%
figure; hold on
offset = 0*ones(size(shotlist));
offset([1:2]) = -0*0.5;
offset([3]) = -0;

for s = 1 :length(shotlist)
h(s) = plot(t{s}(1:length(PG{s})),(PG{s}-1*PG0a{1})*f + offset(s));
end
set(h(1:2),'LineWidth',2)
% Long reference pulse
plot(t0L{1}(1:length(PG0L{1})),(PG0L{1}-1*mean(PG0L{1}(1:200)))*f + offset(1),'g','LineWidth',2)
%plot(t0L{1}(1:length(PG0L{1})),(PG0L{1}-1*PG0m)*f,'g','LineWidth',2)

legend(h,{num2str(shotlist')},'location','NorthWest')
ylabel('mTorr')
ylim([-0.5,5])
xlim([3.5,5.5])
grid on
