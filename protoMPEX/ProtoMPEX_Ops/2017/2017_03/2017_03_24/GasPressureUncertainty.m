% Gas pressure uncertainty on March 24th 2017

close all
clear all

shot    = 13512; % PS1/PS2 = 4000 A, 200 ms RF pulse 
% shot    = 13480; % PS1/PS2 = 4000 A, 200 ms RF pulse 
shotCal = 13500 + [18:21];

mdsconnect('mpexserver')

% SP6.5:
address{1} = '\MPEX::TOP.MACHOPS1:PG3'; % 6.5
[a,ta] = my_mdsvalue_v2(shot  ,address(1));
[b,~] = my_mdsvalue_v2(shotCal,address(1));
% SP4.5:
address{1} = '\MPEX::TOP.MACHOPS1:PG4'; % 4.5
[c,tc] = my_mdsvalue_v2(shot,address(1));
[d,~] = my_mdsvalue_v2(shotCal,address(1));

%%
for sc = 1:length(shotCal)
    B(:,sc) = b{sc};
    D(:,sc) = d{sc};
end
PG6_0m = 2*mean(B,2);
d_PG6_0 = 2*std(B,1,2);
PG6   = 2*a{1}-PG6_0m;
PG6_u = PG6 + 0.5*d_PG6_0;
PG6_l = PG6 - 0.5*d_PG6_0;
t6 = ta{1}(1:length(PG6));

PG4_0m = 10*mean(D,2);
d_PG4_0 = 10*std(D,1,2);
PG4 = 10*c{1}-PG4_0m;
PG4_u = PG4 + 0.5*d_PG4_0;
PG4_l = PG4 - 0.5*d_PG4_0;
t4 = tc{1}(1:length(PG4));

close all
figure; 
subplot(2,1,1)
hold on
plot(t4,PG4,'k','LineWidth',2)
plot(t4,PG4_u,'r')
plot(t4,PG4_l,'r')
xlim([4,4.6])
ylim([-0.2,8])
box on
ylabel('[mTorr]')
grid on

subplot(2,1,2)
hold on
plot(t6,PG6,'k','LineWidth',1)
plot(t6,PG6_u,'r')
plot(t6,PG6_l,'r')
xlim([4,4.6])
ylim([-0.1,1.5])
box on
ylabel('[mTorr]')
grid on

set(gcf,'color','w')
box on




% PG6 = (a{1}-mean(a{1}(1:100)))*2;
% PG4 = (a{1}-mean(a{1}(1:100)))*10;