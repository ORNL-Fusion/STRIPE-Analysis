close all
clear all

mdsconnect('mpexserver')
Shots = 13200 + [80,99,100,101,102,104];

% address{1} = '\MPEX::TOP.FSCOPE:TUBE08:PMT_VOLT';
% address{2} = '\MPEX::TOP.MPEX1:SPARE'; % 28 GHz signal
% address{3} = '\MPEX::TOP.MACHOPS1:TRANS_I'; % Helicon signal
% address{4} = '\MPEX::TOP.T_ZERO'; % Tzero

address{1} = '\MPEX::TOP.MACHOPS1:PG4';


[PG,t] = my_mdsvalue_v2(Shots,address(1));

figure; hold on
for s = 1 :length(Shots)
h(s) = plot(t{s},PG{s}*10)
end
legend(h,{num2str(Shots')})
ylabel('mTorr')
