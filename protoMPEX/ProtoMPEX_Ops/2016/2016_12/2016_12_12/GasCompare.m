close all
clear all

% Gather RF power data from the MDSserver
mdsconnect('mpexserver')
Shots = 13000 + [38,46,51,156];
Shots = 13000 + [6,8:29];
Shots = 00000 + [10907,10919,12991];
Shots = 9800 + [53,74:89];
Shots = 11422;
Shots = 13200 + [59,60]
Shots = [13275,13276,13277,13278,13279,13280,13281,13284];
Shots = 13200 + [75,80,81,88,89,90,91,98];
Shots = 13200 + [75,80,90,91,98,99,100];
Shots = 13200 + [80,99,100,101,102,104];
Shots = 12873;

Shots = 13300 + [19,20,21,23,26];
Shots = 13300 + [41];

% address{1} = '\MPEX::TOP.FSCOPE:TUBE08:PMT_VOLT';
% address{2} = '\MPEX::TOP.MPEX1:SPARE'; % 28 GHz signal
% address{3} = '\MPEX::TOP.MACHOPS1:TRANS_I'; % Helicon signal
% address{4} = '\MPEX::TOP.T_ZERO'; % Tzero
address{1} = '\MPEX::TOP.MACHOPS1:PG3'; % Helicon signal
f = 2;

[PG,t] = my_mdsvalue_v2(Shots,address(1));

%%
figure; hold on
for s = 1 :length(Shots)
h(s) = plot(t{s}(1:length(PG{s})),PG{s}*f)
end
plot(t{1}(1:length(PG{1})),(PG{1}-PG{2})*f + 0.48,'g')
legend(h,{num2str(Shots')})
ylabel('mTorr')
%xlim([4.0,4.7])
