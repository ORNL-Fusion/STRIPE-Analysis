% the purpose of this script is to plot the gas flow from the MFC for shots
% associated with the paper on highly ionized plasmas in PRotoMPEX

clear all
close all

%% Gas flow, no plasma, no magnets
% MFC gas only shot
shot = 13969;
% the gas puff is done upstream at spool 2.5

% Mass flow controller voltage ============================================
address{1} = '\MPEX::TOP.MACHOPS1:MFC_FLOW_D2';
[s_MFC,t_MFC] = my_mdsvalue_v2(shot,address(1));
% Postprocess data
s_SLM = (s_MFC{1}*2000 -100)/1000; % in SLM
MFC{1} = sgolay_t(s_SLM,3,11);
MFC_offset{1} = mean(MFC{1}(find(t_MFC{1}<4)));


% Gas pressure at spool 2.5
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
[s_B2,t_B2] = my_mdsvalue_v2(shot,address(2));
[s0_B2,t0_B2] = my_mdsvalue_v2(13968,address(2));

PG2{1} = (s_B2{1}-s0_B2{1})*2*0.13333; % in Pascals

figure;
hold on;
h(1) = plot(t_MFC{1}(1:end-1),MFC{1}-MFC_offset{1})
h(2) = plot(t_B2{1}(1:end-1),PG2{1})
xlim([3.5,6])
ylim([-0.2,3])
xlabel('time [s]')
grid on
box on
set(gcf,'color','w')
legend(h,'MFC in SLM','PG2 in Pa')

% NOTE:
% in software we request the gas puff to begin at t = 3.85s. The mass flow
% controller begins to register gas flow at t = 3.94s that is 90 ms later.
% the gas pressure at spool 2.5 begins to rise at t = 4.04 s

%% Gas injection relative to RF pulse
% Gas pressure at spool 2.5
shot = 13962;
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
[s_B2,t] = my_mdsvalue_v2(shot,address(2));
t_B2{2} = t{1};
[s0_B2,t0_B2] = my_mdsvalue_v2(13968,address(2));
PG2{2} = (s_B2{1}-s0_B2{1})*2*0.13333; % in Pascals

% RF power trace
address{3} = '\MPEX::TOP.MACHOPS1:RF_FWD_PWR'; % RF power trace
[s_RF,t] = my_mdsvalue_v2(shot,address(3));
RF{2} = -s_RF{1};
t_rf{2} = t{1};

figure;
hold on;
h(1) = plot(t_MFC{1}(1:end-1),MFC{1}-MFC_offset{1})
h(2) = plot(t_B2{2}(1:end-1),PG2{2})
h(3) = plot(t_B2{1}(1:end-1),PG2{1})
h(4) = plot(t_rf{2}(1:end-1),3*RF{2})

xlim([3.5,6])
ylim([-0.2,3])
xlabel('time [s]')
grid on
box on
set(gcf,'color','w')
legend(h,'MFC in SLM','PG2 in Pa')
