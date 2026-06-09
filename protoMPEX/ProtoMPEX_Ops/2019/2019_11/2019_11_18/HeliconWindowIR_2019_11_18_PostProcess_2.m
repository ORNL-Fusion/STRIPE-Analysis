% Postprocessor 2
% This scripts performs the following tasks:

clear all 
close all

% =========================================================================
% Select dataset
% =========================================================================
xp = 1;
SaveFig = 1;

switch xp
    case 1
        magConfig{1} = 'Window limit'
    case 2
        magConfig{2} = 'MPEX limit'
end

% =========================================================================
% Load data
% =========================================================================
load(['preprocessData_xp_',num2str(xp),'.mat'],'shot','rfPwrNet','viewType')

%%
% -------------------------------------------------------------------------
% FP data

try
    hf = findobj('Tag','fpdata');
    close(hf)
end

removeOffset = 1;
% To convert voltage to temperature, multiply by 20
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

% Gather data
addressFP{1} = [RootAddress,'FLUOROPT_1']; 
addressFP{2} = [RootAddress,'FLUOROPT_2']; 
addressFP{3} = [RootAddress,'FLUOROPT_3']; 
addressFP{4} = [RootAddress,'FLUOROPT_4']; 

for ii = 1:4
    [dataFP{ii},t_dataFP{ii}]   = my_mdsvalue_v2(shot,addressFP(ii));
end

for ii = 1:4
    for s = 1:length(shot)
        temperature_fp{ii}{s} = sgolay_t( ( dataFP{ii}{s}- removeOffset*mean(dataFP{ii}{s}(1:200)) )*20 ,3,51);
        temperature_fp_ss{ii}(s) = mean(temperature_fp{ii}{s}( (end-300):end) );
    end
end

figure
set(gcf,'Tag','dataFP','color','w')
hold on
[~,b] = sort(rfPwrNet);

for ii = 1:4
    subplot(2,2,ii); hold on
    for s = 1:length(shot)
        h_fp{ii}(s) = plot(t_dataFP{ii}{s}(1:end-1),temperature_fp{ii}{s},'LineWidth',3);
    end
    legend(h_fp{ii},num2str(rfPwrNet'))
    ylim([0,25])
end

figure; 
set(gcf,'Tag','dataFP','color','w')
for ii = 1:4
    subplot(2,2,ii);
    plot(rfPwrNet(b),temperature_fp_ss{ii}(b),'LineWidth',3)
    ylim([0,25])
    xlim([0,200])
    h_t{ii} = text(20,15,{magConfig{xp};['FP',num2str(ii)]});
    set(h_t{ii},'Interpreter','latex')
    ylabel('${\Delta}T$ [C]','Interpreter','latex')
    xlabel('RF [kW]','Interpreter','latex')
    grid on
end

subplot(2,2,2)
title('Ground side','Interpreter','latex')

subplot(2,2,3)
title('HV side','Interpreter','latex')

if SaveFig
    saveas(gcf,['FP vs RF ',magConfig{xp}],'tiffn')
end

%%
% -------------------------------------------------------------------------
% Gas flow and pressure data

try
    hf = findobj('Tag','gasData');
    close(hf)
end

addressGasFlow{1} = [RootAddress,'MFC_FLOW_D2'];
addressGasFlow{2} = [RootAddress,'MFC_FLOW_HE']; 

for ii = 1:2
    [dataGasFlow{ii},t_dataGasFlow{ii}]   = my_mdsvalue_v2(shot,addressGasFlow(ii));
end

addressGasPressure{1} = [RootAddress,'PG1']; 
addressGasPressure{2} = [RootAddress,'PG2']; 
addressGasPressure{3} = [RootAddress,'PG3']; 
addressGasPressure{4} = [RootAddress,'PG4']; 

shot_cal = 20908; % 2018_04_05, 2.0 kA case, 1 sec RF pulse
fctr = 0.1333;
gainPG = [2,2,2,10];
pressureGaugeLabel = {'12.5','8.5','6.5','2.5'};

for ii = 1:4
    [dataGasPressure{ii}    ,t_dataGasPressure{ii}    ]   = my_mdsvalue_v2(shot,addressGasPressure(ii));
    [dataGasPressure_cal{ii},t_dataGasPressure_cal{ii}]   = my_mdsvalue_v2(shot_cal,addressGasPressure(ii));

    for s = 1:length(shot)
            L1 = length(dataGasPressure{ii}{s}    );
            L2 = length(dataGasPressure_cal{ii}{1});
            if L1>L2
                Ldata = L2;
            else
                Ldata = L1;
            end
        P{ii}{s} = ( dataGasPressure{ii}{s}(1:Ldata)-dataGasPressure_cal{ii}{1}(1:Ldata) )*gainPG(ii)*fctr;
        t_P{ii}{s} = t_dataGasPressure{ii}{s}(1:Ldata);
    end
end
    
% figure
% for ii = 1:2
%     [dataGasFlow{ii},t_dataGasFlow{ii}]   = my_mdsvalue_v2(shot,addressGasFlow(ii));
% end
s = 7;
figure;
set(gcf,'Tag','gasData','color','w')
subplot(2,1,1)
hold on
for ii = 1:2
    plot(t_dataGasFlow{ii}{s}(1:end-1),dataGasFlow{ii}{s},'LineWidth',3)
end
xlim([3.8,5.5])
ylim([0,6])
ylabel('Piezo voltage [V]','Interpreter','latex','Fontsize',12)
box on
grid on
h_t = text(5.08,4.1,{magConfig{xp};['shot: ',num2str(shot(s))]});
set(h_t,'interpreter','Latex','fontsize',13)

subplot(2,1,2)
hold on
for ii = 1:4
    t_rng = find(t_P{ii}{s}>3.9 & t_P{ii}{s}<5.4);
    plot(t_P{ii}{s}(t_rng),P{ii}{s}(t_rng),'LineWidth',3)
end
xlim([3.8,5.5])
ylim([0,1.5])
ylabel('Pressure [Pa]','Interpreter','latex','Fontsize',12)
box on
grid on

if SaveFig
    saveas(gcf,['gas data ',magConfig{xp}],'tiffn')
end
%%
% -------------------------------------------------------------------------
% ne and Te data

try
    hf = findobj('Tag','dlpData');
    close(hf)
end

Stem = '\MPEX::TOP.'; Branch = 'ANALYZED.DLP:'; RootAddress = [Stem,Branch];

addressDLP{1} = [RootAddress,'NE']; 
addressDLP{2} = [RootAddress,'KTE']; 

for ii = 1:2
    [dataDLP{ii},~]   = my_mdsvalue_v2(shot,addressDLP(ii));
%     t_dataDLP{ii} = 1:
    t_rng = 120:180;
    for s = 1:length(shot)
        dataDLP_ss{ii}(s) = mean(dataDLP{ii}{s}(t_rng));
    end
end

figure;
set(gcf,'Tag','dlpData','color','w')
subplot(1,2,1)
hold on
for s = 1:length(shot)
    h_ne(s) = plot(dataDLP{1}{s},'LineWidth',3);
end
grid on
box on
ylabel('$n_e$ [m$^{-3}$]','Interpreter','latex','Fontsize',13)
xlim([0,350])
ylim([0,1.2e20])
l_ne = legend(h_ne,num2str(rfPwrNet'));
axis('square')
set(gca,'FontName','times','FontSize',11)

subplot(1,2,2)
hold on
for s = 1:length(shot)
    h_te(s) = plot(sgolay_t( dataDLP{2}{s},3,11) ,'LineWidth',3);
end
grid on
box on
ylabel('$T_e$ [eV]','Interpreter','latex','Fontsize',13)
h_t = text(50,7,{magConfig{xp}});
set(h_t,'interpreter','Latex','fontsize',13)
xlim([0,200])
ylim([0,8])
axis('square')
set(gca,'FontName','times','FontSize',11)

if SaveFig
    saveas(gcf,['dlp vs time data ',magConfig{xp}],'tiffn')
end


figure;
[~,b] = sort(rfPwrNet);
subplot(1,2,1)
set(gcf,'Tag','dlpData','color','w')
plot(rfPwrNet(b),dataDLP_ss{1}(b),'LineWidth',3)
grid on
box on
ylabel('$n_e$ [m$^{-3}$]','Interpreter','latex','Fontsize',13)
xlim([0,200])
ylim([0,1.2e20])
axis('square')
set(gca,'FontName','times','FontSize',11)

subplot(1,2,2)
set(gcf,'Tag','dlpData','color','w')
plot(rfPwrNet(b),dataDLP_ss{2}(b),'LineWidth',3)
grid on
box on
ylabel('$n_e$ [m$^{-3}$]','Interpreter','latex','Fontsize',13)
xlim([0,200])
ylim([0,7])
axis('square')
h_t = text(90,6,{magConfig{xp}});
set(h_t,'interpreter','Latex','fontsize',13)
set(gca,'FontName','times','FontSize',11)

if SaveFig
    saveas(gcf,['dlp vs RF data ',magConfig{xp}],'tiffn')
end
