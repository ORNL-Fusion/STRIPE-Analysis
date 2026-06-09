% Fetch and load data from MPEX server:
% =========================================================================

clear all
close all

fetchDataFromServer = 0;

if fetchDataFromServer
    % Radial scan:
    % =====================================================================
    % Entire list:
    shotlist =  [31000 + [821,822,823,824,825,826]];
    r        =           [5.5,4.5,3.5,2.5,1.5,0  ] ;
    att      =           [10 ,10 ,10 ,10 ,10 ,10 ] ;
    
    % Address for the bdot probe data:
    % =====================================================================
    % Capacitive:
    DataAddress{1} = ['\MPEX::TOP.MACHOPS1:BDOT_10'];
    
    % DLP 12.5 ion sat
    DataAddress{2} = ['\MPEX::TOP.MACHOPS1:LP_1'];

    % Fwd RF power trace:
    DataAddress{3} = ['\MPEX::TOP.MACHOPS1:RF_FWD_PWR'];
    
    % reflected RF power trace:
    DataAddress{4} = ['\MPEX::TOP.MACHOPS1:RF_REF_PWR'];
    
    % Load data from server:
    % =====================================================================
    for ii = 1:numel(shotlist)
            for ch = 1:4
                [f{ch},t_f{ch}] = my_mdsvalue_v2(shotlist(ii),DataAddress(ch));   
            end

            % Capacitive probe data:
            vcap{ii} = -f{1}{1};
            for jj = 1:numel(vcap{ii})
                if vcap{ii}(jj)>0
                    vrms{ii}(jj) = +40.95.*(vcap{ii}(jj)).^(0.503);
                else
                    vrms{ii}(jj) = -40.95.*(abs(vcap{ii}(jj))).^(0.503);
                end
            end
            t_vcap{ii} = t_f{1}{1};
            t_vrms{ii} = t_vcap{ii}(1:end-1);
            
            % Isat:
            isat{ii}   = f{2}{1};
            t_isat{ii} = t_f{2}{1}(1:end-1);
            
            % RF_FWD:
            FWD{ii}   = f{3}{1};
            t_FWD{ii} = t_f{3}{1}(1:end-1);
            
            % RF_REF:
            REF{ii}   = f{4}{1};
            t_REF{ii} = t_f{4}{1}(1:end-1);
    end
    
    varList = {'shotlist','r','att',...
               'DataAddress',...
               'vcap','t_vcap',...
               'vrms','t_vrms',...
               'isat','t_isat',...
               'FWD','t_FWD',...
               'REF','t_REF'};

           
    save('Step_1_GetData_CapProbe_2021_04_15.mat',varList{:})
else
    load('Step_1_GetData_CapProbe_2021_04_15.mat')
end

%%

n = 1; figure; 
hold on
plot(t_vrms{n},sgolay_t(vrms{n},3,101))
plot(t_FWD{n},FWD{n}*50)
xlim([4.15,4.7])

figure('color','w')
hold on
for ii = 1:numel(vrms)
    plot3(t_vrms{ii},r(ii)*ones(size(vrms{ii})),vrms{ii})
end
xlim([4.15,4.7])
xlabel('time [s]')
ylabel('R [cm]')
zlabel('V_{rms} [V]')
view([30,30])
title('Capacitive probe, central chamber')
