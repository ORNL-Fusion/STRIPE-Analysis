% Fetch and load data from MPEX server:
% =========================================================================

clear all
close all

fetchDataFromServer = 1;

if fetchDataFromServer
    % Radial scan:
    % =====================================================================
    % Entire list:    
    shotlist = 30000 + [224 ,225,227,229,230,231,232,233,234,235,238,236,239 ,240 ,241,244,245,246,247 ,248 ,249 ,250 ];
    r        =         [10.0,5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.5,0.0,3.75,3.25,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0];
    att      =         [10  ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10  ,10  ,10 ,10 ,10 ,10 ,10  ,10  ,10  ,10  ];

    
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

           
    save('Step_1_GetData_CapProbe_2020_06_17.mat',varList{:})
else
    load('Step_1_GetData_CapProbe_2020_06_17.mat')
end

%%

n = 8; figure; 
hold on
plot(t_vrms{n},sgolay_t(vrms{n},3,101))
plot(t_FWD{n},FWD{n}*50)
xlim([4.15,4.7])

figure
hold on
for ii = 1:numel(vrms)
    plot3(t_vrms{ii},r(ii)*ones(size(vrms{ii})),vrms{ii})
end
xlim([4.15,4.7])

