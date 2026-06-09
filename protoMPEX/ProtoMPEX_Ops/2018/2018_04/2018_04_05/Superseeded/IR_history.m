function [name, diameter, thickness, emissivity, density, pxx_per_cm, ...
    pxy_per_cm, spec_heat, ther_cond, heat_diff, temp_range, ...
    target_loc, comsol_file] = IR_history(shot)
% This function returns the historical records for the target as a function
% of mdsplus shot. All units in mks
%   Inputs
%         shot           shot number
%   Outputs
%         too many...    all the stuff you need for the target
    mdsopen('mpex', shot);
    string = mdsvalue('\MPEX::TOP:T_ZERO');
    mdsclose;   
    C = strsplit(string, {' ', '-'});
    shot_sum = strcat('\\mpexserver\Shot Summaries\20', C{3}, '_', C{1}, ...
         '\20', C{3}, '_',C{1},'_',C{2},'_OperShotSummary.xlsx');
    switch shot
        case num2cell(15000:15713) %stainless steel plate
            name = '304 SS';
            diameter = .1143;
            thickness = .000254;
            emissivity = .55;
            density = 8030;
            temp_range = [0,500,1000]; 
            spec_heat = [502.5, 597.1,696.3];
            ther_cond = [14.9,21.8,28.8];
            heat_diff = [3.68e-6,4.55e-6,5.15e-6];
            target_loc = 4.2071;
            pxx_per_cm = 24.14;
            pxy_per_cm = pxx_per_cm;
            comsol_file = ['C:\Users\Public\', ...
                'IRcamera_stainless_201707.mph'];
        case num2cell(15000:18814) %Graphite (ATJ) plate
            name = 'Graphite (ATJ)';
            diameter = .1143; 
            thickness = .003175;
            emissivity = .75;
            density = 1760;
            temp_range = [0,500,1000]; 
            spec_heat = [700,1600,2000];
            ther_cond = [100,61,46];
            heat_diff = [8.1e-5,2.2e-5,1.3e-5];
            target_loc = 4.2071;
            pxx_per_cm = 24.14;
            pxy_per_cm = pxx_per_cm;
            comsol_file = ['C:\Users\Public\', ...
                'IRcamera_graphite_201708.mph'];
        case num2cell(18815:21275) %304SS with embedded Thermocoax cable
            name = '304 SS';
            diameter = .1143;
            thickness = .001524;
            emissivity = .55;
            density = 8030;
            temp_range = [0,500,1000]; 
            spec_heat = [502.5, 597.1,696.3];
            ther_cond = [14.9,21.8,28.8];
            heat_diff = [3.68e-6,4.55e-6,5.15e-6];
            comsol_file = ['C:\Users\Public\', ...
                'IRcamera_stainless_201801.mph'];
            if nargin ==2 
                fileinfo = xlsread(shot_sum,'ShotSummary', ['A8:BE300']);
                index1 = find(fileinfo(:,1) == shot);
                if isempty(index1) 
                    %cannot find shot, so find latest shot to use.
                    index1 = numel(num_data(:,1));
                end
                target_loc = fileinfo(index1,57);
                switch target_loc
                    case 2
                        pxx_per_cm = 23.61;
                    case 0
                        pxx_per_cm = 24.46;
                    case -1
                        pxx_per_cm = 24.3;
                    case -1.25
                        pxx_per_cm = 24.2;
                    case -1.5
                        pxx_per_cm = 24.1;
                    case -1.75
                        pxx_per_cm = 24.;
                    case -2
                        pxx_per_cm = 23.79;
                    case -3
                        pxx_per_cm = 23.69;
                    case -4
                        pxx_per_cm = 23.38;
                    case -5
                        pxx_per_cm = 22.55;
                    case -6
                        pxx_per_cm = 22.55;
                    otherwise % assume it's 0
                        pxx_per_cm = 23.5;  
                end
                target_loc = target_loc*.01 + 4.2071;
            else
                target_loc = 4.2071;
                pxx_per_cm = 24.;
            end
            pxy_per_cm = pxx_per_cm;
        case num2cell(21276:21910) %SiC plate
            name = 'SiC';
            diameter = .1143;
            thickness = .003175;
            emissivity = .75;
            density = 3210;
            temp_range = [300,300];
            spec_heat = [750,750];
            ther_cond = [300,300];
            heat_diff = [1.2e-4,1.2e-4];
            comsol_file = ['C:\Users\Public\IRcamera_SiC_201804.mph'];
            if nargin == 2
                fileinfo = xlsread(shot_sum,'ShotSummary', ['A8:BE300']);
                index1 = find(fileinfo(:,1) == shot);
                if isempty(index1) 
                    %cannot find shot, so find latest shot to use.
                    index1 = numel(num_data(:,1));
                end
                target_loc = fileinfo(index1,57)*.01 + 4.2071;
                pxx_per_cm = 23.65;
            else
                target_loc = 4.2071;
                pxx_per_cm = 23.65;
            end    
            pxy_per_cm = pxx_per_cm;
        case num2cell(21911:22115)
            name = '304 SS';
            diameter = .05;
            thickness = .000254;
            emissivity = .4;
            density = 8030;
            temp_range = [0,500,1000]; 
            spec_heat = [502.5, 597.1,696.3];
            ther_cond = [14.9,21.8,28.8];
            heat_diff = [3.68e-6,4.55e-6,5.15e-6];
            target_loc = 4.1316;
            pxx_per_cm = 50;
            pxy_per_cm = 27.;
            comsol_file = ['C:\Users\Public\', ...
                'IRcamera_stainless_201808.mph'];
        case num2cell([22116:27000])
            name = '304 SS';
            diameter = .1;
            thickness = .001524;
            emissivity = .86;
            density = 8030;
            temp_range = [0,500,1000]; 
            spec_heat = [502.5, 597.1,696.3];
            ther_cond = [14.9,21.8,28.8];
            heat_diff = [3.68e-6,4.55e-6,5.15e-6];
            target_loc = 4.1316;
            if shot < 23100
                pxx_per_cm = 50.;
                pxy_per_cm = 27.;
            else
                pxx_per_cm = 63;
                pxy_per_cm = 32.;
            end
            comsol_file = ['C:\Users\Public\', ...
                'IRcamera_stainless_201809.mph'];
        otherwise
            error(['Did not recognize shot: ',num2str(shot)]);
    end
end