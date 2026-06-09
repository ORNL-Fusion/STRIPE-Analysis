function IR_mdsplus_oneshot(shot, nowrite, mesh_size, time, ...
Px_per_cm, emissivity)
% this function runs the IR heat flux routine for a range of shots, and
% writes it to mdsplus
% Inputs
%      shot      shot number
% Optional inputs
% this can be used, or if not, historical records will be used
%      nowrite         1 (default) means to not write to mdsplus, 0 to
%      mesh_size       optional input for mesh_size (default to -1, which means use COMSOL file)
%      time            optional input for time vector in case default
%                      helicon finder routine does not work.
%      Px_per_cm       pixels per cm
%      Frame_rate      frame_rate
%      Emissivity      emissivity
    if nargin < 1
        shot = mdsvalue('current_shot("mpex")') - 1;
    end
    if nargin < 2
        nowrite = 1;
    end
    if nargin < 3
        mesh_size = -1;
    end
    % JFC ===========
    mdsconnect('mpexserver');
    % ================
    message = mdsopen('mpex', shot);  
    if ismember(shot, [20272,21629, 19625,19626, 21314:21850, 22949:23314, 24061:24258]) ...
            | (message ~= shot) | (message == 'Failed') %22949 to 23314 is dump plate
        disp('Bad shot numbers'); % known bad shot numbers, do nothing and skip
        mdsclose;
    else
        string = mdsvalue('\MPEX::TOP:T_ZERO');
        mdsclose;
        C = strsplit(string, {' ', '-'});
        folder = strcat('\\mpexserver\ProtoMPEX_Data\IR_Camera\','20',...
            C{3},'_',C{1},'_',C{2},'\');
        file_name = strcat(num2str(shot),'.seq');
        file_name = strjoin({'Shot', file_name});
        if exist(strcat(folder,file_name), 'file')
            [name, diameter, thickness, emissivity, density, ...
                pxx_per_cm, pxy_per_cm, spec_heat, ther_cond, ...
                heat_diff, temp_range,target_loc, comsol_file] = IR_history(shot);
            if nargin > 4
                pxx_per_cm = Px_per_cm;
                pxy_per_cm = Py_per_cm;
                emissivity = Emissivity;
            end
            [RawData,frames,seq, frame_rate] = ExtractRawData(folder,file_name,shot);
            [RawTemp] = IntensityTempConv(emissivity,RawData,seq);
            disp(['Max Temp is ', num2str(max(max(max(RawTemp))))]);
            if  nargin < 4
                [time,ind_s,ind_e] = CreateTimeVector ...
                    (RawTemp,shot,frame_rate,frames);
            else
                time = linspace(4.1, 4.6, 50);
                ind_s = 1;
                ind_e = 50;
            end
            if isempty(time)
                mdsopen('mpex', shot);
                mdsput('\MPEX::TOP.ANALYZED.IR:COMMENTS','Low helicon power');
                mdsclose;
            else
                if shot < 21911
                    rad_peak = 2.5;
                elseif shot > 21911 & shot < 24341
                    rad_peak = 1.75;
                else
                    rad_peak  = 1.75;
                end
                [xx_c,xx_s,xx_e,yy_c,yy_s,yy_e] = ...
                FindCenter(RawTemp,floor((ind_s+ind_e)/2),pxx_per_cm, pxy_per_cm, rad_peak);
                xx_s = max(1,ceil((xx_s - xx_c)*1.5) + xx_c);
                xx_e = min(floor((xx_e-xx_c)*1.5) +xx_c, numel(RawTemp(:,1,1)));
                yy_s = max(1,ceil((yy_s - yy_c)*1.5) + yy_c);
                yy_e = min(floor((yy_e-yy_c)*1.5) +yy_c, numel(RawTemp(1,:,1)));
                if time(ind_e) - time(ind_s) > .1
                    ind_s = max(1,ind_s-8); 
                    ind_e = min(numel(time),ind_e + 12);
                else
                    ind_s = max(1,ind_s-3); 
                    ind_e = min(numel(time),ind_e + 3);
                end
                time = time(ind_s:ind_e);
                Temp = RawTemp(xx_s:xx_e,yy_s:yy_e,ind_s:ind_e) + 273.15;
%                 Temp = fillmissing(Temp,'nearest');
                x = (xx_s-xx_c:xx_e-xx_c) / pxx_per_cm / 100.;
                y = (yy_s-yy_c:yy_e-yy_c) / pxy_per_cm / 100.;
                x_start = min(x); x_end = max(x); x_size = numel(x);
                y_start = min(y); y_end = max(y); y_size = numel(y);
                [y_grid,x_grid, t_grid] = meshgrid(y,x,time-min(time));
                M = [y_grid(:), x_grid(:), t_grid(:), Temp(:)];
                %csvwrite('C:\Users\icl\Documents\MATLAB\MPEX\automatic_analysis\IR_analysis\front_temperature.csv', M);
%                 mex_WriteMatrix('C:\Users\Public\front_temperature.csv' ...
%                     ,M,'%4.4f',',','w+');
                [x, y, Heat_Flux] = IR_comsol(mesh_size, time, x_start, ...
                    x_end, x_size, y_start, y_end, y_size, ...
                    Temp,frame_rate,comsol_file);
                Heat_Flux(isnan(Heat_Flux)) = 0;
                Heat_power = squeeze(sum(sum(Heat_Flux))) ...
                 * (x(2) - x(1)) * (y(2) - y(1));
                [Heat_peak, Heat_peak_temp] = max(Heat_Flux);
                [Heat_peak, Heat_peak_y] = max(Heat_peak);
                Heat_peak = squeeze(Heat_peak);
                Heat_peak_y = squeeze(Heat_peak_y);
                Heat_peak_x = Heat_peak_y;
                for ii=1:numel(Heat_peak_x)
                    Heat_peak_x(ii) = Heat_peak_temp(1,Heat_peak_y(ii), ii);
                end
                Heat_peak_y = y(Heat_peak_y);
                Heat_peak_x = x(Heat_peak_x);
            
                clear Heat_peak_temp;
                [xx_c,xx_s,xx_e,yy_c,yy_s,yy_e] = ... 
                FindCenter(Heat_Flux,ceil(numel(time)/2), pxx_per_cm, pxy_per_cm, rad_peak);
                x = x - x(xx_c);
                y = y - y(yy_c);
                if shot > 22631
                   x = x - .006; 
                end
                rad = FindRadius(x,y,Heat_Flux,numel(time)); rad = rad(:,1);
                [Heat_peak_ce, Heat_peak_ed] = FindPeak(x,y,Heat_Flux, numel(time));
                mdsopen('mpex', shot);
                t_on = mdsvalue('ANALYZED.POWER:TSTART_HEL');
                mdsclose;
                if max(Heat_power) > 500
                    t_IR_on = time(find(Heat_power > 50, 1, 'first')) - .02;
                elseif max(Heat_power) > 100
                    t_IR_on = time(find(Heat_power > 20, 1, 'first')) - .02;
                else
                    t_IR_on = t_on;
                end
                time = time - t_IR_on + t_on;
                if ~nowrite
                    mdsopen('mpex', shot); 
                    mdsput('\MPEX::TOP.ANALYZED.IR:COMMENTS', '$', '');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_FLUX', '$', Heat_Flux);
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_PEAK', '$', Heat_peak');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_PEAK_X', '$', Heat_peak_x');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_PEAK_Y', '$', Heat_peak_y');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_POWER', '$', Heat_power');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_PEAK_CE', '$', Heat_peak_ce');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_PEAK_ED', '$', Heat_peak_ed');
                    mdsput('\MPEX::TOP.ANALYZED.IR:HEAT_RADIUS', '$', rad');
                    mdsput('\MPEX::TOP.ANALYZED.IR:T', '$', time');
                    mdsput('\MPEX::TOP.ANALYZED.IR:X', '$', x');
                    mdsput('\MPEX::TOP.ANALYZED.IR:Y', '$', y');
                    mdsput('\MPEX::TOP.ANALYZED.IR:SETUP:EMISSIVITY', '$', emissivity);
                    mdsput('\MPEX::TOP.ANALYZED.IR:SETUP:PX_PER_CM', '$', pxx_per_cm);
                    mdsput('\MPEX::TOP.ANALYZED.IR:SETUP:FRAME_RATE', '$', frame_rate);
                    mdsput('\MPEX::TOP.ANALYZED.IR:SETUP:SURF_TEMP', '$', Temp);
                    mdsclose;
                else
                    plot(time, Heat_peak_ce); title('Heat Peak'); hold on;
                    plot(time, Heat_peak_ed);
                    plot(time, Heat_peak);
                    disp('write option was not specified.');
                end
            end
        else
            comments = 'No IR file';
            mdsopen('mpex', shot);
            mdsput('\MPEX::TOP.ANALYZED.IR:COMMENTS',comments);
            mdsclose;
        end
    end
end