function [x, y, Heat_Flux] = IR_comsol(mesh_size, time, x_start, x_end, ...
    x_size, y_start, y_end, y_size, Temp,frame_rate,comsol_file);
% this function runs the specified comsol file to get the heat flux fropm
% the surface temperature
% Inputs
%       time                   1-D time array
%       Temp                   3-D surface temperature array in X,Y,T
%       frame_rate             frame rate of IR camera
%       comsol_file            string for which comsol file to run
% Outputs
%       Heat_Flux    3-D heat flux array in X,Y,T

    t_start = 0; t_end = max(time)-min(time); t_size = numel(time);  
    ambient_temp = median(median(Temp(:,:,1)));
    try
        comsol_model = mphload(comsol_file);
    catch
        import com.comsol.model.*
        import com.comsol.model.util.*
        try
            mphstart
        catch
        end
        comsol_model = mphload(comsol_file);
    end
    comsol_model.param.set('t_start', t_start);
    comsol_model.param.set('t_end', t_end);
    comsol_model.param.set('t_size', t_size);
    comsol_model.param.set('x_start', x_start);
    comsol_model.param.set('x_end', x_end);
    comsol_model.param.set('x_size', x_size);
    comsol_model.param.set('y_start', y_start);
    comsol_model.param.set('y_end', y_end);
    comsol_model.param.set('y_size', y_size);
    comsol_model.param.set('frame_rate', frame_rate);
    comsol_model.param.set('time_step', 1 / frame_rate);
    if mesh_size ~= -1
        comsol_model.param.set('mesh_size', mesh_size);
    else
        % do nothing
    end
    comsol_model.param.set('ambient_temperature', ambient_temp);
    comsol_model.study('std1').run;
    M = csvread('C:\Users\Public\heatflux.csv', 9);    
    x = unique(M(:,1));
    y = unique(M(:,2));
    M = M(:,3:end);
    Heat_Flux = reshape(M, numel(x), numel(y), t_size);
end