function [heat_flux, heat_peak, heat_peak_x, heat_peak_y, heat_power ...
    heat_radius, T, X, Y, comments] = get_IR(shot)
% Inputs
%     shot               single shot number
% Outputs
%     heat_flux          heat flux 3-D array in X,Y,T
%     heat_peak          peak heat flux, 1-D array in T
%     heat_peak_x        x location of peak heat flux, 1-D array in T
%     heat_peak_y        y location of peak heat flux, 1-D array in T
%     heat_power         power to target 1-D array in T
%     heat_radius        radius of heat flux 1-D array in T
%     T                  time 1-D array
%     X                  x 1-D array
%     Y                  Y 1-D array
%     comments           comments

    if nargin < 1
        shot = mdsvalue('current_shot("mpex")') - 1;
    end
    message = mdsopen('mpex', shot);
    T = mdsvalue('ANALYZED:IR:T');
    X = mdsvalue('ANALYZED:IR:X');
    Y = mdsvalue('ANALYZED:IR:Y');
    heat_flux = mdsvalue('ANALYZED:IR:HEAT_FLUX');
    heat_flux = permute(heat_flux, [2,1,3]);
    heat_peak = mdsvalue('ANALYZED:IR:HEAT_PEAK');
    heat_peak_x = mdsvalue('ANALYZED:IR:HEAT_PEAK_X');
    heat_peak_y = mdsvalue('ANALYZED:IR:HEAT_PEAK_Y');
    heat_power = mdsvalue('ANALYZED:IR:HEAT_POWER');
    heat_radius = mdsvalue('ANALYZED:IR:HEAT_RADIUS');
    comments = mdsvalue('ANALYZED:IR:COMMENTS');
    mdsclose;
end