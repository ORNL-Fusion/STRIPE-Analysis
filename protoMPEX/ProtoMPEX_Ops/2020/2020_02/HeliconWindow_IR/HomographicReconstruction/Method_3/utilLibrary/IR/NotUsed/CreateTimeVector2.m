function [time,ind_s,ind_e] = CreateTimeVector2(RawData,shot,framerate,frames)
% Creates time vector that corresponds to helicon pulse time vector
% Note this function is critically important in getting many diagnostics to
% work. This is because a lot of diagnostics are not on the mdsplus
% timebase and need to sync up appropiately to do so. 
%  Inputs
%       RawData     3-D array (X,Y, T)usually heat flux or raw data
%       shot        shot number
%       framerate   frame rate of IR camera
%       frames      # of frames in IR camera
%  Output
%       time        1-D synced up time array to MDSplus
%       ind_s       starting index when helicon turns on
%       ind_e       ending index when helicon turns off
    [t_h_on, t_h_off] = get_helicon_time(shot);
    if ~isempty(t_h_on) | ~isempty(t_h_off)
        tme = movmean(diff(diff(squeeze(max(max(RawData))))),5);
        
% =========================================================================
% JF Caneses 2019_01_10 The function "isoutlier" is not available on my
% Matlab version, so I had to modify the existing code to make it work in
% my comptuer
% The following is the original code
% =========================================================================
%         ind_s = find(isoutlier(tme, 'ThresholdFactor', 5), 1, 'First'); 
% =========================================================================
% The following is the new addition which uses "peekseek" which is a code
% that I got from the Matlab central
% =========================================================================
        [a,~] = peakseek(tme,40);
        ind_s = a(1);
% =========================================================================
% =========================================================================

        time = (0.:1. /framerate:(frames-1)/framerate);
        if ~isempty(ind_s) & isnumeric(t_h_on) & ~isempty(time)
            time = time + t_h_on  - time(ind_s);
            [test, ind_e] = min(abs(time-t_h_off));
        else
            time = []; ind_s = []; ind_e = [];
        end
    else
        time = []; ind_s = []; ind_e = [];
end

