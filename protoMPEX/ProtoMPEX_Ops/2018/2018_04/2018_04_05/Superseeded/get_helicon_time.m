function [t_h_on, t_h_off] = get_helicon_time(shot)
    % JFC =========
    mdsconnect('mpexserver');
    % ====================
    mdsopen('mpex', shot);
    t_h_on = mdsvalue('ANALYZED.POWER:TSTART_HEL');
    t_h_off = mdsvalue('ANALYZED.POWER:TEND_HEL');
    mdsclose;
end