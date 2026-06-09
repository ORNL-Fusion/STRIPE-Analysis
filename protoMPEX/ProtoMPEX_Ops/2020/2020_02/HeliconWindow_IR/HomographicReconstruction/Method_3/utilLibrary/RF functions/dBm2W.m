function [W] = dBm2W(dBm)
% this function converts, dBm into W
% see this reference: http://www.giangrandi.ch/electronics/anttool/decibel.shtml
W = 10^((dBm-30)/10);


end

