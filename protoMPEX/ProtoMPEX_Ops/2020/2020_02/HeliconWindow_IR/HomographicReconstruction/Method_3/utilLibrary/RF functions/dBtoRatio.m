function [P_ratio,V_ratio] = dBtoRatio(dB)
% Decibels to power and voltage ratio

P_ratio = 10.^(dB/10)
V_ratio = 10.^(dB/20)

end

