function [ Vr ] = Vratio(CalFile, Vmag)

cal = csvread(CalFile,6,0);
P = polyfit(cal(:,2), cal(:,1), 1);

A = polyval(P,Vmag);
Vr = 10.^(A./20);
end

