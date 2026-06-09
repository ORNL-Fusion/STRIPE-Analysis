% Plot the times of the shots on 2017_04_11
clear all
close all


shotlist = 13000 + [794:870];

% shotlist = 13000 + [794:800];

for ii = 1:numel(shotlist)
    d = t_zero(shotlist(ii));
    hour(ii) = str2num(d{1}(10:11));
    min(ii)  = str2num(d{1}(13:14));
end

for ii = 1:numel(shotlist)
    dt(ii) = hour(ii)*60 + min(ii)
end

figure
plot(dt,'ko')