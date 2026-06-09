%%
% This script was sent to me by Cornwall
mdsconnect('mpexserver');

shots = [22601,22596,22595,22593,22592,22591,22576,22590,22575, ...
    22579,22574,22580,22581,22582,22583,22584,22585,22586,22587,22588,22589];
rad = zeros(numel(shots),1);
den_hel = rad;
temp_hel = rad;
den_ech = rad;
temp_ech = rad;
den_hel_err = rad;
temp_hel_err = rad;
den_ech_err = rad;
temp_ech_err = rad;
for ii = 1:numel(shots)
   mdsopen('mpex', shots(ii));
   rad(ii) = mdsvalue('ANALYZED.DLP.SETUP:RAD_LOC');
   temp_time = mdsvalue('ANALYZED.DLP:TIME');
   temp_den = mdsvalue('ANALYZED.DLP:NE');
   temp_temp = mdsvalue('ANALYZED.DLP:KTE'); 
   subs = find( (temp_time < 4.45 & temp_time > 4.4) |  ...
       (temp_time < 4.58 & temp_time > 4.53));
   den_hel(ii) = median(temp_den(subs));
   den_hel_err(ii) = std(temp_den(subs));
   temp_hel(ii) = median(temp_temp(subs));
   temp_hel_err(ii) = std(temp_temp(subs));
   subs = find( (temp_time < 4.505 & temp_time > 4.48));
   den_ech(ii) = median(temp_den(subs));
   den_ech_err(ii) = std(temp_den(subs));
   temp_ech(ii) = median(temp_temp(subs));
   temp_ech_err(ii) = std(temp_temp(subs));
end
%%
figure;
subplot(3,1,1); hold on;
errorbar(rad, den_hel, den_hel_err, '-ro');
errorbar(rad, den_ech, den_ech_err, '--b^');
legend('Density during helicon (m^-^3)', 'Density during ECH (m^-^3)');
subplot(3,1,2); hold on;
errorbar(rad, temp_hel, temp_hel_err, '-ro');
errorbar(rad, temp_ech, temp_ech_err, '--b^');
legend('Temperature during helicon (eV)', 'Temperature during ECH (eV)');
subplot(3,1,3); hold on;
plot(rad, temp_hel .* den_hel, '-ro');
plot(rad, temp_ech .* den_ech, '--b^');
legend('Pressure during helicon', 'Pressure during ECH');
xlabel('Radial distance (m)');