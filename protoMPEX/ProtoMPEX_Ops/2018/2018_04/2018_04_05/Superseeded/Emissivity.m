clear all
close all
emissivity = linspace(0.5,1,5);

maxIntensity = 64726;
Res = 10000;
Inten = linspace(0,(maxIntensity),Res);
Tinterp = zeros(1,Res);
for ii=1:Res
    Tinterp1(ii) = seq.ThermalImage.GetValueFromEmissivity(emissivity(1),Inten(ii));
        Tinterp2(ii) = seq.ThermalImage.GetValueFromEmissivity(emissivity(2),Inten(ii));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp  = interp1(Inten,Tinterp,RawData);

figure; plot(Inten,Tinterp1)
hold on; plot(Inten,Tinterp2)