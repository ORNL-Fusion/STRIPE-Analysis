function [Temp_mean,dTemp] = IntensityTempConv2(emissivity,emissivity_std,RawData,seq)
%Conversion from Intensity to Temperature in Celsius takes place here.
%First an interpolation table is made to speed up this process.
%%%%%%%%Create Interpolation Table%%%%%%%%%%%%%%%%%
maxIntensity = 64726;
Res = 10000;
Inten = linspace(0,(maxIntensity),Res);
Tinterp_mean = zeros(1,Res);
for ii=1:Res
    Tinterp_mean(ii) = seq.ThermalImage.GetValueFromEmissivity(emissivity,Inten(ii));
    Tinterp_U(ii) = seq.ThermalImage.GetValueFromEmissivity(emissivity-emissivity_std,Inten(ii));
    Tinterp_L(ii) = seq.ThermalImage.GetValueFromEmissivity(emissivity+emissivity_std,Inten(ii));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp_mean  = interp1(Inten,Tinterp_mean,RawData);
Temp_U  = interp1(Inten,Tinterp_U,RawData);
Temp_L  = interp1(Inten,Tinterp_L,RawData);
dTemp = abs(Temp_U - Temp_L);
end

