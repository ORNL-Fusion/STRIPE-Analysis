function [Temp] = IntensityTempConv(emissivity,RawData,seq)
%Conversion from Intensity to Temperature in Celsius takes place here.
%First an interpolation table is made to speed up this process.
%%%%%%%%Create Interpolation Table%%%%%%%%%%%%%%%%%
maxIntensity = 64726;
Res = 1000;
Inten = linspace(0,(maxIntensity),Res);
Tinterp = zeros(1,Res);
for ii=1:Res
    Tinterp(ii) = seq.ThermalImage.GetValueFromEmissivity(emissivity,Inten(ii));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp  = interp1(Inten,Tinterp,RawData);
end

