close all;
clear all;

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;
Z = 10;

% % Ionization
% file_inz = 'scd89_ne.dat';
% % Recombination
% file_rcmb = 'acd89_ne.dat';
% % Charge Exchange
% file_cx = 'ccd89_ne.dat';

% Ionization
file_inz = 'scd89_ne.dat';
% Recombination
file_rcmb = 'acd89_ne.dat';
% Charge Exchange
file_cx = 'ccd89_ne.dat';

[IonizationTemp, IonizationDensity, IonizationRateCoeff, IonizationChargeState] = ADF11s(file_inz);
[RecombinationTemp, RecombinationDensity, RecombinationRateCoeff, RecombinationChargeState] = ADF11a(file_rcmb);
[CXTemp, CXDensity, CXRateCoeff, CXChargeState] = ADF11_CX(file_cx);  % ✅ Read Charge Exchange data

% ✅ Convert Data to Log Scale (like Ionization & Recombination)
IonizationData.Temp = IonizationTemp;
IonizationData.Density = log10(10.^IonizationDensity.*1e6);
IonizationData.RateCoeff = log10(10.^IonizationRateCoeff./1e6);
IonizationData.ChargeState = IonizationChargeState;

RecombinationData.Temp = RecombinationTemp;
RecombinationData.Density = log10(10.^RecombinationDensity.*1e6);
RecombinationData.RateCoeff = log10(10.^RecombinationRateCoeff./1e6);
RecombinationData.ChargeState = RecombinationChargeState;

CXData.Temp = CXTemp;
CXData.Density = log10(10.^CXDensity.*1e6);
CXData.RateCoeff = log10(10.^CXRateCoeff./1e6);  % ✅ Normalize Charge Exchange Rates
CXData.ChargeState = CXChargeState;

% Open the file
% ncid = netcdf.create('ADAS_Rates_Ne.nc','NC_WRITE');
ncid = netcdf.create('ADAS_Rates_Ne.nc','NC_WRITE');

% Define the dimensions
dimScalar = netcdf.defDim(ncid,'scalar',1);
dimPair = netcdf.defDim(ncid,'pair',2);
dimTemp_Ionize = netcdf.defDim(ncid,'n_Temperatures_Ionize',length(IonizationData.Temp));
dimDensity_Ionize = netcdf.defDim(ncid,'n_Densities_Ionize',length(IonizationData.Density));
dimTemp_Recombine = netcdf.defDim(ncid,'n_Temperatures_Recombine',length(RecombinationData.Temp));
dimDensity_Recombine = netcdf.defDim(ncid,'n_Densities_Recombine',length(RecombinationData.Density));
dimTemp_CX = netcdf.defDim(ncid,'n_Temperatures_CX',length(CXData.Temp));
dimDensity_CX = netcdf.defDim(ncid,'n_Densities_CX',length(CXData.Density));

dimChargeState_Ionize = netcdf.defDim(ncid,'n_ChargeStates_Ionize',length(IonizationData.ChargeState));
dimChargeState_Recombine = netcdf.defDim(ncid,'n_ChargeStates_Recombine',length(RecombinationData.ChargeState)); 
dimChargeState_CX = netcdf.defDim(ncid,'n_ChargeStates_CX',length(CXData.ChargeState));  % ✅ Add CX charge states

% Define IDs for the dimension variables
Z_ID = netcdf.defVar(ncid,'Atomic_Number','int',[dimScalar]);
TempGridIonization = netcdf.defVar(ncid,'gridTemperature_Ionization','double',[dimTemp_Ionize]);
DensityGridIonization = netcdf.defVar(ncid,'gridDensity_Ionization','double',[dimDensity_Ionize]);
TempGridRecombination = netcdf.defVar(ncid,'gridTemperature_Recombination','double',[dimTemp_Recombine]);
DensityGridRecombination = netcdf.defVar(ncid,'gridDensity_Recombine','double',[dimDensity_Recombine]);
TempGridCX = netcdf.defVar(ncid,'gridTemperature_CX','double',[dimTemp_CX]);  % ✅ CX Grid
DensityGridCX = netcdf.defVar(ncid,'gridDensity_CX','double',[dimDensity_CX]);  % ✅ CX Density Grid

ChargeStateGridIonization = netcdf.defVar(ncid,'gridChargeState_Ionization','double',[dimChargeState_Ionize dimPair]);
ChargeStateGridRecombination = netcdf.defVar(ncid,'gridChargeState_Recombination','double',[dimChargeState_Recombine dimPair]);
ChargeStateGridCX = netcdf.defVar(ncid,'gridChargeState_CX','double',[dimChargeState_CX dimPair]);  % ✅ CX Charge States

% Define the main variable
IonizeCoeff = netcdf.defVar(ncid,'IonizationRateCoeff','double',[dimDensity_Ionize dimTemp_Ionize dimChargeState_Ionize]);
RecombineCoeff = netcdf.defVar(ncid,'RecombinationRateCoeff','double',[dimDensity_Recombine dimTemp_Recombine dimChargeState_Recombine]);
CXCoeff = netcdf.defVar(ncid,'ChargeExchangeRateCoeff','double',[dimDensity_CX dimTemp_CX dimChargeState_CX]);  % ✅ CX Coefficients

% We are done defining the NetCDF
netcdf.endDef(ncid);

% Store the dimension variables
netcdf.putVar(ncid,Z_ID,Z);
netcdf.putVar(ncid,TempGridIonization,IonizationData.Temp);
netcdf.putVar(ncid,DensityGridIonization,IonizationData.Density);
netcdf.putVar(ncid,TempGridRecombination,RecombinationData.Temp);
netcdf.putVar(ncid,DensityGridRecombination,RecombinationData.Density);
netcdf.putVar(ncid,TempGridCX,CXData.Temp);  % ✅ CX Temperature Grid
netcdf.putVar(ncid,DensityGridCX,CXData.Density);  % ✅ CX Density Grid

netcdf.putVar(ncid,ChargeStateGridIonization,IonizationData.ChargeState);
netcdf.putVar(ncid,ChargeStateGridRecombination,RecombinationData.ChargeState);
netcdf.putVar(ncid,ChargeStateGridCX,CXData.ChargeState);  % ✅ CX Charge States

% Store the main variables
netcdf.putVar(ncid,IonizeCoeff,IonizationData.RateCoeff);
netcdf.putVar(ncid,RecombineCoeff,RecombinationData.RateCoeff);
netcdf.putVar(ncid,CXCoeff,CXData.RateCoeff);  % ✅ Write CX Rate Coefficients

% Close the NetCDF
netcdf.close(ncid);

% Interpolation Example for Ionization Rate Coefficients
Coeff = interpn(IonizationData.Density,IonizationData.Temp,IonizationData.RateCoeff(:,:,1),18,log10(1:1:20),'linear',0);
te = IonizationData.Temp;
ne = IonizationData.Density;

% Plot Ionization Rate Coefficients
figure;
plot(te, IonizationRateCoeff(15,:,1));
title('Ionization Rate Coefficients vs Temperature');

figure;
plot(ne, IonizationRateCoeff(:,14,1));
title('Ionization Rate Coefficients vs Density');

% ✅ Plot Charge Exchange Rate Coefficients
figure;
plot(te, CXRateCoeff(15,:,1));
title('Charge Exchange Rate Coefficients vs Temperature');

figure;
plot(ne, CXRateCoeff(:,14,1));
title('Charge Exchange Rate Coefficients vs Density');