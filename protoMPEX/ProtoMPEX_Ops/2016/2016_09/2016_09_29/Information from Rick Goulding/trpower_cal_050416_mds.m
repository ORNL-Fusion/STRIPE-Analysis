function [t,p,volts]=trpower_cal_050416_mds(shotlist)
% R. H. Goulding 5/25/11  Function to calculate |B| at midplane of antenna
% on axis for
% helicon device with magnet locations corresponding to
% "symmetric_large_spacing_feb10_layout.mat" and
% "helicon_ldrd_full_assembly_asm_9_13_10.ar"

%b = magnetic field (tesla)
%cur1 = (equal) currents in outer magnet pair
%cur2 = (equal) currents in inner magnet pair
coeffs=[-6.146e-6,.0017362,.012256,-.0014283,.00080496];
%coeffs=[-6.146e-6,.0017362,.012256,-.0014283,.00080496]; %8471A detector CPT.004
%coeffs=[4.9e-5,7.96e-6,.017444,.010148,.0038132];



%fit to pwr detector characteristic
%attendb=38. ; %pwr attenuation in dB including all attenuation between coax and detector
coupler_atten=67.3; %dB
cable_atten=1.73/2; %dB
splitter_atten=3.6; %dB
atten_correct=-0.66;  % 5/4/16 Additional attenuation to make result agree with scope traces
attendb=coupler_atten+cable_atten+splitter_atten+atten_correct;
pchan='GEN_RF_PWR';
fchan='TRANS_I';
rchan='RF_REF_PWR';
ichan='MN_CURRENT';
[t,array]=read_mdso(pchan,shotlist);
volts=squeeze(array(:,1,:));


sv=size(volts);


p=zeros(sv(1),sv(2))+coeffs(1);



for iord=2:5, p=p+coeffs(iord)*volts.^(iord-1);end;

p=p*10^(attendb/10);


