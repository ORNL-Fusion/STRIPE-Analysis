function [Term] = EmissivityErrorTerm(eps_mean,eps_delta,RawData,seq,ind_s)

% Create TempSpace given the range of Intensities and emissivities
[TempSpace,Inten_vec,eps_vec] = TempEmissivityIntensitySpace(eps_mean,eps_delta,RawData,seq);

%


end

