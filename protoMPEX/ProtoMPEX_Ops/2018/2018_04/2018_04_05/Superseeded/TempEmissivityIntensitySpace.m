function [T1space,Inten_vec,eps_vec] = TempEmissivityIntensitySpace(eps_mean,eps_delta,RawData,seq,ind_s)
% Given the raw data and information on the emissitivy uncertainty,
% calculate the Temp for the entire space of possible emissivity and
% intensity data

% Range of Intensities
Inten_min = min(min(min(RawData)));
Inten_max = max(max(max(RawData)));
% Range of emissivities
eps_min = eps_mean - eps_delta;
eps_max = eps_mean + eps_delta;

% Create space
Ne = 1e2;
eps_vec = linspace(eps_min,eps_max,Ne);
Ni = 2e2;
Inten_vec = linspace(Inten_min,Inten_max,Ni);

% Calculate T1 space
for neps = 1:Ne
    for nInten = 1:Ni
        f(nInten,neps) = seq.ThermalImage.GetValueFromEmissivity(eps_vec(neps),Inten_vec(nInten));
    end
end
T1space = f;

% Calculat T0 space
n = find(eps_vec >= eps_mean,1);
ind_s = 34;
T0space = interp1(Inten_vec,T1space(:,n),RawData(:,:,ind_s));


end