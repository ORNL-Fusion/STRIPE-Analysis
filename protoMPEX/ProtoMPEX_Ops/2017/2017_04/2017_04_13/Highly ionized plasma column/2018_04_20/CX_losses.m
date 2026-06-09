% CX collision frequency

clear all 
close all

Ti = 5; % V
Tg = 300; % K
P0 = linspace(0.01,5); % mTorr
n0 = P0*0.1333/(k_B*Tg);
ne = 5e19;
L = 1;
V  = L*pi*(2.5/100)^2;

h = 3;
switch h
    case 1
    cs  = 2e-19; % based on mom transfer H+ + H2(0)
    vti = sqrt(2*e_c*Ti/(2*m_p));
    case 2       % CX cross section at ion at 5 eV and brennan's solution to the integral
    cs  = 1e-20;
    vti = sqrt(10*2*e_c*Ti/(2*m_p));
    case 3 % CX cross section and assume neutrals are stationary relative to ions
    cs  = 1e-20;
    vti = sqrt(10*2*e_c*Ti/(2*m_p));
end

P_ICH = 12e3;
nu_CX = n0*cs*vti;

P_CX = (5/2)*e_c*Ti*ne*nu_CX*V;

P_ion = P_ICH - P_CX; 

eta_ion = P_ion/P_ICH;

figure 
plot(P0,P_CX/P_ICH)
ylim([0,1])
xlabel('P_0 [mTorr]')
ylabel('CX loss %')


figure 
semilogy(P0,vti./nu_CX)
ylim([0,2])
ylabel('L_{cx} [m]')
xlabel('P_0 [mTorr]')

%% surface plot Ti and n0
Ti = linspace(1,20); % V
Tg = 300; % K
P0 = linspace(0.01,5); % mTorr
n0 = P0*0.1333/(k_B*Tg);
ne = 5e19;
L = 1;
V  = L*pi*(2.5/100)^2;

h = 3;
switch h
    case 1
    cs  = 2e-19; % based on mom transfer H+ + H2(0)
    vti = sqrt(2*e_c*Ti/(2*m_p));
    case 2       % CX cross section at ion at 5 eV and brennan's solution to the integral
    cs  = 1e-20;
    vti = sqrt(10*2*e_c*Ti/(2*m_p));
    case 3 % CX cross section and assume neutrals are stationary relative to ions
    cs  = 1e-20;
    vti = sqrt(10*2*e_c*Ti/(2*m_p));
end

for k = 1:length(Ti)
    for p = 1:length(n0)
    nu_CX(k,p) = n0(p).*cs.*vti(k);
    P_CX(k,p) = (5/2)*e_c*ne*nu_CX(k,p).*Ti(k)*V;
    end
end

figure; [C,h] = contourf(Ti,P0,P_CX*1e-3,10); clabel(C,h)