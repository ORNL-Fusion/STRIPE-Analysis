clear all
close all

% =========================================================================
% Connect to Server: 
% myMDSconnect(1)
mdsconnect('mpexserver') 

% Extract data:
[SHOT,~]=mdsopen('MPEX',8888); % see 8636

Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
%Branch = 'MPEX1:';
RA = [Stem,Branch];

% Data{1} = [RA,'RF_REF_PWR'];
Data{1} = [RA,'RF_REF_PWR'];
Data{2} = [RA,'TARGET_LP'];
Data{3} = [Stem,'FSCOPE:','TUBE04'];
Data{4} = [Stem,'SHOT_NOTE'];

[Vm,~]  = mdsvalue(Data{1});
[Im,~]  = mdsvalue(Data{2});
[Pg,~]  = mdsvalue(Data{3});
[NOTE,~]= mdsvalue(Data{4});
 
[t_1,~]  = mdsvalue(['DIM_OF(',Data{1},')']);
[t_2,~]  = mdsvalue(['DIM_OF(',Data{2},')']);
[t_Pg,~] = mdsvalue(['DIM_OF(',Data{3},')']);

% Create a function where we only need to give the name of the Stem, Branch
% and data name and it returns both the x and y arrays of the data.
% In addition, we could create a function that evaluates all the data or at
% least all the probes and tells you which datas are populated
 
mdsclose;
mdsdisconnect;
% =========================================================================
% DLP data analysis: 
L_tip = 1.2/1000;  % [m]
D_tip = 2*1.127/1000; % [m]
% D_tip = 0.25/1000; % [m]
rp = D_tip/2; % [m]

AreaType = 2;
switch AreaType
    case 1
        Area = L_tip*pi*D_tip + 0.25*pi*D_tip^2;
    case 2
        Area = 2*L_tip*D_tip;
end

e_c = 1.602e-19; 
m_p = 1.6726e-27;
e_0 = 8.854187817*1e-12; % F/m
AMU = 2;

% Break data into cycles: -------------------------------------------------

tStart = 4.17; % [s]
tEnd = 4.32;
rng = find(t_1>=tStart & t_1<=tEnd);
t = t_1(rng);

% Apply correct scale factors to data:
Vm = 1*10*Vm(rng); 
if 1
    Im = sgolay_t(Im,3,69);
end
Im = -10*Im(rng)/150;
tm = 2;

[locs,~] = peakseek(abs(Vm),211);
for s = 1:(length(locs)-1)
    I{s} = Im(locs(s):locs(s+1));
    V{s} = Vm(locs(s):locs(s+1));

    switch tm
        case 1
            time(s) = ( t(locs(s)) + t(locs(s+1)) )/2;
        case 2
            time(s) = t(locs(s));
    end
    
end

x0 = [30e-3,2,0,0,0];

% Fitting routine: --------------------------------------------------------
for s = 1:length(V)
    A1 = V{s};
    A2 = I{s};
    res = @(x)   ( x(1)*tanh((A1-x(4))/x(2)) + (x(3)*(A1-x(4))) + x(5) ) - A2;
    Ifit{s} = @(x)   x(1)*tanh((A1-x(4))/x(2)) + (x(3)*(A1-x(4))) + x(5); 
    
    [x1(:,s),ssq,cnt] = LMFnlsq(res,x0);
end

Te = real(x1(2,:));
Isat = real(x1(1,:));
Cs = sqrt(e_c*Te/(AMU*m_p));
Ni = real(Isat./(0.61*e_c*Area*Cs));
Ld = real(sqrt(e_0*Te./(Ni*e_c)));% Debye length
Xi = rp./Ld;

% =========================================================================
% Plot data:

figure; hold on
plot(t,Vm)
plot(t,Im*1000)

figure; hold on
subplot(3,1,1)
plot(time,Ni,'ko-')
title('N_i')
ylim([0,1e19])

subplot(3,1,2)
plot(time,Te,'ko-')
title('Te')
ylim([0,10])
xlim([tStart,tEnd])

subplot(3,1,3)
plot(time,Xi,'ko-')
title('Xi = r_p/{\lambda_d}')
%set(gcf,'position',[467   188   572   587])
ylim([0,100])

figure; 
for s = 1:16
subplot(4,4,s); hold on
plot(V{s},I{s},'k')
plot(V{s},Ifit{s}(x1(:,s)),'r')
end
% title(['n_e = ',num2str(Ni),'m^{-3} ','  ,T_e = ',num2str(Te),' eV'])


% Pulses today
% 8634 full shot
% 8635 FS
% 8636 FS

% 
% [Pwr,~]=mdsvalue('\MPEX::TOP.MACHOPS1:GEN_RF_PWR')
