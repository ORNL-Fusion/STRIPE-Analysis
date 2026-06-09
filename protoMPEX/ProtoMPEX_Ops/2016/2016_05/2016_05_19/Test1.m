clear all

% =========================================================================
% Connect to Server: 
% myMDSconnect(1)
mdsconnect('mpexserver') 

% Open shot:
[SHOT,~]=mdsopen('MPEX',8882);

% Data location:
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RA = [Stem,Branch];

Data{1} = [RA,'LP_V_RAMP'];
Data{2} = [RA,'TARGET_LP'];
if 0
    Data{3} = [RA,'PWR_28GHZ']; % on May 19th 2016, this channel was not digitized
else
    Data{3} = [Stem,'FSCOPE:','TUBE07:','PMT_VOLT'];
end
Data{4} = [Stem,'SHOT_NOTE'];

% Extract data:
[Vm,~]  = mdsvalue(Data{1});
[Im,~]  = mdsvalue(Data{2});
[ECH,~]  = mdsvalue(Data{3});
[NOTE,~]= mdsvalue(Data{4});

[t_1,~]  = mdsvalue(['DIM_OF(',Data{1},')']);
[t_2,~]  = mdsvalue(['DIM_OF(',Data{2},')']);
[t_ECH,~] = mdsvalue(['DIM_OF(',Data{3},')']);

% Create a function where we only need to give the name of the Stem, Branch
% and data name and it returns both the x and y arrays of the data.
% In addition, we could create a function that evaluates all the data or at
% least all the probes and tells you which datas are populated
 
mdsclose;
mdsdisconnect;
% =========================================================================
% DLP data analysis: 
L_tip = 1.2/1000;  % [m]
D_tip = 0.25/1000; % [m]
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
Vm = 2*11*Vm(rng); 
if 1
    Im = sgolay_t(Im,3,7);
end
Im = -2*Im(rng)/150;
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
%% Plot data:
close all

figure; hold on
plot(t,Vm)
plot(t,Im*1000)

figure; hold on
subplot(3,1,1)
plot(time,Ni,'k-')
title('n_i')
ylim([0,8e19])
ax(1) = gca;

subplot(3,1,2)
plot(time(1:end-3),sgolay_t(Te(1:end-3),3,7),'k-')
%plot(time(1:end-3),Te(1:end-3),'k-')
title('Te')
ylim([0,20])
xlim([tStart,tEnd])
ax(2) = gca;

subplot(3,1,3)
plot(t_ECH,sgolay_t(ECH,3,5),'k-')
title('28 GHz')
%ylim([0,0.1])
ax(3) = gca;
set(ax,'Xlim',[tStart,tEnd])

set(gcf,'color','w')

figure; 
n = 36;
for s = 1:16
subplot(4,4,s); hold on
plot(V{s+n},I{s+n},'k')
plot(V{s+n},Ifit{s+n}(x1(:,s+n)),'r')
title(num2str(time(s+n)))
end
