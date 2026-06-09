% In this script I attempt to extract heat fluxes from the front facing
% surface using IR image and the semi-infinite slan model

clear all
% close all

%%%%%%%%%%%%FIND FILEPATH WITH RAW DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Define path and filenames
% 2018_04_05
shot = 20958; % all probes removed  
PATHNAME = '\\mpexserver\ProtoMPEX_DATA\IR_Camera\2018_04_05\';
% 2018_03_23
shot = 20399; % all probes removed  
% shot = 20400; 
PATHNAME = '\\mpexserver\ProtoMPEX_DATA\IR_Camera\2018_03_23\';
%%%%%%%%%%%%%%%%%%%READ RAW INTENSITY DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[RawData,frames,seq] = ExtractRawData(PATHNAME,FILENAME,shot);
[RawData,frames,seq] = ExtractRawData(PATHNAME,strjoin({'Shot', strcat(num2str(shot), '.seq')}),shot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%Get Temperature From Emissivity%%%%%%%%%%%%%%%%%%%%%%%%%%
emissivity = 1*0.69;
[Temp,dT] = IntensityTempConv2(emissivity,0.2*emissivity,RawData,seq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%Create Time Vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framerate = 100;
[time,ind_s,ind_e] = CreateTimeVector2(RawData,shot,framerate,frames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Create Geometric Coordinates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pxx_per_cm = 24.14;
pxy_per_cm = 24.14;
rad_peak = 2.5;
[xx_c,xx_s,xx_e,yy_c,yy_s,yy_e] = FindCenter(Temp,floor((ind_s+ind_e)/2),pxx_per_cm, pxy_per_cm, rad_peak);


%%%%%%%%Create x and y vectors%%%%%%%%%%%%%
dim = size(Temp);
[xx_vec,yy_vec] = Pix2XY(xx_c,yy_c,pxx_per_cm,dim);
T0 = Temp(:,:,ind_s);

for s = 1:size(Temp,3)
    T(:,:,s) = Temp(:,:,s) - T0;
end

% [TempSpace,Inten_vec,eps_vec] = TempEmissivityIntensitySpace(0.69,0.2*0.69,RawData,seq,ind_s,end_s);



%%
if 0
    figure;
    for s = ind_s:ind_e
        surf(Temp(:,:,s)- T0,'LineStyle','none')
        view([30,60])
        caxis([0,160])
        zlim([0,160])
        ylim([100,250])
        xlim([250,450])
        drawnow
    end
end

%%
for s = ind_s:size(Temp,3)
     Ts1(:,:,s) = (Temp(:,:,s) - T0);
     Ts2(:,:,s) = Ts1(:,:,s).^2;
end

dt = 1/framerate;

for s = 1:size(Ts2,3)-1
    dTs1_dt(:,:,s) = (Ts1(:,:,s+1) - Ts1(:,:,s))/dt;
    dTs2_dt(:,:,s) = (Ts2(:,:,s+1) - Ts2(:,:,s))/dt;
end
%%
if 0 
% dT_dt
    figure;
    for s = ind_s:ind_e
        surf(dTs1_dt(:,:,s),'LineStyle','none')
        view([30,60])
        caxis([0,4e2])
        zlim([0,4e2])
        ylim([100,250])
        xlim([250,450])
        pause(0.001);
    end
end
%%
mdsconnect('mpexserver');
[heat_flux, heat_peak, heat_peak_x, heat_peak_y, heat_power ...
    heat_radius, T, X, Y, comments] = get_IR(shot);

%%
rho = 8000; 
K = 15;
cp = 500;
b = sqrt(rho*K*cp);
if 1
    t_q0(1) = 4.04;
    for s = 1:size(Ts2,3)-1
        q0(:,:,s) = 0.5*sqrt(pi)*b*sqrt(dTs2_dt(:,:,s));
        q0center(s) = max(max(max(q0(:,:,s))));
        t_q0(s+1) = t_q0(s) + dt;
    end
end
%%
% close all
figure; hold on 
h(1) = plot(T,heat_peak*1e-6,'k','LineWidth',2)
h(2) = plot(time(1:end-1)-0.01,real(q0center)*1e-6,'r','LineWidth',2)

xlim([4,5])
ylim([0,2])
legend(h,'COMSOL','InfiniteSlabModel')
box on
set(gcf,'color','w')
set(gca,'PlotBoxAspectRatio',[1 0.4 1])
ylabel('Heat flux [MWm^-2]')
xlabel('t [s]')
title(['shot: ',num2str(shot)])

return

figure; hold on
box on
set(gcf,'color','w')
set(gca,'PlotBoxAspectRatio',[1 0.4 1])
ylabel('Heat flux [MWm^-2]')
title(['shot: ',num2str(shot)])

t_sample = linspace(4.18,4.56);
for s = 1:length(t_sample)
    hold on
t_rng = find(T>t_sample(s),1);
[M,I] = max(max(heat_flux(:,:,t_rng),[],2));
plot(X,heat_flux(I,:,t_rng)*1e-6,'k','LineWidth',2)


t_rng = find(time>t_sample(s)-0.01,1);
[M,I] = max(max(q0(:,:,t_rng),[],2));
plot(yy_vec,real(q0(I,:,t_rng))*1e-6,'r','LineWidth',2)
ylim([0,1.5])
xlim([-4,4]/100)
title(['time: ',num2str(t_sample(s))])
pause(0.01)
if s ~= (length(t_sample))
cla
end
end