% In this script I attempt to extract heat fluxes from the front facing
% surface using IR image and the semi-infinite slab model

clear all
close all

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

% Material properties:
rho = 8000; 
K = 15;
cp = 500;
b = sqrt(rho*K*cp);

% Error terms 
dK = 0.1*K;
drho = 0.1*rho;
dcp = 0.1*cp;

Term_dK = (dK/(2*K))^2;
Term_drho = (drho/(2*rho))^2;
Term_dcp = (dcp/(2*cp))^2;
Term_dT = (dT./(2*T)).^2;

dq_q = sqrt( Term_dK + Term_drho + Term_dcp + Term_dT);

for i = 1:size(T,1)
    for j = 1:size(T,2)
        for tt = 1:size(time,2)
         qSlab(i,j,tt) = sqrt(pi/4)*b*T(i,j,tt)./sqrt(time(tt)-3.84);
        end
    end
end

for tt = 1:size(time,2)
        qSlab_0(tt) = max(max(qSlab(i,j,tt)));
end

%%
for s = ind_s:size(Temp,3)
     T2(:,:,s) = T(:,:,s).^2;
end

dt = 1/framerate;

for s = 1:size(T2,3)-1
    dTs1_dt(:,:,s) = (T(:,:,s+1) - T(:,:,s))/dt;
    dTs2_dt(:,:,s) = (T2(:,:,s+1) - T2(:,:,s))/dt;
end

%%
mdsconnect('mpexserver');
[heat_flux, heat_peak, heat_peak_x, heat_peak_y, heat_power ...
    heat_radius, T, X, Y, comments] = get_IR(shot);

%%
if 1
    t_q(1) = 4.04;
    for s = 1:size(T2,3)-1
        q(:,:,s) = 0.5*sqrt(pi)*b*sqrt(dTs2_dt(:,:,s));
        dq(:,:,s) = q(:,:,s).*dq_q(:,:,s);
        
        q_0(s) = max(max(max(q(:,:,s))));
       
        t_q(s+1) = t_q(s) + dt;
    end
end
 
%%
 close all
figure; hold on 
h(1) = plot(T,heat_peak*1e-6,'k','LineWidth',2)
h(2) = plot(time(1:end-1)-0.01,real(q_0)*1e-6,'r','LineWidth',2)
h(3) = plot(time(1:end-0)-0.01,real(qSlab_0)*1e-6,'g','LineWidth',2)

xlim([4,5])
ylim([0,2])
legend(h,'COMSOL','InfiniteSlabModel')
box on
set(gcf,'color','w')
set(gca,'PlotBoxAspectRatio',[1 0.4 1])
ylabel('Heat flux [MWm^-2]')


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
h(1) = plot(X,heat_flux(I,:,t_rng)*1e-6,'k','LineWidth',2)


t_rng = find(time>t_sample(s)-0.01,1);
[M,I] = max(max(q(:,:,t_rng),[],2));
% plot(yy_vec,real(q(I,:,t_rng))*1e-6,'r','LineWidth',2)
h(2) = errorbar(yy_vec,real(q(I,:,t_rng))*1e-6,real(dq(I,:,t_rng))*1e-6,'r','LineWidth',0.5)
% legend(h,'COMSOL','InfiniteSlabModel')

ylim([0,1.8])
xlim([-4,4]/100)
title(['time: ',num2str(t_sample(s))])
pause(0.01)
if s ~= (length(t_sample))
cla

xlabel('t [s]')
title(['shot: ',num2str(shot)])

end
end