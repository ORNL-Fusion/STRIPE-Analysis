
hold on;

clearvars -except Data
clear all
close all
clc
% load('Mask.mat')

OxyPercent=0.0;


%% Extended voltage range with updated yields not using angles
erosion_rate=1.6073e+17; %20KW 


%%
nP=1e6;
erosionPP=erosion_rate/nP; %atoms/s/particle

% M = csvread('iterGeom.csv');
% r = M(:,1);
% z = M(:,2);
% rGeom = r;
% zGeom = z;
total_redeposition_fraction = [];
local_redep_fraction = [];
prompt_redep_fraction = [];
net_erosion_fraction = [];
angles = 0:10:90;
% for angle=0:10:90
    
plot_tracks = 0;
file = strcat(pwd,'/positions.nc');
hitWall = ncread(file,'hitWall');
nHit = length(find(hitWall));
hasHit = find(hitWall);
notHit = find(hitWall==0);
x0 = ncread(file,'x');
y0 = ncread(file,'y');
z0 = ncread(file,'z');
vx0 = ncread(file,'vx');
vy0 = ncread(file,'vy');
vz0 = ncread(file,'vz');
time0=ncread(file,'time');
distTraveled = ncread(file,'distTraveled');
charge0 = ncread(file,'charge');
weight0 = ncread(file,'weight');
vtot = sqrt(vx0.^2 +vy0.^2 + vz0.^2);
E = 0.5*27*1.66e-27*vtot.^2/1.602e-19;
figure(11)
histogram(E)

figure()
scatter3(x0(hasHit),y0(hasHit),z0(hasHit))
hold on
scatter3(x0(notHit),y0(notHit),z0(notHit))
xlim([-.1 .1])
ylim([-.1 .1])
zlim([0 5])

r0=sqrt(x0.^2+y0.^2);

qR=find(r0>0.0625);

figure; scatter3(x0(qR),y0(qR),z0(qR))

local_length = 8e-4;
prompt_length = 31e-3;
redep = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1));
redep1 = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & (charge0 == 1));
redep2 = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & (charge0 == 2));
redep_local = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & distTraveled < local_length);
redep_prompt = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & distTraveled > local_length & distTraveled < prompt_length);

if plot_tracks
file = strcat(pwd,'/history.nc');
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
sizeArray = size(x);
nP = sizeArray(2);
hit = find(weight(end,:) < 1);

r = sqrt(x.^2 + y.^2);


 figure(1)
 hold on
%  scatter3(x0(redep1),y0(redep1),z0(redep1))
% % i=581
% % plot(r(:,i),z(:,i))
% %41242
for i=1:1:length(x0)
 
% plot(r(:,i),z(:,i))
plot3(x(:,i),y(:,i),z(:,i),'k')
end
end
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% % plot(rGeom,zGeom)
% % hold on
% % scatter(rGeom,zGeom)
% axis equal
hold on
if (exist('x1') == 0)
    fid = fopen(strcat(pwd,'/gitrGeometryPointPlane3d.cfg'));
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:21
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end

%%
surface = find(Zsurface);
nSurfaces = length(a);
%Section for finding a subset of planes to plot
r = sqrt(x1.^2 + y1.^2);
subset = 1:length(x1);%find(r<0.07 & z1> 0.001 & z1 < .20);
% subset = find(r>0.049 & z1 > 1.5 & z1<2);
figure(1)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
% patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface)
title('Sample Particle Tracks for Helicon Erosion Source')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% axis([xMin xMax yMin yMax zMin zMax])
% cb1 = colorbar
set(gca,'fontsize',18)
set(gcf,'color','w')
% set(cb1,'fontsize',14);
%axis equal
%%

file = strcat(pwd,'/surface.nc');
grossDep0 = ncread(file,'grossDeposition');
grossEro0 = ncread(file,'grossErosion');

subset = find(Zsurface);
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
% THETA = atan2d(Y,X);
% ss=sum(sign(THETA),2);
% swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
% THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );

%%
THETA = 0*Z;
for ii=1:length(THETA)
    %     if ii==2457
%     if ii==4859
%         ii
%     end
    THETA(ii,:) = atan2d(Y(ii,:),X(ii,:));
    diff21 = abs(THETA(ii,2) - THETA(ii,1));
    diff31 = abs(THETA(ii,3) - THETA(ii,1));
    diff32 = abs(THETA(ii,3) - THETA(ii,2));
    
    maxdiff = max([diff21,diff31,diff32]);
    if maxdiff > 90
        if THETA(ii,1)>0
            % THETA(ii,1)=THETA(ii,1)+360;
            %     end
            if THETA(ii,2)<0
                THETA(ii,2)=THETA(ii,2)+360;
            end
            if THETA(ii,3)<0
                THETA(ii,3)=THETA(ii,3)+360;
            end
        end
        if THETA(ii,1)<0
            % THETA(ii,1)=THETA(ii,1)+360;
            %     end
            if THETA(ii,2)>0
                THETA(ii,2)=THETA(ii,2)-360;
            end
            if THETA(ii,3)>0
                THETA(ii,3)=THETA(ii,3)-360;
            end
        end
    end
end
THETA=THETA;

for ii=1:length(THETA)
    if THETA(ii,1)  < 0
        THETA(ii,1)=THETA(ii,1)+ 360;
        THETA(ii,2)=THETA(ii,2)+360;
        THETA(ii,3)=THETA(ii,3)+360;
    end
end

%
ElementArea=area(surface)';
pointsize=10;
% Mask(:,1)=Mask(:,1)+1.7462;

GrossErosion=(erosionPP.*grossEro0./ElementArea);
%GrossErosion=GrossErosion*OxyPercent+((1-OxyPercent)*Data.GrossErosion);

figure;
patch(transpose(Z),(transpose(THETA)),0*transpose(Z),GrossErosion,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max(GrossErosion)*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
% patch(transpose(X),transpose(Y),transpose(Z),grossEro0,'FaceAlpha',.3,'EdgeAlpha', 0.3)
c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',15)
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
set(gca, 'FontSize',18)
ylim([0 360])
xlim([0.97 1.3456])
caxis([0 0.5E19])
set(gcf,'color','w')
title('Gross Erosion','FontSize',18);
hold off
% set(gca,'ColorScale','log')
 % caxis([1e22 1e23])

%view(90,90)
%%
% ncdisp(file)
% edist = ncread(file,'surfEDist');
% 
% nE=100;
% eDtotal = reshape(sum(edist,3),[90,nE]);
% 
% figure(105)
% p4 = imagesc(linspace(0,90,90),linspace(0,1000,nE),eDtotal');
% colorbar
% meanEn=mean(eDtotal,1);
% figure; plot(linspace(0,1000,nE),movmean(meanEn,10)');
% 
% 
% % caxis([0 1000])
% set(gca, 'YDir', 'normal');
% 
% numSurfaces=1484
% 
% meanE = zeros(1,length(numSurfaces));
% 
% energy_iead_gitr = linspace(5,995,nE);
% 
% for i=1:numSurfaces
% 
%     iead = reshape(edist(:,:,i),[90,nE]);
% 
%     edist_i = sum(iead);
% 
%     meanE(i)= sum(energy_iead_gitr.*edist_i)/sum(edist_i);
% 
% end
% 
% figure
% 
%  p4 = imagesc(linspace(0,90,90),linspace(0,10540,nE),iead');
%  set(gca,'YDir','normal')
% 
% colorbar
% % caxis([0 0.2])




%%

GrossDeposition=erosionPP.*grossDep0./ElementArea;
%GrossDeposition=GrossDeposition*OxyPercent+((1-OxyPercent)*Data.GrossErosion);

figure(3)
patch(transpose(Z),(transpose(THETA)),0*transpose(Z),GrossDeposition,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max(GrossDeposition)*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

% patch(transpose(X),transpose(Y),transpose(Z),grossDep0,'FaceAlpha',.3,'EdgeAlpha', 0.3)
c=colorbar;
ylabel(c, 'Flux [Al ion Particles/m^2/s]','FontSize',15)
title('Net Erosion','FontSize',18);
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
ylim([0 360])
xlim([0.97 1.3456])
set(gca, 'FontSize',18)
caxis([0 0.5E19])
set(gcf,'color','w')
title('Gross Deposition','FontSize',18);
hold off

% THETA = atan2d(Y,X);
% ss=sum(sign(THETA),2);
% swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
% THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );
% 
% for ii=1:length(THETA)
%     if THETA(ii,1)<0
% THETA(ii,1)=THETA(ii,1)+360;
%     end
%     if THETA(ii,2)<0
% THETA(ii,2)=THETA(ii,2)+360;
%     end
%     if THETA(ii,3)<0
% THETA(ii,3)=THETA(ii,3)+360;
%     end
% end

NetErosion=erosionPP.*(grossEro0-grossDep0)./ElementArea;
%NetErosion=NetErosion*OxyPercent+((1-OxyPercent)*Data.NetErosion);

figure(4)
patch(transpose(Z),(transpose(THETA)),0*transpose(Z),NetErosion,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max(NetErosion)*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',15)
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
ylim([0 360])
xlim([0.97 1.3456])
caxis([0 0.6E19])
set(gca, 'FontSize',18)
set(gcf,'color','w')
title('Net Erosion of Al','FontSize',18);
hold off

netEro = (grossEro0-grossDep0);
netEro(find(netEro >= 0)) = nan;

% figure(5)
% patch(transpose(Z),transpose(THETA),0*transpose(Z),erosionPP.*netEro./ElementArea,'FaceAlpha',1,'EdgeAlpha', 0.3)
% hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netEro./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
% 
% c=colorbar;
% ylabel(c, 'Particles/m^2/s','FontSize',18);
% ylabel('Angle [Deg]','FontSize',18);
% xlabel('Z [m]','FontSize',18);
% set(gcf,'color','w')
% caxis([0 inf])
% title({'Net Erosion of Al','Deposition Regions Only'},'FontSize',18);
% hold off

netDep=abs(netEro);

NetDep=erosionPP.*netDep./ElementArea;
%NetDep=NetDep*OxyPercent+((1-OxyPercent)*Data.NetDeposition);

for ii=1:length(NetDep)
   if isnan(NetDep(ii))==1
       NetDep(ii)=0;
   end
end
Z_cyl = Z;
%%
figure(6)
patch(transpose(Z),(transpose(THETA)),0*transpose(Z),NetDep,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',15)
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
ylim([0 360])
xlim([0.97 1.3456])
caxis([0 0.1E19])
set(gca, 'FontSize',18)
set(gcf,'color','w')
title({'Net Deposition of Al'},'FontSize',18);
hold off

%% End cap physics
nExited = sum(weight0(find(z0>= 0.64 & z0<=0.83)));

ExitedParticles=nExited*erosionPP;
PercentExited=nExited/nP*100;

%%% histogram of particles which hit end cap
nBins = 20;
% edges = linspace(-0.07,0.07,nBins+1);
% centers = linspace(-0.07,0.07,nBins);

edges = linspace(0,0.08,nBins+1);
centers = 0.5*(edges(1:end-1) + edges(2:end))
hist_bins = zeros(1,nBins);
areas = pi*(edges(2:end).^2 - edges(1:end-1).^2);

for i=1:nBins

    hist_bins(i) = sum(weight0(find(z0>= 0.65 & z0<=0.9 & sqrt(x0.^2+y0.^2) > edges(i) & sqrt(x0.^2+y0.^2) <= edges(i+1))));

end
radialprofile=hist_bins*erosionPP;
flux=radialprofile./areas;

figure(321)
plot(centers,movmean(flux,10))
xlabel('X [cm]','FontSize',18);
ylabel('Particles/s','FontSize',18);
set(gcf,'color','w')
title({'Paritlces leaving the helicon'},'FontSize',18);

%%

%For total 
% hitCap1 = find(z0<=-0.199999);
% hitCap2 = find(z0>=0.1999);

%For nuetrals 
% hitCap1 = find(z0<=-0.199999 & charge0==0);
% hitCap2 = find(z0>=0.1999 & charge0==0);

%For ionized
hitCap1 = find(z0>= 0.65 & z0<=0.9 & charge0>0);
hitCap2 = find(z0>=3.0899& z0 <=3.0901 & charge0>0);
% hitCap1 = find(z0<=0.3801 & charge0>0);
% hitCap2 = find(z0>=3.0999& charge0>0);
%  hitCap1 = find(z0<=-0.1499 & charge0>0);
%  hitCap2 = find(z0>=0.1499 & charge0>0);
% hitCap1 = find(z0<=-0.2999 & charge0>0);
% hitCap2 = find(z0>=0.2999 & charge0>0);


figure(322)
s1=scatter3(x0(hitCap1),y0(hitCap1),z0(hitCap1));
hold on
s2=scatter3(x0(hitCap2),y0(hitCap2),z0(hitCap2));
xlabel('X [cm]','FontSize',18);
ylabel('Y [cm]','FontSize',18);
zlabel('Z [cm]','FontSize',18);
set(gcf,'color','w')
title({'Paritlces leaving the helicon'},'FontSize',18);
hold off

% end cap flux 

%Upstream based on geometry
cap1 = find(Zsurface ==0 & z1 < 0.8);

cap1_bins = zeros(1,length(cap1));
for i=1:length(cap1)
  in = inpolygon(x0(hitCap1),y0(hitCap1),[x1(cap1(i)),x2(cap1(i)),x3(cap1(i))],[y1(cap1(i)),y2(cap1(i)),y3(cap1(i))]);
  cap1_bins(i) = sum(weight0(hitCap1(in)));
end
subset = cap1;%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(333)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),erosionPP*cap1_bins./area(cap1),'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Upstream End Cap Flux ','FontSize',18);
xlabel('X [m]','FontSize',18);
ylabel('Y [m]','FontSize',18);
zlabel('Z [m]','FontSize',18);
set(gcf,'color','w')
% axis([xMin xMax yMin yMax zMin zMax])
% cb1 = colorbar
set(gca,'fontsize',18)
% set(cb1,'fontsize',14);
caxis([0 0.8E17])
axis equal
c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',18)



%%
%Downstream based on geometry
cap2 = find(Zsurface ==0 & z1 > 3.05);
cap2_bins = zeros(1,length(cap2));
for i=1:length(cap2)
  in = inpolygon(x0(hitCap2),y0(hitCap2),[x1(cap2(i)),x2(cap2(i)),x3(cap2(i))],[y1(cap2(i)),y2(cap2(i)),y3(cap2(i))]);
  cap2_bins(i) = sum(weight0(hitCap2(in)));
end
subset = cap2;%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)

X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];

%%
ii=0.001;
jj=1;
while ii<=0.0603
    %ii=0.012
UpstreamRadius=ii;
Cap2radius=sqrt(x0(hitCap2).^2+y0(hitCap2).^2);
Cap2radius=find(Cap2radius<=UpstreamRadius );
Sumhitcap2=sum(weight0(hitCap2(Cap2radius)));
DownstreamAverageFlux(jj)=erosionPP*Sumhitcap2./(pi*UpstreamRadius^2); %0.020m is radius that maps downstream to target
radiusx(jj)=ii;

ii=ii+0.001;
jj=jj+1;
end

DownstreamAverageFlux(12) % Flux over different areas
%TargetImpurityFlux=erosionPP*Sumhitcap2./(pi*0.012^2)
DownstreamFlux=erosionPP*cap2_bins./area(cap2);


xflux=linspace(0,0.07,200); %35 is 0.012m
DownstreamRadialFlux=zeros(1,length(xflux)-1);

for ii=1:(length(xflux)-1)
    
Cap2radius=sqrt(x0(hitCap2).^2+y0(hitCap2).^2);
Cap2radius=find(Cap2radius>=xflux(ii) & Cap2radius<=xflux(ii+1));
Sumhitcap2=sum(weight0(hitCap2(Cap2radius)));
DownstreamRadialFlux(ii)=erosionPP*Sumhitcap2./((pi*xflux(ii+1)^2)-(pi*xflux(ii)^2)); %0.020m is radius that maps downstream to target

end

figure
plot(xflux(1:end-1),movmean(DownstreamRadialFlux,10))
title({' Al flux at z=4.14 m for the', 'High Density Magnetized Case'})
xlabel('x [m]')
ylabel('Al ion flux [Particles/m^2/s]')
set(gcf,'color','w')
set(gca,'fontsize',18)


figure
plot(radiusx,DownstreamAverageFlux,'k')

title({' Al flux at z=3.09 m'})
xlabel('Plasma Radius [m]')
ylabel('Al ion flux [Particles/m^2/s]')
set(gcf,'color','w')
set(gca,'fontsize',18)


%DownstreamFluxTot=sum(DownstreamFlux)
DownstreamFluxTot=erosionPP*sum(cap2_bins)/(pi*0.0625^2);


%%
%DownstreamFlux=DownstreamFlux*OxyPercent+((1-OxyPercent)*Data.DownstreamFlux);
figure
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),DownstreamFlux,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Downstream End Cap Flux','FontSize',18)
xlabel('X [m]','FontSize',18)
ylabel('Y [m]','FontSize',18)
zlabel('Z [m]','FontSize',18)
% axis([xMin xMax yMin yMax zMin zMax])
% cb1 = colorbar
set(gcf,'color','w')
set(gca,'fontsize',18)
% set(cb1,'fontsize',14);
axis equal
c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',18)
caxis([0 8E18])

% figure
% histogram(charge0(hitCap1))
% figure
% histogram(charge0(hitCap2))


%InnerCap2;


%%
% figure
% s3=scatter(x0(hitCap2),y0(hitCap2));
% s3.MarkerFaceAlpha =weight0(hitCap1);

%%

Data.GrossErosion=(erosionPP.*grossEro0./ElementArea);
Data.GrossDeposition=erosionPP.*grossDep0./ElementArea;
Data.NetErosion=erosionPP.*(grossEro0-grossDep0)./ElementArea;
Data.NetDeposition=erosionPP.*netDep./ElementArea;
Data.DownstreamFlux=erosionPP*cap2_bins./area(cap2);

%%
% figure(4)
% histogram(charge(end,:))

% figure(123)
% histogram(distTraveled(redep),0:0.0002:0.5)
% % hold on
% % histogram(distTraveled(redep2),0:0.0002:0.5)
% total_redeposition_fraction = [total_redeposition_fraction; length(redep)/length(x0)];
% local_redep_fraction = [local_redep_fraction; length(redep_local)/length(x0)];
% prompt_redep_fraction = [prompt_redep_fraction; length(redep_prompt)/length(x0) ];
% net_erosion_fraction = [net_erosion_fraction; (length(x0) - length(redep))/length(x0) ];
% % end
% close all
% figure(1)
% plot(angles,total_redeposition_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Total Re-deposited Fraction of Eroded W')
% axis([0 90 0.7 0.9])
% set(gca,'fontsize',16)
% 
% figure(2)
% plot(angles,local_redep_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Locally Re-deposited Fraction of Eroded W')
% axis([0 90 0 0.04])
% set(gca,'fontsize',16)
% 
% figure(3)
% plot(angles,prompt_redep_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Promptly Re-deposited Fraction of Eroded W')
% axis([0 90 0.25 0.7])
% set(gca,'fontsize',16)
% 
% figure(4)
% plot(angles,net_erosion_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Net Eroded Fraction of Eroded W')
% axis([0 90 0.1 0.25])
% set(gca,'fontsize',16)

aveSpyl = ncread(file,'aveSpyl');
spylCounts = cast(ncread(file,'spylCounts'),'single');
sumWeightStrike = ncread(file,'sumWeightStrike');
surfSputtDist = ncread(file,'surfSputtDist');

eDist = ncread(file,'surfEDist');
e_grid = 10:10:1000;
s1=reshape(sum(eDist,1),100,length(surface));
s2 = s1.*e_grid'./sum(s1);
meanE = sum(s2);

sp1=reshape(sum(surfSputtDist,1),100,length(surface));
sp2 = sp1.*e_grid'./sum(sp1);
mean_spE = sum(sp2);

yield =aveSpyl./spylCounts;

figure(201)
patch(transpose(Z_cyl),transpose(THETA),0*transpose(Z_cyl),meanE,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
set(gca, 'FontSize',18)
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
ylim([0 360])
xlim([0.97 1.3456])
set(gca, 'FontSize',18)
set(gcf,'color','w')
title({'meanE'},'FontSize',18);
hold off

figure(202)
patch(transpose(Z_cyl),transpose(THETA),0*transpose(Z_cyl),yield,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
set(gcf,'color','w')
ylim([0 360])
xlim([0.97 1.3456])
set(gca, 'FontSize',18)
title({'yield'},'FontSize',18);
hold off

figure(203)
patch(transpose(Z_cyl),transpose(THETA),0*transpose(Z_cyl),mean_spE,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

set(gca, 'FontSize',18)
c=colorbar;
ylabel('Angle [Deg]','FontSize',18);
xlabel('Z [m]','FontSize',18);
set(gcf,'color','w')
title({'mean sputtered energy'},'FontSize',18);
ylim([0 360])
xlim([0.97 1.3456])
set(gca, 'FontSize',18)
hold off

%% Density of charge states

specFile = strcat(pwd,'/spec.nc');
Chargedens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
gridY = ncread(specFile,'gridY');
gridZ = ncread(specFile,'gridZ');

figure
slice1 = sum(Chargedens(:,:,:,1:1),4);
slice1 = reshape(slice1,length(gridR),length(gridY),length(gridZ));
RadialSlice1=slice1(:,:,5999); % Target
% RadialSlice1=slice1(:,:,5000); % ICRF Antenna
RadialSlice1=slice1(:,:,749); % Dump piscesRF
%RadialSlice1=sum(slice1(:,:,5875:5905)); %5905
% RadialSlice1=sum(slice1(:,:,1815:1884)); %5905
% RadialSlice1=reshape(RadialSlice1,1884);
slice1 = sum(slice1,2);
slice1 = reshape(slice1,length(gridR),length(gridZ));
dV = (gridR(2) - gridR(1))*(gridZ(2) - gridZ(1))*(gridY(2) - gridY(1));
slice1 = 1/1e5*1e-10/dV*slice1./length(gridY);
p = pcolor(gridR,gridZ,(slice1'));
p.EdgeColor = 'none';
hold on 
set(gca, 'ColorScale', 'linear', 'FontSize',18)
set(gcf,'color','w')
xlabel('X [cm]')
ylabel('Z [m]')
c=colorbar;
ylabel(c,'Al^{0+} density [arb. units]','FontSize',18)
% % plot([-0.05, 0.005],[-0.00866025403784439,0.00866025403784439],'w','lineWidth',2)
% % axis equal
% % axis([-0.015 0.015 -0.01 0.02])
% % title({'W Impurity Density (All Charges)', 'Averaged in y-direction [m^{-3}]'})
% % xlabel('x [m]') % x-axis label
% % ylabel('z [m]') % y-axis label
% % set(gca,'fontsize',16)
% % colorbar
% % set(gca,'ColorScale','log')
% % caxis([0 200000])

% plots radial
RadialDensity1=RadialSlice1*1e-9*erosion_rate/dV/nP;
AverageRadialDensity1=movmean(RadialDensity1(35,:),1);
% AverageRadialDensity1=mean(RadialDensity1,2);


% for ii=1:(length(xflux)-1)
%     
% Cap2radius=sqrt(x0(hitCap2).^2+y0(hitCap2).^2);
% Cap2radius=find(Cap2radius>=xflux(ii) & Cap2radius<=xflux(ii+1));
% Sumhitcap2=sum(weight0(hitCap2(Cap2radius)));
% AverageRadialDensity1(ii)=slice1./((pi*xflux(ii+1)^2)-(pi*xflux(ii)^2)); %0.020m is radius that maps downstream to target
% 
% end


figure;
plot(gridR,AverageRadialDensity1, 'LineWidth',2)
xlabel('x [m]')
ylabel('Density [Particles/m^3]','FontSize',18, 'LineWidth',2)
% title({' Al density at z=4.14 m for the', 'High Density Magnetized Case'})
set(gca, 'LineWidth',2, 'FontSize',18)
% title({' Al density at z=1.75 m for the', 'High Density Magnetized Case'})
title({' Radial Al density at z=0.64 m'})
xlim([0 0.06])

figure
h=pcolor(gridR,gridY,RadialSlice1'); %counts at axial point
h.EdgeColor = 'none';
title('Counts')
xlim([-0.06 0.06])
ylim([-0.06 0.06])
set(gcf,'color','w')
title({' Radial Al density at z=0.64 m'})
colorbar

figure
h=pcolor(gridR,gridY, RadialDensity1'); %counts at axial point
c=colorbar;
h.EdgeColor = 'none';
% title('Al ions Impurity Density at Z=4.14 m')
title('Al ions Impurity Density at Z=3.09 m')
xlim([-0.06 0.06])
ylim([-0.06 0.06])
xlabel('X [cm]')
ylabel('Y [cm]')
ylabel(c, 'Density [Particles/m^3]','FontSize',18,'LineWidth',2)
set(gca, 'LineWidth',2, 'FontSize',18)
set(gcf,'color','w')

