close all
m=readmatrix('Shot195196.csv');
idx=find(m(:,1)~=0);
m=m(idx,:);
%plot(m(:,1),m(:,2))
R = 1:0.01:3;
z = -1:0.01:1;
[R,z] = meshgrid(R,z);
ad = ((0.30912-0.52748)*R-(2.306-2.2309)*z+2.306*0.52748-0.30912*2.2309)/sqrt((0.30912-0.52748)^2+(2.306-2.2309)^2);
pcolor(R,z,ad)
shading interp
title('distance to antenna [m]')
colorbar;
xlabel('R')
ylabel('z')
vq=interp1(m(:,1),m(:,2),ad,'pchip',nan);
figure;
pcolor(R,z,vq)
shading interp
title('density parametrized by distance to antenna [m^{-3}]')
colorbar;
xlabel('R')
ylabel('z')
