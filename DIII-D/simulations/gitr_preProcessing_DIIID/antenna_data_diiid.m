
ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;


%% Read profiles from SOLPS
R=ncread('profilesDIIID.nc','x');
z=ncread('profilesDIIID.nc','z');

% % B-field
% bz=ncread('profilesDIIID.nc','bz');
% bt=ncread('profilesDIIID.nc','bt');
% br=ncread('profilesDIIID.nc','br');

% br(isnan(br))=0;bz(isnan(bz))=0;bt(isnan(bt))=0;

% Densities
ne=ncread('profilesDIIID.nc','ne');
ni=ncread('profilesDIIID.nc','ni');



% Temperatures
te=ncread('profilesDIIID.nc','te');
ti=ncread('profilesDIIID.nc','ti');

% ne(isnan(ne))=0;ni(isnan(ni))=0;te(isnan(te))=0;ti(isnan(ti))=0;

% Velocity
vt=ncread('profilesDIIID.nc','vt');
vr=ncread('profilesDIIID.nc','vr');
vz=ncread('profilesDIIID.nc','vz');




%Centroids and surface normals
r_centroid= sqrt(centroid(:,1).^2+centroid(:,2).^2);
% norm_vec=readmatrix('norm_vec.csv');

figure; 
%%
% 2D density
figure; imagesc(R, z, ne);
set(gca,'YDir','normal','FontName','Times','FontSize',24);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
cb = colorbar; ylabel(cb,'$n_e\,[\mathrm{m}^{-3}]$','Interpreter','latex');
title('Masked $n_e$ (\psi_N<0.8)','Interpreter','latex');
pbaspect([2 1 1]);
hold on;
scatter(r_centroid,centroid(:,3))
plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);
    % contour(X, Y, psiN, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2)
%%
%Interpolations of profiles onto the geometry
ne_surf=interpn(R,z,ne,r_centroid, centroid(:,3));
ni_surf=interpn(R,z,ni,r_centroid, centroid(:,3));
te_surf=interpn(R,z,te,r_centroid, centroid(:,3));
ti_surf=interpn(R,z,ti,r_centroid, centroid(:,3));


% br_surf=interpn(R,z,br,r_centroid, centroid(:,3));
% bt_surf=interpn(R,z,bt,r_centroid, centroid(:,3));
% bz_surf=interpn(R,z,bz,r_centroid, centroid(:,3));
vp=sqrt(vt.^2+vz.^2+vr.^2);

vr_surf=interpn(R,z,vr,r_centroid, centroid(:,3));
vt_surf=interpn(R,z,vt,r_centroid, centroid(:,3));
vz_surf=interpn(R,z,vz,r_centroid, centroid(:,3));




vp_surf=interpn(R,z,vz,r_centroid, centroid(:,3));
% vp_surf=interpn(R,z,vp,r_centroid, centroid(:,3));

figure;
patch(transpose(X),transpose(Y),transpose(Z),ne_surf,'FaceAlpha',1,'EdgeAlpha', 0.3 );
title('ITER Limiter surface density')        
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')


%Interpolations of profiles onto the geometry
r_efit=r_efit;
z_efit=z_efit;
% figure; imagesc(z_efit,r_efit,Br)
% hold on;
% scatter( centroid(:,3), r_centroid)

br_surf=interpn(z_efit,r_efit,Br, centroid(:,3),r_centroid);
bt_surf=interpn(z_efit,r_efit,Bt, centroid(:,3),r_centroid);
bz_surf=interpn(z_efit,r_efit,Bz, centroid(:,3),r_centroid);

phi_centroid = atan2(centroid(:,2),centroid(:,1));


bx=double(br_surf.*cos(phi_centroid)-bt_surf.*sin(phi_centroid));
by=double(br_surf.*sin(phi_centroid)+bt_surf.*cos(phi_centroid));
bz=double(bz_surf);


b_mag = sqrt(bx.^2 + by.^2 + bz.^2);


ubx = bx./b_mag;
uby = by./b_mag;
ubz = bz./b_mag;


norm_vec_mag = sqrt(norm_vec(:,1).^2+ norm_vec(:,2).^2+norm_vec(:,3).^2);

unorm_vec = norm_vec./norm_vec_mag;

theta = acos(unorm_vec(:,1).*ubx + unorm_vec(:,2).*uby + unorm_vec(:,3).*ubz);



ii=find(theta>pi/2);
theta(ii)=abs(theta(ii)-pi);


% figure;hold on; patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','none')%'none')
figure; quiver3(centroid(:,1),centroid(:,2),centroid(:,3),norm_vec(:,1),norm_vec(:,2),norm_vec(:,3)); hold on;
quiver3(centroid(:,1),centroid(:,2),centroid(:,3),bx,by,bz);



figure;histogram(theta)



