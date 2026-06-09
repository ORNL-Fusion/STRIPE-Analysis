file = 'profilesITER.nc';
% file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profilesProtoMPEX.nc';
x = ncread(file,'x');
z = ncread(file,'z');
ni = ncread(file,'ni');
ne = ncread(file,'ne');
ti = ncread(file,'ti');
te = ncread(file,'te');
vr = ncread(file,'vr');
vt = ncread(file,'vt');
vz = ncread(file,'vz');
br = ncread(file,'br');
bt = ncread(file,'bt');
bz = ncread(file,'bz');
% vz1= ncread(file1,'vz');
figure; imagesc(z,x,br)

% Physical constants:
% =========================================================================

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 20.18;





%Centroids and surface normals
r_centroid= sqrt(centroid(:,1).^2+centroid(:,2).^2);


%Interpolations of profiles onto the geometry
r_efit=x;
z_efit=z;


figure; imagesc(x,z,ne)
hold on;
plot(g.lim(1,:), g.lim(2,:), 'r');
hold on;
scatter(  r_centroid, centroid(:,3))

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



