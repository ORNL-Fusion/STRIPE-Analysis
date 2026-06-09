ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;


%% Read profiles from SOLEDGE
% Grid
R=ncread('profilesWEST.nc','x');
z=ncread('profilesWEST.nc','z');

% B-field
bz=ncread('profilesWEST.nc','bz');
bt=ncread('profilesWEST.nc','bt');
br=ncread('profilesWEST.nc','br');

% br(isnan(br))=0;bz(isnan(bz))=0;bt(isnan(bt))=0;

% Densities
ne=ncread('profilesWEST.nc','ne');
ni=ncread('profilesWEST.nc','no');
% no=ncread('profilesWEST.nc','no');




% Temperatures
te=ncread('profilesWEST.nc','te');
ti=ncread('profilesWEST.nc','to');
% to=ncread('profilesWEST.nc','to');

% ne(isnan(ne))=0;ni(isnan(ni))=0;te(isnan(te))=0;ti(isnan(ti))=0;

% Velocity
vto=ncread('profilesWEST.nc','vto');
vro=ncread('profilesWEST.nc','vro');
vzo=ncread('profilesWEST.nc','vzo');

vto(isnan(vto))=0;

vti=ncread('profilesWEST.nc','vti');
vri=ncread('profilesWEST.nc','vri');
vzi=ncread('profilesWEST.nc','vzi');


%Centroids and surface normals
r_centroid= sqrt(centroid(:,1).^2+centroid(:,2).^2);
% norm_vec=readmatrix('norm_vec.csv');

%Interpolations of profiles onto the geometry
ne_surf=interpn(R,z,ne,r_centroid, centroid(:,3));
ni_surf=interpn(R,z,ni,r_centroid, centroid(:,3));
% no_surf=interpn(R,z,no,r_centroid, centroid(:,3));
te_surf=interpn(R,z,te,r_centroid, centroid(:,3));
ti_surf=interpn(R,z,ti,r_centroid, centroid(:,3));
% to_surf=interpn(R,z,to,r_centroid, centroid(:,3));


br_surf=interpn(R,z,br,r_centroid, centroid(:,3));
bt_surf=interpn(R,z,bt,r_centroid, centroid(:,3));
bz_surf=interpn(R,z,bz,r_centroid, centroid(:,3));

vri_surf=interpn(R,z,vri,r_centroid, centroid(:,3));
vti_surf=interpn(R,z,vti,r_centroid, centroid(:,3));
vzi_surf=interpn(R,z,vzi,r_centroid, centroid(:,3));

vro_surf=interpn(R,z,vro,r_centroid, centroid(:,3));
vto_surf=interpn(R,z,vto,r_centroid, centroid(:,3));
vzo_surf=interpn(R,z,vzo,r_centroid, centroid(:,3));



vp_surf=interpn(R,z,vto,r_centroid, centroid(:,3));
% vp_surf=interpn(R,z,vp',r_centroid, centroid(:,3));





phi_centroid = atan2(centroid(:,2),centroid(:,1));


bx=double(br_surf.*cos(phi_centroid)-bt_surf.*sin(phi_centroid));
by=double(br_surf.*sin(phi_centroid)+bt_surf.*cos(phi_centroid));
bz=double(bz_surf);


b_mag = sqrt(bx.^2 + by.^2 + bz.^2);

% writematrix(b_mag, "b_mag.csv");

ubx = bx./b_mag;
uby = by./b_mag;
ubz = bz./b_mag;


vx=double(vro_surf.*cos(phi_centroid)-vto_surf.*sin(phi_centroid));
vy=double(vro_surf.*sin(phi_centroid)+vto_surf.*cos(phi_centroid));
vz=double(vzo_surf);

v_mag=sqrt(vx.^2+vy.^2+vz.^2);
uvx=vx./v_mag;
uvy=vy./v_mag;
uvz=vy./v_mag;

o8plus_flux_surf=ni_surf.*v_mag;


norm_vec_mag = sqrt(norm_vec(:,1).^2+ norm_vec(:,2).^2+norm_vec(:,3).^2);

unorm_vec = norm_vec./norm_vec_mag;

theta = acos(unorm_vec(:,1).*ubx + unorm_vec(:,2).*uby + unorm_vec(:,3).*ubz);



ii=find(theta>pi/2);
theta(ii)=abs(theta(ii)-pi);

% % Write profiles data on the geometry
% writematrix(ne_surf, 'ne_surf.csv');
% writematrix(te_surf, 'te_surf.csv');
% writematrix(ti_surf, 'ti_surf.csv');
% writematrix(theta,   'theta.csv'  );
writematrix(o8plus_flux_surf, 'o8plus_flux_surf.csv');


% figure;hold on; patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','none')%'none')
figure; quiver3(centroid(:,1),centroid(:,2),centroid(:,3),norm_vec(:,1),norm_vec(:,2),norm_vec(:,3)); hold on;
quiver3(centroid(:,1),centroid(:,2),centroid(:,3),bx,by,bz);



figure;histogram(theta)



