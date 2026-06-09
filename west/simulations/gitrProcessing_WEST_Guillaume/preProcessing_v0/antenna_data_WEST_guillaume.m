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
ni=ncread('profilesWEST.nc','ni');

% Temperatures
te=ncread('profilesWEST.nc','te');
ti=ncread('profilesWEST.nc','ti');

% ne(isnan(ne))=0;ni(isnan(ni))=0;te(isnan(te))=0;ti(isnan(ti))=0;

% Velocity
vt=ncread('profilesWEST.nc','vt');
vr=ncread('profilesWEST.nc','vr');
vz=ncread('profilesWEST.nc','vz');


%Centroids and surface normals
r_centroid= sqrt(centroid(:,1).^2+centroid(:,2).^2);
% norm_vec=readmatrix('norm_vec.csv');

% %Interpolations of profiles onto the geometry
% ne_surf=interpn(R,z,ne,r_centroid, centroid(:,3));
% ni_surf=interpn(R,z,ni,r_centroid, centroid(:,3));
% te_surf=interpn(R,z,te,r_centroid, centroid(:,3));
% ti_surf=interpn(R,z,ti,r_centroid, centroid(:,3));
% 
% 
% br_surf=interpn(R,z,br,r_centroid, centroid(:,3));
% bt_surf=interpn(R,z,bt,r_centroid, centroid(:,3));
% bz_surf=interpn(R,z,bz,r_centroid, centroid(:,3));
% 
% vr_surf=interpn(R,z,vr,r_centroid, centroid(:,3));
% vt_surf=interpn(R,z,vt,r_centroid, centroid(:,3));
% vz_surf=interpn(R,z,vz,r_centroid, centroid(:,3));
% 
% vp_surf=interpn(R,z,vp',r_centroid, centroid(:,3));
% 
% 



phi_centroid = atan2(centroid(:,2),centroid(:,1));



%% GITR geometry

% Sheath calculations
disp('>>> Extrapolating B-field')


        Bx_data=readmatrix('../Guillaume_cases/ne_lim_1p1e16/Bx_out.txt');
        By_data=readmatrix('../Guillaume_cases/ne_lim_1p1e16/By_out.txt');
        Bz_data=readmatrix('../Guillaume_cases/ne_lim_1p1e16/Bz_out.txt');

        xx=Bx_data(:,1);
        yy=By_data(:,2);
        zz=Bz_data(:,3);

        Bx_g=Bx_data(:,4);
        By_g=By_data(:,4);
        Bz_g=Bz_data(:,4);

    
        for i=1:length(centroid)
            distance = sqrt((centroid(i,1) - xx).^2 + (centroid(i,2) - yy).^2 + (centroid(i,3) - zz).^2);

            [M I] = min(distance);
            bx(i)= Bx_g(I);
            by(i)= By_g(I);
            bz(i)= Bz_g(I);
         end


b_mag = sqrt(bx.^2 + by.^2 + bz.^2);

% writematrix(b_mag, "b_mag.csv");

ubx = bx./b_mag;
uby = by./b_mag;
ubz = bz./b_mag;


vx=double(vr_surf.*cos(phi_centroid)-vt_surf.*sin(phi_centroid));
vy=double(vr_surf.*sin(phi_centroid)+vt_surf.*cos(phi_centroid));
vz=double(vz_surf);

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



