ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;



%Centroids and surface normals
r_centroid= sqrt(centroid(:,1).^2+centroid(:,2).^2);

phi_centroid = atan2(centroid(:,2),centroid(:,1));



%% GITR geometry

% Sheath calculations
disp('>>> Extrapolating B-field')


        Bx_data=readmatrix('../Guillaume_cases/ne_lim_7p3e18/Bx_out.txt');
        By_data=readmatrix('../Guillaume_cases/ne_lim_7p3e18/By_out.txt');
        Bz_data=readmatrix('../Guillaume_cases/ne_lim_7p3e18/Bz_out.txt');

        xx=Bx_data(:,1);
        yy=Bx_data(:,2);
        zz=Bx_data(:,3);

        Bx_g=Bx_data(:,4);
        By_g=By_data(:,4);
        Bz_g=Bz_data(:,4);

    
        for i=1:length(centroid)
            distance = sqrt((centroid(i,1) - xx).^2 + (centroid(i,2) - yy).^2 + (centroid(i,3) - zz).^2);

            [M I] = min(distance);
            bx1(i)= Bx_g(I);
            by1(i)= By_g(I);
            bz1(i)= Bz_g(I);
         end

bx_surf=bx1';
by_surf=by1';
bz_surf=bz1';
b_mag = sqrt(bx_surf.^2 + by_surf.^2 + bz_surf.^2);


ubx = bx_surf./b_mag;
uby = by_surf./b_mag;
ubz = bz_surf./b_mag;


norm_vec_mag = sqrt(norm_vec(:,1).^2+ norm_vec(:,2).^2+norm_vec(:,3).^2);

unorm_vec = norm_vec./norm_vec_mag;

theta = acos(unorm_vec(:,1).*ubx + unorm_vec(:,2).*uby + unorm_vec(:,3).*ubz);



ii=find(theta>pi/2);
theta(ii)=abs(theta(ii)-pi);


% figure;hold on; patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','none')%'none')
figure; quiver3(centroid(:,1),centroid(:,2),centroid(:,3),norm_vec(:,1),norm_vec(:,2),norm_vec(:,3)); hold on;
quiver3(centroid(:,1),centroid(:,2),centroid(:,3),bx_surf,by_surf,bz_surf);

figure;histogram(theta)



