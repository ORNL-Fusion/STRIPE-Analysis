
% clc;
% close all;
% clear all;

disp('>>>> Reading g-eqdsk file');


% Read an g-eqdsk file

filename= 'g900003.00230_ITER_15MA_eqdsk16VVHR.txt';
% filename= 'MOB-348s_eqdsk.txt';
g = readg_g3d(filename);

figure; plot_gfile(g)

set(gca,'FontName','times','fontSize',24);
ylabel('$z$ [m$^{-3}$]','interpreter','Latex','fontSize',24)
xlabel('$r$ [m]','interpreter','latex','fontSize',24)
%% Calculate B-field and write it in .csv and .nc format
disp('>>>> Writing the B-field profiles');

nR = 400; %for ITER
nZ = 400;
offset = 3;
r_efit=linspace(g.r(1+offset),g.r(end-offset),nR);
z_efit=linspace(g.z(1+offset),g.z(end-offset),nZ);
[r2D,z2D] = meshgrid(r_efit,z_efit);
[Bout,ierr] = bfield_geq_bicub(g,r2D(:),z2D(:));

[psiN,psi] = calc_psiN(g,r2D(:),z2D(:),[]);
psi=reshape(psi,[nZ,nR]);
psiN=reshape(psiN,[nZ,nR]);

Br = reshape(Bout.br,[nZ,nR]);
Bt = reshape(Bout.bphi,[nZ,nR]);
Bz = reshape(Bout.bz,[nZ,nR]);


ncid = netcdf.create(('./bfield_iter.nc'),'NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);

dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);

gridZnc = netcdf.defVar(ncid,'z','float',dimZ);

brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,r_efit);
netcdf.putVar(ncid,gridZnc,z_efit);


netcdf.putVar(ncid,brnc,Br');
netcdf.putVar(ncid,btnc,Bt');
netcdf.putVar(ncid,bznc,Bz');

netcdf.close(ncid);

disp('>>>> Calculated the B-field profile');


%% === Find OUTER (LFS) WALL INTERSECTION at the MIDPLANE (Z = 0) ==========
disp('>>>> Determining outer wall intersection at Z = 0 (midplane) ...');

if isfield(g,'lim') && ~isempty(g.lim)
    wall_r = g.lim(1,:); 
    wall_z = g.lim(2,:);
    good = isfinite(wall_r) & isfinite(wall_z);
    wall_r = wall_r(good);
    wall_z = wall_z(good);

    % Identify segments where Z changes sign around 0 → crossing midplane
    sgn = sign(wall_z);
    icross = find(diff(sgn) ~= 0);

    if ~isempty(icross)
        R_cross = [];
        for ii = icross
            z1 = wall_z(ii);  z2 = wall_z(ii+1);
            r1 = wall_r(ii);  r2 = wall_r(ii+1);
            % Interpolate R where Z=0
            r_cross = interp1([z1 z2],[r1 r2],0,'linear','extrap');
            R_cross = [R_cross, r_cross];
        end
        % Outer wall → largest R crossing
        R_wall_mid = max(R_cross);
        fprintf('Outer wall intersection at midplane: R = %.3f m, Z = 0.000 m\n', R_wall_mid);

        % Plot and label it
        hold on;
        plot(R_wall_mid, 0, 'ro', 'MarkerFaceColor','r','MarkerSize',10);
        text(R_wall_mid+0.03, 0, sprintf('R=%.3f m', R_wall_mid), ...
             'FontSize',14,'Color','r','FontWeight','bold');
    else
        warning('Wall contour does not cross Z=0. Using global max-R as fallback.');
        [R_wall_mid, iMaxR] = max(wall_r);
        fprintf('Fallback: R = %.3f m (Z = %.3f m)\n', R_wall_mid, wall_z(iMaxR));
        hold on;
        plot(R_wall_mid, wall_z(iMaxR), 'ro','MarkerFaceColor','r','MarkerSize',10);
    end
else
    warning('No g.lim field found in gfile structure.');
end