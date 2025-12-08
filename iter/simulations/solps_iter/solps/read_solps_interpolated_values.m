% close all
% clear all
do_plot = 1;

filename = 'interpolated_values.nc'
atomic_number = ncread(filename,'atomic_number');
charge_number = ncread(filename,'charge_number');
mass_number = ncread(filename,'mass_number');

r_inner_target = ncread(filename,'r_inner_target');
z_inner_target = ncread(filename,'z_inner_target');
r_inner_target_midpoints = ncread(filename,'r_inner_target_midpoints');
z_inner_target_midpoints = ncread(filename,'z_inner_target_midpoints');
rmrs_inner_target = ncread(filename,'rmrs_inner_target');
rmrs_inner_target_midpoints = ncread(filename,'rmrs_inner_target_midpoints');
Bmag_inner_target = ncread(filename,'Bmag_inner_target');
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,Bmag_inner_target,'R-R_{sep}','B [T]',{'B-field Magnitude','Inner Divertor'})
Bangle_inner_target = ncread(filename,'Bangle_inner_target');
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,Bangle_inner_target,'R-R_{sep}','Angle [degrees]',{'B-field Angle','Inner Divertor'})
te_inner_target = ncread(filename,'te_inner_target');
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,te_inner_target,'R-R_{sep}','T_e [eV]',{'Electron Temperature','Inner Divertor'})
ti_inner_target = ncread(filename,'ti_inner_target');
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,ti_inner_target,'R-R_{sep}','T_i [eV]',{'Ion Temperature','Inner Divertor'})
ne_inner_target = ncread(filename,'ne_inner_target');
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,ne_inner_target,'R-R_{sep}','n_e [m^{-3}]',{'Electron Density','Inner Divertor'})
flux_inner_target = abs(ncread(filename,'flux_inner_target'));
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,flux_inner_target,'R-R_{sep}','\Gamma_i [m^{-2}s^{-1}]',{'Ion Flux','Inner Divertor'})
ni_inner_target = ncread(filename,'ni_inner_target');
    ret = line_plot(do_plot,rmrs_inner_target_midpoints,ni_inner_target,'R-R_{sep}','n_i [m^{-3}]',{'Ion Density','Inner Divertor'})
r_outer_target = ncread(filename,'r_outer_target');
z_outer_target = ncread(filename,'z_outer_target');
r_outer_target_midpoints = ncread(filename,'r_outer_target_midpoints');
z_outer_target_midpoints = ncread(filename,'z_outer_target_midpoints');
rmrs_outer_target = ncread(filename,'rmrs_outer_target');
rmrs_outer_target_midpoints = ncread(filename,'rmrs_outer_target_midpoints');
Bmag_outer_target = ncread(filename,'Bmag_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,Bmag_outer_target,'R-R_{sep}','B [T]',{'B-field Magnitude','outer Divertor'})
Bangle_outer_target = ncread(filename,'Bangle_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,Bangle_outer_target,'R-R_{sep}','Angle [degrees]',{'B-field Angle','outer Divertor'})
te_outer_target = ncread(filename,'te_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,te_outer_target,'R-R_{sep}','T_e [eV]',{'Electron Temperature','outer Divertor'})
ti_outer_target = ncread(filename,'ti_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,ti_outer_target,'R-R_{sep}','T_i [eV]',{'Ion Temperature','outer Divertor'})
ne_outer_target = ncread(filename,'ne_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,ne_outer_target,'R-R_{sep}','n_e [m^{-3}]',{'Electron Density','outer Divertor'})
flux_outer_target = ncread(filename,'flux_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,flux_outer_target,'R-R_{sep}','\Gamma_i [m^{-2}s^{-1}]',{'Ion Flux','outer Divertor'})
ni_outer_target = ncread(filename,'ni_outer_target');
    ret = line_plot(do_plot,rmrs_outer_target_midpoints,ni_outer_target,'R-R_{sep}','n_i [m^{-3}]',{'Ion Density','outer Divertor'})

gridr = ncread(filename,'gridr');
gridz = ncread(filename,'gridz');

values = ncread(filename,'values');
    ret = pcolor_plot(do_plot,gridr,gridz,values,'r [m]','z [m]',{'Potential [V]','He Discharge'})
te = ncread(filename,'te');
    ret = pcolor_plot(do_plot,gridr,gridz,te,'r [m]','z [m]',{'T_e [eV]','He Discharge'})
ti = ncread(filename,'ti');
    ret = pcolor_plot(do_plot,gridr,gridz,ti,'r [m]','z [m]',{'T_i [eV]','He Discharge'})
ne = ncread(filename,'ne');
    ret = pcolor_plot(do_plot,gridr,gridz,ne,'r [m]','z [m]',{'n_e [m^{-3}]','He Discharge'})
ni = ncread(filename,'ni');
    ret = pcolor_plot(do_plot,gridr,gridz,ni,'r [m]','z [m]',{'n_i [m^{-3}]','He Discharge'})
flux = ncread(filename,'flux');

mass = ncread(filename,'mass');
    ret = pcolor_plot(do_plot,gridr,gridz,mass,'r [m]','z [m]',{'Mass [amu]','He Discharge'})
charge = ncread(filename,'charge');
    ret = pcolor_plot(do_plot,gridr,gridz,charge,'r [m]','z [m]',{'Charge [#]','He Discharge'})
Br = ncread(filename,'Br');
    ret = pcolor_plot(do_plot,gridr,gridz,Br,'r [m]','z [m]',{'B_r [T]','He Discharge'})
Bt = ncread(filename,'Bt');
    ret = pcolor_plot(do_plot,gridr,gridz,Bt,'r [m]','z [m]',{'B_t [T]','He Discharge'})
Bz = ncread(filename,'Bz');
    ret = pcolor_plot(do_plot,gridr,gridz,Bz,'r [m]','z [m]',{'B_z [T]','He Discharge'})

Bmag = ncread(filename,'Bmag');
vr = ncread(filename,'vr');
    ret = pcolor_plot(do_plot,gridr,gridz,vr,'r [m]','z [m]',{'v_r [m/s]','He Discharge'})
vt = ncread(filename,'vt');
    ret = pcolor_plot(do_plot,gridr,gridz,vt,'r [m]','z [m]',{'v_t [m/s]','He Discharge'})
vz = ncread(filename,'vz');
    ret = pcolor_plot(do_plot,gridr,gridz,vz,'r [m]','z [m]',{'v_z [m/s]','He Discharge'})
Er = ncread(filename,'Er');
    ret = pcolor_plot(do_plot,gridr,gridz,vr,'r [m]','z [m]',{'E_r [V/m]','He Discharge'})
Ez = ncread(filename,'Ez');
    ret = pcolor_plot(do_plot,gridr,gridz,vz,'r [m]','z [m]',{'E_z [V/m]','He Discharge'})

gradTer = ncread(filename,'gradTer');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTer,'r [m]','z [m]',{'B dot \nabla T_{er} [eV/m]','He Discharge'})
    hl = title({'$(\hat{B} \cdot{} \nabla T_{e}) \hat{B}_r [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTe = ncread(filename,'gradTe');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTe,'r [m]','z [m]',{'\hat{B} \dot{} \nabla T_{e} [eV/m]','He Discharge'})
    hl = title({'$\hat{B} \cdot{} \nabla T_{e} [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTet = ncread(filename,'gradTet');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTet,'r [m]','z [m]',{'B dot \nabla T_{et} [eV/m]','He Discharge'})
    hl = title({'$(\hat{B} \cdot{} \nabla T_{e}) \hat{B}_t [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTez = ncread(filename,'gradTez');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTez,'r [m]','z [m]',{'B dot \nabla T_{ez} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{e}) \hat{B}_z [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTi = ncread(filename,'gradTi');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTi,'r [m]','z [m]',{'\nabla T_{i} [eV/m]','He Discharge'})
    hl = title({'$\hat{B} \cdot{} \nabla T_{i} [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTir = ncread(filename,'gradTir');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTir,'r [m]','z [m]',{'B dot \nabla T_{ir} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{i}) \hat{B}_r [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTit = ncread(filename,'gradTit');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTit,'r [m]','z [m]',{'B dot \nabla T_{it} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{i}) \hat{B}_t [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTiz = ncread(filename,'gradTiz');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTiz,'r [m]','z [m]',{'B dot \nabla T_{iz} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{i}) \hat{B}_z [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');

% figure(1)
% h = pcolor(gridr,gridz,te')
% h.EdgeColor = 'none';
% colorbar
% set(gca, 'ColorScale', 'log')

function ret = line_plot(do_plot,x,y,xl,yl,t)
if do_plot
    figure
    plot(x,y,'LineWidth',2)
    title(t)
    ylabel(yl)
    xlabel(xl)
    set(gca,'FontSize',14)
    ret = 1;
else
    ret = 0;
end
end

function ret = pcolor_plot(do_plot,r,z,y,xl,yl,t)
if do_plot
    figure
h = pcolor(r,z,y')
h.EdgeColor = 'none';
colorbar

    title(t)
    ylabel(yl)
    xlabel(xl)
    set(gca,'FontSize',14)
%     axis equal
%     axis([min(r) max(r) min(z) max(z)])
    ret = 1;
else
    ret = 0;
end
end
