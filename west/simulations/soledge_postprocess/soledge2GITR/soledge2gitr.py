import os
import h5py
from shapely.geometry import Point, Polygon
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.interpolate import griddata
from matplotlib.path import Path
from read_wall import Surface


def interpolate_and_save_plasma_field(ref_file, mesh_file, bfield_file, data_file, wall_file, openedge_file, s3x_file, nR=1000, nZ=200):
    # Load reference parameters
    with h5py.File(ref_file, 'r') as ref:
        n0 = ref['/n0'][...]
        T0 = ref['/T0'][...]
        c0 = ref['/c0'][...]
        W0 = ref['/W0'][...]
        rho0 = ref['/rho0'][...]
        R0 = ref['/R0'][...]
        a0 = ref['/a0'][...]
        B0 = ref['/B0'][...]

    # Load wall information geometry
    path = '../'
    wall = Surface(wall_file, "2D")
    domain = wall.polygon
    Rwall, Zwall = domain.exterior.xy

    # Load mesh data
    with h5py.File(mesh_file, 'r') as meshEIRENE:
        tri_knots = meshEIRENE['/triangles/tri_knots'][...] - 1  # Convert to 0-based indexing
        Rk = meshEIRENE['/knots/R'][...] / 100
        Zk = meshEIRENE['/knots/Z'][...] / 100

    # Set a mask on the core boundary: soledge inner boundary
    with h5py.File(s3x_file, 'r') as meshS3X:
        psi = meshS3X['config/psi'][...]
        psisep1 = meshS3X['/config/psisep1'][...]
        psisep2 = meshS3X['/config/psisep2'][...]
        psicore = meshS3X['/config/psicore'][...]
        r2D = meshS3X['/config/r'][...]
        z2D = meshS3X['/config/z'][...]

    print(psi.shape)

    # Load field data
    with h5py.File(data_file, 'r') as data:
        temp_e = data['/triangles/spec0/T'][...].flatten() * T0
        dens_e = data['/triangles/spec0/n'][...].flatten() * n0
        parr_flow_e = data['/triangles/spec0/G'][...].flatten() * c0 * n0

        temp_i = data['/triangles/spec1/T'][...].flatten() * T0
        dens_i = data['/triangles/spec1/n'][...].flatten() * n0
        parr_flow_i = data['/triangles/spec1/G'][...].flatten() * c0 * n0

        # Load data for species o_1 to o_8
        species_data = {}
        for i in range(1, 9):
            temp = data[f'/triangles/spec{i+1}/T'][...].flatten() * T0
            dens = data[f'/triangles/spec{i+1}/n'][...].flatten() * n0
            parr_flow = data[f'/triangles/spec{i+1}/G'][...].flatten() * c0 * n0
            species_data[f'o_{i}'] = {'temp': temp, 'dens': dens, 'parr_flow': parr_flow}

    # Load bfield
    with h5py.File(bfield_file, 'r') as data:
        b_r = data['/triangles/Br'][...].flatten() * B0
        b_z = data['/triangles/Bz'][...].flatten() * B0
        b_phi = data['/triangles/Bphi'][...].flatten() * B0

    # Combine R and Z coordinates for Delaunay triangulation
    points = np.vstack((Rk, Zk)).T
    triangle_points = np.array([points[tri] for tri in tri_knots.T])
    centroids = triangle_points.mean(axis=1)

    # Interpolation grid setup
    r = np.linspace(1.8, 3.2, nR)
    z = np.linspace(-1, 1, nZ)
    grid_r, grid_z = np.meshgrid(r, z)
    grid_points = np.vstack((grid_r.flatten(), grid_z.flatten())).T

    # Function to interpolate field data
    def interpolate_field(field_data):
        interpolated_values_linear = griddata(centroids, field_data, grid_points, method='linear')
        interpolated_values_nearest = griddata(centroids, field_data, grid_points, method='nearest')
        interpolated_values = np.where(np.isnan(interpolated_values_linear), interpolated_values_nearest, interpolated_values_linear)
        return interpolated_values.reshape(nZ, nR)

    # Interpolate psi onto the new grid
    psi_interpolated = griddata((r2D.flatten(), z2D.flatten()), psi.flatten(), (grid_r, grid_z), method='linear')
    psi_interpolated_nearest = griddata((r2D.flatten(), z2D.flatten()), psi.flatten(), (grid_r, grid_z), method='nearest')
    psi_grid = np.where(np.isnan(psi_interpolated), psi_interpolated_nearest, psi_interpolated)

    # Core and wall masks
    psi_core_threshold = psicore[0]
    mask_inside_core = psi_grid < psi_core_threshold
    wall_path = Path(np.vstack((Rwall, Zwall)).T)
    mask_outside_wall = ~wall_path.contains_points(grid_points).reshape(nZ, nR)

    # Interpolate fields and apply masks
    temp_e_grid = interpolate_field(temp_e)
    dens_e_grid = interpolate_field(dens_e)
    parr_flow_e_grid = interpolate_field(parr_flow_e)
    temp_i_grid = interpolate_field(temp_i)
    dens_i_grid = interpolate_field(dens_i)
    parr_flow_i_grid = interpolate_field(parr_flow_i)
    
    b_r_grid = interpolate_field(b_r)
    b_z_grid = interpolate_field(b_z)
    b_phi_grid = interpolate_field(b_phi)
    
    temp_e_grid[mask_inside_core] = 0
    dens_e_grid[mask_inside_core] = 0
    parr_flow_e_grid[mask_inside_core] = 0
    temp_i_grid[mask_inside_core] = 0
    dens_i_grid[mask_inside_core] = 0
    parr_flow_i_grid[mask_inside_core] = 0
    
    b_r_grid[mask_inside_core] = 0
    b_z_grid[mask_inside_core] = 0
    b_phi_grid[mask_inside_core] = 0
    
     # Set interpolated values to zero for points outside the wall and inside the core
    fields = [temp_e_grid, dens_e_grid, temp_i_grid, dens_i_grid, parr_flow_e_grid, parr_flow_i_grid, b_r_grid, b_z_grid, b_phi_grid]
    for field in fields:
        field[mask_outside_wall.reshape(nZ, nR)] = 0
#        field[mask_inside_core_reshaped] = 0

    species_grids = {}
    for name, data in species_data.items():
        temp_grid = interpolate_field(data['temp'])
        dens_grid = interpolate_field(data['dens'])
        parr_flow_grid = interpolate_field(data['parr_flow'])

        # Apply core and wall masks
        temp_grid[mask_inside_core] = 0
        dens_grid[mask_inside_core] = 0
        parr_flow_grid[mask_inside_core] = 0
        temp_grid[mask_outside_wall] = 0
        dens_grid[mask_outside_wall] = 0
        parr_flow_grid[mask_outside_wall] = 0

        species_grids[name] = {'temp': temp_grid, 'dens': dens_grid, 'parr_flow': parr_flow_grid}

   

    Bmag = np.sqrt(b_r_grid**2 + b_z_grid**2 + b_phi_grid**2)
    [grad_te_z, grad_te_r] = np.gradient(temp_e_grid, z[1] - z[0], r[1] - r[0])
    [grad_ti_z, grad_ti_r] = np.gradient(temp_i_grid, z[1] - z[0], r[1] - r[0])

    epsilon = 1e-10
    grad_te_parallel = np.where(Bmag < epsilon, 0, (grad_te_r * b_r_grid + grad_te_z * b_z_grid) / (Bmag + epsilon))
    grad_ti_parallel = np.where(Bmag < epsilon, 0, (grad_ti_r * b_r_grid + grad_ti_z * b_z_grid) / (Bmag + epsilon))
    
    parr_flow_o_8_grid = species_grids['o_8']['parr_flow']
    dens_o_8_grid = species_grids['o_8']['dens']
    temp_o_8_grid = species_grids['o_8']['temp']
    
    # # For oxygen 8+ species
    # min_flux_value = -5e20 # Set to the desired minimum
    # max_flux_value = 5e20 #np.max(temp_e_grid)  # Set to the desired maximum
    
    # min_dens_value = 1e14 # Set to the desired minimum
    # max_dens_value = 5e17 #np.max(temp_e_grid)  # Set to the desired maximum
    
    # min_temp_value = 1 # Set to the desired minimum
    # max_temp_value = 1000 #np.max(temp_e_grid)  # Set to the desired maximum
    
    # For D+ species
    min_flux_value = -1e23 # Set to the desired minimum
    max_flux_value = 1e23 #np.max(temp_e_grid)  # Set to the desired maximum
    
    min_dens_value = 1e14 # Set to the desired minimum
    max_dens_value = 5e19 #np.max(temp_e_grid)  # Set to the desired maximum
    
    min_temp_value = 1 # Set to the desired minimum
    max_temp_value = 1000 #np.max(temp_e_grid)  # Set to the desired maximum
    
    max_Br_value = 0.15
    min_Br_value = -0.15
    
    max_Bz_value = 0.15
    min_Bz_value = -0.15
    
    max_Bphi_value = 5.0
    min_Bphi_value = 0.0

    plt.contour(r2D, z2D, psi, linestyles='--', levels=[psicore[0]], label='Y')
    plt.contour(r2D, z2D, psi, linestyles='--',levels=[psisep1[0]],  label='Y')
    plt.contour(r2D, z2D, psi, linestyles='--',levels=[psisep2[0]],  label='Y')
    # plt.pcolormesh(grid_r, grid_z, parr_flow_i_grid, cmap='coolwarm', norm=LogNorm(vmin=min_temp_value)) # Log scale;plasma, jet, inferno, viridis(default), coolwarm, RdBu, Blues
    # plt.pcolormesh(grid_r, grid_z, parr_flow_i_grid, cmap='coolwarm') #Linear scale; plasma, jet, inferno, viridis(default), coolwarm, RdBu, Blues
    # plt.pcolormesh(grid_r, grid_z, parr_flow_o_8_grid, cmap='coolwarm', vmin=min_flux_value, vmax=max_flux_value)
    # plt.pcolormesh(grid_r, grid_z, dens_o_8_grid, cmap='coolwarm',norm=LogNorm(vmin=min_dens_value, vmax=max_dens_value))
    # plt.pcolormesh(grid_r, grid_z, temp_o_8_grid, cmap='coolwarm', norm=LogNorm(vmin=min_temp_value, vmax=max_temp_value))
    
    # plt.pcolormesh(grid_r, grid_z, parr_flow_i_grid, cmap='coolwarm', vmin=min_flux_value, vmax=max_flux_value)
    # plt.pcolormesh(grid_r, grid_z, dens_i_grid, cmap='coolwarm',norm=LogNorm(vmin=min_dens_value, vmax=max_dens_value))
    # plt.pcolormesh(grid_r, grid_z, temp_i_grid, cmap='coolwarm', norm=LogNorm(vmin=min_temp_value, vmax=max_temp_value))
   
    # plt.pcolormesh(grid_r, grid_z, parr_flow_e_grid, cmap='coolwarm', vmin=min_flux_value, vmax=max_flux_value)
    # plt.pcolormesh(grid_r, grid_z, dens_e_grid, cmap='coolwarm',norm=LogNorm(vmin=min_dens_value, vmax=max_dens_value))
    # plt.pcolormesh(grid_r, grid_z, temp_e_grid, cmap='coolwarm', norm=LogNorm(vmin=min_temp_value, vmax=max_temp_value)) 
    
    
    plt.pcolormesh(grid_r, grid_z, b_r_grid, cmap='coolwarm', vmin=min_Br_value, vmax=max_Br_value)
    # plt.pcolormesh(grid_r, grid_z, b_z_grid, cmap='coolwarm', vmin=min_Bz_value, vmax=max_Bz_value)
    # plt.pcolormesh(grid_r, grid_z, b_phi_grid, cmap='coolwarm', vmin=min_Bphi_value, vmax=max_Bphi_value)
   
    
    plt.plot(Rwall, Zwall, 'k', lw=2.5)
    # Color bar with label
    cbar = plt.colorbar()
    # cbar.set_label('Temperature (eV)', fontsize=18)
    cbar.ax.tick_params(labelsize=12)
    
    np.savetxt("wall_coordinates.txt",
           np.column_stack((Rwall, Zwall)),
           header="Rwall   Zwall",
           comments='')

    # Setting axis limits and aspect ratio
    plt.xlim([np.min(Rwall), np.max(Rwall)])
    plt.ylim([np.min(Zwall), np.max(Zwall)])
    plt.gca().set_aspect('equal', adjustable='box')  # Ensuring equal aspect ratio

    # Adding labels and title with increased font size
    # plt.xlabel('Major Radius (R)', fontsize=18)
    # plt.ylabel('Vertical Position (Z)', fontsize=18)
    
    plt.xlabel('R [m]', fontsize=14)
    plt.ylabel('Z [m]', fontsize=14)
    # plt.title(r'$\Gamma_{O^{8+}}$ [m$^{-2}$s$^{-1}$]', fontsize=18)
    # plt.title(r'$n_{O^{8+}}$ [m$^{-3}$]', fontsize=18)
    # plt.title(r'$T_{O^{8+}}$ [eV]', fontsize=18) 
    
    # plt.title(r'$\Gamma_{D^{+}}$[m$^{-2}$s$^{-1}$]', fontsize=18)
    # plt.title(r'$n_{D^{+}}$[m$^{-3}$]', fontsize=18)
    # plt.title(r'$T_{D^{+}}$ [eV]', fontsize=18) 
    
    # plt.title(r'$\Gamma_{e^{-}}$ [m$^{-2}$s$^{-1}$]', fontsize=18)
    # plt.title(r'$n_{e^{-}}$ [m$^{-3}$]', fontsize=18)
    # plt.title(r'$T_{e^{-}}$ [eV]', fontsize=18) 
    
    plt.title(r'$B_r$ [T]', fontsize=18)
    # plt.title(r'$B_z$ [T]', fontsize=18)
    # plt.title(r'$B_\phi$ [T]', fontsize=18) # Setting tick font size
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)


    # Show legend and plot
    plt.legend()
    plt.show()
    # exit()

    # Write the data to an HDF5 file
    try:
        with h5py.File(openedge_file, 'w') as file:
            file.create_dataset('solps_like/r', data=r)
            file.create_dataset('solps_like/z', data=z)
            file.create_dataset('n_e/temp', data=temp_e_grid)
            file.create_dataset('n_i/temp', data=temp_i_grid)
            file.create_dataset('n_e/dens', data=dens_e_grid)
            file.create_dataset('n_i/dens', data=dens_i_grid)
            file.create_dataset('n_e/parr_flow', data=parr_flow_e_grid)
            file.create_dataset('n_i/parr_flow', data=parr_flow_i_grid)
            for name, grids in species_grids.items():
                file.create_dataset(f'{name}/temp', data=grids['temp'])
                file.create_dataset(f'{name}/dens', data=grids['dens'])
                file.create_dataset(f'{name}/parr_flow', data=grids['parr_flow'])
            file.create_dataset('bfield/b_r', data=b_r_grid)
            file.create_dataset('bfield/b_z', data=b_z_grid)
            file.create_dataset('bfield/b_phi', data=b_phi_grid)
            file.create_dataset('n_e/grad_te_parallel', data=grad_te_parallel)
            file.create_dataset('n_i/grad_ti_parallel', data=grad_ti_parallel)
            print(f"Data written to {openedge_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


# create plasma file for SPARTA-PMI
cases = ['1p5MW']
for case in cases:
    base_dir = f'/Users/78k/Desktop/soledge_postprocess/sim_O_4mw_main/run_dir'
    plasma_dir = f'/Users/78k/Desktop/soledge_postprocess/sim_O_4mw_main/run_dir'
    ref_file = os.path.join(base_dir, 'refParam_raptorX.h5')
    mesh_file = os.path.join(base_dir, 'meshEIRENE.h5')
    s3x_file = os.path.join(base_dir, 'mesh.h5')
    data_file = os.path.join(base_dir, 'plasmaFinal.h5')
    bfield_file = os.path.join(base_dir, 'mesh_raptorX.h5')
    wall_file = os.path.join('/Users/78k/Desktop/soledge_postprocess/soledge2GITR', 'new_wall.txt')
    openedge_file = os.path.join(plasma_dir, f'plasma_background_{case}_o1_to_o8.h5')

    interpolate_and_save_plasma_field(ref_file, mesh_file, bfield_file, data_file, wall_file, openedge_file, s3x_file)
