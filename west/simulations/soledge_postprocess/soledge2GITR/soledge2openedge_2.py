import os
import h5py
from shapely.geometry import Point, Polygon
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.interpolate import griddata
from matplotlib.path import Path
from read_wall import Surface


def interpolate_and_save_plasma_field(ref_file, mesh_file, bfield_file, data_file, wall_file, openedge_file, s3x_file,  nR=1000, nZ=200):
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
    # Load wall information  geometry
    path = '../'
    wall = Surface(wall_file, "2D")
    domain = wall.polygon
    Rwall, Zwall = domain.exterior.xy

    # Load mesh data
    with h5py.File(mesh_file, 'r') as meshEIRENE:
        tri_knots = meshEIRENE['/triangles/tri_knots'][...] - 1  # Convert to 0-based indexing
        Rk = meshEIRENE['/knots/R'][...] / 100
        Zk = meshEIRENE['/knots/Z'][...] / 100


    # set a mask on the core boundary: soledge inner boundary
    
    with h5py.File(s3x_file, 'r') as meshS3X:
        psi = meshS3X['config/psi'][...]
        psisep1=meshS3X['/config/psisep1'][...]
        psisep2=meshS3X['/config/psisep2'][...]
        psicore=meshS3X['/config/psicore'][...]
        r2D=meshS3X['/config/r'][...]
        z2D=meshS3X['/config/z'][...]
        
    print(psi.shape)

#
    # Load field data
    with h5py.File(data_file, 'r') as data:
        temp_e = data['/triangles/spec0/T'][...].flatten() * T0
        dens_e = data['/triangles/spec0/n'][...].flatten() * n0
        parr_flow_e = data['/triangles/spec0/G'][...].flatten() * c0 *n0
        
        temp_i = data['/triangles/spec1/T'][...].flatten() * T0
        dens_i = data['/triangles/spec1/n'][...].flatten() * n0
        parr_flow_i = data['/triangles/spec1/G'][...].flatten() * c0 *n0
        
        temp_o_8 = data['/triangles/spec9/T'][...].flatten() * T0
        dens_o_8 = data['/triangles/spec9/n'][...].flatten() * n0
        parr_flow_o8 = data['/triangles/spec9/G'][...].flatten() * c0 *n0


        # Load bfield bfield_file
    with h5py.File(bfield_file, 'r') as data:
        b_r = data['/triangles/Br'][...].flatten() * B0
        b_z = data['/triangles/Bz'][...].flatten() * B0
        b_phi = data['/triangles/Bphi'][...].flatten() * B0


    # Combine R and Z coordinates for Delaunay triangulation
    points = np.vstack((Rk, Zk)).T
    
    # Create a list of points for each triangle
    triangle_points = np.array([points[tri] for tri in tri_knots.T])
    
    # Calculate the centroid of each triangle to represent it
    centroids = triangle_points.mean(axis=1)
    
    # Given points for interpolation
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

    psi_core_threshold = psicore[0]  # Adjust this index or method to get the scalar value
    mask_inside_core = psi_grid < psi_core_threshold

    # Interpolate fields
    temp_e_grid = interpolate_field(temp_e)
    dens_e_grid = interpolate_field(dens_e)
    parr_flow_e_grid = interpolate_field(parr_flow_e)
    temp_i_grid = interpolate_field(temp_i)
    dens_i_grid = interpolate_field(dens_i)
    parr_flow_i_grid = interpolate_field(parr_flow_i)
    temp_o_8_grid = interpolate_field(temp_o_8)
    dens_o_8_grid = interpolate_field(dens_o_8)
#    psi_grid  = interpolate_field(psi)
    parr_flow_o8_grid = interpolate_field(parr_flow_o8)
    b_r_grid = interpolate_field(b_r)
    b_z_grid = interpolate_field(b_z)
    b_phi_grid = interpolate_field(b_phi)

    # Create a Path object for the wall polygon
    # (Assuming Rwall and Zwall are defined elsewhere in the script)
    wall_path = Path(np.vstack((Rwall, Zwall)).T)
    
    # Create a mask for points outside the wall
    mask_outside_wall = ~wall_path.contains_points(grid_points)

    # Apply mask to interpolated fields to set values to zero inside the core
    temp_e_grid[mask_inside_core] = 0
    dens_e_grid[mask_inside_core] = 0
    parr_flow_e_grid[mask_inside_core] = 0
    temp_i_grid[mask_inside_core] = 0
    dens_i_grid[mask_inside_core] = 0
    parr_flow_i_grid[mask_inside_core] = 0
    temp_o_8_grid[mask_inside_core] = 0
    dens_o_8_grid[mask_inside_core] = 0
    parr_flow_o8_grid[mask_inside_core] = 0
    b_r_grid[mask_inside_core] = 0
    b_z_grid[mask_inside_core] = 0
    b_phi_grid[mask_inside_core] = 0

        # Set interpolated values to zero for points outside the wall and inside the core
    fields = [temp_e_grid, dens_e_grid, temp_i_grid, dens_i_grid, temp_o_8_grid, dens_o_8_grid, parr_flow_e_grid, parr_flow_i_grid, parr_flow_o8_grid ]
    for field in fields:
        field[mask_outside_wall.reshape(nZ, nR)] = 0
#        field[mask_inside_core_reshaped] = 0
    
    Bmag = np.sqrt(b_r_grid*b_r_grid + b_z_grid * b_z_grid +  b_phi_grid * b_phi_grid)
    [grad_ti_z,grad_ti_r] = np.gradient(np.array(temp_e_grid),z[1]-z[0],r[1]-r[0])
    [grad_te_z,grad_te_r] = np.gradient(np.array(temp_i_grid),z[1]-z[0],r[1]-r[0])

    #    Small epsilon value to avoid division by zero
    epsilon = 1e-10
    # Compute parallel gradients with condition
    grad_te_parallel = np.where(Bmag < epsilon, 0, (grad_te_r * b_r_grid + grad_te_z * b_z_grid) / (Bmag + epsilon))
    grad_ti_parallel = np.where(Bmag < epsilon, 0, (grad_ti_r * b_r_grid + grad_ti_z * b_z_grid) / (Bmag + epsilon))
    
    min_temp_value = 1e0  # Adjust this as needed
    
    min_value = -9e20 # Set to the desired minimum
    max_value = 9e20 #np.max(temp_e_grid)  # Set to the desired maximum

    plt.contour(r2D, z2D, psi, linestyles='--', levels=[psicore[0]], label='Y')
    plt.contour(r2D, z2D, psi, linestyles='--',levels=[psisep1[0]],  label='Y')
    plt.contour(r2D, z2D, psi, linestyles='--',levels=[psisep2[0]],  label='Y')
    #plt.pcolormesh(grid_r, grid_z, parr_flow_i_grid, cmap='coolwarm', norm=LogNorm(vmin=min_temp_value)) # plasma, jet, inferno, viridis(default), coolwarm, RdBu, Blues
    #plt.pcolormesh(grid_r, grid_z, parr_flow_i_grid, cmap='coolwarm') # plasma, jet, inferno, viridis(default), coolwarm, RdBu, Blues
    plt.pcolormesh(grid_r, grid_z, parr_flow_o8_grid, cmap='coolwarm', vmin=min_value, vmax=max_value)
    plt.plot(Rwall, Zwall, 'k', lw=2.5)
    # Color bar with label
    cbar = plt.colorbar()
    cbar.set_label('Temperature (eV)', fontsize=18)
    cbar.ax.tick_params(labelsize=12)

    # Setting axis limits and aspect ratio
    plt.xlim([np.min(Rwall), np.max(Rwall)])
    plt.ylim([np.min(Zwall), np.max(Zwall)])
    plt.gca().set_aspect('equal', adjustable='box')  # Ensuring equal aspect ratio

    # Adding labels and title with increased font size
    plt.xlabel('Major Radius (R)', fontsize=18)
    plt.ylabel('Vertical Position (Z)', fontsize=18)
    #plt.title(r'$ Plasma Temperature  Distribution$', fontsize=18)
    
    # Setting tick font size
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)


    # Show legend and plot
    plt.legend()
    plt.show()



    # Write the data to an HDF5 file
    try:
        with h5py.File(openedge_file, 'w') as file:
            print("Creating HDF5 file...")
            file.create_dataset('solps_like/r', data=r)
            file.create_dataset('solps_like/z', data=z)
            file.create_dataset('n_e/temp', data=temp_e_grid)
            file.create_dataset('n_i/temp', data=temp_i_grid)
            file.create_dataset('o_8/temp', data=temp_o_8_grid)
            file.create_dataset('o_8/dens', data=dens_o_8_grid)
            file.create_dataset('n_e/dens', data=dens_e_grid)
            file.create_dataset('n_i/dens', data=dens_i_grid)
            file.create_dataset('n_e/parr_flow', data=parr_flow_e_grid)
            file.create_dataset('n_i/parr_flow', data=parr_flow_i_grid)
            file.create_dataset('o_8/parr_flow', data=parr_flow_o8_grid)
            file.create_dataset('bfield/b_r', data=b_r_grid)  # Uncomment if b_r_grid is defined
            file.create_dataset('bfield/b_z', data=b_z_grid)  # Uncomment if b_z_grid is defined
            file.create_dataset('bfield/b_phi', data=b_phi_grid)  # Uncomment if b_phi_grid is defined
            file.create_dataset('n_e/grad_te_parallel', data=grad_te_parallel)  # Uncomment if grad_te_parallel is defined
            file.create_dataset('n_i/grad_ti_parallel', data=grad_ti_parallel)  # Uncomment if grad_ti_parallel is defined
            print(f"Data written to {openedge_file}")
    
    except Exception as e:
        print(f"An error occurred: {e}")


# create plasma file for SPARTA-PMI

#cases = ['0p75MW', '1p5MW', '3MW']
cases = ['1p5MW']
for case in cases:
    base_dir = f'/home/78k/Desktop/sim_O_4mw_main/run_dir'
    plasma_dir = '/home/78k/Desktop/sim_O_4mw_main/run_dir'
    # Use os.path.join to create full file paths
    ref_file = os.path.join(base_dir, 'refParam_raptorX.h5')
    mesh_file = os.path.join(base_dir, 'meshEIRENE.h5')
    s3x_file = os.path.join(base_dir, 'mesh.h5')
    data_file = os.path.join(base_dir, 'plasmaFinal.h5')
    bfield_file = os.path.join(base_dir, 'mesh_raptorX.h5')
    wall_file = os.path.join('/home/78k/Desktop/soledge2GITR', 'new_wall.txt')
    openedge_file = os.path.join(plasma_dir, f'plasma_background_{case}_o1.h5')

    # Call the function with the file paths
    interpolate_and_save_plasma_field(ref_file, mesh_file, bfield_file, data_file, wall_file,openedge_file, s3x_file)
