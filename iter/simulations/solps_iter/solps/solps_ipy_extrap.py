#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.io import loadmat, savemat
from scipy.interpolate import RegularGridInterpolator, griddata

def main():
    #%% Load SOLPS + EFIT Data from .mat file
    data = loadmat('solps_iter_1540.mat')

    # Extract data
    rS  = data['rS'].squeeze()    
    zS  = data['zS'].squeeze()    
    neS = data['neS']             
    TeS = data['teS']             
    g   = data['g']               
    g_lim = g[0,0]['lim']         

    #%% Data Cleaning
    T_min = 10     
    neS[(neS <= 0) | (~np.isfinite(neS))] = np.nan
    TeS[(TeS <= 0) | (~np.isfinite(TeS))] = np.nan
    TeS[TeS < T_min] = T_min

    #%% Grid Setup
    Nr, Nz = len(rS), len(zS)
    R_min, R_max = np.min(g_lim[0,:]), np.max(g_lim[0,:])
    Z_min, Z_max = np.min(g_lim[1,:]), np.max(g_lim[1,:])
    num_points = 20 

    X, Y = np.meshgrid(np.linspace(R_min, R_max, num_points),
                       np.linspace(Z_min, Z_max, num_points))

    #%% Identify valid (R, Z) points for interpolation
    rkron = np.repeat(rS, Nz)
    zkron = np.tile(zS, Nr)

    neS_T = neS.T  
    TeS_T = TeS.T  

    valid_mask_ne = ~np.isnan(neS_T)
    valid_mask_Te = ~np.isnan(TeS_T)

    coords_ne = np.column_stack((rkron[valid_mask_ne.flatten()], zkron[valid_mask_ne.flatten()]))
    coords_Te = np.column_stack((rkron[valid_mask_Te.flatten()], zkron[valid_mask_Te.flatten()]))

    ne_valid = neS_T[valid_mask_ne]
    Te_valid = TeS_T[valid_mask_Te]

    #%% Midplane Density & Temperature Extrapolation (Polynomial Fit)
    mpfx = np.linspace(8.2, 8.25, 1000)  
    mpfy = np.zeros_like(mpfx)           

    ne_interp = RegularGridInterpolator((zS, rS), neS_T, bounds_error=False, fill_value=np.nan)
    Te_interp = RegularGridInterpolator((zS, rS), TeS_T, bounds_error=False, fill_value=np.nan)

    pts_midplane = np.column_stack((mpfy, mpfx))
    fitDensityAtMidplane = ne_interp(pts_midplane)
    fitTeAtMidplane = Te_interp(pts_midplane)

    p_ne = np.polyfit(mpfx, np.log(fitDensityAtMidplane), 1)
    p_Te = np.polyfit(mpfx, np.log(fitTeAtMidplane), 1)

    mpx = np.linspace(8.1, 9, 1000)
    pts_midplane_full = np.column_stack((np.zeros_like(mpx), mpx))
    densityAtMidplane = ne_interp(pts_midplane_full)
    densityAtMidplane = np.where(np.isnan(densityAtMidplane), 0, densityAtMidplane)
    TeAtMidplane = Te_interp(pts_midplane_full)
    TeAtMidplane = np.where(np.isnan(TeAtMidplane), T_min, TeAtMidplane)

    interpfn = np.clip((mpx - 8.2) / (8.25 - 8.2), 0, 1)
    extrapolatedne1d = interpfn * np.exp(p_ne[1] + mpx * p_ne[0]) + (1 - interpfn) * densityAtMidplane
    extrapolatedTe1d = np.maximum(interpfn * np.exp(p_Te[1] + mpx * p_Te[0]) + (1 - interpfn) * TeAtMidplane, T_min)

    #%% 2D Interpolation for Smooth Extrapolation
    method = "cubic"  # Can also use "linear" or "nearest"

    val_ne = griddata(points=coords_ne,
                      values=ne_valid,
                      xi=(X, Y),
                      method=method,
                      fill_value=np.nan)

    val_Te = griddata(points=coords_Te,
                      values=Te_valid,
                      xi=(X, Y),
                      method=method,
                      fill_value=np.nan)

    val_Te[val_Te < T_min] = T_min  

    #%% Visualization of Extrapolated Electron Density
    plt.figure(figsize=(6,8))
    plt.imshow(val_ne, extent=[X.min(), X.max(), Y.min(), Y.max()],
               origin='lower', aspect='auto',
               norm=colors.LogNorm(vmin=1e10, vmax=1e20))
    plt.colorbar(label='n_e')
    plt.title('Extrapolated Electron Density (n_e)')
    plt.plot(g_lim[0, :], g_lim[1, :], 'r')
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    plt.show()

    #%% Visualization of Extrapolated Electron Temperature
    plt.figure(figsize=(6,8))
    plt.imshow(val_Te, extent=[X.min(), X.max(), Y.min(), Y.max()],
               origin='lower', aspect='auto',
               norm=colors.LogNorm(vmin=T_min, vmax=np.nanmax(val_Te)))
    plt.colorbar(label='T_e')
    plt.title('Extrapolated Electron Temperature (Te)')
    plt.plot(g_lim[0, :], g_lim[1, :], 'r')
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    plt.show()

    #%% Save Extrapolated Data as .mat + CSV
    savemat('extrapolated_data.mat', {
        'rS': rS, 'zS': zS, 'neS': neS, 'TeS': TeS,
        'X': X, 'Y': Y, 'val_ne': val_ne, 'val_Te': val_Te,
        'mpx': mpx,
        'extrapolatedne1d': extrapolatedne1d,
        'extrapolatedTe1d': extrapolatedTe1d,
        'g': g
    })
    np.savetxt('extrapolatedR.csv', X, delimiter=',')
    np.savetxt('extrapolatedZ.csv', Y, delimiter=',')
    np.savetxt('extrapolatedne.csv', val_ne, delimiter=',')
    np.savetxt('extrapolatedTe.csv', val_Te, delimiter=',')

if __name__ == '__main__':
    main()
