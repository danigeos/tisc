#!/usr/bin/env python3
"""
TISC plotting script: 
Generates a 4-panel summary of Topography & Drainage, Erosion/Sedimentation Rates, 
Top Sediment Grain Size, and Salt (Halite & Gypsum) Precipitation.

Usage: python3 tisc.plot.topo+eros+grain+salt.py <projectname>
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {os.path.basename(__file__)} <projectname>")
        sys.exit(1)

    project = sys.argv[1]

    # File names
    f_st    = f"{project}.st"
    f_xyw   = f"{project}.xyw"
    f_grain = f"{project}.grainsize"
    f_salt  = f"{project}.thicksalt"

    # 1. Base Geometry & Erosion (.st file is always generated)
    if not os.path.exists(f_st):
        print(f"Error: Could not find base file '{f_st}'. Run TISC first!")
        sys.exit(1)

    print(f"Reading {f_st}...")
    st_data = np.loadtxt(f_st, comments='#')
    x_1d = st_data[:, 0]
    y_1d = st_data[:, 1]
    
    # Safely determine grid dimensions (TISC iterates y outer, x inner)
    Nx = 0
    for j in range(1, len(x_1d)):
        if y_1d[j] != y_1d[0]:
            Nx = j
            break
    if Nx == 0: 
        Nx = len(x_1d)
    Ny = len(x_1d) // Nx

    x = x_1d.reshape(Ny, Nx)
    y = y_1d.reshape(Ny, Nx)
    topo = st_data[:, 2].reshape(Ny, Nx)
    eros_rate = st_data[:, 4].reshape(Ny, Nx)

    # 2. Drainage Network (.xyw file)
    river_lines = []
    river_widths = []
    lake_grid = None
    precip_grid = None
    if os.path.exists(f_xyw):
        print(f"Reading {f_xyw}...")
        try:
            xyw_num = np.loadtxt(f_xyw, comments='#', usecols=(0, 1, 2, 6, 7, 9))
            x_from, y_from, water_flux, x_to, y_to, precip = xyw_num.T
        except Exception:
            # Fallback if the precipitation column is missing
            xyw_num = np.loadtxt(f_xyw, comments='#', usecols=(0, 1, 2, 6, 7))
            x_from, y_from, water_flux, x_to, y_to = xyw_num.T
            precip = np.zeros_like(water_flux)
            
        xyw_str = np.loadtxt(f_xyw, comments='#', usecols=(4,), dtype=str)
        
        # Lakes (Types 'L' for Lake, 'E' for Exit/Outlet, 'S' for Sea)
        lake_mask = (xyw_str == 'L') | (xyw_str == 'E') | (xyw_str == 'S')
        
        # Create a 2D bitmap for lakes based on the linear mask
        lake_grid = np.zeros(Ny * Nx)
        precip_grid = np.zeros(Ny * Nx)
        
        min_len = min(len(lake_grid), len(lake_mask))
        lake_grid[:min_len] = lake_mask[:min_len]
        lake_grid = lake_grid.reshape(Ny, Nx)
        
        precip_grid[:min_len] = precip[:min_len]
        precip_grid = precip_grid.reshape(Ny, Nx)
        
        # Rivers (draw segments for cells with significant water, exclude internal lake segments)
        r_mask = (water_flux > 0.1) & (xyw_str != 'L') & (xyw_str != 'S')
        river_lines = np.stack([np.column_stack([x_from[r_mask], y_from[r_mask]]), 
                                np.column_stack([x_to[r_mask], y_to[r_mask]])], axis=1)
        river_widths = np.log10(water_flux[r_mask] + 1) * 0.5
    else:
        print(f"Notice: {f_xyw} missing (hydro_model=0?). Drainage will be empty.")

    # 3. Top Sediment Grain Size (.grainsize file)
    grain = np.zeros((Ny, Nx))
    if os.path.exists(f_grain):
        print(f"Reading {f_grain}...")
        grain_data = np.loadtxt(f_grain, comments='#')
        # The last column represents the most recently deposited sedimentary block
        if grain_data.shape[1] > 2:
            grain = grain_data[:, -1].reshape(Ny, Nx)
    else:
        print(f"Notice: {f_grain} missing. Grain size will be 0.")

    # 4. Salt Thicknesses (.thicksalt file)
    gypsum = np.zeros((Ny, Nx))
    halite = np.zeros((Ny, Nx))
    if os.path.exists(f_salt):
        print(f"Reading {f_salt}...")
        salt_data = np.loadtxt(f_salt, comments='#')
        if salt_data.shape[1] >= 4:
            gypsum = salt_data[:, 2].reshape(Ny, Nx)
            halite = salt_data[:, 3].reshape(Ny, Nx)
    else:
        print(f"Notice: {f_salt} missing. Salt plots will be empty.")

    print("Generating plots...")
    # Create Figure
    fig, axs = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
    axs = axs.flatten()

    # ---- Panel 0: Topography + Drainage ----
    topo_vmin = np.nanmin(topo)
    topo_vmax = np.nanmax(topo)
    if topo_vmin >= 0:
        topo_vmin = -1.0
    if topo_vmax <= 0:
        topo_vmax = 1.0

    bath_cmap = plt.get_cmap("Blues_r")(np.linspace(0.2, 1.0, 128))
    land_cmap = plt.get_cmap("terrain")(np.linspace(0.25, 1.0, 128))
    bathy_land_cmap = mcolors.ListedColormap(np.vstack((bath_cmap, land_cmap)))
    topo_norm = mcolors.TwoSlopeNorm(vmin=topo_vmin, vcenter=0.0, vmax=topo_vmax)
    
    im0 = axs[0].pcolormesh(x, y, topo, cmap=bathy_land_cmap, norm=topo_norm, shading='nearest')
    div0 = make_axes_locatable(axs[0])
    cax0 = div0.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im0, cax=cax0, label='Topography (m)')
    
    # Overlay rain shade (top 10%)
    if precip_grid is not None and np.nanmax(precip_grid) > np.nanmin(precip_grid):
        p90 = np.nanpercentile(precip_grid, 90)
        pmax = np.nanmax(precip_grid)
        if pmax > p90:
            # Use contourf to create a smooth, gradient cloud
            levels = np.linspace(p90, pmax, 15)
            axs[0].contourf(x, y, precip_grid, levels=levels, cmap='Blues', alpha=0.6, extend='max')
            # Draw a contour line around the edge for better visibility
            axs[0].contour(x, y, precip_grid, levels=[p90], colors='dodgerblue', linewidths=1.2, alpha=0.8)

    # Overlay lakes as a bitmap
    if lake_grid is not None and np.any(lake_grid):
        lake_masked = np.ma.masked_where(lake_grid == 0, lake_grid)
        axs[0].pcolormesh(x, y, lake_masked, cmap=mcolors.ListedColormap(['#2244C8']), shading='nearest', alpha=0.8)
        
    # Overlay river network
    if len(river_lines) > 0:
        lc = LineCollection(river_lines, linewidths=river_widths, colors='blue', alpha=0.7)
        axs[0].add_collection(lc)
    axs[0].set_title('Topography, Drainage & High Rain (Top 10%)')

    # ---- Panel 1: Erosion / Sedimentation Rate ----
    # Use a centered norm to force white at 0 (Red=Erosion, Blue=Sedimentation)
    max_e = np.max(np.abs(eros_rate))
    if max_e == 0: max_e = 1
    im1 = axs[1].pcolormesh(x, y, eros_rate, cmap='RdBu_r', 
                            norm=mcolors.CenteredNorm(vcenter=0), shading='nearest')
    div1 = make_axes_locatable(axs[1])
    cax1 = div1.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im1, cax=cax1, label='Erosion (+) / Sedim (-) Rate (m/My)')
    axs[1].set_title('Surface Transport Rates')

    # ---- Panel 2: Grain Size ----
    im2 = axs[2].pcolormesh(x, y, grain, cmap='YlOrBr', shading='nearest')
    div2 = make_axes_locatable(axs[2])
    cax2 = div2.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im2, cax=cax2, label='Surface Grain Size (m)')
    axs[2].set_title('Top Sediment Grain Size')

    # ---- Panel 3: Evaporites (Halite + Gypsum) ----
    # Plot topography as a grayscale background for context
    axs[3].pcolormesh(x, y, topo, cmap='Greys', shading='nearest', alpha=0.4)
    
    T_salt = halite + gypsum
    mask_salt = T_salt > 1e-3
    rgba = np.zeros((Ny, Nx, 4)) # Transparent background
    
    # Setup custom axes for the 2D Colorbar Legend
    divider = make_axes_locatable(axs[3])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    
    if np.any(mask_salt):
        T_max = np.percentile(T_salt[mask_salt], 98)
        if T_max < 1e-3: T_max = np.max(T_salt)
        
        T_norm = np.clip(T_salt[mask_salt] / T_max, 0, 1)
        
        Yellow = np.array([1.0, 0.85, 0.0]) # Halite
        Purple = np.array([0.6, 0.1, 0.8])  # Gypsum
        White  = np.array([1.0, 1.0, 1.0])
        
        f_H = halite[mask_salt] / T_salt[mask_salt]
        f_G = gypsum[mask_salt] / T_salt[mask_salt]
        base_color = f_H[:, None] * Yellow + f_G[:, None] * Purple
        
        # Saturation increases quickly, Lightness decreases to make thicker = darker
        S = np.clip(T_norm * 2.5, 0, 1)
        V = np.clip(1.0 - T_norm * 0.6, 0, 1)
        
        rgba[mask_salt, :3] = (White * (1 - S[:, None]) + base_color * S[:, None]) * V[:, None]
        rgba[mask_salt, 3] = 1.0 # Opaque where salt exists
        axs[3].pcolormesh(x, y, rgba, shading='nearest')
        
        # Render the 2D Colorbar Legend
        leg_H, leg_T = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
        leg_S = np.clip(leg_T * 2.5, 0, 1)
        leg_V = np.clip(1.0 - leg_T * 0.6, 0, 1)
        leg_base = leg_H[..., None] * Yellow + (1 - leg_H)[..., None] * Purple
        leg_rgb = (White * (1 - leg_S[..., None]) + leg_base * leg_S[..., None]) * leg_V[..., None]
        
        cax.imshow(leg_rgb, origin='lower', aspect='auto', extent=[0, 1, 0, T_max])
        cax.set_xlabel('Fraction of Halite', fontsize=8)
        cax.set_ylabel('Total Thickness (m)', fontsize=8)
        cax.set_xticks([0, 0.5, 1])
        cax.set_xticklabels(['0\n(Gyp)', '0.5', '1\n(Hal)'])
        cax.tick_params(labelsize=8, pad=2)
        cax.yaxis.set_label_position("right")
        cax.yaxis.tick_right()
    else:
        cax.axis('off')
        
    axs[3].set_title('Evaporite Deposition')

    # General formatting
    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')

    # ---- Add Time Stamp ----
    t_cur, dt_val = None, 1.0
    prm_file = f"{project}.PRM"
    if os.path.exists(prm_file):
        with open(prm_file, 'r', errors='ignore') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == 'dt':
                    try:
                        dt_val = float(parts[1])
                    except ValueError:
                        pass

    for f_name in [f"{project}.hrz", f_grain, f_salt]:
        if os.path.exists(f_name):
            with open(f_name, 'r', errors='ignore') as f:
                header = f.readline()
                match = re.search(r'\(t\s*=\s*([0-9\.\+\-eE]+)', header)
                if match:
                    t_cur = float(match.group(1))
                    break
    
    if t_cur is not None:
        decimals = max(0, int(np.ceil(-np.log10(dt_val)))) if dt_val > 0 else 2
        fig.suptitle(f"Time: {t_cur:.{decimals}f} My", fontsize=18, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96], w_pad=-15.)
    else:
        plt.tight_layout(w_pad=-15.)

    svg_file = f"{project}.svg"
    jpg_file = f"{project}.jpg"
    plt.savefig(svg_file, dpi=200, bbox_inches='tight', format='svg')
    plt.savefig(jpg_file, dpi=200, bbox_inches='tight', format='jpg')
    print(f"Success! Output saved to: {svg_file} and {jpg_file}")

if __name__ == '__main__':
    main()