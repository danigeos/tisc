#!/usr/bin/env python3
"""
 tisc.3Dplot.py — Robust 3D conical-perspective plotter for TISC outputs

 Features
 --------
 • VISIBILITY: Rivers and Lakes are forced ON TOP (zorder=10000+) to prevent hiding.
 • LAKES: Drawn as filled rectangles (River Blue) with Darker Blue outlines.
 • COLORS: Land is Green/Brown; Water is Blue (#1f77b4).
 • FIXED: HRZ water layer stripping (prevents Blue Blanket).
 • DYNAMIC AZIMUTH: Rotates from 30 deg to 150 deg based on simulation time.

 Usage
 -----
 ./tisc.3Dplot.py sinusoidal --ve 30
"""

import os
import re
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.colors import LightSource, Normalize, LinearSegmentedColormap, to_rgba
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ---------- tolerant numeric token ----------
_FLOAT = re.compile(r'^[+\-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+\-]?\d+)?$')

def _flt_tokens(parts):
    out = []
    for p in parts:
        if _FLOAT.match(p):
            out.append(float(p.replace('D', 'E').replace('d', 'e')))
    return out

# ---------- Metadata Extraction ----------

def get_time_info(project):
    prm_file = f"{project}.PRM"
    hrz_file = f"{project}.hrz"
    t_ini, t_fin, t_cur = 0.0, 10.0, 0.0
    
    if os.path.exists(prm_file):
        try:
            with open(prm_file, 'r', errors='ignore') as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 2: continue
                    if parts[0] == 'Timeini': t_ini = float(parts[1])
                    elif parts[0] == 'Timefinal': t_fin = float(parts[1])
        except Exception: pass

    if os.path.exists(hrz_file):
        try:
            with open(hrz_file, 'r', errors='ignore') as f:
                header = f.readline()
                match = re.search(r'\(t\s*=\s*([0-9\.\+\-eE]+)', header)
                if match: t_cur = float(match.group(1))
        except Exception: pass
            
    return t_ini, t_fin, t_cur

# ---------- Loaders ----------

def load_hrz_robust(path):
    if not os.path.exists(path): raise FileNotFoundError(path)
    
    # Check Header for "water"
    has_water_header = False
    with open(path, 'r', errors='ignore') as f:
        for _ in range(20):
            line = f.readline()
            if not line: break
            if line.startswith('#'):
                if 'water' in line.lower():
                    has_water_header = True
            else:
                break 
                
    # Read Data
    xs, ys, zrows, ncols = [], [], [], None
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'): continue
            nums = _flt_tokens(s.replace(',', ' ').split())
            if len(nums) < 3: continue
            if ncols is None: ncols = len(nums)
            elif len(nums) != ncols: nums = (nums + [np.nan] * (ncols - len(nums)))[:ncols]
            xs.append(nums[0]); ys.append(nums[1]); zrows.append(nums[2:])
            
    X, Y, Z = np.asarray(xs), np.asarray(ys), np.asarray(zrows)
    if X.size == 0: raise ValueError("Empty HRZ")
    
    # DROP WATER LAYER IF PRESENT
    if has_water_header and Z.shape[1] > 1:
        print("[Info] Dropping global water layer found in HRZ file.")
        Z = Z[:, :-1]

    xu, yu = np.unique(X), np.unique(Y)
    nx, ny = len(xu), len(yu)
    
    Xkm, Ykm = np.meshgrid(xu, yu)
    Zlist = []
    for j in range(Z.shape[1]):
        G = np.full((ny, nx), np.nan)
        for xv, yv, zv in zip(X, Y, Z[:, j]):
            xi = np.searchsorted(xu, xv)
            yi = np.searchsorted(yu, yv)
            if xi < nx and yi < ny:
                G[yi, xi] = zv
        Zlist.append(G)
    names = [f'z{i+1}' for i in range(Z.shape[1])]
    
    return Xkm, Ykm, Zlist, names

def load_st_topo(path, nodata=-32767.9):
    if not os.path.exists(path): raise FileNotFoundError(path)
    xs, ys, zs = [], [], []
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'): continue
            nums = _flt_tokens(s.replace(',', ' ').split())
            if len(nums) == 5:
                xs.append(nums[0]); ys.append(nums[1]); zs.append(nums[2])
    X, Y, Z = np.asarray(xs), np.asarray(ys), np.asarray(zs)
    Z = np.where(np.isclose(Z, nodata), np.nan, Z)
    xu, yu = np.unique(X), np.unique(Y)
    nx, ny = len(xu), len(yu)
    Xkm, Ykm = np.meshgrid(xu, yu)
    G = np.full((ny, nx), np.nan)
    for xv, yv, zv in zip(X, Y, Z): 
        xi = np.searchsorted(xu, xv)
        yi = np.searchsorted(yu, yv)
        if xi < nx and yi < ny:
            G[yi, xi] = zv
    return Xkm, Ykm, [G], ['topo']

def load_pfl(path):
    if not os.path.exists(path): return None
    pts = []
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'): continue
            nums = _flt_tokens(s.replace(',', ' ').split())
            if len(nums) >= 2: pts.append(nums[:2])
            if len(pts) == 2: break
    return np.array(pts) if len(pts) == 2 else None

def load_lak_squares(path):
    """
    Parses .lak file for distinct lakes and their elevations.
    Returns: [{'z': elevation, 'cells': [[x,y], ...]}, ...]
    """
    if not os.path.exists(path): return []
    
    lakes = []
    current_lake = None
    
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s: continue
            
            if s.startswith('>'):
                if current_lake and current_lake['cells']:
                    lakes.append(current_lake)
                
                current_lake = {'z': 0.0, 'cells': []}
                
                parts = s.split()
                found_z = False
                try:
                    for i, p in enumerate(parts):
                        if p == 'm' and i > 0:
                            current_lake['z'] = float(parts[i-1])
                            found_z = True
                            break
                    if not found_z:
                        match = re.search(r'([0-9\.]+)\s*m', s)
                        if match: current_lake['z'] = float(match.group(1))
                except: pass
            
            elif not s.startswith('#') and current_lake is not None:
                nums = _flt_tokens(s.replace(',', ' ').split())
                if len(nums) >= 2:
                    current_lake['cells'].append([nums[0], nums[1]])
    
    if current_lake and current_lake['cells']:
        lakes.append(current_lake)
    
    print(f"[Info] Found {len(lakes)} lakes in .lak file.")
    return lakes

def process_xyw(path, Xkm, Ykm):
    if not os.path.exists(path): return [], [], None
    rows_xyz, rows_w = [], []
    with open(path,'r',errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'): continue
            parts = s.replace(',', ' ').split()
            if len(parts) < 8: continue
            nums = _flt_tokens(parts)
            if len(nums) < 8: continue
            
            rows_xyz.append(np.array([[nums[0],nums[1],nums[4]],[nums[5],nums[6],nums[7]]], float))
            rows_w.append(float(nums[2]))
    return rows_xyz, rows_w, None

# ---------- Rendering ----------

def create_land_colormap():
    colors = [
        (0.0, "#228B22"), # Forest Green
        (0.3, "#32CD32"), # Lime Green
        (0.6, "#F4A460"), # Sandy Brown
        (1.0, "#FFFFFF")  # White
    ]
    return LinearSegmentedColormap.from_list('land_only', colors)

def bilinear(X, Y, Z, xp, yp):
    xv, yv = X[0, :], Y[:, 0]
    xp = np.clip(xp, xv[0], xv[-1]); yp = np.clip(yp, yv[0], yv[-1])
    ix = np.clip(np.searchsorted(xv, xp) - 1, 0, len(xv) - 2)
    iy = np.clip(np.searchsorted(yv, yp) - 1, 0, len(yv) - 2)
    wx = (xp - xv[ix]) / np.maximum(xv[ix+1] - xv[ix], 1e-12)
    wy = (yp - yv[iy]) / np.maximum(yv[iy+1] - yv[iy], 1e-12)
    z0 = Z[iy, ix] * (1-wx) + Z[iy, ix+1] * wx
    z1 = Z[iy+1, ix] * (1-wx) + Z[iy+1, ix+1] * wx
    return z0 * (1-wy) + z1 * wy

def nice_ticks(vmin, vmax, n=5):
    if not np.isfinite([vmin, vmax]).all() or vmin == vmax: return np.array([vmin])
    return np.linspace(vmin, vmax, n)

def render(project, Xkm, Ykm, Zm_list, lake_list=[], pfl=None,
           azim=-60.0, elev=25.0, ve=10.0,
           xyw_geoms=None, xyw_linewidths=None,
           proj='persp', camdist=None, fov=None, time_val=None):
    
    Z_land = Zm_list[-1]
    Zothers_m = Zm_list[:-1]
    
    # Global limits
    all_z = [Z_land] + Zothers_m
    if lake_list:
        all_z.extend([lake['z'] for lake in lake_list])

    g_zmin = np.nanmin([np.nanmin(z) for z in all_z])
    g_zmax = np.nanmax([np.nanmax(z) for z in all_z])
    
    Zland_sc = (Z_land / 1000.0) * ve
    Zothers_sc = [(z / 1000.0) * ve for z in Zothers_m]
    zmin_p = (g_zmin / 1000.0) * ve
    zmax_p = (g_zmax / 1000.0) * ve

    fig = plt.figure(figsize=(14.0, 10.0), dpi=200)
    ax = fig.add_subplot(111, projection='3d')
    
    dx_grid = float(np.nanmax(Xkm) - np.nanmin(Xkm))
    dy_grid = float(np.nanmax(Ykm) - np.nanmin(Ykm))
    dz_grid = float(zmax_p - zmin_p)
    ax.set_box_aspect([dx_grid, dy_grid, dz_grid if dz_grid > 0 else 1.0])
    ax.set_zlim(zmin_p, zmax_p)
    
    if time_val is not None:
        ax.set_title(f"Time: {time_val:.2f} Myr", fontsize=16)

    if proj == 'ortho':
        ax.set_proj_type('ortho')
    else:
        dist_val = 10.0
        focal_val = 1.0
        if fov is not None:
            fov_val = float(fov)
            focal_val = 1.0 / np.tan(np.deg2rad(fov_val) / 2.0)
            if fov_val >= 130: dist_val = 2.5 
            elif fov_val >= 90: dist_val = 5.5
            elif fov_val >= 60: dist_val = 8.0
            if camdist is not None: dist_val = float(camdist)
        elif camdist is not None:
            dist_val = float(camdist)
        try: ax.set_proj_type('persp', focal_length=focal_val)
        except: ax.set_proj_type('persp')
        ax.dist = dist_val

    # 1. Plot LAND
    ls = LightSource(azdeg=315, altdeg=45)
    land_cmap = create_land_colormap()
    norm = Normalize(vmin=np.nanmin(Z_land), vmax=np.nanmax(Z_land))
    land_fc = ls.shade(Z_land, cmap=land_cmap, norm=norm, vert_exag=ve/1000.0, blend_mode='soft')
    
    ax.plot_surface(Xkm, Ykm, Zland_sc, rstride=1, cstride=1,
                    facecolors=land_fc, linewidth=0, antialiased=True, shade=False, zorder=1)

    # 2. Plot LAKES (Rectangles) - Forced ON TOP
    if lake_list:
        # horizontal size of each grid cell in km
        cell_dx = np.abs(np.mean(np.diff(Xkm[0, :])))
        cell_dy = np.abs(np.mean(np.diff(Ykm[:, 0])))
        hx, hy = cell_dx / 2.0, cell_dy / 2.0

        # same blue as rivers
        river_blue = '#1f77b4'
        edge_col   = '#104060'   # optional darker outline

        for lake in lake_list:
            z_lake_m = lake['z']

            # small vertical lift so the lake surface sits clearly on top
            z_lift = 0.005 * dz_grid
            z_lake_plot = (z_lake_m / 1000.0) * ve + z_lift

            for (cx, cy) in lake['cells']:
                # Build a tiny 2×2 surface patch for each lake cell
                Xsq = np.array([[cx - hx, cx + hx],
                                [cx - hx, cx + hx]])
                Ysq = np.array([[cy - hy, cy - hy],
                                [cy + hy, cy + hy]])
                Zsq = np.full_like(Xsq, z_lake_plot, dtype=float)

                # Filled blue rectangle on top of the terrain
                ax.plot_surface(
                    Xsq, Ysq, Zsq,
                    color=river_blue,
                    linewidth=0,
                    antialiased=False,
                    shade=False,
                    zorder=10001,
                )

                # Optional outline so you still see the lake borders
                x_sq = [cx - hx, cx + hx, cx + hx, cx - hx, cx - hx]
                y_sq = [cy - hy, cy - hy, cy + hy, cy + hy, cy - hy]
                z_sq = [z_lake_plot] * 5
                ax.plot(
                    x_sq, y_sq, z_sq,
                    color=river_blue,
                    linewidth=2.2,
                    zorder=10002,
                )

    # 3. Plot Subsurface
    for i, Zp in enumerate(Zothers_sc[::-1], start=1):
        ax.plot_surface(Xkm, Ykm, Zp, linewidth=0, antialiased=True, shade=True, 
                        color='#808080', alpha=0.55, zorder=0)

    # 4. Plot RIVERS (Z-Order 10003 -> On Top of Lakes)
    if xyw_geoms:
        for geom, lw in zip(xyw_geoms, xyw_linewidths):
            ax.plot(geom[:,0], geom[:,1], (geom[:,2]/1000.0)*ve, color='#1f77b4', linewidth=lw, zorder=10003)

    # 5. Profile
    if pfl is not None:
        (x1, y1), (x2, y2) = pfl
        t = np.linspace(0, 1, 700)
        xp = x1 + (x2-x1)*t; yp = y1 + (y2-y1)*t
        zline = (bilinear(Xkm, Ykm, Z_land, xp, yp) / 1000.0) * ve
        ax.plot(xp, yp, zline, color='k', linewidth=3.0, zorder=10004)
        step = max(1, len(xp)//150)
        for i in range(0, len(xp), step):
            ax.plot([xp[i], xp[i]], [yp[i], yp[i]], [zmin_p, zline[i]], color='k', linewidth=1.0, zorder=999)

    ax.set_xlabel('X (km)'); ax.set_ylabel('Y (km)')
    ax.set_zlabel('Elevation (m)')
    z_ticks = nice_ticks(g_zmin, g_zmax, n=5)
    ax.set_zticks((z_ticks/1000.0)*ve)
    ax.set_zticklabels([f'{int(v)}' for v in z_ticks])
    ax.view_init(elev=elev, azim=azim)

    out = f'{project}.jpg'
    plt.savefig(out, bbox_inches='tight', dpi=300, pad_inches=0.25)
    return out

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description='TISC 3D Plotter: Visualizes HRZ/ST topography, XYW rivers, and Lakes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument('project', help='Project base name.')
    ap.add_argument('--azim', type=float, default=None, 
                    help='Camera Azimuth in degrees. If None, calculated dynamically from Timeini(30) to Timefinal(150).')
    ap.add_argument('--elev', type=float, default=25.0, 
                    help='Camera Elevation in degrees.')
    ap.add_argument('--ve', type=float, default=10.0, 
                    help='Vertical Exaggeration.')
    ap.add_argument('--proj', choices=['persp','ortho'], default='persp', 
                    help='Projection Type.')
    ap.add_argument('--camdist', type=float, default=None, 
                    help='Camera Distance (Legacy).')
    ap.add_argument('--fov', type=float, default=150.0, 
                    help='Field of View in degrees. Controls the "Wide Angle" effect.')
    ap.add_argument('--n-cells-river-plot', type=float, default=0.5, 
                    help='River threshold multiplier.')

    args = ap.parse_args()
    proj = args.project
    
    hrz = f'{proj}.hrz'
    st  = f'{proj}.st'
    pfl = f'{proj}.pfl'
    xyw = f'{proj}.xyw'
    lak = f'{proj}.lak'
    
    t_ini, t_fin, t_cur = get_time_info(proj)
    
    if args.azim is not None:
        azim_val = args.azim
    else:
        if t_fin > t_ini:
            progress = (t_cur - t_ini) / (t_fin - t_ini)
            progress = max(0.0, min(1.0, progress))
            azim_val = 30.0 + (progress * (150.0 - 30.0))
        else:
            azim_val = 30.0
            
    print(f"[Info] Time: {t_cur:.2f} Myr (Range: {t_ini}-{t_fin}). Calculated Azim: {azim_val:.1f}")

    if os.path.exists(hrz):
        Xkm, Ykm, Zm_list, names = load_hrz_robust(hrz)
    elif os.path.exists(st):
        Xkm, Ykm, Zm_list, names = load_st_topo(st)
    else:
        sys.exit(f"Error: {hrz} or {st} not found.")

    dx = np.abs(np.mean(np.diff(Xkm[0, :])))
    dy = np.abs(np.mean(np.diff(Ykm[:, 0])))
    max_segment_len_km = np.sqrt(dx**2 + dy**2) * 1.5

    raw_segs, raw_ws, _ = process_xyw(xyw, Xkm, Ykm)
    
    # Load Lakes (Explicit Squares)
    lake_list = load_lak_squares(lak)

    draw_geoms, draw_widths = [], []
    if raw_segs:
        valid_ws = [w for w in raw_ws if w > 1e-6]
        threshold = 0.0
        if valid_ws:
            hist, bins = np.histogram(valid_ws, bins=50)
            peak_idx = np.argmax(hist)
            mode_q = (bins[peak_idx] + bins[peak_idx+1]) / 2.0
            threshold = mode_q * args.n_cells_river_plot

        kept_ws = []
        for seg, w in zip(raw_segs, raw_ws):
            dist_2d = np.sqrt((seg[1,0]-seg[0,0])**2 + (seg[1,1]-seg[0,1])**2)
            if dist_2d > max_segment_len_km: continue
            if w < threshold: continue
            draw_geoms.append(seg)
            kept_ws.append(w)
        
        if kept_ws:
            kept_ws = np.array(kept_ws)
            w_min, w_max = np.min(kept_ws), np.max(kept_ws)
            for w in kept_ws:
                width = 0.8
                if w_max > w_min and w > 0:
                    try:
                        num = np.log10(w) - np.log10(w_min)
                        den = np.log10(w_max) - np.log10(w_min)
                        width = 0.8 + (num / den) * 2.7
                    except:
                        pass
                draw_widths.append(width)

    jpg = render(proj, Xkm, Ykm, Zm_list, lake_list=lake_list, pfl=load_pfl(pfl),
                 azim=azim_val, elev=args.elev, ve=args.ve,
                 xyw_geoms=draw_geoms, xyw_linewidths=draw_widths,
                 proj=args.proj, camdist=args.camdist, fov=args.fov, time_val=t_cur)
    print('WROTE:', jpg)

if __name__ == '__main__':
    main()