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
import glob
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.colors import LightSource, Normalize, LinearSegmentedColormap, to_rgba
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.path import Path as MplPath
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

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
    zmin_prm, zmax_prm = None, None
    
    if os.path.exists(prm_file):
        try:
            with open(prm_file, 'r', errors='ignore') as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 2: continue
                    if parts[0] == 'Timeini': t_ini = float(parts[1])
                    elif parts[0] == 'Timefinal': t_fin = float(parts[1])
                    elif parts[0] == 'zmin': zmin_prm = float(parts[1])
                    elif parts[0] == 'zmax': zmax_prm = float(parts[1])
        except Exception: pass

    if os.path.exists(hrz_file):
        try:
            with open(hrz_file, 'r', errors='ignore') as f:
                header = f.readline()
                match = re.search(r'\(t\s*=\s*([0-9\.\+\-eE]+)', header)
                if match: t_cur = float(match.group(1))
        except Exception: pass
            
    return t_ini, t_fin, t_cur, zmin_prm, zmax_prm

# ---------- Loaders ----------

def load_hrz_robust(path):
    if not os.path.exists(path): raise FileNotFoundError(path)
    
    # Check Header for "water" and count density tokens
    has_water_header = False
    dens_tokens = []
    with open(path, 'r', errors='ignore') as f:
        for _ in range(20):
            line = f.readline()
            if not line: break
            if line.startswith('#'):
                low = line.lower()
                if 'water' in low:
                    has_water_header = True
                if 'dens' in low:
                    dens_tokens = _flt_tokens(line.replace(',', ' ').split())
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
    
    # Only drop a trailing water column if there is actually an extra column beyond the rock/sediment densities
    if has_water_header and Z.shape[1] > 1 and (not dens_tokens or Z.shape[1] > len(dens_tokens)):
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
    if not os.path.exists(path):
        return None
    densities = []
    ages = []
    rows = []
    ncols = None
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                low = s.lower()
                if 'densities' in low:
                    densities = _flt_tokens(s.replace(',', ' ').split())
                elif 'ages' in low:
                    ages = _flt_tokens(s.replace(',', ' ').split())
                continue
            nums = _flt_tokens(s.replace(',', ' ').split())
            if len(nums) < 4:
                continue
            if ncols is None:
                ncols = len(nums)
            elif len(nums) != ncols:
                nums = (nums + [np.nan] * (ncols - len(nums)))[:ncols]
            rows.append(nums)
    if len(rows) < 2:
        return None
    data = np.asarray(rows, dtype=float)
    x = data[:, 0]
    y = data[:, 1]
    l = data[:, 2]
    z = data[:, 3:]

    # Drop water/ice layer at the top if present (density <= 1100), matching GMT col < numcols
    while z.shape[1] > 1 and len(densities) >= z.shape[1] and densities[z.shape[1] - 1] <= 1100:
        z = z[:, :-1]
        densities = densities[:z.shape[1]]

    return {
        'x': x,
        'y': y,
        'l': l,
        'z': z,
        'densities': densities,
        'ages': ages,
    }

def load_lak_squares(path):
    """
    Parses .lak file for distinct continental lakes and their elevations.
    Skips sea cells ('S') and tiny discharge puddles, matching tisc.topo3D.gmt.job.
    Returns: [{'z': elevation, 'endorheic': bool, 'cells': [[x,y], ...]}, ...]
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
                
                current_lake = {'z': 0.0, 'disch': None, 'outl': 1, 'cells': []}
                
                parts = s.split()
                try:
                    for i, p in enumerate(parts):
                        if p == 'm' and i > 0:
                            current_lake['z'] = float(parts[i-1])
                        elif p.startswith('outl') and i > 0:
                            current_lake['outl'] = int(float(parts[i-1]))
                        elif p.startswith('inp') and i > 1:
                            current_lake['disch'] = float(parts[i-2])
                except Exception:
                    pass
            
            elif not s.startswith('#') and current_lake is not None:
                parts = s.replace(',', ' ').split()
                if len(parts) >= 3 and parts[2] == 'S':
                    continue
                if current_lake['disch'] is not None and current_lake['disch'] <= 2.5:
                    continue
                nums = _flt_tokens(parts[:2])
                if len(nums) >= 2:
                    current_lake['cells'].append([nums[0], nums[1]])
    
    if current_lake and current_lake['cells']:
        lakes.append(current_lake)
    
    print(f"[Info] Found {len(lakes)} continental lakes in .lak file.")
    return lakes

def process_xyw(path, Xkm, Ykm):
    if not os.path.exists(path): return [], [], 0.0
    rows_xyz, rows_w = [], []
    sea_level = 0.0
    with open(path,'r',errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s: continue
            if s.startswith('#'):
                m = re.search(r'sea_level:\s*([0-9\.\+\-eE]+)', s)
                if m:
                    try: sea_level = float(m.group(1))
                    except Exception: pass
                continue
            parts = s.replace(',', ' ').split()
            if len(parts) < 9: continue
            # Skip sea ('S') and lake ('L') nodes, matching tisc.topo3D.gmt.job ($5 != "S" && $5 != "L")
            dr_type = parts[4]
            if dr_type in ('S', 'L'):
                continue
            nums = _flt_tokens(parts[:4] + parts[5:9])
            if len(nums) < 8: continue
            
            rows_xyz.append(np.array([[nums[0],nums[1],nums[4]],[nums[5],nums[6],nums[7]]], float))
            rows_w.append(float(nums[2]))
    return rows_xyz, rows_w, sea_level

# ---------- UNIT File Parser & Tectonic Arrows ----------

def parse_unit_file(path, t_ini=0.0):
    t_load = t_ini
    t_stop = 9999.0
    vx, vy = 0.0, 0.0
    z_default = -1.0e6
    mode_interp = 4
    polys = []
    cur_poly = None
    in_polys = False

    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.split('#')[0].strip()
            if not s:
                continue
            parts = s.split()
            if not in_polys and not _FLOAT.match(parts[0]):
                key = parts[0]
                if len(parts) >= 2 and _FLOAT.match(parts[1]):
                    val = float(parts[1].replace('D', 'E').replace('d', 'e'))
                    if key in ('time', 'time_load'):
                        t_load = val
                    elif key == 'time_stop':
                        t_stop = val
                    elif key == 'vel_x':
                        vx = val
                    elif key == 'vel_y':
                        vy = val
                    elif key == 'z_default':
                        z_default = val
                    elif key == 'mode_interp':
                        mode_interp = int(val)
                continue
            if _FLOAT.match(parts[0]):
                if len(parts) == 1 or (len(parts) >= 2 and not _FLOAT.match(parts[1])):
                    in_polys = True
                    z_val = float(parts[0].replace('D', 'E').replace('d', 'e'))
                    ptype = parts[1][0].lower() if len(parts) >= 2 else 'n'
                    cur_poly = {'z': z_val, 'type': ptype, 'pts': []}
                    polys.append(cur_poly)
                elif len(parts) >= 2 and _FLOAT.match(parts[1]):
                    in_polys = True
                    if cur_poly is not None:
                        px = float(parts[0].replace('D', 'E').replace('d', 'e'))
                        py = float(parts[1].replace('D', 'E').replace('d', 'e'))
                        cur_poly['pts'].append((px, py))

    # Convert polygon coordinates from meters to km if given in meters (matching libreria.c)
    all_coords = [c for p in polys for pt in p['pts'] for c in pt]
    scale = 1.0 / 1000.0 if (all_coords and max(abs(c) for c in all_coords) > 2000.0) else 1.0
    for p in polys:
        p['pts'] = np.array([[pt[0] * scale, pt[1] * scale] for pt in p['pts']], dtype=float)

    return {
        'path': path,
        't_load': t_load,
        't_stop': t_stop,
        'vx': vx,
        'vy': vy,
        'z_default': z_default,
        'mode_interp': mode_interp,
        'polys': [p for p in polys if len(p['pts']) >= 3],
    }

def _dist_point_to_poly_edges(xp, yp, pts):
    n = len(pts)
    min_d = 1.0e30
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        dx, dy = x2 - x1, y2 - y1
        l2 = dx * dx + dy * dy
        if l2 < 1e-12:
            d = float(np.hypot(xp - x1, yp - y1))
        else:
            t = max(0.0, min(1.0, ((xp - x1) * dx + (yp - y1) * dy) / l2))
            d = float(np.hypot(xp - (x1 + t * dx), yp - (y1 + t * dy)))
        if d < min_d:
            min_d = d
    return min_d

def eval_unit_polygons(unit, xp, yp):
    """Evaluates TISC mode-4 nested polygons at point (xp, yp) in km."""
    polys = unit['polys']
    z = unit['z_default']
    if not polys:
        return z
    path0 = MplPath(polys[0]['pts'])
    inside_prev = path0.contains_point((xp, yp), radius=-1e-5) or (
        _dist_point_to_poly_edges(xp, yp, polys[0]['pts']) < 1e-3
    )
    if inside_prev:
        d_prev = _dist_point_to_poly_edges(xp, yp, polys[0]['pts'])
        z = polys[0]['z'] if polys[0]['type'] == 'c' else (
            unit['z_default'] + (polys[0]['z'] - unit['z_default']) * min(d_prev, 1.0)
        )
    else:
        d_prev = 0.0

    for k in range(1, len(polys)):
        pk = polys[k]
        pathk = MplPath(pk['pts'])
        inside_k = pathk.contains_point((xp, yp), radius=-1e-5) or (
            _dist_point_to_poly_edges(xp, yp, pk['pts']) < 1e-3
        )
        dk = _dist_point_to_poly_edges(xp, yp, pk['pts'])
        if inside_k:
            z = pk['z']
        elif inside_prev:
            if polys[k - 1]['type'] != 'c':
                denom = d_prev + dk
                z = polys[k - 1]['z'] + (pk['z'] - polys[k - 1]['z']) * (
                    d_prev / denom if denom > 1e-9 else 0.0
                )
            break
        inside_prev = inside_k
        d_prev = dk
    return z

def get_unit_arrows(project, Xkm, Ykm, t_cur, t_ini=0.0):
    """
    Finds all <project>*.UNIT files active at t_cur with non-zero velocity
    and computes boundary (or interior) arrow anchors and velocity vectors.
    """
    if t_cur is None:
        return []
    pattern = f"{project}*.UNIT"
    candidates = glob.glob(pattern)
    # Sort naturally by unit index
    def _unit_sort_key(p):
        m = re.search(r'(\d+)\.UNIT$', p)
        return (0, int(m.group(1))) if m else (1, p)
    unit_files = sorted(candidates, key=_unit_sort_key)
    if not unit_files:
        return []

    x0, xf = float(np.nanmin(Xkm)), float(np.nanmax(Xkm))
    y0, yf = float(np.nanmin(Ykm)), float(np.nanmax(Ykm))
    Lx, Ly = xf - x0, yf - y0
    max_L = max(Lx, Ly)

    # Choose a clean round coordinate step along boundaries (~100 km for a 500-750 km domain)
    raw_step = max_L / 7.5
    mag = 10.0 ** np.floor(np.log10(max(raw_step, 1e-6)))
    for mult in (1.0, 2.0, 2.5, 5.0, 10.0):
        if mult * mag >= raw_step * 0.85:
            step = mult * mag
            break
    else:
        step = 10.0 * mag

    tx_min = np.ceil((x0 - 1e-6) / step) * step
    tx_max = np.floor((xf + 1e-6) / step) * step
    x_ticks = np.arange(tx_min, tx_max + 0.1 * step, step)

    ty_min = np.ceil((y0 - 1e-6) / step) * step
    ty_max = np.floor((yf + 1e-6) / step) * step
    y_ticks = np.arange(ty_min, ty_max + 0.1 * step, step)

    boundaries = [
        ('N', [(float(x), yf) for x in x_ticks], (0.0, 1.0)),
        ('S', [(float(x), y0) for x in x_ticks], (0.0, -1.0)),
        ('W', [(x0, float(y)) for y in y_ticks], (-1.0, 0.0)),
        ('E', [(xf, float(y)) for y in y_ticks], (1.0, 0.0)),
    ]

    arrows_by_pos = {}
    for uf in unit_files:
        u = parse_unit_file(uf, t_ini=t_ini)
        vx, vy = u['vx'], u['vy']
        speed = float(np.hypot(vx, vy))
        if speed <= 1e-6:
            continue
        if not (u['t_load'] - 1e-4 <= t_cur <= u['t_stop'] + 1e-4):
            continue
        if not u['polys']:
            continue

        zmax_u = max(p['z'] for p in u['polys'])
        if zmax_u <= 0:
            continue
        max_polys = [p for p in u['polys'] if np.isclose(p['z'], zmax_u)]
        is_const = any(p['type'] == 'c' for p in max_polys)
        z_thresh = (0.95 if is_const else 0.75) * zmax_u

        unit_arrows = []
        for bname, bpts, (nx_out, ny_out) in boundaries:
            # Skip boundaries parallel to the velocity vector
            v_norm_comp = abs(vx * nx_out + vy * ny_out) / speed
            if v_norm_comp < 0.20:
                continue
            inward = (vx * nx_out + vy * ny_out) <= 0.0

            for bx, by in bpts:
                z0_val = eval_unit_polygons(u, bx, by)
                if bname in ('N', 'S'):
                    z_minus = eval_unit_polygons(u, max(x0, bx - 0.20 * step), by)
                    z_plus  = eval_unit_polygons(u, min(xf, bx + 0.20 * step), by)
                else:
                    z_minus = eval_unit_polygons(u, bx, max(y0, by - 0.20 * step))
                    z_plus  = eval_unit_polygons(u, bx, min(yf, by + 0.20 * step))

                if z0_val >= z_thresh and z_minus >= z_thresh and z_plus >= z_thresh:
                    key = (bname, round(bx, 2), round(by, 2))
                    arrows_by_pos[key] = {
                        'x': bx,
                        'y': by,
                        'vx': vx,
                        'vy': vy,
                        'speed': speed,
                        'inward': inward,
                        'bname': bname,
                    }
                    unit_arrows.append(key)

        # Fallback for an active moving unit completely inside the domain
        if not unit_arrows and max_polys:
            pts_in = max_polys[-1]['pts']
            cx = float(np.clip(np.mean(pts_in[:, 0]), x0, xf))
            cy = float(np.clip(np.mean(pts_in[:, 1]), y0, yf))
            key = ('INT', round(cx, 2), round(cy, 2))
            arrows_by_pos[key] = {
                'x': cx,
                'y': cy,
                'vx': vx,
                'vy': vy,
                'speed': speed,
                'inward': False,
                'bname': 'INT',
            }

    return list(arrows_by_pos.values())

# ---------- Rendering ----------

def color_and_illuminate_topo(Xkm, Ykm, Z_land, sea_level=0.0):
    """
    Colors and illuminates Z_land using the TISC GMT bathymetry + continental palette
    (generate_topo_cpt.job) and directional arctan slope shading (grdgradient -A45).
    """
    z_rel = Z_land - sea_level
    stops_z = np.array([
        -5000.0, -3000.0, -1500.0, -800.0, -200.0, -0.01,
        0.0, 200.0, 600.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0
    ], dtype=float)
    stops_rgb = np.array([
        [130,  60, 170],
        [ 90,  90, 180],
        [ 90, 150, 190],
        [150, 215, 225],
        [210, 240, 240],
        [210, 240, 240],
        [100, 200, 100],
        [160, 255, 100],
        [255, 255, 120],
        [255, 220,  60],
        [210, 180, 110],
        [210, 150,  80],
        [190, 120,  50],
        [205, 205, 205],
        [245, 245, 245],
    ], dtype=float) / 255.0

    rgb = np.empty(Z_land.shape + (3,), dtype=float)
    z_flat = np.nan_to_num(z_rel, nan=0.0)
    for c in range(3):
        rgb[..., c] = np.interp(z_flat, stops_z, stops_rgb[:, c])

    # Compute directional gradient from NE (azimuth 45 deg, matching GMT grdgradient -A45)
    dx_m = float(np.abs(np.mean(np.diff(Xkm[0, :])))) * 1000.0
    dy_m = float(np.abs(np.mean(np.diff(Ykm[:, 0])))) * 1000.0
    if Xkm[0, -1] < Xkm[0, 0]:
        dx_m = -dx_m
    if Ykm[-1, 0] < Ykm[0, 0]:
        dy_m = -dy_m

    dz_dy, dz_dx = np.gradient(np.nan_to_num(Z_land, nan=0.0), dy_m, dx_m)
    az_rad = np.deg2rad(45.0)
    # Positive when slope faces toward the light source (NE)
    dir_slope = -(dz_dx * np.sin(az_rad) + dz_dy * np.cos(az_rad))

    # Nonlinear arctan intensity in [-1, 1], matching GMT's ATAN hillshading
    intens = (2.0 / np.pi) * np.arctan(dir_slope * 28.0)

    # Apply GMT-like illumination: darken shadows, brighten highlights
    pos = np.maximum(intens, 0.0)[..., None]
    neg = np.minimum(intens, 0.0)[..., None]
    rgb_shaded = rgb * (1.0 + 0.65 * neg)
    rgb_shaded = rgb_shaded + (1.0 - rgb_shaded) * (0.38 * pos)
    rgb_shaded = np.clip(rgb_shaded, 0.0, 1.0)

    rgba = np.ones(Z_land.shape + (4,), dtype=float)
    rgba[..., :3] = rgb_shaded
    return rgba

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

def _draw_2d_block_arrow(ax, fig, tip_or_tail_in, u_dir_2d, inward=True, speed=2.0, zorder=20000):
    """
    Draws a GMT-style red block arrow sized in physical inches and added to ax in transData.
    If inward=True (-Svh), tip_or_tail_in is the arrow tip (head).
    If inward=False (-Sv), tip_or_tail_in is the arrow tail.
    """
    ux, uy = float(u_dir_2d[0]), float(u_dir_2d[1])
    norm = np.hypot(ux, uy)
    if norm < 1e-9:
        return
    ux, uy = ux / norm, uy / norm
    # Perpendicular unit vector
    px, py = -uy, ux

    # Scale arrow dimensions gently with velocity magnitude (matching GMT 0.6 / 1.0 / 1.4 cm)
    scale = float(np.clip((max(speed, 0.2) / 2.5) ** 0.22, 0.75, 1.25))
    L = 0.46 * scale          # total arrow length in inches
    head_l = 0.22 * scale     # head length in inches
    shaft_hw = 0.070 * scale  # half-width of shaft in inches
    head_hw = 0.155 * scale   # half-width of head in inches

    x_a, y_a = float(tip_or_tail_in[0]), float(tip_or_tail_in[1])
    if inward:
        tip_x, tip_y = x_a, y_a
        tail_x, tail_y = tip_x - L * ux, tip_y - L * uy
    else:
        tail_x, tail_y = x_a, y_a
        tip_x, tip_y = tail_x + L * ux, tail_y + L * uy

    neck_x, neck_y = tip_x - head_l * ux, tip_y - head_l * uy

    verts_in = np.array([
        [tail_x + shaft_hw * px, tail_y + shaft_hw * py],
        [neck_x + shaft_hw * px, neck_y + shaft_hw * py],
        [neck_x + head_hw * px,  neck_y + head_hw * py],
        [tip_x,                  tip_y],
        [neck_x - head_hw * px,  neck_y - head_hw * py],
        [neck_x - shaft_hw * px, neck_y - shaft_hw * py],
        [tail_x - shaft_hw * px, tail_y - shaft_hw * py],
    ], dtype=float)

    # Convert physical inch vertices into ax.transData so they are invariant to savefig dpi & tight bbox
    verts_disp = fig.dpi_scale_trans.transform(verts_in)
    verts_data = ax.transData.inverted().transform(verts_disp)

    poly = MplPolygon(
        verts_data,
        closed=True,
        facecolor='#ff0000',
        edgecolor='black',
        linewidth=1.5,
        transform=ax.transData,
        zorder=zorder,
        clip_on=False,
    )
    fig.add_artist(poly)

def render(project, Xkm, Ykm, Zm_list, lake_list=[], pfl=None,
           azim=-60.0, elev=25.0, ve=10.0,
           xyw_geoms=None, xyw_linewidths=None,
           proj='persp', camdist=None, fov=None, time_val=None,
           t_ini=0.0, sea_level=0.0,
           zmin_override=None, zmax_override=None):
    
    Z_land = Zm_list[-1]
    Zothers_m = Zm_list[:-1]
    
    # Global limits
    all_z = [Z_land] + Zothers_m
    if lake_list:
        all_z.extend([lake['z'] for lake in lake_list])
    if pfl is not None and isinstance(pfl, dict) and 'z' in pfl:
        all_z.append(pfl['z'])

    g_zmin = float(np.nanmin([np.nanmin(z) for z in all_z]))
    g_zmax = float(np.nanmax([np.nanmax(z) for z in all_z]))

    if pfl is not None and isinstance(pfl, dict) and 'z' in pfl:
        pfl_z0_min = float(np.nanmin(pfl['z'][:, 0]))
        span = max(g_zmax - g_zmin, 1000.0)
        if g_zmin > pfl_z0_min - 0.10 * span:
            g_zmin = pfl_z0_min - 0.15 * span

    if zmin_override is not None:
        g_zmin = float(zmin_override)
    if zmax_override is not None:
        g_zmax = float(zmax_override)
    
    Zland_sc = (Z_land / 1000.0) * ve
    Zothers_sc = [(z / 1000.0) * ve for z in Zothers_m]
    zmin_p = (g_zmin / 1000.0) * ve
    zmax_p = (g_zmax / 1000.0) * ve

    fig = plt.figure(figsize=(14.0, 10.0), dpi=300)
    try:
        ax = fig.add_subplot(111, projection='3d', computed_zorder=False)
    except TypeError:
        ax = fig.add_subplot(111, projection='3d')
        ax.computed_zorder = False
    
    dx_grid = float(np.nanmax(Xkm) - np.nanmin(Xkm))
    dy_grid = float(np.nanmax(Ykm) - np.nanmin(Ykm))
    dz_grid = float(zmax_p - zmin_p)
    ax.set_box_aspect([dx_grid, dy_grid, dz_grid if dz_grid > 0 else 1.0])
    ax.set_zlim(zmin_p, zmax_p)

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

    # 1. Plot LAND with GMT bathymetry/topography palette and directional hillshading
    land_fc = color_and_illuminate_topo(Xkm, Ykm, Z_land, sea_level=sea_level)
    ax.plot_surface(Xkm, Ykm, Zland_sc, rstride=1, cstride=1,
                    facecolors=land_fc, linewidth=0, antialiased=True, shade=False, zorder=1)

    # 2. Plot Continental LAKES (batched Poly3DCollection)
    cell_dx = float(np.abs(np.mean(np.diff(Xkm[0, :]))))
    cell_dy = float(np.abs(np.mean(np.diff(Ykm[:, 0]))))
    if lake_list:
        hx, hy = cell_dx * 0.52, cell_dy * 0.52
        col_open = (20.0 / 255.0, 80.0 / 255.0, 255.0 / 255.0, 1.0)
        col_endo = (80.0 / 255.0, 70.0 / 255.0, 255.0 / 255.0, 1.0)
        lake_quads = []
        lake_colors = []
        z_lift = 0.003 * dz_grid

        for lake in lake_list:
            z_lake_m = lake['z']
            z_lake_plot = (z_lake_m / 1000.0) * ve + z_lift
            lcol = col_endo if lake.get('outl', 1) == 0 else col_open
            for (cx, cy) in lake['cells']:
                # Ensure lake quad sits at least slightly above local terrain node
                z_loc = float(bilinear(Xkm, Ykm, Zland_sc, cx, cy))
                zp = max(z_lake_plot, z_loc + z_lift)
                lake_quads.append(np.array([
                    [cx - hx, cy - hy, zp],
                    [cx + hx, cy - hy, zp],
                    [cx + hx, cy + hy, zp],
                    [cx - hx, cy + hy, zp],
                ], dtype=float))
                lake_colors.append(lcol)

        if lake_quads:
            poly_lakes = Poly3DCollection(lake_quads, facecolors=lake_colors,
                                          edgecolors='none', linewidths=0, zorder=10001)
            ax.add_collection3d(poly_lakes)

    # 3. Plot Subsurface (only if no vertical cross-section wall is plotted)
    if pfl is None:
        for i, Zp in enumerate(Zothers_sc[::-1], start=1):
            ax.plot_surface(Xkm, Ykm, Zp, linewidth=0, antialiased=True, shade=True, 
                            color='#808080', alpha=0.55, zorder=0)

    # 4. Plot RIVERS (batched Line3DCollection with line-of-sight ridge occlusion check)
    if xyw_geoms:
        az_rad = np.deg2rad(azim)
        el_rad = np.deg2rad(elev)
        ux_cam = float(np.cos(az_rad))
        uy_cam = float(np.sin(az_rad))
        tan_el = float(np.tan(el_rad))
        cell_d = np.hypot(cell_dx, cell_dy)
        x0_dom, xf_dom = float(np.nanmin(Xkm)), float(np.nanmax(Xkm))
        y0_dom, yf_dom = float(np.nanmin(Ykm)), float(np.nanmax(Ykm))

        geoms_arr = np.asarray(xyw_geoms, dtype=float)
        xm = 0.5 * (geoms_arr[:, 0, 0] + geoms_arr[:, 1, 0])
        ym = 0.5 * (geoms_arr[:, 0, 1] + geoms_arr[:, 1, 1])
        zm_sc = (0.5 * (geoms_arr[:, 0, 2] + geoms_arr[:, 1, 2]) / 1000.0) * ve

        occluded = np.zeros(len(xyw_geoms), dtype=bool)
        for step_mult in (1.5, 3.0, 5.0, 8.0, 12.0, 18.0):
            s = step_mult * cell_d
            xs = xm + s * ux_cam
            ys = ym + s * uy_cam
            in_dom = (xs >= x0_dom) & (xs <= xf_dom) & (ys >= y0_dom) & (ys <= yf_dom)
            if not np.any(in_dom):
                continue
            z_terr = bilinear(Xkm, Ykm, Zland_sc, xs, ys)
            z_ray = zm_sc + s * tan_el
            occluded |= in_dom & (z_terr > z_ray + 0.012 * dz_grid)

        vis_segs = []
        vis_lws = []
        z_riv_lift = 0.0015 * dz_grid
        for keep, geom, lw in zip(~occluded, xyw_geoms, xyw_linewidths):
            if not keep:
                continue
            s3 = np.column_stack([geom[:, 0], geom[:, 1], (geom[:, 2] / 1000.0) * ve + z_riv_lift])
            vis_segs.append(s3)
            vis_lws.append(lw)

        if vis_segs:
            riv_coll = Line3DCollection(vis_segs, colors='#0025dd', linewidths=vis_lws, zorder=10003)
            ax.add_collection3d(riv_coll)

    # 5. 2D Cross Section (Vertical Walls, matching tisc.topo3D.gmt.job)
    x0, xf = float(np.nanmin(Xkm)), float(np.nanmax(Xkm))
    y0, yf = float(np.nanmin(Ykm)), float(np.nanmax(Ykm))
    if pfl is not None and isinstance(pfl, dict) and 'z' in pfl:
        Lx = xf - x0
        Ly = yf - y0
        crop_limit1 = Ly
        crop_limit2 = Ly + Lx

        px = pfl['x']
        py = pfl['y']
        pl = pfl['l']
        pz = pfl['z'].copy()
        dens = pfl.get('densities', [])

        # Re-evaluate horizons from HRZ grid using exact bilinear interpolation when available
        # to eliminate inverse-distance sawtooth oscillation along boundary walls
        if len(Zm_list) == pz.shape[1]:
            for k in range(pz.shape[1]):
                pz[:, k] = bilinear(Xkm, Ykm, Zm_list[k], px, py)
        else:
            z_surf_wall = bilinear(Xkm, Ykm, Z_land, px, py)
            pz[:, -1] = z_surf_wall
            for k in range(pz.shape[1] - 1):
                pz[:, k] = np.minimum(pz[:, k], z_surf_wall)

        # Enforce monotonic non-decreasing horizons from base to top
        for k in range(1, pz.shape[1]):
            pz[:, k] = np.maximum(pz[:, k], pz[:, k - 1])

        # In GMT: view_azimut < 180 views from SE (East wall + South wall visible: l > crop_limit1)
        #         view_azimut >= 180 views from SW (West wall + South wall visible: l < crop_limit2)
        right_view = (np.cos(np.deg2rad(azim)) <= 0.0)

        is_boundary_box = (np.nanmax(pl) > crop_limit2 * 0.95)
        if is_boundary_box:
            if not right_view:
                mask = pl > crop_limit1
            else:
                mask = pl < crop_limit2
            if np.count_nonzero(mask) >= 2:
                px = px[mask]
                py = py[mask]
                pl = pl[mask]
                pz = pz[mask, :]

        # Colors from tisc.common.gmt.job and tisc.topo3D.gmt.job
        color_basement = (165.0 / 255.0, 120.0 / 255.0, 0.0 / 255.0)    # 165/120/0
        color_thrust   = (170.0 / 255.0, 110.0 / 255.0, 10.0 / 255.0)   # 170/110/10
        color_sediment = (250.0 / 255.0, 210.0 / 255.0, 0.0 / 255.0)    # 250/210/0

        if is_boundary_box and g_zmin <= sea_level <= g_zmax:
            sl_sc = (sea_level / 1000.0) * ve
            ax.plot([x0, xf], [y0, y0], [sl_sc, sl_sc], color='k', linewidth=2.0, zorder=10004)
            if right_view:
                ax.plot([x0, x0], [y0, y0 + 0.75 * Ly], [sl_sc, sl_sc], color='k', linewidth=2.0, zorder=10004)
            else:
                ax.plot([xf, xf], [y0, y0 + 0.55 * Ly], [sl_sc, sl_sc], color='k', linewidth=2.0, zorder=10004)

        pz_clamped = np.maximum(pz, g_zmin)
        pz_sc = (pz_clamped / 1000.0) * ve

        # Basement polygon (from zmin_p up to horizon 0)
        z_base_top = pz_sc[:, 0]
        z_base_bot = np.full_like(px, zmin_p)
        verts_base = np.column_stack([
            np.concatenate([px, px[::-1]]),
            np.concatenate([py, py[::-1]]),
            np.concatenate([z_base_top, z_base_bot[::-1]])
        ])
        poly_base = Poly3DCollection([verts_base], facecolors=[color_basement],
                                     edgecolors='k', linewidths=1.0, zorder=10005)
        ax.add_collection3d(poly_base)

        # Overlying units (col = 4 .. numcols-1 in GMT -> k = 1 .. num_horiz-1 in 0-based)
        num_horiz = pz_sc.shape[1]
        for k in range(1, num_horiz):
            z_bot = pz_sc[:, k - 1]
            z_top = pz_sc[:, k]
            verts_unit = np.column_stack([
                np.concatenate([px, px[::-1]]),
                np.concatenate([py, py[::-1]]),
                np.concatenate([z_top, z_bot[::-1]])
            ])
            d_val = dens[k] if k < len(dens) else 2150.0
            u_color = color_thrust if d_val >= 2400.0 else color_sediment
            poly_unit = Poly3DCollection([verts_unit], facecolors=[u_color],
                                         edgecolors='k', linewidths=1.0, zorder=10005 + k)
            ax.add_collection3d(poly_unit)

        if is_boundary_box:
            corner_x = x0 if right_view else xf
            ic = int(np.argmin((px - corner_x)**2 + (py - y0)**2))
            for idx in (0, ic, len(px) - 1):
                ax.plot([px[idx], px[idx]], [py[idx], py[idx]],
                        [zmin_p, pz_sc[idx, -1]], color='k', linewidth=2.2, zorder=10006 + num_horiz)

    ax.set_xlabel('X (km)'); ax.set_ylabel('Y (km)')
    ax.set_zlabel('Elevation (m)')
    z_ticks = nice_ticks(g_zmin, g_zmax, n=5)
    ax.set_zticks((z_ticks/1000.0)*ve)
    ax.set_zticklabels([f'{int(v)}' for v in z_ticks])
    ax.view_init(elev=elev, azim=azim)

    # Force canvas layout update so 3D-to-2D projection transforms are exact
    fig.canvas.draw()
    proj_mat = ax.get_proj()

    def _pt3d_to_inches(x3, y3, z3):
        x2, y2, _ = proj3d.proj_transform(x3, y3, z3, proj_mat)
        disp = ax.transData.transform((x2, y2))
        return fig.dpi_scale_trans.inverted().transform(disp)

    # 6. Plot Red Tectonic Velocity Arrows from <project>*.UNIT files
    unit_arrows = get_unit_arrows(project, Xkm, Ykm, time_val, t_ini=t_ini)
    if unit_arrows:
        d_km = 0.05 * max(dx_grid, dy_grid, 1.0)
        for arr in unit_arrows:
            bx, by = arr['x'], arr['y']
            vx, vy = arr['vx'], arr['vy']
            spd = arr['speed']
            ux, uy = vx / spd, vy / spd
            # Anchor at z=0 on boundaries (matching tisc.gmt.job), or at surface for interior units
            if arr['bname'] == 'INT':
                z_anch = max(float(bilinear(Xkm, Ykm, Zland_sc, bx, by)), 0.0)
            elif arr['bname'] == 'N':
                z_anch = max(float(bilinear(Xkm, Ykm, Zland_sc, bx, by)), 0.0)
            else:
                z_anch = 0.0

            p0_in = _pt3d_to_inches(bx, by, z_anch)
            p1_in = _pt3d_to_inches(bx + ux * d_km, by + uy * d_km, z_anch)
            u_dir_2d = p1_in - p0_in
            _draw_2d_block_arrow(ax, fig, p0_in, u_dir_2d, inward=arr['inward'], speed=spd, zorder=20000)

    # 7. Plot Time in My above the North edge of the 3D block (matching tisc.topo3D.gmt.job)
    if time_val is not None:
        if abs(time_val - round(time_val, 1)) > 1e-3:
            time_str = f"{time_val:.2f} My"
        else:
            time_str = f"{time_val:.1f} My"
        x_mid = 0.5 * (x0 + xf)
        p_top_in = _pt3d_to_inches(x_mid, yf, zmax_p)
        p_txt_disp = fig.dpi_scale_trans.transform((p_top_in[0], p_top_in[1] + 0.35))
        p_txt_data = ax.transData.inverted().transform(p_txt_disp)
        ax.text2D(
            p_txt_data[0], p_txt_data[1],
            time_str,
            transform=ax.transData,
            fontsize=20,
            fontweight='bold',
            ha='center',
            va='bottom',
            color='black',
            path_effects=[pe.withStroke(linewidth=3.0, foreground='white')],
            zorder=25000,
        )

    out = f'{project}.jpg'
    plt.savefig(out, bbox_inches='tight', dpi=300, pad_inches=0.25)
    plt.close(fig)
    return out

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description='TISC 3D Plotter: Visualizes HRZ/ST topography, XYW rivers, and Lakes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument('project', help='Project base name.')
    ap.add_argument('--azim', type=float, default=None, 
                    help='Camera Azimuth in degrees. If None, calculated dynamically from Timeini(-150) to Timefinal(-60).')
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

    args, unknown = ap.parse_known_args()

    if unknown:
        print(f"Warning: Unknown parameters passed and ignored: {unknown}")

    proj = args.project
    
    hrz = f'{proj}.hrz'
    st  = f'{proj}.st'
    pfl = f'{proj}.pfl'
    xyw = f'{proj}.xyw'
    lak = f'{proj}.lak'
    
    t_ini, t_fin, t_cur, zmin_prm, zmax_prm = get_time_info(proj)
    
    if args.azim is not None:
        azim_val = args.azim
    else:
        # Matches GMT view_ang_ini=240 (azim=-150) to view_ang_final=150 (azim=-60)
        if t_fin > t_ini:
            progress = (t_cur - t_ini) / (t_fin - t_ini)
            progress = max(0.0, min(1.0, progress))
            azim_val = -150.0 + (progress * (-60.0 - (-150.0)))
        else:
            azim_val = -150.0
            
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

    raw_segs, raw_ws, sea_level = process_xyw(xyw, Xkm, Ykm)
    
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
                 proj=args.proj, camdist=args.camdist, fov=args.fov, time_val=t_cur,
                 t_ini=t_ini, sea_level=sea_level,
                 zmin_override=zmin_prm, zmax_override=zmax_prm)
    print('WROTE:', jpg)

if __name__ == '__main__':
    main()