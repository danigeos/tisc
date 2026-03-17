#!/usr/bin/env python3
"""
General tAo plotting script with hydro-climate panel.

Usage:
    python tao.plot.py <projectname>

Required files:
    <project>.PRM   : domain & density parameters
    <project>.pfl   : geometry (interfaces) + densities + ages
    <project>.eros  : erosion, erosion rate, topography, discharge, precip...

Optional files:
    <project>.xzt   : deflection w + surface at multiple times (or one time)
    <project>.rain / .prec : (legacy) precipitation profile (x, precip)

Output:
    <project>.png
    <project>.svg
"""

import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

SEDIMENT_DENSITY_THRESHOLD = 2500.0  # density threshold for sediments


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def find_file(project, exts):
    """Return first existing file with one of the given extensions."""
    for ext in exts:
        path = project + ext
        if os.path.exists(path):
            return path
    return None


# ---------------------------------------------------------------------
# PRM
# ---------------------------------------------------------------------
def parse_prm(prm_path):
    """Extract key domain parameters from a tAo .PRM file."""
    params = {}
    with open(prm_path, encoding="latin-1") as f:
        text = f.read()

    def get(key, cast=float, default=None):
        # case-insensitive, skip comments
        pattern = rf"^\s*{key}\s+([^\n#]+)"
        m = re.search(pattern, text, flags=re.MULTILINE | re.IGNORECASE)
        if not m:
            return default
        token = m.group(1).split()[0]
        try:
            return cast(token)
        except ValueError:
            return default

    params["Nx"] = int(get("Nx", float, 0) or 0)
    for key in ("x0", "xf", "xmin", "xmax", "zmin", "zmax"):
        params[key] = get(key, float, None)

    params["denscrust"] = get("denscrust", float, None)
    params["densmantle"] = get("densmantle", float, None)

    return params


# ---------------------------------------------------------------------
# PFL
# ---------------------------------------------------------------------
def parse_pfl(pfl_path):
    """
    Parse a tAo *.pfl file.

    Returns
    -------
    x : (N,) array (km)
    z : (N, M) array (m)  interfaces along profile
    densities : (M,) array (kg/m3)
    ages : (M,) array (Ma)
    """
    with open(pfl_path, encoding="latin-1") as f:
        lines = f.read().splitlines()

    dens_line = None
    ages_line = None
    for line in lines:
        if "#" in line and "densit" in line.lower():
            dens_line = line
        if "#" in line and "x(km" in line.lower() and "age" in line.lower():
            ages_line = line
        if dens_line is not None and ages_line is not None:
            break

    if dens_line is None or ages_line is None:
        raise RuntimeError(
            f"Could not find densities / ages header in {os.path.basename(pfl_path)}"
        )

    densities = np.array(
        [float(s) for s in re.findall(r"[-+]?\d+(?:\.\d+)?", dens_line)],
        dtype=float,
    )
    ages = np.array(
        [float(s) for s in re.findall(r"[-+]?\d+(?:\.\d+)?", ages_line)],
        dtype=float,
    )

    data_rows = []
    for line in lines:
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            continue
        nums = [float(s) for s in re.findall(r"[-+]?\d+(?:\.\d+)?", line)]
        if not nums:
            continue
        data_rows.append(nums)

    data = np.array(data_rows, dtype=float)

    # First column = x, remaining columns = horizons
    x = data[:, 0]
    z = data[:, 1:]

    M_data = z.shape[1]
    M_header = len(ages)

    # If header & data horizons mismatch, trim to the minimum so we don't crash
    if M_header != M_data:
        M = min(M_header, M_data)
        print(
            f"Warning: {os.path.basename(pfl_path)} has {M_data} horizons in data "
            f"but {M_header} ages/densities in header; truncating to {M}."
        )
        z = z[:, :M]
        densities = densities[:M]
        ages = ages[:M]

    return x, z, densities, ages


# ---------------------------------------------------------------------
# EROS (robust, now pulls discharge / precip / evap)
# ---------------------------------------------------------------------
def read_eros(eros_path):
    """
    Read project.eros file in a robust way.

    We only require the first 3 numeric columns:
        x, erosion, eros_rate

    If present on a given line, we interpret extra numeric columns as:
        4th numeric: topography (m)
        5th numeric: discharge
        6th numeric: precipitation
        7th numeric: evaporation

    Any text columns (e.g. 'R') are ignored.
    Extra numeric columns beyond the 7th are ignored.

    Returns
    -------
    x : (N,) km
    erosion : (N,) m
    erosion_rate : (N,) m/My
    topo : (N,) m or None
    discharge : (N,) or None
    precip : (N,) or None
    evap : (N,) or None
    """
    xs = []
    erosions = []
    rates = []
    topos = []
    qs = []
    ps = []
    es = []

    with open(eros_path, encoding="latin-1") as f:
        for line in f:
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                continue
            # extract all floats on the line
            nums = [float(s) for s in re.findall(r"[-+]?\d+(?:\.\d+)?", line)]
            if len(nums) < 3:
                continue

            x = nums[0]
            erosion = nums[1]
            rate = nums[2]
            topo = nums[3] if len(nums) >= 4 else np.nan
            q = nums[4] if len(nums) >= 5 else np.nan
            # nums[5] is sedload, which we ignore here
            p = nums[6] if len(nums) >= 7 else np.nan
            e = nums[7] if len(nums) >= 8 else np.nan

            xs.append(x)
            erosions.append(erosion)
            rates.append(rate)
            topos.append(topo)
            qs.append(q)
            ps.append(p)
            es.append(e)

    if not xs:
        raise RuntimeError(
            f"No usable numeric data found in {os.path.basename(eros_path)}"
        )

    x_arr = np.array(xs, dtype=float)
    erosion_arr = np.array(erosions, dtype=float)
    rate_arr = np.array(rates, dtype=float)
    topo_arr = np.array(topos, dtype=float)
    q_arr = np.array(qs, dtype=float)
    p_arr = np.array(ps, dtype=float)
    e_arr = np.array(es, dtype=float)

    def clean_optional(arr):
        if np.all(np.isnan(arr)):
            return None
        return arr

    topo_arr = clean_optional(topo_arr)
    q_arr = clean_optional(q_arr)
    p_arr = clean_optional(p_arr)
    e_arr = clean_optional(e_arr)

    return x_arr, erosion_arr, rate_arr, topo_arr, q_arr, p_arr, e_arr


# ---------------------------------------------------------------------
# XZT (optional)
# ---------------------------------------------------------------------
def read_xzt(xzt_path):
    """
    Robust reader for project.xzt file.

    Supported layouts per row:

    (1) Single time:
        x, w, h

    (2) Multiple times:
        x, w(t0), h(t0), w(t1), h(t1), ..., w(tN), h(tN)

    In both cases we return the *final* (w, h) along x.

    Returns
    -------
    x : (N,) km
    w : (N,) m   (final time deflection)
    h : (N,) m or None   (final time surface)
    """
    rows = []
    with open(xzt_path, encoding="latin-1") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            nums = [float(s) for s in re.findall(r"[-+]?\d+(?:\.\d+)?", line)]
            if len(nums) < 2:
                continue
            rows.append(nums)

    if not rows:
        raise RuntimeError(f"No numeric data found in {os.path.basename(xzt_path)}")

    xs, ws, hs = [], [], []
    for nums in rows:
        x = nums[0]
        rest = nums[1:]

        if len(rest) == 2:
            # x, w, h (single time)
            w_last, h_last = rest[0], rest[1]
        elif len(rest) >= 2:
            # x, w0, h0, w1, h1, ..., wN, hN
            nt = len(rest) // 2
            w_last = rest[2 * (nt - 1)]
            h_last = rest[2 * (nt - 1) + 1] if 2 * nt - 1 < len(rest) else None
        else:
            continue

        xs.append(x)
        ws.append(w_last)
        hs.append(h_last if h_last is not None else np.nan)

    x_arr = np.array(xs, dtype=float)
    w_arr = np.array(ws, dtype=float)
    h_arr = np.array(hs, dtype=float)
    if np.all(np.isnan(h_arr)):
        h_arr = None
    return x_arr, w_arr, h_arr


# ---------------------------------------------------------------------
# Geometry plotting
# ---------------------------------------------------------------------
def plot_geometry(ax, x, z, densities, ages, denscrust=None, zmin=None, zmax=None):
    """Plot 2D section of geometry from pfl.

    - Area below the first z-horizon filled with denscrust
    - Inter-horizon intervals colored by density (non-sediments)
      or by age (sediments)
    - Ultra-thin black contours for non-sedimentary unit boundaries
    """
    nsurf = z.shape[1]
    nx = z.shape[0]

    # custom colour maps
    cmap_dens = LinearSegmentedColormap.from_list(
        "dens_cm", ["#4b2e19", "#7b1a1a"]  # dark brown -> dark red
    )
    cmap_age = LinearSegmentedColormap.from_list(
        "age_cm", ["#ffff66", "#ff8c00"]  # yellow -> dark orange
    )

    crust_dens = densities[densities > SEDIMENT_DENSITY_THRESHOLD]
    if crust_dens.size == 0:
        crust_dens = densities
    if denscrust is not None:
        crust_dens = np.append(crust_dens, denscrust)

    dens_norm = plt.Normalize(vmin=float(crust_dens.min()),
                              vmax=float(crust_dens.max()))
    age_norm = plt.Normalize(vmin=float(ages.min()),
                             vmax=float(ages.max()))

    dens_mappable = None
    age_mappable = None

    # --- region below first horizon = denscrust ---
    if zmin is not None:
        Xb = np.vstack([x, x])
        Z_lower = np.full_like(x, zmin)
        Z_upper = z[:, 0]
        Zmin = np.minimum(Z_lower, Z_upper)
        Zmax = np.maximum(Z_lower, Z_upper)
        Yb = np.vstack([Zmin, Zmax])

        val = denscrust if denscrust is not None else densities[0]
        vals = np.full((1, nx - 1), val)

        pm = ax.pcolormesh(
            Xb, Yb, vals,
            cmap=cmap_dens,
            norm=dens_norm,
            shading="auto",
        )
        dens_mappable = pm

        # top of crust contour
        ax.plot(x, z[:, 0], color="k", linewidth=0.2)

    # --- intervals between horizons ---
    for j in range(nsurf - 1):
        dens_j = densities[j]
        age_j = ages[j]

        Z_lower = z[:, j]
        Z_upper = z[:, j + 1]
        Zmin = np.minimum(Z_lower, Z_upper)
        Zmax = np.maximum(Z_lower, Z_upper)

        X = np.vstack([x, x])          # km
        Y = np.vstack([Zmin, Zmax])    # m

        vals = np.full(
            (1, nx - 1),
            dens_j if dens_j > SEDIMENT_DENSITY_THRESHOLD else age_j,
        )

        if dens_j > SEDIMENT_DENSITY_THRESHOLD:
            pm = ax.pcolormesh(
                X, Y, vals,
                cmap=cmap_dens,
                norm=dens_norm,
                shading="auto",
            )
            dens_mappable = pm

            # ultra-thin black contours for unit boundaries
            ax.plot(x, Zmin, color="k", linewidth=0.2)
            ax.plot(x, Zmax, color="k", linewidth=0.2)
        else:
            pm = ax.pcolormesh(
                X, Y, vals,
                cmap=cmap_age,
                norm=age_norm,
                shading="auto",
            )
            age_mappable = pm

    ax.set_ylabel("Elevation (m)")
    ax.grid(True, linestyle=":", linewidth=0.5)

    return dens_mappable, age_mappable, dens_norm, age_norm


# ---------------------------------------------------------------------
# Main

# ---------------------------------------------------------------------
def main():
    if len(sys.argv) != 2:
        print("Usage: python tao.plot.py <projectname>", file=sys.stderr)
        sys.exit(1)

    project = sys.argv[1]

    prm_path = find_file(project, [".PRM", ".prm"])
    pfl_path = find_file(project, [".pfl", ".PFL"])
    eros_path = find_file(project, [".eros", ".EROS"])
    xzt_path = find_file(project, [".xzt", ".XZT"])  # optional

    if prm_path is None:
        print(f"ERROR: Could not find {project}.PRM (or .prm)", file=sys.stderr)
        sys.exit(1)
    if pfl_path is None:
        print(f"ERROR: Could not find {project}.pfl (or .PFL)", file=sys.stderr)
        sys.exit(1)
    if eros_path is None:
        print(f"ERROR: Could not find {project}.eros (or .EROS)", file=sys.stderr)
        sys.exit(1)

    prm = parse_prm(prm_path)
    x_pfl, z_pfl, densities, ages = parse_pfl(pfl_path)

    # 3 panels: top = erosion / -w, middle = Q/P/E, bottom = geometry
    fig, (ax_hydro, ax_top, ax_geom) = plt.subplots(
        3, 1,
        figsize=(10, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 2]},
    )

    surface_x = None
    surface_z = None

    # -------------------------
    # top panel: erosion, erosion rate (and optionally -w)
    # -------------------------
    x_e, erosion, erosion_rate, topo_e, q_e, p_e, e_e = read_eros(eros_path)
    # convert precipitation & evaporation to m/yr (assuming original units are mm/yr)
    if p_e is not None:
        p_e = p_e / 1e3
    if e_e is not None:
        e_e = e_e / 1e3
    ax_top.plot(x_e, erosion_rate, label="erosion rate (m/Ma)")
    ax_top.plot(x_e, erosion, label="erosion (m)")
    if topo_e is not None:
        surface_x, surface_z = x_e, topo_e

    # deflection axis (right), only if we have xzt
    ax_defl = None
    if xzt_path is not None:
        ax_defl = ax_top.twinx()
        x_x, w, h = read_xzt(xzt_path)

        print(f"Deflection w range: min = {w.min():.3f} m, max = {w.max():.3f} m")

        w_plot = -w  # plot -w as requested
        ax_defl.plot(x_x, w_plot, label="-deflection w (m)", linestyle="--")
        ax_defl.set_ylabel("-deflection w (m)")

        if surface_x is None and h is not None:
            surface_x, surface_z = x_x, h

    # combined legend (left + right if deflection exists)
    lines1, labels1 = ax_top.get_legend_handles_labels()
    if ax_defl is not None:
        lines2, labels2 = ax_defl.get_legend_handles_labels()
        lines = lines1 + lines2
        labels = labels1 + labels2
    else:
        lines, labels = lines1, labels1
    if lines:
        ax_top.legend(lines, labels, loc="upper right")

    # no left axis title; units are in legend
    ax_top.set_ylabel("")
    ax_top.grid(True, linestyle=":", linewidth=0.5)

    # -------------------------
    # hydro panel (now top): precipitation & evaporation (left), discharge (right)
    # -------------------------
    ax_q = None
    # left axis: P & E (up to 5 m/yr)
    plotted_left = False
    if p_e is not None:
        # precipitation: solid blue
        ax_hydro.plot(x_e, p_e, label="precipitation", color="tab:blue", linestyle="-")
        plotted_left = True
    if e_e is not None:
        # evaporation: dashed blue
        ax_hydro.plot(x_e, e_e, label="evaporation", color="tab:blue", linestyle="--")
        plotted_left = True

    if plotted_left:
        ax_hydro.set_ylabel("P, E (m/yr)")

    ax_hydro.grid(True, linestyle=":", linewidth=0.5)

    # right axis: discharge (up to 2000 m³/s)
    if q_e is not None:
        ax_q = ax_hydro.twinx()
        # discharge: dark blue
        ax_q.plot(x_e, q_e, label="discharge", color="navy", linestyle="-")
        ax_q.set_ylabel("Discharge (m³/s)")

    # combined legend for hydro panel
    lines_left, labels_left = ax_hydro.get_legend_handles_labels()
    if ax_q is not None:
        lines_right, labels_right = ax_q.get_legend_handles_labels()
        hydro_lines = lines_left + lines_right
        hydro_labels = labels_left + labels_right
    else:
        hydro_lines, hydro_labels = lines_left, labels_left

    if hydro_lines:
        ax_hydro.legend(hydro_lines, hydro_labels, loc="upper right")

    # -------------------------
    # bottom panel: geometry
    # -------------------------
    dens_map, age_map, dens_norm, age_norm = plot_geometry(
        ax_geom,
        x_pfl,
        z_pfl,
        densities,
        ages,
        denscrust=prm.get("denscrust"),
        zmin=prm.get("zmin"),
        zmax=prm.get("zmax"),
    )

    # overlay present-day surface if known
    if surface_x is not None and surface_z is not None:
        ax_geom.plot(surface_x, surface_z, color="k", linewidth=0.8, zorder=5)

    ax_geom.set_xlabel("Distance (km)")



    # -------------------------
    # layout & colorbars at bottom
    # -------------------------
    if dens_map is not None or age_map is not None:
        fig.subplots_adjust(bottom=0.18)

    fig.suptitle(f"tAo project: {project}", fontsize=12)
    fig.tight_layout(rect=[0.0, 0.18, 1.0, 0.95])

    # then add small bottom colorbars
    if dens_map is not None:
        sm_dens = plt.cm.ScalarMappable(norm=dens_norm, cmap=dens_map.cmap)
        sm_dens.set_array([])
        cax1 = fig.add_axes([0.15, 0.07, 0.30, 0.02])
        cbar1 = fig.colorbar(sm_dens, cax=cax1, orientation="horizontal")
        cbar1.set_label("Density (kg/m³)")

    if age_map is not None:
        sm_age = plt.cm.ScalarMappable(norm=age_norm, cmap=age_map.cmap)
        sm_age.set_array([])
        cax2 = fig.add_axes([0.55, 0.07, 0.30, 0.02])
        cbar2 = fig.colorbar(sm_age, cax=cax2, orientation="horizontal")
        cbar2.set_label("Sediment age (Ma)")

    png_path = f"{project}.png"
    svg_path = f"{project}.svg"
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)

    print(f"Wrote {png_path} and {svg_path}")


if __name__ == "__main__":
    main()
