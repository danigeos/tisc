#!/usr/bin/env python3
"""
TISC multi-panel plotting script (single page).

Usage
-----
    python tisc.plot.py LagoMare
    python tisc.plot.py sinusoidal

Inputs for project P (if present)
---------------------------------
P.xyw   : drainage output (required)
P.st    : erosion output (topo, accum_eros, eros_rate)
P.xyzt  : isostatic subsidence/deflection (3rd column = w)
P.bas   : river basins (blocks separated by '>', area in '# END river basin ... of XXX km2: ...')

Output (single multi-panel figure)
----------------------------------
P.svg
P.jpg

Panels
------
1) Top row, full width:    Topography + drainage (exorheic/endorheic, lakes)
2) Middle row, full width: Erosion rate + isostatic subsidence (w) contours
3) Bottom row, left:       Length–elevation profile (largest basin)
4) Bottom row, right:      Chi–elevation profile (largest basin)
"""

import argparse
import io
import os
import re
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection


# ----------------------------- CONFIG -----------------------------
# Visual tuning defaults (you can adjust these quickly)
EXO_COLOR = "#4d8fd1"
ENDO_COLOR = "#5f7ad6"
HSPACE     = 0.28          # vertical space between rows
WSPACE     = 0.10          # horizontal space between main and cbar columns
WIDTH_RATIOS  = [20.0, 1.5]# left: plots, right: colorbars
HEIGHT_RATIOS = [1.0, 1.0, 0.45]  # topo, erosion, profiles (shorter row)
FIG_SIZE  = (12, 14)
DPI       = 250

# Final layout safety margins (figure-fraction)
LEFT_MIN, RIGHT_MAX = 0.05, 0.98
BOTTOM_MIN, TOP_MAX = 0.05, 0.97

# Final layout tweaks
PROFILE_SCALE = 1.20   # 20% larger (width & height)
PROFILE_DY    = -0.05  # 5% lower (figure coords)
CBAR_SHRINK   = 0.70   # 70% of original width (right edge fixed)


# ----------------------------- I/O -----------------------------

def load_drainage_xyw(path: str) -> pd.DataFrame:
    """
    Load drainage (.xyw) with at least:
      x y water sed type topo x_to y_to ...
    Extra cols are kept as generic col_i.
    """
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None)
    base_cols = [
        "x_km", "y_km",
        "water_m3s", "sed_kgs",
        "type", "topo_m",
        "x_to", "y_to",
    ]
    extra_cols = [f"col_{i}" for i in range(df.shape[1] - len(base_cols))]
    df.columns = base_cols + extra_cols
    return df


def load_erosion_st(path: str) -> pd.DataFrame:
    """
    Load erosion (.st):
      x y topo accum_erosion eros_rate
    """
    st_df = pd.read_csv(path, comment="#", sep=r"\s+", header=None)
    st_df.columns = [
        "x_km", "y_km",
        "topo_m_st",
        "accum_eros_m",
        "eros_rate_m_per_My",
    ]
    return st_df


def load_deflection_xyzt(path: str) -> pd.DataFrame:
    """
    Load deflection (.xyzt). We need col 2 -> 'w'.
    """
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None)
    ncol = df.shape[1]
    names = ["x_km", "y_km", "w"] + [f"col{i}" for i in range(3, ncol)]
    df.columns = names[:ncol]
    return df


def load_basins(path: str) -> List[Tuple[float, pd.DataFrame]]:
    """
    Parse .bas into list of (area_km2, basin_df).
    Each block starts with '>' and ends with a line like:
      '# END river basin N of AREA km2: ...'
    Data lines contain 17 columns:
      x_km y_km dischg_m3s masstr_kgs type topo_m length_km chi_km
      x_to_km y_to_km topo_to_m length_to_km chi_to_km eros_rate_m_per_My
      accum_eros_m level nodes_above
    """
    basins: List[Tuple[float, pd.DataFrame]] = []
    cur_lines: List[str] = []
    cur_area: Optional[float] = None

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#TISC output") or s.startswith("#x["):
                continue
            if s.startswith(">"):
                cur_lines = []
                cur_area = None
                continue
            if s.startswith("# END river basin"):
                m = re.search(r"of\s+([\d\.]+)\s+km2", s)
                if m:
                    cur_area = float(m.group(1))
                if cur_lines:
                    buf = io.StringIO("\n".join(cur_lines))
                    bdf = pd.read_csv(buf, sep=r"\s+", header=None)
                    bdf.columns = [
                        "x_km", "y_km",
                        "dischg_m3s", "masstr_kgs",
                        "type", "topo_m",
                        "length_km", "chi_km",
                        "x_to_km", "y_to_km",
                        "topo_to_m", "length_to_km",
                        "chi_to_km", "eros_rate_m_per_My",
                        "accum_eros_m", "level", "nodes_above",
                    ]
                    basins.append((cur_area if cur_area is not None else np.nan, bdf))
                cur_lines = []
                cur_area = None
                continue
            if s.startswith("#"):
                continue
            cur_lines.append(s)

    return basins


# ------------------------ Helpers / Grid / Routing ------------------------

def estimate_grid_spacing(df: pd.DataFrame) -> Tuple[float, float]:
    xs = np.sort(df["x_km"].unique())
    ys = np.sort(df["y_km"].unique())
    dxs = np.diff(xs)
    dys = np.diff(ys)
    dx = float(np.min(dxs[dxs > 0])) if np.any(dxs > 0) else 1.0
    dy = float(np.min(dys[dys > 0])) if np.any(dys > 0) else 1.0
    return dx, dy


def build_grid(df: pd.DataFrame, value_col: str):
    grid_df = df.pivot(index="y_km", columns="x_km", values=value_col)
    xs = grid_df.columns.to_numpy()
    ys = grid_df.index.to_numpy()
    return xs, ys, grid_df.to_numpy()


def classify_routing(df: pd.DataFrame) -> pd.Series:
    """
    Label rows 'exo' if they reach sea / sub-sea or exit domain;
    otherwise 'endo' (loops / internal sinks).
    """
    coords = list(zip(df["x_km"].values, df["y_km"].values))
    index_from_coord = {c: i for i, c in enumerate(coords)}
    topo_clean = df["topo_m"].copy()
    topo_clean[topo_clean <= -3e4] = np.nan

    status: dict[int, str] = {}

    def step(idx: int, max_steps: int = 10000) -> str:
        if idx in status:
            return status[idx]
        visited = set()
        cur = idx
        for _ in range(max_steps):
            if cur in status:
                res = status[cur]
                break
            if cur in visited:
                res = "endo"
                break
            visited.add(cur)

            row = df.iloc[cur]
            topo_val = topo_clean.iloc[cur]
            if (row["type"] == "S") or (np.isfinite(topo_val) and topo_val <= 0):
                res = "exo"
                break

            tgt_coord = (row["x_to"], row["y_to"])
            if tgt_coord == (row["x_km"], row["y_km"]):
                res = "endo"
                break

            nxt = index_from_coord.get(tgt_coord)
            if nxt is None:
                res = "exo"
                break

            cur = nxt
        else:
            res = "endo"

        for v in visited:
            status[v] = res
        return res

    return pd.Series([step(i) for i in range(len(df))], index=df.index, name="drain_type")


def nice_ticks(vmin: float, vmax: float, target: int = 8) -> np.ndarray:
    span = vmax - vmin
    if span <= 0:
        return np.array([vmin])
    raw_step = span / target
    mag = 10 ** np.floor(np.log10(raw_step))
    candidates = np.array([1, 2, 2.5, 5, 10]) * mag
    step = candidates[np.argmin(np.abs(candidates - raw_step))]
    start = np.ceil(vmin / step) * step
    end = np.floor(vmax / step) * step
    return np.arange(start, end + 0.5 * step, step)


# ----------------------------- Main Plotter -----------------------------

def plot_tisc_project(project: str,
                      min_discharge_quantile: float = 0.3,
                      max_link_factor: float = 3.0):

    drainage_path = f"{project}.xyw"
    erosion_path  = f"{project}.st"
    xyzt_path     = f"{project}.xyzt"
    bas_path      = f"{project}.bas"

    if not os.path.exists(drainage_path):
        raise FileNotFoundError(f"Drainage file not found: {drainage_path}")

    have_erosion = os.path.exists(erosion_path)
    have_xyzt    = os.path.exists(xyzt_path)
    have_bas     = os.path.exists(bas_path)

    # -- Topography & drainage --
    df = load_drainage_xyw(drainage_path)
    df["topo_m_clean"] = df["topo_m"]
    df.loc[df["topo_m_clean"] <= -3e4, "topo_m_clean"] = np.nan

    xs, ys, topo_grid = build_grid(df, "topo_m_clean")
    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()

    valid_topo = topo_grid[np.isfinite(topo_grid)]
    sea_vals   = valid_topo[valid_topo < 0]
    land_vals  = valid_topo[valid_topo > 0]
    topo_vmin  = np.nanpercentile(sea_vals, 25.0) if sea_vals.size > 0 else np.nanpercentile(valid_topo, 25.0)
    topo_vmax  = np.nanmax(land_vals) if land_vals.size > 0 else np.nanmax(valid_topo)
    topo_vmin  = min(topo_vmin, -1.0)
    topo_vmax  = max(topo_vmax,  1.0)

    bath_cmap = mpl.colormaps.get_cmap("Blues_r").resampled(128)
    land_cmap = mpl.colormaps.get_cmap("terrain").resampled(128)
    colors = np.vstack((bath_cmap(np.linspace(0.2, 1.0, 128)),
                        land_cmap(np.linspace(0.25, 1.0, 128))))
    bathy_land_cmap = mpl.colors.ListedColormap(colors, name="bathy_land")
    topo_norm = mpl.colors.TwoSlopeNorm(vmin=topo_vmin, vcenter=0.0, vmax=topo_vmax)

    df["drain_type"] = classify_routing(df)

    dx_grid, dy_grid = estimate_grid_spacing(df)
    max_link_len = max_link_factor * max(dx_grid, dy_grid)

    flow_df = df[(df["water_m3s"] > 0) & (df["type"] != "S")].copy()
    dx = flow_df["x_to"] - flow_df["x_km"]
    dy = flow_df["y_to"] - flow_df["y_km"]
    dist = np.sqrt(dx * dx + dy * dy)
    flow_df = flow_df[dist <= max_link_len].copy()
    if not flow_df.empty:
        q_thresh = flow_df["water_m3s"].quantile(min_discharge_quantile)
        flow_df = flow_df[flow_df["water_m3s"] >= q_thresh].copy()

    exo_riv  = flow_df[flow_df["drain_type"] == "exo"].copy()
    endo_riv = flow_df[flow_df["drain_type"] == "endo"].copy()

    def make_segments(subdf: pd.DataFrame):
        if subdf.empty:
            return [], []
        w = subdf["water_m3s"].values
        wmin, wmax = np.percentile(w, [1, 99])
        w = np.clip(w, wmin, wmax)
        wn = (w - wmin) / (wmax - wmin + 1e-9)
        widths = 0.2 + 2.5 * wn
        segs = list(zip(subdf[["x_km", "y_km"]].to_numpy(),
                        subdf[["x_to", "y_to"]].to_numpy()))
        return segs, widths

    exo_segs,  exo_w  = make_segments(exo_riv)
    endo_segs, endo_w = make_segments(endo_riv)

    # Lakes (remove single-node)
    lake_df = df[df["type"] == "L"].copy()
    lake_coords = list(zip(lake_df["x_km"].values, lake_df["y_km"].values))
    lake_set = set(lake_coords)
    multi_mask = []
    for (x, y) in lake_coords:
        has_neighbor = False
        for dxn in (-dx_grid, 0.0, dx_grid):
            for dyn in (-dy_grid, 0.0, dy_grid):
                if dxn == 0.0 and dyn == 0.0:
                    continue
                if (x + dxn, y + dyn) in lake_set:
                    has_neighbor = True
                    break
            if has_neighbor:
                break
        multi_mask.append(has_neighbor)
    lake_df_multi = lake_df[np.array(multi_mask)] if len(lake_df) else lake_df
    exo_lakes  = lake_df_multi[lake_df_multi["drain_type"] == "exo"]
    endo_lakes = lake_df_multi[lake_df_multi["drain_type"] == "endo"]

    # -- Erosion & norm --
    eros_grid = None
    eros_cmap = None
    eros_norm = None
    if have_erosion:
        st_df = load_erosion_st(erosion_path)
        eros_grid_df = st_df.pivot(index="y_km", columns="x_km", values="eros_rate_m_per_My")
        eros_grid = eros_grid_df.reindex(index=ys, columns=xs).to_numpy()
        eros_vals = eros_grid[np.isfinite(eros_grid)]
        eros_max  = np.nanpercentile(np.abs(eros_vals), 99.0) if eros_vals.size > 0 else 1.0
        if eros_max == 0: eros_max = 1.0
        eros_norm = mpl.colors.TwoSlopeNorm(vmin=-eros_max, vcenter=0.0, vmax=eros_max)
        eros_cmap = mpl.colormaps.get_cmap("seismic")

    # -- Deflection (w) --
    w_grid = None
    if have_xyzt:
        xyzt_df = load_deflection_xyzt(xyzt_path)
        w_grid_df = xyzt_df.pivot(index="y_km", columns="x_km", values="w")
        w_grid = w_grid_df.reindex(index=ys, columns=xs).to_numpy()

    # -- Prepare profiles (largest basin) --
    profiles = None
    if have_bas and eros_grid is not None and eros_norm is not None:
        basins = load_basins(bas_path)
        if basins:
            largest_area, bas = max(basins, key=lambda t: t[0])
            bas = bas.copy()
            mask_len = np.isfinite(bas["length_km"]) & np.isfinite(bas["length_to_km"]) \
                       & np.isfinite(bas["topo_m"]) & np.isfinite(bas["topo_to_m"])
            mask_chi = np.isfinite(bas["chi_km"]) & np.isfinite(bas["chi_to_km"]) \
                       & np.isfinite(bas["topo_m"]) & np.isfinite(bas["topo_to_m"])
            bas_len = bas[mask_len]
            bas_chi = bas[mask_chi]

            # widths by discharge
            q = bas["dischg_m3s"].values
            if q.size > 0:
                qmin, qmax = np.percentile(q, [1, 99])
                q_clip = np.clip(q, qmin, qmax)
                q_norm = (q_clip - qmin) / (qmax - qmin + 1e-9)
                widths = 0.5 + 2.5 * q_norm
            else:
                widths = np.full(len(bas), 1.0)

            # colors by eros_rate (same norm/cmap as erosion map)
            colors_profile = eros_cmap(eros_norm(bas["eros_rate_m_per_My"].values))

            len_segs = np.stack([
                bas_len[["length_km", "topo_m"]].to_numpy(),
                bas_len[["length_to_km", "topo_to_m"]].to_numpy()
            ], axis=1) if len(bas_len) > 0 else np.empty((0, 2, 2))

            chi_segs = np.stack([
                bas_chi[["chi_km", "topo_m"]].to_numpy(),
                bas_chi[["chi_to_km", "topo_to_m"]].to_numpy()
            ], axis=1) if len(bas_chi) > 0 else np.empty((0, 2, 2))

            profiles = dict(
                largest_area=largest_area,
                bas_len=bas_len, bas_chi=bas_chi,
                len_segs=len_segs, chi_segs=chi_segs,
                widths=widths, colors_profile=colors_profile
            )

    # ------------------ Figure layout (single page) ------------------
    fig = plt.figure(figsize=FIG_SIZE, dpi=DPI)

    # 3 rows x 2 columns (right col for colorbars)
    gs = fig.add_gridspec(
        3, 2,
        width_ratios=WIDTH_RATIOS,
        height_ratios=[
            HEIGHT_RATIOS[0],
            HEIGHT_RATIOS[1] if eros_grid is not None else 0.001,
            HEIGHT_RATIOS[2] if profiles is not None else 0.001,
        ],
        wspace=WSPACE,
        hspace=HSPACE,
    )

    # Axes
    ax1 = fig.add_subplot(gs[0, 0])  # topo
    ax2 = fig.add_subplot(gs[1, 0]) if eros_grid is not None else None  # erosion
    if profiles is not None:
        bottom_gs = gs[2, 0].subgridspec(1, 2, wspace=0.28)
        ax_len = fig.add_subplot(bottom_gs[0, 0])
        ax_chi = fig.add_subplot(bottom_gs[0, 1])
    else:
        ax_len = ax_chi = None

    # ---- Panel 1: Topography + drainage ----
    im1 = ax1.imshow(
        topo_grid,
        extent=[x_min, x_max, y_min, y_max],
        origin="lower",
        aspect="equal",
        cmap=bathy_land_cmap,
        norm=topo_norm,
    )
    if exo_segs:
        ax1.add_collection(LineCollection(exo_segs, linewidths=exo_w, color=EXO_COLOR, alpha=0.8, label="Exorheic rivers"))
    if endo_segs:
        ax1.add_collection(LineCollection(endo_segs, linewidths=endo_w, color=ENDO_COLOR, alpha=0.8, label="Endorheic rivers"))
    if not exo_lakes.empty:
        ax1.scatter(exo_lakes["x_km"], exo_lakes["y_km"], s=10, edgecolors="none", c=EXO_COLOR, alpha=0.9, label="Exorheic lakes")
    if not endo_lakes.empty:
        ax1.scatter(endo_lakes["x_km"], endo_lakes["y_km"], s=12, edgecolors="none", c=ENDO_COLOR, alpha=0.9, label="Endorheic lakes")
    ax1.set_xlim(x_min, x_max); ax1.set_ylim(y_min, y_max)
    ax1.set_ylabel("y (km)")
    ax1.set_title(f"{project}: Topography & Drainage")
    ax1.set_xticks(nice_ticks(x_min, x_max)); ax1.set_yticks(nice_ticks(y_min, y_max))
    ax1.legend(loc="lower left", fontsize=8)
    cax1 = fig.add_subplot(gs[0, 1])
    cbar1 = fig.colorbar(im1, cax=cax1); cbar1.set_label("Topography (m)")

    # ---- Panel 2: Erosion + w contours ----
    if eros_grid is not None and ax2 is not None:
        im2 = ax2.imshow(
            eros_grid,
            extent=[x_min, x_max, y_min, y_max],
            origin="lower",
            aspect="equal",
            cmap=eros_cmap,
            norm=eros_norm,
        )
        ax2.set_xlabel("x (km)"); ax2.set_ylabel("y (km)")
        title2 = f"{project}: Erosion rate (m/My)"
        if w_grid is not None:
            title2 += " & isostatic subsidence (w) contours"
        ax2.set_title(title2)
        ax2.set_xticks(nice_ticks(x_min, x_max)); ax2.set_yticks(nice_ticks(y_min, y_max))
        cax2 = fig.add_subplot(gs[1, 1])
        cbar2 = fig.colorbar(im2, cax=cax2); cbar2.set_label("Erosion rate (m/My)")

        if w_grid is not None:
            w_valid = w_grid[np.isfinite(w_grid)]
            if w_valid.size > 0:
                levels = np.linspace(np.nanmin(w_valid), np.nanmax(w_valid), 10)
                cs = ax2.contour(xs, ys, w_grid, levels=levels, colors="k", linewidths=0.5, alpha=0.8)
                ax2.clabel(cs, inline=True, fontsize=6, fmt="%.1f")

    # ---- Panel 3 & 4: Profiles (largest basin) ----
    if profiles is not None and ax_len is not None and ax_chi is not None:
        la = profiles["largest_area"]
        bas_len = profiles["bas_len"]; bas_chi = profiles["bas_chi"]
        len_segs = profiles["len_segs"]; chi_segs = profiles["chi_segs"]
        widths = profiles["widths"]; colors_profile = profiles["colors_profile"]

        if len(len_segs) > 0:
            lc_len = LineCollection(len_segs, linewidths=widths[bas_len.index], colors=colors_profile[bas_len.index], alpha=0.9)
            ax_len.add_collection(lc_len)
        ax_len.set_xlabel("Length along river (km)")
        ax_len.set_ylabel("Elevation (m)")
        ax_len.set_title(f"{project}: Length–elevation profile\nLargest basin ({la:.0f} km²)")
        ax_len.autoscale()

        if len(chi_segs) > 0:
            lc_chi = LineCollection(chi_segs, linewidths=widths[bas_chi.index], colors=colors_profile[bas_chi.index], alpha=0.9)
            ax_chi.add_collection(lc_chi)
        ax_chi.set_xlabel("Chi (km)")
        ax_chi.set_title(f"{project}: Chi–elevation profile")
        ax_chi.autoscale()

    # ------------------ Final layout tweaks (safe & robust) ------------------
    # 1) Profiles: scale and lower a bit, but clamp to margins
    for ax in [ax_len, ax_chi]:
        if ax is None:
            continue
        left, bottom, width, height = ax.get_position().bounds
        width  *= PROFILE_SCALE
        height *= PROFILE_SCALE
        bottom += PROFILE_DY
        left   = max(LEFT_MIN,  min(left,   RIGHT_MAX - width))
        bottom = max(BOTTOM_MIN, min(bottom, TOP_MAX  - height))
        ax.set_position([left, bottom, width, height])

    # 2) Colorbars: shrink width to CBAR_SHRINK, keep right edge fixed and clamp
    for cax in [locals().get("cax1"), locals().get("cax2")]:
        if cax is None:
            continue
        left, bottom, width, height = cax.get_position().bounds
        new_w   = width * CBAR_SHRINK
        new_l   = left + (width - new_w)  # right edge fixed
        new_l   = max(LEFT_MIN, min(new_l, RIGHT_MAX - new_w))
        bottom  = max(BOTTOM_MIN, min(bottom, TOP_MAX - height))
        cax.set_position([new_l, bottom, new_w, height])

    # ------------------ Save ------------------
    svg_path = f"{project}.svg"
    jpg_path = f"{project}.jpg"
    fig.savefig(svg_path)
    fig.savefig(jpg_path, dpi=300)
    print(f"Saved: {svg_path}\nSaved: {jpg_path}")


# ----------------------------- CLI -----------------------------

def main():
    parser = argparse.ArgumentParser(description="Plot TISC outputs for a project into one multi-panel figure.")
    parser.add_argument("project", help="Project name (expects project.xyw and optionally .st, .xyzt, .bas)")
    parser.add_argument("--min-discharge-quantile", type=float, default=0.3,
                        help="Quantile below which rivers are not plotted (default 0.3).")
    parser.add_argument("--max-link-factor", type=float, default=3.0,
                        help="Max allowed routing-link length as factor of grid spacing.")
    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Warning: Unknown parameters passed and ignored: {unknown}")

    plot_tisc_project(
        project=args.project,
        min_discharge_quantile=args.min_discharge_quantile,
        max_link_factor=args.max_link_factor,
    )


if __name__ == "__main__":
    main()
