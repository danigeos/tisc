#!/usr/bin/env python3
"""
tisc.UNITdrawer.py — km units, inline Z input (type & Enter),
user-order ring interpolation, live 2D+3D (Shift+Left-drag tilt),
save to <project>.UNIT with Z line (no 'Z' prefix) + coords (6 sig figs),
snap-to-grid inside domain, accept clicks anywhere (even beyond axes).
"""

import argparse, re, sys
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

# Disable toolbar & default keymaps that interfere
import matplotlib as mpl
mpl.rcParams['toolbar'] = 'none'
for k in ('keymap.fullscreen','keymap.home','keymap.back','keymap.forward',
          'keymap.pan','keymap.zoom','keymap.grid','keymap.yscale',
          'keymap.xscale','keymap.save','keymap.copy','keymap.quit'):
    mpl.rcParams[k] = []

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.path import Path as MplPath
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ---------------- PRM parsing ---------------- #

PRM_FIELD_REGEX = r'(?im)^\s*{name}\s+([+\-]?\d+(?:\.\d+)?(?:e[+\-]?\d+)?)'

def parse_prm(prm_path: Path):
    txt = prm_path.read_text()
    def get(name, required=True):
        m = re.search(PRM_FIELD_REGEX.format(name=re.escape(name)), txt)
        if not m:
            if required: raise ValueError(f"Parameter '{name}' not found in {prm_path}")
            return None
        return float(m.group(1))
    data = {
        "xmin": get("xmin"), "xmax": get("xmax"),
        "ymin": get("ymin"), "ymax": get("ymax"),
        "Nx":   int(get("Nx", required=False)) if re.search(PRM_FIELD_REGEX.format(name="Nx"), txt) else None,
        "Ny":   int(get("Ny", required=False)) if re.search(PRM_FIELD_REGEX.format(name="Ny"), txt) else None,
    }
    if data["Ny"] is None: data["Ny"] = data["Nx"]
    return data

# ---------------- Geometry helpers ---------------- #

def point_in_poly_grid(xs, ys, poly: List[Tuple[float,float]]):
    path = MplPath(poly, closed=True)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel()])
    return path.contains_points(pts).reshape(Y.shape)

def dist_to_poly_boundary(X, Y, poly: List[Tuple[float,float]]):
    px = np.asarray([p[0] for p in poly]); py = np.asarray([p[1] for p in poly])
    xi, yi = px[:-1], py[:-1]; xj, yj = px[1:], py[1:]
    vx, vy = xj - xi, yj - yi
    Xe, Ye = X[..., None], Y[..., None]
    wx, wy = Xe - xi, Ye - yi
    seglen2 = vx*vx + vy*vy + 1e-30
    t = np.clip((wx*vx + wy*vy)/seglen2, 0.0, 1.0)
    projx, projy = xi + t*vx, yi + t*vy
    dx, dy = Xe - projx, Ye - projy
    return np.sqrt(np.min(dx*dx + dy*dy, axis=-1))

# ---------------- App ---------------- #

class PolyApp:
    def __init__(self, dom_m, project: str, out_dir: Path):
        self.project = project
        self.out_dir = out_dir

        m2km = lambda v: None if v is None else v/1000.0
        self.dom_km = {
            "xmin": m2km(dom_m.get("xmin")), "xmax": m2km(dom_m.get("xmax")),
            "ymin": m2km(dom_m.get("ymin")), "ymax": m2km(dom_m.get("ymax")),
            "Nx": dom_m.get("Nx") or 81, "Ny": dom_m.get("Ny") or (dom_m.get("Nx") or 81),
        }

        # Build grid (km) for snapping + interpolation
        self.xs, self.ys, (self.Xg, self.Yg) = self._ensure_grid()

        # 2D editor
        self.fig, self.ax = plt.subplots(num=f"TISC UNIT drawer (km) — {project}", facecolor='white')
        w, h = self.fig.get_size_inches(); self.fig.set_size_inches(w*1.5, h*1.5, forward=True)
        self.ax.set_title(f"Draw polygons for '{project}'.UNIT (km)")
        self.ax.set_xlabel("x (km)"); self.ax.set_ylabel("y (km)")
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_xlim(self.xs.min(), self.xs.max())
        self.ax.set_ylim(self.ys.min(), self.ys.max())
        self.ax.grid(True, alpha=0.3)

        # Shortcuts as figure subtitle (outside plot)
        shortcuts = ("Shortcuts: Left-click add (snaps inside domain) | Right-click/Enter close→Z | "
                     "Shift+Left-drag tilt 3D | z undo | d delete | n new | s save | q/Esc quit")
        self.fig.subplots_adjust(top=0.88)
        self.fig.suptitle(shortcuts, y=0.98, fontsize=9)

        # State
        self.current_pts: List[Tuple[float,float]] = []
        self.current_line = Line2D([], [], marker='o', linestyle='-', lw=1.5, alpha=0.9)
        self.ax.add_line(self.current_line)
        self.polygons: List[List[Tuple[float,float]]] = []
        self.poly_patches: List[MplPolygon] = []
        self.poly_z: List[float] = []

        # Inline Z input (no widget focus needed)
        self.awaiting_z = False
        self._z_buffer = ""
        self._z_ax = self.fig.add_axes([0.15, 0.02, 0.35, 0.06]); self._z_ax.set_axis_off()

        # Elevation: 2D + 3D
        self.im = None; self.cbar = None
        self.fig3d = plt.figure(num=f"Elevation 3D — {self.project}", figsize=(8, 6))
        self.fig3d.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)  # tighter
        self.ax3d = self.fig3d.add_subplot(111, projection='3d')
        # Start positive Z up with slight tilt
        self.azim = -60.0; self.elev = 30.0
        self.ax3d.view_init(elev=self.elev, azim=self.azim)
        self.ax3d.set_xlabel("x (km)"); self.ax3d.set_ylabel("y (km)"); self.ax3d.set_zlabel("Z")
        self.cbar3d = None  # single 3D colorbar

        # Shift+drag tracking
        self.shift_down = False
        self.tilt_dragging = False
        self.tilt_last_xy = (0.0, 0.0)

        # Status (small inside axes)
        self.status = self.ax.text(
            0.01, 0.01, "0 polygons | 0 points (km)", transform=self.ax.transAxes,
            va='bottom', ha='left', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, lw=0.5),
        )

        # Events
        self.fig.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_button_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.fig.canvas.mpl_connect('key_release_event', self.on_key_release)

        # Prime empty Z plots
        self.compute_and_plot_z()

    # ---------- Grid ---------- #
    def _ensure_grid(self):
        Nx, Ny = int(self.dom_km["Nx"]), int(self.dom_km["Ny"])
        xmin, xmax = self.dom_km["xmin"], self.dom_km["xmax"]
        ymin, ymax = self.dom_km["ymin"], self.dom_km["ymax"]
        xs = np.linspace(xmin, xmax, Nx); ys = np.linspace(ymin, ymax, Ny)
        return xs, ys, np.meshgrid(xs, ys)

    def _snap_to_grid(self, x, y):
        ix = np.abs(self.xs - x).argmin()
        iy = np.abs(self.ys - y).argmin()
        return float(self.xs[ix]), float(self.ys[iy])

    def _in_domain(self, x: float, y: float) -> bool:
        return (self.xs.min() <= x <= self.xs.max()) and (self.ys.min() <= y <= self.ys.max())

    def _display_to_data(self, x_disp: float, y_disp: float):
        """Map figure/window pixel coords → data coords, even when pointer is outside axes."""
        inv = self.ax.transData.inverted()
        x_data, y_data = inv.transform((x_disp, y_disp))
        return float(x_data), float(y_data)

    # ---------- Events ---------- #

    def on_button_press(self, event):
        # Shift + Left-drag for 3D tilt can start anywhere over the figure
        if self.shift_down and event.button == 1 and event.x is not None and event.y is not None:
            self.tilt_dragging = True
            self.tilt_last_xy = (event.x, event.y)
            return

        # While waiting for Z entry, ignore editing clicks
        if self.awaiting_z:
            return

        # Accept clicks ANYWHERE in the window (even outside axes)
        if event.button in (1, 3) and event.x is not None and event.y is not None:
            x_d, y_d = self._display_to_data(event.x, event.y)

            if event.button == 1:
                # Snap inside domain; keep raw coords outside domain
                if self._in_domain(x_d, y_d):
                    xg, yg = self._snap_to_grid(x_d, y_d)
                    self.current_pts.append((xg, yg))
                else:
                    self.current_pts.append((x_d, y_d))
                self._redraw_current()

            elif event.button == 3:
                self.close_current_polygon()

    def on_button_release(self, event):
        self.tilt_dragging = False

    def on_motion(self, event):
        if not self.tilt_dragging: return
        if event.x is None or event.y is None: return
        dx = event.x - self.tilt_last_xy[0]
        dy = event.y - self.tilt_last_xy[1]
        self.azim = (self.azim + dx * 0.3) % 360.0
        self.elev = float(np.clip(self.elev - dy * 0.2, 5.0, 90.0))
        self.ax3d.view_init(elev=self.elev, azim=self.azim)
        self.fig3d.canvas.draw_idle()
        self.tilt_last_xy = (event.x, event.y)

    def on_key_press(self, event):
        # Track shift state
        if event.key in ('shift', 'shiftleft', 'shiftright'):
            self.shift_down = True
            return

        # Inline Z input: capture keys globally until Enter
        if self.awaiting_z:
            if event.key in ('enter','return'):
                self._submit_z_from_buffer(); return
            elif event.key == 'backspace':
                self._z_buffer = self._z_buffer[:-1]
            elif event.key == 'delete':
                self._z_buffer = ""
            elif len(event.key) == 1 and event.key in "0123456789.-+eE":
                self._z_buffer += event.key
            self._render_z_overlay()
            return

        # Normal editing keys
        if event.key in ('enter','return'):
            self.close_current_polygon()
        elif event.key == 'z' and self.current_pts:
            self.current_pts.pop(); self._redraw_current()
        elif event.key == 'd':
            self.delete_last_polygon()
        elif event.key == 'n' and self.current_pts:
            self.close_current_polygon()
        elif event.key in ('q','escape'):
            plt.close(self.fig); plt.close(self.fig3d)
        elif event.key == 's':
            self.save_and_exit()

    def on_key_release(self, event):
        if event.key in ('shift', 'shiftleft', 'shiftright'):
            self.shift_down = False

    # ---------- Helpers ---------- #

    def _redraw_current(self):
        xs = [p[0] for p in self.current_pts]; ys = [p[1] for p in self.current_pts]
        self.current_line.set_data(xs, ys); self._update_status(); self.fig.canvas.draw_idle()

    def _update_status(self):
        self.status.set_text(f"{len(self.polygons)} polygons | {len(self.current_pts)} points (km)")

    # ---------- Inline Z input overlay ---------- #

    def _open_z_prompt(self):
        self.awaiting_z = True
        self._z_buffer = ""
        self._render_z_overlay()

    def _render_z_overlay(self):
        self._z_ax.cla(); self._z_ax.set_axis_off()
        msg = f"Z value: {self._z_buffer}  (press Enter)"
        self._z_ax.text(0, 0.5, msg, va='center', ha='left', fontsize=10)
        self.fig.canvas.draw_idle()

    def _submit_z_from_buffer(self):
        try: z_val = float(self._z_buffer) if self._z_buffer.strip() != "" else 0.0
        except Exception: z_val = 0.0
        self.poly_z.append(z_val)
        self.awaiting_z = False
        self._z_ax.cla(); self._z_ax.set_axis_off()
        self._update_status(); self.fig.canvas.draw_idle()
        self.compute_and_plot_z()

    # ---------- Polygon ops ---------- #

    def close_current_polygon(self):
        if len(self.current_pts) < 3: return
        if self.current_pts[0] != self.current_pts[-1]:
            self.current_pts.append(self.current_pts[0])
        self.polygons.append(self.current_pts.copy())
        patch = MplPolygon(self.current_pts, closed=True, fill=False, edgecolor='C1', lw=2)
        self.poly_patches.append(patch); self.ax.add_patch(patch)
        self.current_pts.clear(); self.current_line.set_data([], [])
        self._update_status(); self.fig.canvas.draw_idle()
        self._open_z_prompt()

    def delete_last_polygon(self):
        if not self.polygons: return
        self.polygons.pop()
        if self.poly_patches: self.poly_patches.pop().remove()
        if self.poly_z: self.poly_z.pop()
        self._update_status(); self.fig.canvas.draw_idle()
        self.compute_and_plot_z()

    # ---------- Grid / Interpolation / Plot ---------- #

    def compute_Z(self):
        xs, ys, X, Y = self.xs, self.ys, self.Xg, self.Yg
        Z = np.zeros_like(X)
        if len(self.polygons) == 0:
            return xs, ys, Z
        polys = self.polygons  # user order
        zvals = [self.poly_z[i] if i < len(self.poly_z) else 0.0 for i in range(len(polys))]
        inside = [point_in_poly_grid(xs, ys, poly) for poly in polys]
        dists  = [dist_to_poly_boundary(X, Y, poly) for poly in polys]
        for i in range(len(polys)):
            mask_i = inside[i]
            if i < len(polys) - 1:
                mask_next = inside[i+1]
                # inside i but outside next → blend
                ring = mask_i & (~mask_next)
                if ring.any():
                    d_outer = dists[i][ring]
                    d_inner = dists[i+1][ring]
                    t = d_outer / (d_outer + d_inner + 1e-12)
                    Z[ring] = (1.0 - t) * zvals[i] + t * zvals[i+1]
                # inside i AND inside next → Z_i
                nested = mask_i & mask_next
                if nested.any():
                    Z[nested] = zvals[i]
            else:
                # last polygon: inside → Z_last
                if mask_i.any():
                    Z[mask_i] = zvals[i]
        return xs, ys, Z

    def compute_and_plot_z(self):
        xs, ys, Z = self.compute_Z()
        extent = [xs.min(), xs.max(), ys.min(), ys.max()]

        # 2D heatmap (single colorbar)
        if self.im is None:
            self.im = self.ax.imshow(Z, origin='lower', extent=extent, interpolation='nearest', cmap='terrain')
            vmin, vmax = float(np.nanmin(Z)), float(np.nanmax(Z)); vmax = vmax if vmax != 0 else 1.0
            self.im.set_norm(Normalize(vmin=vmin, vmax=vmax))
            self.cbar = self.fig.colorbar(self.im, ax=self.ax, label='Z')
        else:
            self.im.set_data(Z); self.im.set_extent(extent)
            vmin, vmax = float(np.nanmin(Z)), float(np.nanmax(Z)); vmax = vmax if vmax != 0 else 1.0
            self.im.set_norm(Normalize(vmin=vmin, vmax=vmax))
            if self.cbar: self.cbar.update_normal(self.im)
        self.fig.canvas.draw_idle()

        # 3D surface (single colorbar; Z positive upwards; tight layout)
        X, Y = np.meshgrid(xs, ys)
        self.ax3d.cla()
        surf = self.ax3d.plot_surface(X, Y, Z, cmap='terrain', linewidth=0, antialiased=False, rstride=1, cstride=1)
        self.ax3d.set_xlabel("x (km)"); self.ax3d.set_ylabel("y (km)"); self.ax3d.set_zlabel("Z")
        # Axis limits + box aspect to reduce empty space
        self.ax3d.set_xlim(xs.min(), xs.max())
        self.ax3d.set_ylim(ys.min(), ys.max())
        zmin, zmax = float(np.nanmin(Z)), float(np.nanmax(Z))
        if zmax == zmin: zmax = zmin + 1.0
        self.ax3d.set_zlim(zmin, zmax)
        try:
            self.ax3d.zaxis.set_inverted(False)  # ensure positive-upwards
            xr = xs.max()-xs.min(); yr = ys.max()-ys.min(); zr = max(zmax - zmin, 1e-9)
            self.ax3d.set_box_aspect((xr, yr, zr))
        except Exception:
            pass
        self.ax3d.view_init(elev=self.elev, azim=self.azim)

        # Single 3D colorbar
        if self.cbar3d is not None:
            try: self.cbar3d.remove()
            except Exception: pass
            self.cbar3d = None
        self.cbar3d = self.fig3d.colorbar(surf, ax=self.ax3d, shrink=0.72, aspect=18, pad=0.02, label='Z')

        self.fig3d.canvas.draw_idle()

    # ---------- Save (Z line + coords, 6 sig figs) ---------- #

    def save_and_exit(self):
        out_path = self.out_dir / f"{self.project}.UNIT"
        write_unit_with_z_and_coords_km(out_path, self.polygons, self.poly_z)
        print(f"Saved {len(self.polygons)} polygons to {out_path} (units: km; Z line + coords, 6 sig figs)")
        plt.close(self.fig); plt.close(self.fig3d)

# ---------------- File I/O ---------------- #

def write_unit_with_z_and_coords_km(path: Path, polygons: List[List[Tuple[float,float]]], poly_z: List[float]) -> None:
    """
    Output (km):
      # UNIT polygons generated <timestamp>
      <Z_value>              <-- number only, no leading 'Z'
      x y                    <-- 6 significant digits
      ...
    (blank line between polygons)
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(path, "w") as f:
        f.write(f"# UNIT polygons generated {ts}\n")
        for i, poly in enumerate(polygons):
            z_val = poly_z[i] if i < len(poly_z) else 0.0
            f.write(f"{z_val:.6g}\n")  # Z line, no 'Z ' prefix
            for x, y in poly:
                f.write(f"{x:.6g} {y:.6g}\n")
            if i != len(polygons) - 1:
                f.write("\n")

# ---------------- CLI ---------------- #

def main():
    ap = argparse.ArgumentParser(description="Draw polygons for TISC UNIT files (km units) + live 3D")
    ap.add_argument("project_or_prm", help="Project name (expects <name>.PRM) or a path to a .PRM file")
    ap.add_argument("--outdir", default=".", help="Output directory (default: current folder)")
    args = ap.parse_args()

    prm_input = Path(args.project_or_prm)
    if prm_input.suffix.lower() == ".prm" and prm_input.exists():
        prm_path = prm_input; project = prm_path.stem
    else:
        project = args.project_or_prm; prm_path = Path(f"{project}.PRM")

    if not prm_path.exists():
        print(f"Warning: PRM file not found: {prm_path} — using default km view.", file=sys.stderr)
        dom_m = {"xmin": -1_000.0, "xmax": 1_000.0, "ymin": -1_000.0, "ymax": 1_000.0, "Nx": 81, "Ny": 81}
    else:
        dom_m = parse_prm(prm_path)
        if dom_m.get("Nx") is None: dom_m["Nx"], dom_m["Ny"] = 81, 81

    km = lambda v: None if v is None else v/1000.0
    print(f"Loaded domain from {prm_path if prm_path.exists() else '(none)'} (km): "
          f"x in [{km(dom_m['xmin'])}, {km(dom_m['xmax'])}], y in [{km(dom_m['ymin'])}, {km(dom_m['ymax'])}] | "
          f"Nx={dom_m.get('Nx')} Ny={dom_m.get('Ny')}")

    app = PolyApp(dom_m, project, Path(args.outdir))
    plt.show()

if __name__ == "__main__":
    main()
