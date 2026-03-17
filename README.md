# tAo & TISC: Coupled Tectonic and Surface Process Models

This repository contains the source code for two closely related scientific modeling softwares: **tAo** and **TISC**. They are managed together in this monorepo to facilitate the sharing of common libraries and ensure consistency.

---

### TISC: 2D/3D Planform Modeling

TISC (Tectonics, Isostasy, Surface Processes, and Climate) is a 2D/3D planform model for simulating the interaction between tectonic processes, lithospheric flexure, and surface processes like erosion and sedimentation.

*   **Source Code:** `./tisc/`
*   **Documentation:** `./tisc/doc/`

### tAo: 1D/2D Vertical Cross-Section Modeling

tAo is a 1D (pseudo-2D) model for simulating the evolution of a vertical cross-section of the lithosphere, including flexure, faulting, and surface processes. It is ideal for modeling foreland basin systems.

*   **Source Code:** `./tao/`
*   **Documentation:** `./tao/doc/`

---

## Compilation

To compile both projects, simply run `make` from this root directory:

```sh
make
```

This will place all executables in the root `bin/` directory and all object files in the root `build/` directory.