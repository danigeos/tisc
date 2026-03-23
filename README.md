# TISC & tAo: Coupled Tectonic and Surface Processes Models

This repository contains the source code for two closely related scientific modeling softwares: **TISC** and **tAo**. They are managed together in this monorepo to facilitate the sharing of common libraries and ensure consistency.

![Demo](./images/demo_LagoMare.gif)

---

### TISC: 2D (planform modeling, pseudo 3D)

TISC (Tectonics, Isostasy, Surface Processes, and Climate) is a 2D/3D planform model for simulating the interaction between tectonic processes, lithospheric flexure, and surface processes like erosion and sedimentation.

*   **Source Code:** `./tisc/`
*   **Documentation:** `./tisc/doc/`

### tAo: (vertical cross-section modeling, pseudo 2D)

tAo is a 1D (pseudo-2D) model for simulating the evolution of a vertical cross-section of the lithosphere, including flexure, faulting, and surface processes. It is ideal for modeling foreland basin systems.

*   **Source Code:** `./tao/`
*   **Documentation:** `./tao/doc/`

---

## Compilation

To compile both projects, simply run `make` from this root directory:

```sh
make
```

This will place all executables in the root `bin/` directory.