# tAo Limitations

* **Discrete Movement**: Blocks move only an entire number of nodes at each time step. Block deformation during motion follows a vertical shear approach (preserving vertical thickness). Rotation of Blocks is not allowed.
* **Crust/Mantle Deformation**: tAo doesn't account for deformation in lower crust or mantle, and the deformation in the upper crust (by shortening) is kinematically introduced in the model. Subduction is not considered either. Different velocities of shortening for each detachment fault are (still) not implemented.
* **Compaction**: Gradual sediment compaction is not considered, but this shouldn't affect the sediment onlap/offlap geometry in the external side of the foreland basin, but only the thickness of the Blocks.
* **Facies**: Grain size or sediment facies distribution of the deposits are not differentiated.
