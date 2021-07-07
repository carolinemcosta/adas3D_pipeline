# adas3D_pipeline
Repository for pipeline prototype to create labeled volumetric mesh from Adas3D data. Reads in surfaces describing the different LV regions, generates a point cloud with these, generates a tetrahedral mesh, and interpolates the point cloud data onto the mesh.

## Parameters
Path to data folder containing "VTK_seperate_layers_3VTK/" folder with layers 10 to 90. Each layer subfolder must contain "BZ", "CO", "HE", and "EX" VTK polydata files. These files represent surfaces for border zone, scar core, healthy tissue, and excluded regions. The folder must also contain a subfolder "Other" with the LV surface mesh called "LV_M-DE-MRI.vtk".

Also needs the meshtool binary name.

## Requirements
Python 3.6
Numpy 1.18.5
Meshtool (https://bitbucket.org/aneic/meshtool/src/master/)
