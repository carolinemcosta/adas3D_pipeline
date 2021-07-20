# adas3D_pipeline
Repository for pipeline prototype to create labeled volumetric meshes from Adas3D data. The meshes can be created from surfaces describing the different LV regions or surfaces plus scalar fields for different ventricular layers. The code reads the surfaces and tags as a point cloud with data, generates a tetrahedral mesh, and interpolates the point cloud data onto the mesh.

## Examples

A mesh can be generated using the surfaces defining scar, BZ, and healthy regions. In this case, the mesh tags are generated based on each surface. An example of usage is shown below. A coarse mesh (resolution=0.8mm), appropriate for Eikonal simulations is generated.

**Usage:**

`python3 adas_to_mesh.py --adas_folder="VTK_seperate_layers_3VTK" --mesh_res=0.8 --mesh_name="lv_mesh_eikonal"`

A mesh can be also generated using the LV surface layers defining with a scalar field representing the normalized LGE intensity. In this case, the mesh tags for  scar, BZ, and healthy are generated based set thresholds for BZ ans scar. An example of usage is shown below. A fine mesh (resolution=0.35mm), appropriate for monodomain simulations is generated.

**Usage:**

`python3 adas_to_mesh.py --adas_folder="VTK_seperate_layers" --scar_threshold=0.8 --bz_threshold=0.6 --mesh_res=0.35 --mesh_name="lv_mesh_monodomain"`

Example mesh generated from an Adas3D dataset using threshold method. Blue is healthy (tag=1), light read is border zone (tag=2), and dark red is scar (tag=3)


![mesh800um](https://user-images.githubusercontent.com/81109384/126347516-59d72b2f-9444-412f-b035-765e9ac31e32.png)


## Requirements
- Python 3.6
- Numpy 1.18.5
- Meshtool (https://bitbucket.org/aneic/meshtool/src/master/)
- PyVista (https://docs.pyvista.org/)

