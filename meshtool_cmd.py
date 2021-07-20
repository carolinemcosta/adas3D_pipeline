import os

def surf2mesh(meshtool_bin, surf_name, mesh_name):
  # generate tet mesh from LV surface
  if not os.path.isfile(mesh_name+".elem"):
    cmd = "%s generate mesh -surf=%s -ifmt=vtk -outmsh=%s -ofmt=carp_txt"%(meshtool_bin, surf_name, mesh_name)
    print(cmd)
    os.system(cmd)
  else:
    print("%s already exists. Mesh not generated."%mesh_name)
  

def refine_mesh(meshtool_bin, mesh_name, mesh_res, mesh_name_ref):
  # refine mesh for Eikonal or monodomain simulation
  if not os.path.isfile(mesh_name_ref+".elem"):
    cmd = "%s resample mesh -msh=%s -avrg=%1.2f -outmsh=%s -postsmth=0 -ofmt=carp_txt"%(meshtool_bin, mesh_name, mesh_res, mesh_name_ref)
    print(cmd)
    os.system(cmd)
  else:
    print("%s already exists. No mesh refinement done."%mesh_name_ref)    
    

def interpolate_cloud(meshtool_bin, mesh_name, pts_name, tags_name, out_tags_name):
  # interpolate point cloud data onto tet mesh
  if not os.path.isfile(out_tags_name):
    cmd = "%s interpolate clouddata -omsh=%s -pts=%s -idat=%s -odat=%s -mode=1"%(meshtool_bin, mesh_name, pts_name, tags_name, out_tags_name)
    print(cmd)
    os.system(cmd)
  else:
    print("%s already exists. No interpolation done."%out_tags_name)
    

def node2elem(meshtool_bin, mesh_name, tags_name_pts, tags_name_elem):
  # interpolate tags from nodes to elements
  if not os.path.isfile(tags_name_elem):
    cmd = "%s interpolate node2elem -omsh=%s -idat=%s -odat=%s"%(meshtool_bin, mesh_name, tags_name_pts, tags_name_elem)
    print(cmd)
    os.system(cmd)
  else:
    print("%s already exists. No interpolation done."%tags_name_elem)
