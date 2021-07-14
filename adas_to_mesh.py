import argparse
import numpy as np
import os
import pyvista as pv

def read_vtk_polydata(mesh_name):
  mesh = pv.read(mesh_name+".vtk")
  data = mesh.active_scalars
  
  return data


def get_layer(meshtool_bin, region_name, excluded_name, point_cloud, tags_data, region_tag, excluded_tag):
  # convert to CARP format
  if not os.path.isfile(region_name+".pts"):
    cmd = "%s convert -imsh=%s.vtk -omsh=%s  -ofmt=carp_txt"%(meshtool_bin, region_name, region_name)
    os.system(cmd)
  
  if not os.path.isfile(excluded_name+".pts"):
    cmd = "%s convert -imsh=%s.vtk -omsh=%s  -ofmt=carp_txt"%(meshtool_bin, excluded_name, excluded_name)
    os.system(cmd)

  # read surface points
  region_pts = np.loadtxt(region_name+".pts", skiprows=1, dtype=float)
  excluded_pts = np.loadtxt(excluded_name+".pts", skiprows=1, dtype=float)
  
  # create tags 
  region_tags = np.ones(len(region_pts))*region_tag
  
  # find closest region point and change tag - this is to remove the "excluded layer" from other regions
  npts = len(excluded_pts)
  for i in range(npts):    
    idx = np.sum((region_pts - excluded_pts[i,:])**2, axis=1).argmin() # min distance
    region_tags[idx] = excluded_tag    
    
  # concatenate points and tags
  point_cloud = np.concatenate((point_cloud, region_pts), axis=0)
  tags_data = np.concatenate((tags_data.T,region_tags), axis=0)      
  
  return point_cloud, tags_data
    
def get_point_cloud_tags(meshtool_bin, mesh_dir, full_pts_name, full_tags_name, healthy_tag, bz_tag, core_tag, excluded_tag):
  
  if not os.path.isfile(full_pts_name):
    print("Creating point cloud and tags with all layers...")
    
    point_cloud = np.empty([1,3])
    tags_data = np.empty([1,])
    
    # read in all surfaces as points clouds and concatate points and tags, excluding "excluded" surfaces
    for layer in range(10,100,10):
      # excluded
      excluded_name = "%s/Layer%d/EX_%d"%(mesh_dir, layer, layer)

      # healthy
      healthy_name = "%s/Layer%d/HE_%d"%(mesh_dir, layer, layer)
      if os.path.isfile(healthy_name+".vtk") and os.path.isfile(excluded_name+".vtk"):
        point_cloud, tags_data = get_layer(meshtool_bin, healthy_name, excluded_name, point_cloud, tags_data, healthy_tag, excluded_tag)
        
      # BZ
      bz_name = "%s/Layer%d/BZ_%d"%(mesh_dir, layer, layer)
      if os.path.isfile(bz_name+".vtk") and os.path.isfile(excluded_name+".vtk"):
        point_cloud, tags_data = get_layer(meshtool_bin, bz_name, excluded_name, point_cloud, tags_data, bz_tag, excluded_tag)
        
      # core
      core_name = "%s/Layer%d/CO_%d"%(mesh_dir, layer, layer)
      if os.path.isfile(core_name+".vtk") and os.path.isfile(excluded_name+".vtk"):
        point_cloud, tags_data = get_layer(meshtool_bin, core_name, excluded_name, point_cloud, tags_data, core_tag, excluded_tag)


    # write out full point cloud and data 
    npoints = len(point_cloud)
    with open(full_pts_name, 'w') as f:
      f.write("%d\n"%(npoints-1))
      for i in range(1,npoints):
        f.write("%f %f %f\n"%(point_cloud[i][0],point_cloud[i][1],point_cloud[i][2]))

    with open(full_tags_name, 'w') as f:
      for i in range(1,npoints):
        f.write("%d\n"%(tags_data[i]))
      
    
def assign_int_tags(full_tags_name_elem, full_tags_name_elem_int, lv_mesh_name, lv_mesh_name_tag):
  # round tag to nearest integer
  if not os.path.isfile(full_tags_name_elem_int):
    print("Formatting tags...")    
    tag_data = np.loadtxt(full_tags_name_elem, dtype=float)
    tag_data_int = np.rint(tag_data)    
    np.savetxt(full_tags_name_elem_int, tag_data_int, fmt='%d')  
  else:
    tag_data_int = np.loadtxt(full_tags_name_elem_int, dtype=int)
  
  # assign tags to new mesh
  if not os.path.isfile(lv_mesh_name_tag+".elem"):
    print("Assigning tags...")
    elems = np.loadtxt(lv_mesh_name+".elem", skiprows=1, usecols=(1,2,3,4,5), dtype=int)  
    esize = len(elems)
    
    with open(lv_mesh_name_tag+".elem", "w") as f:
      f.write("%d\n"%esize)
      for i in range(esize):
        f.write("Tt %d %d %d %d %d\n"%(elems[i,0], elems[i,1], elems[i,2], elems[i,3], tag_data_int[i]))

    cmd="ln -s %s.lon %s.lon"%(lv_mesh_name, lv_mesh_name_tag)
    os.system(cmd)
    cmd="ln -s %s.pts %s.pts"%(lv_mesh_name, lv_mesh_name_tag)
    os.system(cmd)
  
  
def main(args):
  # set arguments
  mesh_dir = args.adas_folder + "/VTK_seperate_layers_3VTK/"
  meshtool_bin = args.meshtool_bin

  full_pts_name = "%s/point_cloud.pts"%(mesh_dir)
  full_tags_name = "%s/tags.dat"%(mesh_dir)
  
  healthy_tag = 1
  bz_tag = 2
  core_tag = 3
  excluded_tag = healthy_tag # same as healthy
  
  # create point cloud with all layers and generate tags
  get_point_cloud_tags(meshtool_bin, mesh_dir, full_pts_name, full_tags_name, healthy_tag, bz_tag, core_tag, excluded_tag)

  # generate tet mesh from LV surface
  lv_surf_name = "%s/Other/LV_M-DE-MRI.vtk"%mesh_dir
  lv_mesh_name = "%s/lv_mesh"%mesh_dir
  if not os.path.isfile(lv_mesh_name+".elem"):
    cmd = "%s generate mesh -surf=%s -ifmt=vtk -outmsh=%s -ofmt=carp_txt"%(meshtool_bin, lv_surf_name, lv_mesh_name)
    print(cmd)
    os.system(cmd)

  # interpolate point cloud data onto tet mesh
  full_tags_name_pts = "%s/lv_mesh_tags_pts.dat"%(mesh_dir)
  if not os.path.isfile(full_tags_name_pts):
    cmd = "%s interpolate clouddata -omsh=%s -pts=%s -idat=%s -odat=%s -mode=1"%(meshtool_bin, lv_mesh_name, full_pts_name, full_tags_name, full_tags_name_pts)
    print(cmd)
    os.system(cmd)
  
  # interpolate tags from nodes to elements
  full_tags_name_elem = "%s/lv_mesh_tags_elem.dat"%(mesh_dir)
  if not os.path.isfile(full_tags_name_elem):
    cmd = "%s interpolate node2elem -omsh=%s -idat=%s -odat=%s"%(meshtool_bin, lv_mesh_name, full_tags_name_pts, full_tags_name_elem)
    print(cmd)
    os.system(cmd)
  
  # round tags to nearest integer and assign to mesh
  full_tags_name_elem_int = "%s/lv_mesh_tags_elem_int.dat"%(mesh_dir)
  lv_mesh_name_tag = "%s/lv_mesh_tagged"%mesh_dir
  assign_int_tags(full_tags_name_elem, full_tags_name_elem_int, lv_mesh_name, lv_mesh_name_tag)
  
      
if __name__== "__main__":
  
  parser = argparse.ArgumentParser()
  parser.formtter_class = argparse.ArgumentDefaultsHelpFormatter
  
  parser.add_argument('--adas_folder', type=str, default=None, help='Provide the directory for the Adas3D data. Must contain the VTK_seperate_layers_3VTK folder')
  parser.add_argument('--meshtool_bin', type=str, default="meshtool", help='Provide name of meshtool binary')
  args = parser.parse_args()
  
  main(args)
  
