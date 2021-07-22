import argparse
import numpy as np
import os
import pyvista as pv
import glob

import meshtool_cmd as mc


def read_vtk_polydata(mesh_name):
  mesh = pv.read(mesh_name)
  pts = mesh.points
  data = mesh.active_scalars
  
  #print(mesh_name)
  #print("Pts size:", pts)
  #if not type(data) == None:
    #print("Data size:", data)
  
  return pts, data


def write_point_cloud_data(full_pts_name, full_tags_name, point_cloud, tags_data):
  npoints = len(point_cloud)
  fp = open(full_pts_name, 'w')
  fd = open(full_tags_name, 'w')

  fp.write("%d\n"%(npoints-1))
  for i in range(1,npoints): # skip first line of concatated array
    fp.write("%f %f %f\n"%(point_cloud[i][0],point_cloud[i][1],point_cloud[i][2]))
    fd.write("%d\n"%(tags_data[i]))

  fp.close()
  fd.close()


def get_excluded_ids(region_pts, excluded_pts):
  # find closest region point and change tag - this is to remove the "excluded layer" from other regions
  npts = len(excluded_pts)
  ids = np.ones(npts, dtype=int)
  for i in range(npts):    
    ids[i] = np.sum((region_pts - excluded_pts[i,:])**2, axis=1).argmin() # min distance
  
  return ids

    
def get_layer(layer_name, excluded_name, point_cloud, tags_data, excluded_tag, region_tags=[], bz_scar_thresholds=[]):

  # read surface points
  layer_pts, layer_data = read_vtk_polydata(layer_name+".vtk")
  excluded_pts, _ = read_vtk_polydata(excluded_name+".vtk")

  layer_tags = np.ones(len(layer_pts), dtype=int)*region_tags[0] # healthy or region tag

  # using LV layers instead of scar and BZ surfaces
  if len(bz_scar_thresholds) > 0:
    # apply BZ threshold
    bz = layer_data > bz_scar_thresholds[0]
    layer_tags[bz] = region_tags[1] # BZ tag
    
    # apply scar threshold
    scar = layer_data > bz_scar_thresholds[1]
    layer_tags[scar] = region_tags[2] # scar tag

  # find change tag in excluded regions
  ids = get_excluded_ids(layer_pts, excluded_pts)
  layer_tags[ids] = excluded_tag
    
  # concatenate points and tags
  point_cloud = np.concatenate((point_cloud, layer_pts), axis=0)
  tags_data = np.concatenate((tags_data.T,layer_tags), axis=0)
  
  return point_cloud, tags_data


def get_point_cloud_tags_from_regions(mesh_dir, full_pts_name, full_tags_name, healthy_tag, bz_tag, core_tag, excluded_tag):  
  if not os.path.isfile(full_pts_name):
    print("Creating point cloud and tags using region layers...")
    
    point_cloud = np.empty([1,3])
    tags_data = np.empty([1,])
    
    # read in all surfaces as points clouds and concatate points and tags, excluding "excluded" surfaces
    files_list = os.listdir(mesh_dir)

    for file_name in files_list:
      if file_name.startswith("Layer"):
        layer = file_name.split("Layer")[1]
        # excluded
        excluded_name = "%s/Layer%s/EX_%s"%(mesh_dir, layer, layer)

        # healthy
        healthy_name = "%s/Layer%s/HE_%s"%(mesh_dir, layer, layer)
        if os.path.isfile(healthy_name+".vtk") and os.path.isfile(excluded_name+".vtk"):
          #get_layer(layer_name, excluded_name, point_cloud, tags_data, excluded_tag, region_tags=[], bz_scar_thresholds=[])
          point_cloud, tags_data = get_layer(healthy_name, excluded_name, point_cloud, tags_data, excluded_tag, [healthy_tag])
          
        # BZ
        bz_name = "%s/Layer%s/BZ_%s"%(mesh_dir, layer, layer)
        if os.path.isfile(bz_name+".vtk") and os.path.isfile(excluded_name+".vtk"):
          point_cloud, tags_data = get_layer(bz_name, excluded_name, point_cloud, tags_data, excluded_tag, [bz_tag])
          
        # core
        core_name = "%s/Layer%s/CO_%s"%(mesh_dir, layer, layer)
        if os.path.isfile(core_name+".vtk") and os.path.isfile(excluded_name+".vtk"):
          point_cloud, tags_data = get_layer(core_name, excluded_name, point_cloud, tags_data, excluded_tag, [core_tag])

    # write out full point cloud and data 
    write_point_cloud_data(full_pts_name, full_tags_name, point_cloud, tags_data)
    

def get_point_cloud_tags_from_lv(mesh_dir, full_pts_name, full_tags_name, healthy_tag, bz_tag, core_tag, excluded_tag, bz_scar_thresholds):  
  if not os.path.isfile(full_pts_name):
    print("Creating point cloud and tags using LV layers ...")
    
    point_cloud = np.empty([1,3])
    tags_data = np.empty([1,])
    
    tags_list = [healthy_tag, bz_tag, core_tag]
    
    files_list = os.listdir(mesh_dir)

    for file_name in files_list:      
      if file_name.endswith("_Excluded-DE-MRI.vtk"):
        print(file_name)
        base_name = file_name.split("_Excluded")[0]
        layer_name = mesh_dir + "/" + base_name + "-DE-MRI"
        print(layer_name)
        excluded_name = mesh_dir + "/" + base_name + "_Excluded-DE-MRI"
        print(excluded_name)
      
        if os.path.isfile(layer_name+".vtk"):
          point_cloud, tags_data = get_layer(layer_name, excluded_name, point_cloud, tags_data, excluded_tag, tags_list, bz_scar_thresholds)          

    # write out full point cloud and data 
    write_point_cloud_data(full_pts_name, full_tags_name, point_cloud, tags_data)

    
def assign_int_tags(full_tags_name_elem, full_tags_name_elem_int, lv_mesh_name, lv_mesh_name_tag):
  # assign tags to new mesh
  if not os.path.isfile(lv_mesh_name_tag+".elem"):
    # round tag to nearest integer
    if not os.path.isfile(full_tags_name_elem_int):
      print("Formatting tags...")    
      tag_data = np.loadtxt(full_tags_name_elem, dtype=float)
      tag_data_int = np.rint(tag_data)    
      np.savetxt(full_tags_name_elem_int, tag_data_int, fmt='%d')  
    else:
      tag_data_int = np.loadtxt(full_tags_name_elem_int, dtype=int)
      
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
  else:
    print("%s already exists. No tags assigned."%lv_mesh_name_tag)
  
  
def main(args):
  
  # set arguments
  mesh_dir = args.adas_folder[0]
  if not os.path.isdir(mesh_dir):
    raise NameError("ADAS folder not found")

  scar_threshold = args.scar_threshold
  bz_threshold = args.bz_threshold
  if bz_threshold > scar_threshold:
    raise ValueError("The border zone thrshold must be smaller than the scar threshold")

  meshtool_bin = args.meshtool_bin   
  if not os.path.isdir(mesh_dir):
    raise NameError("Cannot find meshtool. If binary not \"meshtool\", please define \"meshtool_bin\" parameter")
  
  mesh_res = args.mesh_res
  base_name = args.mesh_name[0]

  # define file names
  lv_mesh_name = mesh_dir + "/" + base_name
  full_pts_name = "%s/point_cloud_test.pts"%(mesh_dir)
  full_tags_name = "%s/tags_test.dat"%(mesh_dir)

  if scar_threshold and bz_threshold:
    lv_surf_name = glob.glob("%s/*Left*Ventricle-DE-MRI.vtk"%mesh_dir)[0]
  else:
    lv_surf_name = "%s/Other/LV_M-DE-MRI.vtk"%mesh_dir  

  if not os.path.isfile(lv_surf_name):
    raise NameError("LV surface file not found")
  
  # define tags
  healthy_tag = 1
  bz_tag = 2
  core_tag = 3
  excluded_tag = healthy_tag # same as healthy

  # create point cloud with all layers and generate tags
  if scar_threshold and bz_threshold: # Using LV layer and thresholds
    get_point_cloud_tags_from_lv(mesh_dir, full_pts_name, full_tags_name, healthy_tag, bz_tag, core_tag, excluded_tag, [bz_threshold, scar_threshold])
  else: # Using regions surfaces
    get_point_cloud_tags_from_regions(mesh_dir, full_pts_name, full_tags_name, healthy_tag, bz_tag, core_tag, excluded_tag)  

  # generate tet mesh from LV surface
  mc.surf2mesh(meshtool_bin, lv_surf_name, lv_mesh_name)  
    
  # refine mesh for Eikonal or monodomain simulation
  lv_mesh_name_ref = lv_mesh_name + "_%dum"%int(mesh_res*1000)
  mc.refine_mesh(meshtool_bin, lv_mesh_name, mesh_res, lv_mesh_name_ref)

  # interpolate point cloud data onto tet mesh
  full_tags_name_pts = lv_mesh_name_ref + "_tags_pts.dat"
  mc.interpolate_cloud(meshtool_bin, lv_mesh_name_ref, full_pts_name, full_tags_name, full_tags_name_pts)
  
  # interpolate tags from nodes to elements
  full_tags_name_elem = lv_mesh_name_ref + "_tags_elem.dat"
  mc.node2elem(meshtool_bin, lv_mesh_name_ref, full_tags_name_pts, full_tags_name_elem)
  
  # round tags to nearest integer and assign to mesh
  full_tags_name_elem_int = lv_mesh_name_ref + "_tags_elem_int.dat"
  lv_mesh_name_tag = lv_mesh_name_ref + "_tagged"
  assign_int_tags(full_tags_name_elem, full_tags_name_elem_int, lv_mesh_name_ref, lv_mesh_name_tag)
  
  # convert from mm to um
  mc.to_mm(meshtool_bin, lv_mesh_name_tag)
  
      
if __name__== "__main__":
  
  parser = argparse.ArgumentParser()
  parser.formtter_class = argparse.ArgumentDefaultsHelpFormatter
  
  parser.add_argument('--adas_folder', type=str, default=None, nargs='+', help='Provide the directory for the Adas3D data.')
  parser.add_argument('--mesh_name', type=str, default="lv_mesh", nargs='+', help='Define name of output mesh.')
  parser.add_argument('--mesh_res', type=float, default=0.8, nargs='?', help='Mean resolution to refine final mesh.')
  parser.add_argument('--meshtool_bin', type=str, default="meshtool", nargs='?', help='Provide name of meshtool binary')
  parser.add_argument('--scar_threshold', type=float, default=None, nargs='?', help='For a threshold-based mesh define the scar threshold')
  parser.add_argument('--bz_threshold', type=float, default=None, nargs='?', help='For a threshold-based mesh define the bz threshold')
  
  args = parser.parse_args()
  
  main(args)
  
