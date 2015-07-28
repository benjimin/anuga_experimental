"""

This scenario is a floodplain containing an array of buildings to induce anisotropic flow.

This test is based from (and aims to replicate):

Sanders, Schubert, Gallegos, Integral formulation of shallow-water equations 
with anisotropic porosity for urban flood modeling, J.Hydrol. (2008) 362, 19-38.


"""

# configurable parameters

method = 'DE1' # method could be DE1 or subgrid.

mesh_size = None

breaklines = None # whether or not to place break lines around buildings

building_height = 2.5 # metres (should exceed 1.15m)

display_figure = True # whether to show a diagram of the scenario
display_points = 15000





# scenario constants:  (refer to pg.30-31 of citation, do not change)

floodplain_size = [1500,500] # metres
floodplain_slope = [-0.0001,0] # constant slope is entirely in longitudinal direction
manning_coefficient = 0.035 # m^(-1/3) s
building_array = (5,5) # number of buildings in each direction
building_size = [10,20] # m
building_spacing = [30,30] # m (periodicity; subtract size to obtain gaps) 
building_centroid = [750,250] #m
building_array_angle = 45 # degrees (anti-clockwise)
discharge = 400 # cubic metres per second at upstream (left) boundary
simulation_time = 4000 # seconds (i.e. to approximately reach steady state)


"""

Need to transfer the urban centre (i.e. a collection of polygons) onto the
(vectorised) elevation function. That is, given any list of point
coordinates, we want to test each for overlap with the polygons.
This is basically a problem from vector image drawing, or from GIS.
A simpler (raster) image drawing problem would apply if filling values
directly onto a grid rather than workding with coordinates of the 
unstructured mesh.

Note: older Matplotlib versions (before 1.3) had a bug in "contains_points"
which could be worked-around using "contains_point" instead.

"""

import numpy as np

# convert some lists to vectors for convenience
building_array = np.array(building_array)
building_size = np.array(building_size)
building_spacing = np.array(building_spacing)
building_centroid = np.array(building_centroid)
floodplain_slope = np.array(floodplain_slope)

# use libraries for polygons
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
  
# Construct one building
dx,dy = 0.5*building_size
polygon_vertices = [(-dx,dy),(dx,dy),(dx,-dy),(-dx,-dy),(-dx,dy)]
building = Path(polygon_vertices,closed=True)
  
# Calculate the position of each building (relative to the town centroid)
izero = 0.5*(building_array-1)
locations = (building_spacing*(indices-izero) for indices in np.ndindex(tuple(building_array)))
  
# Construct the urban centre
translations = (Affine2D().translate(tx,ty) for tx,ty in locations)
buildings = [building.transformed(new_position) for new_position in translations]
town = Path.make_compound_path(*buildings)
  
# rotate anticlockwise then (chaining) translate urban centroid.
placement = Affine2D().rotate_deg(building_array_angle).translate(*building_centroid)
town = town.transformed(placement)
  
def elevation(x,y): # test: does this break if x,y are NOT vectors?  
  points = np.vstack((x,y)).T # convert pair-of-lists (input) to list-of-pairs
  
  indoor_boolean = town.contains_points(points)
  urban_thickness = np.where(indoor_boolean, building_height, 0)
  
  hillside = np.dot(points, floodplain_slope)
  
  return hillside + urban_thickness
  # convert boolean to elevations


if display_figure:
  import matplotlib.pyplot as plt 
  sx,sy = floodplain_size
  x = np.random.rand(display_points)*sx
  y = np.random.rand(display_points)*sy
  plt.figure().add_subplot(111,aspect='equal')
  plt.scatter(x,y,c=elevation(x,y),alpha=0.5,cmap='Paired',linewidth=0)
  plt.colorbar()
  plt.show()


"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig=plt.figure()
ax = fig.add_subplot(111)
patch=patches.PathPatch(building,facecolor='orange',lw=2)
ax.add_patch(patch)
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
plt.show()






def building(x,y):
    from numpy import piecewise as pw
    return 

# building -- a step function
# array -- sum (with translation)
# rotate -- this is just a linear change of variables.
# slope -- simplest to just add with building elevation

# test that building works:
assert building(0.9*building_size[0]/2.0,0.9*building_size[1]/2.0)==building_height, "Building interiors failed"
assert building(1.1*building_size[0]/2.0,1.1*building_size[1]/2.0)==0, "Building exteriors failed"


import unittest
class Tests(unittest.TestCase):
  def test_building(self):
    self.assertEqual(building(
    
    
    
x = np.array([32,0.5,-1])
y = np.array([12,0.5,-1])
z = np.vstack((x,y)).T

print elevation(x,y) + 1
    
def elevation(x,y):
  points = np.vstack((x,y)).T # convert pair-of-lists to list-of-pairs
  
  dx,dy = np.asarray(building_size)*0.5
  polygon_vertices = [(-dx,dy),(dx,dy),(dx,-dy),(-dx,-dy)]
  #polygon_codes = [path.Path.MOVETO] + [path.Path.LINETO]*3 + [path.Path.CLOSEPOLY]
  building = path.Path(polygon_vertices)#,codes=polygon_codes)
  
  #print polygon_vertices
  #print x
  
  #print points
  
  indoor = building.contains_points(points) # obtain boolean
  #indoor = building.contains_points(np.array([[32. ,  12. ], [  0.5  , 0.5], [ -1.  , -1. ]])) # obtain boolean
  return np.where(indoor,building_height,0.0) # convert boolean to elevations    
    
"""
