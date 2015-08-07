"""

This scenario is a floodplain containing an array of buildings to induce anisotropy in the flow.

This test is based from (and aims to replicate):

Sanders, Schubert, Gallegos, Integral formulation of shallow-water equations 
with anisotropic porosity for urban flood modeling, J.Hydrol. (2008) 362, 19-38.


"""

# testing options:

display_figure = True # whether to show a diagram of the scenario
display_points = 50000 # number of elevation samples



# configurable parameters of the solver: 

method = 'DE1' # method might be DE1 or subgrid.

mesh_resolution = 0.5*(100)**2   # m^2, minimum triangle area

breaklines_around_buildings = True # Boolean, whether the mesh should conform with urban planning




# extra scenario constants:

building_height = 2.5 # metres (should exceed 1.15m)

cliff_fall = -10 # metres (substitute for over-fall boundary-condition)
bucket_width = 50 # m (extends right edge of domain)



# scenario constants:  (refer to pg.30-31 of citation, do not change)

floodplain_size = [1500,500] # metres
floodplain_slope = [-0.001,0] # constant slope gradient is entirely in longitudinal direction
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
which could be worked-around using "contains_point" or by breaking compound
paths into polygons defined without explicit closure.

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
  
# rotate anticlockwise then (chaining) shift urban centroid into position.
placement = Affine2D().rotate_deg(building_array_angle).translate(*building_centroid)
town = town.transformed(placement)
  
def elevation(x,y): 
  points = np.vstack((x,y)).T # convert pair-of-lists (input) to list-of-pairs
  
  indoor_boolean = town.contains_points(points)
  urban_thickness = np.where(indoor_boolean, building_height, 0)
    
  hillside = np.dot(points, floodplain_slope)
  
  cliff = np.where(x > floodplain_size[0], cliff_fall, 0) # fall past right side of plain
  
  return hillside + urban_thickness + cliff
  

"""

Breaklines:

Sanders, et al, used a free-overfall boundary-condition at the downstream (right) boundary. 
For the purposes of testing solver algorithms in which not all boundary types have been
implemented yet (e.g. using only reflective boundary walls), we can simulate free-overfall 
by implementing a cliff edge (and downstream bucket). We can improve this simulation by
aligning the mesh with the cliff edge, i.e., by using anuga's riverwall or breakline feature.

Moreover, we may wish to align mesh to building edges, at least for comparison.

"""

# create cliff-edge polyline
cx,cy = floodplain_size
cliff = [(cx,0),(cx,cy)]

breaklines = [cliff]

if breaklines_around_buildings:
  breaklines += town.to_polygons() # Note all already explicitly closed
  

"""

Mesh generation.

"""

import anuga

# implicitly-closed polygon enclosing the floodplain plus the basin below the cliff
floodplain =  [(0,0),(0,cy),(cx,cy),(cx,0)]
cx += bucket_width
full_extent = [(0,0),(0,cy),(cx,cy),(cx,0)]

# TODO: ought to use interior regions to reduce bucket resolution


meshname = 'urban.msh' # if None then anuga won't generate_mesh

mesh = anuga.create_mesh_from_regions(full_extent,boundary_tags={'outer boundary':[0,1,2,3]},
  maximum_triangle_area=mesh_resolution, breaklines=breaklines, filename=meshname)
  
"""

Domain creation

Could start out with the floodplain dry (stage equal to elevation perhaps),
although it might converge faster if adding the steady-state water depth
(known analytically if we can neglect the buildings). 
Alternatively, could set a constant value (to check what happens when the water reaches
a bucket that isn't quite empy).



"""

domain = anuga.create_domain_from_file(meshname)

domain.set_flow_algorithm(method)

domain.set_quantity('stage',cliff_fall) # start out dry? Oh wait..
domain.set_quantity('friction',manning_coefficient)
domain.set_quantity('elevation',elevation)

domain.set_name('urban') # .sww output file
domain.set_datadir('.') # current working directory

domain.set_boundary({'outer boundary': anuga.Reflective_boundary(domain)}) # hard walls

from anuga import Inlet_operator
left_edge = [(0.0,cy),(0,0)]
Inlet_operator(domain,left_edge,discharge) # set up constant water inflow

"""

Run the simulation

"""

t1 = simulation_time 
t1 = 400

for t in domain.evolve(yieldstep=20,finaltime=t1):
  print domain.timestepping_statistics()



# Although the output is already saved (in .sww using NetCDF) there is the option to 
# reformat as .tif for a quick check (using "display" on UNIX) of the final status:
#from anuga import plot_utils as util
#util.Make_Geotif('urban.sww',output_quantities=['depth'],myTimeStep='last',
#  CellSize=1.0,EPSG_CODE=32756,bounding_polygon=full_extent,k_nearest_neighbours=1)

"""

Can compare plot with fig. 8 of Sanders et al. to check everything looks sensible

"""

if display_figure:
  import matplotlib.pyplot as plt
  from matplotlib.patches import PathPatch
  fig = plt.figure()
  padding = 30
  xl,yl = [-padding,cx+padding],[-padding,cy+padding] # plot extent
  
  
  # first subplot tests the setup prior to the domain existing.
  ax1 = fig.add_subplot(311,aspect='equal') # set aspect ratio
  plt.title('Set-up')
  plt.xlim(xl)
  plt.ylim(yl)
  
  # randomly test elevation function
  sx,sy = floodplain_size
  sx += bucket_width
  x = np.random.rand(display_points)*sx
  y = np.random.rand(display_points)*sy  
  plt.scatter(x,y,c=elevation(x,y),alpha=0.5,cmap='Paired',linewidth=0, marker='.')
  plt.colorbar() # elevation is mapped to colour
  
  # check the mesh looks ok
  pts = mesh.getMeshVertices()
  paths = (Path([pts[i],pts[j],pts[k],pts[i]]) for i,j,k in mesh.getTriangulation())
  ax1.add_patch(PathPatch(Path.make_compound_path(*paths),fill=None,alpha=0.2))
  
  # check the buildings (and cliff) align with expectations
  shapes = Path.make_compound_path(town,Path(cliff)) # town + cliff
  ax1.add_patch(PathPatch(shapes,fill=None,linewidth=0.5,color='black',alpha=0.8,linestyle='dashed'))
  



  # second subplot tests the results from the final time-step
  ax2 = fig.add_subplot(312,aspect='equal') # set aspect ratio
  plt.title('Results')
  plt.xlim(xl)
  plt.ylim(yl)
  
  # obtain water depth
  water_depth = domain.get_quantity('height')
  x,y,values,triangles = water_depth.get_vertex_values(smooth=False)
  
  # interpolate separately within different triangles
  plt.tripcolor(x,y,triangles,values,shading='gouraud',alpha=0.8)
  plt.colorbar()
  
  # get centroid coordinates
  centx,centy = np.hsplit(domain.get_centroid_coordinates(),2)
  
  u = domain.get_quantity('xmomentum').get_values(location='centroids')
  v = domain.get_quantity('ymomentum').get_values(location='centroids')
  
  plt.quiver(centx,centy,u,v,alpha=0.6)
  


  # why does interpolation become unstable for large (not small) time periods?


  # third subplot 
  ax2 = fig.add_subplot(313,aspect='equal') # set aspect ratio
  plt.title('Interpolation')
  plt.xlim(xl)
  plt.ylim(yl)
  
  
  # interpolate water depth, from centroid values, using radial basis functions
  water_depth = domain.get_quantity('stage') - domain.get_quantity('elevation')
  from scipy.interpolate import Rbf
  depth = Rbf(centx,centy,water_depth.get_values(location='centroids'))
  
  # produce image of smoothly-interpolated water depth
  i = 5
  X,Y = np.meshgrid(np.arange(0,sx,i),np.arange(0,sy,i))
  Z = depth(X,Y)
  plt.pcolor(X,Y,Z)
  plt.colorbar()
  
  # lightly indicate centroid positions
  plt.scatter(centx,centy,s=5,marker='x',alpha=0.2)
  
  # interpolate the velocities also.
  # note: 

  xvel = domain.get_quantity('xmomentum')/water_depth
  yvel = domain.get_quantity('ymomentum')/water_depth
  u = xvel.get_values(location='centroids')
  v = yvel.get_values(location='centroids')
  
  i = 20 # interval at which to place quivers
  X,Y = np.meshgrid(np.arange(0,sx,i),np.arange(0,sy,i)) # grid for placing quivers
  U = Rbf(centx,centy,u)
  V = Rbf(centx,centy,v) # interpolate onto grid using radial basis functions
  plt.quiver(X,Y,U(X,Y),V(X,Y),alpha=0.6)
  
  # would probably prefer a streamline plot?
  

  
  plt.show()
  


  """
  # interpolate separately within different triangles
  plt.tripcolor(x,y,triangles,values,alpha=0.8)
  plt.colorbar()
  
  # keep mesh visible
  plt.triplot(x,y,triangles,'-',color='k',alpha=0.1) 

  """

  """
  # keep mesh visible
  plt.triplot(x,y,triangles,'-',color='k',alpha=0.1) 
  
  # quiver plot of flow field
  
  # obtain water velocity
  xvel = domain.get_quantity('xmomentum')/water_depth
  yvel = domain.get_quantity('ymomentum')/water_depth
  u = xvel.get_values()
  v = xvel.get_values()
  
  # need to get both components and smoothly interpolate first.
  
  i = 20 # interval at which to place quivers
  X,Y = np.meshgrid(np.arange(0,sx,i),np.arange(0,sy,i)) # grid for placing quivers
  U = Rbf(x,y,u)(X,Y)
  V = Rbf(x,y,v)(X,Y) # interpolate onto grid using radial basis functions
  plt.quiver(X,Y,U,V)

  
  
    
  """
