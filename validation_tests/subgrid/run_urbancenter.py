"""

This scenario is a floodplain containing an array of buildings to induce anisotropy in the flow.

This test is based from (and aims to replicate):

Sanders, Schubert, Gallegos, Integral formulation of shallow-water equations 
with anisotropic porosity for urban flood modeling, J.Hydrol. (2008) 362, 19-38.



Usage examples:

python run_urbancenter.py --res 80 --breaklines
python run_urbancenter.py --res 80 --breaklines --hires 10
python run_urbancenter.py --res 80 --subgrid
python 

"""

# testing options:

display_figure = False # whether to show a diagram of the scenario
display_points = 5000 # number of elevation samples



# configurable parameters of the solver: 


method = 'DE1' # method might be discontinuous-elevation (DE1) or subgrid (DE_SG).

mesh_resolution = 950 #0.5*(50)**2   # m^2, maximum triangle area

higher_resolution_near_centre = None # None, or value in m^2. Say, 10
higher_resolution_padding_width = 50 # m, width of padding of hi-res region

breaklines_around_buildings = False # Boolean, whether the mesh should conform with urban planning



# NB: Sanders et. al. use max area constraints of 80m2 or 950m2,
# resulting in 1257 to 14640 triangles for the plain.
# Elsewhere (fig.4) they consider 25, 500, or 10000 m2 with the same-size buildings.


# extra scenario constants:

building_height = 2.5 # metres (should exceed 1.15m)

cliff_fall = -10 # metres (substitute for over-fall boundary-condition)
bucket_width = 500 # m (extends right edge of domain)
# nb: must fit a couple million cubic metres, i.e. width*500m*fall >> discharge*time


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
Permit setting options for the solver from the command line,
to override defaults chosen above.

e.g. python run_urbancenter.py --res 50 --subgrid

This is to facilitate batch comparisons of solvers.

"""

import sys
import getopt
opts,args = getopt.getopt(sys.argv[1:], "", ["subgrid","res=","breaklines","hires="])
for opt,arg in opts:
  if opt == "--subgrid":
    method = "DE_SG"
  elif opt == "--res":
    mesh_resolution = float(arg) 
  elif opt == "--breaklines":
    breaklines_around_buildings = True
  elif opt == "--hires":
    higher_resolution_near_centre = float(arg)



# Now that settings are finalised, give a label for the run.
name = "urban_"        
name += "SG" if method == "DE_SG" else method
if breaklines_around_buildings: name += "b"
name += "_" + str(int(mesh_resolution))
if higher_resolution_near_centre: name += "_" + str(int(higher_resolution_near_centre))

print name

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

# define hi-res region
dx,dy = 0.5*building_size + building_spacing*izero + higher_resolution_padding_width
polygon_vertices = [(-dx,dy),(dx,dy),(dx,-dy),(-dx,-dy),(-dx,dy)]
higher_resolution = Path(polygon_vertices,closed=True).transformed(placement)

# conditionally declare implicitly-closed hi-res region for mesh algorithm
hires = higher_resolution_near_centre
assert hires is None or hires < mesh_resolution
interior_regions = [(higher_resolution.to_polygons()[0], hires)] if hires else []




# implicitly-closed polygon enclosing the floodplain plus the basin below the cliff
floodplain =  [(0,0),(0,cy),(cx,cy),(cx,0)]
cx += bucket_width
full_extent = [(0,0),(0,cy),(cx,cy),(cx,0)]



import anuga

meshname = name + '.msh' # if None then anuga won't generate_mesh

mesh = anuga.create_mesh_from_regions(full_extent,boundary_tags={'outer boundary':[0,1,2,3]},
  maximum_triangle_area=mesh_resolution, breaklines=breaklines, filename=meshname, interior_regions=interior_regions)

  
"""

Domain creation

Could start out with the floodplain dry (stage equal to elevation perhaps),
although it might converge faster if adding the steady-state water depth
(known analytically if we can neglect the buildings). 
Alternatively, could set a constant value (to see what happens when the flow reaches
a bucket that isn't quite empty).



"""

domain = anuga.create_domain_from_file(meshname)

domain.set_flow_algorithm(method)

domain.set_quantity('stage',cliff_fall) # start out dry? Oh wait..
domain.set_quantity('friction',manning_coefficient)
domain.set_quantity('elevation',elevation, location='centroids')
# Note assigning elev. @ centroids not vertices is critical to avoid
# smoothing, since deliberately aligned vertices onto breaklines
# where discontinuities are intended.

domain.set_name(name) # .sww output file
domain.set_datadir('.') # current working directory

domain.set_store_vertices_uniquely() # store multiple values (extrapolated separately from every adjoining cell) at each vertex

domain.set_boundary({'outer boundary': anuga.Reflective_boundary(domain)}) # hard walls

from anuga import Inlet_operator
left_edge = [(0.0,cy),(0,0)]
Inlet_operator(domain,left_edge,discharge) # set up constant water inflow

"""

Subgrid topography

Precompute tables relating water stage to volume, etc.

The tables are based on a prior constant approximation for the stage gradients.
In this scenario, the steady-state stage gradients should roughly match the bed-slope
(which in turn may crudely be approximated as zero). However, the average bed-slope  
would make a poor approximation for cells that straddle the boundaries of buildings.

"""

# express the constant as a vectorised function
friction = lambda x,y: manning_coefficient + 0*x

if (domain.flow_algorithm=='DE_SG'):
  domain.subgrid_data.make_subgrid_tables(elevation,friction,
    reference_gradient_type='zero',
    approx_grid_spacing=[15.,15.], # this defines the subgrid resolution
    max_reference_depth=25.0 # exceeds estimated maximum water depth, and "second max ref depth"..
    )
  domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()


"""

Run the simulation

"""


t1 = simulation_time 
#t1 = 30 # quicker debugging

for t in domain.evolve(yieldstep=50,finaltime=t1):
  print domain.timestepping_statistics()
  #print domain.timestep_fluxcalls, domain.max_flux_update_frequency, domain._order_, domain.default_order



"""
---------------------------------------------------------------------------------

Can compare plot with fig. 8 of Sanders et al. to check everything looks sensible

"""

# domain to plot
sx,sy = floodplain_size
sx += bucket_width


def draw_outlines(ax,*args): # draw buildings (and cliff) 
  shapes = Path.make_compound_path(town,Path(cliff),higher_resolution) # town + cliff
  ax.add_patch(PathPatch(shapes,fill=None,*args))


def mesh_elev(ax1): # show the setup prior to the domain existing.
  
  # randomly test elevation function
  x = np.random.rand(display_points)*sx
  y = np.random.rand(display_points)*sy  
  plt.scatter(x,y,c=elevation(x,y),alpha=0.5,cmap='Paired',linewidth=0, marker='.')
  plt.colorbar() # elevation is mapped to colour
  
  # check the mesh looks ok
  pts = mesh.getMeshVertices()
  paths = (Path([pts[i],pts[j],pts[k],pts[i]]) for i,j,k in mesh.getTriangulation())
  ax1.add_patch(PathPatch(Path.make_compound_path(*paths),fill=None,alpha=0.2))
  
  # check the buildings (and cliff) align with expectations
  draw_outlines(ax1) # linewidth=0.5,color='black',alpha=0.8,linestyle='dashed' ?
  


def triangulation(ax2): # show the raw output from the final time-step
  
  # obtain water depth
  water_depth = domain.get_quantity('height')
  x,y,values,triangles = water_depth.get_vertex_values(smooth=False)
  
  # interpolate separately within different triangles
  plt.tripcolor(x,y,triangles,values,shading='gouraud',alpha=0.8)
  plt.colorbar()
  
  # show momenta at centroids
  centx,centy = np.hsplit(domain.get_centroid_coordinates(),2)
  u = domain.get_quantity('xmomentum').get_values(location='centroids')
  v = domain.get_quantity('ymomentum').get_values(location='centroids')
  plt.quiver(centx,centy,u,v,alpha=0.6)
  


def interpolate(quantity,x=None,y=None): # smooth interpolation
    cx,cy = np.hsplit(domain.get_centroid_coordinates(),2) # centroid coordinates
    cz = quantity.get_values(location='centroids') # centroid quantity-values
    from scipy.interpolate import Rbf # Radial-basis-function interpolation method
    return Rbf(cx,cy,cz,function='cubic') # choose cubics for those basis functions

def grid(i): # regular grid over domain with given increment
    return np.meshgrid(np.arange(0,sx,i),np.arange(0,sy,i))


def interp(ax3): # smooth interpolation of results
  

  
  water_depth = domain.get_quantity('stage') - domain.get_quantity('elevation')
  
  # produce image of smoothly-interpolated water depth
  X,Y = grid(5)
  Z = interpolate(water_depth)
  plt.pcolormesh(X,Y,Z(X,Y),shading='gouraud',cmap='gist_rainbow')
  # tries to match colors from paper.
  
  plt.colorbar()
  
  # present the water velocities
  u = interpolate(domain.get_quantity('xmomentum')/water_depth)
  v = interpolate(domain.get_quantity('ymomentum')/water_depth)
  X,Y = grid(30)
  plt.quiver(X,Y,u(X,Y),v(X,Y),alpha=0.6)
  # might prefer a streamline plot?
  
  # superimpose buildings and cliff again
  draw_outlines(ax3) # linewidth=0.5,color='black',alpha=0.8 ?
  
  
def closeup():  # Sanders et. al. fig.9
  
  # construct grids
  i = 2
  j = 20
  x,y = np.meshgrid(np.arange(550,950,i),np.arange(75,425,i))
  X,Y = np.meshgrid(np.arange(550,950,j),np.arange(75,425,j))
  
  # gather data
  water_depth = domain.get_quantity('height')
  u = interpolate(domain.get_quantity('xmomentum')/water_depth)
  v = interpolate(domain.get_quantity('ymomentum')/water_depth)
  Z = interpolate(water_depth)
  
  # configure plot
  fig = plt.figure()
  ax = fig.add_subplot(111,aspect='equal') # set aspect ratio
  ax.set_xlim(550,950)
  ax.set_ylim(75,425)
  
  # draw
  plt.pcolormesh(x,y,Z(x,y),cmap='gist_rainbow_r',vmax=1.15,vmin=0.45)
  plt.colorbar()
  plt.quiver(X,Y,u(X,Y),v(X,Y),alpha=0.6)
  draw_outlines(ax) # buildings
  




if display_figure:
  import matplotlib.pyplot as plt
  from matplotlib.patches import PathPatch
  fig = plt.figure()
  
  def subfig(total,index):
    padding = 30
    xl,yl = [-padding,cx+padding],[-padding,cy+padding] # plot extent
    ax = fig.add_subplot(total,1,index,aspect='equal') # set aspect ratio
    #plt.title(label)
    plt.xlim(xl)
    plt.ylim(yl)
    return ax
  
  def do_plots(*args):
    n = len(args)
    for i,func in enumerate(args,start=1):
      func(subfig(n,i))
    

  #do_plots(mesh_elev,triangulation,interp)
  do_plots(mesh_elev,triangulation)
  
  #closeup()
  
  plt.show() # pause for a look at the outputs
