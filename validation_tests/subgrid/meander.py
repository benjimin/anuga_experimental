"""

This scenario represents a meander in a steady flowing river.

The riverbed is a U-bend shaped like the Panathenaic Stadium,
and the banks are only partially filled by the river.



This test is based from (and aims to replicate):

G.S. Stelling, Quadtree flood simulations with sub-grid digital
elevation models, Proc. Inst. Civil Eng., Water Management 165 567-580 (2012).



A motivating feature of this scenario is that the riverbed elevation will vary 
continuously, but not only in a trivial plane, thus serving as a test for 
within-cell topographical-variation. At least one of the boundaries of the wetted 
region does not correspond to infinitely-steep walls, but is free to shift spatially 
according to the flow.

"""

# scenario constants: (refer to pg.576 of citation, do not change)

bed_slope = 2e-4 # constant downhill slope in flow direction
outer_radius = 320 # metres
incision_width = 213 # m, maximum width of river before flooding onto plain
manning_coefficient = 0.026 # m^(-1/3) s
total_length = 960 # m
cross_section_amplitude = -4 # m
total_drop = 0.30 # m, taken from fig.13 pg.577

# initial values

max_depth = 0.9 # m, at inflow/outflow (e.g., 3.9, 1.9 or 0.9; table 1, pg.576)

# extra scenario constants: 

#pond_incision_depth = -8 # m
pond_width = 1000 # m
pond_length = 1000 # m
delta = True # boolean, whether to try softening the pond transition; overrides depth.

# time-stepping options

time = 5000 # s, total simulation time to produce (an approximation of) convergence
time_between_output = 100 # s, temporal resolution of output

# solver parameters

mesh_resolution = 0.5*(15**2) # m^2, maximum triangle area;
                              # note Stelling uses squares of 20*20, 40*40 and 80*80 m^2,
			      # albeit with 1m^2 subgrid.
method = 'DE1'


"""
The topography is akin to a negative of the Panathenaic Stadium.
Two straight and parallel channels, connected by a semi-circular section
to form the shape of the letter U.

The ground elevation is deepest at the outside (with an abrupt discontinuity) 
and shallows toward the inside (with a flat region in the middle, level with 
the surrounding plain).

The cross-section of the riverbed is sinusoidal (i.e. a quarter-cycle, so the 
transverse or radial gradient vanishes where the bed is deepest, and is steepest
where it is shallowest). 

There is a small constant linear bed slope, in the direction of the river. This means
that the incision slightly deepens as the river approaches the bend, and deepens
further as the river moves away from the bend. (This accumulates as a small
discontinuity at the inside edge of the incision.) 

The original publication is not explicit about whether the slope continues through 
the bend itself. It does however indicate (fig.13) the approximate amount that the
steady-state stage differs between ends of the U-bend. This is a good proxy for
the accumulated elevation offset, because the boundary conditions are described
as imposing the same max water depth at each end of the U-bend, and the figure 
also verifies that (to first order) the water has approximately the same (flat) 
cross-section at either end. (The implication turns out to be that the semicircle 
takes an effective length of 220m.)
"""

inner_radius = outer_radius - incision_width
arm_length = total_length - outer_radius


# Calculate how the bed slope continues around the bend
bed_drop_over_semicircle = total_drop - 2 * arm_length * bed_slope
assert bed_slope>=0 and total_drop>=0 and bed_drop_over_semicircle>=0
pseudo = bed_drop_over_semicircle/bed_slope # pseudo-length of semi-circle

from math import pi
import numpy as np




# between the inside and outside radius spans a quarter period of a sinusoidal curve
angular_freq = pi/(2*incision_width)
def sinusoid(r): 
  return cross_section_amplitude * np.sin(angular_freq * (r - inner_radius))
def cross_section(r,z):
  # note, cannot chain numpy inequalities a<r<b
  return np.where((inner_radius<r)&(r<outer_radius), sinusoid(r) - bed_slope * z, 0.)


# in the upper half-space, cross-section is a function primarily of radial coordinate.
def elevation_for_bend(x,y):
  radius = np.sqrt(x**2+y**2) 
  theta = np.arctan2(y,x) # ranges from 0 (at right) to pi (at left) for y +ve.
  return cross_section(radius, theta * pseudo/pi)


# in the lower half-space, cross-section is a function primarily of absolute lateral displacement.
def elevation_for_straights(x,y):
  abs_x = np.fabs(x) # |x| is symmetric across y-axis
  distance = y * np.sign(x) + arm_length + np.where(x < 0, pseudo, 0.) # dist. from upstream pond
  return np.where(y < -arm_length, boundary_pond(abs_x,distance), cross_section(abs_x,distance))


# incorporate the lakes to buffer the inflow/outflow water levels.
def boundary_pond(offset,*args): 
  return np.where(offset < inner_radius, 0., pond_incision_depth) 
if delta:
  def boundary_pond(offset, distance):
    drop = total_drop * (distance > 0) # binary
    return np.where(offset < outer_radius, sinusoid(offset) - drop, cross_section_amplitude)
  


def elevation(x,y): # overall  
  return np.where(y>0, elevation_for_bend(x,y), elevation_for_straights(x,y))
 


"""

We are interested in the steady-state solution.

Unfortunately, as we are approximating the inflow/outflow boundaries using ponds of finite
size, we are limited in how long the simulation can run (before it evolves too far from the
target state). 

We can only minimise this evolution rate by using ponds with large area (and with large
depths to minimise additional effects related to the physics of the pond mouths). 

We will need to monitor the evolving water level in each pond. This will also be used to
compute the flow rates.

As a simple initial condition, we can set the left and right sides of the scenario to have
different water levels. This is effectively a dam-break (but is the best that a dam-break
can be, in terms of approximating the steady-state solution). 

"""


# initial values

upper_stage = cross_section_amplitude + max_depth
lower_stage = upper_stage - total_drop
#print cross_section_amplitude , lower_stage , upper_stage , 0.
assert cross_section_amplitude < lower_stage < upper_stage < 0.

def initial_stage(x,y): return np.where(x > 0, upper_stage, lower_stage)





import anuga

pond_radius = pond_width + inner_radius
south = -(arm_length + pond_length)
def rectangle(x1,y1,x2,y2): return [(x1,y1),(x1,y2),(x2,y2),(x2,y1)]
u_bend_area = rectangle(-outer_radius, -arm_length, outer_radius, outer_radius)
ponds_area = rectangle(-pond_radius, south, pond_radius, -arm_length)
full_extent = ponds_area[:2] + u_bend_area + ponds_area[2:] # T-shaped domain

# align the mesh to the (inside) pond edges using break-lines, or cut-out part of the region entirely
full_extent += [(inner_radius,south),(inner_radius,-arm_length) , (-inner_radius,-arm_length),(-inner_radius,south)]
#full_extent = map(list,full_extent) if tuples were a problem..



meshfile = 'ubend.msh'
mesh = anuga.create_mesh_from_regions(full_extent,
        boundary_tags={'outer boundary':range(len(full_extent))},
        filename=meshfile,
	interior_regions=[(u_bend_area, mesh_resolution)])
domain = anuga.create_domain_from_file(meshfile)



# ------------- assign quantities ---

# shift to using internal coordinates
# (ironically, losing numerical-precision due to a feature designed to mitigate such)
def translate(f): 
  xll,yll = map(min,zip(*full_extent))
  return lambda x,y: f(x+xll,y+yll)
  
elevation = translate(elevation)
initial_stage = translate(initial_stage)

domain.set_quantity('elevation',elevation,location='centroids') # assign from centroid positions instead of smoothing
domain.set_quantity('friction',manning_coefficient)
domain.set_quantity('stage',initial_stage) # smoothing here will assist convergence



# -------- set up solver -------



domain.set_flow_algorithm(method)

domain.set_name('meander.sww')
domain.set_datadir('.')
domain.set_store_vertices_uniquely()

domain.set_boundary({'outer boundary': anuga.Reflective_boundary(domain)}) # hard walls


if method=='DE_SG': # prepare sub-grid tables
  def friction(x,y): return manning_coefficient + 0*x # vectorised
  kwargs = {'reference_gradient_type':'zero', 'max_reference_depth':10}
  domain.subgrid_data.make_subgrid_tables(elevation,friction,**kwargs)
  domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()







# run simulation
for t in domain.evolve(yieldstep=time_between_output,finaltime=time):
  print domain.timestepping_statistics()


"""
Current idea is not to set any inlet.

Alternatively, could have a inlet and outlet. 

But trying to balance flow rates to maintain surface levels... might be simpler just
to make the operator for constant levels work with subgrid.

So: question: regards the code for anuga.Inlet_operator.Inlet_operator, 
is anything special about getting it to work for subgrid?

Or.. could use flow rates as my independent variable?
"""






"""

The paper plots results in a system of coordinates where the slope has the same effective
length as either arm (i.e. a typical radius of ~203m). 
"""


raise SystemExit # skip plots




import matplotlib.pyplot as plt



fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')

"""from matplotlib.path import Path
from matplotlib.patches import PathPatch
pts = mesh.getMeshVertices()
paths = (Path([pts[i],pts[j],pts[k],pts[i]]) for i,j,k in mesh.getTriangulation())
ax.add_patch(PathPatch(Path.make_compound_path(*paths),fill=None))
x,y = zip(*pts) #x,y = zip(*full_extent) inapplicable due to internal transformation"""

x,y,values,triangles = domain.get_quantity('stage').get_vertex_values(smooth=False)
plt.tripcolor(x,y,triangles,values,shading='gouraud',alpha=0.2,vmax=3)
print 'Output value range',min(values),max(values)
plt.colorbar()

limits = lambda x: (min(x),max(x))
ax.set_xlim(*limits(x))
ax.set_ylim(*limits(y))


x,y = np.hsplit(domain.get_centroid_coordinates(),2)
values = domain.get_quantity('stage').get_values(location='centroids')
h = domain.get_quantity('height').get_values(location='centroids')
def distance(x,y): # so-called "dimensionless distance"
  return np.where(y>0, np.arctan2(y,x)/pi, np.where(x>0,y/arm_length,1-y/arm_length))+1
fig = plt.figure()
X,Y = translate(distance)(x,y)[h!=0], values[h!=0]
plt.scatter(X,Y,alpha=0.2)


plt.show()


