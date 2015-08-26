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

# scenario constants: (refer to pg.576 pf citation, do not change)

bed_slope = 2e-4 # constant slope in flow direction
outer_radius = 320 # metres
incision_width = 213 # m, maximum width of river before flooding onto plain
manning_coefficient = 0.026 # m^(-1/3) s
total_length = 960 # m
cross_section_amplitude = 4 # m

# extra scenario constants: 

bed_drop_over_semicircle = 0 # m. 
pond_incision_depth = -30 # m
pond_width = 10#00 # m
pond_length = 10#00 # m

# initial values

upper_stage = -1.
lower_stage = -2.

# solver parameters

mesh_resolution = 0.5*(50**2) # m^2, maximum triangle area



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

The paper is not explicit about whether the slope continues through the bend itself.
It does utilise a system of coordinates where the slope has the same effective
length as either arm (i.e. a typical radius of ~203m). 

It is not explicit (in the original publication) whether the slope continues through the 
bend itself. It does however indicate (fig.13) that the steady-state stage (a proxy for
elevation) varies by about 0.30m over the course of the U-bend (as though the semi-circle 
had incorporated a couple hundred metres worth of slope). 

The paper does adopt a system of coordinate where the semi-circle has
the same "dimensionless" length as either arm.


"""

inner_radius = outer_radius - incision_width
arm_length = total_length - outer_radius

psuedo = abs(bed_drop_over_semicircle)/bed_slope # calculate psuedo-distance around semi-circle

from math import pi
import numpy as np

# for the cross-section: in the upper half-space, function primarily of radial coordinate.
#                        in the lower half-space, function primarily of absolute lateral displacement.

# between the inside and outside radius spans a quarter period of a sinusoidal curve
angular_freq = pi/(2*incision_width)
sinusoid = lambda r: cross_section_amplitude * np.sin(angular_freq * (r - inner_radius))
def cross_section(r,z):
  # note, cannot chain numpy inequalities a<r<b
  return np.where((inner_radius<r)&(r<outer_radius), sinusoid(r) - bed_slope * z, 0.)
  
def elevation_for_bend(x,y):
  radius = np.sqrt(x**2+y**2) 
  theta = np.arctan2(y,x) # ranges from 0 (at right) to pi (at left) for y +ve.
  return cross_section(radius, theta * psuedo/pi)

boundary_pond = lambda offset: np.where(offset < inner_radius, 0., pond_incision_depth)
def elevation_for_straights(x,y):
  abs_x = np.fabs(x) # |x| is symmetric across y-axis
  distance = y * np.sign(x) + arm_length + np.where(x < 0, psuedo, 0.) # from upstream pond
  return np.where(y < -arm_length, boundary_pond(abs_x), cross_section(abs_x,distance))
  
elev = lambda x,y: np.where(y>0, elevation_for_bend(x,y), elevation_for_straights(x,y))
 


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



initial_stage = lambda x,y: np.where(x > 0, upper_stage, lower_stage)

friction = lambda x,y: manning_coefficient



import anuga

pond_radius = pond_width + inner_radius
south = -(arm_length + pond_length)
rectangle = lambda x1,y1,x2,y2: [(x1,y1),(x1,y2),(x2,y2),(x2,y1)]
u_bend_area = rectangle(-outer_radius, -arm_length, outer_radius, outer_radius)
ponds_area = rectangle(-pond_radius, south, pond_radius, -arm_length)
full_extent = ponds_area[:2] + u_bend_area + ponds_area[2:] # T-shaped domain

# align the mesh to the (inside) pond edges using break-lines, or cut-out part of the region entirely
full_extent += [(inner_radius,south),(inner_radius,-arm_length) , (-inner_radius,-arm_length),(-inner_radius,south)]



meshfile = 'ubend.msh'
mesh = anuga.create_mesh_from_regions(full_extent,
        boundary_tags={'outer boundary':range(len(full_extent))},
        filename=meshfile,
	interior_regions=[(u_bend_area, mesh_resolution)])
domain = anuga.create_domain_from_file(meshfile)


domain.set_quantity('elevation',elev,location='centroids') # assign from centroid positions instead of smoothing
domain.set_quantity('friction',manning_coefficient)
domain.set_quantity('stage',initial_stage)


import re
print filter(re.compile('.*ref.*').match,dir(mesh))
print filter(re.compile('.*ref.*').match,dir(domain))
#print dir(domain)

print domain.geo_reference
print mesh.geo_reference
print mesh.geo_reference==domain.geo_reference
print dir(domain.geo_reference)

quit()


import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')

"""
from matplotlib.path import Path
from matplotlib.patches import PathPatch
pts = mesh.getMeshVertices()
paths = (Path([pts[i],pts[j],pts[k],pts[i]]) for i,j,k in mesh.getTriangulation())
ax.add_patch(PathPatch(Path.make_compound_path(*paths),fill=None))
x,y = zip(*pts) #x,y = zip(*full_extent) inapplicable due to internal transformation
"""

x,y,values,triangles = domain.get_quantity('elevation').get_vertex_values(smooth=False)
plt.tripcolor(x,y,triangles,elev(x,y),shading='gouraud',alpha=0.2)
print min(values),max(values)


limits = lambda x: (min(x),max(x))
ax.set_xlim(*limits(x))
ax.set_ylim(*limits(y))
plt.show()


