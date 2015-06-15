"""

Stationary water over topography (2 basins connected by a narrow channel)

This is used to check well balancing for the subgrid algorithm

Velocities should be 'numerical zero' (e.g. < 1.0e-10) throughout the
simulation

It seems a good idea to run this with stages both below and above the maximum
bed elevation

Also, this case can simply be modified to have different stages in each basin,
which provides a good dynamic test of the subgrid algorithm

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import anuga
import numpy
from anuga.structures.inlet_operator import Inlet_operator
from anuga import create_domain_from_regions

#------------------------------------------------------------------------------
# Useful parameters for controlling this case
#------------------------------------------------------------------------------

floodplain_length = 2000.0 # Model domain length
floodplain_width = 100.0 # Model domain width
floodplain_slope = 1./300.
chan_width = 5. # Bankfull width of the channel
man_n = 0.03 # Manning's n
l0 = 20. # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)

initial_stage = 2.50

flow_algorithm = 'DE_SG' # DE0 # DE1
reference_gradient_type = 'bed-slope' # 'zero'

assert chan_width < floodplain_width, \
        ' ERROR: Channel width is greater than floodplain width'


def topography(x,y):
    # 2 Basins with a narrow gap connecting them
    elev = x*0.0
    elev = elev + 3.0*( (y > (floodplain_length/2.0 - 100.))*\
                        (y < (floodplain_length/2.0 + 100.))*\
                        ((x < (floodplain_width/2.0 - chan_width/2.0)) +\
                         (x > (floodplain_width/2.0 + chan_width/2.0))) )
    
    return elev


# Set friction 
def frictionFun(x,y):
    return man_n + x*0.0


def stageFun(x,y):
    return initial_stage + x*0.0


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

# Define boundary polygon -- in this case clockwise around the proposed boundary
boundary_polygon = [ [0.,0.], 
                     [0., floodplain_length], 
                     [floodplain_width/2. - chan_width/2., floodplain_length], 
                     [floodplain_width/2. + chan_width/2., floodplain_length], 
                     [floodplain_width, floodplain_length], 
                     [floodplain_width, 0.], 
                     [floodplain_width/2. + chan_width/2., 0.], 
                     [floodplain_width/2. - chan_width/2., 0.] 
                     ]

## Define channel polygon, where we can optionally refine the resolution. 
## Note that we set this a distance l0 inside the boundary polygon, so the polygons
## do not overlap. 

breakLines={}

regionPtAreas = [ [ 0.01, 0.01, 0.5*l0*l0]] # x, y, triangle area


if __name__ == '__main__':

    # Define domain with appropriate boundary conditions
    anuga.create_mesh_from_regions(boundary_polygon, 
                             boundary_tags={'left': [0],
                                            'top1': [1],
                                            'chan_out': [2],
                                            'top2': [3],
                                            'right': [4],
                                            'bottom1': [5],
                                            'chan_in': [6],
                                            'bottom2': [7] },
                             maximum_triangle_area = 0.5*l0*l0,
                             minimum_triangle_angle = 28.0,
                             filename = 'channel_floodplain1.msh',
                             breaklines=breakLines.values(),
                             regionPtArea=regionPtAreas,
                             use_cache=False,
                             verbose=False)

    domain=anuga.create_domain_from_file('channel_floodplain1.msh')

    domain.set_flow_algorithm(flow_algorithm)

    domain.set_name('balance_check') # Output name
    domain.set_store_vertices_uniquely(True)

    domain.set_quantity('elevation', topography, location='centroids')
    domain.set_quantity('stage', stageFun)
    domain.set_quantity('xmomentum', 0.0)
    domain.set_quantity('ymomentum', 0.0)
    domain.set_quantity('friction', frictionFun)

    #
    # Set the alphax, alphay quantities
    #
    if(domain.flow_algorithm == 'DE_SG'):
        domain.set_quantity('alphax',0.)
        domain.set_quantity('alphay',0.)

        # Make subgrid lookup functions
        domain.subgrid_data.make_subgrid_tables(
            topography,
            frictionFun,
            approx_grid_spacing = [1.0, 1.0],
            max_reference_depth = 100.0,
            constant_subgrid_values = False,
            reference_gradient_type = reference_gradient_type)

        # Set initial values for volume quantities
        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()


    ## Boundary conditions
    Br = anuga.Reflective_boundary(domain) # Solid reflective wall
    domain.set_boundary({'left': Br, 
                         'right': Br, 
                         'top1': Br, 
                         'top2': Br, 
                         'bottom1': Br, 
                         'bottom2': Br, 
                         'chan_out': Br, 
                         'chan_in': Br})

    domain.CFL = 0.9

    for t in domain.evolve(yieldstep=1.0, finaltime=10*60.0):

        print domain.time, domain.timestep

        vv = domain.quantities['u_vol'].centroid_values**2+domain.quantities['v_vol'].centroid_values**2

        # A velocity^2 like quantity -- avoid division by zero
        velsq = (vv/(domain.quantities['vol'].centroid_values**2 + 1.0e-10))

        print 'vel_max: ', (velsq.max())**0.5

        vloc = velsq.argmax()
        print 'Depth@vel_max: ', domain.quantities['vol'].centroid_values[vloc]/domain.subgrid_data.subgrid_wet_area[vloc]

