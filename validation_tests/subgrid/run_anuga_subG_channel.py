"""

Water flowing down a channel with a floodplain.

This example checks to what extent we can coarsely resolve a channel, and still
get the correct steady-uniform flow answer

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
floodplain_width = 20.0 # Model domain width
floodplain_slope = 1./300.
chan_initial_depth = 0.60 # Initial depth of water in the channel
chan_bankfull_depth = 1.0 # Bankfull depth of the channel
chan_width = 10. # Bankfull width of the channel
bankwidth = 2. # Width of the bank regions -- note that these protrude into the channel
man_n = 0.03 # Manning's n
l0 = 5.0 # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)
Qin = 4.5
flow_in_yval = 5.0
flow_algorithm = 'DE_SG'
reference_gradient_type = 'zero' #'bed-slope'


assert chan_width < floodplain_width, \
        ' ERROR: Channel width is greater than floodplain width'

assert bankwidth < chan_width/2., \
        'ERROR: The bank width must be less than half the channel width'


# Function for topography
def topography(x, y):
    # Longitudinally sloping floodplain with channel in centre
    elev1 = -y*floodplain_slope - chan_bankfull_depth*\
            (x > (floodplain_width/2. - chan_width/2.))*\
            (x < (floodplain_width/2. + chan_width/2.)) 
    # Add banks
    if(bankwidth > 0.0):
        leftbnk = floodplain_width/2. - chan_width/2.
        rightbnk = floodplain_width/2. + chan_width/2.
        # Left bank
        elev2 = elev1 + (chan_bankfull_depth \
                - chan_bankfull_depth/bankwidth*(x - leftbnk))*\
                (x > leftbnk)*(x < leftbnk + bankwidth)
        # Right bank
        elev2 = elev2 + (chan_bankfull_depth \
                + chan_bankfull_depth/bankwidth*(x - rightbnk))*\
                (x > rightbnk - bankwidth)*(x < rightbnk)
    
    if(bankwidth == 0.0):
        elev2 = elev1

    # Add a wall around the domain
    #sides = ( ( floodplain_width - x< l0) + (x < l0) + (floodplain_length -y < l0) + (y<l0)).nonzero()[0]
    #elev2[sides]=3.0

    return elev2

#Function for stage
def stageFun(x,y):
    return -y*floodplain_slope - chan_bankfull_depth + chan_initial_depth 

def frictionFun(x,y):
    return man_n + x*0.


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

breakLines={'n1':[[floodplain_width/2.- chan_width/2., floodplain_length],
                 [floodplain_width/2. - chan_width/2.,0.]],
            'n2':[[floodplain_width/2.+ chan_width/2., floodplain_length],
                  [floodplain_width/2. + chan_width/2.,0.]]}
#            # These ones are inside the channel 
#            'n3':[[floodplain_width/2.- chan_width/2.+bankwidth, floodplain_length],
#                 [floodplain_width/2. - chan_width/2.+bankwidth,0.]],
#            'n4':[[floodplain_width/2.+ chan_width/2.-bankwidth, floodplain_length],
#                  [floodplain_width/2. + chan_width/2-bankwidth,0.]]
#            }

regionPtAreas = [ [ 0.01, 0.01, 0.5*l0*l0],
                  [ floodplain_width/2., 0.01, 0.5*l0*l0], 
                  [ floodplain_length - 0.01, 0.01, 0.5*l0*l0] ]

#regionPtAreas=[[0.01, 0.01, 0.5*l0*l0],
#               [floodplain_width/2.-chan_width/2.+0.01, 0.01, 0.5*l0*l0*0.25],
#               [floodplain_width/2.-chan_width/2.+bankwidth+0.01, 0.01, 0.5*l0*l0],
#               [floodplain_width/2.+chan_width/2.-bankwidth+0.01, 0.01, 0.5*l0*l0*0.25],
#               [floodplain_width/2.+chan_width/2.+0.01, 0.01, 0.5*l0*l0]]



# Main script here

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
    #domain.set_flow_algorithm('DE0')
    #domain.set_flow_algorithm('DE1')

    #domain.set_timestepping_method('euler')

    domain.set_name('channel') # Output name
    domain.set_store_vertices_uniquely(True)



    if domain.flow_algorithm=='DE_SG' and False:
        from anuga.utilities import subGridUtil as sgu

        mean_topography=sgu.make_spatially_averaged_Fun(
            topography, domain, approx_grid_spacing = [0.5, 0.5])

        domain.set_quantity('elevation', mean_topography, location='centroids')

    else:
        domain.set_quantity('elevation', topography, location='centroids')

    domain.set_quantity('stage', stageFun)
    domain.set_quantity('xmomentum', 0.)
    domain.set_quantity('ymomentum', 0.)
    domain.set_quantity('friction', frictionFun)


    if(domain.flow_algorithm == 'DE_SG'):
        domain.set_quantity('alphax', 0.)
        domain.set_quantity('alphay', 0.)

        # Make subgrid lookup functions
        domain.subgrid_data.make_subgrid_tables(
            topography, 
            frictionFun, 
            approx_grid_spacing=[0.5,0.5], 
            max_reference_depth=10.0,
            constant_subgrid_values=False,
            reference_gradient_type=reference_gradient_type)        
            #reference_gradient_type='bed-slope')        

        # Set initial values for volume quantities
        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()

    ######################################################################
    ## DISCHARGE 
    line1 = [ [floodplain_width/2. - chan_width/2., flow_in_yval],\
              [floodplain_width/2. + chan_width/2., flow_in_yval] \
              ]

    Inlet_operator(domain, line1, Qin)

    print 'Discharge in = ', Qin 


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

    for t in domain.evolve(yieldstep=5.0, finaltime=100*60.0):
        print domain.timestepping_statistics()

        domain.report_water_volume_statistics()

        if domain.flow_algorithm == 'DE_SG':
            vv = domain.quantities['u_vol'].centroid_values**2+domain.quantities['v_vol'].centroid_values**2
            velsq = (vv/(domain.quantities['vol'].centroid_values**2+1.0e-10))
            print 'vel_max: ', (velsq.max())**0.5
            vloc = velsq.argmax()
            print 'Depth@vel_max: ', domain.quantities['vol'].centroid_values[vloc]/domain.subgrid_data.subgrid_wet_area[vloc]

