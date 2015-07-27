"""

Very shallow water flowing down a coarsely resolved slope

Many algorithms have trouble with this when the elevation range of each cell
is much larger than the centroid water depth. But it is an important case to
get basically correct.

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Import standard shallow water domain and standard boundaries.
import anuga
import numpy
from anuga.structures.inlet_operator import Inlet_operator
from anuga import create_domain_from_regions
#------------------------------------------------------------------------------
# Useful parameters for controlling this case
#------------------------------------------------------------------------------

stage_rel_to_plain = -0.05 # Target steady state analytic water depth
drain_depth = -1
filename='drainage'

floodplain_length = 1000.0 # Model domain length
floodplain_width = 100.0 # Model domain width
floodplain_slope = 1./10.
chan_initial_depth = 0.0 # Initial depth of water in the channel
chan_bankfull_depth = 0.0 # Bankfull depth of the channel
chan_width = 5. # Bankfull width of the channel
bankwidth = 2. # Width of the bank regions -- note that these protrude into the channel
man_n=0.03 # Manning's n
l0 = 15 # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)
l0 = 2

flow_in_yval = 5.0 # y-value of line along which the input discharge is passed
#Qin = 0.5 # Input discharge

flow_algorithm = 'DE1'
reference_gradient_type = 'bed-slope'

drain_count = 3
drains_proportion_of_width = 0.25
drains_proportion_of_length = 0.8


"""
We know Q = sqrt(abs(slope))/n_manning * integral of depth^(5/3) across transverse axis
And since the integrand is piecewise constant...
"""
#assert stage_rel_to_plain > drain_depth
drain_portion = drains_proportion_of_width * max(0.0,stage_rel_to_plain-drain_depth)**(5.0/3.0)
plain_portion = (1.0-drains_proportion_of_width) * max(0.0,stage_rel_to_plain)**(5.0/3.0)
Qin = abs(floodplain_slope)**(0.5)/man_n * floodplain_width*(drain_portion + plain_portion)
assert Qin > 0

assert int(drain_count)>=1
assert 0<drains_proportion_of_width<1
assert 0<drains_proportion_of_length<1

assert chan_width < floodplain_width, \
        ' ERROR: Channel width is greater than floodplain width'

assert bankwidth < chan_width/2., \
	'ERROR: The bank width must be less than half the channel width'

def topography(x,y):
    elev=-y*floodplain_slope
    # modification to cut some drains into the bed
    drain_width = (floodplain_width/drain_count)
    padding = drain_width * 0.5*(1-drains_proportion_of_width)
    x = x % drain_width
    """
    if float(y)/floodplain_length > (1-drains_proportion_of_length):
        if padding<x<(drain_width-padding):
            elev += drain_depth
    """
    # try to make it vectorisable
    is_drain = ((y/floodplain_length > (1-drains_proportion_of_length))
                  * (padding<x) * (x<(drain_width-padding)) )
    elev += is_drain*drain_depth
    return elev

# Very low friction [frictionless would be ideal for testing, but can't be done]
def frictionFun(x,y):
    return man_n + x*0.

def stageFun(x,y):
    out = topography(x,y)-100. # Dry
    return out

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

#breakLines={}

breakLines = {'n1': [ [0., 11.], [floodplain_width, 11.] ] }

#breakLines={'n1':[[floodplain_width/2.- chan_width/2., floodplain_length],
#                 [floodplain_width/2. - chan_width/2.,0.]],
#            'n2':[[floodplain_width/2.+ chan_width/2., floodplain_length],
#                  [floodplain_width/2. + chan_width/2.,0.]]}
#            # These ones are inside the channel 
#            'n3':[[floodplain_width/2.- chan_width/2.+bankwidth, floodplain_length],
#                 [floodplain_width/2. - chan_width/2.+bankwidth,0.]],
#            'n4':[[floodplain_width/2.+ chan_width/2.-bankwidth, floodplain_length],
#                  [floodplain_width/2. + chan_width/2-bankwidth,0.]]
#            }

regionPtAreas = [ [ 0.01, 0.01, 0.01*l0*l0],
                  [ floodplain_width-0.01, floodplain_length-0.01, 0.5*l0*l0]]

#regionPtAreas=[[0.01, 0.01, 0.5*l0*l0],
#               [floodplain_width/2.-chan_width/2.+0.01, 0.01, 0.5*l0*l0*0.25],
#               [floodplain_width/2.-chan_width/2.+bankwidth+0.01, 0.01, 0.5*l0*l0],
#               [floodplain_width/2.+chan_width/2.-bankwidth+0.01, 0.01, 0.5*l0*l0*0.25],
#               [floodplain_width/2.+chan_width/2.+0.01, 0.01, 0.5*l0*l0]]

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

    domain.set_name(filename) # Output name
    domain.set_store_vertices_uniquely(True)

    domain.quantities_to_be_stored['elevation'] = 2

    #from anuga.utilities import subGridUtil as sgu
    #mean_topography=sgu.make_spatially_averaged_Fun(
    #    topography, domain, approx_grid_spacing = [0.5, 0.5])
    #domain.set_quantity('elevation', mean_topography, location='centroids')
    domain.set_quantity('elevation', topography) #, location='centroids')
    domain.set_quantity('stage', stageFun)
    domain.set_quantity('xmomentum', 0.)
    domain.set_quantity('ymomentum', 0.)
    domain.set_quantity('friction', frictionFun)
    #domain.set_quantity('height',0.1)

    if(domain.flow_algorithm=='DE_SG'):
        domain.set_quantity('alphax',0.)
        domain.set_quantity('alphay',0.)

        # Make subgrid lookup functions
        domain.subgrid_data.make_subgrid_tables(
            topography, frictionFun, 
            approx_grid_spacing=[5.0,5.0], 
            max_reference_depth=10., 
            constant_subgrid_values=False,
            reference_gradient_type='bed-slope') 
            #reference_gradient_type='zero')        

        # Set initial values for volume quantities
        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
    else:
        #domain.beta_w = 0.0
        #domain.beta_uh = 0.0
        #domain.beta_vh = 0.0
        #domain.beta_w_dry = 0.0
        #domain.beta_uh_dry = 0.0
        #domain.beta_vh_dry = 0.0
        pass

    ######################################################################
    ## DISCHARGE
    if True:
        line1 = [ [0., flow_in_yval],\
                  [floodplain_width, flow_in_yval] \
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

    #domain.CFL=0.1

    for t in domain.evolve(yieldstep=5.0,finaltime=100*60.0):
        print domain.time, domain.edge_timestep.min()
        domain.report_water_volume_statistics()
        #domain.report_cells_with_small_local_timestep()
        if domain.flow_algorithm == "DE_SG":
            vv = domain.quantities['u_vol'].centroid_values**2+domain.quantities['v_vol'].centroid_values**2
            velsq = (vv/(domain.quantities['vol'].centroid_values**2+1.0e-10))
            print 'vel_max: ', (velsq.max())**0.5
            vloc = velsq.argmax()
            print 'Depth@vel_max: ', domain.quantities['vol'].centroid_values[vloc]/domain.subgrid_data.subgrid_wet_area[vloc]
            #print domain.timestepping_statistics()
            #print 'vol: ', domain.quantities['vol'].centroid_values.sum()
            #print 'BF: ',  domain.get_boundary_flux_integral()
            #print 'Vol - BF: ',  domain.quantities['vol'].centroid_values.sum() - domain.get_boundary_flux_integral()

    # Analytical solution
    W = floodplain_width
    Sf = floodplain_slope
    n = man_n
    Q = Qin

    u = ( Sf/n**2 *(Q/W)**(4./3.))**(3./10.)
    h = Q/(u*W)

    print 'Analytical u: ', u, ' h: ', h
