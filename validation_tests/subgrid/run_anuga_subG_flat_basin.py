"""

Simple water flow example using ANUGA
Water flowing down a channel with a floodplain

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

floodplain_length = 2000.0 # Model domain length
floodplain_width = 100.0 # Model domain width
floodplain_slope = 1./300.
chan_initial_depth = 0.60 # Initial depth of water in the channel
chan_bankfull_depth = 1.0 # Bankfull depth of the channel
chan_width = 5. # Bankfull width of the channel
bankwidth = 2. # Width of the bank regions -- note that these protrude into the channel
man_n=0.03 # Manning's n
l0 = 30. # Length scale associated with triangle side length in channel (min_triangle area = 0.5*l0^2)

assert chan_width < floodplain_width, \
        ' ERROR: Channel width is greater than floodplain width'

assert bankwidth < chan_width/2., \
        'ERROR: The bank width must be less than half the channel width'

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
#breakLines={'n1':[[floodplain_width/2.- chan_width/2., floodplain_length],
#                 [floodplain_width/2. - chan_width/2.,0.]],
#            'n2':[[floodplain_width/2.+ chan_width/2., floodplain_length],
#                  [floodplain_width/2. + chan_width/2.,0.]]
#            #,
#            }
            ## These ones are inside the channel 
            #'n3':[[floodplain_width/2.- chan_width/2.+bankwidth, floodplain_length],
            #     [floodplain_width/2. - chan_width/2.+bankwidth,0.]],
            #'n4':[[floodplain_width/2.+ chan_width/2.-bankwidth, floodplain_length],
            #      [floodplain_width/2. + chan_width/2-bankwidth,0.]]
            #}

regionPtAreas = [ [ 0.01, 0.01, 0.5*l0*l0]]

#regionPtAreas=[[0.01, 0.01, 0.5*l0*l0],
#               [floodplain_width/2.-chan_width/2.+0.01, 0.01, 0.5*l0*l0*0.25],
#               [floodplain_width/2.-chan_width/2.+bankwidth+0.01, 0.01, 0.5*l0*l0],
#               [floodplain_width/2.+chan_width/2.-bankwidth+0.01, 0.01, 0.5*l0*l0*0.25],
#               [floodplain_width/2.+chan_width/2.+0.01, 0.01, 0.5*l0*l0]]

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

domain.set_flow_algorithm('DE_SG')
#domain.set_flow_algorithm('DE0')
#domain.set_flow_algorithm('DE2')

domain.set_name('flat_basin') # Output name
domain.set_store_vertices_uniquely(True)


## Function for topography
#def topography(x, y):
#    # Longitudinally sloping floodplain with channel in centre
#    elev1= -y*floodplain_slope - chan_bankfull_depth*\
#            (x>(floodplain_width/2. - chan_width/2.))*\
#            (x<(floodplain_width/2. + chan_width/2.)) 
#    # Add banks
#    if(bankwidth>0.0):
#        leftbnk = floodplain_width/2. - chan_width/2.
#        rightbnk = floodplain_width/2. + chan_width/2.
#        # Left bank
#        elev2 = elev1 + (chan_bankfull_depth \
#                - chan_bankfull_depth/bankwidth*(x - leftbnk))*\
#                (x>leftbnk)*(x < leftbnk + bankwidth)
#        # Right bank
#        elev2 = elev2 + (chan_bankfull_depth \
#                + chan_bankfull_depth/bankwidth*(x - rightbnk))*\
#                (x>rightbnk-bankwidth)*(x < rightbnk)
#    
#    if(bankwidth==0.0):
#        elev2 = elev1
#
#    # Add a wall around the domain
#    #sides = ( ( floodplain_width - x< l0) + (x < l0) + (floodplain_length -y < l0) + (y<l0)).nonzero()[0]
#    #elev2[sides]=3.0
#
#    return elev2

def topography(x,y):
    # 2 Basins with a narrow gap connecting them
    elev=x*0
    elev = elev+3.*((y>(floodplain_length/2.-100.))*(y<(floodplain_length/2+100.)))*((x<floodplain_width/2.-chan_width/2.)+(x>floodplain_width/2.+chan_width/2.))
    return elev

# Very low friction [frictionless would be ideal for testing, but can't be done]
def frictionFun(x,y):
    return 3.0e-02+x*0.

def stageFun(x,y):
    return 1.0*(y<floodplain_length/2.) + 0.
    #topo = topography(x,y)
    #return topo+0.001+1.5*((x> (floodplain_width/3.))*(x< (floodplain_width*2./3.))*(y>(floodplain_length/3.))*(y<(floodplain_length*2./3.)))
    #return topo.max()+1.0+x*0.

#
from anuga.utilities import subGridUtil
mean_topography = subGridUtil.make_spatially_averaged_Fun(topography, domain)

#domain.set_quantity('elevation', topography, location='centroids')
domain.set_quantity('elevation', mean_topography, location='centroids')
domain.set_quantity('stage', stageFun)
domain.set_quantity('xmomentum', 0.)
domain.set_quantity('ymomentum', 0.)
domain.set_quantity('friction', frictionFun)
#domain.set_quantity('height',0.1)

if(domain.flow_algorithm=='DE_SG'):
    domain.set_quantity('alphax',0.)
    domain.set_quantity('alphay',0.)

    # Make subgrid lookup functions
    domain.subgrid_data.make_subgrid_tables(topography, frictionFun, 
        approx_grid_spacing=[0.5,0.5], max_reference_depth=10., 
        constant_subgrid_values=False, 
        #constant_subgrid_values=True, 
        reference_gradient_type='zero')        
        #reference_gradient_type='bed-slope')        

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
flow_in_yval=10.0
if True:
    line1 = [ [floodplain_width/2. - chan_width/2., flow_in_yval],\
              [floodplain_width/2. + chan_width/2., flow_in_yval] \
              ]
    Qin=0.0
    
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

#import pdb
#pdb.set_trace()

#if False:
#    print 'Checking volume quantities'
#
#    print domain.quantities['u_vol'].centroid_values[0:20]/domain.quantities['vol'].centroid_values[0:20]
#
#
#    domain.set_quantity('height',0.5)
#    domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
#    print domain.quantities['u_vol'].centroid_values[0:20]/domain.quantities['vol'].centroid_values[0:20]
#
#    domain.set_quantity('height',1.0)
#    domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
#    print domain.quantities['u_vol'].centroid_values[0:20]/domain.quantities['vol'].centroid_values[0:20]
#
#    domain.set_quantity('height',20.0)
#    domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
#    print domain.quantities['u_vol'].centroid_values[0:20]/domain.quantities['vol'].centroid_values[0:20]


domain.CFL=0.9

for t in domain.evolve(yieldstep=1.0,finaltime=10*60.0):
    #import pdb
    #pdb.set_trace()
    print domain.time, domain.timestep
    if domain.flow_algorithm == 'DE_SG':
        vv = domain.quantities['u_vol'].centroid_values**2+domain.quantities['v_vol'].centroid_values**2
        velsq = (vv/(domain.quantities['vol'].centroid_values**2+1.0e-10))
        print 'vel_max: ', (velsq.max())**0.5
        vloc = velsq.argmax()
        print 'Depth@vel_max: ', domain.quantities['vol'].centroid_values[vloc]/domain.subgrid_data.subgrid_wet_area[vloc]
    #print domain.timestepping_statistics()
    #print 'vol: ', domain.quantities['vol'].centroid_values.sum()
    #print 'BF: ',  domain.get_boundary_flux_integral()
    #print 'Vol - BF: ',  domain.quantities['vol'].centroid_values.sum() - domain.get_boundary_flux_integral()

