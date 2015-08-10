"""
Test the subgrid_data code

Gareth Davies, Geoscience Australia, 2014-2015
"""

import unittest


import anuga
import numpy
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.utilities import plot_utils as util
from anuga.config import g
from anuga.utilities import spatialInputUtil as su

import anuga.shallow_water.subgrid_data as subgrid_data

verbose = False

# Various constants which matter
floodplain_slope_val = 1./300.
height_init = 0.1

class Test_subgrid_data(unittest.TestCase):
    """Test the subgrid data functions

    """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def create_subgrid_domain(self, gradtype='zero'):
        """ Make a subgrid domain. Useful for testing it

        """
        #----------------------------------------------------------------------
        # Useful parameters for controlling this case
        #----------------------------------------------------------------------

        floodplain_length = 500.0 # Model domain length
        floodplain_width = 14.0 # Model domain width
        floodplain_slope = floodplain_slope_val
        chan_bankfull_depth = 1.0 # Bankfull depth of the channel
        chan_width = 10. # Bankfull width of the channel
        bankwidth = 2. # Width of the bank regions -- these protrude inward
        man_n = 0.03 # Manning's n
        l0 = 5.000 # Length scale of triangle side length in channel

        assert chan_width < floodplain_width, \
                ' ERROR: Channel width is greater than floodplain width'

        assert bankwidth < chan_width/2., \
                'ERROR: The bank width must be less than half the channel width'

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        # Define boundary polygon -- in this case clockwise around the proposed
        # boundary
        boundary_polygon = [ 
            [0.,0.], 
            [0., floodplain_length], 
            [floodplain_width/2. - chan_width/2., floodplain_length], 
            [floodplain_width/2. + chan_width/2., floodplain_length], 
            [floodplain_width, floodplain_length], 
            [floodplain_width, 0.], 
            [floodplain_width/2. + chan_width/2., 0.], 
            [floodplain_width/2. - chan_width/2., 0.] 
                             ]

        # Define channel polygon, where we can optionally refine the
        # resolution.  Note that we set this a distance l0 inside the boundary
        # polygon, so the polygons do not overlap. 

        breakLines={
            'n1':[[floodplain_width/2.- chan_width/2., floodplain_length],
                 [floodplain_width/2. - chan_width/2.,0.]],
            'n2':[[floodplain_width/2.+ chan_width/2., floodplain_length],
                  [floodplain_width/2. + chan_width/2.,0.]],
            # These ones are inside the channel 
            'n3':[[floodplain_width/2.- chan_width/2.+bankwidth, floodplain_length],
                  [floodplain_width/2. - chan_width/2.+bankwidth,0.]],
            'n4':[[floodplain_width/2.+ chan_width/2.-bankwidth, floodplain_length],
                  [floodplain_width/2. + chan_width/2-bankwidth,0.]]
                   }

        regionPtAreas=[
            [0.01, 0.01, 0.5*l0*l0],
            [floodplain_width/2.-chan_width/2.+0.01, 0.01, 0.5*l0*l0*0.25],
            [floodplain_width/2.-chan_width/2.+bankwidth+0.01, 0.01, 0.5*l0*l0],
            [floodplain_width/2.+chan_width/2.-bankwidth+0.01, 0.01, 0.5*l0*l0*0.25],
            [floodplain_width/2.+chan_width/2.+0.01, 0.01, 0.5*l0*l0]]

        # Define domain with appropriate boundary conditions
        anuga.create_mesh_from_regions(
            boundary_polygon,
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
            verbose=verbose)

        domain=anuga.create_domain_from_file('channel_floodplain1.msh')

        domain.set_flow_algorithm('DE_SG')

        domain.set_name('channel_floodplain1')
        domain.set_store_vertices_uniquely(True)

        # Function for topography
        def topography(x, y):
            # Longitudinally sloping floodplain with channel in centre
            elev1= -y*floodplain_slope - chan_bankfull_depth*\
                    (x>(floodplain_width/2. - chan_width/2.))*\
                    (x<(floodplain_width/2. + chan_width/2.)) 

            # Add banks
            if(bankwidth > 0.0):
                leftbnk = floodplain_width/2. - chan_width/2.
                rightbnk = floodplain_width/2. + chan_width/2.
                # Left bank
                elev2 = elev1 + (chan_bankfull_depth \
                        - chan_bankfull_depth/bankwidth*(x - leftbnk))*\
                        (x>leftbnk)*(x < leftbnk + bankwidth)
                # Right bank
                elev2 = elev2 + (chan_bankfull_depth \
                        + chan_bankfull_depth/bankwidth*(x - rightbnk))*\
                        (x>rightbnk-bankwidth)*(x < rightbnk)

            if(bankwidth == 0.0):
                elev2 = elev1

            return elev2

        def frictionFun(x,y):
            return man_n + x*0.

        def stageFun(x,y):
            return topography(x,y)+height_init

        domain.set_quantity('elevation', topography, location='centroids')
        domain.set_quantity('stage', stageFun)
        domain.set_quantity('xmomentum', 0.)
        domain.set_quantity('ymomentum', 0.)
        domain.set_quantity('friction', frictionFun)
        domain.set_quantity('height', height_init)
        domain.set_quantity('alphay', floodplain_slope**0.5)
        domain.set_quantity('alphax', 0.)

        domain.subgrid_data.make_subgrid_tables(
            topography, frictionFun, approx_grid_spacing=[1.5,1.5],
            reference_gradient_type=gradtype,verbose=verbose)

        return domain


    def table_check(self, vartype = 'edge', reference_gradient_type='zero'):
        """Check we can get the right table values

            The code is written with variable names as though edges
            are being tested -- but it applies also to centroid tables
        """

        # Define properties on a hypothetical edge

        reference_elevation = 1.
        upper_depth_offset = 10.
        max_reference_depth = 100.
        num_table_col = 22
        length_or_area = 100. # This is interpreted as the triangle area if vartype=='centroid'

        real_elevation = numpy.linspace(-5.,5.,num=50)+reference_elevation
        real_friction = real_elevation*0. + 0.03
        # Assume 'flat' free surface
        if reference_gradient_type=='zero':
            stage_at_refdepth0 = reference_elevation + 0.*real_elevation
        elif reference_gradient_type=='bed-slope':
            stage_at_refdepth0 = real_elevation*1.0

        if vartype=='edge':

            edge_table = subgrid_data.compute_subgrid_table(
                max_reference_depth,
                upper_depth_offset,
                reference_elevation,
                num_table_col,
                stage_at_refdepth0,
                real_elevation,
                real_friction,
                length_or_area,
                feature_type='edge')

        elif vartype=='centroid':

            edge_table = subgrid_data.compute_subgrid_table(
                max_reference_depth,
                upper_depth_offset,
                reference_elevation,
                num_table_col,
                stage_at_refdepth0,
                real_elevation,
                real_friction,
                length_or_area,
                feature_type='volume')


        if vartype == 'edge':
            assert edge_table.shape[0] == 6
        elif vartype == 'centroid':
            assert edge_table.shape[0] == 4
        
        assert edge_table.shape[1] == num_table_col

        # Check stage is increasing
        stage_copy = edge_table[0,:]*1.0
        stage_copy.sort()
        assert numpy.all(stage_copy == edge_table[0,:])

        # Check stage range
        assert numpy.allclose(edge_table[0,0],
            (real_elevation-stage_at_refdepth0).min() + reference_elevation)
        assert numpy.allclose(edge_table[0,-1],
            reference_elevation + max_reference_depth )

        ## FIXME: We have been changing the second-last table value, so I
        ##        comment this out for now, but add it back if a method is
        ##        finalised.
        #assert numpy.allclose(edge_table[0,-2], 
        #    (real_elevation-stage_at_refdepth0).max() +\
        #     reference_elevation + upper_depth_offset)

        # Smallest edge length should have at least one subgrid point wet
        min_length_error = \
            (edge_table[1,0] >= 1./len(real_elevation)*length_or_area)
        assert min_length_error
        assert max(edge_table[1,:]) == length_or_area

        # Everything else should have first column = 0
        assert all(edge_table[2:-1, 0]==0.)

        # Compute the integrals directly and compare with subgrid table values,
        # for some particular cases
        for i in [-1, -2, 5, 10]:
            real_stage = \
                (edge_table[0,i] - reference_elevation) + stage_at_refdepth0
            depths = numpy.maximum( real_stage - real_elevation, 0.)

            # Check edge length
            assert numpy.allclose(edge_table[1,i], 
                length_or_area*(depths>0.).mean())

            # Check int(h)
            assert numpy.allclose(edge_table[2,i],
                length_or_area*depths.mean())

            # Check int( 1/n * h^{5/3} )
            int_val = (depths**(5./3.)/real_friction).mean()*length_or_area 
            assert abs(edge_table[3,i] - int_val) < 1.0e-06

            if vartype=='edge':
                # Check int( 1/n^2 * h^{7/3} )
                int_val = (depths**(7./3.)/real_friction**2).mean()*length_or_area 
                assert abs(edge_table[4,i] - int_val) < 1.0e-06

                # Check int( h^{2} )
                int_val = (depths**(2.)).mean()*length_or_area 
                assert abs(edge_table[5,i] - int_val) < 1.0e-06

        return


    def test_subgrid_edge_table(self):
        """Check we can get the right table values

        """

        self.table_check('edge', 'zero')
        self.table_check('edge', 'bed-slope')
        return


    def test_subgrid_centroid_table(self):
        """Check we can get the right table values

        """
        self.table_check('centroid', 'zero')
        self.table_check('centroid', 'bed-slope')

        return

    def test_make_exponentially_spaced_sequence(self):
        """
        """
        from scipy import log
        from scipy import exp

        x = subgrid_data.make_exponentially_spaced_sequence(-0.1, 50.2, 20)
        assert numpy.allclose(x.max() , 50.2)
        assert numpy.allclose(x.min() , -0.1)
        assert len(x) == 20

        # Check constant log spacing
        log_spacing = log(x - x.min()+1.)
        log_diff = log_spacing[1:] - log_spacing[0:-1]
        assert abs(log_diff.max() - log_diff.min())<1.0e-05

        return


    def test_make_subgrid_table_for_one_cell(self):
        """ Tests that the centroid and edge tables are correctly produced for a
            a single cell

        """

        reference_gradient = [0., 0.]
        max_reference_depth = 10000.

        # Ideal Triangle vertices
        myTri = numpy.array([[0., 0.], [10., 0.], [0., 10.]])
        # Suppose elevation = x+y over this triangle
        reference_elevation = (0.+10.+10.)/3.
        cell_area = 0.5*10.*10.
        centroid_x = (0.+10.+0.)/3.
        centroid_y = (0.+10.+0.)/3.

        edge_lengths = numpy.array([(200.)**0.5, 10., 10.])      

        # Get points inside triangle
        pts = su.gridPointsInPolygon(myTri,approx_grid_spacing=[0.5,0.5])
        x = pts[:,0]
        y = pts[:,1]
        pType = numpy.repeat(-1,len(x))

        # Get edge points
        ne = int(numpy.ceil(edge_lengths[0]/2))
        e0x = numpy.linspace(myTri[1,0],myTri[2,0],ne)
        e0y = numpy.linspace(myTri[1,1],myTri[2,1],ne)
        e0z = e0x+e0y
        
        ne = int(numpy.ceil(edge_lengths[1]/2))
        e1x = numpy.linspace(myTri[0,0],myTri[2,0],ne)
        e1y = numpy.linspace(myTri[0,1],myTri[2,1],ne)
        e1z = e1x+e1y
        
        ne = int(numpy.ceil(edge_lengths[2]/2))
        e2x = numpy.linspace(myTri[0,0],myTri[1,0],ne)
        e2y =numpy.linspace(myTri[0,1],myTri[1,1],ne)
        e2z = e2x+e2y

        # Collate all points
        x = numpy.hstack([x,e0x,e1x,e2x])
        y = numpy.hstack([y,e0y,e1y,e2y])
        pType = numpy.hstack([pType, numpy.repeat(0,len(e0x)), 
            numpy.repeat(1,len(e1x)), numpy.repeat(2,len(e2x))])

        # Friction + topography on the triangle z=x+y
        all_topo = x+y
        all_friction = all_topo*0.+0.03

        centroid_table, edge_tables = \
            subgrid_data.make_subgrid_table_for_one_cell(
                reference_elevation,reference_gradient, 
                max_reference_depth,\
                centroid_x,centroid_y,cell_area, edge_lengths,\
                x,y, all_topo, all_friction, pType)

        # Check ranges of reference depth
        assert centroid_table[0,:].max() == max_reference_depth +\
             reference_elevation

        # Check that peak area = cell area
        assert numpy.allclose(centroid_table[1,-1], cell_area)
        assert numpy.allclose(edge_tables[0][1,-1], edge_lengths[0])
        assert numpy.allclose(edge_tables[1][1,-1], edge_lengths[1])
        assert numpy.allclose(edge_tables[2][1,-1], edge_lengths[2])

        # Edge 1 vars should be the same as edge2 vars
        assert numpy.allclose(edge_tables[1], edge_tables[2])

        # Check the minum stage on each edge
        # (Min elevation = 10 on edge 0)
        assert numpy.allclose(edge_tables[0][0,0], 10.)
        # (Min elevation = 0 on edge 1)
        assert numpy.allclose(edge_tables[1][0,0], 0.)
        # (Min elevation = 0 on edge 2)
        assert numpy.allclose(edge_tables[2][0,0], 0.)

        #######################################################################

        # Now run with constant_subgrid_values, and check that all edges are
        # the same
        centroid_table, edge_tables = \
            subgrid_data.make_subgrid_table_for_one_cell(
                reference_elevation,reference_gradient, 
                max_reference_depth,\
                centroid_x,centroid_y,cell_area, edge_lengths,\
                x,y, all_topo, all_friction, pType,
                constant_subgrid_values=True)

        # All stage curves should now be the same
        assert numpy.allclose(edge_tables[0][0,:], edge_tables[1][0,:])
        assert numpy.allclose(edge_tables[0][0,:], edge_tables[2][0,:])
        assert numpy.allclose(edge_tables[0][0,:], centroid_table[0,:])

        # Edge 1 2 vars should all be the same
        assert numpy.allclose(edge_tables[1], edge_tables[2])
        # Edge 0 should be like the other 2 edges, rescaled by sqrt(2)
        # (except for stage, of course -- and except for the shallowest
        # 'wet-length' which is set based on the number of points on the edge)
        assert numpy.allclose(edge_tables[0][2:,:]/2.**0.5, 
                              edge_tables[1][2:,:])

        # The centroid table should be like the first few entries of the
        # edge tables, rescaled (except of course for stage -- and not
        # area/length because of the 'near dry' treatment)
        assert numpy.allclose(centroid_table[2:,:]/cell_area,
                              edge_tables[0][2:4,:]/edge_lengths[0])
        # The subgrid tables have been tested elsewhere so this routine
        # does not need to check micro-details of table values
        return


    def test_make_subgrid_tables(self):
        """ Test the class to make the subgrid tables

        """

        for gradtype in ['zero', 'bed-slope']:

            # Flag for whether edge neighbours have the same lookup table
            shared_edge_neighbours = (gradtype=='zero')

            domain = self.create_subgrid_domain(gradtype=gradtype)

            # Useful shorthands
            centroid_table = domain.subgrid_data.subgrid_centroid_table
            centroid_i_ncol = domain.subgrid_data.subgrid_centroid_i_ncol
            centroid_starting_index = \
                domain.subgrid_data.subgrid_centroid_starting_index
            centroid_last_lookup_col = \
                domain.subgrid_data.subgrid_centroid_last_lookup_col
            edge_table = domain.subgrid_data.subgrid_edge_table
            edge_starting_index = \
                domain.subgrid_data.subgrid_edge_starting_index
            edge_i_ncol = domain.subgrid_data.subgrid_edge_i_ncol
            edge_last_lookup_col = \
                domain.subgrid_data.subgrid_edge_last_lookup_col

            # Check that the indexes into the subgrid tables make sense
            assert 4*centroid_i_ncol.sum() == centroid_table.shape[0]

            if shared_edge_neighbours:
                # Must account for repetition of edges
                double_edges = (domain.neighbours.flatten()>=0).nonzero()[0]
                assert 6*(edge_i_ncol.sum() - \
                        edge_i_ncol[double_edges].sum()/2) == \
                    edge_table.shape[0]
            else:
                assert 6*edge_i_ncol.sum() == edge_table.shape[0]

            if shared_edge_neighbours:
                # Check that neighbouring edges point to the same part of the
                # table
                neighbours = domain.neighbours.flatten()
                neighbour_edges = domain.neighbour_edges.flatten()
                for i in range(len(neighbours)):
                    ni = neighbours[i]
                    if ni>=0:
                        nei = ni*3 + neighbour_edges[i]
                        assert edge_i_ncol[i] == edge_i_ncol[nei]
                        assert edge_starting_index[i] == edge_starting_index[nei]
                        assert edge_last_lookup_col[i] == edge_last_lookup_col[nei]

            # Check that all the starting columns have all variables 0 except
            # for stage and area/length
            for i in range(len(centroid_starting_index)):
                for j in range(2,4):
                    start = centroid_starting_index[i] + j*centroid_i_ncol[i]
                    assert numpy.allclose(centroid_table[start], 0.)
                # Check that the wet area is always >=0.
                start = centroid_starting_index[i] + 1*centroid_i_ncol[i]
                end = start + centroid_i_ncol[i]
                assert (centroid_table[start:end] >= 0.).all()
        
            # As above for edges
            for i in range(len(edge_starting_index)):
                for j in range(2,6):
                    start = edge_starting_index[i] + j*edge_i_ncol[i]
                    assert numpy.allclose(edge_table[start], 0.)
                # Check that the wet length is always >=0.
                start = edge_starting_index[i] + 1*edge_i_ncol[i]
                end = start + edge_i_ncol[i]
                assert (edge_table[start:end] >= 0.).all()

            # Check that all subgrid tables are nondecreasing
            for i in range(len(centroid_starting_index)):
                for j in range(4):
                    start = centroid_starting_index[i]+ centroid_i_ncol[i]*j
                    end = centroid_starting_index[i] + centroid_i_ncol[i]
                    difftab = centroid_table[(start+1):end] - \
                              centroid_table[start:(end-1)]
                    if not (difftab>=0.).all():
                        print difftab.min()
                        raise Exception('values below 0')

            # Check that all subgrid tables are nondecreasing
            for i in range(len(edge_starting_index)):
                for j in range(4):
                    start = edge_starting_index[i]+ edge_i_ncol[i]*j
                    end = edge_starting_index[i] + edge_i_ncol[i]
                    difftab = edge_table[(start+1):end] - \
                              edge_table[start:(end-1)]
                    if not (difftab>=0.).all():
                        print difftab.min()
                        raise Exception('values below 0')

            # Crude check that there are more columns in the table than
            # cells/edges Should be true with current implementation
            assert edge_table.shape[0] > len(edge_i_ncol)
            assert centroid_table.shape[0] > len(centroid_i_ncol)

        return


    def test_set_subgrid_volume_quantities_from_reference_quantities(self):
        """Check that we can set the volume quantities based on reference quantities
        """

        domain = self.create_subgrid_domain(gradtype='bed-slope')

        # Set a large height [for this scenario] so vol $\simeq$ height*area
        # The relation won't be exact because the reference elevation != cell mean elevation
        # (because of finite sampling of elevation values)
        #domain.set_quantity('height', 4.0)
        domain.add_quantity('stage', 3.9) # Stage was initially 0.1 above bed
        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
        meanDepth = domain.quantities['vol'].centroid_values/domain.areas

        assert (meanDepth.min() > 3.9)
        assert (meanDepth.max() < 4.1)

        # Try another depth 
        # domain.set_quantity('height', 10.0)
        domain.add_quantity('stage', 6.0)
        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
        meanDepth=domain.quantities['vol'].centroid_values/domain.areas

        assert (meanDepth.min() > 9.9)
        assert (meanDepth.max() < 10.1)

        domain = self.create_subgrid_domain(gradtype='zero')
        # Try again with a different extrapolation
        #domain.set_quantity('height', 10.0)
        domain.add_quantity('stage', 9.9)
        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
        meanDepth = domain.quantities['vol'].centroid_values/domain.areas

        assert (meanDepth.min() > 9.9)
        assert (meanDepth.max() < 10.1)

        return

    def test_set_reference_quantities_from_subgrid_volume_quantities(self):
        """
            Check that we can set the reference quantities based on the volume
        """
        domain = self.create_subgrid_domain(gradtype='bed-slope')

        domain.subgrid_data.set_subgrid_volume_quantities_from_reference_quantities()
        domain.subgrid_data.set_quantities_from_subgrid_volume_quantities()

        dh0 = domain.quantities['height'].centroid_values - height_init
        assert numpy.allclose(dh0, 0.)
        dalphay = domain.quantities['alphay'].centroid_values - (floodplain_slope_val)**0.5
        assert numpy.allclose(dalphay, 0.)
        assert numpy.allclose(domain.quantities['alphax'].centroid_values, 0.)
        depth = domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values
        assert numpy.allclose(depth, height_init)

        return


# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_subgrid_data, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
    
