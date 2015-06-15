"""

Class to make subgrid lookup tables

Gareth Davies, Geoscience Australia 2014

"""

import os
import scipy
from scipy import log
from scipy import exp
import anuga
import anuga.utilities.spatialInputUtil as su


class SubgridData:

    """
        Hold look-up-table data + interpolation functions to get the sub-grid
        topography information

        -------------------------------------------------------------------
        1) For integrated centroid quantities the class has variables

        - subgrid_centroid_table (a large 1D numpy array -- note each
          variable in the diagram below is a column vector, all are
          concatenated together):

            stage_0
            wet_area_0
            int (h)_0
            int (1/n h^(5/3))_0
            stage_1
            wet_area_1
            int (h)_1
            int (1/n h^(5/3))_1
            ... 
            ... 
            ... 
            ... 
            stage_(nc-1)
            wet_area(nc-1)
            int (h)_(nc-1)
            int (1/n h^(5/3))_(nc-1)

          where nc = number of centroids.

          The vectors stage_i, wet_area_i, int(h)_i, int(1/n h^(5/3))_i
          define the lookup table for centroid i. 

          They are all of the same length (for fixed i), but their lengths
          may differ to the lengths of other cells lookup table data.

          KEY IDEA: ALL THE _i VECTORS CAN BE THOUGHT OF AS REPRESENTING
          SEPARATE ROWS OF A TABLE FOR cell_i -- ALTHOUGH FOR EFFICIENT
          MEMORY ACCESS WE STORE THE WHOLE LOOKUP TABLE AS A 1D ARRAY
          Explanation: The only alternative I see is to store as a table
          with all stages on first row, all wet_areas on 2nd row, etc.
          This has bad memory access properties -- I tried it at first and
          the code was too slow


        - subgrid_centroid_var_names - list of length 4, descriptors of
          variables above


        - subgrid_centroid_starting_index = 1d array of length nc
          Gives the first index in the subgrid_centroid_table corresponding
          to cell_i. (Will be the index of smallest stage value for cell_i)

        - subgrid_centroid_last_lookup_col = 1d array of length nc.
          The last stage index we looked up for each cell, offset from the
          first index for the cell (used to speed up the table search)

             last_stage_index_looked_up = subgrid_centroid_starting_index +
                                          subgrid_centroid_last_lookup_col

          The term 'col' here makes sense if we imagine the data for cell_i
          is a table with stage_i, wet_area_i, .. each on a different row


        - subgrid_centroid_i_ncol = 1d array of length nc. The length of
          stage_i (or any of the other cell i variables)

          NOTE: This can be used to efficiently look up other variables in 
          the table -- e.g. 
          last_wet_area_index_looked_up = subgrid_centroid_starting_index +
                                          subgrid_centroid_last_lookup_col 
                                          + 1*subgrid_centroid_i_ncol
          last_int(h)_index_looked_up = subgrid_centroid_starting_index +
                                          subgrid_centroid_last_lookup_col 
                                          + 2*subgrid_centroid_i_ncol
          last_int(1/n h^(5/3))_index_looked_up = 
                                          subgrid_centroid_starting_index +
                                          subgrid_centroid_last_lookup_col 
                                          + 3*subgrid_centroid_i_ncol

          The term 'col' here makes sense if we imagine the data for cell_i
          is a table with stage_i, wet_area_i, .. each on a different row


        - subgrid_wet_area = 1d array of length nc. Stores the
          wet areas (useful to avoid making into a quantity)



        ######################################################################
        ##
        ## DETAILS ON EDGE LOOKUP TABLE QUANTITIES BELOW HERE
        ##
        ######################################################################


        2) For integrated edge quantities we have (similar to centroids)

        - subgrid_edge_table (large 1d array -- note each variable is a
          column vector, and there are 2 new variables compared with the 
          centroid case):

          stage_0,              
          wet_length_0,         
          int (h)_0 ,           
          int (1/n h^(5/3))_0 , 
          int (1/n^2 h^(7/3))_0 , 
          int (h^(2))_0 ,       
          stage_1,             
          wet_length_1,        
          int (h)_1,           
          int (1/n h^(5/3))_1, 
          int (1/n^2 h^(7/3))_1, 
          int (h^(2))_1,       
          ... 
          ... 
          ... 
          ... 
          ... 
          ... 
          stage_(neu-1)
          wet_length(neu-1)
          int (h)_(neu-1)
          int (1/n h^(5/3))_(neu-1)
          int (1/n^2 h^(7/3))_(neu-1)
          int (h^(2))_(neu-1)

          where neu = number of UNIQUE edge tables
             [ neu can be = 3 x (number of centroids)
              if edge tables are different when viewed from the left/right 
              triangle,

              OR, 

              neu can be: (number of unique edges) if the same edge table 
              can be used for the left/right triangles ]

          The vectors stage_i, wet_length_i, int(h)_i, etc define the lookup
          table for edge_i.
          They are all of the same length (for fixed i), but their lengths
          may differ to the lengths of other edge lookup table data.

          KEY IDEA: ALL THE _i VECTORS CAN BE THOUGHT OF AS REPRESENTING
          SEPARATE ROWS OF A TABLE FOR edge_i -- ALTHOUGH FOR EFFICIENT
          MEMORY ACCESS WE STORE THE WHOLE LOOKUP TABLE AS A 1D ARRAY


        - subgrid_edge_var_names -- list of length 6, descriptors of
          variables above


        - subgrid_edge_starting_index = array of length 3*nc
          First index in edge table corresponding to edge i
          Note there will be repeated values if edges are shared


        - subgrid_edge_last_lookup_col = array of length 3*nc
          Entry of the last stage index we looked up for this edge, as
          offset from the first column for the edge

          last_stage_index_looked_up = subgrid_edge_starting_index +
                                       subgrid_edge_last_lookup_col


        - subgrid_edge_i_ncol = array of length 3*nc
          The number of columns of the lookup table related to edge i

          last_wet_length_index_looked_up = subgrid_edge_starting_index +
                                            subgrid_edge_last_lookup_col +
                                            1*subgrid_edge_i_ncol
          last_int(h)_index_looked_up = subgrid_edge_starting_index +
                                        subgrid_edge_last_lookup_col +
                                        2*subgrid_edge_i_ncol
          last_int(1/n h^(5/3))_index_looked_up = 
                                          subgrid_edge_starting_index +
                                          subgrid_edge_last_lookup_col +
                                          3*subgrid_edge_i_ncol
          last_int(1/n^2 h^(7/3))_index_looked_up = 
                                          subgrid_edge_starting_index +
                                          subgrid_edge_last_lookup_col +
                                          4*subgrid_edge_i_ncol
          last_int(h^(2))_index_looked_up = 
                                          subgrid_edge_starting_index +
                                          subgrid_edge_last_lookup_col +
                                          5*subgrid_edge_i_ncol

    """

    def __init__(self, domain):
        """Create dummy look-up tables:


        """

        self.domain = domain

        nc = len(domain.centroid_coordinates[:, 0])

        # To count unique edges, note that shared edges
        # have domain.neighbours > = 0, while boundary edges
        # do not
        neu = (domain.neighbours >= 0).sum() / 2. +\
              (domain.neighbours < 0).sum()
        assert round(neu) == neu, 'Problem counting edges'

        self.num_centroid_var = 4
        self.num_edge_var = 6

        # These are the variables in the lookup table rows
        self.subgrid_centroid_var_names =\
            ['stage', 'wet_area', 'int[h]', 'int[1/n h^(5/3)]']

        self.subgrid_edge_var_names =\
            ['stage', 'wet_length', 'int[h]', 'int[1/n h^(5/3)]',
             'int[1/n^2 h^(7/3)]', 'int[h^2]']

        # The lookup tables can be large, and must be contiguous to go to C
        # Here we use dummy values
        self.subgrid_centroid_table = scipy.ascontiguousarray(
            scipy.zeros(shape=0))
        self.subgrid_edge_table = scipy.ascontiguousarray(
            scipy.zeros(shape=0))

        # Initialise these to negative (unset flag)
        self.subgrid_centroid_starting_index = scipy.zeros(nc).astype(int) - 1
        self.subgrid_edge_starting_index = scipy.zeros(3 * nc).astype(int) - 1

        # Reasonable default values
        self.subgrid_centroid_last_lookup_col = scipy.zeros(nc).astype(int)
        self.subgrid_edge_last_lookup_col = scipy.zeros(3 * nc).astype(int)

        # Initialise these to negative (unset flag)
        self.subgrid_centroid_i_ncol = scipy.zeros(nc).astype(int) - 1
        self.subgrid_edge_i_ncol = scipy.zeros(3 * nc).astype(int) - 1

        # Defaults
        self.areas = domain.areas
        self.subgrid_wet_area = 1.0 * (self.areas)

        return

###############################################################################

    def set_subgrid_volume_quantities_from_reference_quantities(self):
        """
        Update vol, u_vol, v_vol using the reference quantities + interpolation
        tables

        """
        from swSG1_domain_ext import \
            set_subgrid_volume_quantities_from_reference_quantities as ssvqfrq

        ssvqfrq(self.domain)


        # Check sensible
        assert scipy.all(self.domain.quantities['vol'].centroid_values >= 0.)

        return

###############################################################################

    def set_quantities_from_subgrid_volume_quantities(self):
        """
        Update centroid values of stage / height / alphax / alphay from the
        subgrid volume quantities
        """
        from swSG1_domain_ext import \
            set_quantities_from_subgrid_volume_quantities as sqfsvq

        sqfsvq(self.domain)

        return

###############################################################################

    def make_subgrid_tables(
        self,
        elevation_function,
        friction_function,
        approx_grid_spacing=[1., 1.],
        chunk_size=1e+04,
        verbose=True,
        max_reference_depth=2.0e+04,
        reference_gradient_type='zero',
        constant_subgrid_values=False,
    ):
        """Given a function elevation_function to define the elevations, return
        set of functions describing sub-grid topography

        It does this by generating a grid near the mesh triangle, with points
        spaced ~ approx_grid_spacing, then finding the points inside the mesh
        triangle, computing elevation_function at each, and using these to make
        the sub-grid functions

        INPUTS: @param self -- a subgrid_data object
                @param elevation_function -- the function f(x,y) defining the
                    elevation
                @param friction_function -- the function f(x,y) defining manning's n
                @param approx_grid_spacing -- averaging is computed from points
                    in each triangle, generated by
                    anuga.utilities.SpatialInputUtil.gridPointsInPolygon
                    , with the value of approx_grid_spacing passed there
                @param chunk_size -- Number of mesh triangles to work on in each
                    call to elevation_function. A suitably large chunk_size can
                    reduce function call overhead for some elevation_functions, 
                    but might have consume lots of memory if there are many 
                    grid-points in each triangle
                @param verbose -- print information
                @param max_reference_depth -- extend the lookup table so the
                    'largest reference depth = max_reference_depth'. Default
                    value should be ok for water flows on earth
                @param reference_gradient_type -- how to compute the
                    'reference' water surface gradient -- values are 'zero' or
                    'bed-slope'
                @param constant_subgrid_values -- If TRUE, make each large 
                    cell have constant elevation/friction. This removes the
                    subgrid effect, and can be useful for debugging

        OUTPUT:
            Builds the subgrid lookup table information

        """

        domain = self.domain

        chunk_size = int(chunk_size)

        xc = domain.centroid_coordinates[:, 0]
        yc = domain.centroid_coordinates[:, 1]

        reference_elevation = domain.quantities['elevation'].centroid_values

        # Set the assumed stage gradient used for subgrid computations
        if(reference_gradient_type == 'zero'):

            reference_gradient = domain.centroid_coordinates * 0.
            # If the assumed stage gradient is zero then
            # each edge can have a single lookup table
            share_edge_neighbours = True

        elif(reference_gradient_type == 'bed-slope'):

            domain.quantities['elevation'].compute_gradients()
            # Reshape it
            reference_gradient = scipy.vstack(
                [domain.quantities['elevation'].x_gradient,
                 domain.quantities['elevation'].y_gradient]
            ).transpose()

            # If there is a nonzero reference gradient, then each edge needs a
            # separate lookup table for each neighbouring triangle
            share_edge_neighbours = False
        else:

            raise Exception, 'reference gradient type not recognized'

        # Clear existing subgrid data
        self.subgrid_centroid_table = scipy.zeros(shape=(0))
        self.subgrid_edge_table = scipy.zeros(shape=(0))

        # Set starting_index to negative [flags that cell is unset]
        self.subgrid_centroid_starting_index = \
            0 * self.subgrid_centroid_starting_index - 1
        self.subgrid_edge_starting_index = \
            0 * self.subgrid_edge_starting_index - 1

        self.subgrid_centroid_last_lookup_col *= 0
        self.subgrid_edge_last_lookup_col *= 0

        self.subgrid_centroid_i_ncol *= 0
        self.subgrid_edge_i_ncol *= 0

        # Build subgrid tables
        # Process triangles in chunks to reduce function call overhead
        lx = len(xc)
        chunk_count = scipy.ceil(lx * 1. / (1. * chunk_size)).astype(int)
        for i in range(chunk_count):

            # Compute subgrid tables in triangles lb:ub
            lb = i * chunk_size
            ub = min((i + 1) * chunk_size, lx)

            if verbose:
                print 'Computing functions in triangles ', lb, '-', ub - 1

            # Store x, y, triangle_index, edgeflag for values inside the cell
            px = []
            py = []
            p_index = []
            # p_type contains flag denoting the point is associated with a 
            # centroid (-1), or edge0 (0), or edge1 (1), or edge2 (2)
            p_type = []

            # Loop over this 'chunk' of cells to get the coordinates -- then
            # pass them to the elevation / friction functions in one hit
            for j in range(lb, ub):

                # Compute x,y, triangleindex for values INSIDE the cell
                # We will deal with values on the edges later
                vertex_coords = \
                    domain.mesh.vertex_coordinates[range(3 * j, 3 * j + 3), :]

                # Make grid of points in the triangle
                pts = su.gridPointsInPolygon(
                    vertex_coords, approx_grid_spacing=approx_grid_spacing)

                px = px + [pts[:, 0]]
                py = py + [pts[:, 1]]
                p_index = p_index + [scipy.repeat(j, len(pts[:, 0]))]
                # 'Centre' values get type=-1
                p_type = p_type + [scipy.repeat(-1, len(pts[:, 0]))]

                # Do the same for points along each edge -- these are required
                # for computing both edge AND centroid tables
                for e in [0, 1, 2]:
                    edge_vertices = [0, 1, 2]
                    # Edge 0, uses vertices 1,2
                    # Edge 1, uses vertices 0,2
                    # Edge 2, uses vertices 0,1
                    edge_vertices.remove(e)

                    nx = int(max(scipy.ceil(
                        domain.edgelengths[j, e] / approx_grid_spacing[0]), 3))

                    p_e_x = scipy.linspace(vertex_coords[edge_vertices[0], 0],
                                           vertex_coords[edge_vertices[1], 0],
                                           num=nx)

                    p_e_y = scipy.linspace(vertex_coords[edge_vertices[0], 1],
                                           vertex_coords[edge_vertices[1], 1],
                                           num=nx)
                    px = px + [p_e_x]
                    py = py + [p_e_y]
                    p_index = p_index + [scipy.repeat(j, len(p_e_x))]

                    # Edge e gets type = e
                    p_type = p_type + [scipy.repeat(e, len(p_e_x))]

            # Convert px, py, p_index, p_type to scipy arrays
            px = scipy.hstack(px)
            py = scipy.hstack(py)
            p_index = scipy.hstack(p_index)
            p_type = scipy.hstack(p_type)

            # Get function values at all px,py
            if verbose:
                print '  Evaluating function at ', len(px), ' points'
            all_elevation = elevation_function(px, py)
            all_friction = friction_function(px, py)

            # Can have problems if e.g. the friction function returns a
            # constant
            msg = 'Output of topography function and friction function are' +\
                ' not the same length. Ensure both return a vector'
            assert len(all_elevation) == len(all_friction), msg

            # Compute the lookup table for each triangle and adjust the subgrid
            # data appropriately
            for j in range(lb, ub):

                cell_indices = (p_index == j).nonzero()[0]

                assert len(cell_indices) > 0

                (centroid_table, edge_tables) = \
                    make_subgrid_table_for_one_cell(
                        reference_elevation[j],
                        reference_gradient[j, :],
                        max_reference_depth,
                        domain.centroid_coordinates[j, 0],
                        domain.centroid_coordinates[j, 1],
                        self.areas[j],
                        domain.edgelengths[j, :],
                        px[cell_indices],
                        py[cell_indices],
                        all_elevation[cell_indices],
                        all_friction[cell_indices],
                        p_type[cell_indices],
                        constant_subgrid_values=constant_subgrid_values)

                # Append data to centroid table + index arrays
                self.subgrid_centroid_starting_index[j] = \
                    int(self.subgrid_centroid_table.shape[0])
                self.subgrid_centroid_i_ncol[j] = int(centroid_table.shape[1])

                # Store as 1d vector
                self.subgrid_centroid_table = \
                    scipy.hstack([self.subgrid_centroid_table,
                                  centroid_table.flatten()])

                # Append data to edge table + index arrays
                for e in [0, 1, 2]:
                    edge_index = 3 * j + e
                    if self.subgrid_edge_starting_index[edge_index] >= 0:
                        # Edge is already set
                        continue
                    edge_table = edge_tables[e]
                    self.subgrid_edge_starting_index[edge_index] = \
                        int(self.subgrid_edge_table.shape[0])
                    self.subgrid_edge_i_ncol[edge_index] = \
                        int(edge_table.shape[1])

                    # Append data to edge table as 1d array
                    self.subgrid_edge_table = \
                        scipy.hstack([self.subgrid_edge_table,
                                      edge_table.flatten()])

                    # Check for neighbour, and point it to the same part of the
                    # table
                    n = domain.neighbours.flatten()[edge_index]
                    if ((n >= 0) and share_edge_neighbours):
                        # Find the neighbour edge index
                        ne = 3 * n + \
                            domain.neighbour_edges.flatten()[edge_index]
                        self.subgrid_edge_starting_index[ne] = \
                            self.subgrid_edge_starting_index[edge_index]
                        self.subgrid_edge_i_ncol[ne] = \
                            self.subgrid_edge_i_ncol[edge_index]

        if share_edge_neighbours is False:
            if verbose:
                print 'Defragmenting edge table'

            # Defragment table so that neighbouring edges are next to each other
            # This means the code doesn't need to make large jumps often when
            # accessing the subgrid tables.

            new_subgrid_edge_table = self.subgrid_edge_table * 0.

            new_subgrid_edge_starting_index = \
                self.subgrid_edge_starting_index * 0 - 1

            # new_subgrid_edge_i_ncol = \
            #    self.subgrid_edge_i_ncol*0 - 1

            # Loop over all edges and populate table in order we will access
            new_edge_table_counter = 0
            edge_table_nrow = edge_tables[0].shape[0]
            for k in range(len(xc)):
                for e in [0, 1, 2]:

                    edge_index = 3 * k + e

                    if new_subgrid_edge_starting_index[edge_index] == -1:
                        old_si = self.subgrid_edge_starting_index[edge_index]
                        table_ncol = self.subgrid_edge_i_ncol[edge_index]
                        old_ei = old_si + table_ncol * edge_table_nrow
                        edge_table = self.subgrid_edge_table[old_si:old_ei]

                        # Add to new table
                        new_si = new_edge_table_counter
                        new_ei = new_si + table_ncol * edge_table_nrow
                        new_subgrid_edge_table[new_si:new_ei] = edge_table
                        new_subgrid_edge_starting_index[edge_index] = new_si
                        new_edge_table_counter += table_ncol * edge_table_nrow

                        # Add neighbour if it exists
                        n = domain.neighbours.flatten()[edge_index]

                        if(n >= 0):
                            edge_index = 3 * n + \
                                domain.neighbour_edges.flatten()[edge_index]
                            assert new_subgrid_edge_starting_index[
                                edge_index] == -1

                            old_si = self.subgrid_edge_starting_index[
                                edge_index]
                            table_ncol = self.subgrid_edge_i_ncol[edge_index]
                            old_ei = old_si + table_ncol * edge_table_nrow
                            edge_table = self.subgrid_edge_table[old_si:old_ei]

                            # Add to new table
                            new_si = new_edge_table_counter
                            new_ei = new_si + table_ncol * edge_table_nrow
                            new_subgrid_edge_table[new_si:new_ei] = edge_table
                            new_subgrid_edge_starting_index[
                                edge_index] = new_si
                            new_edge_table_counter += table_ncol * \
                                edge_table_nrow

            self.subgrid_edge_table = new_subgrid_edge_table
            self.subgrid_edge_starting_index = new_subgrid_edge_starting_index
        # End Defragmentation

        # Force large arrays to be contiguous (needed for C)
        self.subgrid_centroid_table = \
            scipy.ascontiguousarray(self.subgrid_centroid_table)
        self.subgrid_edge_table = \
            scipy.ascontiguousarray(self.subgrid_edge_table)

        return

#
# CLASS ENDS HERE
###############################################################################


def make_subgrid_table_for_one_cell(
        reference_elevation,
        reference_gradient,
        max_reference_depth,
        centroid_x,
        centroid_y,
        cell_area,
        edge_lengths,
        x, y,
        all_elevation, all_friction,
        p_type,
        constant_subgrid_values=False):
    """

    Make the subgrid lookup table for a single cell, given a reference
    elevation, the x,y,topography,friction data

    @param reference_elevation = The 'reference centroid elevation' in the cell
        [probably the average or median elevation?]
    @param reference_gradient = The nominal stage gradient in the cell used for
        volume lookup computations [probably 0 or some local approximation
        of the bed slope]
    @param max_reference_depth = The maximum depth of flooding in the cell 
        supported by the lookup table [suggest a large number]
    @param centroid_x
    @param centroid_y
    @param cell_area
    @param edge_lengths = lengths of edges 0,1,2 in order
    @param x = The subgrid point x-coordinates
    @param y = "" y coordinates
    @param all_elevation = The elevation values at the subgrid coordinates
    @param all_friction = The friction values at the subgrid coordinates
    @param p_type = Flag saying whether each subgrid point is an interior point
        (-1) or edge 0 [0] or edge 1 [1] or edge 2 [2]
    @param constant_subgrid_values = If True, assume all subgrid data has the
        same value (reference_elevation) -- useful for debugging (mimics the
        non-subgrid case)
    """

    # Find points inside the cell, and on each edge
    interior_pts = (p_type == -1).nonzero()[0]
    edge0_pts = (p_type == 0).nonzero()[0]
    edge1_pts = (p_type == 1).nonzero()[0]
    edge2_pts = (p_type == 2).nonzero()[0]

    # Avoid modifying inputs
    subgrid_elevation = all_elevation
    subgrid_friction = all_friction

    # constant_subgrid_values can be used to 'not' apply subgrid topography
    # (useful for testing)
    if constant_subgrid_values:
        # reference_elevation
        subgrid_elevation = all_elevation * 0. + all_elevation.mean()
        # reference_friction
        subgrid_friction = all_friction * 0. + all_friction.mean()

    # Each cell has a single reference stage, elevation, and depth, where:
    #
    # reference_depth = reference_stage - reference_elevation; and
    #
    # reference_stage is the 'cell stage' and varies during the simulation; and
    #
    # reference_elevation is fixed (probably the mean cell elevation or similar)
    #
    # stage(x,y) = reference_depth + stage_at_refdepth0 [can go below
    #              subgrid_elevation]
    # depth(x,y) = max( stage(x,y) - subgrid_elevation(x,y), 0.0)

    # Compute stage at x,y when "reference_stage = reference elevation" (so
    # reference_depth=0), and "reference_stage_gradient = reference_gradient"
    # Useful since we can convert from reference_depth to 'real stage'
    stage_at_refdepth0 = reference_elevation +\
        (x - centroid_x) * reference_gradient[0] +\
        (y - centroid_y) * reference_gradient[1]

    # CENTROID TABLE HERE

    # For numerical stability it is important that the 'first cell to get wet',
    # on EITHER the interior or an edge, is included in the lookup tables for
    # the volume. We can 'just' include the smallest subgrid point, or we can
    # include all edge points. The former should be more accurate but the
    # latter might give better numerical behaviour (less 'jumps' in the lookup
    # table)
    include_all_edge_points = True

    if not include_all_edge_points:

        first_flooded_index = (stage_at_refdepth0 - subgrid_elevation).argmax()
        if not (first_flooded_index in interior_pts.tolist()):
            interior_pts = scipy.array(interior_pts.tolist() +
                                       [first_flooded_index])

    else:

        interior_pts = scipy.array(interior_pts.tolist() +
                                   edge0_pts.tolist() +
                                   edge1_pts.tolist() +
                                   edge2_pts.tolist())

    num_centroid_table_col = 10
    upper_depth_offset = 5.  # FIXME: Ad-hoc

    centroid_table = compute_subgrid_table(
        max_reference_depth,
        upper_depth_offset,
        reference_elevation,
        num_centroid_table_col,
        stage_at_refdepth0[interior_pts],
        subgrid_elevation[interior_pts],
        subgrid_friction[interior_pts],
        cell_area, 
        feature_type='volume')

    # EDGE TABLES HERE

    num_edge_table_col = 10

    # upper_depth_offset = 5. # FIXME: Ad-hoc

    edge_pts = [edge0_pts, edge1_pts, edge2_pts]
    edge_tables = [[], [], []]
    for i in range(3):
        edgei_pts = edge_pts[i]
        edgei_length = edge_lengths[i]
        edge_tables[i] = compute_subgrid_table(
            max_reference_depth,
            upper_depth_offset,
            reference_elevation,
            num_edge_table_col,
            stage_at_refdepth0[edgei_pts],
            subgrid_elevation[edgei_pts],
            subgrid_friction[edgei_pts],
            edgei_length, 
            feature_type='edge')

    return centroid_table, edge_tables

###############################################################################


def make_exponentially_spaced_sequence(min_value, max_value, length):
    """Make a sequence from min_value to max_value with exponential spacing and
       specified length

    """
    assert min_value < max_value
    assert length >= 2

    # Make exponentially spaced table depths
    log_spacing = \
        scipy.linspace(0., log(max_value - min_value + 1.),
                       num=length)

    output_sequence = min_value + exp(log_spacing) - 1.

    return output_sequence

###############################################################################


def compute_subgrid_table(
        max_reference_depth,
        upper_depth_offset,
        reference_elevation,
        num_table_col,
        stage_at_refdepth0,
        subgrid_elevation,
        subgrid_friction,
        length_or_area,
        feature_type):
    """ Make the subgrid table for a single edge or volume

        @param max_reference_depth = the largest reference depth to 
            feature in the table. Usually a very large number since the table
            can't be looked up beyond this
        @param upper_depth_offset = the second largest reference depth
            in the table will ensure all points are flooded by at least
            upper_depth_offset. Intuitively, we would not often want
            the cell to be flooded above this (since the table may not
            be accurate above this)
        @param reference_elevation = the cells reference elevation
        @param num_table_col = how many columns in the subgrid table
        @param stage_at_refdepth0 = vector of values for x,y describing
            the subgrid topography and assumed water surface elevation, with 
            stage(x,y) = reference_depth(x,y) + stage_at_refdepth0(x,y)
        @param subgrid_elevation = vector of subgrid elevation values
        @param subgrid_friction = vector of subgrid friction values
        @param length_or_area = edge_length if feature_type = 'edge', or
            triangle_area if feature_type = 'volume'
        @param feature_type = 'edge' or 'volume'

    """

    if feature_type == 'edge':
        num_subgrid_var = 6
    elif feature_type == 'volume':
        num_subgrid_var = 4
    else:
        raise Exception('feature_type must be "edge" or "volume"')

    min_reference_depth = (subgrid_elevation - stage_at_refdepth0).min()
    min_interior_reference_depth = min_reference_depth

    second_max_reference_depth = upper_depth_offset +\
        (subgrid_elevation - stage_at_refdepth0).max()

    assert max_reference_depth > second_max_reference_depth,\
        'max_reference_depth is <= the second max reference depth'

    table_depths = make_exponentially_spaced_sequence(
        min_reference_depth, second_max_reference_depth,
        length=num_table_col - 1)

    # Append the maximum_depth
    table_depths = scipy.hstack([table_depths,
                                 scipy.array([max_reference_depth])])
    table_stages = table_depths + reference_elevation

    subgrid_table = scipy.zeros(shape=(num_subgrid_var, len(table_stages)))

    for i in range(len(table_stages)):

        # (cell depths can be negative)
        cell_depths = table_depths[i] + stage_at_refdepth0 - subgrid_elevation

        # Stage
        subgrid_table[0, i] = table_stages[i]

        # Wetted length = length * wet fraction
        # Include area of cells that are 'almost' wet to make behaviour more
        # similar to non-subgrid in the case of a constant subgrid elevation
        subgrid_table[1, i] = length_or_area *\
            max((cell_depths > -1.0e-06).mean(), 1. / len(cell_depths))

        # Now ensure depth > 0.
        cell_depths = scipy.maximum(cell_depths, 0.)
        # Wetted cross-section
        subgrid_table[2, i] = cell_depths.mean() * length_or_area

        # int (1/n h^(5/3))
        subgrid_table[3, i] = length_or_area *\
            (cell_depths ** (5. / 3.) / subgrid_friction).mean()

        if feature_type == 'edge':
            # int (1/n^2 h^(7/3))
            subgrid_table[4, i] = length_or_area *\
                (cell_depths ** (7. / 3.) / subgrid_friction ** 2).mean()

            # int ( h^(2))
            subgrid_table[5, i] = length_or_area *\
                (cell_depths ** (2.)).mean()

    return subgrid_table
