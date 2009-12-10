#!/usr/bin/env python

import unittest
import sys
import os
from math import sqrt





from anuga.interface import Domain
from anuga.interface import rectangular_cross

from anuga_parallel.distribute_mesh import pmesh_divide_metis
from anuga_parallel.distribute_mesh import build_submesh
from anuga_parallel.distribute_mesh import submesh_full, submesh_ghost, submesh_quantities
from anuga_parallel.distribute_mesh import extract_hostmesh, rec_submesh, send_submesh

from anuga_parallel.interface import myid, numprocs, barrier

import numpy as num


import exceptions
class ParallelException(exceptions.Exception):
    # Call to pypar.abort() somewhere here
    
    def __init__(self, value):
        import pypar
        print 'Terminating parallel processes'
        pypar.abort()



def topography(x,y): 
    return -x/2


def xcoord(x,y):
    return x

def ycoord(x,y):
    return y








###############################################################

class Test_parallel_distribute_mesh(unittest.TestCase):


    def test_distibute_three_processors(self):
        """
        Do a parallel test of distributing a rectangle onto 3 processors
        """

        import pypar

        myid = pypar.rank()
        numprocs = pypar.size()

        assert(numprocs == 3)

        print numprocs

        barrier()

        if myid == 0:

            points, vertices, boundary = rectangular_cross(2,2)

            domain = Domain(points, vertices, boundary)


            domain.set_quantity('elevation', topography) # Use function for elevation
            domain.set_quantity('friction', 0.0)         # Constant friction 
            domain.set_quantity('stage', expression='elevation') # Dry initial stage
            domain.set_quantity('xmomentum', expression='friction + 2.0') # 
            domain.set_quantity('ymomentum', ycoord) #

            #----------------------------------------------------------------------------------
            # Test pmesh_divide_metis
            #----------------------------------------------------------------------------------
            nodes, triangles, boundary, triangles_per_proc, quantities = pmesh_divide_metis(domain,numprocs)

            assert(num.allclose(nodes,points))

            true_vertices = [[0, 9, 1], [3, 9, 0], [4, 9, 3], [1, 9, 4], [1, 10, 2], [4, 10, 1], [5, 10, 4], [2, 10, 5], [3, 11, 4], [6, 11, 3], [7, 11, 6], [4, 11, 7], [4, 12, 5], [7, 12, 4], [8, 12, 7], [5, 12, 8]]

            true_triangles = [[4, 9, 3], [4, 12, 5], [7, 12, 4], [8, 12, 7], [5, 12, 8], [0, 9, 1], [1, 9, 4], [1, 10, 2], [4, 10, 1], [5, 10, 4], [2, 10, 5], [3, 9, 0], [3, 11, 4], [6, 11, 3], [7, 11, 6], [4, 11, 7]]

            assert(num.allclose(vertices,true_vertices))
            assert(num.allclose(triangles,true_triangles))

            assert(num.allclose(triangles_per_proc,[5,6,5]))


            #----------------------------------------------------------------------------------
            # Test build_submesh
            #----------------------------------------------------------------------------------
            submesh = build_submesh(nodes, triangles, boundary, quantities, triangles_per_proc)


            assert(num.allclose(submesh['full_nodes'][0],[[3.0, 0.5, 0.0], [4.0, 0.5, 0.5], [5.0, 0.5, 1.0], [7.0, 1.0, 0.5], [8.0, 1.0, 1.0], [9.0, 0.25, 0.25], [12.0, 0.75, 0.75]]))
            assert(num.allclose(submesh['full_nodes'][1],[[0.0, 0.0, 0.0], [1.0, 0.0, 0.5], [2.0, 0.0, 1.0], [4.0, 0.5, 0.5], [5.0, 0.5, 1.0], [9.0, 0.25, 0.25], [10.0, 0.25, 0.75]]))
            assert(num.allclose(submesh['full_nodes'][2],[[0.0, 0.0, 0.0], [3.0, 0.5, 0.0], [4.0, 0.5, 0.5], [6.0, 1.0, 0.0], [7.0, 1.0, 0.5], [9.0, 0.25, 0.25], [11.0, 0.75, 0.25]]))


            assert(num.allclose(submesh['ghost_nodes'][0],[[0.0, 0.0, 0.0], [1.0, 0.0, 0.5], [2.0, 0.0, 1.0], [6.0, 1.0, 0.0], [10.0, 0.25, 0.75], [11.0, 0.75, 0.25]]))
            assert(num.allclose(submesh['ghost_nodes'][1],[[3.0, 0.5, 0.0], [7.0, 1.0, 0.5], [8.0, 1.0, 1.0], [11.0, 0.75, 0.25], [12.0, 0.75, 0.75]]))
            assert(num.allclose(submesh['ghost_nodes'][2],[[1.0, 0.0, 0.5], [5.0, 0.5, 1.0], [8.0, 1.0, 1.0], [12.0, 0.75, 0.75]]))



            true_full_triangles = [num.array([[ 4,  9,  3],
                                              [ 4, 12,  5],
                                              [ 7, 12,  4],
                                              [ 8, 12,  7],
                                              [ 5, 12,  8]]),
                                   num.array([[ 0,  9,  1],
                                              [ 1,  9,  4],
                                              [ 1, 10,  2],
                                              [ 4, 10,  1],
                                              [ 5, 10,  4],
                                              [ 2, 10,  5]]),
                                   num.array([[ 3,  9,  0],
                                              [ 3, 11,  4],
                                              [ 6, 11,  3],
                                              [ 7, 11,  6],
                                              [ 4, 11,  7]])]


            assert(num.allclose(submesh['full_triangles'][0],true_full_triangles[0]))
            assert(num.allclose(submesh['full_triangles'][1],true_full_triangles[1]))
            assert(num.allclose(submesh['full_triangles'][2],true_full_triangles[2]))

            true_ghost_triangles = [num.array([[ 5,  0,  9,  1],
                                               [ 6,  1,  9,  4],
                                               [ 8,  4, 10,  1],
                                               [ 9,  5, 10,  4],
                                               [10,  2, 10,  5],
                                               [11,  3,  9,  0],
                                               [12,  3, 11,  4],
                                               [13,  6, 11,  3],
                                               [14,  7, 11,  6],
                                               [15,  4, 11,  7]]),
                                    num.array([[ 0,  4,  9,  3],
                                               [ 1,  4, 12,  5],
                                               [ 2,  7, 12,  4],
                                               [ 4,  5, 12,  8],
                                               [11,  3,  9,  0],
                                               [12,  3, 11,  4]]),
                                    num.array([[ 0,  4,  9,  3],
                                               [ 1,  4, 12,  5],
                                               [ 2,  7, 12,  4],
                                               [ 3,  8, 12,  7],
                                               [ 5,  0,  9,  1],
                                               [ 6,  2,  9,  4]])]



            assert(num.allclose(submesh['ghost_triangles'][0],true_ghost_triangles[0]))
            assert(num.allclose(submesh['ghost_triangles'][1],true_ghost_triangles[1]))
            assert num.allclose(submesh['ghost_triangles'][2],true_ghost_triangles[2]), ParallelException('X')

            

            true_full_commun = [{0: [1, 2], 1: [1, 2], 2: [1, 2], 3: [2], 4: [1]}, {5: [0, 2], 6: [0, 2], 7: [], 8: [0], 9: [0], 10: [0]}, {11: [0, 1], 12: [0, 1], 13: [0], 14: [0], 15: [0]}]

            assert(true_full_commun == submesh['full_commun'])


            true_ghost_commun = [num.array([[ 5,  1],
                                            [ 6,  1],
                                            [ 8,  1],
                                            [ 9,  1],
                                            [10,  1],
                                            [11,  2],
                                            [12,  2],
                                            [13,  2],
                                            [14,  2],
                                            [15,  2]]),
                                 num.array([[ 0,  0],
                                            [ 1,  0],
                                            [ 2,  0],
                                            [ 4,  0],
                                            [11,  2],
                                            [12,  2]]),
                                 num.array([[0, 0],
                                            [1, 0],
                                            [2, 0],
                                            [3, 0],
                                            [5, 1],
                                            [6, 1]])]

            assert(num.allclose(submesh['ghost_commun'][0],true_ghost_commun[0]))
            assert(num.allclose(submesh['ghost_commun'][1],true_ghost_commun[1]))
            assert(num.allclose(submesh['ghost_commun'][2],true_ghost_commun[2]))



        barrier()
        #--------------------------------
        # Now do the comunnication part
        #--------------------------------


        if myid == 0:
            #----------------------------------------------------------------------------------
            # Test send_submesh
            #----------------------------------------------------------------------------------
            for p in range(1, numprocs):
                send_submesh(submesh, triangles_per_proc, p, verbose=False)

            #----------------------------------------------------------------------------------
            # Test extract_hostmesh
            #----------------------------------------------------------------------------------
            points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict  =\
            extract_hostmesh(submesh, triangles_per_proc)


            true_points =  [[0.5, 0.0], [0.5, 0.5], [0.5, 1.0], [1.0, 0.5], [1.0, 1.0], [0.25, 0.25], [0.75, 0.75], [0.0, 0.0], [0.0, 0.5], [0.0, 1.0], [1.0, 0.0], [0.25, 0.75], [0.75, 0.25]]

            true_vertices = [[1, 5, 0], [1, 6, 2], [3, 6, 1], [4, 6, 3], [2, 6, 4], [7, 5, 8], [8, 5, 1], [1, 11, 8], [2, 11, 1], [9, 11, 2], [0, 5, 7], [0, 12, 1], [10, 12, 0], [3, 12, 10], [1, 12, 3]]


            true_ghost_recv = {1: [num.array([5, 6, 7, 8, 9]), num.array([ 5,  6,  8,  9, 10])], 2: [num.array([10, 11, 12, 13, 14]), num.array([11, 12, 13, 14, 15])]}


            true_full_send = {1: [num.array([0, 1, 2, 4]), num.array([0, 1, 2, 4])], 2: [num.array([0, 1, 2, 3]), num.array([0, 1, 2, 3])]}

            assert(num.allclose(points,   true_points))
            assert(num.allclose(vertices, true_vertices))
            assert(num.allclose(ghost_recv_dict[1],true_ghost_recv[1]))
            assert(num.allclose(ghost_recv_dict[2],true_ghost_recv[2]))
            assert(num.allclose(full_send_dict[1],true_full_send[1]))
            assert(num.allclose(full_send_dict[2],true_full_send[2]))

            #print triangles_per_proc

        else:
            #----------------------------------------------------------------------------------
            # Test rec_submesh
            #----------------------------------------------------------------------------------
            points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict, no_full_nodes, no_full_trigs = rec_submesh(0, verbose=False)    

            if myid == 1:


                true_points =  [[0.0, 0.0], [0.0, 0.5], [0.0, 1.0], [0.5, 0.5], [0.5, 1.0], [0.25, 0.25], [0.25, 0.75], [0.5, 0.0], [1.0, 0.5], [1.0, 1.0], [0.75, 0.25], [0.75, 0.75]] 

                true_vertices =  [[0, 5, 1], [1, 5, 3], [1, 6, 2], [3, 6, 1], [4, 6, 3], [2, 6, 4], [3, 5, 7], [3, 11, 4], [8, 11, 3], [4, 11, 9], [7, 5, 0], [7, 10, 3]]

                true_ghost_recv =  {0: [num.array([6, 7, 8, 9]), num.array([0, 1, 2, 4])], 2: [num.array([10, 11]), num.array([11, 12])]}

                true_full_send =  {0: [num.array([0, 1, 3, 4, 5]), num.array([ 5,  6,  8,  9, 10])], 2: [num.array([0, 1]), num.array([5, 6])]}

                assert(num.allclose(points,   true_points))
                assert(num.allclose(vertices, true_vertices))
                assert(num.allclose(ghost_recv_dict[0],true_ghost_recv[0]))
                assert(num.allclose(ghost_recv_dict[2],true_ghost_recv[2]))
                assert(num.allclose(full_send_dict[0],true_full_send[0]))
                assert(num.allclose(full_send_dict[2],true_full_send[2])) 


            if myid == 2:

                true_points =   [[0.0, 0.0], [0.5, 0.0], [0.5, 0.5], [1.0, 0.0], [1.0, 0.5], [0.25, 0.25], [0.75, 0.25], [0.0, 0.5], [0.5, 1.0], [1.0, 1.0], [0.75, 0.75]]

                true_vertices =  [[1, 5, 0], [1, 6, 2], [3, 6, 1], [4, 6, 3], [2, 6, 4], [2, 5, 1], [2, 10, 8], [4, 10, 2], [9, 10, 4], [0, 5, 7], [7, 5, 2]]


                true_ghost_recv =   {0: [num.array([5, 6, 7, 8]), num.array([0, 1, 2, 3])], 1: [num.array([ 9, 10]), num.array([5, 6])]}

                true_full_send =   {0: [num.array([0, 1, 2, 3, 4]), num.array([11, 12, 13, 14, 15])], 1: [num.array([0, 1]), num.array([11, 12])]}


                assert(num.allclose(points,   true_points))
                assert(num.allclose(vertices, true_vertices))
                assert(num.allclose(ghost_recv_dict[0],true_ghost_recv[0]))
                assert(num.allclose(ghost_recv_dict[1],true_ghost_recv[1]))
                assert(num.allclose(full_send_dict[0],true_full_send[0]))
                assert(num.allclose(full_send_dict[1],true_full_send[1])) 

#-------------------------------------------------------------

if __name__=="__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_distribute_mesh, 'test')
    runner.run(suite)



        

        