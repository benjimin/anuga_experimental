path = "~/short/jobs/urban_centre/all/"

names = ['urban_DE1_160.sww',
'urban_DE1_300.sww',
'urban_DE1_500.sww',
'urban_DE1_750.sww',
'urban_DE1_80.sww',
'urban_DE1_80_10.sww',
'urban_DE1_80_20.sww',
'urban_DE1_80_22.sww',
'urban_DE1_950.sww',
'urban_DE1b_160.sww',
'urban_DE1b_300.sww',
'urban_DE1b_500.sww',
'urban_DE1b_750.sww',
'urban_DE1b_80.sww',
'urban_DE1b_80_10.sww',
'urban_DE1b_80_20.sww',
'urban_DE1b_80_22.sww',
'urban_DE1b_950.sww',
'urban_SG_160.sww',
'urban_SG_300.sww',
'urban_SG_500.sww',
'urban_SG_750.sww',
'urban_SG_80.sww',
'urban_SG_80_10.sww',
'urban_SG_80_20.sww',
'urban_SG_80_22.sww',
'urban_SG_950.sww']


fns = [path+n for n in names]


import numpy as np
import matplotlib.pyplot as plt
from math import pi
from anuga import plot_utils as util

"""
Utilities
"""

def smooth_centroids_to_vertices(centroid_values,triangles):
  nverts = np.asarray(triangles).max()+1 # get potential number of vertices
  triangles_at_vertices = [[] for element in range(nverts)] # declare list of lists
  for tri,vertices in enumerate(triangles):
    for vert in vertices: 
      triangles_at_vertices[vert].append(tri)  # populate
  return [np.mean(centroid_values[at_vertex]) for at_vertex in triangles_at_vertices] # apply

def filter_duplicate_vertices(vertex_data):
  raw_verts = zip(vertex_data.x, vertex_data.y) # list all vertex points
  uniq_verts = list(set(raw_verts)) # filter out duplicates
  map_pts_to_uniq_ind = dict((pt,i) for (i,pt) in enumerate(uniq_verts)) # map points into new index
  convert = [map_pts_to_uniq_ind[pt] for pt in raw_verts] # map old indices to new indices
  triplets = [[convert[i] for i in tri] for tri in vertex_data.vols] # remove vertex-duplicates from triangles list
  return uniq_verts,triplets

def interpolate(vertex_data, centroid_values=None, vertex_values=None, method='linear'): # stepwise, linear, or cubic interpolation
  uniq_verts,triplets = filter_duplicate_vertices(vertex_data)
  from matplotlib import tri as mtri
  T = mtri.Triangulation(*zip(*uniq_verts),triangles=triplets) # utilise clean triangulation
  if method=='stepwise': # step-wise 
    return lambda x,y,c=np.asarray(centroid_values),f=mtri.TrapezoidMapTriFinder(T): c[f(x,y)]
  else: # interpolation 
    if vertex_values is None: vertex_values = smooth_centroids_to_vertices(centroid_values,triplets)
    f = mtri.CubicTriInterpolator if method=='cubic' else mtri.LinearTriInterpolator
    return lambda x,y,f=f(T,vertex_values): f(x,y)
    
def transect(ax,filename,npoints=10**5,method='stepwise'):
  """
  This looks at values along a circle enclosing the town
  
  The town is centred at (750,250), and consists of a 5x5 array
  of 10x20 buildings placed with 30m periodicity. 
  
  This implies that the smallest bounding rectangle for the
  town has an area of ((5-1)*30+10)x((5-1)*30+20).
  
  The radius of the smallest enclosing circle is given by the 
  distance from the town centroid to a corner of the enclosing
  rectangle.
  
  Let the circle have padding akin to the minimum gap between
  the buildings, which is 30-max(10,20).
  """
  v = util.get_output(filename, timeSlices='last') # vertex data
  c = util.get_centroids(filename, timeSlices='last') # centroid data  
  
  padding = 10
  w,h = (5-1)*30+10, (5-1)*30+20
  radius = 0.5*(w**2+h**2)**0.5 + padding
  interp = interpolate(v,centroid_values=c.elev,method=method) 
 
  theta = np.linspace(0,2*pi,npoints)
  X = radius*np.cos(theta)+750
  Y = radius*np.sin(theta)+250
  Z = interp(X,Y)
  
  return theta,Z    
  
"""
Orchestrate plotting
"""


fig,axes = plt.subplots(nrows=3,ncols=9)


#for ax,f in zip([i for ii in axes for i in ii], fns):
ax = axes[1]
f = fns[1]
x,y = transect(ax,f)
ax.plot(x,y)
  
  
plt.show()
