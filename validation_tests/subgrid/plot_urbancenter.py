import numpy as np
import matplotlib.pyplot as plt


filename = 'urban_DE1_800_50.sww'


"""
Load prepared data
"""


from anuga import plot_utils as util
v = util.get_output(filename, timeSlices='last') # vertex data
c = util.get_centroids(filename, timeSlices='last') # centroid data
# Note, x/y/elev are already squeezed, unlike stage/height.

print "Simulated (final) time is {0} seconds".format(c.time[-1])

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



"""
Design specific plots
"""



def elevations(ax):
  plt.title('Elevation')
  plt.scatter(v.x,v.y)
  plt.tripcolor(v.x,v.y,v.vols,v.elev,shading='gouraud',alpha=0.8)
  
def triangulation(ax,**opts):
  plt.title('Final water depth')
  plt.tripcolor(v.x,v.y,v.vols,v.stage[0],shading='gouraud',alpha=0.8,**opts)
  plt.colorbar()
  plt.quiver(c.x,c.y,c.xvel,c.yvel,alpha=0.5)
  
def closeup(ax):
  triangulation(ax,cmap='gist_rainbow_r',vmax=1.15,vmin=0.65)
  ax.set_xlim(550,950)
  ax.set_ylim(75,425)
  
def transect(ax):
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
  padding = 10
  w,h = (5-1)*30+10, (5-1)*30+20
  radius = 0.5*(w**2+h**2)**0.5 + padding
  npoints = 10**5
  interp = interpolate(v,centroid_values=c.stage,method='stepwise') # <---
  from math import pi
  theta = np.linspace(0,2*pi,npoints)
  X = radius*np.cos(theta)+750
  Y = radius*np.sin(theta)+250
  Z = interp(X,Y)
  
  elevations(ax) # do basemap
  plt.plot(X,Y)
  plt.title('transect')
  
  fig=plt.figure()
  plt.plot(theta,Z)
  plt.xlabel('bearing')
  plt.ylabel('stage')
  
  


"""
Orchestrate plots
"""

def plot():
  fig = plt.figure()
  ax = fig.add_subplot(111,aspect='equal')
  limits = lambda x: (min(x),max(x))
  ax.set_xlim(*limits(v.x))
  ax.set_ylim(*limits(v.y))
  return ax 
  
def plots(*args):
  for func in args:
     func(plot())
  plt.show()


#plots(elevations,triangulation,closeup,transect)
plots(transect)

