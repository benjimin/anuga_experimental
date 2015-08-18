import numpy as np
import matplotlib.pyplot as plt


filename = 'urban.sww'


"""
Load prepared data
"""


from anuga import plot_utils as util
v = util.get_output(filename, timeSlices='last') # vertex data
c = util.get_centroids(filename, timeSlices='last') # centroid data
# Note, x/y/elev are already squeezed, unlike stage/height.



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


plots(elevations,triangulation,closeup)
