import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

filename = 'meander.sww'

from anuga import plot_utils as util
v = util.get_output(filename) # vertex data
c = util.get_centroids(filename) # centroid data

vx,vy = v.x+v.xllcorner,v.y+v.yllcorner # transforms
cx,cy = c.x+c.xllcorner,c.y+c.yllcorner


"""fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')
limits = lambda x: (min(x),max(x))
ax.set_xlim(*limits(vx))
ax.set_ylim(*limits(vy))"""

fig,ax = plt.subplots()
i = -1 # final timestep
w = [c.height[i]!=0] # wetted subset of cells
plt.triplot(vx,vy,v.vols,alpha=0.1) # mesh background
#plt.tripcolor(vx,vy,v.vols[w],v.height[i],shading='gouraud',alpha=0.9,vmax=4)
plt.tripcolor(vx,vy,v.vols[w],v.stage[i],shading='gouraud',alpha=0.9,vmax=max(c.stage[i][w]),vmin=min(c.stage[i][w]))
plt.colorbar()
plt.quiver(cx[w],cy[w],c.xmom[i][w],c.ymom[i][w],alpha=0.35)
ax.axis('tight')
plt.show()

"""
fig,ax = plt.subplots()
limits = lambda x: (min(x),max(x))
ax.set_xlim(*limits(vx))
ax.set_ylim(*limits(vy))
def start():
  plt.triplot(vx,vy,v.vols,alpha=0.1) # mesh background 
  seq=update(0)
  #plt.colorbar(seq[0])
  #ax.axis('equal')
  return seq
def update(i):
  w = [c.height[i]!=0]
  tri=plt.tripcolor(vx,vy,v.vols[w],v.height[i],shading='gouraud',alpha=0.9,vmax=4,animated=True)
  #vec=plt.quiver(cx[w],cy[w],c.xvel[i][w],c.yvel[i][w],alpha=0.35,animated=True)
  vec = plt.quiver([],[],[],[],animated=True)
  return tri,vec
hdl = ani.FuncAnimation(fig,update,frames=len(v.timeSlices),
  init_func=start,blit=True,repeat=False,interval=500)
plt.show()
"""




# -- watch convergence of Stelling fig.13

a = 640 # arm length
from math import pi
def distance(x,y): # so-called "dimensionless distance"
  return np.where(y>0, np.arctan2(y,x)/pi, np.where(x>0,y/a,1-y/a))+1


fig,ax = plt.subplots()
def setup():
  seq = update(0)
  ax.axis('tight')
  return seq
def update(i):
  X = distance(cx,cy)  [c.height[i]!=0]
  Y = c.stage[i]       [c.height[i]!=0]
  plot = ax.scatter(X,Y,animated=True)
  return plot, # as sequence

reference_handle = ani.FuncAnimation(fig,update,frames=len(v.timeSlices),
  init_func=setup,blit=True,repeat=False,interval=50)
  
plt.show()


#"""
