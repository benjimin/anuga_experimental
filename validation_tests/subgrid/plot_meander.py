import numpy as np
import matplotlib.pyplot as plt


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


plt.tripcolor(vx,vy,v.vols,v.height[-1],shading='gouraud',alpha=0.9,vmax=4)
plt.colorbar()

plt.quiver(cx,cy,c.xvel[2],c.yvel[2],alpha=0.4)

ax.axis('equal')



a = 640 # arm length
from math import pi
def distance(x,y): # so-called "dimensionless distance"
  return np.where(y>0, np.arctan2(y,x)/pi, np.where(x>0,y/a,1-y/a))+1

"""

for i in [0,-1]:
  fig = plt.figure()
  X = distance(cx,cy)  [c.height[i]!=0]
  Y = c.stage[i]       [c.height[i]!=0]
  plt.scatter(X,Y,alpha=0.4)
  
  
fig = plt.figure()  
for i in [0,-1]:  
  X = distance(cx,cy)  [c.height[i]!=0]
  Y = c.stage[i]       [c.height[i]!=0]
  plt.scatter(X,Y,alpha=0.4)

"""
plt.show()

#-------------------------
#"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np



#--------------------------
fig,ax = plt.subplots()

sc = None
def setup():
  global sc
  print "get ready!"
  i = 0
  X = distance(cx,cy)  [c.height[i]!=0]
  Y = c.stage[i]       [c.height[i]!=0]
  sc = ax.scatter(X,Y,animated=True)
  ax.axis('tight')
  return sc,
def update(i):
  global sc
  print "oooh",i
  X = distance(cx,cy)  [c.height[i]!=0]
  Y = c.stage[i]       [c.height[i]!=0]
  sc = ax.scatter(X,Y,animated=True)
  return sc,

import matplotlib.animation as ani
reference_handle = ani.FuncAnimation(fig,update,interval=5,init_func=setup,blit=True,repeat=10)
plt.show()


#"""
