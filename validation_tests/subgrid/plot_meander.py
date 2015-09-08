import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

filename = 'meander.sww'

from anuga import plot_utils as util
v = util.get_output(filename) # vertex data
c = util.get_centroids(filename) # centroid data

vx,vy = v.x+v.xllcorner,v.y+v.yllcorner # transforms
cx,cy = c.x+c.xllcorner,c.y+c.yllcorner



# --- basemap

fig,ax = plt.subplots()
i = -1 # final timestep
w = [c.height[i]!=0] # wetted subset of cells
plt.triplot(vx,vy,v.vols,alpha=0.1) # mesh background
#plt.tripcolor(vx,vy,v.vols,v.elev,shading='gouraud',alpha=0.9) # elevation
plt.tripcolor(vx,vy,v.vols[w],v.stage[i],shading='gouraud',alpha=0.9,vmax=max(c.stage[i][w]),vmin=min(c.stage[i][w])) # wet stage
plt.colorbar()
plt.quiver(cx[w],cy[w],c.xmom[i][w],c.ymom[i][w],alpha=0.3,color='white',scale_units='x',scale=0.1)
# scale units can relate to the data, or pixels, or window. The scale is inverse of length per data.
# ideally, I would want to automatically look at the distribution of (nonzero) speeds, and relate this to the mean triangle width..
# whereas existing autoscale seems to struggle with wide distributions..
ax.axis('equal')
plt.title('Final water stage (with momentum vectors)')
plt.show()

"""
def start():
  plt.triplot(vx,vy,v.vols,alpha=0.1) # mesh background 
  seq=update(0)
  plt.colorbar()#seq[0])
  #ax.axis('equal')
  return seq
def update(i):
  w = [c.height[i]!=0]
  #tri=plt.tripcolor(vx,vy,v.vols[w],v.height[i],shading='gouraud',alpha=0.9,vmax=4,animated=True)
  tri=plt.tripcolor(vx,vy,v.vols,froude[i],alpha=0.9,animated=True,vmin=np.nanmin(froude[i]),vmax=np.nanmax(froude[i]))
  #vec=plt.quiver(cx[w],cy[w],c.xvel[i][w],c.yvel[i][w],alpha=0.35,animated=True)
  #vec = plt.quiver([],[],[],[],animated=True)
  return tri,#vec
hdl = ani.FuncAnimation(fig,update,frames=len(v.timeSlices),
  init_func=start,blit=True,repeat=False,interval=500)
plt.show()
"""

# -- check flow rate (over time)



def get_area(tri): # get area of any triangle cell
  p1,p2,p3 = ((vx[i],vy[i]) for i in tri) # triangle vertices
  ds = lambda (x,y),(x2,y2): ((x-x2)**2+(y-y2)**2)**(0.5) # distance
  a,b,c = ds(p1,p2),ds(p2,p3),ds(p3,p1) # side lengths
  s = 0.5*(a+b+c) # semiperimeter
  return (s*(s-a)*(s-b)*(s-c))**(0.5) # Heron's formula
areas = np.array(map(get_area,v.vols)) # all cells

# areas = util.triangle_areas(v)

volume = c.height * areas

total = np.diff(volume.sum(axis=1)) # differentiate total volume
upper = np.diff(volume[:,cx>0].sum(axis=1))
lower = np.diff(volume[:,cx<0].sum(axis=1))
t = (c.time[1:]+c.time[:-1])/2. # mid-timesteps
fig,ax = plt.subplots()
for y,txt in [(-total,'total leakage rate'),(-upper,'flow from east'),(lower,'flow to west')]:
  ax.plot(t,y,label=txt,alpha=0.6)
plt.legend()
ax.set_ylabel('Flow rate (m^3/s)')
ax.set_xlabel('Time (s)')
plt.show()



# -- watch convergence of Stelling fig.13

a = 640 # arm length
from math import pi
def distance(x,y): # so-called "dimensionless distance"
  return np.where(y>0, np.arctan2(y,x)/pi, np.where(x>0,y/a,1-y/a))+1


fig,ax = plt.subplots()
def setup():
  seq = update(0)
  ax.set_xlabel('Dimensionless distance')
  ax.set_ylabel('Water stage level: m')
  ax.set_xlim(-1.5,4.5) # don't lose focus onto the furthest reaches of the ponds
  return seq
def update(i):
  X = distance(cx,cy)  [c.height[i]!=0]
  Y = c.stage[i]       [c.height[i]!=0]
  plot = ax.scatter(X,Y,alpha=0.3,animated=True)
  return plot, # as sequence
reference_handle = ani.FuncAnimation(fig,update,frames=len(v.timeSlices),
  init_func=setup,blit=True,repeat=False,interval=50)
plt.show()


#"""
