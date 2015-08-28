import numpy as np
import matplotlib.pyplot as plt

filename = 'meander.sww'

from anuga import plot_utils as util
v = util.get_output(filename) # vertex data
c = util.get_centroids(filename) # centroid data


print v.x.shape

print v.elev.shape

print v.height.shape

print v.stage.shape

print dir(v)




fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')
limits = lambda x: (min(x),max(x))
ax.set_xlim(*limits(v.x))
ax.set_ylim(*limits(v.y))

#plt.tripcolor(v.x,v.y,v.vols,v.stage[0],shading='gouraud',alpha=0.85)
plt.tripcolor(v.x,v.y,v.vols,v.elev,shading='gouraud',alpha=0.85)
plt.colorbar()

plt.show()

