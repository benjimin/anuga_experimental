from anuga import plot_utils as util
from matplotlib import pyplot
import numpy

from run_anuga_subG_drainage import stage_rel_to_plain as h
from run_anuga_subG_drainage import filename,drain_depth

#p = util.get_centroids(filename+'.sww')

fig,axes = pyplot.subplots(nrows=2)

for ax,sww,label in zip(axes,['drainage-new.sww.de1','drainage-new.sww.sg'],['DE1','Subgrid']):
  #pyplot.ion()
  p = util.get_centroids(sww)

  ax.scatter(p.y,p.elev+0.1*p.y,color='brown',alpha=0.4)
  ax.scatter(p.y,p.stage[-1,:]+0.1*p.y,s=10,alpha=0.7)
  ax.scatter(p.y,p.height[-1,:],s=80,facecolor='none',edgecolor='b',alpha=0.35)

  #pyplot.plot(numpy.array([0., 1000.]), numpy.array([0.01, 0.01]),'-')
  ax.plot([0,1000], [h,h],'--')
  ax.plot([0,1000], [0,0],'--',color='black')
  ax.plot([0,1000], [drain_depth,drain_depth],'--',color='black')

  ax.set_xlim(150,550)
  ax.set_ylim(-0.6,0.8)
  ax.set_xlabel('Distance (m)')
  ax.set_ylabel(label+' water depth (m)')


pyplot.show()
