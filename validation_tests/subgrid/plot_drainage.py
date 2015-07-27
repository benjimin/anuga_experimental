from anuga import plot_utils as util
from matplotlib import pyplot
import numpy

from run_anuga_subG_drainage import stage_rel_to_plain as h
from run_anuga_subG_drainage import filename,drain_depth

p = util.get_centroids(filename+'.sww')

pyplot.ion()

pyplot.scatter(p.y,p.elev+0.1*p.y,color='brown')

pyplot.scatter(p.y,p.stage[-1,:]+0.1*p.y)


#pyplot.plot(numpy.array([0., 1000.]), numpy.array([0.01, 0.01]),'-')
pyplot.plot([0,1000], [h,h],'--')
pyplot.plot([0,1000], [0,0],'--',color='black')
pyplot.plot([0,1000], [drain_depth,drain_depth],'--',color='black')
