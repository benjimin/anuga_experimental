from anuga import plot_utils as util
from matplotlib import pyplot
import numpy

p = util.get_centroids('shallow_slope.sww')

pyplot.ion()

pyplot.scatter(p.y,p.stage[p.stage.shape[0]-1,:]-p.elev)
pyplot.plot(numpy.array([0., 1000.]), numpy.array([0.01, 0.01]),'-')

