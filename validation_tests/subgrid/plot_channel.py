from anuga import plot_utils as util
from matplotlib import pyplot
import numpy
import run_anuga_subG_channel


# Compute analytical solution, which follows from the
# constant friction slope of steady, uniform flow

topography = run_anuga_subG_channel.topography
Qin = run_anuga_subG_channel.Qin

floodplain_length = run_anuga_subG_channel.floodplain_length
floodplain_width = run_anuga_subG_channel.floodplain_width

channel_width = run_anuga_subG_channel.chan_width

slp = run_anuga_subG_channel.floodplain_slope 
man_n = run_anuga_subG_channel.man_n 

k = (slp*(1./man_n)**2)**0.5 # At any point, the analytical solution says U = k*d^(2/3)

y_loc = floodplain_length*0.4 # Compare with the analytical solution here (avoiding boundary effects)

# Convenient vectors
xs = numpy.linspace(0.0, floodplain_width, num=600)
ys = numpy.linspace(0.0, floodplain_length, num=1200)

# Function to calculate the discharge, given the channel centre depth dc, assuming
# steady uniform flow
def discharge_su(stage):
    # Simple numerical integration
    depths = stage - topography(xs, y_loc + xs*0.)
    depths = depths * (depths > 0.)

    # Discharge = integral (U * depth) over the cross-section
    # U = k * depth **(2/3)
    discharge = k * (xs[1]-xs[0]) * (depths**(5.0/3.0) ).sum()
    return discharge

# Function that will be minimized to find the depth associated with discharge
# Qin
def minme(stage):
    q1 = discharge_su(stage)
    return (q1-Qin)**2.0

# Minimise the function mimne, to find the centre depth.
import scipy.optimize
stage_analytical = scipy.optimize.fmin(minme, x0=0.0)[0]

print 'stage_analytical ',stage_analytical



#########################################
# Plot the results
#########################################

p = util.get_centroids('drainage.sww')

# Get central channel indices
keep = (abs(p.x - floodplain_width/2.0) < 0.5*(channel_width/2.0)).nonzero()
#keep = range(len(p.x))

pyplot.ion()

pyplot.scatter(p.y[keep], p.elev[keep], color='brown')
pyplot.scatter(p.y[keep], p.stage[-1,keep], color='red')
pyplot.plot([y_loc, y_loc], [p.elev.min(), p.elev.max()], '-')
pyplot.plot([0.0, floodplain_length], [stage_analytical, stage_analytical], '-')

central_depth_analytical = stage_analytical - topography(floodplain_width/2.0, y_loc)
pyplot.plot(ys, topography(floodplain_width/2.0, ys) + central_depth_analytical, '-', color='black')
