
from pylab import *
from clawpack.visclaw.data import ClawPlotData
from setrun import setrun
from scipy import interpolate
from numpy import sqrt, sin, exp
import os

close('all')

def savetxt_verbose(fname,d):
    print "Created ",fname
    savetxt(fname,d,fmt='%12.3e')

def savefig_verbose(fname):
    print "Created ",fname
    savefig(fname)

plotdata = ClawPlotData()
rundata = setrun()

plotdata.outdir = os.getenv('PWD') # set to the current output directory
ngauges = 516
trans_halfwidth = rundata.probdata.trans_halfwidth
weight_method = 1 # 0 for constant, 1 for gaussian, 2 for super-gaussian
g = plotdata.getgauge(0)

# Interpolate to equally spaced time points, 
tend = g.t[-5]
tstart = g.t[0]
dt = 25e-9   # should evenly divide 500e-9 (500 ns)
tt = arange(tstart, tend, dt)
tt_length = len(tt)
p_average = zeros((tt_length,))
u_average = zeros((tt_length,))

figure(1,figsize=(12,6))
axp = axes()
figure(2,figsize=(12,6))
axu = axes()

# compute average p and u values
for gaugeno in range(ngauges):
    g = plotdata.getgauge(gaugeno)
    y = g.location[0]*sin(g.location[1])
    if (weight_method == 0):
        weight = 1.0
    elif (weight_method == 1):
        weight = exp(-(y/trans_halfwidth)**2.0)
    elif (weight_method == 2):
        weight = exp(-(y/trans_halfwidth)**8.0)
    p = -(g.q[0] + g.q[1] + g.q[2])/3.0
    pp = interpolate.interp1d(g.t, p)
    p_average += weight*pp(tt)
    u = g.q[6,:]
    uu = interpolate.interp1d(g.t, u)
    u_average += weight*uu(tt)
    cy = -y/(trans_halfwidth + 1.0e-10) # to avoid negative machine precision values
    cz = g.location[2]/trans_halfwidth
    c = [(1-cy)*(1-cz), cz, cy]
    figure(1)
    axp.plot(g.t, p, color=c)
    figure(2)
    axu.plot(g.t, u, color=c)

figure(1)
axp.set_title("p at gauges on transducer (red-center, blue-y_edge, green-z_edge)")
#savefig_verbose('p_gauges.png')

figure(2)
axu.set_title("u at gauges on transducer (red-center, blue-y_edge, green-z_edge)")
savefig_verbose('u_gauges.png')


p_average = p_average / ngauges
u_average = u_average / ngauges

figure(figsize=(12,6))
plot(tt, p_average, 'r',label='b=1.25a')  # interpolated
title("Average pressure on transducer")
#savefig_verbose('p_average.png')

figure(figsize=(12,6))
plot(tt, u_average, 'r',label='b=1.25a')  # interpolated
title("Average u velocity on transducer")
savefig_verbose('u_average.png')

d = vstack((tt, p_average, u_average)).T
suffix = os.getenv('Z')
savetxt_verbose('t_p_u_average_' + suffix + '.txt', d)

show()
