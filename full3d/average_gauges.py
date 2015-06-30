
from pylab import *
from clawpack.visclaw.data import ClawPlotData
from scipy import interpolate
import os

close('all')

def savetxt_verbose(fname,d):
    print "Created ",fname
    savetxt(fname,d,fmt='%12.3e')

def savefig_verbose(fname):
    print "Created ",fname
    savefig(fname)

plotdata = ClawPlotData()

plotdata.outdir = os.getenv('PWD') # set to the current output directory
ngauges = 1961
g = plotdata.getgauge(1)
p_average = zeros(g.t.shape)
u_average = zeros(g.t.shape)

figure(1,figsize=(12,6))
axp = axes()
figure(2,figsize=(12,6))
axu = axes()

# compute average p and u values
for gaugeno in range(ngauges):
    g = plotdata.getgauge(gaugeno)
    p = -(g.q[0] + g.q[1] + g.q[2])/3.0
    p_average += p
    u = g.q[6,:]
    u_average += u
    cy = abs(gaugeno/float(ngauges))  # for one-sided gauges only above y=0
    c = [1-cy, 0, cy]  # for upper half domain
    figure(1)
    axp.plot(g.t, p, color=c)
    figure(2)
    axu.plot(g.t, u, color=c)

figure(1)
axp.set_title("p at gauges on transducer (red near center, blue near edge)")
#savefig_verbose('p_gauges.png')

figure(2)
axu.set_title("u at gauges on transducer (red near center, blue near edge)")
#savefig_verbose('u_gauges.png')


p_average = p_average / ngauges
u_average = u_average / ngauges

# Interpolate to equally spaced time points, 
tend = g.t[-1]
tstart = g.t[0]
dt = 25e-9   # should evenly divide 500e-9 (500 ns)
tt = arange(tstart, tend, dt)

pp = interpolate.interp1d(g.t, p_average)
pt = pp(tt)

uu = interpolate.interp1d(g.t, u_average)
ut = uu(tt)

figure(figsize=(12,6))
plot(g.t, p_average, 'b')  # original
plot(tt, pt, 'r',label='b=1.25a')  # interpolated
title("Average pressure on transducer")
#savefig_verbose('p_average.png')

figure(figsize=(12,6))
plot(g.t, u_average, 'b')  # original
plot(tt, ut, 'r',label='b=1.25a')  # interpolated
title("Average u velocity on transducer")
savefig_verbose('u_average.png')

d = vstack((tt, pt)).T
#savetxt_verbose('p_average.txt', d)

d = vstack((tt, ut)).T
#savetxt_verbose('u_average.txt', d)

d = vstack((tt, pt, ut)).T
suffix = os.getenv('Z')
savetxt_verbose('t_p_u_average_' + suffix + '.txt', d)

show()
