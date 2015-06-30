import numpy

from clawpack.clawutil.data import ClawData
probdata = ClawData()
probdata.read('setprob.data',force=True)

def mapc2p(xc,yc):
    
    # Note this will be removed once bug fix is implemented in slices
    xp = xc
    yp = yc

    return xp, yp
