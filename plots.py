import numpy as np
from hydrogen_wave_func import hydrogen_wave_func
from mayavi import mlab

print "starting.."
for n in [4]:
    for l in range(0, n):
        for m in range(0, l+1):
            print "n = %d, l = %d, m = %d" % (n, l, m)
            x,y,z,density = hydrogen_wave_func(n, l, m,10,20,20,20)
            print x,y,z,density
            density.dump('density.dat')
            density = np.load('density.dat')

            figure = mlab.figure('DensityPlot')
            mag = density / np.amax(density)
            pts = mlab.points3d(mag, opacity=0.5, transparent=True)
            mlab.colorbar(orientation='vertical')
            mlab.axes()
            mlab.savefig("images/%s_%s_%s.png" % (n,l,m))
            mlab.close()