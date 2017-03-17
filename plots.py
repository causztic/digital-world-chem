import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from hydrogen_wave_func import hydrogen_wave_func

x,y,z,mag=hydrogen_wave_func(3,1,0,10,20,20,20)
x.dump('xdata310.dat')
y.dump('ydata310.dat')
z.dump('zdata310.dat')
mag.dump('density310.dat')

x = np.load('xdata310.dat')
y = np.load('ydata310.dat')
z = np.load('zdata310.dat')
mag = np.load('density310.dat')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for a in range(0,len(mag)):
    for b in range(0,len(mag)):
        for c in range(0,len(mag)):
            ax.scatter(x[a][b][c],y[a][b][c],z[a][b][c], marker='o', alpha=(mag[a][b][c]/np.amax(mag)))

plt.show()