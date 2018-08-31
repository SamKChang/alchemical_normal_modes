from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
img = read_png('mo_N2.png')
x, y = ogrid[0:img.shape[0], 0:img.shape[1]]
ax = gca(projection='3d')
ax.plot_surface(x, y, 10, rstride=5, cstride=5, facecolors=img)
show()
