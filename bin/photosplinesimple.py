import numpy as np
from photospline import glam_fit, ndsparse, bspline

numpts = 500
x1 = np.logspace(-1,5,numpts)
z = 1e-30 * x1**-1 * np.array(x1>1e2,dtype=float)
z = np.where(z<=0,1e-50,z)
print(z)

w = np.ones_like(z)

order = 2
penalty_order = order
knots = [np.logspace(-2,6,65)]
smooth = 1#3.14159e3

zs, w = ndsparse.from_data(z, w)
spline = glam_fit(zs,w,[x1],knots,[order],[smooth],[penalty_order])

import pylab
# Plot individual basis splines 
xfine = np.logspace(np.log10(knots[0][0]), np.log10(knots[0][-1]), 10001)
# splines = [np.array([bspline(knots[0], x, n, order) for x in xfine]) for n in range(0,len(knots[0])-2-1)]
# for n in range(len(splines)):
#     pylab.plot(xfine, spline.coefficients[n]*splines[n])
# Plot the spline surface (sum of all the basis functions)
pylab.plot(xfine, spline.grideval([xfine]), label='Spline fit', color='k')
pylab.scatter(x1, z, label='Data')
pylab.legend(loc='upper left')
pylab.loglog()
pylab.savefig("photosimple.png")