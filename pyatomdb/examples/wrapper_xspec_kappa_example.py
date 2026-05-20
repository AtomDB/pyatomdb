# temporarily add the wrappers folder to the path
import sys
sys.path.append('../wrappers')

import matplotlib.pyplot as plt

from kappa_xspec import *

# declare a new model
# inital import creates pykappa, pyvkappa, pyvvkappa, analagous to apec, vapec, vvapec

m = xspec.Model('pykappa')

m.show()

m.pykappa.kT=3.0

m.pykappa.kappa = 2.5

# let's plot a spectrum
xspec.Plot.device='/xs'


#set dummy response
xspec.AllData.dummyrsp(0.1, 10, 10001,'lin')

xspec.Plot('model')

x1 = xspec.Plot.x(1)
m1 = xspec.Plot.model(1)


m.pykappa.Velocity = 1000.0


xspec.Plot('model')

x2 = xspec.Plot.x(1)
m2 = xspec.Plot.model(1)
fig = plt.figure()
fig.show()
ax = fig.add_subplot(111)
ax.plot(x1,m1, label='0km/s')
ax.plot(x2,m2, label='1000km/s')
ax.legend(loc=0)
plt.draw()

zzz=input('Press enter')
