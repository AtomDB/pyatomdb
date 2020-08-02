from apec_xspec import *

# declare a new model
# inital import creates pyapec, pyvapec, pyvvapec, analagous to apec, vapec, vvapec

m1 = xspec.Model('pyapec')

m1.show()

m1.pyapec.kT=5.0


# let's plot a spectrum
xspec.Plot.device='/xs'

xspec.Plot('model')

# you can fiddle with some aspects of the model directly
cie.set_eebrems(False) # turn off electron-electron bremsstrahlung

xspec.Plot('model')

# it is also possible to go in and change broadening, etc. The interface can be
# adjust quite easily by looking at apec_xspec.py to add extra parameters



