import pyatomdb
import numpy
import pylab


# declare the Collisional Ionization Equilibrium session
sess = pyatomdb.spectrum.CIESession()

# create a set of energy bins (in keV) for the response. Note these are
# the n edges of the n-1 bins.
ebins = numpy.linspace(0.6,12.4,1000)

# set the response (raw keyword tells pyatomdb it is not a real response file)
sess.set_response(ebins, raw=True)

kT = 0.4 # temperature in keV
spec = sess.return_spectrum(kT)


# alternatively, could set to real RMF, ARF files using:
# sess.set_response(rmf, arf=arf)
# where rmf and arf are the file names

# Turn on (or off) thermal broadening, and if desired velocity broadening
# in this case turn on thermal broadening, and set velocity broadening to 500 km/s

sess.set_broadening(True, velocity_broadening=500.0)

# return spectrum at given temperature
kT = 0.8 # temperature in keV
spec = sess.return_spectrum(kT)

# note returned spectrum has units of photons cm^5 s^-1 bin^-1, and  has 1 less
# value than the energy bin grid. Prepend a 0 to it for plotting purposes
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)
ax.plot(sess.ebins_out, numpy.append(0, spec), drawstyle='steps', label='dummy response')

#label the figure
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^3$ s$^{-1}$ bin$^{-1}$)')

# now do the same thing but with an instrument response
sess.set_response('../tests/testdata/aciss_heg1_cy19.grmf', \
                   arf='../tests/testdata/aciss_heg1_cy19.garf')

# return the spectrum
spec = sess.return_spectrum(kT)
ax.plot(sess.ebins_out, numpy.append(0, spec), drawstyle='steps', label='HEG')

# change abundance set
sess.set_abundset('Feldman')
spec = sess.return_spectrum(kT)
ax.plot(sess.ebins_out, numpy.append(0, spec), drawstyle='steps', label='HEG+Feldman')
ax.legend(loc=0)
pylab.draw()

zzz=input("Press enter to exit")




