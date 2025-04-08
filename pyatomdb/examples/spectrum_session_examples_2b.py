import pyatomdb
import numpy
import pylab

# line broadening of weak features

# declare an CIESession

sess = pyatomdb.spectrum.CIESession()

# create a set of energy bins (in keV) for the response. Note these are
# the n edges of the n-1 bins.
ebins = numpy.linspace(6.9,7.0,501)

# set the response (raw keyword tells pyatomdb it is not a real response file)
sess.set_response(ebins, raw=True)

kT = 10.0 # temperature in keV
vel = 400.0 # velocity to broaden with when using velocity broadening (km/s)

spec = sess.return_spectrum(kT)

# we now have a spectrum and the energy bins,

# prepend a zero to the spectrum so that energy bins and spectrum have
# the same length. If you use a different plotting system you may
# need to add this to the end.
spec = numpy.append(0, spec)

# Returned spectrum has units of photons cm^5 s^-1 bin^-1
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)
ax.plot(sess.ebins_out, spec, drawstyle='steps', label='dummy response')

# label the figure
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')
ax.set_yscale('log')


sess.set_broadening(True, velocity_broadening=vel)
spec = sess.return_spectrum(kT)
spec = numpy.append(0, spec)
ax.plot(sess.ebins_out, spec, drawstyle='steps', label='%.0f km/s'%(vel))

# Some features from weak lines aren't broadened. Let's fix that
# by setting the minimum emissivity limit to a lower value 
# (default is 1e-18 ph cm^3 s^-1)

sess.set_broadening(True, velocity_broadening=vel, broaden_limit=1e-30)
spec = sess.return_spectrum(kT)
spec = numpy.append(0, spec)
ax.plot(sess.ebins_out, spec, drawstyle='steps', label='Low Limit, %.0f km/s'%(vel))

# make plot appear

# now repeat the process with a real response
# add legends
ax.legend(loc=0)

# zoom in on small sections of the spectrum
#ax.set_xlim([1.150, 1.180])

# draw graphs
pylab.draw()

zzz=input("Press enter to continue")

# save image files
fig.savefig('spectrum_session_examples_2b_1.pdf')
fig.savefig('spectrum_session_examples_2b_1.svg')

