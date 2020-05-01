import pyatomdb
import numpy
import pylab


# declare the Collisional Ionization Equilibrium session
sess = pyatomdb.spectrum.CIESession()

# create a set of energy bins (in keV) for the response. Note these are
# the n edges of the n-1 bins.
ebins = numpy.linspace(0.3,1.0,10000)

# set the response (raw keyword tells pyatomdb it is not a real response file)
sess.set_response(ebins, raw=True)

kT = 0.4 # temperature in keV
spec = sess.return_spectrum(kT)

# we now have a spectrum and the energy bins,

# prepend a zero to the spectrum so that energy bins and spectrum have
# the same length. If you use a different plotting system you may
# need to add this to the end.
spec = numpy.append(0, spec)

# Returned spectrum has units of photons cm^5 s^-1 bin^-1
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(211)
ax.plot(sess.ebins_out, spec, drawstyle='steps', label='dummy response')

# label the figure
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')

# zoom in a bit
#ax.set_xlim([0.8,0.85])

# make plot appear

# now repeat the process with a real response
# set response
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

# get spectrum, prepend zero
spec = sess.return_spectrum(kT)
spec = numpy.append(0, spec)

ax2 = fig.add_subplot(212)
ax2.plot(sess.ebins_out, spec, drawstyle='steps', label='Chandra MEG')

# add legends
ax.legend(loc=0)
ax2.legend(loc=0)

# zoom in on small sections of the spectrum
ax.set_xlim([0.7, 0.9])
ax2.set_xlim([0.7, 0.9])

# set axes
ax2.set_xlabel('Energy (keV)')
ax2.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')

# adjust plot spacing so labels are visible
pylab.matplotlib.pyplot.subplots_adjust(hspace = 0.34)

# draw graphs
pylab.draw()

zzz=input("Press enter to continue")

# save image files
fig.savefig('spectrum_session_examples_1_1.pdf')
fig.savefig('spectrum_session_examples_1_1.svg')

