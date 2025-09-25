import pyatomdb, numpy
import pylab # for demo plot

# make a kappa model
a=pyatomdb.spectrum.KappaSession()

# set the response (here a dummy, can use rmf and arf filenames too)
ebins = numpy.linspace(0.01,5,1001)
a.set_response(ebins, raw=True)

# define the kappa parameters
kp = 5.0 # kappa parameter, >2
T = 2e5  # Temperature, K

# fetch the spectrum, in photons bin-1 cm^3 s^-1
spectrum = a.return_spectrum(T, kp, teunit='K')

# Additional options:
# turn on and off line emission, continuum emission and pseudocontinuum
# using the dolines, docont and dopseudo properties
# (all three are True by default)

a.dolines=False
spectrum2 = a.return_spectrum(T, kp, teunit='K')

# test plot
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

ax.loglog(ebins, numpy.append(spectrum[0], spectrum), drawstyle='steps', label='lines')
ax.loglog(ebins, numpy.append(spectrum2[0],spectrum2), drawstyle='steps', label='nolines')

ax.legend(loc=0)
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')
ax.set_title("T=%.1e K, kappa=%.1f"%(T, kp))

pylab.draw()

zzz = input("Press enter to continue")

pylab.savefig("spectrum_kappa_examples_1.pdf")
pylab.savefig("spectrum_kappa_examples_1.svg")
