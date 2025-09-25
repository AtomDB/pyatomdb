import pyatomdb, numpy
import pylab # for demo plot

# make a kappa model
a=pyatomdb.spectrum.KappaSession()

# set the response (here a dummy, can use rmf and arf filenames too)
ebins = numpy.linspace(0.1,5,1001)
a.set_response(ebins, raw=True)

# define the kappa parameters
kp = 5.0 # kappa parameter, >2
T = 1e6  # Temperature, K

# fetch the spectrum, in photons bin-1 cm^3 s^-1
spectrum = a.return_spectrum(T, kp, teunit='K')
spectrum = numpy.append(spectrum[0],spectrum)

kp_hi = 50
# make a maxwellian equivalent (kp -> infinity)
spectrum_kmaxw = a.return_spectrum(T, kp_hi, teunit='K')
spectrum_kmaxw = numpy.append(spectrum_kmaxw[0],spectrum_kmaxw)

# now make a true maxwellian plasma using a CIE session
b = pyatomdb.spectrum.CIESession()
b.set_response(ebins, raw=True)
spectrum_maxw = b.return_spectrum(T, teunit='K')
spectrum_maxw = numpy.append(spectrum_maxw[0],spectrum_maxw)


# test plot
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

ax.loglog(ebins, spectrum, drawstyle='steps', label='$\kappa =%.1f$'%(kp))
ax.loglog(ebins, spectrum_kmaxw, drawstyle='steps', label='$\kappa = %.1f$'%(kp_hi))
ax.loglog(ebins, spectrum_maxw, drawstyle='steps', label='Maxw')

ax.legend(loc=0)
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')
ax.set_title("T=%.1e K"%(T))
ax.set_ylim(ymin=1e-22)

pylab.draw()

zzz = input("Press enter to continue")

pylab.savefig("spectrum_kappa_examples_2.pdf")
pylab.savefig("spectrum_kappa_examples_2.svg")
