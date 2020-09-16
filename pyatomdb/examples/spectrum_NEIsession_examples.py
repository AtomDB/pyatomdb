import pyatomdb, numpy, pylab


# declare the Non-Equilibrium Ionization session
sess = pyatomdb.spectrum.NEISession()

# set the response (here, a dummy response)
ebins = numpy.linspace(0.6,12.4,1000)
sess.set_response(ebins, raw=True)

# alternatively, could set to real RMF, ARF files using:
# sess.set_response(rmf, arf=arf)
# where rmf and arf are the file names

# Turn on (or off) thermal broadening, and if desired velocity broadening
# in this case turn on thermal broadening, and set velocity broadening to 500 km/s

sess.set_broadening(True, velocity_broadening=500.0)

# return spectrum at given temperature
kT = 0.8 # temperature in keV
tau = 1e11 # n_e*t cm^-3 s
kT_init=0.01 # initial temperature (for initial ion distribution)

spec = sess.return_spectrum(kT, tau, Te_init=kT_init)

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
spec = sess.return_spectrum(kT, tau, Te_init=kT_init)
ax.plot(sess.ebins_out, numpy.append(0, spec), drawstyle='steps', label='HEG Ioniz')

# same again, but for a recombining spectrum - set initial temperature > kT
kT_init=10.0
spec = sess.return_spectrum(kT, tau, Te_init=kT_init)

# plot the results
ax.plot(sess.ebins_out, numpy.append(0, spec), drawstyle='steps', label='HEG Recomb')
ax.legend(loc=0)
pylab.draw()

zzz=input("Press enter to exit")




