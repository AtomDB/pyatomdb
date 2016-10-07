import pyatomdb, numpy, pylab

# set up a grid of energy bins to model the spectrum on:
ebins=numpy.linspace(0.3,10,1000)

# define a broadening, in keV, for the lines
de = 0.01

# define the temperature at which to plot (keV)
te = 3.0

# find the index which is closest to this temperature
ite = pyatomdb.spectrum.get_index( te, teunits='keV', logscale=False)

# create both a broadened and an unbroadened spectrum
a = pyatomdb.spectrum.make_spectrum(ebins, ite,dummyfirst=True)
b = pyatomdb.spectrum.make_spectrum(ebins, ite, broadening=de, \
                                    broadenunits='kev',dummyfirst=True)
# The dummyfirst argument adds an extra 0 at teh beginning of the
# returned array so it is the same length as ebins. It allows
# accurate plotting using the "drawstyle='steps'" flag to plot.

# plot the results
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

ax.loglog(ebins, a, drawstyle='steps', label='Unbroadened')
ax.loglog(ebins, b, drawstyle='steps', label='sigma = %.2f'%(de))
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Emissivity (ph cm$^{3}$ s$^{-1}$ bin$^{-1}$)')
ax.legend(loc=0)
pylab.draw()
zzz = raw_input("Press enter to continue")

print "Listing lines between 1 and 2 A"
# now list the lines in a wavelength region
llist = pyatomdb.spectrum.list_lines([1,2.0], index=ite)
# print these to screen
pyatomdb.spectrum.print_lines(llist)
# print to screen, listing the energy, not the wavelength
print "Listing lines between 1 and 2 A, using keV."

pyatomdb.spectrum.print_lines(llist, specunits = 'keV')

