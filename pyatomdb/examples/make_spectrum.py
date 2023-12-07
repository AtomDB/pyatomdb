import pyatomdb, numpy, pylab

# set up a grid of energy bins to model the spectrum on:
ebins=numpy.linspace(0.3,10,10000)

# define a broadening, in keV, for the lines
de = 0.01

# define the temperature at which to plot (keV)
te = 3.0

# declare a CIESession
CIE = pyatomdb.spectrum.CIESession()

# set the response. Raw keyword makes it a diagonal response
# (no line spread or effective area)

CIE.set_response(ebins, raw=True)

# create both a broadened and an unbroadened spectrum
a = CIE.return_spectrum(te)

# turn on broadening
CIE.set_broadening(True)
b = CIE.return_spectrum(te)

# plot the results
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

ax.loglog(ebins, numpy.append(a[0], a), drawstyle='steps', label='Unbroadened')
ax.loglog(ebins, numpy.append(b[0],b), drawstyle='steps', label='Thermal')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Emissivity (ph cm$^{3}$ s$^{-1}$ bin$^{-1}$)')
ax.legend(loc=0)
pylab.draw()
zzz = input("Press enter to continue")

print("Listing lines between 1 and 2 A")
# now list the lines in a wavelength region
llist = CIE.return_linelist(te, [1,2.0])
# print these to screen
# print to screen, listing the energy, not the wavelength
print("Listing lines between 1 and 2 A, first 10 lines only.")

print(CIE.format_linelist(llist[:10]))

