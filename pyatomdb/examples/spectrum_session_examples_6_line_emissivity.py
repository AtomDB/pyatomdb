import pyatomdb
import numpy
import pylab


# declare the Collisional Ionization Equilibrium session
sess = pyatomdb.spectrum.CIESession()

# set response
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

# get the emissivity of the resonance line of O VII at 5 different temperatures

kTlist = numpy.array([0.1,0.3,0.5,0.7,0.9])
Z=8
z1=7
up=7
lo = 1
ldata = sess.return_line_emissivity(kTlist, Z, z1, up, lo)

fig = pylab.figure()
fig.show()
ax= fig.add_subplot(111)

ax.semilogy(ldata['Te'], ldata['epsilon'])

ax.set_xlabel("Temperature (keV)")
ax.set_ylabel("Epsilon ph cm$^3$ s$^{-1}$")

pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('spectrum_session_examples_6_1.pdf')
fig.savefig('spectrum_session_examples_6_1.svg')
