import pyatomdb
import numpy
import pylab


# declare the Non-Equilibrium Ionization session. I am only using
# oxygen and neon for simplicity
sess = pyatomdb.spectrum.NEISession(elements=[8,10])

# set the response (raw keyword tells pyatomdb it is not a real response file)
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

kT = 0.4 # temperature in keV
tau = 1e10 # n_e * t in cm^-3 s
spec = sess.return_spectrum(kT, tau, init_pop=1.0)

# Returned spectrum has units of photons cm^5 s^-1 bin^-1
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

offset = max(spec)
ax.plot(sess.ebins_out, numpy.append(0,spec), drawstyle='steps', label='mild recombining')


spec = sess.return_spectrum(kT, tau, init_pop='ionizing')
ax.plot(sess.ebins_out, offset+numpy.append(0,spec), drawstyle='steps', label='ionizing')

spec = sess.return_spectrum(kT, tau, init_pop='recombining')
ax.plot(sess.ebins_out, 2*offset+numpy.append(0,spec), drawstyle='steps', label='recombining')



# manually specified initial population
init_pop_dict = {}
init_pop_dict[8] = numpy.zeros(9) # 9 stages of oxygen (include neutral an fully stripped)
init_pop_dict[10] = numpy.zeros(11)
#put everything in the H-like stage
init_pop_dict[8][8] = 1.0
init_pop_dict[10][10] = 1.0


spec = sess.return_spectrum(kT, tau, init_pop=init_pop_dict)
ax.plot(sess.ebins_out, 3*offset+numpy.append(0,spec), drawstyle='steps', label='all H-like')


# label the figure
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')

# zoom in a bit
ax.set_xlim([0.6,1.2])
ax.legend(loc=0)
# make plot appear
pylab.draw()

zzz=input("Press enter to continue")

# save image files
fig.savefig('spectrum_session_examples_7_1.pdf')
fig.savefig('spectrum_session_examples_7_1.svg')

