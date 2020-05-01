import pyatomdb
import numpy
import pylab


# declare the Collisional Ionization Equilibrium session
sess = pyatomdb.spectrum.CIESession()

# set response
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

# for plotting
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

# make a spectrum at 1keV
s1 = sess.return_spectrum(1.0)

# change the abundance of Iron
sess.set_abund(26,2.0)
s2 = sess.return_spectrum(1.0)
# change the abundance of Magnesium and Silicon
sess.set_abund([12,14],0.5)
s3 = sess.return_spectrum(1.0)
# change the abundance of Neon and Iron separately
sess.set_abund([10,26],[2.0,0.5])
s4 = sess.return_spectrum(1.0)
# plot all this
offset = max(s1)
ax.plot(sess.ebins_out, numpy.append(0,s1), drawstyle='steps', label='Original')
ax.plot(sess.ebins_out, numpy.append(0,s2)+offset, drawstyle='steps', label='Fex2')
ax.plot(sess.ebins_out, numpy.append(0,s3)+offset*2, drawstyle='steps', label='Mg, Six0.5')
ax.plot(sess.ebins_out, numpy.append(0,s4)+offset*3, drawstyle='steps', label='Nex2, Fex0.5')
ax.legend(ncol=4, loc=0)
ax.set_xlim([0.8,1.5])
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')
pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('spectrum_session_examples_3_1.pdf')
fig.savefig('spectrum_session_examples_3_1.svg')

