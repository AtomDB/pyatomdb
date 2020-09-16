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

# get a list of the available abundance sets
sess.set_abundset()
# this prints them to screen.



for iabset, abset in enumerate(['Allen', 'AG89', 'GA88', 'Feldman', 'GA10', 'Lodd09', 'AE82']):

  # set abundset
  sess.set_abundset(abset)

  # get spectrum at 1keV
  s1 = sess.return_spectrum(1.0)
  if iabset == 0:
    baseoffset = max(s1)

  offset = iabset*baseoffset


  #plot
  ax.plot(sess.ebins_out, numpy.append(0,s1)+offset, drawstyle='steps', label=abset)

ax.set_xlim([0.8,2.1])

ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')
ax.legend(loc=0, ncol=4)
pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('spectrum_session_examples_4_1.pdf')
fig.savefig('spectrum_session_examples_4_1.svg')

