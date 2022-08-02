import pyatomdb
import numpy
import pylab


# declare the Collisional Ionization Equilibrium session
sess = pyatomdb.spectrum.CIESession()


kT = 0.4 # temperature in keV

# Returned spectrum has units of photons cm^5 s^-1 bin^-1

# create figure
gs = pylab.matplotlib.gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 
fig = pylab.figure()
fig.show()
ax1 = fig.add_subplot(gs[0]) # for data
ax2 = fig.add_subplot(gs[1], sharex=ax1) # for ratio

# label the figure
ax2.set_xlabel('Energy (keV)')
ax1.set_ylabel('Intensity (ph cm$^5$ s$^{-1}$ bin$^{-1}$)')
ax2.set_ylabel('Ratio')

# set response - using a regular matrix
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

# get spectrum, prepend zero
spec = sess.return_spectrum(kT)
spec = numpy.append(spec[0], spec)

ax1.plot(sess.ebins_out, spec, drawstyle='steps', label='Regular')

# change the response to a sparse one
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf',
                  sparse=True)
spec_sparse = sess.return_spectrum(kT)
spec_sparse = numpy.append(spec_sparse[0], spec_sparse)

ax1.plot(sess.ebins_out, spec_sparse, '--', drawstyle='steps', label='Sparse')

# plot the ratio
ax2.plot(sess.ebins_out, spec/spec_sparse, drawstyle='steps',\
         label='Regular/Sparse')

# add legend
ax1.legend(loc=0)
ax2.legend(loc=0)

# zoom in on small sections of the spectrum
ax1.set_xlim([0.5, 2.0])

# adjust plot spacing so labels are visible
pylab.matplotlib.pyplot.subplots_adjust(hspace = 0)

# draw graphs
pylab.draw()

zzz=input("Press enter to continue")

# save image files
fig.savefig('spectrum_session_examples_1b_1.pdf')
fig.savefig('spectrum_session_examples_1b_1.svg')

