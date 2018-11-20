import pyatomdb, pylab, numpy, time


rmf = '/export1/projects/atomdb_308/hitomi/resp_100041010sxs.rmf'
arf = '/export1/projects/atomdb_308/hitomi/arf_100041010sxs.arf'


# create a session object. This contains the apec files, response files,
# and any previously calculated spectrs so that simple multiplication
# can be used to get the results without recalculating everything from scratch

data = pyatomdb.spectrum.Session()

# If you want to specify custom energy bins:
ebins = numpy.linspace(1,2,1001)
data.set_specbins(ebins, specunits='A')

# alternative method: just load the response and use its binning. Note
# that this will always be in keV currently, because reasons.

data.set_response(rmf, arf=arf)
ebins = data.ebins_response

# now get the spectrum at 4keV. This calculates (and stores) the 
# spectrum at each temperature nearby (~3.7, 4.3 keV)
# then linearly interpolates between the result
#
# vector is stored for each element at each temperature
# so if you change temperature/abundance, it's a simple multiplication and
# interpolation instead of a total recalculation

t0 = time.time()
s=data.return_spectra(4.0, teunit='keV')
t1 = time.time()
# let's change the abundance
data.set_abund([1,2,3,4,26],0.5)

# and see how fast this goes this time, changing temperature and abund
s2=data.return_spectra(4.1, teunit='keV')
t2 = time.time()

print("first spectrum took %g seconds" %(t1-t0))
print("second spectrum took %g seconds" %(t2-t1))
print("note how much faster the second one was as I didn't recalculate everything from scratch!")


#linedata = pyatomdb.pyfits.open('/export1/atomdb_latest/apec_v3.0.8_line.fits')


#spec = speclo*rlo + specup*rup

# some plotting of things

fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)
#ax2 = fig.add_subplot(212, sharex=ax)
s = numpy.append(0,s)
s += 1e-40

s2 = numpy.append(0,s2)
s2 += 1e-40
#spec = numpy.append(0,spec)
#spec += 1e-40

ax.plot(ebins, s, drawstyle='steps')
ax.plot(ebins, s2, drawstyle='steps')
#ax.plot(elo, spec, drawstyle='steps')


#ax2.plot(elo, s/spec, drawstyle='steps')


zzz=input('Press enter to exit')
