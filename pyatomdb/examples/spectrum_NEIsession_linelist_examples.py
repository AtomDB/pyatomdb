import pyatomdb, numpy, pylab


# declare the Non-Equilibrium Ionization session
sess = pyatomdb.spectrum.NEISession()

kT = 0.8 # temperature in keV
tau = 1e11 # n_e*t cm^-3 s
kT_init=0.01 # initial temperature (for initial ion distribution)

waverange = [21,27]

# find the lines in a certain wavelength range
llist = sess.return_linelist(kT,tau,  waverange)

# filter out to get only lines within 10% of the strongest
llist = llist[ llist['Epsilon'] > max(llist['Epsilon'])*0.1]

# print the results:

# column names
print(llist.dtype.names)

# data
print(llist)


# now let's see how one of the lines varies as  a function of temperature. Let's look
# at the O VII resonance line (7 -> 1)

# list of temperatures to look at
kTlist = numpy.linspace(0.05,3.0,20)

# get the emissivity as a function of temperature (returned in a dict)
ldata = sess.return_line_emissivity(kTlist, 8, 7, 7, 1)

# print out the data
for i in range(len(kTlist)):
  print("%.3f keV %e ph cm^3 s^-1"%(ldata['Te'][i], ldata['epsilon'][i]))
