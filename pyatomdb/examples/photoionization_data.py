import pyatomdb, numpy,os, pylab
try:
  import astropy.io.fits as pyfits
except:
  import pyfits

# This is a sample routine that reads in the photoionization data
# It also demonstrates using get_data, which should download the data you 
# need automagically from the AtomDB site.
#
# It also shows how to get the raw XSTAR PI cross sections. 

# going to get PI cross section from iron 16+ to 17+ (Fe XVII-XVIII)
Z = 26
z1 = 17

# get the AtomDB level data
lvdata = pyatomdb.atomdb.get_data(Z, z1, 'LV')

# get the XSTAR PI data from AtomDB
pidata = pyatomdb.atomdb.get_data(Z, z1, 'PI')

# set up the figure
fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)


# to calculate the cross section (in cm^2) at a given single energy E (in keV)
# does not currently work with vector input, so have to call in a loop if you 
# want multiple energies [I will fix this]

E = 10.

# get the ground level (the 0th entry in LV file) data
lvd = lvdata[1].data[0]

# This is the syntax for calculating the PI cross section of a given line
# This will work for non XSTAR data too.
sigma = pyatomdb.atomdb.sigma_photoion(E, Z, z1,lvd['phot_type'], lvd['phot_par'], \
               xstardata=pidata, xstarfinallev=1)


# To get the raw XSTAR cross sections (units: energy = keV, cross sections = Mb)               
# for level 1 -> 1 (ground to ground)
pixsec = pyatomdb.atomdb.sort_pi_data(pidata, 1,1)
ax.loglog(pixsec['energy'], pixsec['pi_param']*1e-18, label='raw xstar data')

# label the plot
ax.set_title('Plotting raw XSTAR PI cross sections. Fe XVII gnd to Fe XVIII gnd')
ax.set_xlabel("Energy (keV)")
ax.set_ylabel("PI cross section (cm$^{2}$)")

pylab.draw()
zzz=raw_input('press enter to continue')
