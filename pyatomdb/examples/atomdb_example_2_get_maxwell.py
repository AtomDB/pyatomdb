import pyatomdb, numpy, pylab

# We will get some maxwellian rates for a O 6+ 7->1

Z = 8
z1 = 7
up = 7
lo = 1

# excitation

Te = numpy.logspace(6, 7, 15)

#  >>> # (1) Get excitation rates for row 12 of an Fe XVII file
#  >>> colldata = pyatomdb.atomdb.get_data(26,17,'EC')
#  >>> exc, dex = get_maxwell_rate(Te, colldata=colldata, index=12)

#  >>> # (2) Get excitation rates for row 12 of an Fe XVII file
#  >>> exc, dex = get_maxwell_rate(Te, Z=26,z1=17, index=12)

#  >>>  (3) Get excitation rates for transitions from level 1 to 15 of FE XVII
#  >>> exc, dex = get_maxwell_rate(Te, Z=26, z1=17, dtype='EC', finallev=15, initlev=1)

datacache = {}

# get data by specifying ion, upper and lower levels
exc, dexc = pyatomdb.atomdb.get_maxwell_rate(Te, Z=Z, z1=z1, initlev = lo, finallev=up, dtype='EC', datacache=datacache)

fig= pylab.figure()
fig.show()
ax = fig.add_subplot(111)
ax.loglog(Te, exc, label='excitation')
ax.loglog(Te, dexc, label = 'deexcitation')


# preload data and find a specific transition
ecdat = pyatomdb.atomdb.get_data(8,7,'EC', datacache=datacache)

i = numpy.where( (ecdat[1].data['Upper_Lev']==up) &\
                 (ecdat[1].data['Lower_Lev']==lo))[0][0]

exc, dex = pyatomdb.atomdb.get_maxwell_rate(Te,colldata=ecdat, index=i, datacache=datacache)

ax.loglog(Te, exc, 'o', label='excitation')
ax.loglog(Te, dexc, 'o', label = 'deexcitation')

ax.legend(loc=0)
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Rate Coefficient (cm$^3$ s$^{-1}$)")

pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('atomdb_examples_2_1.pdf')
fig.savefig('atomdb_examples_2_1.svg')

# you can also obtain ionization or recombination rates (see get_maxwell_rates writeup for options)
# in most cases you need to specify initlev and finallev as 1 to get the ion to ion rates.

ion= pyatomdb.atomdb.get_maxwell_rate(Te,Z=Z, z1=z1, initlev=1, finallev=1, dtype='CI', datacache=datacache)
print(ion)
#
# However it is recommended that you use get_ionrec_rate instead - it gives you everythign in one go

# combining the different types
ion, rec = pyatomdb.atomdb.get_ionrec_rate(Te, Z=Z, z1=z1, datacache=datacache)

# as separate entities
CI, EA, RR, DR=pyatomdb.atomdb.get_ionrec_rate(Te, Z=Z, z1=z1, datacache=datacache, separate=True)

print("Collisional ionization: ",CI)
print("Excitation Autoionization: ",EA)
print("Radiative Recombination: ",RR)
print("Dielectroni Recombination: ",DR)
ax.cla()

ax.loglog(Te, ion, label='Ionization')
ax.loglog(Te, rec, label='Recombination')

ax.loglog(Te, CI, '--',label='CI')
ax.loglog(Te, EA, '--',label='EA')
ax.loglog(Te, RR, '--',label='RR')
ax.loglog(Te, DR, '--',label='DR')

ax.legend(loc=0)
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Rate Coefficient (cm$^3$ s$^{-1}$)")

pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('atomdb_examples_2_2.pdf')
fig.savefig('atomdb_examples_2_2.svg')


#get_maxwell_rate(Te, colldata=False, index=-1, lvdata=False, Te_unit='K', \
#                     lvdatap1=False, ionpot = False, \
#                     force_extrap=False, silent=True,\
#                     finallev=False, initlev=False,\
#                     Z=-1, z1=-1, dtype=False, exconly=False,\
#                     datacache=False, settings=False, ladat=False):
