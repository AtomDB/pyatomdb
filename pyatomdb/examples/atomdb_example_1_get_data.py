import pyatomdb, time

# Example of using get data. I am going to get a range of different atomic data
#
# This routine returns the raw opened FITS data. No other processing is done.

Z = 8
z1 = 7

# get the energy levels of O VII

a = pyatomdb.atomdb.get_data(Z,z1,'LV')

# print out the first few lines
# note that for all the atomic data files, the useful data is in HDU 1.

print(a[1].data.names)

print(a[1].data[:5])

# get radiative transition information
b = pyatomdb.atomdb.get_data(Z,z1,'LA')

# find the transitions from level 7

bb = b[1].data[ (b[1].data['UPPER_LEV']==7)]
print()
print(bb.names)
print(bb)

# and so on.

# For non ion-specific data, Z and z1 are ignored.

abund = pyatomdb.atomdb.get_data(0,0,'ABUND')
print()
print(abund[1].data.names)
print(abund[1].data[:2])


# How the datacache works
# The datacache variable stores opened files in memory, preventing
# having to go and re-open them each time. This reduces disk access.

# make sure file is downloaded for a fair test:
a = pyatomdb.atomdb.get_data(Z,z1,'EC')

# case 1: no datacache
t1 = time.time()
for i in range(10):
  a = pyatomdb.atomdb.get_data(Z,z1,'EC')
t2 = time.time()

# case 2: using datacache
datacache = {}
t3 = time.time()
for i in range(10):
  a = pyatomdb.atomdb.get_data(Z,z1,'EC', datacache=datacache)
t4 = time.time()

print("Time without / with datacache: %f / %f seconds"%(t2-t1, t4-t3))
