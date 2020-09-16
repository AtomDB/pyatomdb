import pyatomdb, numpy, pylab

# Calculate an equilibrium ionization balance for an element
kT = 0.5
Z = 8

# going to declare a datacache to speed up the calculations as we
# will be repeating things
datacache = {}

ionfrac_fast = pyatomdb.apec.return_ionbal(Z, kT, teunit='keV', datacache=datacache)

# The "fast" keyword makes the code open up a precalculated set of
# ionization balance data, a grid os 1251 values from 10^4K <= T <= 10^9K.
# By default it is true. Here, we turn it off for comparison
ionfrac_slow = pyatomdb.apec.return_ionbal(Z, kT, teunit='keV', datacache=datacache,
                                      fast=False)

fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

ax.semilogy(range(0,Z+1), ionfrac_fast, 'o', label='fast')
ax.semilogy(range(0,Z+1), ionfrac_slow, '^', label='slow')
ax.set_ylim([1e-6,1])
ax.set_xlabel('Ion Charge')
ax.set_ylabel('Fraction')
ax.legend(loc=0)

pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('apec_examples_1_1.pdf')
fig.savefig('apec_examples_1_1.svg')


# now to do the same for several different ways of defining non-equilbrium:
tau = 1e11
# 1: from neutral
init_pop='ionizing' # this is the default
nei1 = pyatomdb.apec.return_ionbal(Z, kT, tau = tau, init_pop = init_pop,
                                   teunit='keV', datacache=datacache)


# 2: from fully stripped
init_pop='recombining'
nei2 = pyatomdb.apec.return_ionbal(Z, kT, tau = tau, init_pop = init_pop,
                                   teunit='keV', datacache=datacache)


# 3: from a given Temperature
init_pop=0.1
nei3 = pyatomdb.apec.return_ionbal(Z, kT, tau = tau, init_pop = init_pop,
                                   teunit='keV', datacache=datacache)


# 4: from a given array
init_pop=numpy.array([0.1, #O0+
                      0.1, #O1+
                      0.1, #O2+
                      0.1, #etc
                      0.3,
                      0.1,
                      0.1,
                      0.1,
                      0.0])
nei4 = pyatomdb.apec.return_ionbal(Z, kT, tau = tau, init_pop = init_pop,
                                   teunit='keV', datacache=datacache)

# 5: from a dict of arrays
# (this is useful if you have lots of elements to do)
init_pop = {}
init_pop[8] = numpy.array([0.5, #O0+
                      0.0, #O1+
                      0.0, #O2+
                      0.0, #etc
                      0.0,
                      0.0,
                      0.0,
                      0.0,
                      0.5])
nei5 = pyatomdb.apec.return_ionbal(Z, kT, tau = tau, init_pop = init_pop,
                                   teunit='keV', datacache=datacache)

ax.cla()
ax.semilogy(range(Z+1), nei1, 'o', label='ionizing')
ax.semilogy(range(Z+1), nei2, '<', label='recombining')
ax.semilogy(range(Z+1), nei3, 's', label='0.1keV')
ax.semilogy(range(Z+1), nei4, '^', label='array')
ax.semilogy(range(Z+1), nei5, 'x', label='dict')

ax.set_ylim([1e-6,1])
ax.set_xlabel('Ion Charge')
ax.set_ylabel('Fraction')
ax.legend(loc=0)

pylab.draw()
zzz=input("Press enter to continue")
# save image files
fig.savefig('apec_examples_1_2.pdf')
fig.savefig('apec_examples_1_2.svg')
