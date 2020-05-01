import pyatomdb
import numpy
import pylab


# declare the Collisional Ionization Equilibrium session
sess = pyatomdb.spectrum.CIESession()

# set response
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

# get a linelist between 1.85 and 1.88 Angstroms for a kT=4keV plasma

llist = sess.return_linelist(4.0,[1.85,1.88])

# select only the strongest 10 lines
llist.sort(order='Epsilon')
llist = llist[-10:]

# show the column names
print(llist.dtype.names)

# expected output:
# ('Lambda', 'Lambda_Err', 'Epsilon', 'Epsilon_Err', 'Element', 'Ion', 'UpperLev', 'LowerLev')

# print the data
print(llist)

# expected output:
#[(1.8635   , nan, 8.8202775e-19, nan, 26, 24,    47, 1)
# (1.8532   , nan, 1.0113579e-18, nan, 26, 24, 10089, 6)
# (1.8622   , nan, 1.0681659e-18, nan, 26, 24, 10242, 3)
# (1.863    , nan, 2.5785468e-18, nan, 26, 24, 10247, 2)
# (1.8611   , nan, 2.8751725e-18, nan, 26, 24,    48, 1)
# (1.8659   , nan, 3.8414809e-18, nan, 26, 24, 10248, 3)
# (1.8554125, nan, 7.2097873e-18, nan, 26, 25,     6, 1)
# (1.8595169, nan, 7.7184372e-18, nan, 26, 25,     5, 1)
# (1.8681941, nan, 1.2001933e-17, nan, 26, 25,     2, 1)
# (1.8503995, nan, 3.5732330e-17, nan, 26, 25,     7, 1)]

# we can then do the same thing, but including the effective area of the instrument.
# So the Epsilon columnw will be multiplied by the effective area at the line's energy
llist = sess.return_linelist(4.0,[1.85,1.88], apply_aeff=True)

llist.sort(order='Epsilon')
llist = llist[-10:]

print('After apply effective area:')
print(llist)
#[(1.8635   , nan, 1.6238816e-18, nan, 26, 24,    47, 1)
# (1.8532   , nan, 1.7191134e-18, nan, 26, 24, 10089, 6)
# (1.8622   , nan, 1.9665764e-18, nan, 26, 24, 10242, 3)
# (1.863    , nan, 4.7473049e-18, nan, 26, 24, 10247, 2)
# (1.8611   , nan, 5.2934158e-18, nan, 26, 24,    48, 1)
# (1.8659   , nan, 7.3193854e-18, nan, 26, 24, 10248, 3)
# (1.8554125, nan, 1.2779108e-17, nan, 26, 25,     6, 1)
# (1.8595169, nan, 1.3680673e-17, nan, 26, 25,     5, 1)
# (1.8681941, nan, 2.2867946e-17, nan, 26, 25,     2, 1)
# (1.8503995, nan, 6.0738072e-17, nan, 26, 25,     7, 1)]
