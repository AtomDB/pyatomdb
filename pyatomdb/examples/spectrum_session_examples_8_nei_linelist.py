import pyatomdb
import numpy
import pylab


# declare the Non-Equilibrium Ionization session. I am only using
# oxygen and neon for simplicity
sess = pyatomdb.spectrum.NEISession()

# set the response (raw keyword tells pyatomdb it is not a real response file)
sess.set_response('aciss_meg1_cy22.grmf', arf = 'aciss_meg1_cy22.garf')

kT = 4.0 # temperature in keV
tau = 1e11 # n_e * t in cm^-3 s
llist = sess.return_linelist(kT, tau, [1.7,1.9],init_pop=1.0)

# filter lines to strongest 20
llist.sort(order='Epsilon')
llist= llist[-20:]

print("strongest 20 lines")
print(llist.dtype.names)
print(llist)

llist = sess.return_linelist(kT, tau, [1.7,1.9],by_ion_drv=True,init_pop=1.0)

# filter lines to strongest 20
llist.sort(order='Epsilon')
llist= llist[-20:]

print("Strongest 20 lines, separated by driving ion")
print(llist.dtype.names)
print(llist)


# Expected output
# strongest 20 lines
# ('Lambda', 'Lambda_Err', 'Epsilon', 'Epsilon_Err', 'Element', 'Ion', 'UpperLev', 'LowerLev')
# [(1.8827   , nan, 4.3892664e-19, nan, 26, 22, 20009, 1)
#  (1.8622   , nan, 4.5014698e-19, nan, 26, 24, 10242, 3)
#  (1.8824   , nan, 4.8421002e-19, nan, 26, 22,   483, 1)
#  (1.8825   , nan, 5.8088691e-19, nan, 26, 22,   482, 1)
#  (1.8851   , nan, 5.9669826e-19, nan, 26, 22, 20008, 2)
#  (1.8747   , nan, 6.6686722e-19, nan, 26, 24,    44, 1)
#  (1.8737   , nan, 6.9597891e-19, nan, 26, 23, 20126, 3)
#  (1.863    , nan, 1.0866564e-18, nan, 26, 24, 10247, 2)
#  (1.8756001, nan, 1.0874650e-18, nan, 26, 23, 20122, 4)
#  (1.8706   , nan, 1.3853237e-18, nan, 26, 24,    46, 1)
#  (1.8738   , nan, 1.4878691e-18, nan, 26, 24,    45, 1)
#  (1.8659   , nan, 1.6188786e-18, nan, 26, 24, 10248, 3)
#  (1.857    , nan, 1.9076071e-18, nan, 26, 24,    50, 1)
#  (1.8554125, nan, 2.3919818e-18, nan, 26, 25,     6, 1)
#  (1.8595169, nan, 2.8080335e-18, nan, 26, 25,     5, 1)
#  (1.8635   , nan, 3.6176967e-18, nan, 26, 24,    47, 1)
#  (1.8704   , nan, 7.2666667e-18, nan, 26, 23,   309, 1)
#  (1.8681941, nan, 8.8379641e-18, nan, 26, 25,     2, 1)
#  (1.8611   , nan, 1.1892281e-17, nan, 26, 24,    48, 1)
#  (1.8503995, nan, 1.4463351e-17, nan, 26, 25,     7, 1)]
# Strongest 20 lines, separated by driving ion
# ('Lambda', 'Lambda_Err', 'Epsilon', 'Epsilon_Err', 'Element', 'Elem_drv', 'Ion', 'Ion_drv', 'UpperLev', 'LowerLev')
# [(1.8622   , nan, 4.5014698e-19, nan, 26, 26, 24, 25, 10242, 3)
#  (1.8824   , nan, 4.8421002e-19, nan, 26, 26, 22, 22,   483, 1)
#  (1.8825   , nan, 5.6233666e-19, nan, 26, 26, 22, 22,   482, 1)
#  (1.8851   , nan, 5.9669826e-19, nan, 26, 26, 22, 23, 20008, 2)
#  (1.8747   , nan, 6.5364458e-19, nan, 26, 26, 24, 24,    44, 1)
#  (1.8737   , nan, 6.9597891e-19, nan, 26, 26, 23, 24, 20126, 3)
#  (1.863    , nan, 1.0866564e-18, nan, 26, 26, 24, 25, 10247, 2)
#  (1.8756001, nan, 1.0874650e-18, nan, 26, 26, 23, 24, 20122, 4)
#  (1.8706   , nan, 1.3512313e-18, nan, 26, 26, 24, 24,    46, 1)
#  (1.8738   , nan, 1.4619381e-18, nan, 26, 26, 24, 24,    45, 1)
#  (1.8659   , nan, 1.6188786e-18, nan, 26, 26, 24, 25, 10248, 3)
#  (1.857    , nan, 1.8978882e-18, nan, 26, 26, 24, 24,    50, 1)
#  (1.8554125, nan, 2.3303760e-18, nan, 26, 26, 25, 25,     6, 1)
#  (1.8595169, nan, 2.7651213e-18, nan, 26, 26, 25, 25,     5, 1)
#  (1.8681941, nan, 3.3107048e-18, nan, 26, 26, 25, 25,     2, 1)
#  (1.8635   , nan, 3.6059715e-18, nan, 26, 26, 24, 24,    47, 1)
#  (1.8681941, nan, 5.3853784e-18, nan, 26, 26, 25, 24,     2, 1)
#  (1.8704   , nan, 7.2068723e-18, nan, 26, 26, 23, 23,   309, 1)
#  (1.8611   , nan, 1.1865456e-17, nan, 26, 26, 24, 24,    48, 1)
#  (1.8503995, nan, 1.4405425e-17, nan, 26, 26, 25, 25,     7, 1)]
# 
