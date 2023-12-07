import pyatomdb, numpy, os, pylab
import astropy.io.fits as pyfits

"""
This code is an example of generating a cooling curve: the total power
radiated in keV cm3 s-1 by each element at each temperature. It will
generate a text file with the emission per element at each temperature
from 1e4 to 1e9K. This has been modified for a non-equilibrium plasma.

This is similar to the atomdb.lorentz_power function, but with a few
normalizations removed to run a little quicker.

Note that the Anders and Grevesse (1989) abundances are built in to
this. These can be looked up using atomdb.get_abundance(abundset='AG89'),
or the 'angr' column of the table at
https://heasarc.nasa.gov/xanadu/xspec/xspec11/manual/node33.html#SECTION00631000000000000000

Adjustable parameters (energy range, element choice) are in the block
marked #### ADJUST THINGS HERE

Usage: python3 calc_power_nei.py
"""


def calc_power(Zlist, nei, Tlist, Taulist):
  init_pop = 'ionizing'
  res = {}
  res['power'] = numpy.zeros([len(Tlist), len(Taulist)])
  res['temperature'] = []
  res['tau'] = []
  
  # get photon energy of each bin
  en = nei.ebins_out
  en = (en[1:]+en[:-1])/2
  
  for i, T in enumerate(Tlist):
    kT = T/(1000*11604.5)
    res['temperature'].append(kT)
    for j, Tau in enumerate(Taulist):
      print("Doing temperature iteration %i of %i"%(i, len(Tlist)))
      
      res['tau'].append(Tau)


      # This code would go through each element 1 by 1.
      # We will do all elements at the same time to simplify
      
#      for Z in Zlist:
#        if Z==0:
#          # This is the electron-electron bremstrahlung component alone
#
#          #set all abundances to 1 (I need a full census of electrons in the plasma for e-e brems)
#          nei.set_abund(Zlist[1:], 1.0)
#          # turn on e-e bremsstrahlung
#          nei.set_eebrems(True)
#          spec = nei.return_spectrum(kT, Tau, init_pop=init_pop, \
#                                     dolines=False, docont=False, dopseudo=False)
#        else:
#          # This is everything else, element by element.
#
#          # turn off all the elements
#          nei.set_abund(Zlist[1:], 0.0)
#          # turn back on this element
#          nei.set_abund(Z, 1.0)
#          # turn off e-e bremsstrahlung (avoid double counting)
#          nei.set_eebrems(False)
#
#          spec = nei.return_spectrum(kT, Tau, init_pop=init_pop)
#

      spec = nei.return_spectrum(kT, Tau, init_pop=init_pop)
      E = spec*en
      res['power'][i][j] = sum(E)
  return(res)

if __name__=='__main__':


  ############################
  #### ADJUST THINGS HERE ####
  ############################

  # Elements to include
  #Zlist = range(31) <- all the elements
  Zlist = [0,1,2,6,7,8,10,12,13,14,16,18,20,26,28] #<- just a few
  # Note that for this purpose, Z=0 is the electron-electron bremsstrahlung
  # continuum. This is not a general AtomDB convention, just what I've done here
  # to make this work.

  # specify photon energy range you want to integrate over (min = 0.001keV, max=100keV)
  Elo = 0.001 #keV
  Ehi = 100.0 #

  # temperatures at which to calculate curve (K)
  Tlist = numpy.logspace(4,9,51)

  # specify output file name (default output.txt)
  outfile = 'output.txt'

  ################################
  #### END ADJUST THINGS HERE ####
  ################################

  # set up the spectrum
  Zlist = numpy.array([0,1,2,6,7,8,10,12,13,14,16,18,20,26,28]) #<- just a few
  
  nei = pyatomdb.spectrum.NEISession(elements=Zlist[Zlist>0])
  if 0 in Zlist:
    nei.set_eebrems(True)
  ebins = numpy.linspace(Elo, Ehi, 10001)
  nei.set_response(ebins, raw=True)
  nei.set_eebrems(True)
  # crunch the numbers
  Taulist = [1e8, 1e10, 1e14]
  k = calc_power(Zlist, nei, Tlist, Taulist)


  fig = pylab.figure()
  fig.show()
  ax = fig.add_subplot(111)

  for i, tau in enumerate(Taulist):
    ax.loglog(k['temperature'], k['power'][:,i]*pyatomdb.const.ERG_KEV, label=repr(int(numpy.log10(tau))))
  ax.legend(loc=0)
  pylab.draw()
  
  ax.set_xlabel("Temperature (keV)")
  ax.set_ylabel("Radiated Power 0.001-100keV\n(erg cm$^{3}$ s$^{-1}$)")
  

  zzz=input("Press enter to continue")
  fig.savefig('calc_power_nei_examples_1_1.pdf')
  fig.savefig('calc_power_nei_examples_1_1.svg')
