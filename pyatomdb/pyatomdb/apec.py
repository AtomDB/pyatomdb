"""
The apec module contains routines crucial for the APEC code. This also
includes some interfaces to external C libraries (or will, eventually).

Version 0.1 - initial release
Adam Foster September 16th 2015

"""

import numpy, copy
import const, atomdb
import scipy

def calc_full_ionbal(Te, tau, init_pop=False, Te_init=False, Zlist=False, teunit='K'):
  """
  Calculate the ionization balance for all the elements in Zlist.
  
  One of init_pop or Te_init should be set. If neither is set, assume 
  all elements start from neutral.
  
  
  Parameters
  ----------
  Te : float
    electron temperature in keV or K (default K)
  tau : float
    N_e * t for the non-equilibrium ioniziation
  init_pop : dict of float arrays, indexed by Z
    initial populations. E.g. init_pop[6]=[0.1,0.2,0.3,0.2,0.2,0.0,0.0]
  Te_init : float
    initial ionization balance temperature, same units as Te
  Zlist : int array
    array of nuclear charges to include in calculation (e.g. [8,26] for 
    oxygen and iron)
  teunit : {'K' , 'keV'}
    units of temperatures (default K) 
  
  Returns
  -------
  final_pop : dict of float arrays, indexed by Z
    final populations. E.g. final_pop[6]=[0.1,0.2,0.3,0.2,0.2,0.0,0.0]
  
  """
  
  # input checking
  if teunit.lower() == 'kev':
    kT = Te*1.0
  elif teunit.lower() == 'ev':
    kT = Te/1000.0
  elif teunit.lower() == 'k':
    kT = Te*const.KBOLTZ
  else:
    print "*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunits)

  # input checking
  if Te_init != False:
    if teunit.lower() == 'kev':
      kT_init = Te_init*1.0
    elif teunit.lower() == 'ev':
      kT_init = Te_init/1000.0
    elif teunit.lower() == 'k':
      kT_init = Te_init*const.KBOLTZ
    else:
      print "*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
            (teunits)
  if not Zlist:
    Zlist = range(1,29)

  if (not Te_init) & (not init_pop):
    print "Warning: you have not specified an initial temperature or "+\
          "ion population: assuming everything is neutral"
    init_pop={}
    for Z in Zlist:
      init_pop[Z] = numpy.zeros(Z+1, dtype=float)
      init_pop[Z][0] = 1.0

  
  if (Te_init!=False) & (init_pop!=False):
    print "Warning: you have specified both an initial temperature and "+\
          "ion population: using ion population"
  
  datacache={}
  
  if not init_pop:
    init_pop = {}
    for Z in Zlist:
      ionrate = numpy.zeros(Z, dtype=float)
      recrate = numpy.zeros(Z, dtype=float)
      
      for z1 in range(1,Z+1):
        ionrate[z1-1], recrate[z1-1] = \
          atomdb.get_ionrec_rate(kT_init, False,  Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache)
      # now solve
      init_pop[Z] = solve_ionbal(ionrate, recrate)
      
  # now solve the actual ionization balance we want.
  pop = {}
  for Z in Zlist:
    ionrate = numpy.zeros(Z, dtype=float)
    recrate = numpy.zeros(Z, dtype=float)
    for z1 in range(1,Z+1):
      ionrate[z1-1], recrate[z1-1] = \
        atomdb.get_ionrec_rate(kT, False,  Te_unit='keV', \
                   Z=Z, z1=z1, datacache=datacache)
    
    # now solve
    pop[Z] = solve_ionbal(ionrate, recrate, init_pop=init_pop[Z], tau=tau)
    
  return pop
      

  
  
  
def solve_ionbal(ionrate, recrate, init_pop=False, tau=False):
  """
  solve_ionbal: given a set of ionization and recombination rates, find
  the equilibrium ionization balance. If init_pop and tau are set, do an
  non-equilibrium calculation starting from init_pop and evolving for 
  n_e * t = tau (cm^-3 s)
  
  Parameters
  ----------
  ionrate : float array
    the ionization rates, starting with neutral ionizing to +1
  recrate : float array
    the recombination rates, starting with singly ionized recombining to neutral
  init_pop : float array
    initial population of ions for non-equlibrium calculations. Will be renormalised to 1.
  tau : float
    N_e * t for the non-equilibrium ioniziation
  
  Returns
  -------
  final_pop : float array
    final populations.
  
  Notes
  -----
  Note that init_pop & final_pop will have 1 more element than ionrate and recrate.
  
  """
#  
#  Version 0.1 Initial Release
#  Adam Foster 16th September 2015
#
  # first, calculate the equilibrium solution
  
#  if (init_pop==False) & (tau==False): do_equilib=True
  try:
    if init_pop==False:
      do_equilib=True
    else:
      do_equilib=False
  except ValueError:
    do_equilib=False
  Z = len(ionrate)
  b = numpy.zeros(Z+1, dtype=numpy.float32)
  a = numpy.zeros((Z+1,Z+1), dtype=numpy.float32)
  
  
  for iZ in range(0,Z):
    a[iZ,iZ] -= (ionrate[iZ])
    a[iZ+1,iZ] += (ionrate[iZ])

    a[iZ,iZ+1] += (recrate[iZ])
    a[iZ+1,iZ+1] -= (recrate[iZ])

  # conservation of population
  for iZ in range(0,Z+1):
    a[0,iZ]=1.0
  b[0]=1.0
  
  lu,piv,eqpop,info=scipy.linalg.lapack.dgesv(a,b)

  eqpop[eqpop<0] = 0.0
  eqpop[0] = max([1.0-sum(eqpop[1:]), 0])
  
  if do_equilib == True:
    return eqpop
  
  # now the NEI part
  #remake matrix a
  Z=len(ionrate)+1
  ndim=Z
  AA = numpy.zeros((ndim-1,ndim-1), dtype=numpy.float32)
  # populate with stuff

  for iCol in range(ndim-1):
    for iRow in range(ndim-1):
      
      if (iRow==0):
        if (iCol==0):
          if (Z>1):
            AA[0,iCol] = -(ionrate[0] + ionrate[1] + recrate[0])
          else:
            AA[0,iCol] = -(ionrate[0] + recrate[0])
          
        
        if (iCol==1): AA[0,iCol] = recrate[1] - ionrate[0]
        if (iCol>1):
          AA[0,iCol] = -ionrate[0]
      else:
        if (iRow==iCol+1):  AA[iRow,iCol]= ionrate[iRow]
        if (iRow==iCol):
          if (iRow+2<ndim):
            
            AA[iRow,iCol]=-(ionrate[iRow+1]+recrate[iRow])
          else:
            AA[iRow,iCol]=-recrate[iRow]
          
        
        if (iRow==iCol-1):
           AA[iRow,iCol]= recrate[iRow+1]
  
  
  wr_la,wi_la,vl_la,vr_la, info=scipy.linalg.lapack.dgeev(AA)

  
  # now copy VR to AA
  AA=vr_la
  
  # b_in is difference betwen initial and final populations
  b_in=numpy.float32(init_pop[1:])-eqpop[1:]
  
  # solve for initial conditions
  lu,piv,b_out,info=scipy.linalg.lapack.dgesv(AA,b_in)

  # include exponential decay term
  C = b_out*numpy.exp( wr_la*numpy.float32(tau) )
  
  # get change
  G = numpy.dot(vr_la,C)

  # solve for population
  Ion_pop=numpy.zeros(Z)
  Ion_pop[1:]= G + eqpop[1:]
  
  # make sure everything is > 0
  Ion_pop[Ion_pop<0] = 0.0
  
  # set neutral to be the residual (or zero)
  Ion_pop[0]=max([1-sum(Ion_pop[1:]), 0.0])

  # return the data
  return Ion_pop
  
