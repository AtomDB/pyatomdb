"""
The apec module contains routines crucial for the APEC code. This also
includes some interfaces to external C libraries (or will, eventually).

Version 0.1 - initial release
Adam Foster September 16th 2015

"""

import numpy, copy

import scipy

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
  
