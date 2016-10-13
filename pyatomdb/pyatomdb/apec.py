"""
The apec module contains routines crucial for the APEC code. This also
includes some interfaces to external C libraries (or will, eventually).

Version 0.1 - initial release
Adam Foster September 16th 2015

"""

import numpy, copy, pickle
import util, atomdb, const, os, atomic, time
import scipy, ctypes
import astropy.io.fits as pyfits
from joblib import Parallel, delayed

import pylab

def calc_full_ionbal(Te, tau=1e14, init_pop=False, Te_init=False, Zlist=False, teunit='K',\
                    extrap=False, cie=True, settings=False):
  """
  Calculate the ionization balance for all the elements in Zlist.

  One of init_pop or Te_init should be set. If neither is set, assume
  all elements start from neutral.


  Parameters
  ----------
  Te : float
    electron temperature in keV or K (default K)
  tau : float
    N_e * t for the non-equilibrium ioniziation (default 1e14)
  init_pop : dict of float arrays, indexed by Z
    initial populations. E.g. init_pop[6]=[0.1,0.2,0.3,0.2,0.2,0.0,0.0]
  Te_init : float
    initial ionization balance temperature, same units as Te
  Zlist : int array
    array of nuclear charges to include in calculation (e.g. [8,26] for
    oxygen and iron)
  teunit : {'K' , 'keV'}
    units of temperatures (default K)
  extrap : bool
    Extrappolate rates to values outside their given range. (default False)
  cie : bool
    If true, collisional ionization equilbrium calculation
    (tau, init_pop, Te_init all ignored)

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
  if util.keyword_check(Te_init):
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

  if (not util.keyword_check(Te_init)) & (not util.keyword_check(init_pop)) &\
     (not util.keyword_check(cie)):
    print "Warning: you have not specified an initial temperature or "+\
          "ion population: assuming everything is neutral"
    init_pop={}
    for Z in Zlist:
      init_pop[Z] = numpy.zeros(Z+1, dtype=float)
      init_pop[Z][0] = 1.0


  if (util.keyword_check(Te_init)!=False) & (util.keyword_check(init_pop)!=False):
    print "Warning: you have specified both an initial temperature and "+\
          "ion population: using ion population"

  datacache={}

  if not init_pop:
    init_pop = {}
    for Z in Zlist:
      ionrate = numpy.zeros(Z, dtype=float)
      recrate = numpy.zeros(Z, dtype=float)
      if cie:
        kT_init=kT
      for z1 in range(1,Z+1):
        tmp = \
          atomdb.get_ionrec_rate(kT_init, False,  Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)
        
        
        ionrate[z1-1], recrate[z1-1]=tmp
      # now solve
      
      init_pop[Z] = solve_ionbal(ionrate, recrate)
      print "initial pop: ", init_pop
  if cie:
    return init_pop

  # now solve the actual ionization balance we want.
  pop = {}
  for Z in Zlist:
    ionrate = numpy.zeros(Z, dtype=float)
    recrate = numpy.zeros(Z, dtype=float)
    for z1 in range(1,Z+1):
      ionrate[z1-1], recrate[z1-1] = \
        atomdb.get_ionrec_rate(kT, False,  Te_unit='keV', \
                   Z=Z, z1=z1, datacache=datacache, extrap=extrap)

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
          if (Z>2):
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

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calc_brems_gaunt(E, T, z1, brems_type, datacache=False, \
                                          settings=False):
  """
  calculate the bremstrahulung free-free gaunt factor

  Parameters
  ----------
  E : float
    Energy (in keV) to calculate gaunt factor
  T : float
    Temperature (in K) of plasma
  z1 : int
    Ion charge +1 of ion (e.g. 6 for C VI)
  brems_type : int
    Type of bremstrahlung requested:
    1 = HUMMER = Non-relativistic: 1988ApJ...327..477H
    2 = KELLOGG = Semi-Relativistic: 1975ApJ...199..299K
    3 = RELATIVISTIC = Relativistic: 1998ApJ...507..530N
    4 = BREMS_NONE = no bremstrahlung
  settings : dict
    See description in atomdb.get_data
  datacache : dict
    Used for caching the data. See description in atomdb.get_data


  Returns
  -------
  gaunt_ff : float
    The gaunt factor for the free-free process.

  """

  Evec, Eisvec = util.make_vec(E)
  gaunt_ff = numpy.zeros(len(Evec), dtype=float)
  NUM_J=8
  NUM_I=11
  z0 = z1-1.0 # we need to use the actual ion charge
  if brems_type==const.HUMMER:
    # read in the hummer data
    hdat = atomdb.get_data(False, False,'hbrems', settings = settings,\
                                               datacache=datacache)
    gaunt_D=hdat['BR_GAUNT'].data['COEFFICIENT']

    gamma2 = z0**2 * const.RYDBERG/(const.KBOLTZ*T)
    if ((gamma2 < 1e-3) | (gamma2 > 1e3)):
      if (z0<10):
        print "brems_hummer: Warning, gamma^2 = %e is out of range."%(gamma2)
      gaunt_ff[:]=1.0

    else:
      u = Evec/(const.KBOLTZ*T)
      j =  numpy.where(u <1.e-4)[0]
      if len(j) != 0:
        print "brems_hummer: Warning, u is out of range: ", u[j]
        gaunt_ff[j]=1.0

      #j = numpy.where(u>31.6227766)[0]
      gaunt_ff[u>31.6227766]=1.0

      iii = numpy.where(gaunt_ff > 0.99)[0]
      if len(iii) > 0:
#        print "u too large at %i locations"%(len(iii))
#        print Evec[iii]
         pass
        #zzz=raw_input()
      j = numpy.where(gaunt_ff<1.0)[0]
      if len(j) > 0:
        x_u = (2*numpy.log10(u[j])+2.5)/5.5
        x_g = numpy.log10(gamma2)/3.0

        c_j = numpy.zeros(NUM_J, dtype=float)
        for jj in xrange(NUM_J):
          # We then sum the Chebyshev series, but only 0.5* the first value
          tmpgaunt = gaunt_D[jj*NUM_I:(jj+1)*NUM_I]
          tmpgaunt[0] *=0.5
          c_j[jj]=numpy.polynomial.chebyshev.chebval(x_g, tmpgaunt)
        # again, sum chebshev polynomial with first index * 0.5
        c_j[0]*=0.5
        for ij, jj in enumerate(j):
          gaunt_ff[jj] = numpy.polynomial.chebyshev.chebval(x_u[ij], c_j)
    if not Eisvec:
      gaunt_ff = gaunt_ff[0]
    return gaunt_ff
  elif brems_type==const.KELLOGG:

    gam2 = numpy.array([.7783,1.2217,2.6234, 4.3766,20.,70.])
    gam3 = numpy.array([1.,1.7783,3.,5.6234, 10.,30. ])
    a    = numpy.array([1.001, 1.004, 1.017, 1.036, 1.056, 1.121, 1.001, 1.005,\
                        1.017, 1.046, 1.073, 1.115, .9991, 1.005, 1.03, 1.055,\
                        1.102, 1.176, .997, 1.005, 1.035,  1.069, 1.134, 1.186,\
                        .9962, 1.004, 1.042, 1.1, 1.193, 1.306,.9874, .9962,\
                        1.047, 1.156, 1.327, 1.485, .9681, .9755, 1.02, 1.208,\
                        1.525, 1.965, .3029, .1616, .04757, .013, .0049, -.0032,\
                        .4905, .2155, .08357, .02041, .00739, 2.9e-4, .654, \
                        .2833, .08057, .03257, .00759, -.00151, 1.029, .391, \
                        .1266, .05149, .01274, .00324, .9569, .4891, .1764, \
                        .05914, .01407, -2.4e-4, 1.236, .7579, .326, .1077, \
                        .028, .00548, 1.327, 1.017,.6017, .205, .0605, .00187, \
                        -1.323, -.254, -.01571, -.001, -1.84e-4, 8e-5, -4.762, \
                        -.3386, -.03571, -.001786, -3e-4, 1e-5, -6.349, -.4206, \
                        -.02571, -.003429, -2.34e-4, 5e-5, -13.231, -.59, \
                        -.04571, -.005714, -4.45e-4, -4e-5, -7.672, -.6852, \
                        -.0643, -.005857,  -4.2e-4, 4e-5, -7.143, -.9947, \
                        -.12, -.01007, -8.51e-4, -4e-5, -3.175, -1.116, -.227, \
                        -.01821, -.001729, 2.3e-4])


    kT = const.KBOLTZ*T
    gam = z0 * z0 * const.RYDBERG / kT
    gam1 = min(gam*1.e3,100.)
    if (gam > .1):
      gaunt_ff = kurucz(Evec/kT, gam)
    elif (kT==0.0):
      print "brems_kellog: Zero temperature!"
    else:
      u=Evec/kT
      gaunt_ff[(u>50) | (u==0)]=0.0
      i = numpy.where( (u<50)& (u !=0.0))[0]
      g=1.0
      u2=u[i]**2
      u1=0.5*u[i]
      t=u1/3.75
      u12=u1*0.5

      ak = numpy.zeros(len(i), dtype=float)

      ii = numpy.where(u12<=1.0)[0]
      t2=t[ii]**2
      t4=t2**2
      t8=t4**2
      ai = t2 * 3.5156229 + 1. + t4 * 3.089942 + t2 * 1.2067492 * t4 +\
           t8 * .2659732 + t8 * .0360768 * t2 + t8 * .0045813 * t4
      u122 = u12[ii] * u12[ii]

      ak[ii] = -numpy.log(u12[ii]) * ai - .57721566 + u122 *\
             (u122 * \
               (u122 * \
                 (u122 * \
                   (u122 * \
                     (u122 * 7.4e-6 + 1.075e-4) + \
                    .00262698) +\
                  .0348859) +\
                .23069756) +\
              .4227842)

      ii = numpy.where(u12>1.0)[0]
      type2=ii*1
      u12im = -1. / u12[ii]
      ak[ii] = u12im*\
                (u12im*\
                  (u12im*\
                    (u12im*\
                      (u12im*\
                        (u12im*5.3208e-4 + .0025154) +\
                       .00587872) +\
                     .01062446) +\
                   .02189568) +\
                 .07832358) + \
               1.25331414
      ak[ii] /= numpy.exp(u1[ii]) * numpy.sqrt(u1[ii])




      # see if born approx works
      born = numpy.exp(u1) * .5513 * ak

      if gam1<1.0:
        gaunt_ff[i] = born
        print "FORCE BORN"
      else:
        # go to do polynomial expansion

        u1 =u[i]
        u1[u1<0.003]=0.003
        u2=u1**2


        n = numpy.zeros(len(i),dtype=int)
        #m = numpy.zeros(len(i),dtype=int)

        n[(u1<=0.03)]=1
        n[(u1>0.03)&(u1<=0.3)]=2
        n[(u1>0.3)&(u1<=1.0)]=3
        n[(u1>1.0)&(u1<=5.0)]=4
        n[(u1>5.0)&(u1<=15.0)]=5
        n[(u1>15.0)]=6
        if (gam1<=1.773):
          m=1
        elif ((gam1>1.773)&(gam1<=3.0)):
          m=2
        elif ((gam1>3.0)&(gam1<=5.6234)):
          m=3
        elif ((gam1>5.6234)&(gam1<=10.0)):
          m=4
        elif ((gam1>10.0)&(gam1<=30.0)):
          m=5
        else:
          m=6

        m1=m+1

        g1= (a[n+(m+7)*6-49] +\
             a[n+(m+14)*6-49]*u1+\
             a[n+(m+21)*6-49]*u2)*born
        g2 = (a[n+(m1+7)*6-49] +\
              a[n+(m1+14)*6-49]*u1+\
              a[n+(m1+21)*6-49]*u2)*born
        p=(gam1-gam3[m-1])/gam2[m-1]
        gaunt_ff[i]= (1.0-p)*g1 + p*g2

    #aktype=numpy.ones(len(ak), dtype=int)
    #aktype[type2]=2
#    for i in range(len(E)):
#      print "E[%i]=%e ak type %i, m=%i, n=%i, u=%e, g_ff=%e"%\
#            (i,E[i], aktype[i],m,n[i], u[i],gaunt_ff[i])
#
#    zzz=raw_input()
    if not Eisvec:
      gaunt_ff = gaunt_ff[0]
    return gaunt_ff
  elif brems_type==const.RELATIVISTIC:

    #initialize
    rdat = atomdb.get_data(False, False,'rbrems', settings = settings,\
                                               datacache=datacache)
#    gaunt_D=rdat['BR_GAUNT'].data['COEFFICIENT']
    gaunt_U = rdat['GAUNT_FF'].data['LOG10_U']
    gaunt_Z = rdat['GAUNT_FF'].data['Z']
    gaunt_Ng = rdat['GAUNT_FF'].data['N_GAMMA2']
    gaunt_g2 = rdat['GAUNT_FF'].data['LOG10_GAMMA2']
    gaunt_gf = rdat['GAUNT_FF'].data['GAUNT_FF']


    gamma2 = numpy.log10(z0*z0*const.RYDBERG/(const.KBOLTZ*T))

    # extract the gaunt factors
    if z0 in gaunt_Z:
      # exact element is in file
      Uvec, GauntFFvec = extract_gauntff(z0, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
#      print "we are here"
#      for iuvec in xrange(len(Uvec)):
#        print "Uvec[%i]=%e, GauntFFvec[%i]=%e"%(iuvec, Uvec[iuvec],\
#                                                iuvec, GauntFFvec[iuvec])
      #zzz=raw_input()
    else:
      # find nearest elements in the file
      zlo = gaunt_Z[numpy.where(gaunt_Z < z0)[0]]
      if len(zlo) ==0:
        zlo=0
      else:
        zlo = max(zlo)
      zup = gaunt_Z[numpy.where(gaunt_Z > z0)[0]]
      if len(zup)==0:
        zup=100
      else:
        zup=min(zup)

      # check for various cases:
      if (zlo==0) & (zup<100):
        Uvec, GauntFFvec = extract_gauntff(zup, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
      if (zlo>0) & (zup==100):
        Uvec, GauntFFvec = extract_gauntff(zlo, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
      if (zlo==0) & (zup==100):
          #no match found
          print "brems relativistic: we should never be here"
      if (zlo>0) & (zup<100):
        # we are going to interpolate between the 2 of these
        Uvecl, GauntFFvecl = extract_gauntff(zlo, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
        Uvecu, GauntFFvecu = extract_gauntff(zup, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)

        if len(numpy.where(numpy.abs(Uvecl-Uvecu)>0.001)[0]) != 0 :
          print "Error: brems_relativistic: U vector mismatch ",  Uvecl, Uvecu

        Uvec = Uvecl
        GauntFFvec = ((zup-z0)*GauntFFvecl + (z0-zlo)*GauntFFvecu)/(zup-zlo)
#    print GauntFFvec
    # ok, so now we have Uvec and GauntFFvec assigned
    u = numpy.log10(Evec/(const.KBOLTZ*T))
    #for (iBin=0;iBin < N; iBin++) {

    # First, some special cases on u
    # Long wavelength.  Use approximation from Scheuer 1960, MNRAS, 120,231
    # As noted in Hummer, 1988, ApJ, 325, 477

    # low energy
    i = numpy.where(u<Uvec[0])[0]
    gaunt_ff[i]=-0.55133*(0.5*numpy.log(10.)*gamma2+numpy.log(10.)*u+0.056745)

    #high energy
    i = numpy.where(u>Uvec[-1])[0]
    gaunt_ff[i]=1.0

    #other energy
    i = numpy.where((u>=Uvec[0]) & (u <=Uvec[-1]))[0]
    gaunt_ff[i] = numpy.interp(u,Uvec, GauntFFvec)

    if not Eisvec:
      gaunt_ff = gaunt_ff[0]
    return gaunt_ff
  else:
    print "UNKNOWN BREMS TYPE: ", brems_type
    return -1


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


def extract_gauntff(Z, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf):
  """
  Extract the appropriate Gaunt free-free factor from the relativistic
  data tables of Nozawa, Itoh, & Kohyama, 1998 ApJ, 507,530

  Parameters
  ----------
  Z : int
    Z for which result is required
  gamma2 : array(float)
    gamma^2 in units of Z^2 Rydbergs/kT
  gaunt_U : array(float)
    u=E/kT
  gaunt_Z : array(int)
    nuclear charge
  gaunt_Ng : array(int)
    number of gamma^2 factors
  gaunt_g2 : array(float)
    gamma^2 factors
  gaunt_gf : array(float)
    ff factors

  Returns
  -------
  array(float)
    Gaunt factors.

  References
  ----------
  Nozawa, Itoh, & Kohyama, 1998 ApJ, 507,530
  """


    # look for one which matches the supplied Z
  i = numpy.where(gaunt_Z == Z)[0]
  Uvec = gaunt_U[i]
  GauntFFvec = numpy.zeros(len(i),dtype=float)

    # find the :"goo", interpolable values
  ii = numpy.where((gamma2>gaunt_g2[i,0]) &\
                   (gamma2<gaunt_g2[i,gaunt_Ng[i]-1]))[0]
#    print "i=",i, "len(i)=", len(i)
#    print "ii=",ii, "len(ii)=", len(ii)
#    print "i[ii]=",i[ii]
#    print "into the interpolation:"
  for iii in ii:
#      print "---"
#      print "iii=",iii
#      print "i[iii]=",i[iii]
#      print "gaunt_g2[i[iii]]", gaunt_g2[i[iii]]
#      print "gaunt_gf[i[iii]]", gaunt_gf[i[iii]]
#      print "gamma2", gamma2
#
#      print "len(GauntFFvec)", len(GauntFFvec)
#      print "gamma2 = %e, "%(gamma2)
    ng = gaunt_Ng[i[iii]]

#      for iv in xrange(len(gaunt_g2[i[iii]][:ng])):
#        print "gaunt_g2[%i][%i] = %e, gaunt_gf[%i][%i]=%e"%(i[iii],iv,gaunt_g2[i[iii]][iv],\
#                                                            i[iii],iv,gaunt_gf[i[iii]][iv])

      #, gaunt_g2[i[iii]], gaunt_gf[i[iii]]
    GauntFFvec[iii] = 10**(numpy.interp(gamma2,gaunt_g2[i[iii]][:ng], \
                                        numpy.log10(gaunt_gf[i[iii]][:ng])))
      #print "iii",iii, "GauntFFvec[iii]", GauntFFvec[iii]
    #zzz=raw_input()
  ii = numpy.where(gamma2<gaunt_g2[i,0])[0]
  GauntFFvec[ii]=gaunt_gf[i[ii],0]

  ii = numpy.where(gamma2>gaunt_g2[i,gaunt_Ng[i]-1])[0]
  GauntFFvec[ii]=gaunt_gf[i[ii],gaunt_Ng[i[ii]]-1]
  return Uvec, GauntFFvec



#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def do_brems(Z, z1, T, abund, brems_type, eedges):
  """
  Calculate the bremstrahlung emission in units of photon cm^3 s^-1 bin^-1

  Parameters
  ----------
  Z : int
    nuclear charge for which result is required
  z1 : int
    ion charge +1
  T : float
    temperture (Kelvin)
  abund : float
    elemental abundance (should be between 1.0 and 0.0)
  brems_type : int
    Type of bremstrahlung requested:
    1 = HUMMER = Non-relativistic: 1988ApJ...327..477H
    2 = KELLOGG = Semi-Relativistic: 1975ApJ...199..299K
    3 = RELATIVISTIC = Relativistic: 1998ApJ...507..530N
    4 = BREMS_NONE = no bremstrahlung
  eedges : array(float)
    The energy bin edges for the spectrum (keV)

  Returns
  -------
  array(float)
     bremstrahlung emission in units of photon cm^3 s^-1 bin^-1

  """

  kT = T*const.KBOLTZ # convert to keV

  E = (eedges[1:]+eedges[:-1])/2.0

  dE= (eedges[1:]-eedges[:-1])

  EkT=E/kT

  gaunt_ff = calc_brems_gaunt(E, T, z1, brems_type)

  emission = const.BREMS_COEFF*abund*\
      (z1-1)*(z1-1)*numpy.exp(-EkT)*(dE/numpy.sqrt(T))*gaunt_ff/(E*const.ERG_KEV)

  return emission




#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def kurucz(uin, gam):
  """
  Correction factors to Kellogg bremstrahlung calculation by Bob Kurucz

  Parameters
  ----------
  uin : array(float)
    energy grid, units of E/kT (both in keV)
  gam : array(float)
    Z**2/T, in units of Rydbergs

  Returns
  -------
  array(float)
    gaunt factors at high gam (> 0.1)
  """

  ya = numpy.array([5.4,5.25, 5.,4.69,4.48,4.16,3.85,4.77,4.63,4.4,\
                    4.13,3.87,3.52, 3.27,4.15,4.02,3.8,3.57,3.27,\
                    2.98,2.7,3.54,3.41,3.22,2.97,2.7,2.45,2.2,2.94,\
                    2.81,2.65,2.44,2.21,2.01,1.81,2.41,2.32,2.19,2.02,\
                    1.84,1.67,1.5,1.95,1.9,1.8,1.68,1.52,1.41,1.3,\
                    1.55,1.56,1.51,1.42,1.33,1.25,1.17,1.17,1.3,1.32,\
                    1.3,1.2,1.15,1.11,.86,1.,1.15,1.18,1.15,1.11,1.08,\
                    .59,.76,.97,1.09,1.13,1.1,1.08,.38,.53,.76,.96,\
                    1.08,1.09,1.09])
  gaunt = numpy.ones(len(uin), dtype=float)


  rj = numpy.log10(gam) * 2. + 3.
  j = numpy.array(numpy.floor(rj), dtype=int)
  rj =  numpy.array(j, dtype=float)

  rk = numpy.log10(uin) * 2. + 9.
  k = numpy.array(numpy.floor(rk), dtype=int)
  k[k<1] = 1
  rk =  numpy.array(k, dtype=float)



  gaunt[(j >= 7) | (j < 1) | (k>=12)]=0.0

  i = numpy.where(gaunt>0)[0]


  t = (numpy.log10(gam) - (rj - 3.) / 2.) / .5
  u = (numpy.log10(uin[i]) - (rk[i] - 9.) / 2.) / .5
  # so t is a scalar
  # u is a vector

  gaunt[i] = (1. - t) * (1. - u) * ya[j + k[i] * 7 - 8] + \
    t * (1.-u)* ya[j + 1 + k[i]*7 - 8] + t*u*ya[j + 1 + (k[i] + 1) * 7 - 8] +\
    (1. - t) * u * ya[j + (k[i] + 1) * 7 - 8]

  return gaunt

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


def calc_ee_brems(E, T, N):
  """
  calculate the electron-electron bremsstrahlung.

  Parameters
  ----------
  E : array (float)
    energy grid (keV)
  T : float
    Electron temperature (keV)
  N : float
    electron density (cm^-3)

  Returns
  -------
  array(float)
    ee_brems in photons cm^s s^-1 keV-1 at each point E.
    This should be multiplied by the bin width to get flux per bin.

  References
  ----------
  Need to check this!
  """
#
#  T is the electron temperature (in keV)
#  N is the electron density (in cm^-3)
#  E is the list of energies (in keV)
#

# series of data constants
# Region I, k_BT<=1 keV
  #print "T=%f"%(T)

  aI1 = numpy.array([(3.15847E+0, -2.52430E+0, 4.04877E-1, 6.13466E-1, 6.28867E-1, 3.29441E-1),
             (2.46819E-2, 1.03924E-1, 1.98935E-1, 2.18843E-1, 1.20482E-1, -4.82390E-2),
             (-2.11118E-2, -8.53821E-2, -1.52444E-1, -1.45660E-1, -4.63705E-2, 8.16592E-2),
             (1.24009E-2, 4.73623E-2, 7.51656E-2, 5.07201E-2, -2.25247E-2, -8.17151E-2),
             (-5.41633E-3, -1.91406E-2, -2.58034E-2, -2.23048E-3, 5.07325E-2, 5.94414E-2),
             (1.70070E-3, 5.39773E-3, 4.13361E-3, -1.14273E-2, -3.23280E-2, -2.19399E-2),
             (-3.05111E-4, -7.26681E-4, 4.67015E-3, 1.24789E-2, -1.16976E-2, -1.13488E-2),
             (-1.21721E-4, -7.47266E-4, -2.20675E-3, -2.74351E-3, -1.00402E-3, -2.38863E-3),
             (1.77611E-4, 8.73517E-4, -2.67582E-3, -4.57871E-3, 2.96622E-2, 1.89850E-2),
             (-2.05480E-5, -6.92284E-5, 2.95254E-5, -1.70374E-4, -5.43191E-4, 2.50978E-3),
             (-3.58754E-5, -1.80305E-4, 1.40751E-3, 2.06757E-3, -1.23098E-2, -8.81767E-3)])

# column j=0-5; line i=0-10
  aI2 = numpy.array([(-1.71486E-1, -3.68685E-1, -7.59200E-2, 1.60187E-1, 8.37729E-2),
             (-1.20811E-1, -4.46133E-4, 8.88749E-2, 2.50320E-2, -1.28900E-2),
             (9.87296E-2, -3.24743E-2, -8.82637E-2, -7.52221E-3, 1.99419E-2),
             (-4.59297E-2, 5.05096E-2, 5.58818E-2, -9.11885E-3, -1.71348E-2),
             (-2.11247E-2, -5.05387E-2, 9.20453E-3, 1.67321E-2, -3.47663E-3),
             (1.76310E-2, 2.23352E-2, -4.59817E-3, -8.24286E-3, -3.90032E-4),
             (6.31446E-2, 1.33830E-2, -8.54735E-2, -6.47349E-3, 3.72266E-2),
             (-2.28987E-3, 7.79323E-3, 7.98332E-3, -3.80435E-3, -4.25035E-3),
             (-8.84093E-2, -2.93629E-2, 1.02966E-1, 1.38957E-2, -4.22093E-2),
             (4.45570E-3, -2.80083E-3, -5.68093E-3, 1.10618E-3, 2.33625E-3),
             (3.46210E-2, 1.23727E-2, -4.04801E-2, -5.68689E-3, 1.66733E-2)])
# column j=6-10, line 0-10

# Region II (1 keV<=k_B<=300 keV)
  abII = numpy.array([(-3.7369800E+1, -9.3647000E+0, 9.2170000E-1, -1.1628100E+1, -8.6991000E+0),
              (3.8036590E+2, 9.5918600E+1, -1.3498800E+1, 1.2560660E+2, 6.3383000E+1),
              (-1.4898014E+3, -3.9701720E+2, 7.6453900E+1, -5.3274890E+2, -1.2889390E+2),
              (2.8614150E+3, 8.4293760E+2, -2.1783010E+2, 1.1423873E+3, -1.3503120E+2),
              (-2.3263704E+3, -9.0730760E+2, 3.2097530E+2, -1.1568545E+3, 9.7758380E+2),
              (-6.9161180E+2, 3.0688020E+2, -1.8806670E+2, 7.5010200E+1, -1.6499529E+3),
              (2.8537893E+3, 2.9129830E+2, -8.2416100E+1, 9.9681140E+2, 1.2586812E+3),
              (-2.0407952E+3, -2.9902530E+2, 1.6371910E+2, -8.8818950E+2, -4.0474610E+2),
              (4.9259810E+2, 7.6346100E+1, -6.0024800E+1, 2.5013860E+2, 2.7335400E+1)])
# column a_0j-a2j, b0j,b1j; line j=0-8

  cII = numpy.array([(-5.7752000E+0, 3.0558600E+1, -5.4327200E+1, 3.6262500E+1, -8.4082000E+0),
             (4.6209700E+1, -2.4821770E+2, 4.5096760E+2, -3.1009720E+2, 7.4792500E+1),
             (-1.6072800E+2, 8.7419640E+2, -1.6165987E+3, 1.1380531E+3, -2.8295400E+2),
             (3.0500700E+2, -1.6769028E+3, 3.1481061E+3, -2.2608347E+3, 5.7639300E+2),
             (-3.2954200E+2, 1.8288677E+3, -3.4783930E+3, 2.5419361E+3, -6.6193900E+2),
             (1.9107700E+2, -1.0689366E+3, 2.0556693E+3, -1.5252058E+3, 4.0429300E+2),
             (-4.6271800E+1, 2.6056560E+2, -5.0567890E+2, 3.8008520E+2, -1.0223300E+2)])
# column CII_2j-cII_6j; line j=0-6

# Region III (300 keV<= k_BT<=7 MeV)
  abIII = numpy.array([(5.2163300E+1, 4.9713900E+1, 6.4751200E+1, -8.5862000E+0, 3.7643220E+2),
               (-2.5703130E+2, -1.8977460E+2, -2.1389560E+2, 3.4134800E+1, -1.2233635E+3),
               (4.4681610E+2, 2.7102980E+2, 1.7414320E+2, -1.1632870E+2, 6.2867870E+2),
               (-2.9305850E+2, -2.6978070E+2, 1.3650880E+2, 2.9654510E+2, 2.2373946E+3),
               (0.0000000E+0, 4.2048120E+2, -2.7148990E+2, -3.9342070E+2, -3.8288387E+3),
               (7.7047400E+1, -5.7662470E+2, 8.9321000E+1, 2.3754970E+2, 2.1217933E+3),
               (-2.3871800E+1, 4.3277900E+2, 5.8258400E+1, -3.0600000E+1, -5.5166700E+1),
               (0.0000000E+0, -1.6053650E+2, -4.6080700E+1, -2.7617000E+1, -3.4943210E+2),
               (1.9970000E-1, 2.3392500E+1, 8.7301000E+0, 8.8453000E+0, 9.2205900E+1)])
# column aII_0j-aIII_2j, bIII_0j, bIII_1j; line j=0-8
  #print "E:",E
  aI = numpy.hstack((aI1,aI2))

  inum1 = numpy.arange(11)
  jnum = numpy.resize(inum1,(11,11))
  inum = numpy.transpose(jnum)
  aII = numpy.zeros((9,3))
  bII = numpy.zeros((9,2))
  [aII, bII] = numpy.hsplit(abII,numpy.array([3]))

  numII = numpy.arange(9)
  aIIj = numpy.transpose(numpy.resize(numII,(3,9)))
  aIIi = numpy.resize(numpy.arange(3),(9,3))
  bIIj = numpy.transpose(numpy.resize(numpy.arange(9),(2,9)))
  bIIi = numpy.resize(numpy.arange(2),(9,2))
  cIIj = numpy.transpose(numpy.resize(numpy.arange(7),(5,7)))
  cIIi = numpy.resize(numpy.arange(2,7),(7,5))

  aIII = numpy.zeros((9,3))
  bIII = numpy.zeros((9,2))

  kTIII = 1.e3 #1 MeV, in unit of keV
  [aIII, bIII] = numpy.hsplit(abIII,numpy.array([3]))


  tao = T/const.ME_KEV
  # find the length of the input
  if isinstance(E, (collections.Sequence, numpy.ndarray)):
    x = E/T
    (numx,) = x.shape
  else:
    Earray = numpy.array([E])
    x = Earray/T
    numx=1
  #print "x", x
  #print "numx:", numx
  GI = numpy.zeros((numx,))
  AIIr = numpy.zeros((numx,))
  BIIr = numpy.zeros((numx,))
  GpwII = numpy.zeros((numx,))
  GII = numpy.zeros((numx,))
  FCCII = numpy.zeros((numx,))
  Ei0 = numpy.zeros((numx,))
  GpwIII = numpy.zeros((numx,))
  if T<0.05:
    ret = numpy.zeros(len(x), dtype=float)
  elif 0.05<=T<70.:
    # hmm
    # ARF REDO
    GI=numpy.zeros(len(x))
    theta = (1/1.35) * ( numpy.log10(tao) + 2.65)
    bigx = (1/2.5) * (numpy.log10(x) + 1.50)
    for i in range(11):
      for j in range(11):
         GI += aI[i,j]*(theta**i)*(bigx**j)
    GI *= numpy.sqrt(8/(3*numpy.pi))
    ret = 1.455e-16*N**2*numpy.exp(-x)/(x*numpy.sqrt(tao))*GI
  elif 70.<=T<300.:
    taoII = tao
    for k in range(numx):
      def integrand(t):
        return numpy.exp(-1.0*t)/t
      [Ei0[k,],error] = scipy.integrate.quad(integrand,x[k,],\
                                             numpy.Inf,args=())
      AIIr[k,] = numpy.sum(aII*taoII**(aIIj/8.)*x[k,]**(aIIi))
      BIIr[k,] = numpy.sum(bII*taoII**(bIIj/8.)*x[k,]**(bIIi))
      FCCII[k,] = 1.+numpy.sum(cII*taoII**(cIIj/6.)*x[k,]**(cIIi/8.))
      GpwII[k,] = numpy.sum(aII*taoII**(aIIj/8.)*x[k,]**(aIIi))-\
                      numpy.exp(x[k,])*(-1.0)*Ei0[k,]*\
                      numpy.sum(bII*taoII**(bIIj/8.)*x[k,]**(bIIi))
      GII[k,] = GpwII[k,]*FCCII[k,]
    ret = 1.455e-16*N**2*numpy.exp(-x)/(x*numpy.sqrt(taoII))*GII
  elif 300.<=T<7000.:
    taoIII = taoII
    for k in range(numx):
      GpwIII[k,] = numpy.sum(aIII*taoIII**(aIIj/8.)*x[k,]**(aIIi))-\
                   numpy.exp(x[k,])*(-1.0)*Ei0[k,]*\
                   numpy.sum(bIII*taoIII**(bIIj/8.)*x[k,]**(bIIi))
    ret = 1.455e-16*N**2*numpy.exp(-x)/\
           (x*numpy.sqrt(taoIII))*GpwIII
  else:
    taoIV = tao
    GIV = 3./(4.*numpy.pi*numpy.sqrt(taoIV))*\
          (28./3.+2.*x+x**2/2.+2.*(8./3.+4.*x/3.+x**2)*\
          (numpy.log(2.*taoIV)-0.57721)-numpy.exp(x) \
          *(-1.0*Ei0)*(8./3.-4.*x/3.+x**2))

    ret = 1.455e-16*N**2*numpy.exp(-x)/(x*numpy.sqrt(taoIV))*GIV
  return ret/const.ME_KEV


def parse_par_file(fname):
  """
  Parse the apec.par input file for controlling APEC

  Parameters
  ----------
  fname : string
    file name

  Returns
  -------
  dict
    The settings in "key:value" pairs.
  """


  # Adam Foster 2016-07-27

  # check the par file exists
  f = open(fname,'r')

  data = {}
  rawdata = {}
  for line in f:
    # check for comment lines
    if line[0] == '#': continue

    lall = line.split('"')
    
    l = []
    for ill, ll in enumerate(lall):
      if ill % 2 == 0:
        # even, not a quote
        ltmp = ll.split(',')
        if ill > 0:
          ltmp=ltmp[1:]
        for iltmp in ltmp:
          l.append(iltmp)
      else:
        l=l[:-1]
        l.append(ll)
          

        
        

      
        
        
        
        
    name = l[0]
    dtype = l[1]
    hidden = l[2]
    value = l[3]
    minval = l[4]
    maxval= l[5]
    descr = l[6]

    rawdata[name]={}
    rawdata[name]['dtype']=dtype
    rawdata[name]['hidden']=hidden
    rawdata[name]['value']=value
    rawdata[name]['minval']=minval
    rawdata[name]['maxval']=maxval
    rawdata[name]['descr']=descr

    # some checking
    if dtype=='s':
      value = value.strip('"')
      minval = minval.strip('"')

      #string
      if '|' in minval:
        # check the value is one of the options
        options = minval.split('|')

        if not value in options:
          print "Error in paramater file %s. Item %s: %s is not in allowed parameter list %s"%\
                (fname, name, value, minval)
          print options
        else:
          data[name] = value
      else:
        data[name] = value
    elif dtype =='b':
      #boolean, yes or no.
      if not value in ['yes','no']:
        print "Error in paramater file %s. Item %s: should be boolean yes or no, supplied value is %s"%\
              (fname, name, value)
      else:
        if value=='yes':
          data[name]=True
        else:
          data[name]=False
    elif dtype=='r':
      # real
      value = float(value)

      if len(minval) > 0:
        minval = float(minval)
        if value < minval:
          print "Error in paramater file %s. Item %s: should be > %e, supplied value is %e"%\
              (fname, name, minval, value)

      if len(maxval) > 0:
        maxval = float(maxval)
        if value > maxval:
          print "Error in paramater file %s. Item %s: should be < %e, supplied value is %e"%\
              (fname, name, maxval, value)

      data[name] = value
    elif dtype=='i':
      # integer
      value = int(value)
      if len(minval) > 0:
        minval = int(minval)
        if value < minval:
          print "Error in paramater file %s. Item %s: should be > %i, supplied value is %i"%\
              (fname, name, minval, value)

      if len(maxval) > 0:
        maxval = int(maxval)
        if value > maxval:
          print "Error in paramater file %s. Item %s: should be < %i, supplied value is %i"%\
              (fname, name, maxval, value)
      data[name] = value
    else:
      print "Error in paramater file %s. Item %s: unknown data type %s"%\
          (fname, name, dtype)
  
    # now some massaging of the parameters
  if 'BremsType' in data.keys():
    if data['BremsType']=='Hummer':
      data['BremsType']=const.HUMMER
    elif data['BremsType']=='Kellogg':
      data['BremsType']=const.KELLOGG
    elif data['BremsType']=='Relativistic':
      data['BremsType']=const.RELATIVISTIC
    else:
      print "UNKNOWN BREMS TYPE: %s" %(data['BremsType'])
  
  if 'NEIMinEpsilon' in data.keys():
      tmp = data['NEIMinEpsilon'].split(',')
      tmp2 = []
      for i in tmp:
        tmp2.append(float(i))
      data['NEIMinEpsilon']=tmp2

  if 'NEIMinFrac' in data.keys():
      tmp = data['NEIMinFrac'].split(',')
      tmp2 = []
      for i in tmp:
        tmp2.append(float(i))
      data['NEIMinFrac']=tmp2

  if 'IncAtoms' in data.keys():
    elsymblist = data['IncAtoms'].split(',')
    Zlist=[]
    for iel in elsymblist:
      Zlist.append(atomic.elsymb_to_Z(iel))
    data['Zlist']=Zlist
    
  if not ('WriteIon' in data.keys()):
    data['WriteIon'] = False

        
  return data


def make_vector(linear, minval, step, nstep):
  """
  Create a vector from the given inputs

  Parameters
  ----------

  linear: boolean
    Whether the array should be linear or log spaced

  minval: float
    initial value of the array. In dex if linear==False

  step: float
    step between points on the array. In dex if linear==False

  nstep: int
    number of steps

  Returns
  -------

  array(float)
    array of values spaced out use the above parameters

  """

  arr = numpy.arange(nstep)*step+minval

  if linear==False:
    arr=10**arr

  return arr


def make_vector_nbins(linear, minval, maxval, nstep):
  """
  Create a vector from the given inputs

  Parameters
  ----------

  linear: boolean
    Whether the array should be linear or log spaced

  minval: float
    initial value of the array. In dex if linear==False

  maxval: float
    maximum value of the array. In dex if linear==False

  nstep: int
    number of steps

  Returns
  -------

  array(float)
    array of values spaced out use the above parameters

  """

  arr = numpy.linspace(minval, maxval, nstep+1)

  if linear==False:
    arr=10**arr

  return arr


def run_apec(fname):
  """
  Run the entire APEC code using the data in the parameter file fname

  Parameters
  ----------
  fname : string
    file name

  Returns
  -------
  None

  """

  # get the input settings
  settings = parse_par_file(fname)

  # need to parse the IncAtoms parameter - this defines which atoms
  # we will include
  #
  # we will transfer this to a "Zlist" parameter, in the form of
  # nuclear charges

  Zlist = settings['Zlist']
  print "I will be running Z=", Zlist
  # run for each element, temperature, density
  
  lhdulist = []
  chdulist = []
  
  seclHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=generate_datatypes('lineparams'))
  seccHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=generate_datatypes('cocoparams'))
  
  
  for iTe in range(settings['NumTemp']):

    te = make_vector(settings['LinearTemp'], \
                     settings['TempStart'], \
                     settings['TempStep'], \
                     settings['NumTemp'])[iTe]

    if settings['TempUnits']=='keV':
      te /= const.KBOLTZ


    for iDens in range(settings['NumDens']):
      dens = make_vector(settings['LinearDens'], \
                         settings['DensStart'], \
                         settings['DensStep'], \
                         settings['NumDens'])[iDens]

      # AT THIS POINT, GENERATE SHELL WHICH WILL GO IN THE HDU OF CHOICE
      
      if settings['Ionization']=='CIE':
        linedata = numpy.zeros(0,dtype=generate_datatypes('linelist_cie'))
        cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))
      if settings['Ionization']=='NEI':
        linedata = numpy.zeros(0,dtype=generate_datatypes('linetype'))
        cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))

      for Z in Zlist:
        print "Calling run_apec_element for Z=%i Te=%e dens=%e at %s"%(Z, te, dens, time.asctime())
        dat = run_apec_element(settings, te, dens, Z)
        # append this data to the output
        linedata = numpy.append(linedata, dat['lines'])
        cocodata = continuum_append(cocodata, dat['cont'])


      
      # now make an HDU for all of this
      if settings['Ionization']=='CIE':
        LHDUdat = create_lhdu_cie(linedata)
      elif settings['Ionization']=='NEI':
        LHDUdat = create_lhdu_nei(linedata)
      # now update the headers
      iseclHDUdat=iDens+iTe*settings['NumDens']
      LHDUdat.header['EXTNAME']=("EMISSIVITY","name of this binary table extension")
      LHDUdat.header['EXTVER']=(iseclHDUdat+1,"Index for this EMISSIVITY extension")

      LHDUdat.header['HDUNAME'] = ("L%.2f_%.2f"%(numpy.log10(te), numpy.log10(dens)),\
                             'Spectral emission data')
      LHDUdat.header['HDUCLASS'] = ("Proposed OGIP",\
                             'Proposed OGIP standard')
      LHDUdat.header['HDUCLAS1']=("LINE MODEL",\
                             'Line emission spectral model')
      LHDUdat.header['HDUCLAS2']=("LINE",\
                             'Emission line data')
      if settings['Ionization']=='CIE':
        LHDUdat.header['HDUVERS1']=("1.0.0",\
                               'version of format')
      elif settings['Ionization']=='NEI':
        LHDUdat.header['HDUVERS1']=("2.0.0",\
                               'version of format')

      LHDUdat.header['TEMPERATURE']=(te,\
                             'Electron temperature')
      LHDUdat.header['DENSITY']=(dens,\
                             'Electron density')
      LHDUdat.header['TIME']=(0,\
                             'IN EQUILIBRIUM')
      if settings['Ionization']=='CIE':
        tot_emiss = sum(linedata['epsilon']*const.HC_IN_ERG_A/linedata['lambda'])
      else:
        tot_emiss=0.0
      LHDUdat.header['TOT_LINE']=(tot_emiss,\
                             'Total Line Emission (erg cm^3 s^-1)')
      LHDUdat.header['N_LINES']=(len(linedata),\
                             'Number of emission lines')

      seclHDUdat['kT'][iseclHDUdat]=te*const.KBOLTZ
      seclHDUdat['EDensity'][iseclHDUdat]= dens
      seclHDUdat['Time'][iseclHDUdat]= 0.0
      seclHDUdat['Nelement'][iseclHDUdat]= len(Zlist)
      seclHDUdat['Nline'][iseclHDUdat]= len(linedata)
      
      lhdulist.append(LHDUdat)


# continuum data
      # now make an HDU for all of this
      CHDUdat = create_chdu_cie(cocodata)
      
      # now update the headers
      iseccHDUdat=iDens+iTe*settings['NumDens']

      CHDUdat.header['EXTNAME']=("EMISSIVITY","name of this binary table extension")
      CHDUdat.header['EXTVER']=(iseccHDUdat+1,"Index for this EMISSIVITY extension")

      CHDUdat.header['HDUNAME'] = ("C%.2f_%.2f"%(numpy.log10(te), numpy.log10(dens)),\
                             'Spectral emission data')
      CHDUdat.header['HDUCLASS'] = ("Proposed OGIP",\
                             'Proposed OGIP standard')
      CHDUdat.header['HDUCLAS1']=("COMP CONT MODEL",\
                             'Compressed continua spectra')
      CHDUdat.header['HDUCLAS2']=("COCO",\
                             'Compressed continuum data')
      CHDUdat.header['HDUVERS1']=("1.0.0",\
                             'version of format')
      CHDUdat.header['TEMPERATURE']=(te,\
                             'Electron temperature')
      CHDUdat.header['DENSITY']=(dens,\
                             'Electron density')
      CHDUdat.header['TIME']=("%.2e"%(0),\
                             'IN EQUILIBRIUM')
      tot_emiss = calc_total_coco(cocodata, settings)
      
      if settings['Ionization']=='CIE':
        CHDUdat.header['TOT_COCO']=(tot_emiss,\
                               'Total Emission (erg cm^3 s^-1)')
      else:
        CHDUdat.header['TOT_COCO']=(0.0,\
                               'Total Emission (erg cm^3 s^-1)')
        
      seccHDUdat['kT'][iseccHDUdat]=te*const.KBOLTZ
      seccHDUdat['EDensity'][iseccHDUdat]= dens
      seccHDUdat['Time'][iseccHDUdat]= 0.0
      seccHDUdat['NElement'][iseccHDUdat]= len(cocodata)
      seccHDUdat['NCont'][iseccHDUdat]= max(cocodata['N_Cont'])
      seccHDUdat['NPseudo'][iseccHDUdat]= max(cocodata['N_Pseudo'])
      
      chdulist.append(CHDUdat)










    # make secHDUdat into a fits table
  seclHDU = create_lparamhdu_cie(seclHDUdat)
  seclHDU.header['EXTNAME']=('PARAMETERS','name of this binary table extension')
  seclHDU.header['HDUCLASS']=('Proposed OGIP','Proposed OGIP standard')
  seclHDU.header['HDUCLAS1']=('LINE MODEL','line emission spectra model')
  seclHDU.header['HDUCLAS2']=('Parameters','extension containing parameter info')
  seclHDU.header['HDUVERS1']=('1.0.0','version of format')


  seccHDU = create_cparamhdu_cie(seccHDUdat)
  seccHDU.header['EXTNAME']=('PARAMETERS','name of this binary table extension')
  seccHDU.header['HDUCLASS']=('Proposed OGIP','Proposed OGIP standard')
  seccHDU.header['HDUCLAS1']=('COMP CONT MODEL','compressed continua spectra')
  seccHDU.header['HDUCLAS2']=('Parameters','extension containing parameter info')
  seccHDU.header['HDUVERS1']=('1.0.0','version of format')


  fileroot = settings['OutputFileStem']
  if settings['Ionization']=='NEI':
    fileroot+='_nei'
  # create the Primary HDU
  PrilHDU = pyfits.PrimaryHDU()
  lhdulist.insert(0,PrilHDU)
  lhdulist.insert(1,seclHDU)
  tmplhdulist = pyfits.HDUList(lhdulist)
  tmplhdulist.writeto('%s_line.fits'%(fileroot), clobber=True, checksum=True)


  PricHDU = pyfits.PrimaryHDU()
  chdulist.insert(0,PricHDU)
  chdulist.insert(1,seccHDU)
  tmpchdulist = pyfits.HDUList(chdulist)
  if settings['Ionization']=='CIE':
    tmpchdulist.writeto('%s_coco.fits'%(fileroot), clobber=True, checksum=True)
  elif settings['Ionization']=='NEI':
    tmpchdulist.writeto('%s_comp.fits'%(fileroot), clobber=True, checksum=True)

def calc_total_coco(cocodata, settings):
  """
  Calculate the total emission in erg cm^3 s^-1
  
  """
  emiss = 0.0
  
  ebins = make_vector_nbins(settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])

  ecent = (ebins[1:]+ebins[:-1])/2
  
  width = ebins[1:]-ebins[:-1]
  
  for i in range(len(cocodata)):
    E_Cont = cocodata[i]['E_Cont'][:cocodata[i]['N_Cont']]
    Cont = cocodata[i]['Continuum'][:cocodata[i]['N_Cont']]
    Cont_erg = Cont*E_Cont*const.ERG_KEV
    
    e = numpy.interp(ecent, E_Cont, Cont_erg)
    emiss += sum(e*width)

    E_Cont = cocodata[i]['E_Pseudo'][:cocodata[i]['N_Pseudo']]
    Cont = cocodata[i]['Pseudo'][:cocodata[i]['N_Pseudo']]
    Cont_erg = Cont*E_Cont*const.ERG_KEV
    
    e = numpy.interp(ecent, E_Cont, Cont_erg)
    emiss += sum(e*width)
    

  return emiss

def continuum_append(a,b):
  """
  Join two continuum arrays together, expanding arrays as necessary
  
  Parameters
  ----------
  a: numpy.array(dtype=continuum)
    The first array
  b: numpy.array(dtype=continuum)
    The second array
  
  Returns
  c: numpy.array(dtype=continuum)
    The two arrays combined, with continuum arrays resized as required.
  """
  if len(a) == 0:
    npseudomax = max(b['N_Pseudo'])
    ncontmax = max(b['N_Cont'])
    nlines = len(b)
  
  else:
    npseudomax = max([max(a['N_Pseudo']), max(b['N_Pseudo'])])
    ncontmax = max([max(a['N_Cont']), max(b['N_Cont'])])
    nlines = len(a) + len(b)
  
  c = numpy.zeros(nlines, dtype=generate_datatypes('continuum',npseudo=npseudomax, ncontinuum=ncontmax))
  ic = 0
  if len(a) > 0:
    for ia in range(len(a)):
      c['Z'][ic] = a['Z'][ia]
      c['rmJ'][ic] = a['rmJ'][ia]
      c['N_Cont'][ic] = a['N_Cont'][ia]
      c['E_Cont'][ic][:c['N_Cont'][ic]] = a['E_Cont'][ia][:a['N_Cont'][ia]]
      c['Continuum'][ic][:c['N_Cont'][ic]] = a['Continuum'][ia][:a['N_Cont'][ia]]
      c['N_Pseudo'][ic] = a['N_Pseudo'][ia]
      c['E_Pseudo'][ic][:c['N_Pseudo'][ic]] = a['E_Pseudo'][ia][:a['N_Pseudo'][ia]]
      c['Pseudo'][ic][:c['N_Pseudo'][ic]] = a['Pseudo'][ia][:a['N_Pseudo'][ia]]
      ic +=1
  for ib in range(len(b)):
    c['Z'][ic] = b['Z'][ib]
    c['rmJ'][ic] = b['rmJ'][ib]
    c['N_Cont'][ic] = b['N_Cont'][ib]
    c['E_Cont'][ic][:c['N_Cont'][ic]] = b['E_Cont'][ib][:b['N_Cont'][ib]]
    c['Continuum'][ic][:c['N_Cont'][ic]] = b['Continuum'][ib][:b['N_Cont'][ib]]
    c['N_Pseudo'][ic] = b['N_Pseudo'][ib]
    c['E_Pseudo'][ic][:c['N_Pseudo'][ic]] = b['E_Pseudo'][ib][:b['N_Pseudo'][ib]]
    c['Pseudo'][ic][:c['N_Pseudo'][ic]] = b['Pseudo'][ib][:b['N_Pseudo'][ib]]
    ic +=1
  
  return c

def create_lhdu_cie(linedata):

  # sort the data
  linedata.sort(order=['lambda'])
  linedata = linedata[::-1]
  linedata.sort(order=['element','ion'], kind='mergesort')
  

  
  cols = []
  cols.append(pyfits.Column(name='Lambda', format='1E', unit="A", array=linedata['lambda']))
  cols.append(pyfits.Column(name='Lambda_Err', format='1E', unit="A", array=linedata['lambda_err']))
  cols.append(pyfits.Column(name='Epsilon', format='1E', unit="photons cm^3 s^-1", array=linedata['epsilon']))
  cols.append(pyfits.Column(name='Epsilon_Err', format='1E', unit="photons cm^3 s^-1", array=linedata['epsilon_err']))
  cols.append(pyfits.Column(name='Element', format='1J',  array=linedata['element']))
  cols.append(pyfits.Column(name='Ion', format='1J',  array=linedata['ion']))
  cols.append(pyfits.Column(name='UpperLev', format='1J', array=linedata['upperlev']))
  cols.append(pyfits.Column(name='LowerLev', format='1J',  array=linedata['lowerlev']))
  
  coldefs = pyfits.ColDefs(cols)
  tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
  print "Created a linelist HDU with %i lines"%(len(tbhdu.data))
  return tbhdu

def create_lhdu_nei(linedata):

  # sort the data
  linedata.sort(order=['lambda'])
  linedata = linedata[::-1]
  linedata.sort(order=['element','ion'], kind='mergesort')

  
  cols = []
  cols.append(pyfits.Column(name='Lambda', format='1E', unit="A", array=linedata['lambda']))
  cols.append(pyfits.Column(name='Lambda_Err', format='1E', unit="A", array=linedata['lambda_err']))
  cols.append(pyfits.Column(name='Epsilon', format='1E', unit="photons cm^3 s^-1", array=linedata['epsilon']))
  cols.append(pyfits.Column(name='Epsilon_Err', format='1E', unit="photons cm^3 s^-1", array=linedata['epsilon_err']))
  cols.append(pyfits.Column(name='Element', format='1J',  array=linedata['element']))
  cols.append(pyfits.Column(name='Element_drv', format='1J',  array=linedata['element']))
  cols.append(pyfits.Column(name='Ion', format='1J',  array=linedata['ion']))
  cols.append(pyfits.Column(name='Ion_drv', format='1J',  array=linedata['ion_drv']))
  cols.append(pyfits.Column(name='UpperLev', format='1J', array=linedata['upperlev']))
  cols.append(pyfits.Column(name='LowerLev', format='1J',  array=linedata['lowerlev']))

  coldefs = pyfits.ColDefs(cols)
  tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
  return tbhdu


def create_chdu_cie(cocodata):
  
  ncont = max(cocodata['N_Cont'])
  npseudo = max(cocodata['N_Pseudo'])
  
  cols = []
  cols.append(pyfits.Column(name='Z', format='1J', array=cocodata['Z']))
  cols.append(pyfits.Column(name='rmJ', format='1J', array=cocodata['rmJ']))
  cols.append(pyfits.Column(name='N_Cont', format='1J', array=cocodata['N_Cont']))
  cols.append(pyfits.Column(name='E_Cont', format='%iE'%(ncont), unit='keV', array=cocodata['E_Cont']))
  cols.append(pyfits.Column(name='Continuum', format='%iE'%(ncont), unit='photons cm^3 s^-1 keV^-1', array=cocodata['Continuum']))
  cols.append(pyfits.Column(name='Cont_Err', format='%iE'%(ncont), unit='photons cm^3 s^-1 keV^-1', array=cocodata['Cont_Err']))
  cols.append(pyfits.Column(name='N_Pseudo', format='1J', array=cocodata['N_Pseudo']))
  cols.append(pyfits.Column(name='E_Pseudo', format='%iE'%(npseudo), unit='keV', array=cocodata['E_Pseudo']))
  cols.append(pyfits.Column(name='Pseudo', format='%iE'%(npseudo), unit='photons cm^3 s^-1 keV^-1', array=cocodata['Pseudo']))
  cols.append(pyfits.Column(name='Pseudo_Err', format='%iE'%(npseudo), unit='photons cm^3 s^-1 keV^-1', array=cocodata['Pseudo_Err']))
  
  coldefs = pyfits.ColDefs(cols)
  tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
  return tbhdu



def create_lparamhdu_cie(linedata):
  
  cols = []
  cols.append(pyfits.Column(name='kT', format='1E', unit="keV", array=linedata['kT']))
  cols.append(pyfits.Column(name='EDensity', format='1E', unit="cm**-3", array=linedata['EDensity']))
  cols.append(pyfits.Column(name='Time', format='1E', unit="s", array=linedata['Time']))
  cols.append(pyfits.Column(name='Nelement', format='1J', array=linedata['Nelement']))
  cols.append(pyfits.Column(name='Nline', format='1J',  array=linedata['Nline']))
  
  coldefs = pyfits.ColDefs(cols)
  tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
  return tbhdu

def create_cparamhdu_cie(cocodata):
  cols = []
  cols.append(pyfits.Column(name='kT', format='1E', unit="keV", array=cocodata['kT']))
  cols.append(pyfits.Column(name='EDensity', format='1E', unit="cm**-3", array=cocodata['EDensity']))
  cols.append(pyfits.Column(name='Time', format='1E', unit="s", array=cocodata['Time']))
  cols.append(pyfits.Column(name='NElement', format='1J', array=cocodata['NElement']))
  cols.append(pyfits.Column(name='NCont', format='1J', array=cocodata['NCont']))
  cols.append(pyfits.Column(name='NPseudo', format='1J', array=cocodata['NPseudo']))

  coldefs = pyfits.ColDefs(cols)
  tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
  return tbhdu

  
  
def run_apec_element(settings, te, dens, Z):
  """
  Run the APEC code using the settings provided for one element

  Parameters
  ----------
  settings: dictionary
    The settings read from the apec.par file by parse_par_file

  te: float
    The electron temperature (K)

  dens: float
    The electron density (cm^-3)

  Z: int
    The nuclear charge of the element

  Returns
  -------
  None
  """

  if settings['Ionization']=='NEI':
    z1list = range(1, Z+2)
    ionfrac = numpy.ones(len(z1list), dtype=float)

  elif settings['Ionization']=='CIE':
    z1list = range(1, Z+2)

    # time to calculate the ionization balance
    if settings['UseIonBalanceTable']:
      # read the ionization balance table
      ionfrac = atomdb.get_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, te)
    else:
      # calculate the ionization balance
      ionftmp = calc_full_ionbal(te, 1e14, Te_init=te, Zlist=[Z], extrap=True)
      ionfrac = ionftmp[Z]

  else:
    print "ERROR: settings['Ionization'] must be CIE or NEI, not %s"%(settings['Ionization'])


  abundfile = atomdb.get_filemap_file('abund',\
                                      Z,\
                                      False,\
                                      fmapfile=settings['FileMap'],\
                                      atomdbroot=os.path.expandvars('$ATOMDB'),\
                                      misc=True)

  abundances = atomdb.get_abundance(abundfile, settings['Abundances'])
  abund=abundances[Z]

  # create placeholders for all the data
  
  linelist = numpy.zeros(0, dtype=generate_datatypes('linetype'))
  contlist = {}
  pseudolist = {}
  
  # now go through each ion and assemble the data
  ebins = numpy.linspace(0.01,100,100001)
  ecent = (ebins[1:]+ebins[:-1])/2

  z1_drv_list = numpy.arange(1,Z+2, dtype=int)
  for z1_drv in range(1, Z+2):
    
    tmplinelist, tmpcontinuum, tmppseudocont = run_apec_ion(settings, te, dens, Z, z1_drv, ionfrac, abund)
    linelist = numpy.append(linelist, tmplinelist)
    contlist[z1_drv] = tmpcontinuum
    pseudolist[z1_drv] = tmppseudocont
#  out = Parallel(n_jobs=4, max_nbytes=1e2)(delayed(run_apec_ion)(settings, te, dens,Z,z1_drv,ionfrac,abund) for z1_drv in range(1,Z+2))
  #for ii,i in enumerate(out):
    #linelist = numpy.append(linelist, i[0])
    #contlist[ii+1] = i[1]
    #pseudolist[ii+1] = i[2]

  # now merge these together.
  if settings['Ionization']=='CIE':
    
    cieout = generate_cie_outputs(settings, Z, linelist, contlist, pseudolist)
    return cieout
  elif settings['Ionization']=='NEI':
    ionftmp= calc_full_ionbal(te, 1e14, Te_init=te, Zlist=[Z], extrap=True)
    ionfrac_nei = ionftmp[Z]
    neiout = generate_nei_outputs(settings, Z, linelist, contlist, pseudolist, ionfrac_nei)
    return neiout

def generate_nei_outputs(settings, Z, linelist, contlist, pseudolist, ionfrac_nei):
  """
  Convert a linelist and continuum values into a non-equilibrium AtomDB fits output
  
  Parameters
  ----------
  settings: dictionary
    The settings read from the apec.par file by parse_par_file

  Z: int
    The nuclear charge of the element
  
  linelist: numpy.array(dtype=linelisttype)
    The list of lines, separated by ion
  
  contlist: dict
    Dictionary with the different continuum contributions from each ion.
    Each is an array of ph cm^3 s^-1 bin^-1
  
  pseudolist: dict
    Dictionary with the different pseudocontinuum contributions from each ion.
    Each is an array of ph cm^3 s^-1 bin^-1

  Returns
  -------
  None
  """
  
  # ok, so we are first up going to combine the lines from the different ions
  linelist.sort(order=['element','ion','upperlev','lowerlev'])
  igood = numpy.ones(len(linelist), dtype=bool)

  linelist= linelist[igood]
  

  

  ebins =  make_vector_nbins(settings['LinearGrid'], \
                             settings['GridMinimum'], \
                             settings['GridMaximum'], \
                             settings['NumGrid'])

  pseudo = {}
  cont = {}  
  # now do some weak line filtering
  igood = numpy.ones(len(linelist), dtype=bool)
  print "initially we have %i lines for Z =%i"%(len(linelist), Z)
  for z1 in range(1, Z+2):
    ionfrac = ionfrac_nei[z1-1]
    mineps = settings['NEIMinEpsilon'][0]
    pseudotmp = numpy.zeros(len(ebins)-1, dtype=float)

    pseudotmp += pseudolist[z1]
    
    for i in range(len(settings['NEIMinFrac'])):
      if ionfrac < settings['NEIMinFrac'][i]:
        mineps = settings['NEIMinEpsilon'][i+1]
    print "z1 = %i. Ionfrac = %e. mineps = %e"%(z1,ionfrac, mineps)
    
    
    weaklines = linelist[(linelist['element']==Z) &\
                (linelist['ion_drv']==z1) &\
                (linelist['epsilon']<mineps) &\
                (linelist['lambda']>const.HC_IN_KEV_A /settings['GridMaximum']) &\
                (linelist['lambda']<const.HC_IN_KEV_A /settings['GridMinimum'])]
    print "identified %i weak lines"%(len(weaklines))
    
    for line in weaklines:
      e = const.HC_IN_KEV_A /line['lambda']
      ibin = numpy.where(ebins>e)[0][0] - 1
      pseudotmp[ibin]+=line['epsilon']

    igood[(linelist['element']==Z) &\
          (linelist['ion_drv']==z1) &\
          (linelist['epsilon']<mineps)] = False
    print "Filtered by mineps %e: from %i to %i lines"%(mineps, len(igood), sum(igood))
    conttmp = contlist[z1]['rrc']+contlist[z1]['twophot']+contlist[z1]['brems']

  
    econt, contin = compress_continuum(ebins, conttmp, const.TOLERANCE, minval=1e-38)
    
    epseudo, pseudocont = compress_continuum(ebins, pseudotmp, const.TOLERANCE, minval=1e-38)
    
    cont[z1] = {}
    cont[z1]['E_Cont'] = econt
    cont[z1]['Cont'] = contin
    cont[z1]['E_Pseudo'] = epseudo
    cont[z1]['Pseudo'] = pseudocont
  
  igood[(linelist['lambda']<const.HC_IN_KEV_A /settings['GridMaximum']) |\
        (linelist['lambda']>const.HC_IN_KEV_A /settings['GridMinimum'])] = False
         
  print "Filtering lines on wavelength: keeping %i of %i lines"%\
        (sum(igood), len(igood))
  linelist = linelist[igood]

  ret={}
  ret['lines']=linelist
  maxnpseudo = 0
  maxncont = 0
  
  for i in cont.keys():
    if len(cont[i]['E_Cont'])>maxncont:
      maxncont= len(cont[i]['E_Cont'])
    if len(cont[i]['E_Pseudo'])>maxnpseudo:
      maxnpseudo= len(cont[i]['E_Pseudo'])

  ret['cont'] = numpy.zeros(len(cont.keys()), dtype=generate_datatypes('continuum', npseudo=maxnpseudo, ncontinuum=maxncont))

  for iz1, z1 in enumerate(cont.keys()):
    ret['cont']['Z'][iz1] = Z
    ret['cont']['rmJ'][iz1] = z1
    ret['cont']['N_Cont'][iz1] = len(cont[z1]['E_Cont'])
    ret['cont']['E_Cont'][iz1][:ret['cont']['N_Cont'][iz1]] = cont[z1]['E_Cont']
    ret['cont']['Continuum'][iz1][:ret['cont']['N_Cont'][iz1]] = cont[z1]['Cont']
    ret['cont']['N_Pseudo'][iz1] = len(cont[z1]['E_Pseudo'])
    ret['cont']['E_Pseudo'][iz1][:ret['cont']['N_Pseudo'][iz1]] = cont[z1]['E_Pseudo']
    ret['cont']['Pseudo'][iz1][:ret['cont']['N_Pseudo'][iz1]] = cont[z1]['Pseudo']
    
  print "returning ret['lines'] with length %i"%(len(ret['lines']))
  return ret
  
  
  

def generate_cie_outputs(settings, Z, linelist, contlist, pseudolist):
  """
  Convert a linelist and continuum values into an equilibrium AtomDB fits output
  
  Parameters
  ----------
  settings: dictionary
    The settings read from the apec.par file by parse_par_file

  Z: int
    The nuclear charge of the element
  
  linelist: numpy.array(dtype=linelisttype)
    The list of lines, separated by ion
  
  contlist: dict
    Dictionary with the different continuum contributions from each ion.
    Each is an array of ph cm^3 s^-1 bin^-1
  
  pseudolist: dict
    Dictionary with the different pseudocontinuum contributions from each ion.
    Each is an array of ph cm^3 s^-1 bin^-1

  Returns
  -------
  None
  """
  
  # ok, so we are first up going to combine the lines from the different ions
  linelist.sort(order=['element','ion','upperlev','lowerlev'])
  igood = numpy.ones(len(linelist), dtype=bool)

#  pickle.dump(linelist, open('whole_linedebugZ%i.pkl'%(Z), 'wb'))
  
  for i in range(1,len(linelist)):
    if ((linelist['element'][i]==linelist['element'][i-1]) &\
        (linelist['ion'][i]==linelist['ion'][i-1]) &\
        (linelist['upperlev'][i]==linelist['upperlev'][i-1]) &\
        (linelist['lowerlev'][i]==linelist['lowerlev'][i-1])):
      linelist['epsilon'][i] += linelist['epsilon'][i-1]
      igood[i-1] = False
  
  linelist= linelist[igood]
  
  linelist_equ = numpy.zeros(len(linelist), dtype=generate_datatypes('linelist_cie'))

  for i in linelist_equ.dtype.names:
    linelist_equ[i] = linelist[i]


  ebins =  make_vector_nbins(settings['LinearGrid'], \
                             settings['GridMinimum'], \
                             settings['GridMaximum'], \
                             settings['NumGrid'])
  pseudo   = numpy.zeros(len(ebins)-1, dtype=float)

  # now do some weak line filtering
  MinEpsilon = settings['MinEpsilon']
  weaklines = linelist_equ[(linelist_equ['epsilon']< MinEpsilon) &\
                       (linelist_equ['lambda']>const.HC_IN_KEV_A /settings['GridMaximum']) &\
                       (linelist_equ['lambda']<const.HC_IN_KEV_A /settings['GridMinimum'])]

  for line in weaklines:
    e = const.HC_IN_KEV_A /line['lambda']
    ibin = numpy.where(ebins>e)[0][0] - 1
    pseudo[ibin]+=line['epsilon']

  linelist_equ = linelist_equ[linelist_equ['epsilon'] > MinEpsilon]


  # now do the same for the continuum:

  cont_rrc = numpy.zeros(len(ebins)-1, dtype=float)
  cont_2ph = numpy.zeros(len(ebins)-1, dtype=float)
  cont_bre = numpy.zeros(len(ebins)-1, dtype=float)
  
  for z1 in  contlist.keys():
    cont_rrc += contlist[z1]['rrc']
    cont_2ph += contlist[z1]['twophot']
    cont_bre += contlist[z1]['brems']
    
  for z1 in  pseudolist.keys():
    pseudo += pseudolist[z1]
  
  # compress things

#  ret = {}
#  ret['lines']= linelist_equ
#  ret['cont_rrc'] = cont_rrc
#  ret['cont_2ph'] = cont_2ph
#  ret['cont_bre'] = cont_bre
#  ret['pseudo'] = pseudo
#  ret['ebins'] = ebins
#  ret['Z'] = Z
  
  # create the real outputs
  
  econt, cont = compress_continuum(ebins, cont_rrc+cont_2ph+cont_bre, const.TOLERANCE, minval=1e-38)
  epseudo, pseudo = compress_continuum(ebins, pseudo, const.TOLERANCE, minval=1e-38)
  
  # make adjustments for "zero"
  


  pseudo[pseudo < 0] = 0
  igood = numpy.ones(len(pseudo), dtype=bool)
  for i in range(1, len(pseudo)-1):
    if ((pseudo[i]<=0) &\
        (pseudo[i-1]<=0) &\
        (pseudo[i+1]<=0)):
      igood[i] = False
  pseudo = pseudo[igood]
  epseudo = epseudo[igood]
        
  
  
  ret={}
  ret['lines']=linelist_equ
  maxnpseudo = len(pseudo)
  maxncont = len(cont)

  ret['cont'] = numpy.zeros(1, dtype=generate_datatypes('continuum', npseudo=maxnpseudo, ncontinuum=maxncont))

  ret['cont']['Z'][0] = Z
  ret['cont']['rmJ'][0] = 0
  ret['cont']['N_Cont'][0] = maxncont
  ret['cont']['E_Cont'][0][:maxncont] = econt
  ret['cont']['Continuum'][0][:maxncont] = cont
  ret['cont']['N_Pseudo'][0] = maxnpseudo
  ret['cont']['E_Pseudo'][0][:maxnpseudo] = epseudo
  ret['cont']['Pseudo'][0][:maxnpseudo] = pseudo
  
  return ret
  
    


def gather_rates(Z, z1, te, dens, datacache=False, settings=False,\
                 do_la=True, do_ai=True, do_ec=True, do_pc=True,\
                 do_ir=True):
  """
  fetch the rates for all the levels of Z, z1

  Parameters
  ----------
  Z: int
    The nuclear charge of the element
  z1 : int
    ion charge +1
  te : float
    temperture (Kelvin)
  dens: float
    electron density (cm^-3)
  settings : dict
    See description in atomdb.get_data
  datacache : dict
    Used for caching the data. See description in atomdb.get_data

  Returns
  -------
  up: numpy.array(float)
    Initial level of each transition
  lo: numpy.array(float)
    Final level of each transition
  rate: numpy.array(float)
    Rate for each transition (in s-1)
  """
  Te_arr, dummy = util.make_vec(te)

  lvdat = atomdb.get_data(Z, z1, 'LV', datacache=datacache, \
                            settings = settings)
  nlev = len(lvdat[1].data)
  
  diagterms = numpy.zeros(nlev)
  # get the LA data:
  if ['AAUT_REF'] in lvdat[1].data.names:
    has_sum_lv = True
  else:
    has_sum_lv = False

  if has_sum_lv:
    diagterms+= lvdat[1].data['AAUT_TOT']+lvdat[1].data['ARAD_TOT']

  laup = numpy.zeros(0, dtype=int)
  lalo = numpy.zeros(0, dtype=int)
  larate = numpy.zeros(0, dtype=float)

  if do_la:
    ladat = atomdb.get_data(Z, z1, 'LA', datacache=datacache, \
                            settings = settings)

    if ladat != False:
      laup = numpy.zeros(len(ladat[1].data), dtype=int)
      lalo = numpy.zeros(len(ladat[1].data), dtype=int)
      larate = numpy.zeros(len(ladat[1].data), dtype=float)

      laup[:] = ladat[1].data['Upper_Lev'][:] - 1
      lalo[:] = ladat[1].data['Lower_Lev'][:] - 1
      larate[:] = ladat[1].data['Einstein_A'][:]
      
      if not(has_sum_lv):
        for i in range(len(laup)):
          diagterms[laup[i]] +=larate[i]


      # create dummy results

  # get the AI data:
  aiup = numpy.zeros(0, dtype=int)
  ailo = numpy.zeros(0, dtype=int)
  airate = numpy.zeros(0, dtype=float)

  if do_ai:
    aidat = atomdb.get_data(Z, z1, 'AI', datacache=datacache, \
                            settings = settings)
    if aidat != False:
      aiup = numpy.zeros(len(aidat[1].data), dtype=int)
      ailo = numpy.zeros(len(aidat[1].data), dtype=int)
      airate = numpy.zeros(len(aidat[1].data), dtype=float)

      aiup[:] = aidat[1].data['Level_Init'][:] - 1
      ailo[:] = 0 # punt all of this into the ground state for maths reasons
      airate[:] = aidat[1].data['Auto_Rate'][:]
      if not(has_sum_lv):
        for i in range(len(aiup)):
          diagterms[aiup[i]] +=airate[i]


  # get the EC data:

  ecup = numpy.zeros(0, dtype=int)
  eclo = numpy.zeros(0, dtype=int)
  ecrate = numpy.zeros(0, dtype=float)

  if do_ec:
    ecdat = atomdb.get_data(Z, z1, 'EC', datacache=datacache, \
                            settings = settings)
    if ecdat != False:

      lvdat = atomdb.get_data(Z, z1, 'LV', datacache=datacache, \
                            settings = settings)

      ecup = numpy.zeros(len(ecdat[1].data), dtype=int)
      eclo = numpy.zeros(len(ecdat[1].data), dtype=int)
      ecrate = numpy.zeros(len(ecdat[1].data), dtype=float)

      decup = numpy.zeros(len(ecdat[1].data), dtype=int)
      declo = numpy.zeros(len(ecdat[1].data), dtype=int)
      decrate = numpy.zeros(len(ecdat[1].data), dtype=float)
      idex = 0
    # need to loop and calculate each result

#    Te_arr = numpy.array([te])
      deglarr = lvdat[1].data['LEV_DEG'][ecdat[1].data['lower_lev']-1]
      deguarr = lvdat[1].data['LEV_DEG'][ecdat[1].data['upper_lev']-1]
      deltaearr = lvdat[1].data['ENERGY'][ecdat[1].data['upper_lev']-1]-\
                  lvdat[1].data['ENERGY'][ecdat[1].data['lower_lev']-1]
  
  
      for i in range(len(ecdat[1].data)):
        # check if we need to swap things
        if deltaearr[i] >=0:
          degl = deglarr[i]
          degu = deguarr[i]
          ilo = ecdat[1].data['Lower_Lev'][i] - 1
          iup = ecdat[1].data['Upper_Lev'][i] - 1
        else:
          degl = deguarr[i]
          degu = deglarr[i]
          ilo = ecdat[1].data['Upper_Lev'][i] - 1
          iup = ecdat[1].data['Lower_Lev'][i] - 1
          deltaearr[i] *= -1
          
        
        ladat = atomdb.get_data(Z, z1, 'LA', datacache=datacache, \
                            settings = settings)

        exc,dex, tmp = atomdb.calc_maxwell_rates(ecdat[1].data['coeff_type'][i],\
                                     ecdat[1].data['min_temp'][i],\
                                     ecdat[1].data['max_temp'][i],\
                                     ecdat[1].data['temperature'][i],\
                                     ecdat[1].data['effcollstrpar'][i],\
                                     deltaearr[i]/1e3, Te_arr, Z, \
                                     degl, degu,\
                                     force_extrap=True,\
                                     levdat=lvdat,ladat=ladat, \
                                     lolev=ilo+1, uplev=iup+1)

        ecup[i] = ilo
        eclo[i] = iup
        ecrate[i] = exc*dens

        if dex > 0:
          decup[idex] = iup
          declo[idex] = ilo
          decrate[idex] = dex*dens
          idex += 1

    # now merge the de-excitation data
      decup = decup[:idex]
      declo = declo[:idex]
      decrate = decrate[:idex]

      ecup = numpy.append(ecup, decup)
      eclo = numpy.append(eclo, declo)
      ecrate = numpy.append(ecrate, decrate)

    # create dummy results
#      if not(has_sum_lv):
      for i in range(len(ecup)):
        diagterms[ecup[i]] +=ecrate[i]


  # get the PC data:
  pcup = numpy.zeros(0, dtype=int)
  pclo = numpy.zeros(0, dtype=int)
  pcrate = numpy.zeros(0, dtype=float)

  if do_pc:
    pcdat = atomdb.get_data(Z, z1, 'PC', datacache=datacache, \
                            settings = settings)
    if pcdat != False:
      pcup = numpy.zeros(len(pcdat[1].data), dtype=int)
      pclo = numpy.zeros(len(pcdat[1].data), dtype=int)
      pcrate = numpy.zeros(len(pcdat[1].data), dtype=float)

      dpcup = numpy.zeros(len(pcdat[1].data), dtype=int)
      dpclo = numpy.zeros(len(pcdat[1].data), dtype=int)
      dpcrate = numpy.zeros(len(pcdat[1].data), dtype=float)
      idex = 0
    # need to loop and calculate each result

#    Te_arr = numpy.array([te])
      deglarr = lvdat[1].data['LEV_DEG'][pcdat[1].data['lower_lev']-1]
      deguarr = lvdat[1].data['LEV_DEG'][pcdat[1].data['upper_lev']-1]
      deltaearr = lvdat[1].data['ENERGY'][pcdat[1].data['upper_lev']-1]-\
                  lvdat[1].data['ENERGY'][pcdat[1].data['lower_lev']-1]



      for i in range(len(pcdat[1].data)):


        exc,dex, tmp = atomdb.calc_maxwell_rates(pcdat[1].data['coeff_type'][i],\
                                     pcdat[1].data['min_temp'][i],\
                                     pcdat[1].data['max_temp'][i],\
                                     pcdat[1].data['temperature'][i],\
                                     pcdat[1].data['effcollstrpar'][i],\
                                     deltaearr[i]/1e3, Te_arr, Z, \
                                     deglarr[i], deguarr[i],\
                                     force_extrap=True)


        pcup[i] = pcdat[1].data['Lower_Lev'][i] - 1
        pclo[i] = pcdat[1].data['Upper_Lev'][i] - 1
        pcrate[i] = exc*dens

        if dex > 0:
          dpcup[idex] = pcdat[1].data['Upper_Lev'][i] - 1
          dpclo[idex] = pcdat[1].data['Lower_Lev'][i] - 1
          dpcrate[idex] = dex*dens
          idex += 1

    # now merge the de-excitation data
      dpcup = dpcup[:idex]
      dpclo = dpclo[:idex]
      dpcrate = dpcrate[:idex]

      pcup = numpy.append(pcup, dpcup)
      pclo = numpy.append(pclo, dpclo)
      pcrate = numpy.append(pcrate, dpcrate)

#      if not(has_sum_lv):
      for i in range(len(pcup)):
        diagterms[pcup[i]] +=pcrate[i]

  # get the IR data for colln ionization:
  irup = numpy.zeros(0, dtype=int)
  irlo = numpy.zeros(0, dtype=int)
  irrate = numpy.zeros(0, dtype=float)
  if do_ir:
    irdat = atomdb.get_data(Z, z1, 'IR', datacache=datacache, \
                            settings = settings)
    if irdat != False:
      irup = numpy.zeros(len(irdat[1].data), dtype=int)
      irlo = numpy.zeros(len(irdat[1].data), dtype=int)
      irrate = numpy.zeros(len(irdat[1].data), dtype=float)
  
      iir = 0
      # need to loop and calculate each result
  
      Te_arr = numpy.array(te)
      irtmp = irdat[1].data[(irdat[1].data['TR_TYPE']=='XI') | \
                            (irdat[1].data['TR_TYPE']=='CI')]
  
      deglarr = lvdat[1].data['LEV_DEG'][irtmp['level_init']-1]
  
      
      lvdatp1 = atomdb.get_data(Z,z1+1,'LV', settings=settings, datacache=datacache)
      
      if not (lvdatp1==False) :
  
        deglarr = lvdatp1[1].data['LEV_DEG'][irtmp['level_final']-1]
      else:
        deglarr = numpy.ones(len(irtmp['level_final']), dtype=int)
  
      ionpot = float(irdat[1].header['ionpot'])
  
      for i in range(len(irdat[1].data)):
        if not irdat[1].data['TR_TYPE'][i] in ['XI','CI']:
          continue
        else:
          
          rate=atomdb.get_maxwell_rate(Te_arr, irdat, i, lvdat, \
                                     lvdatap1=lvdatp1, ionpot=ionpot, \
                                     exconly=True)
  
  
          irup[iir] = irdat[1].data['level_init'][i] - 1
          irlo[iir] = 0
          irrate[iir] = rate*dens
          iir += 1
  
  
      # now merge the de-excitation data
      irup = irup[:iir]
      irlo = irlo[:iir]
      irrate = irrate[:iir]
#      for i in range(len(irup)):
#        print "IR %i %i = %e"%(irup[i], irlo[i], irrate[i])
#      if not(has_sum_lv):
      for i in range(len(irup)):
        diagterms[irup[i]] +=irrate[i]


  
  up_out = numpy.append(laup, numpy.append(aiup, numpy.append(ecup, numpy.append(pcup, irup))))
  lo_out = numpy.append(lalo, numpy.append(ailo, numpy.append(eclo, numpy.append(pclo, irlo))))
  rate_out = numpy.append(larate, numpy.append(airate, numpy.append(ecrate, numpy.append(pcrate, irrate))))

  up_out = numpy.append(up_out, numpy.arange(nlev, dtype=int))
  lo_out = numpy.append(lo_out, numpy.arange(nlev, dtype=int))
  rate_out = numpy.append(rate_out, diagterms*-1)
  
  
  return up_out, lo_out, rate_out

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def generate_datatypes(dtype, npseudo=0, ncontinuum=0):
  """
  returns the various data types needed by apec

  Parameters
  ----------
  dtype : string
    One of "linetype", "cielinetype", "continuum"
  npseudo : int (default=0)
    Number of pseudocontinuum points for "continuum" type
  ncontinuum : int (default=0)
    Number of continuum points for "continuum" type
    
  
  Returns
  -------
  numpy.dtype
    The data dtype in question
  """
  if dtype == 'linetype':
    ret = numpy.dtype({'names':['lambda',\
                                 'lambda_err',\
                                 'epsilon',\
                                 'epsilon_err',\
                                 'element',\
                                 'ion', \
                                 'element_drv',\
                                 'ion_drv', \
                                 'upperlev',\
                                 'lowerlev'],\
                        'formats':[numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int]})
                                   
  elif dtype =='linelist_cie':
    ret = numpy.dtype({'names':['lambda',\
                                 'lambda_err',\
                                 'epsilon',\
                                 'epsilon_err',\
                                 'element',\
                                 'ion', \
                                 'upperlev',\
                                 'lowerlev'],\
                        'formats':[numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int]})

  elif dtype == 'linetype_cap':
    ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Ion', \
                                 'Element_Drv',\
                                 'Ion_Drv', \
                                 'UpperLev',\
                                 'LowerLev'],\
                        'formats':[numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int]})

  elif dtype == 'linelist_cie_cap':
    ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Ion', \
                                 'UpperLev',\
                                 'LowerLev'],\
                        'formats':[numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.float,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int,\
                                   numpy.int]})


  elif dtype == 'continuum':
    if ncontinuum==0:
      ncontinuum+=1
    if npseudo==0:
      npseudo+=1
      
    ret = numpy.dtype({'names':['Z','rmJ','N_Cont','E_Cont','Continuum','Cont_Err','N_Pseudo','E_Pseudo','Pseudo','Pseudo_Err'],\
                       'formats':[int, int, \
                                  int, (float, ncontinuum), (float, ncontinuum),(float, ncontinuum),\
                                  int, (float, npseudo), (float, npseudo),(float, npseudo)]})
  elif dtype == 'lineparams':
    ret = numpy.dtype({'names':['kT','EDensity','Time','Nelement','Nline'],\
                       'formats':[float, float, float, int, int]})
  elif dtype == 'cocoparams':
    ret = numpy.dtype({'names':['kT','EDensity','Time','NElement','NCont', 'NPseudo'],\
                       'formats':[float, float, float, int, int, int]})
    
  else:
    print "Unknown dtype %s in generate_datatypes"%(dtype)
  return ret

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def solve_level_pop(init,final,rates,settings):
  """
  Solve the level population
  
  Parameters
  ----------
  init : array(int)
    The initial level for each transition
  final : array(int)
    The initial level for each transition
  rate : array(float)
    The rate for each transition
    
  settings: dictionary
    The settings read from the apec.par file by parse_par_file

  Returns
  -------
  array(float)
    The level population
  """

  import scipy.sparse as sparse
  from scipy.sparse.linalg import spsolve
  #stuff = {}
  #stuff['init'] = init
  #stuff['final'] = final
  #stuff['rates'] = rates
  
#  pickle.dump(stuff, open('tmpA.pkl','wb'))
#  pickle.dump(stuff, open('tmpB.pkl','wb'))
#  zzz=raw_input()
  # creat the arrays
  nlev = max([max(init), max(final)])+1

  if nlev <= const.NLEV_NOSPARSE:
    # convert to a regular solver
    matrixA = numpy.zeros([nlev,nlev], dtype=float)
    matrixB = numpy.zeros([nlev], dtype=float)
    
    for i in range(len(init)):
      matrixA[final[i], init[i]] += rates[i]
      #matrixA[init[i], init[i]] -= rates[i]

    # popn conservation
    matrixB[0] = 1.0
    matrixA[0,:] = 1.0


        
    # bug-u-fix
    for i in range(1, len(matrixB)):
      if matrixA[i,i] >= 0:
        matrixA[i,i]=-1e10
        print "ATieing level %i to ground with rate 1e10"%(i) 
    
#    a = {}
#    a['A'] = matrixA
#    a['B'] = matrixB
#    pickle.dump(a, open('dump_a.pkl','wb'))
    
    
    
    try:
      popn = numpy.linalg.solve(matrixA, matrixB)
    except numpy.linalg.linalg.LinAlgError:
      raise


  else:
    # add into sparse solver

#    tmp={}
#    tmp['init'] = init
#    tmp['final'] = final
#    tmp['rate'] = rates
#    pickle.dump(tmp, open('dump_tmp.pkl','wb'))


    matrixA={}
    matrixB = numpy.zeros(nlev, dtype=float)
#    matrixA['init'] = numpy.append(init, init)
#    matrixA['final'] = numpy.append(final, init)
#    matrixA['rate'] = numpy.append(rates, -1.0*rates)
    matrixA['init'] = init
    matrixA['final'] = final
    matrixA['rate'] = rates
    
    # filter for the ground state levels
    i = matrixA['final']>0
    matrixA['final'] = matrixA['final'][i]
    matrixA['init'] = matrixA['init'][i]
    matrixA['rate'] = matrixA['rate'][i]
    # add in continuity eqn

    matrixA['final'] = numpy.append(matrixA['final'], \
                                 numpy.zeros(nlev, dtype=int))
    matrixA['init'] = numpy.append(matrixA['init'], \
                                 numpy.arange(nlev, dtype=int))
    matrixA['rate'] = numpy.append(matrixA['rate'], \
                                 numpy.ones(nlev, dtype=float))

    
    

    A = sparse.coo_matrix((matrixA['rate'],\
                          (matrixA['final'],matrixA['init'])), \
                          shape=(nlev,nlev)).tocsr()
    
    hasdat = numpy.zeros(len(matrixB), dtype=bool)
    for i in range(1,len(hasdat)):
      if A[i,i]>=0.0:
         A[i,i] = -1e10
         print "BTieing level %i to ground with rate 1e10"%(i) 

    matrixB[0] = 1.0
    tmp={}
    #tmp['A'] = matrixA
    #tmp['B'] = matrixB
#    pickle.dump(tmp, open('dump_tmp.pkl','wb'))
    popn = spsolve(A, matrixB)
    # check the solution
#    soln = numpy.allclose(numpy.dot(A, popn), matrixB)
#    if soln==False:
#      print "ERROR Solving population matrix!"
#    for i in range(len(popn)):
      #print i, popn[i]
    #print popn
    #zzz=raw_input('mycheck1')  
  
  return popn

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def do_lines(Z, z1, lev_pop, N_e, datacache=False, settings=False, z1_drv_in=-1):
  """
  Convert level populations into line lists
  
  Parameters
  ----------
  Z: int
    The nuclear charge of the element
  z1 : int
    Ion charge +1 of ion (e.g. 6 for C VI)
  lev_pop : array(float)
    The level population for the ion. Should already have elemental abundance
    and ion fraction multiplied in.
  N_e : float
    Electron Density (cm^-3)
  datacache : dict
    Used for caching the data. See description in atomdb.get_data
  settings : dict
    See description in atomdb.get_data
  z1_drv_in : int
    the driving ion for this calculation, if not z1 (defaults to z1)
  
  Returns
  -------
  linelist: numpy.dtype(linetype)
    The list of lines and their emissivities. see generate_datatypes
  twophot: array(float)
    The two-photon continuum on the grid specified by the settings
    If settings['TwoPhoton'] is False, then returns a grid of zeros.
  """
  
  ladat = atomdb.get_data(Z,z1,'LA', datacache=datacache, settings=settings)
  lvdat = atomdb.get_data(Z,z1,'LV', datacache=datacache, settings=settings)
  ebins = make_vector_nbins(settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])
  twoph = numpy.zeros(settings['NumGrid'], dtype=float)  
  linelist = numpy.zeros(len(ladat[1].data), \
             dtype= generate_datatypes('linetype'))
  goodlines = numpy.ones(len(linelist), dtype=bool)
  
  if z1_drv_in < 0:
    z1_drv = z1
  else:
    z1_drv = z1_drv_in
    
  # now do this
  linelist['epsilon'] = ladat[1].data.field('einstein_a') * \
                        lev_pop[ladat[1].data.field('upper_lev')-1]/N_e



  for iline, line in enumerate(ladat[1].data):
    if numpy.isfinite(line['wave_obs']):
      if line['wave_obs'] > 0:
        linelist[iline]['lambda'] = line['wave_obs']
      else:
        linelist[iline]['lambda'] = line['wavelen']
      
    else:
      linelist[iline]['lambda'] = line['wavelen']
  
  linelist['lambda_err'][:] = numpy.nan
  linelist['epsilon_err'][:] = numpy.nan
  linelist['element'][:] = Z
  linelist['element_drv'][:] = Z
  linelist['ion'][:] = z1
  linelist['ion_drv'][:] = z1_drv
  linelist['upperlev'][:] = ladat[1].data['upper_lev']
  linelist['lowerlev'][:] = ladat[1].data['lower_lev']

  # I have a linelist. Yay.
  
  # now check for 2 photon transitions
  
  if (Z-z1==1):
    # He-like:
    nup = lvdat[1].data['n_quan'][linelist['upperlev']-1]
    lup = lvdat[1].data['l_quan'][linelist['upperlev']-1]
    degup = lvdat[1].data['lev_deg'][linelist['upperlev']-1]
    deggd = lvdat[1].data['lev_deg'][0]
    ila = numpy.where((nup==2) &\
                      (lup==0) &\
                      (degup==deggd) &\
                      (linelist['lowerlev']==1))[0]
    if len(ila)>0:

      # do the 2 photon fun
      twoph = numpy.zeros(len(ebins)-1, dtype=float)  
      if settings['TwoPhoton']:
#        if len(ila) > 1:
        
        for iila in ila:    
          tmp2ph={}
          tmp2ph['lambda'] = linelist['lambda'][iila]
          tmp2ph['lev_pop'] = lev_pop[ladat[1].data.field('upper_lev')[iila]-1]/N_e
          tmp2ph['einstein_a'] = ladat[1].data.field('einstein_a')[iila]
          twoph += atomdb.calc_two_phot(tmp2ph['lambda'],\
                                       tmp2ph['einstein_a'],\
                                       tmp2ph['lev_pop'], ebins)    
      goodlines[ila]=False


  elif (Z-z1==0):
    # H-like
    nup = lvdat[1].data['n_quan'][linelist['upperlev']-1]
    lup = lvdat[1].data['l_quan'][linelist['upperlev']-1]
    degup = lvdat[1].data['lev_deg'][linelist['upperlev']-1]
    deggd = lvdat[1].data['lev_deg'][0]
    ila = numpy.where((nup==2) &\
                      (lup==0) &\
                      (degup==deggd)&\
                      (linelist['lowerlev']==1))[0]
    if len(ila)>0:
      twoph = numpy.zeros(len(ebins)-1, dtype=float)  
      if settings['TwoPhoton']:
        for iila in ila:

          tmp2ph={}
          tmp2ph['lambda'] = linelist['lambda'][iila]
          tmp2ph['lev_pop'] = lev_pop[ladat[1].data.field('upper_lev')[iila]-1]/N_e
          tmp2ph['einstein_a'] = ladat[1].data.field('einstein_a')[iila]
          
          twoph += atomdb.calc_two_phot(tmp2ph['lambda'],\
                                     tmp2ph['einstein_a'],\
                                     tmp2ph['lev_pop'], ebins)
      goodlines[ila]=False
    
  elif (Z-z1==3):
    nup = lvdat[1].data['n_quan'][linelist['upperlev']-1]
    lup = lvdat[1].data['l_quan'][linelist['upperlev']-1]
    degup = lvdat[1].data['lev_deg'][linelist['upperlev']-1]
    deggd = lvdat[1].data['lev_deg'][0]
    ila = numpy.where((linelist['upperlev']==2) &\
                      (linelist['lowerlev']==1))[0]
    if len(ila)>0:
      twoph = numpy.zeros(len(ebins)-1, dtype=float)  
      if settings['TwoPhoton']:
        for iila in ila:
          tmp2ph={}
          tmp2ph['lambda'] = linelist['lambda'][iila]
          tmp2ph['lev_pop'] = lev_pop[ladat[1].data.field('upper_lev')[iila]-1]/N_e
          tmp2ph['einstein_a'] = ladat[1].data.field('einstein_a')[iila]
      
          twoph += atomdb.calc_two_phot(tmp2ph['lambda'],\
                                     tmp2ph['einstein_a'],\
                                     tmp2ph['lev_pop'], ebins)
      goodlines[ila]=False

  linelist = linelist[goodlines]
  
  return linelist, twoph
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def calc_satellite(Z, z1, T, datacache=False, settings=False):
  """
  Calcaulate DR satellite lines
  
  Parameters
  ----------
  Z: int
    The nuclear charge of the element  Z: int
  z1 : int
    Recombined Ion charge +1 of ion (e.g. 5 for C VI -> C V)
  te: float
    The electron temperature (K)
  settings: dictionary
    The settings read from the apec.par file by parse_par_file

  Returns
  -------
  array(linelist)
    List of DR lines
  """  
  
  # recombining ion charge
  z1_drv=z1+1
  
  drdat = atomdb.get_data(Z,z1_drv,'DR', settings=settings, datacache=datacache)
  lvdat = atomdb.get_data(Z,z1_drv,'LV', settings=settings, datacache=datacache)
  lvdatrec = atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)
  
  
#  pseudocont = numpy.zeros(settings['NumGrid']-1, dtype=float)

#  print "start DR"

  if drdat==False:
    linelist = numpy.zeros(0,dtype=generate_datatypes('linetype'))
    lev_rates_in = 0.0
  else:

    if not(lvdatrec):
      # no level data for recombined ion.
      lomax = max(drdat[1].data['lowerlev'])
      lev_rates_in=numpy.zeros(lomax+1, dtype=float)
    
    else:
      lev_rates_in = numpy.zeros(len(lvdatrec[1].data), dtype=float)


    linelist = numpy.zeros(len(drdat[1].data),dtype=generate_datatypes('linetype'))
    kT = const.KBOLTZ*T
    for iline in range(len(linelist)):
      epsilson = 0.0
      ll = drdat[1].data['lower_lev'][iline]
      lu = drdat[1].data['upper_lev'][iline]
      lam = drdat[1].data['wavelen'][iline]
      if (numpy.isfinite(drdat[1].data['wave_obs'][iline]) &\
          (drdat[1].data['wave_obs'][iline]>0.0)):
        lam = drdat[1].data['wave_obs'][iline]
      lamerr = drdat[1].data['wave_err'][iline]
      e_excite = drdat[1].data['E_excite'][iline]
      
      if drdat[1].data['dr_type'][iline]==const.ROMANIK:
        q_exc = drdat[1].data['SatelInt'][iline]
        epsilon = q_exc*numpy.exp(-e_excite/kT)/(T**1.5)
      elif drdat[1].data['dr_type'][iline]==const.SAFRANOVA:
        gl = lvdat[1].data['lev_deg'][0]
        q_exc = drdat[1].data['SatelInt'][iline]
        epsilon = const.SAF_COEFF* (const.RYDBERG/kT)**1.5 * \
                  (q_exc/gl) * numpy.exp(-e_excite/kT)
      else:
        print "Error in calc_satellite: unknown DR type %i"%\
               (drdat[1].data['type'][iline])
        epsilon = numpy.nan
      if ll >= len(lev_rates_in):
        print "warning: DR satellite line recombining into non existant level %i of ion Z=%i, z1=%i"%\
              (ll, Z, z1)
      else:
        lev_rates_in[ll-1] += epsilon

      linelist[iline]['lambda'] = lam
      linelist[iline]['lambda_err'] = lamerr
      linelist[iline]['epsilon'] = epsilon
      linelist[iline]['epsilon_err'] = numpy.nan
      linelist[iline]['element'] = Z
      linelist[iline]['ion'] = z1
      linelist[iline]['element_drv'] = Z
      linelist[iline]['ion_drv'] = z1_drv
      linelist[iline]['upperlev'] = lu
      linelist[iline]['lowerlev'] = ll
      
    linelist = linelist[numpy.where(linelist['epsilon'] > 0)[0]]

    
  return linelist, lev_rates_in


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calc_recomb_popn(levpop, Z, z1, z1_drv,T, dens, drlevrates, rrlevrates,\
                     settings=False, datacache=False):
  """
  Calculate the level population of a recombined ion
  
  Parameters
  ----------
  levpop: array(float)
    Level populations, already taking into account elemental abundance
    and ion fraction of z1_drv
  Z: int
  z1: int
  z1_drv: int
  te: electron temperature (K)
  dens: electron density (cm^-3)
  drlevrates: array(float)
    Rates into each level from DR calculations
  rrlevrates: array(float)
    Rates into each level from RR calculations
  
  Returns
  -------
  array(float)
    Level population
  """
  import scipy.sparse as sparse
  from scipy.sparse.linalg import spsolve
  
  # levpop at this point should alread have the corrected abundance in 
  # there
  
  lvdat = atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)
  if not lvdat:
    nlev = 1
    levpop = numpy.zeros(1, dtype=float)
    return levpop
  nlev = len(lvdat[1].data)
  Tarr, dummy = util.make_vec(T)
  
  if nlev > const.NLEV_NOSPARSE:
    
  
  # sort the levels
    aidat = atomdb.get_data(Z,z1,'AI', settings=settings, datacache=datacache)
    if aidat:
      ailev = numpy.array(util.unique(aidat[1].data['level_init']))-1
      nailev = len(ailev)
      isbound = numpy.ones(nlev, dtype=bool)
      isbound[ailev]=False
    else:
      nailev = 0
      isbound = numpy.ones(nlev, dtype=bool)

    #sortdat = sort_levels(isbound, nlev, nailev, Z, z1, settings)

    recombrate = numpy.zeros(nlev, dtype=float)
  
  # get the recomb data
  
    irdat = atomdb.get_data(Z, z1, 'IR',  settings=settings, datacache=datacache)
  
    for iir, ir in enumerate(irdat[1].data):
      # check we have the right data types
      if ir['TR_TYPE'] in ['RR','DR','XR']:
        recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)*dens
        recombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]
  
  # so recombrate has all the influxes
    matrixB = recombrate
#    print "Sum recomb rates in level>1:", sum(matrixB[1:])
    
    
    
    matrixA = {}
    matrixA['init'], matrixA['final'], matrixA['rate']=\
    gather_rates(Z, z1, T, dens, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)
#    matrixA['init'] = numpy.array(ladat[1].data['upper_lev'], dtype=int)-1
#    matrixA['final'] = numpy.array(ladat[1].data['lower_lev'], dtype=int)-1
#    matrixA['rate'] = numpy.array(ladat[1].data['einstein_a'], dtype=float)

    ### FIXME HERE

    matrixB *= -1 
    nlev = len(matrixB)
    
    # append extras onto matrixA:
    
#    matrixA['final'] = numpy.append(matrixA['final'],matrixA['init'])
#    matrixA['init'] = numpy.append(matrixA['init'],matrixA['init'])
#    matrixA['rate'] = numpy.append(matrixA['rate'],-1*matrixA['rate'])
    
    # remove all matrix A from and to level 0
    i = matrixA['final']>0
    matrixA['final'] = matrixA['final'][i]
    matrixA['init'] = matrixA['init'][i]
    matrixA['rate'] = matrixA['rate'][i]
    
    # subtract 1 from the levels
    matrixA['init']-= 1
    matrixA['final']-= 1
    
    # solve
    A = sparse.coo_matrix((matrixA['rate'],\
                           (matrixA['final'],matrixA['init'])), \
                           shape=(nlev-1,nlev-1)).tocsr()
    levpop_this = numpy.zeros(len(matrixB))
#    for i in numpy.where(matrixB != 0)[0]:
#      print i, matrixB[i]
    if sum(matrixB[1:] < 0):
      levpop_this[1:] = spsolve(A, matrixB[1:])



#    levpop_this = calc_cascade(recombrate, Z, z1, isbound, sortdat, T, settings, noauto=True)
  else:
    rrrecombrate = numpy.zeros(nlev, dtype=float)
    drrecombrate = numpy.zeros(nlev, dtype=float)
    irdat = atomdb.get_data(Z, z1, 'IR', settings=settings, datacache=datacache)

    havedrrate=False
    haverrrate=False
    for iir, ir in enumerate(irdat[1].data):
      # check we have the right data types
      if ir['TR_TYPE'] in ['RR','XR']:

        recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)
        
        rrrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*dens

        if ((ir['TR_TYPE']=='RR') & (ir['level_final']>1)):
          haverrrate=True
      if ir['TR_TYPE'] in ['DR']:
        recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)
        drrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*dens
        if ((ir['TR_TYPE']=='DR') & (ir['level_final']>1)):
          havedrrate=True

    if havedrrate:
      tmpdrlevrates=0.0
    else:
      tmpdrlevrates=drlevrates
    if haverrrate:
      tmprrlevrates=0.0
    else:
      tmprrlevrates=rrlevrates
    if isinstance(rrlevrates, numpy.ndarray):
      sumrrlevrates = sum(rrlevrates)
    else:
      sumrrlevrates = 0.0
    if isinstance(drlevrates, numpy.ndarray):
      sumdrlevrates = sum(drlevrates)
    else:
      sumdrlevrates = 0.0

    print "DR: sum from satellite lines: %e, sum from IR file: %e" %\
          (sumdrlevrates, sum(drrecombrate))
    print "RR: sum from PI xsections: %e, sum from IR file: %e" %\
          (sumrrlevrates, sum(rrrecombrate))
#    for i in range(len(rrlevrates)):
#      print "%i %e %e %e"%(i+i, rrlevrates[i], rrrecombrate[i], rrlevrates[i]/rrrecombrate[i])
#    zzz=raw_input()
    matrixB = rrrecombrate+drrecombrate+tmpdrlevrates+tmprrlevrates
    matrixA = numpy.zeros([nlev,nlev],dtype=float)
    #irdat.close()

#    print "nlev", nlev
    ladat = atomdb.get_data(Z, z1, 'LA', settings=settings, datacache=datacache)
#    print "max(la['upper_lev']): ",max(ladat[1].data['upper_lev'])

    matrixA_in = {}
    matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
    gather_rates(Z, z1, T, dens, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

    for i in range(len(matrixA_in['init'])):
      matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]

    # solve unless matrixB ==0
    if sum(matrixB[1:])>0:
      matrixB = -1*matrixB
      levpop_this = calc_cascade_population(matrixA, matrixB)
    else:
      levpop_this = numpy.zeros(nlev)
    
    

  
  return levpop_this
  
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


def calc_cascade_population(matrixA, matrixB):
  
  # replace the first line in the matrix with the population 
  # conservation equation (sum of level pops = 1)
  
  matrixB[0] = 1.0
  matrixA[0,:] = 1.0

  mb =matrixB[1:]
  ma =matrixA[1:,1:]
    
  # solve
  try:
    popn = numpy.linalg.solve(ma,mb)
  except numpy.linalg.linalg.LinAlgError:
    if ma[0,0]==0.0:
#      print 'hacking ma[0,0]=1.0'
      ma[0,0]=-1.0
      try:
        popn = numpy.linalg.solve(ma,mb)
      except numpy.linalg.linalg.LinAlgError:
        print "failed again"

        # look for levels with no way to ground. Put in a -1rate
        for i in range(ma.shape[0]):
          if ma[i,i]>= 0.0:
            ma[i,i]=-1.0

        try:
          popn = numpy.linalg.solve(ma,mb)
        except:
          print 'triple fail'
          print ma
          print mb
          raise
#    print 'mb'
#    print mb
#    tmp = open('t.pkl','w')
#    dat = {}
#    dat['mb']=mb
#    dat['ma']=ma
#    pickle.dump(dat,tmp)
#    tmp.close()
    
  #check
  soln = numpy.allclose(numpy.dot(ma, popn), mb)
  
  if soln==False:
    print "ERROR Solving population matrix!"
  popn=numpy.append(numpy.array([0.0]), popn)
  return popn
  

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
    
def calc_ioniz_popn(levpop, Z, z1, z1_drv,T, Ne, settings=False, \
                    datacache=False, do_xi=False):
  """
  Calculate the level population due to ionization into the ion
  
  Parameters
  ----------
  levpop: array(float)
    The level population of the parent ion. Should already have abundance
    and ion fraction built in.
  Z: int
  z1: int
  z1_drv: int
  T: float
  Ne: float
  settings: dict
  datacache: dict
  do_xi: bool
    Include collisional ionization
  
  Returns
  -------
  levpop_out: array(float)
    The level populations of the Z,z1 ion
  """
  # levpop at this point should alread have the corrected abundance in 
  # there
  
  # This calculates the population of levels of z1, given a previous
  # ion's (z1-1) population of levpop.
  import scipy.sparse as sparse
  from scipy.sparse.linalg import spsolve

  lvdat = atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)
  
  # if we have no lv data, ignore.
  if not util.keyword_check(lvdat):
    nlev = 1
    return numpy.array([0.0])
  nlev = len(lvdat[1].data)
  
  # get populating rate from previous ion
  aidat = atomdb.get_data(Z, z1-1, 'AI', settings=settings, datacache=datacache)
  ionizrateai=numpy.zeros(nlev, dtype=float)
  ionizrateir=numpy.zeros(nlev, dtype=float)
  if aidat:
    for ai in aidat[1].data:
      ionizrateai[ai['level_final']-1] += levpop[ai['level_init']-1]*\
                                     ai['auto_rate']
  
    #aidat.close()


  if do_xi:

    irdat = atomdb.get_data(Z, z1-1, 'IR', settings=settings, datacache=datacache)
    ionpot = float(irdat[1].header['ionpot'])
    if z1 >1:
      lvdatm1 = atomdb.get_data(Z, z1-1, 'LV', settings=settings, datacache=datacache)
  # go through each excitation, have fun
  
    for iir, ir in enumerate(irdat[1].data):
      if ir['TR_TYPE'] in ['XI']:
        Te =  numpy.array([T])
        
        ionrate=atomdb.get_maxwell_rate(Te, irdat, iir, lvdatm1, \
                                     lvdatap1=lvdat, ionpot=ionpot)
        ionizrateir[ir['level_final']-1] += levpop[ir['level_init']-1]*\
                                       ionrate
#  print 
  ionizrate=ionizrateir+ionizrateai

  matrixB = ionizrate
  
  # save some time if there is nothing to ionize.
  if sum(matrixB[1:]) ==0:
    levpop_this = numpy.zeros(len(matrixB))
    return levpop_this
 
 
  matrixA_in={}
  matrixA_in['init'], matrixA_in['final'], matrixA_in['rate'] = \
   gather_rates(Z, z1, T, Ne, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=True, do_ec=False, do_pc=False,\
                 do_ir=False)
#  tmp = {}
#  tmp['A']=matrixA_in
#  tmp['B']=matrixB
#  pickle.dump(tmp, open('tmp_%i_%i_%i.pkl'%(Z,z1_drv,z1),'wb'))

#  print "Z = %i, z1 = %i, z1_drv = %i"%(Z, z1, z1_drv)
#  print "MatrixB"
  
#  for i in range(len(matrixB)):
#    print i, matrixB[i]
    
  # fix the rates
  hasdat = numpy.zeros(len(matrixB), dtype=bool)
  for i in range(len(matrixA_in['init'])):
    if matrixA_in['init'][i]==matrixA_in['final'][i]:
      if matrixA_in['rate'][i] >=0.0:
        matrixA_in['rate'][i] -=1e10
        print "CTieing level %i to ground with rate 1e10"%(i) 
#  j = numpy.where(hasdat==False)[0]
#  if len(j) > 0:
#    tmpinit=j
#    tmpfinal=j
#    tmprate = numpy.zeros(len(j))
#    tmprate[:]=-1e-10
#    matrixA_in['init'] = numpy.append(matrixA_in['init'], tmpinit)
#    matrixA_in['final'] = numpy.append(matrixA_in['final'], tmpfinal)
#    matrixA_in['rate'] = numpy.append(matrixA_in['rate'], tmprate)
 
  
  
  if (nlev <= const.NLEV_NOSPARSE):
    # convert to a regular solver
    matrixA = numpy.zeros([nlev,nlev], dtype=float)
    
    for i in range(len(matrixA_in['init'])):
      matrixA[matrixA_in['final'][i], matrixA_in['init'][i]] += matrixA_in['rate'][i]
#      matrixA[matrixA_in['init'][i], matrixA_in['init'][i]] -= matrixA_in['rate'][i]

    # popn conservation
    matrixB[0] = 1.0
    matrixA[0,:] = 1.0

    # bug-u-fix
    for i in range(1, len(matrixB)):
      if matrixA[i,i] >= 0:
        matrixA[i,i]=-1e10
        print "FIXING matrixA[%i,%i] = -1.0"%(i,i)

    popn = numpy.zeros(nlev)
    try:
      popn[1:] = numpy.linalg.solve(matrixA[1:,1:], matrixB[1:]*-1)
    except numpy.linalg.linalg.LinAlgError:
      "EEK ERROR!"
      raise

    #if sum(matrixB[1:])<1e-40:
      #levpop_this = numpy.zeros(len(matrixB))
    #else:
      #matrixB = -1*matrixB
      #levpop_this = numpy.zeros(len(matrixB))
      
      # check for zeros
      #for iLev in range(1,len(matrixB)):
         #if not(matrixA[iLev,iLev]< 0.0):
         #  matrixA[iLev,iLev]=-1e10
      #levpop_this[1:] = numpy.linalg.solve(matrixA[1:,1:], matrixB[1:])
#      calc_cascade_population(matrixA, matrixB)

  else:
    # add into sparse solver
    matrixA={}
    matrixB *= -1 
    nlev = len(matrixB)
    
    if sum(matrixB)>=0:
      return numpy.zeros(len(matrixB))
    
    
    # add in the reverse processes
#    matrixA['init'] = numpy.append(matrixA['init'], matrixA['init'])
#    matrixA['final'] = numpy.append(matrixA['final'], matrixA['init'])
#    matrixA['rate'] = numpy.append(matrixA['rate'], -1*matrixA['rate'])
    
    # remove ground level
    i = (matrixA_in['init']>0) & (matrixA_in['final']>0)
    
    matrixA['init'] = matrixA_in['init'][i]
    matrixA['final'] = matrixA_in['final'][i]
    matrixA['rate'] = matrixA_in['rate'][i]
    
    
    # subtract 1 from the levels
    matrixA['init']-=1
    matrixA['final']-=1
    
    
         
    A = sparse.coo_matrix((matrixA['rate'],\
                           (matrixA['final'],matrixA['init'])), \
                           shape=(nlev-1,nlev-1)).tocsr()

    popn = numpy.zeros(len(matrixB))
    popn[1:] = spsolve(A, matrixB[1:])
    #print levpop_this
    #zzz=raw_input('argh')
  return popn
  
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
 

def run_apec_ion(settings, te, dens, Z, z1, ionfrac, abund):
  """
  Run the APEC code using the settings provided for an individual ion.

  Parameters
  ----------
  settings: dictionary
    The settings read from the apec.par file by parse_par_file
  te: float
    The electron temperature (K)
  dens: float
    The electron density (cm^-3)
  Z: int
    The nuclear charge of the element
  z1: int
    The ion charge +1 of the ion
  ionfrac: float
    The fractional abundance of this ion (between 0 and 1)
  abund: float
    The elemental abundance of the element (normalized to H)

  Returns
  -------
  linelist : numpy array
    List of line details and emissivities
  continuum : array
    Continuum emission in photons bin-1 s-1.
    This is a 3-item dict, with "rrc", "twophot", "brems" entries
    for each continuum source
  pseudocont : array
    Pseudo Continuum emission in photons bin-1 s-1
  """

  # get the data.
  z1_drv=z1*1
  # First, get the elemental abundance
  # make some compatilibty copies of entries in the settings file
  settings['filemap'] = settings['FileMap']
  settings['atomdbroot'] = os.path.expandvars('$ATOMDB')



  # get the output energy bins for the continuum
  ebins = make_vector_nbins(settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])

  ## FIXME CUTOFF FOR MIN IONPOP
  linelist = numpy.zeros(0, dtype=generate_datatypes('linetype'))
  pseudo = numpy.zeros(settings['NumGrid'], dtype=float)
  continuum = {}
  continuum['brems'] = numpy.zeros(settings['NumGrid'], dtype=float)
  continuum['twophot'] = numpy.zeros(settings['NumGrid'], dtype=float)
  continuum['rrc'] = numpy.zeros(settings['NumGrid'], dtype=float)


  ## FIXME CUTOFF FOR MIN IONPOP
  if ionfrac[z1_drv-1] < const.MIN_IONPOP:
    print "OMITTING Z=%i, z1=%i as ionfrac %e is below threshold of %e"%\
           (Z, z1_drv, ionfrac[z1_drv-1], const.MIN_IONPOP)
    return  linelist, continuum, pseudo
  print "NOT OMITTING Z=%i, z1=%i as ionfrac %e is above threshold of %e"%\
           (Z, z1_drv, ionfrac[z1_drv-1], const.MIN_IONPOP)

  # set up the datacache

  datacache = {}

  # Find the number of levels
  lvdat = atomdb.get_data(Z,z1,'LV', datacache=datacache, settings=settings)
  linelist_exc = numpy.zeros(0, dtype= generate_datatypes('linetype'))
  linelist_dr = numpy.zeros(0, dtype= generate_datatypes('linetype'))
  linelist_rec = numpy.zeros(0, dtype= generate_datatypes('linetype'))
  if lvdat!= False:
  
    # find the number of levels
    nlev = len(lvdat[1].data)

    # gather all the level to level rates

    up, lo, rates = gather_rates(Z, z1, te, dens, datacache=datacache, settings=settings)
    
    # solve everything
    lev_pop = solve_level_pop(up,lo,rates, settings)

    # just in case, add zeros to lengthen the lev_pop appropriately
    if len(lev_pop) < nlev:
      lev_pop = numpy.append(lev_pop, numpy.zeros(nlev-len(lev_pop), dtype=float))
  
    # fix any sub-zero level populations
    lev_pop[lev_pop<0] = 0.0  
    
  # now we have the level populations, make a line list for each ion
  # need to make a 2 photon spectrum too
    twoph = numpy.zeros(settings['NumGrid'], dtype=float)
    
  # scale lev_pop by the ion and element abundance.
    print "lev_pop Z=%i, z1=%i,z1_drv=%i, abund*ionfrac=%e, sum(pop)=%e:"%(Z,z1, z1, abund*ionfrac[z1-1], sum(lev_pop)*abund*ionfrac[z1-1])
#    for i in range(len(lev_pop)):
#      print i, lev_pop[i]
    lev_pop *= abund*ionfrac[z1-1]
    for i in range(len(lev_pop)):
      print i, lev_pop[i], lev_pop[i]/(abund*ionfrac[z1-1])
    
    linelist_exc,  continuum['twophot'] = do_lines(Z, z1, lev_pop, dens, datacache=datacache, settings=settings, z1_drv_in=z1_drv)
  
    print "Excitation Z=%i z1=%i z1_drv=%i created %i lines"%(Z, z1, z1_drv, len(linelist_exc))
  # remove this as this conversion now done to lev_pop before calling do_lines
    #linelist_exc['epsilon']*=ionfrac[z1-1]*abund
    #continuum['twophot']*=ionfrac[z1-1]*abund
  else:
    lev_pop=numpy.ones(1, dtype=float)*abund*ionfrac[z1-1]
    print "lev_pop Z=%i, z1=%i,z1_drv=%i, abund*ionfrac=%e, sum(pop)=%e: (ZEROS)"%(Z,z1, z1, abund*ionfrac[z1-1], sum(lev_pop)*abund*ionfrac[z1-1])
    for i in range(len(lev_pop)):
       print i, lev_pop[i]
    
  # calculate some continuum processes: brems
  
  if settings['Bremsstrahlung'] ==True:
    brems = do_brems(Z, z1, te, 1.0, settings['BremsType'], ebins)
    # scale for  ion and element abundance.
    continuum['brems']=brems*abund*ionfrac[z1-1]
  else:
    continuum['brems']=numpy.zeros(len(ebins)-1, dtype=float)

  # now look at the neighbouring ions
  datacache={}
#  z1_drv=z1*1
  if z1_drv>1:
    z1=z1_drv-1
    # do the DR satellite lines
    if settings['DRSatellite']:
      linelist_dr, drlevrates = calc_satellite(Z, z1, te, datacache=datacache, settings=settings)
      linelist_dr['epsilon']*=ionfrac[z1_drv-1]*abund
      drlevrates *=ionfrac[z1_drv-1]*abund
    else:
      linelist_dr = numpy.zeros(0, dtype= generate_datatypes(linetype))
      drlevrates = 0.0

    # Radiative Recombination
    if settings['RRC']:
      rrc, rrlevrates = atomdb.calc_rad_rec_cont(Z, z1, z1_drv, te, ebins, settings=settings, datacache=datacache)
      continuum['rrc'] = rrc*ionfrac[z1_drv-1]*abund    
      rrlevrates*=ionfrac[z1_drv-1]*abund
      
#      print "rrlevrates from Z=%i, z1_drv=%i to z1=%i"%(Z, z1_drv, z1)
#      for ilev in range(len(rrlevrates)):
#        print ilev, rrlevrates[ilev]
    else:
      continuum['rrc'] = numpy.zeros(len(ebins)-1, dtype=float)
      rrlevrates=0.0
#    for i in range(len(rrlevrates)):
#      print i+1, rrlevrates[i]
#    zzz=raw_input()
    # get the level specific recombination rates
#    print len(drlevrates)
    
    # if there is recombination to process:
    tmpdrlevrates,xxx = util.make_vec(drlevrates)
    tmprrlevrates,xxx = util.make_vec(rrlevrates)
    
    
    if sum(tmpdrlevrates) + sum(tmprrlevrates)>0:
      print "initial_lev_pop", lev_pop
      levpop_recomb=calc_recomb_popn(lev_pop, Z, z1,\
                                      z1_drv, te, dens, drlevrates,\
                                      rrlevrates,\
                                      datacache=datacache, settings=settings)
      print "final_lev_pop", lev_pop
                                      
      #zzz=raw_input()
      linelist_rec, tmptwophot = \
               do_lines(Z, z1, levpop_recomb , dens, datacache=datacache, settings=settings, z1_drv_in=z1_drv)
      continuum['twophot']+= tmptwophot
      
      print "linelist_rec Z=%i, z1=%i,z1_drv=%i, nlines=%i:"%(Z,z1, z1_drv, len(linelist_rec))
      

  # now do the ionizing cases
  linelist_ion = numpy.zeros(0,dtype= generate_datatypes('linetype'))
  if z1_drv < Z:
    datacache={}
    z1=z1_drv+1
    lev_pop_parent = lev_pop*1.0
    print "Sum lev_pop_parent[1:] = %e"%(sum(lev_pop_parent[1:]))
    while (sum(lev_pop_parent[1:]) > 1e-40) &\
          (z1 <= Z):
      
      # do we need to include collisional ionzation here?
      # only in first ion
      if z1== z1_drv+1:
        do_xi = True
      else:
        do_xi = False
      lev_pop=calc_ioniz_popn(lev_pop_parent, Z, z1, z1_drv, te, dens, \
                             settings=settings, datacache=datacache, \
                             do_xi=do_xi)
      
      lev_pop[lev_pop<const.MIN_LEVPOP] = 0.0
      if sum(lev_pop[1:]) > 0:
        linelist_ion_tmp, tmptwophot = \
               do_lines(Z, z1, lev_pop, dens, datacache=datacache, settings=settings,   z1_drv_in=z1_drv)
        print "linelist_ion Z=%i, z1=%i,z1_drv=%i, nlines=%i:"%(Z,z1, z1_drv, len(linelist_ion_tmp))

        linelist_ion = numpy.append(linelist_ion, linelist_ion_tmp)
        continuum['twophot']+=tmptwophot

      lev_pop_parent = lev_pop
      z1+=1

  # generate return data
  linelist = numpy.append(linelist_exc, numpy.append(linelist_dr, numpy.append(linelist_ion, linelist_rec)))

  # filter line list
  MinEpsilon = settings['MinEpsilon']
  if settings['Ionization']=='CIE':
    MinEpsilon*=0.001
#  istrong = linelist>=MinEpsilon
  pseudocont = numpy.zeros(len(ebins)-1, dtype=float)
  if len(linelist) > 0:
    weaklines = linelist[(linelist['epsilon']< MinEpsilon) &\
                         (linelist['lambda']>const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                         (linelist['lambda']<const.HC_IN_KEV_A /settings['GridMinimum'])]
#  print "filtered out % i weak lines"%(len(weaklines))
  
    for line in weaklines:
      e = const.HC_IN_KEV_A /line['lambda']
      ibin = numpy.where(ebins>e)[0][0] - 1
      pseudocont[ibin]+=line['epsilon']

    linelist = linelist[linelist['epsilon'] > MinEpsilon]
#  print "kept  % i strong lines"%(len(weaklines))

#  ret = {}
#  ret['lines'] = linelist
#  ret['Z'] = Z
#  ret['z1'] = z1
#  ret['pseudocont'] = pseudocont
#  ret['continuum'] = continuum
#  ret['ionfrac'] = ionfrac
#  fname = 'dump_%i_%i.pkl'%(Z,z1_drv)
#  pickle.dump(ret, open(fname, 'wb'))
#  print "wrote %s"%(fname)
  if settings['WriteIon']==True:
    ret = {}
    ret['lines'] = linelist
    ret['continuum'] = continuum
    ret['pseudocont'] = pseudocont
    ret['ionfrac'] = ionfrac
    ret['te'] = te
    ret['dens'] = dens
    ret['settings'] = settings
    ret['abund'] = abund
    fname = settings['WriteIonFname']+'.pkl'
    pickle.dump(ret, open(fname, 'wb'))
  return linelist, continuum, pseudocont

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def compress_continuum(xin, yin, tolerance, minval = 0.0):
  """
  Compress the continuum into linear interpolatable grids
  
  Parameters
  ----------
  xin : array(float)
    The bin edges (keV)
  yin : array(float)
    The continuum in photons (or ergs) cm^3 s^-1 bin^-1. 
    Should be 1 element shorter then xin
  tolerance : float
    The tolerance of the final result (if 0.01, the result will always
    be within 1% of the original value)
  
  Returns
  -------
  xout : array (float)
    The energy points of the compressed energy grid (keV)
  yout : array (float)
    The continuum, in photons(or ergs) cm^3 s^-1 keV^-1
  
  """
  
  from pyatomdb import liblinapprox
  npts = len(yin)

  xin_tmp = ctypes.c_double * npts
  yin_tmp = ctypes.c_double * npts
  ein_tmp = ctypes.c_double * npts
  xout_tmp = ctypes.c_double * npts
  xout_tmp = ctypes.c_double * npts
  yout_tmp = ctypes.c_double * npts
  wout_tmp = ctypes.c_double * npts

  x = xin_tmp()
  y = yin_tmp()
  e = ein_tmp()
  xout= xout_tmp()
  yout= yout_tmp()
  wout= wout_tmp()

  m = ctypes.c_int(npts)
  
  for i in range(npts):

    x[i] = (xin[i]+xin[i+1])/2.0
    y[i] = yin[i]/(xin[i+1]-xin[i]) # convert to per keV
    e[i] = tolerance*y[i]
  k = ctypes.c_int(0)
  ip = ctypes.c_int(-6)
  
  
  liblinapprox.linear_approx(ctypes.byref(x),\
                            ctypes.byref(y),\
                            ctypes.byref(e),\
                            ctypes.byref(m),\
                            ctypes.byref(xout),\
                            ctypes.byref(yout),\
                            ctypes.byref(wout),\
                            ctypes.byref(k),ctypes.byref(ip))
  
  nout = k.value+1
  xout = numpy.array(xout[:nout], dtype=float)
  yout = numpy.array(yout[:nout], dtype=float)
  
  # now trim if required
  i = numpy.where(yout < minval)[0]
  if len(i) > 0:
    # we have something to do
    yout[i] = 0.0
    isgood = numpy.ones(len(yout), dtype=bool)
    for i in range(1, len(yout)-1):
      if ( (yout[i]<=minval) &\
           (yout[i-1]<=minval) &\
           (yout[i+1]<=minval)):
        isgood[i] = False
    xout = xout[isgood]
    yout = yout[isgood]
           
    
    
  
  return xout, yout


def wrap_ion_directly(fname, ind, Z, z1):
  import time
  # read in the file
  settings = parse_par_file(fname)



  te = make_vector(settings['LinearTemp'], \
                   settings['TempStart'], \
                   settings['TempStep'], \
                   settings['NumTemp'])

  if settings['TempUnits']=='keV':
    te /= const.KBOLTZ


  dens = make_vector(settings['LinearDens'], \
                     settings['DensStart'], \
                     settings['DensStep'], \
                     settings['NumDens'])

  ite = ind /len(dens)
  idens = ind%len(dens)
  print ite, idens
  Te = te[ite]
  Dens = dens[idens]

  if settings['Ionization']=='NEI':
    z1list = range(1, Z+2)
    ionfrac = numpy.ones(len(z1list), dtype=float)

  elif settings['Ionization']=='CIE':
    z1list = range(1, Z+2)

    # time to calculate the ionization balance
    if settings['UseIonBalanceTable']:
      # read the ionization balance table
      ionfrac = atomdb.get_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, te)
    else:
      # calculate the ionization balance
      ionftmp = calc_full_ionbal(Te, 1e14, Te_init=Te, Zlist=[Z], extrap=True)
      ionfrac = ionftmp[Z]

  else:
    print "ERROR: settings['Ionization'] must be CIE or NEI, not %s"%(settings['Ionization'])


  abundfile = atomdb.get_filemap_file('abund',\
                                      Z,\
                                      False,\
                                      fmapfile=settings['FileMap'],\
                                      atomdbroot=os.path.expandvars('$ATOMDB'),\
                                      misc=True)

  #abundances = atomdb.get_abundance(abundfile, settings['Abundances'])
  abundances = atomdb.get_abundance()
  abund=abundances[Z]

  # update the output filename
  settings['WriteIonFname'] ="Z_%i_z1_%i_iT_%iiN_%i.pkl"%(Z,z1,ite,idens)
  settings['HDUIndex'] = ind
  settings['Te_used'] = Te
  settings['Ne_used'] = Dens
  
  x=run_apec_ion(settings, Te, Dens, Z, z1, ionfrac, abund)
  
  ret = {}
  ret['settings'] = settings
  ret['data'] = x
  fname = settings['OutputFileStem']+'_'+settings['WriteIonFname']
  pickle.dump(ret, open(fname,'wb'))
  print "wrote file %s"%(fname)
  print "Finished cleanly at %s"%(time.asctime())



def wrap_run_apec(fname, readpickle=False):
  """
  After running the APEC code ion by ion, use this to combine into 
  FITS files.

  Parameters
  ----------
  fname : string
    file name
  
  readpickle : bool
    Load apec results by element from pickle files, instead of regenerating

  Returns
  -------
  None

  """

  # get the input settings
  settings = parse_par_file(fname)

  # need to parse the IncAtoms parameter - this defines which atoms
  # we will include
  #
  # we will transfer this to a "Zlist" parameter, in the form of
  # nuclear charges

  if len(Zlist)==0:
    Zlist = settings['Zlist']

  print "I will be running Z=", Zlist
  # run for each element, temperature, density
  
  lhdulist = []
  chdulist = []
  
  seclHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=generate_datatypes('lineparams'))
  seccHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=generate_datatypes('cocoparams'))
  
  
  for iTe in range(settings['NumTemp']):

    te = make_vector(settings['LinearTemp'], \
                     settings['TempStart'], \
                     settings['TempStep'], \
                     settings['NumTemp'])[iTe]

    if settings['TempUnits']=='keV':
      te /= const.KBOLTZ


    for iDens in range(settings['NumDens']):
      dens = make_vector(settings['LinearDens'], \
                         settings['DensStart'], \
                         settings['DensStep'], \
                         settings['NumDens'])[iDens]

      # AT THIS POINT, GENERATE SHELL WHICH WILL GO IN THE HDU OF CHOICE
      
      if settings['Ionization']=='CIE':
        linedata = numpy.zeros(0,dtype=generate_datatypes('linelist_cie'))
        cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))
      if settings['Ionization']=='NEI':
        linedata = numpy.zeros(0,dtype=generate_datatypes('linetype'))
        cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))

      for Z in Zlist:
        print "Calling run_apec_element for Z=%i Te=%e dens=%e at %s"%(Z, te, dens, time.asctime())
        dat = wrap_run_apec_element(settings, te, dens, Z,iTe,iDens, readpickle=readpickle)
        # append this data to the output
        #pickle.dump(dat, open('dump_%i.pkl'%(Z),'wb'))
        linedata = numpy.append(linedata, dat['lines'])
        cocodata = continuum_append(cocodata, dat['cont'])
        print "Z=%i, nlines=%i"%(Z, len(dat['lines']))

      
      # now make an HDU for all of this
      if settings['Ionization']=='CIE':
        LHDUdat = create_lhdu_cie(linedata)
      elif settings['Ionization']=='NEI':
        pickle.dump(linedata,open('argh.pkl','wb'))
        LHDUdat = create_lhdu_nei(linedata)
        
      # now update the headers
      iseclHDUdat=iDens+iTe*settings['NumDens']
      LHDUdat.header['EXTNAME']=("EMISSIVITY","name of this binary table extension")
      LHDUdat.header['EXTVER']=(iseclHDUdat+1,"Index for this EMISSIVITY extension")

      LHDUdat.header['HDUNAME'] = ("L%.2f_%.2f"%(numpy.log10(te), numpy.log10(dens)),\
                             'Spectral emission data')
      LHDUdat.header['HDUCLASS'] = ("Proposed OGIP",\
                             'Proposed OGIP standard')
      LHDUdat.header['HDUCLAS1']=("LINE MODEL",\
                             'Line emission spectral model')
      LHDUdat.header['HDUCLAS2']=("LINE",\
                             'Emission line data')
      if settings['Ionization']=='CIE':
        LHDUdat.header['HDUVERS1']=("1.0.0",\
                               'version of format')
      elif settings['Ionization']=='NEI':
        LHDUdat.header['HDUVERS1']=("2.0.0",\
                               'version of format')

      LHDUdat.header['TEMPERATURE']=(te,\
                             'Electron temperature')
      LHDUdat.header['DENSITY']=(dens,\
                             'Electron density')
      LHDUdat.header['TIME']=(0,\
                             'IN EQUILIBRIUM')
      if settings['Ionization']=='CIE':
        tot_emiss = sum(linedata['epsilon']*const.HC_IN_ERG_A/linedata['lambda'])
      else:
        tot_emiss=0.0
      LHDUdat.header['TOT_LINE']=(tot_emiss,\
                             'Total Line Emission (erg cm^3 s^-1)')
      LHDUdat.header['N_LINES']=(len(linedata),\
                             'Number of emission lines')

      seclHDUdat['kT'][iseclHDUdat]=te*const.KBOLTZ
      seclHDUdat['EDensity'][iseclHDUdat]= dens
      seclHDUdat['Time'][iseclHDUdat]= 0.0
      seclHDUdat['Nelement'][iseclHDUdat]= len(Zlist)
      seclHDUdat['Nline'][iseclHDUdat]= len(linedata)
      
      lhdulist.append(LHDUdat)


# continuum data
      # now make an HDU for all of this
      CHDUdat = create_chdu_cie(cocodata)
      
      # now update the headers
      iseccHDUdat=iDens+iTe*settings['NumDens']

      CHDUdat.header['EXTNAME']=("EMISSIVITY","name of this binary table extension")
      CHDUdat.header['EXTVER']=(iseccHDUdat+1,"Index for this EMISSIVITY extension")

      CHDUdat.header['HDUNAME'] = ("C%.2f_%.2f"%(numpy.log10(te), numpy.log10(dens)),\
                             'Spectral emission data')
      CHDUdat.header['HDUCLASS'] = ("Proposed OGIP",\
                             'Proposed OGIP standard')
      CHDUdat.header['HDUCLAS1']=("COMP CONT MODEL",\
                             'Compressed continua spectra')
      CHDUdat.header['HDUCLAS2']=("COCO",\
                             'Compressed continuum data')
      CHDUdat.header['HDUVERS1']=("1.0.0",\
                             'version of format')
      CHDUdat.header['TEMPERATURE']=(te,\
                             'Electron temperature')
      CHDUdat.header['DENSITY']=(dens,\
                             'Electron density')
      CHDUdat.header['TIME']=("%.2e"%(0),\
                             'IN EQUILIBRIUM')
      tot_emiss = calc_total_coco(cocodata, settings)
      
      if  settings['Ionization']=='CIE':
        CHDUdat.header['TOT_COCO']=(tot_emiss,\
                               'Total Emission (erg cm^3 s^-1)')
      else:
        CHDUdat.header['TOT_COCO']=(0.0,\
                               'Total Emission (erg cm^3 s^-1)')

      seccHDUdat['kT'][iseccHDUdat]=te*const.KBOLTZ
      seccHDUdat['EDensity'][iseccHDUdat]= dens
      seccHDUdat['Time'][iseccHDUdat]= 0.0
      seccHDUdat['NElement'][iseccHDUdat]= len(cocodata)
      seccHDUdat['NCont'][iseccHDUdat]= max(cocodata['N_Cont'])
      seccHDUdat['NPseudo'][iseccHDUdat]= max(cocodata['N_Pseudo'])
      
      chdulist.append(CHDUdat)

    # make secHDUdat into a fits table
  seclHDU = create_lparamhdu_cie(seclHDUdat)
  seclHDU.header['EXTNAME']=('PARAMETERS','name of this binary table extension')
  seclHDU.header['HDUCLASS']=('Proposed OGIP','Proposed OGIP standard')
  seclHDU.header['HDUCLAS1']=('LINE MODEL','line emission spectra model')
  seclHDU.header['HDUCLAS2']=('Parameters','extension containing parameter info')
  seclHDU.header['HDUVERS1']=('1.0.0','version of format')


  seccHDU = create_cparamhdu_cie(seccHDUdat)
  seccHDU.header['EXTNAME']=('PARAMETERS','name of this binary table extension')
  seccHDU.header['HDUCLASS']=('Proposed OGIP','Proposed OGIP standard')
  seccHDU.header['HDUCLAS1']=('COMP CONT MODEL','compressed continua spectra')
  seccHDU.header['HDUCLAS2']=('Parameters','extension containing parameter info')
  seccHDU.header['HDUVERS1']=('1.0.0','version of format')


  fileroot = settings['OutputFileStem']
  if settings['Ionization']=='NEI':
    fileroot+='_nei'
  # create the Primary HDU
  PrilHDU = pyfits.PrimaryHDU()
  lhdulist.insert(0,PrilHDU)
  lhdulist.insert(1,seclHDU)
  tmplhdulist = pyfits.HDUList(lhdulist)

  PricHDU = pyfits.PrimaryHDU()
  chdulist.insert(0,PricHDU)
  chdulist.insert(1,seccHDU)
  tmpchdulist = pyfits.HDUList(chdulist)

  # now add header blurb
  generate_apec_headerblurb(settings, tmplhdulist, tmpchdulist)

  # write out results
  tmplhdulist.writeto('%s_line.fits'%(fileroot), clobber=True, checksum=True)

  if settings['Ionization']=='CIE':
    tmpchdulist.writeto('%s_coco.fits'%(fileroot), clobber=True, checksum=True)
  elif settings['Ionization']=='NEI':
    tmpchdulist.writeto('%s_comp.fits'%(fileroot), clobber=True, checksum=True)



def wrap_run_apec_element(settings, te, dens, Z, ite, idens, writepickle=False, readpickle=False):
  """
  Combine wrap_run_apec_ion results for an element

  Parameters
  ----------
  settings: dictionary
    The settings read from the apec.par file by parse_par_file

  te: float
    The electron temperature (K)

  dens: float
    The electron density (cm^-3)

  Z: int
    The nuclear charge of the element

  ite: int
    The temperature index
  
  idens: int
    The density index

  writepickle: bool
    Dump data into a pickle file. Useful for rapidly combining data
    after runs.

  readpickle: bool
    Read data from a pickle file. Useful for rapidly combining data
    after runs. Usually the result of a previous call using
    writepickle=True

  Returns
  -------
  None
  """

  if readpickle:
      setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)
      ret = pickle.load(open(setpicklefname,'rb'))
      print "read %s"%(setpicklefname)
      return ret

  linelist = numpy.zeros(0, dtype=generate_datatypes('linetype'))
  contlist = {}
  pseudolist = {}



  
  # now go through each ion and assemble the data
  ebins = numpy.linspace(0.01,100,100001)
  ecent = (ebins[1:]+ebins[:-1])/2

  for z1_drv in range(1,Z+2):
    setpicklefname = "%s_Z_%i_z1_%i_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,z1_drv,ite,idens)
    print "loading %s"%(setpicklefname)
    if not os.path.exists(setpicklefname):
      print "Warning: no such file: %s"%(setpicklefname)
      contlist[z1_drv]={}
      contlist[z1_drv]['rrc']=numpy.zeros(settings['NumGrid'], dtype=float)
      contlist[z1_drv]['twophot']=numpy.zeros(settings['NumGrid'], dtype=float)
      contlist[z1_drv]['brems']=numpy.zeros(settings['NumGrid'], dtype=float)
    
      pseudolist[z1_drv] = numpy.zeros(settings['NumGrid'], dtype=float)
    else:
      dat= pickle.load(open(setpicklefname,'rb'))    
      tmplinelist, tmpcontinuum, tmppseudocont = dat['data']
      # check for NAN
      nlines = len(tmplinelist)
      print nlines
      
      ngoodlines = sum(numpy.isfinite(tmplinelist['epsilon']))
      if nlines != ngoodlines:
        print "Bad lines found in %s"%(setpicklefname)
      linelist = numpy.append(linelist, tmplinelist)
      for key in tmpcontinuum.keys():
        tmpncont = len(tmpcontinuum[key])
        if tmpncont != sum(numpy.isfinite(tmpcontinuum[key])):
          print "Bad continuum found in %s %s"%(key, setpicklefname),
#          if key=='rrc':
#            tmpcontinuum['rrc'][numpy.isnan(tmpcontinuum['rrc'])]=0.0
          print ""
            
      contlist[z1_drv] = tmpcontinuum

      tmpncont = len(tmppseudocont)
      if tmpncont != sum(numpy.isfinite(tmppseudocont)):
        print "Bad pseudocont found in %s"%( setpicklefname)
      pseudolist[z1_drv] = tmppseudocont

  # now merge these together.
  if settings['Ionization']=='CIE':
    
    cieout = generate_cie_outputs(settings, Z, linelist, contlist, pseudolist)
    if writepickle:
      setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)
      pickle.dump(cieout, open(setpicklefname,'wb'))
      print "wrote %s"%(setpicklefname)
    return cieout
  elif settings['Ionization']=='NEI':
    ionftmp= calc_full_ionbal(te, 1e14, Te_init=te, Zlist=[Z], extrap=True)
    ionfrac_nei = ionftmp[Z]
    neiout = generate_nei_outputs(settings, Z, linelist, contlist, pseudolist, ionfrac_nei)
    if writepickle:
      setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)
      pickle.dump(neiout, open(setpicklefname,'wb'))
      print "wrote %s"%(setpicklefname)
    return neiout



#-------------------------------------------------------------------------------



def run_wrap_run_apec(fname, Z, iTe, iDens):
  """
  After running the APEC code ion by ion, use this to combine into 
  FITS files.

  Parameters
  ----------
  fname : string
    file name of par file
  Z: int
    The atomic numbers
  iTe: int
    The temperature index
  iDens: int
    The density index
    
  Returns
  -------
  None

  """

  # get the input settings
  settings = parse_par_file(fname)

  # need to parse the IncAtoms parameter - this defines which atoms
  # we will include
  #
  # we will transfer this to a "Zlist" parameter, in the form of
  # nuclear charges


  print "I will be running Z=", Z
  # run for each element, temperature, density
  
  lhdulist = []
  chdulist = []
  
  seclHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=generate_datatypes('lineparams'))
  seccHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=generate_datatypes('cocoparams'))
  
  
  te = make_vector(settings['LinearTemp'], \
                   settings['TempStart'], \
                   settings['TempStep'], \
                   settings['NumTemp'])[iTe]

  if settings['TempUnits']=='keV':
    te /= const.KBOLTZ


  dens = make_vector(settings['LinearDens'], \
                     settings['DensStart'], \
                     settings['DensStep'], \
                     settings['NumDens'])[iDens]

  # AT THIS POINT, GENERATE SHELL WHICH WILL GO IN THE HDU OF CHOICE
      
  if settings['Ionization']=='CIE':
    linedata = numpy.zeros(0,dtype=generate_datatypes('linelist_cie'))
    cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))
  if settings['Ionization']=='NEI':
    linedata = numpy.zeros(0,dtype=generate_datatypes('linetype'))
    cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))

  print "Calling run_apec_element for Z=%i Te=%e dens=%e at %s"%(Z, te, dens, time.asctime())
  dat = wrap_run_apec_element(settings, te, dens, Z,iTe,iDens, writepickle=True)
  print "Done safely"

#-------------------------------------------------------------------------------

def generate_apec_headerblurb(settings, linehdulist, cocohdulist):
  """
  Generate all the headers for an apec run, and apply them to the HDUlist.
  
  Parameters
  ----------
  
  settings: dict
    The output of read_apec_parfile
  hdulist : list or array of fits HDUs
    The hdus to have headings added.
  
  Returns
  -------
  None
  """
  
#  import subprocess
#  label = subprocess.check_output(["git","describe"])
  hl = linehdulist[0].header
  hc = cocohdulist[0].header
  
  hl.add_comment("FITS (Flexible Image Transport System) format is defined in 'Astronomy")
  hl.add_comment("and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")
  hl.append(("CONTENT",'SPECTRUM','Emission Line Project output'),end=True)
  hl.append(("FILENAME",'original','Parent File'),end=True)
  hl.append(("ORIGIN",'APEC','origin of FITS file'),end=True)
  hl.append(("HDUVERS1",'1.0.0','version of format'),end=True)


  hc.add_comment("FITS (Flexible Image Transport System) format is defined in 'Astronomy")
  hc.add_comment("and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")
  hc.append(("CONTENT",'SPECTRUM','Emission Line Project output'),end=True)
  hc.append(("FILENAME",'original','Parent File'),end=True)
  hc.append(("ORIGIN",'APEC','origin of FITS file'),end=True)
  hc.append(("HDUVERS1",'1.0.0','version of format'),end=True)

  
  fmap = atomdb.read_filemap(settings['FileMap'])
  ftype_comments={}
  ftype_comments['ir']='ionization/recombination'
  ftype_comments['lv']='energy levels'
  ftype_comments['la']='emission lines'
  ftype_comments['ec']='electron collisions'
  ftype_comments['pc']='proton collisions'
  ftype_comments['dr']='dielectronic satellites'
  ftype_comments['pi']='photoionization data'
  ftype_comments['ai']='autoionization data'
  
  
  for Z in settings['Zlist']:
    for z1 in range(1, Z+1):
      i = numpy.where((fmap['Z']==Z) & (fmap['z1']==z1))[0]
      if len(i)==0: continue
      i = i[0]
      for ftype in ['ir','lv','la','ec','pc','dr','ci','ai','em','pi']:
        if fmap[ftype][i]=='': continue
        if not os.path.exists(fmap[ftype][i]):
          a=atomdb.get_data(Z,z1,ftype)
        datsum=pyfits.getval(fmap[ftype][i], 'DATASUM', 1)
        s='CS%s%02i%02i'%(ftype.upper(), Z,z1)
        c='Checksum for %s: %s'%(atomic.spectroscopic_name(Z,z1),\
                                 ftype_comments[ftype])
        hl.append((s,datsum,c),end=True)
        hc.append((s,datsum,c),end=True)

  hl.add_comment("*** BEGIN INITIALIZATION VALUES ***", after=len(hl)-1)
  hl.add_comment("Input Files",after=len(hl)-1)
  hl.append(('SIONBAL',settings['IonBalanceTable'],'Ionization bal'),end=True)
  hl.append(('SFILEMAP',settings['FileMap'], 'Filemap file name'), end=True)
  hl.append(('SABUND_SOURCE',settings['Abundances'], 'Abundance source'), end=True)
  hl.add_comment('Output Files',after=len(hl)-1)
  hl.append(('OUTPUT_LINE','','Output line list file'), end=True)
  hl.append(('SMODELFILE',settings['OutputFileStem'],'Output file head'), end=True)
  hl.append(('OUTPUT_COCO','','Output compressed continuum file'), end=True)
  hl.append(('SMODELFILE',settings['OutputFileStem'],'Output file head'), end=True)
  hl.append(('SMODELNAME',settings['OutputFileStem'],'XSPEC Model Name'), end=True)

  hc.add_comment("*** BEGIN INITIALIZATION VALUES ***",after=len(hl)-1)
  hc.add_comment("Input Files",after=len(hl)-1)
  hc.append(('SIONBAL',settings['IonBalanceTable'],'Ionization bal'),end=True)
  hc.append(('SFILEMAP',settings['FileMap'], 'Filemap file name'), end=True)
  hc.append(('SABUND_SOURCE',settings['Abundances'], 'Abundance source'), end=True)
  hc.add_comment('Output Files',after=len(hl)-1)
  hc.append(('OUTPUT_LINE','','Output line list file'), end=True)
  hc.append(('SMODELFILE',settings['OutputFileStem'],'Output file head'), end=True)
  hc.append(('OUTPUT_COCO','','Output compressed continuum file'), end=True)
  hc.append(('SMODELFILE',settings['OutputFileStem'],'Output file head'), end=True)
  hc.append(('SMODELNAME',settings['OutputFileStem'],'XSPEC Model Name'), end=True)

  hl.add_comment("Physical Processes",after=len(hl)-1)
  hl.append(('DO_RRC',settings['RRC'],'Include Radiative Recombination Continuum'), end=True)
  hl.append(('DO_SATEL',settings['DRSatellite'],'Include Dielectronic Satellite Lines'), end=True)
  hl.append(('DO_LINES',settings['EmissionLines'],'Include Emission Lines'), end=True)
  hl.append(('DO_TWOPH',settings['TwoPhoton'],'Two Photon Continuum'), end=True)
  hl.append(('DO_BREMS',settings['Bremsstrahlung'],'Include Bremsstrahlung Continuum'), end=True)
  hl.append(('SBREMS_TYPE', settings['BremsType'], 'Bremsstrahlung type'), end=True)

  hc.add_comment("Physical Processes",after=len(hl)-1)
  hc.append(('DO_RRC',settings['RRC'],'Include Radiative Recombination Continuum'), end=True)
  hc.append(('DO_SATEL',settings['DRSatellite'],'Include Dielectronic Satellite Lines'), end=True)
  hc.append(('DO_LINES',settings['EmissionLines'],'Include Emission Lines'), end=True)
  hc.append(('DO_TWOPH',settings['TwoPhoton'],'Two Photon Continuum'), end=True)
  hc.append(('DO_BREMS',settings['Bremsstrahlung'],'Include Bremsstrahlung Continuum'), end=True)
  hc.append(('SBREMS_TYPE', settings['BremsType'], 'Bremsstrahlung type'), end=True)


  hl.add_comment('Spectral Binning',after=len(hl)-1)
  hl.append(('BIN_LIN_ENERGY', settings['LinearGrid'],'Linear energy spacing'), end=True)
  hl.append(('INUM_E',settings['NumGrid'],'Number of energy bins'), end=True)
  hl.append(('DE_START', settings['GridMinimum'],'Starting energy bin'),end=True)
  hl.append(('DE_END', settings['GridMaximum'],'Final energy bin'),end=True)
  hl.append(('DMIN_EPSILON', settings['MinEpsilon'],'[ph cm^3/s] Minimum output line emissivity'),end=True)


  hc.add_comment('Spectral Binning',after=len(hl)-1)
  hc.append(('BIN_LIN_ENERGY', settings['LinearGrid'],'Linear energy spacing'), end=True)
  hc.append(('INUM_E',settings['NumGrid'],'Number of energy bins'), end=True)
  hc.append(('DE_START', settings['GridMinimum'],'Starting energy bin'),end=True)
  hc.append(('DE_END', settings['GridMaximum'],'Final energy bin'),end=True)
  hc.append(('DMIN_EPSILON', settings['MinEpsilon'],'[ph cm^3/s] Minimum output line emissivity'),end=True)

  logdens =True
  if settings['LinearDens']:
    logdens=False

  logtemp =True
  if settings['LinearTemp']:
    logtemp=False

  hl.add_comment('Physical Parameters',after=len(hl)-1)
  hl.append(('LOG_DENS', logdens, 'Logarithmic density spacing'), end=True)
  hl.append(('INUM_DENSITIES', settings['NumDens'], 'Number of densities'), end=True)
  hl.append(('DDENSITY_START',settings['DensStart'],'Starting density'),end=True)
  hl.append(('DDENSITY_STEP', settings['DensStep'],'Density step size'), end=True)
  hl.append(('LOG_TEMP', logtemp, 'Logarithmic temperature spacing'), end=True)
  hl.append(('INUM_TEMP',settings['NumTemp'],'Number of temperatures'), end=True)
  hl.append(('DTEMP_START',settings['TempStart'],'Starting temperature'), end=True)
  hl.append(('DTEMP_STEP', settings['TempStep'], 'Temperature step size'), end=True)
  hl.add_comment('Atoms Included',after=len(hl)-1)
  for Z in settings['Zlist']:
    hl.append(('SATOM',atomic.Ztoelsymb(Z), 'Element name'), end=True)
  hl.add_comment('*** END INITIALIZATION VALUES ***',after=len(hl)-1)
  
  hc.add_comment('Physical Parameters',after=len(hl)-1)
  hc.append(('LOG_DENS', logdens, 'Logarithmic density spacing'), end=True)
  hc.append(('INUM_DENSITIES', settings['NumDens'], 'Number of densities'), end=True)
  hc.append(('DDENSITY_START',settings['DensStart'],'Starting density'),end=True)
  hc.append(('DDENSITY_STEP', settings['DensStep'],'Density step size'), end=True)
  hc.append(('LOG_TEMP', logtemp, 'Logarithmic temperature spacing'), end=True)
  hc.append(('INUM_TEMP',settings['NumTemp'],'Number of temperatures'), end=True)
  hc.append(('DTEMP_START',settings['TempStart'],'Starting temperature'), end=True)
  hc.append(('DTEMP_STEP', settings['TempStep'], 'Temperature step size'), end=True)
  hc.add_comment('Atoms Included',after=len(hl)-1)
  for Z in settings['Zlist']:
    hc.append(('SATOM',atomic.Ztoelsymb(Z), 'Element name'), end=True)
  hc.add_comment('*** END INITIALIZATION VALUES ***',after=len(hl)-1)
  
    


def solve_ionbal_eigen(Z, Te, init_pop=False, tau=False, Te_init=False, \
                       teunit='K', \
                       filename=False):
  """
  Solve the ionization balance for a range of ions using the eigenvector
  approach and files as distributed in XSPEC.

  Parameters
  ----------
  Z : int
    atomic number of element
  Te : float
    electron temperature, default in K
  init_pop : float array
    initial population of ions for non-equlibrium calculations. Will be renormalised to 1.
  tau : float
    N_e * t for the non-equilibrium ioniziation
  Te_init : float
    initial ionization balance temperature, same units as Te
  teunit : {'K' , 'keV'}
    units of temperatures (default K)
  filename : string
    Can optionally point directly to the file in question, i.e. to look at older data
    look at $HEADAS/../spectral/modelData/eigenELSYMB_v3.0.fits.
    If not set, download from AtomDB FTP site.

  Returns
  -------
  final_pop : float array
    final populations.

  """
#
#  Version 0.1 Initial Release
#  Adam Foster 16th September 2015
#
  # first, calculate the equilibrium solution

#  if (init_pop==False) & (tau==False): do_equilib=True

  init_pop_set = util.keyword_check(init_pop)
  tau_set = util.keyword_check(tau)
  Te_init_set = util.keyword_check(Te_init)
  
  if (not tau_set):
    # we need to do equilbirum
    do_equilib = True
    Te_equilib = Te
  else:
    if (not init_pop_set):
    # we need to do equilbirum
      do_equilib = True
      Te_equilib = Te_init
      if not Te_init_set:
        print "ERROR: need to specift Te_init or init_pop for NEI calculation"
        return
    else:
      do_equilib=False
  if (not tau_set):
    do_nei=False
  else:
      do_nei=True
  
  if util.keyword_check(filename):
    # we have a filename specified!
    fname = os.path.expandvars(fname)
    if not os.path.isfile(fname):
      print "Specified file %s does not exist. Exiting"%(fname)
      return
    d = pyfits.open(fname)
  else:
    d = atomdb.get_data(Z, False, 'eigen')

  telist = numpy.logspace(4,9,1251)
  
  if do_equilib:
    itelist = numpy.argsort((telist-Te_equilib)**2)
    ite = [min(itelist[:2]), max(itelist[:2])]
    Tdiff = telist[ite[1]] - telist[ite[0]]
    if Tdiff > 0.0:
      factorlow = (telist[ite[1]]-Te_equilib)/Tdiff
      factorhigh = (Te_equilib-telist[ite[0]])/Tdiff
      equilib = factorlow * d['EIGEN'].data['FEQB'][ite[0]]+\
                factorhigh * d['EIGEN'].data['FEQB'][ite[1]]
    else:
      equilib = d['EIGEN'].data['FEQB'][ite[0]]
      
  
  if do_nei:
    if (not init_pop_set):
      init_pop = equilib
    
    Tindex = numpy.argmin((telist-Te)**2)
    
    lefteigenvec = numpy.zeros([Z,Z], dtype=float)
    righteigenvec = numpy.zeros([Z,Z], dtype=float)
    for i in range(Z):
      for j in range(Z):
        lefteigenvec[i,j] = d['EIGEN'].data['VL'][Tindex][i*Z+j]
        righteigenvec[i,j] = d['EIGEN'].data['VR'][Tindex][i*Z+j]
        
    delt = 1.0/(len(telist)-1.0)
    work = numpy.zeros(Z, dtype=float)
    work = init_pop[1:] - d['EIGEN'].data['FEQB'][Tindex][1:]
    fspec = numpy.zeros(Z)
    for i in range(Z):
      for j in range(Z):
        fspec[i] +=lefteigenvec[i][j] * work[j]
    
    fspectmp = numpy.matrix(lefteigenvec) * numpy.matrix(work).transpose()
    worktmp = fspectmp *numpy.exp(d['EIGEN'].data['EIG'][Tindex] * delt * tau)

    frac = numpy.zeros(Z+1)
    for i in range(Z):
      for j in range(Z):
        frac[i+1] += worktmp[j,0]*righteigenvec[j][i]
      frac[i+1] += d['EIGEN'].data['FEQB'][Tindex][i+1]
    
    frac[frac<0.0] = 0.0
      
    if sum(frac)< 1.0:
      frac[0] = 1.0-sum(frac)
    print "frac:",frac
    
  return
  
