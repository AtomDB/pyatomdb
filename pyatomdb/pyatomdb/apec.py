"""
The apec module contains routines crucial for the APEC code. This also
includes some interfaces to external C libraries (or will, eventually).

Version 0.1 - initial release
Adam Foster September 16th 2015

"""

import numpy, copy
import util, atomdb, const
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
