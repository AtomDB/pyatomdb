#import matplotlib.pyplot as plt
#import pylab
import pyatomdb
import sys
import math
#from pyatomdb import *
import numpy, copy, pickle
import hashlib
import os, time, sys
import util, const, atomic, atomdb
import scipy, ctypes
import astropy.io.fits as pyfits
#from __main__ import *




#Type=str(input("ion/element/allelement?"))











class ion:


    """
    Defining the class ion:

    Parameters
    ----------
     Z: int
     The nuclear charge of the element
     z1 : int
     ion charge 
     Te : float
     temperture (Kelvin)
     N_e: float
     electron density (cm^-3)
     upper_level:
     Electron is making transition from this level
     lower_level:
     Electron is making transition to this level
     factor: float
     A multiplicative factor that can change the atomic data stored in $ATOMDB/APED. 
     This determines which atomic data (e.g. collision strength) is adjusted by what percentage. 
     factor = 1.1 means 10% increase in the atomic     data 
     errortype:str
     Type of the error, e.g. 'LA', 'EC' etc. 
     errortype ='LA' means Einstein A will be multiplied by the factor.  
     settings : dict
     See description in atomdb.get_data

    """

    #Z=int(sys.argv[1])
    #z1=int(sys.argv[2])
    #Te=float(sys.argv[3])
    #N_e=float(sys.argv[4])
    #upper_level = int(sys.argv[5])
    #lower_level = int(sys.argv[6])
    #factor = float(sys.argv[5])
    #errortype= str(sys.argv[6])


    
    




    
    
    
    
    


    def __init__(self, Z, z1, Te=None, N_e=None, ind=None, settings_fname='$ATOMDB/apec.par'):
        
        '''
        Z=self.Z
        z1=self.z1
        Te=self.Te
        N_e=self.N_e
        factor=self.factor
        errortype=self.errortype
        ind=self.ind
        '''
        #ionfrac= element.calc_elem_ionbal_delta(self, Z, Te, factor, errortype)
        #ionfrac= element.calc_elem_ionbal(self, Z, Te)

        #ionfrac=self.ionfrac
        self.Z=Z
        self.z1=z1
        #self.Te=Te
        #self.N_e=N_e
        Abund= atomdb.get_abundance(abundfile=False, \
                abundset='AG89', element=[-1], datacache=False, settings=False, show=False)
        self.Abund=Abund
        #fname = os.path.expandvars('$ATOMDB/apec.par')
        #settings = parse_par_file(fname)
        settings = parse_par_file(os.path.expandvars(settings_fname))
        self.settings = settings
    
        
        


        




    
    
    


    """
    Disable printing.
    """
    def blockPrint():
        sys.stdout = open(os.devnull, 'w')

    """
    Restore printing
    """
    def enablePrint():
        sys.stdout = sys.__stdout__

    
    blockPrint()
   


    def datacache(self, Z, z1, settings=False):


        """
        Creating an empty dictionary datacache
    
        Reading atomic data for:
        ‘IR’ - ionization and recombination
        ‘LV’ - energy levels
        ‘LA’ - radiative transition data (lambda and A-values)
        ‘EC’ - electron collision data
        ‘PC’ - proton collision data
        ‘DR’ - dielectronic recombination satellite line data
        ‘PI’ - XSTAR photoionization data
        ‘AI’ - autoionization data

        Parameters
        ----------
        Z: int
          The nuclear charge of the element
        z1 : int
          ion charge
        datacache : dict
          Used for caching the data. See description in atomdb.get_data
        settings : dict
          See description in atomdb.get_data

        Returns: datacache
        Storing in datacache for later use
        """
       
        datacache={}
        #Z=self.Z
        #z1=self.z1
        
        
        #for z1 in range(1, Z+1):
         # Alldat = atomdb.get_data(Z, z1, 'ALL', datacache=datacache, \
                            #settings = settings)
        Alldat = atomdb.get_data(Z, z1,  'ALL', datacache=datacache, \
                            settings = settings)
        
        return datacache




    
    






    def gather_rates(self, Te, N_e, datacache=True, settings=False,\
                 do_la=True, do_ai=True, do_ec=True, do_pc=True,\
                 do_ir=True):
      """
      fetch the rates for all the levels of Z, z1.  'IR', 'LV', 'LA', 'EC', 'PC', 'DR', 'PI', 'AI' rates are being read from datacache.
      Parameters
      ----------
      Z: int
        The nuclear charge of the element
      z1 : int
        ion charge +1
      Te : float
        temperture (Kelvin)
      N_e: float
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
      Z=self.Z
      z1=self.z1
      #datacache=self.datacache
      Te_arr, dummy = util.make_vec(Te)
      

      Atomic_data =  ion.datacache(self, Z, z1)
      

      lvdat =  Atomic_data['data'][Z][z1]['LV']

      
      
      nlev = len(lvdat[1].data)

      diagterms = numpy.zeros(nlev)
  # get the LA data:
      if 'AAUT_TOT' in lvdat[1].data.names:
        has_sum_lv = True
      else:
        has_sum_lv = False
      

      laup = numpy.zeros(0, dtype=int)
      lalo = numpy.zeros(0, dtype=int)
      larate = numpy.zeros(0, dtype=float)

      if do_la:
        t1=time.time()
        if has_sum_lv:
          diagterms+=lvdat[1].data['ARAD_TOT']

        

        ladat =  Atomic_data['data'][Z][z1]['LA']


        

        if ladat != False:

          laup = ladat[1].data['Upper_Lev'] - 1
          lalo = ladat[1].data['Lower_Lev'] - 1
          larate = ladat[1].data['Einstein_A']

          if not(has_sum_lv):
#        diagterms[laup] +=larate
            for i in range(len(laup)):
              diagterms[laup[i]] +=larate[i]
        t2 = time.time()
        


      # create dummy results

  # get the AI data:
      aiup = numpy.zeros(0, dtype=int)
      ailo = numpy.zeros(0, dtype=int)
      airate = numpy.zeros(0, dtype=float)

      if do_ai:
        if has_sum_lv:
          diagterms+= lvdat[1].data['AAUT_TOT']


        t1=time.time()
        aidat =  Atomic_data['data'][Z][z1]['AI']

        
        if aidat != False:
          aiup = aidat[1].data['Level_Init'][:] - 1
          ailo = numpy.zeros(len(aidat[1].data), dtype=int)
          airate = aidat[1].data['Auto_Rate']

          if not(has_sum_lv):
            for i in range(len(aiup)):
              diagterms[aiup[i]] +=airate[i]
        t2=time.time()
        

  # get the EC data:

      ecup = numpy.zeros(0, dtype=int)
      eclo = numpy.zeros(0, dtype=int)
      ecrate = numpy.zeros(0, dtype=float)

      if do_ec:
        t1 = time.time()
        ecdat =  Atomic_data['data'][Z][z1]['EC']
        
        if ecdat != False:

          lvdat =  Atomic_data['data'][Z][z1]['LV']
          

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

          ladat =  Atomic_data['data'][Z][z1]['LA']

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




            exc,dex, tmp = atomdb._calc_maxwell_rates(ecdat[1].data['coeff_type'][i],\
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
            ecrate[i] = exc*N_e

            if dex > 0:
              decup[idex] = iup
              declo[idex] = ilo
              decrate[idex] = dex*N_e
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
        t2 = time.time()
        

  # get the PC data:
      pcup = numpy.zeros(0, dtype=int)
      pclo = numpy.zeros(0, dtype=int)
      pcrate = numpy.zeros(0, dtype=float)

      if do_pc:
        t1 = time.time()
        pcdat =  Atomic_data['data'][Z][z1]['PC']
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


            exc,dex, tmp = atomdb._calc_maxwell_rates(pcdat[1].data['coeff_type'][i],\
                                     pcdat[1].data['min_temp'][i],\
                                     pcdat[1].data['max_temp'][i],\
                                     pcdat[1].data['temperature'][i],\
                                     pcdat[1].data['effcollstrpar'][i],\
                                     deltaearr[i]/1e3, Te_arr, Z, \
                                     deglarr[i], deguarr[i],\
                                     force_extrap=True)


            pcup[i] = pcdat[1].data['Lower_Lev'][i] - 1
            pclo[i] = pcdat[1].data['Upper_Lev'][i] - 1
            pcrate[i] = exc*N_e

            if dex > 0:
              dpcup[idex] = pcdat[1].data['Upper_Lev'][i] - 1
              dpclo[idex] = pcdat[1].data['Lower_Lev'][i] - 1
              dpcrate[idex] = dex*N_e
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
        t2=time.time()
        

  # get the IR data for colln ionization:
      irup = numpy.zeros(0, dtype=int)
      irlo = numpy.zeros(0, dtype=int)
      irrate = numpy.zeros(0, dtype=float)
      if do_ir:
        t1 = time.time()
        irdat =  Atomic_data['data'][Z][z1]['IR']

        
        if irdat != False:
          irup = numpy.zeros(len(irdat[1].data), dtype=int)
          irlo = numpy.zeros(len(irdat[1].data), dtype=int)
          irrate = numpy.zeros(len(irdat[1].data), dtype=float)

          iir = 0
      # need to loop and calculate each result

          Te_arr = numpy.array(Te)
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
              irrate[iir] = rate*N_e
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
        t2=time.time()
        

      tmp={}
      tmp['up']={}
      tmp['lo']={}
      tmp['rate']={}
      tmp['up']['ec'] = ecup
      tmp['lo']['ec'] = eclo
      tmp['rate']['ec'] = ecrate

      tmp['up']['pc'] = pcup
      tmp['lo']['pc'] = pclo
      tmp['rate']['pc'] = pcrate

      tmp['up']['ir'] = irup
      tmp['lo']['ir'] = irlo
      tmp['rate']['ir'] = irrate

      tmp['up']['ai'] = aiup
      tmp['lo']['ai'] = ailo
      tmp['rate']['ai'] = airate

      tmp['up']['la'] = laup
      tmp['lo']['la'] = lalo
      tmp['rate']['la'] = larate

      tmp['up']['diag'] = numpy.arange(nlev, dtype=int)
      tmp['lo']['diag'] = numpy.arange(nlev, dtype=int)
      tmp['rate']['diag'] = diagterms*-1

#  pickle.dump(tmp, open('tmp_gather_%i_%i.pkl'%(Z,z1),'wb'))

      t1=time.time()
      up_out = numpy.append(laup, numpy.append(aiup, numpy.append(ecup, numpy.append(pcup, irup))))
      lo_out = numpy.append(lalo, numpy.append(ailo, numpy.append(eclo, numpy.append(pclo, irlo))))
      rate_out = numpy.append(larate, numpy.append(airate, numpy.append(ecrate, numpy.append(pcrate, irrate))))
      

      up_out = numpy.append(up_out, numpy.arange(nlev, dtype=int))
      lo_out = numpy.append(lo_out, numpy.arange(nlev, dtype=int))
      rate_out = numpy.append(rate_out, diagterms*-1)
      t2=time.time()
      #self.up_out=up_out
      #self.lo_out=lo_out
      #self.rate_out=rate_out
      
      return up_out, lo_out, rate_out
      






    


    




    def ion_fraction(self, Z, z1, Te):

        """
        Parameters
        ----------
        Z: int
        The nuclear charge of the element
        z1 : int
        ion charge

        Returns: ionfraction for ion charge z1
        """

        
        
        for z1 in range(1,Z+1):
            ionfrac = element.calc_elem_ionbal(self,Z,Te)
 
        return ionfrac[z1-2]

    
    
  


    blockPrint()

    def solve_level_pop(self, init,final,rates,settings):
      """
      Solve the level population

      Parameters
      ----------
      init : array(int)
        The initial level for each transition
      final : array(int)
        The initial level for each transition
      rates : array(float)
        The rate for each transition

      settings: dictionary
        The settings read from the apec.par file by parse_par_file

      Returns
      -------
      array(float)
        The level population
      """

      settings=self.settings
      import scipy.sparse as sparse
      from scipy.sparse.linalg import spsolve

  # creat the arrays
      nlev = max([max(init), max(final)])+1
#--  print "Starting Solve Level Pop at %s"%(time.asctime())
      tstart = time.time()
      k=numpy.where((final==0) & (init==0))[0]
      if len(k)==1:
        rates[k[0]] = 0.0
      k=numpy.isfinite(rates)
      rates=rates[k]
      init=init[k]
      final=final[k]



      nlev_old = max([max(init), max(final)])+1
      irate = numpy.where(rates>1e-40)[0]



      nlev = max([max(init[irate]), max(final[irate])])


      irate = numpy.where((init<=nlev) & (final<=nlev))[0]

      init=init[irate]
      final=final[irate]
      rates=rates[irate]
      nlev = max([max(init), max(final)])+1
      print("nlev = %i, nlev_old =%i"%(nlev, nlev_old))


      if nlev <= const.NLEV_NOSPARSE:
        #print("Using regular solver.")
        #print("Starting generation of matrixA at %s"%(time.asctime()))
        t1=time.time()
    # convert to a regular solver
        matrixA = numpy.zeros([nlev,nlev], dtype=float)
        matrixB = numpy.zeros([nlev], dtype=float)

        t0 = time.time()
        t1=time.time()

        for i in range(len(init)):
          matrixA[final[i], init[i]] += rates[i]
        t2 = time.time()
        #print("time differences: %f vs %f seconds"%(t1-t0, t2-t1))


    # popn conservation
        matrixB[0] = 1.0
        matrixA[0,:] = 1.0

        #print("Starting check of diagonal terms at %s"%(time.asctime()))

    # bug-u-fix
        for i in range(1, len(matrixB)):
          if matrixA[i,i] >= 0:
            matrixA[i,i]=-1e10
            print("ATieing level %i to ground with rate 1e10"%(i))
        #print("Finished check of diagonal terms at %s"%(time.asctime()))


        t2=time.time()
        #print("Finished generation of matrixA at %s: took %f seconds"%(time.asctime(), t2-t1))

        #print("Starting calling solver at %s"%(time.asctime()))
        try:
          popn = numpy.linalg.solve(matrixA, matrixB)
        except numpy.linalg.linalg.LinAlgError:
          raise
        t3 = time.time()
        #print("Finished calling solver at %s: took %f seconds"%(time.asctime(), t3-t2))

      else:

        #print("using sparse solver")
        matrixA={}
        matrixB = numpy.zeros(nlev, dtype=float)

        matrixA['init'] = init
        matrixA['final'] = final
        matrixA['rate'] = rates


        maxlev = max(matrixA['final'][matrixA['rate']>1e-40])
        k = numpy.where(matrixA['rate']<1e-40)[0]

    # filter for the ground state levels

        i = matrixA['final']>0
        matrixA['final'] = matrixA['final'][i]
        matrixA['init'] = matrixA['init'][i]
        matrixA['rate'] = matrixA['rate'][i]
    # add in continuity eqn

        i = (matrixA['final']<=maxlev) & (matrixA['init']<=maxlev)
        matrixA['final'] = matrixA['final'][i]
        matrixA['init'] = matrixA['init'][i]
        matrixA['rate'] = matrixA['rate'][i]


        matrixA['final'] = numpy.append(matrixA['final'], \
                                 numpy.zeros(maxlev, dtype=int))
        matrixA['init'] = numpy.append(matrixA['init'], \
                                 numpy.arange(maxlev, dtype=int))
        matrixA['rate'] = numpy.append(matrixA['rate'], \
                                 numpy.ones(maxlev, dtype=float))

        A = numpy.zeros([maxlev+1,maxlev+1])
        for i in range(len(matrixA['final'])):
          A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]


        matrixB[0] = 1.0
        A[0,:] = 1.0

        hasdat = numpy.zeros(len(matrixB), dtype=bool)
        for i in range(1,maxlev+1):
          if A[i,i]>=0.0:
            A[i,i] = -1e10
            #print("BTieing level %i to ground with rate 1e10"%(i))

        matrixB[0] = 1.0

        popn = numpy.zeros(nlev_old)
        matrixA=0
        popn[:maxlev+1] = numpy.linalg.solve(A, matrixB[:maxlev+1])


    # sparse solver sometimes returns small negative pop. set to 0.

        popn[popn<0] = 0.0
    
      tfinish=time.time()
      #print("Finished Solve Level Pop at %s: took %i seconds"%(time.asctime(), tfinish-tstart))
      return popn


    enablePrint()







    def level_population(self, Te, N_e):


        """
        Calculate level populations for all the levels of Z, z1

        Parameters
        ----------
        Z: int
          The nuclear charge of the element
        z1 : int
          ion charge
        Te : float
          Temperature in K
        N_e : float
          Electron Density (cm^-3)
        settings : dict
        """
        datacache={}
        settings=self.settings
        up, lo, rates = ion.gather_rates(self, Te, N_e)
        lev_pop = ion.solve_level_pop (self, up, lo, rates, settings)
        #lev_pop *= ion.ion_fraction(self, Z, z1)
        return lev_pop


    




        




    
    


    def calc_recomb_popn(self, levpop, z1_drv, Te, N_e, drlevrates, rrlevrates,\
                     settings=False, datacache=False, dronly=False,\
                     rronly=False):
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
      Te: electron temperature (K)
      N_e: electron density (cm^-3)
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
      Z=self.Z
      z1=self.z1
      

      # levpop at this point should alread have the corrected abundance in
      # there

      Atomic_data =  ion.datacache(self, Z, z1)
      lvdat =  Atomic_data['data'][Z][z1]['LV']
      
      if not lvdat:
        nlev = 1
        levpop = numpy.zeros(1, dtype=float)
        return levpop
      nlev = len(lvdat[1].data)
      Tarr, dummy = util.make_vec(Te)

      if nlev > const.NLEV_NOSPARSE:
        print("using sparse solver for recomb")

    
        aidat =  Atomic_data['data'][Z][z1]['AI']
        
        if aidat:
          ailev = numpy.array(util.unique(aidat[1].data['level_init']))-1
          nailev = len(ailev)
          isbound = numpy.ones(nlev, dtype=bool)
          isbound[ailev]=False
        else:
          nailev = 0
          isbound = numpy.ones(nlev, dtype=bool)

    

        recombrate = numpy.zeros(nlev, dtype=float)

        #get the recomb data
        irdat =  Atomic_data['data'][Z][z1]['IR']
        

        for iir, ir in enumerate(irdat[1].data):
        #check we have the right data types
          if ir['TR_TYPE'] in ['RR','DR','XR']:
            recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)*N_e
            if not (numpy.isfinite(recrate)):
              print("iir=%i, recrate is not finite!"%(iir))
            else:
              recombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]

        maxlev = numpy.where(recombrate > 0)[0]
        if len(maxlev) == 0:  # so recombrate has all the influxes
          return numpy.zeros(nlev)
        maxlev=maxlev[-1]
        matrixB = recombrate
        #print "Sum recomb rates in level>1:", sum(matrixB[1:])



        matrixA = {}
        matrixA['init'], matrixA['final'], matrixA['rate']=\
          ion.gather_rates(self, Z, z1, Te, N_e, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)


        matrixB *= -1
        nlev = len(matrixB)



        i = (matrixA['final']>0) & (matrixA['init']>0)
        matrixA['final'] = matrixA['final'][i]
        matrixA['init'] = matrixA['init'][i]
        matrixA['rate'] = matrixA['rate'][i]

    #remove all matrix A from and to high levels
        i = (matrixA['final']<maxlev+1) & (matrixA['init']<maxlev+1)
        matrixA['final'] = matrixA['final'][i]
        matrixA['init'] = matrixA['init'][i]
        matrixA['rate'] = matrixA['rate'][i]

    #subtract 1 from the levels
        matrixA['init']-= 1
        matrixA['final']-= 1

        A  = numpy.zeros([maxlev,maxlev])
        for i in range(len(matrixA['final'])):
          A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]


        levpop_this = numpy.zeros(len(matrixB))

        if sum(matrixB[1:] < 0):
          levpop_this[1:maxlev+1] = numpy.linalg.solve(A, matrixB[1:maxlev+1])



      else:

        #print("using regular solver for recomb")

        rrrecombrate = numpy.zeros(nlev, dtype=float)
        drrecombrate = numpy.zeros(nlev, dtype=float)
        irdat =  Atomic_data['data'][Z][z1]['IR']
        

        havedrrate=False
        haverrrate=False
        for iir, ir in enumerate(irdat[1].data):
      # check we have the right data types
          if ir['TR_TYPE'] in ['RR','XR']:

            recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)

            rrrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*N_e

            if ((ir['TR_TYPE'] in ['RR','XR']) & (ir['level_final']>1)):
              haverrrate=True
          if ir['TR_TYPE'] in ['DR','XD']:
            recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)
            drrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*N_e
            if ((ir['TR_TYPE'] in ['DR','XD']) & (ir['level_final']>1)):
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

        #print("DR: sum from satellite lines: %e, sum from IR file: %e" %\
              #(sumdrlevrates, sum(drrecombrate)))
        #print("RR: sum from PI xsections: %e, sum from IR file: %e" %\
              #(sumrrlevrates, sum(rrrecombrate)))

        matrixB = rrrecombrate+drrecombrate+tmpdrlevrates+tmprrlevrates
        if dronly:
          matrixB = drrecombrate+tmpdrlevrates
        if rronly:
          matrixB = rrrecombrate+tmprrlevrates


        matrixA = numpy.zeros([nlev,nlev],dtype=float)

        ladat =  Atomic_data['data'][Z][z1]['LA']
        

        matrixA_in = {}
        matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
          ion.gather_rates(self, Te, N_e, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

        #datacache={}


        for i in range(len(matrixA_in['init'])):
          matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]

    # solve unless matrixB ==0
        if sum(matrixB[1:])>0:
          matrixB = -1*matrixB
          levpop_this = setup.calc_cascade_population(self, matrixA, matrixB)
        else:
          levpop_this = numpy.zeros(nlev)


      #print("level population for recombination into Z=%i, z1=%i, z1_drv=%i, levpop=%i"%\
            #(Z, z1, z1_drv,levpop))
      #for i in range(len(levpop_this)):
        #print(i, levpop_this[i])
      return levpop_this









    def calc_ioniz_popn(self, levpop, Z, z1, z1_drv, Te, N_e, settings=False, \
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
        factor: float  
        errortype:str
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

        datacache = {}
        #print("Starting calc_ioniz_popn at %s"%(time.asctime()))
        Atomic_data =  ion.datacache(self, Z, z1)
        lvdat =  Atomic_data['data'][Z][z1]['LV']
        
        #lvdat = atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)

  # if we have no lv data, ignore.
        if not util.keyword_check(lvdat):
            nlev = 1
            return numpy.array([0.0])
        nlev = len(lvdat[1].data)

  # get populating rate from previous ion
        Atomic_data_prev =  ion.datacache(self, Z, z1-1)
        aidat = Atomic_data_prev['data'][Z][z1-1]['AI']
        ionizrateai=numpy.zeros(nlev, dtype=float)
        ionizrateir=numpy.zeros(nlev, dtype=float)

        #print("Starting calc_ioniz_popn aidat loop at %s"%(time.asctime()))
        if aidat:
            tmp_pop = levpop[aidat[1].data['level_init']-1]
            for iai in range(len(aidat[1].data)-1):
              ionizrateai[aidat[1].data['level_final'][iai]-1] += \
                  tmp_pop[iai]*aidat[1].data['auto_rate'][iai]

    #aidat.close()
        #print("Finished calc_ioniz_popn aidat loop at %s"%(time.asctime()))


        #print("Starting calc_ioniz_popn xidat loop at %s"%(time.asctime()))
        if do_xi:

            irdat = Atomic_data_prev['data'][Z][z1-1]['IR']
            #irdat = atomdb.get_data(Z, z1-1, 'IR', settings=settings, datacache=datacache)
            ionpot = float(irdat[1].header['ionpot'])
            if z1 >1:
                lvdatm1 = Atomic_data_prev['data'][Z][z1-1]['LV']
  # go through each excitation, have fun

            for iir, ir in enumerate(irdat[1].data):
                if ir['TR_TYPE'] in ['XI']:
                    Tarr =  numpy.array([Te])
                    ionrate=atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdatm1, \
                                     lvdatap1=lvdat, ionpot=ionpot)
                    
                    ionizrateir[ir['level_final']-1] += levpop[ir['level_init']-1]*\
                                       ionrate




        ionizrate=ionizrateir+ionizrateai
        matrixB = ionizrate

  # save some time if there is nothing to ionize.

        if sum(matrixB[1:]) ==0:
            levpop_this = numpy.zeros(len(matrixB))
            return levpop_this

        maxlev = numpy.where(matrixB > 1e-40)[0]
        if len(maxlev)==0:
            popn = numpy.zeros(len(matrixB))
            return popn
        maxlev=maxlev[-1]
  
        matrixA_in={}
        matrixA_in['init'], matrixA_in['final'], matrixA_in['rate'] = \
            ion.gather_rates(self,  Te, N_e, datacache=datacache, settings=settings,\
                        do_la=True, do_ai=True, do_ec=False, do_pc=False,\
                        do_ir=False)

        i = (matrixA_in['init']<=maxlev) & (matrixA_in['final']<=maxlev)
        matrixA_in['init']=matrixA_in['init'][i]
        matrixA_in['final']=matrixA_in['final'][i]
        matrixA_in['rate']=matrixA_in['rate'][i]

  # fix the rates
        for i in range(len(matrixA_in['init'])):
            if matrixA_in['init'][i]==matrixA_in['final'][i]:
                if matrixA_in['rate'][i] >=0.0:
                    matrixA_in['rate'][i] -=1e10
                    #print("CTieing level %i to ground with rate 1e10"%(i))

        if (maxlev <= const.NLEV_NOSPARSE):
    # convert to a regular solver
            #print("regular solver")
            matrixA = numpy.zeros([maxlev+1,maxlev+1], dtype=float)

            for i in range(len(matrixA_in['init'])):
                matrixA[matrixA_in['final'][i], matrixA_in['init'][i]] += matrixA_in['rate'][i]


    # bug-u-fix
            for i in range(1, maxlev):
                if matrixA[i,i] >= 0:
                    matrixA[i,i]=-1e10
                    #print("FIXING matrixA[%i,%i] = -1.0"%(i,i))

            popn = numpy.zeros(nlev)

            matrixB*=-1

            try:
                popn[1:maxlev] = numpy.linalg.solve(matrixA[1:maxlev,1:maxlev], matrixB[1:maxlev])
            except numpy.linalg.linalg.LinAlgError:
                "EEK ERROR!"
                raise

    
        else:
    # add into sparse solver
            #print("Using sparse solver")
            matrixA={}
            matrixB *= -1
            nlev = len(matrixB)

            if sum(matrixB)>=0:
                return numpy.zeros(len(matrixB))

    # remove ground level
            i = (matrixA_in['init']>0) & (matrixA_in['final']>0)

            matrixA['init'] = matrixA_in['init'][i]
            matrixA['final'] = matrixA_in['final'][i]
            matrixA['rate'] = matrixA_in['rate'][i]

            i = (matrixA['init']<=maxlev+1) & (matrixA['final']<=maxlev+1)

            matrixA['init'] = matrixA['init'][i]
            matrixA['final'] = matrixA['final'][i]
            matrixA['rate'] = matrixA['rate'][i]


    # subtract 1 from the levels
            matrixA['init']-=1
            matrixA['final']-=1

#    A = sparse.coo_matrix((matrixA['rate'],\
#                           (matrixA['final'],matrixA['init'])), \
#                           shape=(maxlev,maxlev)).tocsr()
            A = numpy.zeros([maxlev, maxlev])
            for i in range(len(matrixA['final'])):
                A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]

            popn = numpy.zeros(len(matrixB))
#    popn[1:maxlev+1] = spsolve(A, matrixB[1:maxlev+1])
            popn[1:maxlev+1] = numpy.linalg.solve(A, matrixB[1:maxlev+1])

            popn_bak = popn*1.0

            popn[popn<0] = 0.0
    
  
        
        return popn







#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

  





    def kurucz(self, uin, gam):
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






    def calc_brems_gaunt(self, E, T, z1, brems_type, datacache=False, \
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
        Atomic_data =  ion.datacache(self, Z, z1)
        hdat = Atomic_data['data'][False][False]['hbrems'] 

        #hdat = atomdb.get_data(False, False,'hbrems', settings = settings,\
                                               #datacache=datacache)
        gaunt_D=hdat['BR_GAUNT'].data['COEFFICIENT']

        gamma2 = z0**2 * const.RYDBERG/(const.KBOLTZ*T)

        if ((gamma2 < 1e-3) | (gamma2 > 1e3)):
          if (z0<10):
            print("brems_hummer: Warning, gamma^2 = %e is out of range."%(gamma2))
          gaunt_ff[:]=1.0

        else:
          u = Evec/(const.KBOLTZ*T)
          j =  numpy.where(u <1.e-4)[0]
          if len(j) != 0:
            print("brems_hummer: Warning, u is out of range: ", u[j])
            gaunt_ff[j]=1.0


          gaunt_ff[u>31.6227766]=1.0


      # all out of range data points are set to 1.0 now. Do the good stuff.

          j = numpy.where(gaunt_ff<1.0)[0]
          if len(j) > 0:
            x_u = (2*numpy.log10(u[j])+2.5)/5.5
            x_g = numpy.log10(gamma2)/3.0

            c_j = numpy.zeros(NUM_J, dtype=float)
            for jj in range(NUM_J):
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
          gaunt_ff = ion.kurucz(self, Evec/kT, gam)
        elif (kT==0.0):
          print("brems_kellog: Zero temperature!")
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
#        print("FORCE BORN")
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
        Atomic_data =  ion.datacache(self, Z, z1)
        rdat = Atomic_data['data'][False][False]['rbrems'] 
        #rdat = atomdb.get_data(False, False,'rbrems', settings = settings,\
                                               #datacache=datacache)
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
          Uvec, GauntFFvec = element.extract_gauntff(self, z0, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
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
            Uvec, GauntFFvec = element.extract_gauntff(self,zup, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
          if (zlo>0) & (zup==100):
            Uvec, GauntFFvec = element.extract_gauntff(self,zlo, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
          if (zlo==0) & (zup==100):
          #no match found
              print("brems relativistic: we should never be here")
          if (zlo>0) & (zup<100):
        # we are going to interpolate between the 2 of these
            Uvecl, GauntFFvecl = element.extract_gauntff(self,zlo, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)
            Uvecu, GauntFFvecu = element.extract_gauntff(self,zup, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf)

            if len(numpy.where(numpy.abs(Uvecl-Uvecu)>0.001)[0]) != 0 :
              print("Error: brems_relativistic: U vector mismatch ",  Uvecl, Uvecu)

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
        gaunt_ff[i]=-0.55133*(0.5*numpy.log(10.)*gamma2+numpy.log(10.)*u[i]+0.056745)

    #high energy
        i = numpy.where(u>Uvec[-1])[0]
        gaunt_ff[i]=1.0

    #other energy
        i = numpy.where((u>=Uvec[0]) & (u <=Uvec[-1]))[0]
        gaunt_ff[i] = numpy.interp(u[i],Uvec, GauntFFvec)

        if not Eisvec:
          gaunt_ff = gaunt_ff[0]
        return gaunt_ff
      else:
        print("UNKNOWN BREMS TYPE: ", brems_type)
        return -1


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------







    def do_brems(self, Z, z1, T, abund, brems_type, eedges):
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

      gaunt_ff = ion.calc_brems_gaunt(self, E, T, z1, brems_type)

      emission = const.BREMS_COEFF*abund*\
          (z1-1)*(z1-1)*numpy.exp(-EkT)*(dE/numpy.sqrt(T))*gaunt_ff/(E*const.ERG_KEV)

      return emission






 

    


    

    def do_lines(self,  Te, N_e, datacache=False, settings=False, z1_drv_in=-1, *args, **kwargs):
      #lev_pop=None, 
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
      #print("starting do_lines at %s"%(time.asctime()))
      Z=self.Z
      z1=self.z1
      settings=self.settings
      tstart=time.time()
      
      Atomic_data =  ion.datacache(self, Z, z1)
      #ladat = atomdb.get_data(Z,z1,'LA', datacache=datacache, settings=settings)
      ladat = Atomic_data['data'][Z][z1]['LA']
      lvdat = Atomic_data['data'][Z][z1]['LV']

      #ladat = atomdb.get_data(Z,z1,'LA', datacache=datacache, settings=settings)
      # lvdat = atomdb.get_data(Z,z1,'LV', datacache=datacache, settings=settings)
      ebins = setup.make_vector_nbins(self, settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])
      twoph = numpy.zeros(settings['NumGrid'], dtype=float)
      linelist = numpy.zeros(len(ladat[1].data), \
             dtype= setup.generate_datatypes(self,'linetype'))
      goodlines = numpy.ones(len(linelist), dtype=bool)

      if z1_drv_in < 0:
        z1_drv = z1
      else:
        z1_drv = z1_drv_in

  # now do this

      ## checking if lev_pop is given ##
      if not kwargs:
        lev_pop=ion.level_population(self, Te,N_e)
      else:
        lev_pop = kwargs.get('lev_pop')

      linelist['epsilon'] = ladat[1].data.field('einstein_a') * \
                        lev_pop[ladat[1].data.field('upper_lev')-1]/N_e
      linelist['lambda'] = ladat[1].data.field('wavelen')


      igood = numpy.isfinite(ladat[1].data['wave_obs'])

      if sum(igood) > 0 :
        igood = numpy.where(igood==True)[0]

        igood = igood[numpy.where(ladat[1].data['wave_obs'][igood]>0)[0]]

      if len(igood)  >0:
        linelist['lambda'][igood]=ladat[1].data['wave_obs'][igood]
#  for iline, line in enumerate(ladat[1].data):
#    if numpy.isfinite(line['wave_obs']):
#      if line['wave_obs'] > 0:
#        linelist[iline]['lambda'] = line['wave_obs']
#      else:
#        linelist[iline]['lambda'] = line['wavelen']

#    else:
#      linelist[iline]['lambda'] = line['wavelen']

      linelist['lambda_err'] = numpy.nan
      linelist['epsilon_err'] = numpy.nan
      linelist['element'] = Z
      linelist['elem_drv'] = Z
      linelist['ion'] = z1
      linelist['ion_drv']= z1_drv
      linelist['upperlev'] = ladat[1].data['upper_lev']
      linelist['lowerlev'] = ladat[1].data['lower_lev']

  # I have a linelist. Yay.
      t1=time.time()
      #print("finished making linelist at %s: took %f seconds"%(time.asctime(), t1-tstart))

  # now check for 2 photon transitions
      #print("starting check for 2 photon transitions at %s"%(time.asctime()))
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
      t2=time.time()
      #print("finished checking two photon transitions at %s: took %f seconds"%(time.asctime(), t2-t1))

      linelist = linelist[goodlines]

      tfinish=time.time()
      #print("finished do_lines at %s, took %f seconds"%(time.asctime(), tfinish-tstart))

      return linelist, twoph
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


    def calc_satellite(self, Te, Z=None, z1=None, datacache=False, settings=False):
      """
      Calcaulate DR satellite lines

      Parameters
      ----------
      Z: int
        The nuclear charge of the element  Z: int
      z1 : int
        Recombined Ion charge +1 of ion (e.g. 5 for C VI -> C V)
      Te: float
        The electron temperature (K)
      settings: dictionary
        The settings read from the apec.par file by parse_par_file

      Returns
      -------
      array(linelist)
        List of DR lines
      array(levlistin)
        Rates into each lower level, driven by DR
      """

  # recombining ion charge
      Z=self.Z
      z1=self.z1
      z1_drv=z1+1


      Atomic_data =  ion.datacache(self, Z, z1)
      #print(Atomic_data)
      Atomic_data_drv= ion.datacache(self, Z, z1+1)
      
      
      drdat = Atomic_data_drv['data'][Z][z1_drv]['DR']
      lvdat =  Atomic_data_drv['data'][Z][z1_drv]['LV']
      lvdatrec = Atomic_data['data'][Z][z1]['LV']
      

      

#  pseudocont = numpy.zeros(settings['NumGrid']-1, dtype=float)

#  print "start DR"

      if drdat==False:
        linelist = numpy.zeros(0,dtype=setup.generate_datatypes(self, 'linetype'))
        lev_rates_in = 0.0
      else:

        if not(lvdatrec):
      # no level data for recombined ion.
          lomax = max(drdat[1].data['lowerlev'])
          lev_rates_in=numpy.zeros(lomax+1, dtype=float)

        else:
          lev_rates_in = numpy.zeros(len(lvdatrec[1].data), dtype=float)


        linelist = numpy.zeros(len(drdat[1].data),dtype=setup.generate_datatypes(self, 'linetype'))
        kT = const.KBOLTZ*Te
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
            print("Error in calc_satellite: unknown DR type %i"%\
               (drdat[1].data['type'][iline]))
            epsilon = numpy.nan

      # if there is a translation table from DR SATELLITE to AtomDB levels, use it.
          if len(drdat) > 2:
            itmp = drdat[2].data[drdat[2].data['DRLEVID']==ll]['APEDID']
            if len(itmp) ==1:
              if itmp[0] != 0:
                ll = itmp[0]


          if ll >= len(lev_rates_in):
            print("warning: DR satellite line recombining into non existant level %i of ion Z=%i, z1=%i"%\
              (ll, Z, z1))
          else:
            lev_rates_in[ll-1] += epsilon

          linelist[iline]['lambda'] = lam
          linelist[iline]['lambda_err'] = lamerr
          linelist[iline]['epsilon'] = epsilon
          linelist[iline]['epsilon_err'] = numpy.nan
          linelist[iline]['element'] = Z
          linelist[iline]['ion'] = z1
          linelist[iline]['elem_drv'] = Z
          linelist[iline]['ion_drv'] = z1_drv
          linelist[iline]['upperlev'] = lu
          linelist[iline]['lowerlev'] = ll

        linelist = linelist[numpy.where(linelist['epsilon'] > 0)[0]]


      return linelist, lev_rates_in


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------







    def run_apec_ion(self, Te, N_e, Z=None, z1=None):
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
        
        
        Z=self.Z
        z1=self.z1
        settings=self.settings
        Abund=self.Abund
        
        ionfrac = atomdb._get_precalc_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, Te)
        z1_drv=z1*1
        
        
        settings['WriteIonFname'] = "Z_%i_z1_%i_T_%i_N_%i"%(Z,z1,Te,N_e)
        settings['filemap'] = settings['FileMap']
        settings['atomdbroot'] = os.path.expandvars('$ATOMDB')



  # get the output energy bins for the continuum
        ebins = setup.make_vector_nbins(self, settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])

  ## FIXME CUTOFF FOR MIN IONPOP
        linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
        pseudo = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum = {}
        continuum['brems'] = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum['twophot'] = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum['rrc'] = numpy.zeros(settings['NumGrid'], dtype=float)


  ## FIXME CUTOFF FOR MIN 
        #if ionfrac[z1_drv-1] < const.MIN_IONPOP:
        #    return  linelist, continuum, pseudo
        

  # set up the datacache

        datacache = {}

  # Find the number of levels
        Atomic_data =  ion.datacache(self, Z, z1)
        #print(Atomic_data)
        lvdat = Atomic_data['data'][Z][z1]['LV']
        
        
        linelist_exc = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        linelist_dr = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        linelist_rec = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        if lvdat!= False:
    # find the number of levels
            nlev = len(lvdat[1].data)

    # check if we need to do any of the line-related calculations

            if (settings['EmissionLines'] or settings['TwoPhoton']):


    # gather all the level to level rates
                
                
                up, lo, rates = ion.gather_rates(self, Te, N_e, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
                lev_pop = ion.level_population(self, Te, N_e)

                #print(lev_pop)

    # just in case, add zeros to lengthen the lev_pop appropriately
                if len(lev_pop) < nlev:
                    lev_pop = numpy.append(lev_pop, numpy.zeros(nlev-len(lev_pop), dtype=float))

    # fix any sub-zero level populations
                lev_pop[lev_pop<0] = 0.0
                lev_pop *= Abund[Z]*ionfrac[z1-1]
    
      
                linelist_exc,  continuum['twophot'] = ion.do_lines(self,  Te, N_e, datacache=datacache, settings=settings, z1_drv_in=z1_drv, lev_pop=lev_pop)
                
                print(linelist_exc)
                #linelist_excited= ion.linelist(self, Z, z1, Te, N_e)


            else:
                # skipping the exact level calculation, fill it with zeros, ground state with 1.
                lev_pop = numpy.zeros(nlev, dtype=float)
                lev_pop[0] = 1.0*Abund[Z]*ionfrac[z1-1]

        else:
            lev_pop=numpy.ones(1, dtype=float)*Abund[Z]*ionfrac[z1-1]


        if settings['Bremsstrahlung'] ==True:
            brems = ion.do_brems(self, Z, z1, Te, N_e, settings['BremsType'], ebins)
    # scale for  ion and element abundance.
            continuum['brems']=brems*Abund[Z]*ionfrac[z1-1]

    
        else:
            continuum['brems']=numpy.zeros(len(ebins)-1, dtype=float)
    # now look at the neighbouring ions
        if z1_drv>1:
            z1=z1_drv-1
            if settings['DRSatellite']:
                linelist_dr, drlevrates = ion.calc_satellite(self, Z, z1, Te, datacache=datacache, settings=settings)
                linelist_dr['epsilon']*=Abund[Z]*ionfrac[z1_drv-1]
                drlevrates *=Abund[Z]*ionfrac[z1_drv-1]

            else:
                linelist_dr = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
                drlevrates = 0.0

    # Radiative Recombination
            if settings['RRC']:
      
                rrc, rrlevrates = atomdb.calc_rad_rec_cont(Z, z1, z1_drv, Te, ebins, settings=settings, datacache=datacache)
                continuum['rrc'] = rrc*Abund[Z]*ionfrac[z1_drv-1]
                rrlevrates*=Abund[Z]*ionfrac[z1_drv-1]

            else:
                continuum['rrc'] = numpy.zeros(len(ebins)-1, dtype=float)
                rrlevrates=0.0


    # if there is recombination to process:
            tmpdrlevrates,xxx = util.make_vec(drlevrates)
            tmprrlevrates,xxx = util.make_vec(rrlevrates)

            
            if sum(tmpdrlevrates) + sum(tmprrlevrates)>0:
      

                levpop_recomb= ion.calc_recomb_popn(self, lev_pop, z1_drv, Te, N_e,  drlevrates,\
                                     rrlevrates,\
                                     datacache=datacache, settings=settings)

                


                linelist_rec, tmptwophot =ion.do_lines(self, Te , N_e,  datacache=datacache, settings=settings, z1_drv_in=z1_drv, lev_pop=levpop_recomb)  
                print(linelist_rec)

                continuum['twophot']+= tmptwophot
            
        # now do the ionizing cases
        linelist_ion = numpy.zeros(0,dtype= setup.generate_datatypes(self,'linetype'))
        if z1_drv < Z:
    #datacache={}
            z1=z1_drv+1
            lev_pop_parent = lev_pop*1.0

            while (sum(lev_pop_parent[1:]) > 1e-40) &\
                (z1 <= Z):

                if z1== z1_drv+1:
                    do_xi = True
                else:
                    do_xi = False

                lev_pop = ion.calc_ioniz_popn(self, lev_pop_parent, Z, z1, z1_drv, Te, N_e, \
                                settings=settings, datacache=datacache, \
                                do_xi=do_xi)
                print(lev_pop, len(lev_pop))

                lev_pop[lev_pop<const.MIN_LEVPOP] = 0.0
                if sum(lev_pop[1:]) > 0:
                    linelist_ion_tmp, tmptwophot = \
                        ion.do_lines(self, Te, N_e,  datacache=datacache, settings=settings,  z1_drv_in=z1_drv, lev_pop=lev_pop)

                    linelist_ion = numpy.append(linelist_ion, linelist_ion_tmp)
                    continuum['twophot']+=tmptwophot

                lev_pop_parent = lev_pop
                z1+=1

        linelist = numpy.append(linelist_exc, numpy.append(linelist_dr, numpy.append(linelist_ion, linelist_rec)))
        
        MinEpsilon = settings['MinEpsilon']
        if settings['Ionization']=='CIE':
            MinEpsilon*=0.001

        pseudocont = numpy.zeros(len(ebins)-1, dtype=float)

        if len(linelist) > 0:
            weaklines = linelist[(linelist['epsilon']< MinEpsilon) &\
                         (linelist['lambda']>const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                         (linelist['lambda']<const.HC_IN_KEV_A /settings['GridMinimum'])]

            for line in weaklines:
                e = const.HC_IN_KEV_A /line['lambda']
                ibin = numpy.where(ebins>e)[0][0] - 1
                pseudocont[ibin]+=line['epsilon']

            linelist = linelist[linelist['epsilon'] > MinEpsilon]
        
      
        if settings['WriteIon']==True:
        
            
            ret = {}
            ret['lines'] = linelist
            ret['continuum'] = continuum
            ret['pseudocont'] = pseudocont
            ret['ionfrac'] = ionfrac
            ret['te'] = Te
            ret['dens'] = N_e
            ret['settings'] = settings
            ret['abund'] = Abund[Z]
            fname = settings['WriteIonFname']+'.pkl'
        
            pickle.dump(ret, open(fname, 'wb'))

        print(linelist)
        return linelist, continuum, pseudocont



    




    



    
    def wrap_ion_directly(self, ind, Te, N_e, Z, z1):
        """Changes from apec.py
          Use either ind, or Te and N_e for deriving the pickle file. if ind is set to zero, use the values
          of Te and N_e from sys(argv). If ind is non zero, for example, ind=100, Te and Dens can be derived
          from it

          Other changes: No need to calculate ionfrac and abund here, as ion.run_apec_ion already includes the 
          ionfrac and abundance calculation. Line 4045 to 4074 in apec.py are hence not required. Check CEi and 
          NEI!!! Is this true????

        """
  
  # read in the file
        

        fname = os.path.expandvars('$ATOMDB/apec.par')
        settings = parse_par_file(fname)

        datacache={}
        if ind==0:


          if settings['Ionization']=='NEI':
            z1list = list(range(1, Z+2))
            ionfrac = numpy.ones(len(z1list), dtype=float)

          elif settings['Ionization']=='CIE':
            z1list = list(range(1, Z+2))

    # time to calculate the ionization balance
            if settings['UseIonBalanceTable']:
      # read the ionization balance table
                ionfrac = atomdb._get_precalc_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, Te)
                
            else:
      # calculate the ionization balance
                ionftmp = element.calc_elem_ionbal(self, Z, Te)
                ionfrac = ionftmp

          else:
            print("ERROR: settings['Ionization'] must be CIE or NEI, not %s"%(settings['Ionization']))


          abundfile = atomdb.get_filemap_file('abund',\
                                      Z,\
                                      False,\
                                      fmapfile=settings['FileMap'],\
                                      atomdbroot=os.path.expandvars('$ATOMDB'),\
                                      misc=True)

          abundances = atomdb.get_abundance(abundfile=abundfile, abundset=settings['Abundances'])
  #abundances = atomdb.get_abundance()
          Abund=abundances
          

  # update the output filename
          
          settings['HDUIndex'] = ind
          settings['Te_used'] = Te
          settings['Ne_used'] = N_e

          x=ion.run_apec_ion(self, settings, Z, z1, Te, N_e, ionfrac, Abund)

          ret = {}
          ret['settings'] = settings
          ret['data'] = x
          settings['WriteIonFname'] ="Z_%i_z1_%i_iT_%iiN_%i.pkl"%(Z,z1,Te,N_e)
          fname_pickle = settings['OutputFileStem']+'_'+settings['WriteIonFname']
          pickle.dump(ret, open(fname_pickle,'wb'))
          print("wrote file %s"%(fname_pickle))
          print("Finished cleanly at %s"%(time.asctime()))


        
        else:
          te = setup.make_vector(self,settings['LinearTemp'], \
                   settings['TempStart'], \
                   settings['TempStep'], \
                   settings['NumTemp'])

          if settings['TempUnits']=='keV':
            te /= const.KBOLTZ


          dens = setup.make_vector(self,settings['LinearDens'], \
                     settings['DensStart'], \
                     settings['DensStep'], \
                     settings['NumDens'])

          ite = ind //len(dens)
          idens = ind%len(dens)
          
          Te = te[ite]
          Dens = dens[idens]

          if settings['Ionization']=='NEI':
            z1list = list(range(1, Z+2))
            ionfrac = numpy.ones(len(z1list), dtype=float)

          elif settings['Ionization']=='CIE':
            z1list = list(range(1, Z+2))

    # time to calculate the ionization balance
            if settings['UseIonBalanceTable']:
      # read the ionization balance table
              ionfrac = atomdb._get_precalc_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, te)
            else:
      # calculate the ionization balance
              ionftmp = element.calc_elem_ionbal(self, Z, Te)
              ionfrac = ionftmp

          else:
            print("ERROR: settings['Ionization'] must be CIE or NEI, not %s"%(settings['Ionization']))


          abundfile = atomdb.get_filemap_file('abund',\
                                      Z,\
                                      False,\
                                      fmapfile=settings['FileMap'],\
                                      atomdbroot=os.path.expandvars('$ATOMDB'),\
                                      misc=True)

          abundances = atomdb.get_abundance(abundfile=abundfile, abundset=settings['Abundances'])
  #abundances = atomdb.get_abundance()
          Abund=abundances
          print("Z_%i_z1_%i_iT_%iiN_%i.pkl"%(Z,z1,ite,idens))
  # update the output filename
          #settings['WriteIonFname'] ="Z_%i_z1_%i_iT_%iiN_%i.pkl"%(Z,z1,ite,idens)
          settings['HDUIndex'] = ind
          settings['Te_used'] = Te
          settings['Ne_used'] = Dens

          x = ion.run_apec_ion(self, settings, Z, z1, Te, Dens,  ionfrac, Abund)
          
          

          ret = {}
          ret['settings'] = settings
          ret['data'] = x
          settings['WriteIonFname'] ="Z_%i_z1_%i_iT_%iiN_%i.pkl"%(Z,z1,ite,idens)
          fname_pickle = settings['OutputFileStem']+'_'+settings['WriteIonFname']
          pickle.dump(ret, open(fname_pickle,'wb'))
          print("wrote file %s"%(fname_pickle))
          print("Finished cleanly at %s"%(time.asctime()))

  







    


    """writing the level_population in a file
    """
    def writeout_data_levelpop(self,Z, z1, Te, N_e, settings):

        
        f_levpop = open("apec_levpop.txt", "w+")
        #lev_pop_write=ion.level_pop
        lev_pop_write = ion.level_population(self, Z, z1, Te, N_e, settings)
        for element in lev_pop_write:
            f_levpop.write('%e\n'%(element))
        f_levpop.close()

        return



    


    




    """
    Convert level populations into line lists. Similar to do_lines in apec.py. The main differences are"
    1) Atomic data (ladat, lvdat) is being read from the stored datacache
    2) Updated line emissivities that reflect change in emissivity due to change/uncertainty in atomic data

    Parameters
    ----------
    Z: int
    The nuclear charge of the element
    z1 : int
    Ion charge
    Te : float
    Temperature in K
    N_e : float
    Electron Density (cm^-3)
    settings : dict
    See description in atomdb.get_data
    z1_drv_in : int
    the driving ion for this calculation, if not z1 (defaults to z1)

    Returns: linelist
    The list of lines and their emissivities.
    """
    

                                                              
   

    """writing the linelist in a file
    """
    def write_linelist(self,Z, z1, Te, N_e):
        f_linelist = open("apec_linelist.txt", "w+")
        line_emissivities = ion.run_apec_ion(self, Z, z1, Te, N_e)
        for t in line_emissivities:
            line = ',   '.join(str(x) for x in t)
            f_linelist.write((line) +"\n")
        f_linelist.close()

        return






    


    """
    Convert level population into emissivity of invidual lines adjusted with factor   
    Parameters
    ----------
    Z: int
     The nuclear charge of the element
    z1 : int
     Ion charge
    Te : float
     Temperature in K
    N_e : float
     Electron Density (cm^-3)
    upper_level:
     Electron is making transition from this level
    lower)level:
     Electron is making transition to this level
    factor: float
     this determines which atomic data (e.g. collision strength) is adjusted by what percentage. 
     factor = 1.1 means 10% increase in the atomic     data 
    settings : dict
     See description in atomdb.get_data
    z1_drv_in : int
     the driving ion for this calculation, if not z1 (defaults to z1)

    
    Returns: line emissivity
        Line emissivity for an individual line. Sensitive to change in atomic data.
    """
    

    '''
    def single_emissivity(self, Z, z1, Te, N_e, upper_level, lower_level, factor, settings=False, z1_drv_in=-1):

        
        datacache={}
        upper_level = self.upper_level
        lower_level =self.lower_level
        factor =self.factor
        Atomic_data =  ion.datacache(self, Z, z1)

        ladat = Atomic_data['data'][Z][z1]['LA']
        lvdat =  Atomic_data['data'][Z][z1]['LV']

        linelist = numpy.zeros(len(ladat[1].data), \
             dtype= apec.generate_datatypes('linetype'))
        goodlines = numpy.ones(len(linelist), dtype=bool)
        

        if z1_drv_in < 0:
            z1_drv = z1
        else:
            z1_drv = z1_drv_in

        
        lev_pop = ion.level_population(self, Z, z1, Te, N_e,  settings)


        LA =  ladat[1].data.field('einstein_a')
        ULV =  ladat[1].data.field('upper_lev')
        LLV = ladat[1].data.field('lower_lev')
        WV  = ladat[1].data.field('wavelen')
        WV_OBS = ladat[1].data.field('wave_obs')

        #for k in range(len(LLV)):
        #    if type(WV_OBS) == 'nan':
         #       Wavelength = WV[k]
          #  else:
           #     Wavelength = WV_OBS[k]


        #ion.enablePrint()
        for i in range(len(LLV)):
            if ULV[i]==upper_level:
                if LLV[i]==lower_level:
                    if self.errortype=='LA':
                        Emissivity = factor*LA[i]* \
                                lev_pop[ULV[i]-1]/N_e


                        if math.isnan(WV_OBS[i])==True:
                            Wavelength = WV[i]
                        else:
                            Wavelength = WV_OBS[i]

                        Wavelength_err = numpy.nan
                        Emissivity_err = numpy.nan


                    else:
                        Emissivity = LA[i]* \
                                lev_pop[ULV[i]-1]/N_e



                        if math.isnan(WV_OBS[i])==True:
                            Wavelength = WV[i]
                        else:
                            Wavelength = WV_OBS[i]
        
                        Wavelength_err = numpy.nan
                        Emissivity_err = numpy.nan

        return Wavelength, Wavelength_err, Emissivity, Emissivity_err, Z, z1, upper_level, lower_level


        '''






    '''

    """
    Reading emissivity from apec_nei_line.fits and scaling with abundances from "AG89"
    Parameters
    ----------
    Z: int
     The nuclear charge of the element
    z1 : int
     Ion charge
    upper_level:
     Electron is making transition from this level
    lower)level:
     Electron is making transition to this level
    z1_drv : int
     the driving ion for this calculation, if not z1 (defaults to z1)

    
    Returns: Temperatures in kelvin and corresponding emissivities for a single transition
    """
    def return_emissivity(self, Z, z1, upper_level, lower_level, z1_drv=-1):
    
    
        if z1_drv == -1:
            z1_drv=z1


        linefile= os.path.expandvars("$ATOMDB/apec_nei_line.fits")
        a = pyfits.open(linefile)
        kTlist = a[1].data['KT']
        kTlist_Kelvin= kTlist/(const.KBOLTZ)
        emiss = numpy.zeros(len(kTlist))

        Abund = atomdb.get_abundance(abundfile=False, \
                abundset='AG89', element=[-1], datacache=False, settings=False, show=False)

        for i in range(len(kTlist)):
            ihdu = i+2
            d = a[ihdu].data
            l = d[(d['element']== Z) &\
                 (d['ion'] == z1) &\
                 (d['ion_drv'] == z1_drv) &\
                 (d['upperlev'] == upper_level) &\
                 (d['lowerlev'] == lower_level)]

            for ll in l:
                emiss[i] += ll['epsilon']/Abund[self.Z]
                
        return kTlist_Kelvin, emiss



        '''





class variableapec():




    def __init__(self, Z, z1, Te, N_e, uls, lls, factor, errortype, upper_level, lower_level, ind, Abund, settings, ionfrac):
        '''
        Z=self.Z
        z1=self.z1
        Te=self.Te
        N_e=self.N_e
        factor=self.factor
        errortype=self.errortype
        ind=self.ind
        '''
        ionfrac= element.calc_elem_ionbal_delta(self, Z, Te, factor, errortype)
        
        Abund= atomdb.get_abundance(abundfile=False, \
                abundset='AG89', element=[-1], datacache=False, settings=False, show=False)
        fname = os.path.expandvars('$ATOMDB/apec.par')
        settings = parse_par_file(fname)








    def gather_rates_delta(self, Z, z1, Te, N_e, factors, errortypes, upper_levels, lower_levels, datacache=True, settings=False,\
                 do_la=True, do_ai=True, do_ec=True, do_pc=True,\
                 do_ir=True):
      """
      fetch the rates for all the levels of Z, z1.  'IR', 'LV', 'LA', 'EC', 'PC', 'DR', 'PI', 'AI' rates are being read from datacache.
      Parameters
      ----------
      Z: int
        The nuclear charge of the element
      z1 : int
        ion charge +1
      Te : float
        temperture (Kelvin)
      N_e: float
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
      larate_new = []
      
      pcrate_new = []
      airate_new = []
      ecrate_new=[]
      irrate_new=[]
      
      change_num_la=0
      change_num_ec=0
      change_num_pc=0
      change_num_ai=0
      change_num_ir=0
      #factor_AI = 1.0
      #factor_EC = 1.0
      #factor_IR = 1.0
      #factor_LA = 1.0

      Atomic_data =  ion.datacache(self, Z, z1)
      



      #for fac in range(len(factors)):
      # factor = factors[fac]



      for err in range(len(errortypes)):

        errortype = errortypes[err]

        if errortype == 'LA':
          larate_new = []
          factor_LA = factors[err][:]

        if errortype == 'EC':
          ecrate_new = []
          factor_EC = factors[err][:]

        if errortype == 'PC':
          pcrate_new = []
          factor_PC = factors[err][:]

        if errortype == 'AI':
          airate_new = []
          factor_AI = factors[err][:]

        if errortype == 'IR':
          irrate_new = []
          factor_IR = factors[err][:]




        #for num in range(len(upper_levels)):
        #  upper_level = int(upper_levels[num])
        #  lower_level = int(lower_levels[num])

        Te_arr, dummy = util.make_vec(Te)



        lvdat =  Atomic_data['data'][Z][z1]['LV']


        nlev = len(lvdat[1].data)






        diagterms = numpy.zeros(nlev)
  # get the LA data:
        if 'AAUT_TOT' in lvdat[1].data.names:
          has_sum_lv = True
        else:
          has_sum_lv = False
      

        laup = numpy.zeros(0, dtype=int)
        lalo = numpy.zeros(0, dtype=int)
        larate = numpy.zeros(0, dtype=float)

        if do_la:
          t1=time.time()
          if has_sum_lv:
            diagterms+=lvdat[1].data['ARAD_TOT']

        

          ladat =  Atomic_data['data'][Z][z1]['LA']
  

          if ladat != False:

            laup = ladat[1].data['Upper_Lev'] - 1
            lalo = ladat[1].data['Lower_Lev'] - 1
            larate = ladat[1].data['Einstein_A']



            get_indecies = lambda xx, xxs: [iii for (yy, iii) in zip(xxs, range(len(xxs))) if xx == yy]
            index = 0.1     ##### initialize a value for index

            for num in range(len(upper_levels)):
              upper_level = int(upper_levels[num])
              lower_level = int(lower_levels[num])


              
              index_up = get_indecies(upper_level-1,laup)
              index_lo = get_indecies(lower_level-1,lalo)

              for i in range(len(index_up)):
                if index_up[i] in index_lo:
                  index = index_up[i]
                  break
              
              if errortype =='LA':
                if index == 0.1:   ### checking if transition exists
                  larate = larate
                else:
                
                  larate[index] = larate[index] * factor_LA[num]
                  change_num_la = index
            #print(factor_LA)


          #now multiply  larate with factor if errortype = 'LA'


            if not(has_sum_lv):
#         diagterms[laup] +=larate
              for i in range(len(laup)):
                diagterms[laup[i]] +=larate[i]
          t2 = time.time()
        


      # create dummy results

  # get the AI data:
        aiup = numpy.zeros(0, dtype=int)
        ailo = numpy.zeros(0, dtype=int)
        airate = numpy.zeros(0, dtype=float)

        if do_ai:
          if has_sum_lv:
            diagterms+= lvdat[1].data['AAUT_TOT']


          t1=time.time()
          aidat =  Atomic_data['data'][Z][z1]['AI']

        
          if aidat != False:
            aiup = aidat[1].data['Level_Init'][:] - 1
            ailo = numpy.zeros(len(aidat[1].data), dtype=int)
            airate = aidat[1].data['Auto_Rate']



          index = 0.1 

          for num in range(len(upper_levels)):
              upper_level = int(upper_levels[num])
              lower_level = int(lower_levels[num])


              
              index_up = get_indecies(upper_level-1,aiup)
              index_lo = get_indecies(lower_level-1,ailo)

              for i in range(len(index_up)):
                if index_up[i] in index_lo:
                  index = index_up[i]
                  break
              
              if errortype =='AI':
                if index == 0.1:   ### checking if transition exists
                  airate = airate
                else:
                  airate[index] = airate[index] * factor_AI[num]
                  change_num_ai = index
           


          if not(has_sum_lv):
            for i in range(len(aiup)):
                diagterms[aiup[i]] +=airate[i]
          t2=time.time()
        

  # get the EC data:

        ecup = numpy.zeros(0, dtype=int)
        eclo = numpy.zeros(0, dtype=int)
        ecrate = numpy.zeros(0, dtype=float)

        if do_ec:
          t1 = time.time()
          ecdat =  Atomic_data['data'][Z][z1]['EC']
        
          if ecdat != False:

            lvdat =  Atomic_data['data'][Z][z1]['LV']
          

            ecup = numpy.zeros(len(ecdat[1].data), dtype=int)
            eclo = numpy.zeros(len(ecdat[1].data), dtype=int)
            ecrate = numpy.zeros(len(ecdat[1].data), dtype=float)

            decup = numpy.zeros(len(ecdat[1].data), dtype=int)
            declo = numpy.zeros(len(ecdat[1].data), dtype=int)
            decrate = numpy.zeros(len(ecdat[1].data), dtype=float)
            idex = 0
    # need to loop and calculate each result

#     Te_arr = numpy.array([te])
            deglarr = lvdat[1].data['LEV_DEG'][ecdat[1].data['lower_lev']-1]
            deguarr = lvdat[1].data['LEV_DEG'][ecdat[1].data['upper_lev']-1]
            deltaearr = lvdat[1].data['ENERGY'][ecdat[1].data['upper_lev']-1]-\
                  lvdat[1].data['ENERGY'][ecdat[1].data['lower_lev']-1]

            ladat =  Atomic_data['data'][Z][z1]['LA']

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




              exc,dex, tmp = atomdb._calc_maxwell_rates(ecdat[1].data['coeff_type'][i],\
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
              ecrate[i] = exc*N_e

              if dex > 0:
                decup[idex] = iup
                declo[idex] = ilo
                decrate[idex] = dex*N_e
                idex += 1

    

            index = 0.1 
            for num in range(len(upper_levels)):
              upper_level = int(upper_levels[num])
              lower_level = int(lower_levels[num])


              
              index_up = get_indecies(upper_level-1,eclo)
              index_lo = get_indecies(lower_level-1,ecup)

              for i in range(len(index_lo)):
                if index_lo[i] in index_up:
                  index = index_lo[i]
                  break
              
              if errortype =='EC':
                if index == 0.1:   ### checking if transition exists
                  ecrate = ecrate
                else:
                  ecrate[index] = ecrate[index] * factor_EC[num]
                  change_num_ec = index
                
                print(index)

    #now multiply the excitation and deexcitation rate with factor if errortype = 'EC'
    # now merge the de-excitation data
            decup = decup[:idex]
            declo = declo[:idex]
            decrate = decrate[:idex]

            for j in range(len(declo)):
              if decup[j] == upper_level-1:
                if declo[j] == lower_level -1:
                      if errortype=='EC':
                          decrate[j] = decrate[j]
                            #change_num=j

            ecup = numpy.append(ecup, decup)
            eclo = numpy.append(eclo, declo)
            ecrate = numpy.append(ecrate, decrate)

    # create dummy results
#       if not(has_sum_lv):
            for i in range(len(ecup)):
              diagterms[ecup[i]] +=ecrate[i]
          t2 = time.time()
        

  # get the PC data:
        pcup = numpy.zeros(0, dtype=int)
        pclo = numpy.zeros(0, dtype=int)
        pcrate = numpy.zeros(0, dtype=float)

        if do_pc:
          t1 = time.time()
          pcdat =  Atomic_data['data'][Z][z1]['PC']
          if pcdat != False:
            pcup = numpy.zeros(len(pcdat[1].data), dtype=int)
            pclo = numpy.zeros(len(pcdat[1].data), dtype=int)
            pcrate = numpy.zeros(len(pcdat[1].data), dtype=float)

            dpcup = numpy.zeros(len(pcdat[1].data), dtype=int)
            dpclo = numpy.zeros(len(pcdat[1].data), dtype=int)
            dpcrate = numpy.zeros(len(pcdat[1].data), dtype=float)
            idex = 0
    # need to loop and calculate each result

#     Te_arr = numpy.array([te])
            deglarr = lvdat[1].data['LEV_DEG'][pcdat[1].data['lower_lev']-1]
            deguarr = lvdat[1].data['LEV_DEG'][pcdat[1].data['upper_lev']-1]
            deltaearr = lvdat[1].data['ENERGY'][pcdat[1].data['upper_lev']-1]-\
                    lvdat[1].data['ENERGY'][pcdat[1].data['lower_lev']-1]



            for i in range(len(pcdat[1].data)):


              exc,dex, tmp = atomdb._calc_maxwell_rates(pcdat[1].data['coeff_type'][i],\
                                     pcdat[1].data['min_temp'][i],\
                                     pcdat[1].data['max_temp'][i],\
                                     pcdat[1].data['temperature'][i],\
                                     pcdat[1].data['effcollstrpar'][i],\
                                     deltaearr[i]/1e3, Te_arr, Z, \
                                     deglarr[i], deguarr[i],\
                                     force_extrap=True)


              pcup[i] = pcdat[1].data['Lower_Lev'][i] - 1
              pclo[i] = pcdat[1].data['Upper_Lev'][i] - 1
              pcrate[i] = exc*N_e

              if dex > 0:
                dpcup[idex] = pcdat[1].data['Upper_Lev'][i] - 1
                dpclo[idex] = pcdat[1].data['Lower_Lev'][i] - 1
                dpcrate[idex] = dex*N_e
                idex += 1

          # now multiply the excitation and deexcitation rate with factor if errortype = 'PC' 
          #merge the de-excitation data

            index = 0.1 
            for num in range(len(upper_levels)):
              upper_level = int(upper_levels[num])
              lower_level = int(lower_levels[num])


              
              index_up = get_indecies(upper_level-1,pclo)
              index_lo = get_indecies(lower_level-1,pcup)

              for i in range(len(index_lo)):
                if index_lo[i] in index_up:
                  index = index_lo[i]
                  break
              
              if errortype =='PC':
                if index == 0.1:   ### checking if transition exists
                  pcrate = pcrate
                else:
                  pcrate[index] = pcrate[index] * factor_PC[num]
                  change_num_pc = index
          

    # now merge the de-excitation data
            dpcup = dpcup[:idex]
            dpclo = dpclo[:idex]
            dpcrate = dpcrate[:idex]


            for l in range(len(dpclo)):
              if dpcup[l]==upper_level-1:
                if dpclo[l]==lower_level-1:
                      if errortype=='PC':
                          dpcrate[l] = dpcrate[l]



            pcup = numpy.append(pcup, dpcup)
            pclo = numpy.append(pclo, dpclo)
            pcrate = numpy.append(pcrate, dpcrate)

#       if not(has_sum_lv):
            for i in range(len(pcup)):
              diagterms[pcup[i]] +=pcrate[i]
          t2=time.time()
        

  # get the IR data for colln ionization:
        irup = numpy.zeros(0, dtype=int)
        irlo = numpy.zeros(0, dtype=int)
        irrate = numpy.zeros(0, dtype=float)
        if do_ir:
          t1 = time.time()
          irdat =  Atomic_data['data'][Z][z1]['IR']

        
          if irdat != False:
            irup = numpy.zeros(len(irdat[1].data), dtype=int)
            irlo = numpy.zeros(len(irdat[1].data), dtype=int)
            irrate = numpy.zeros(len(irdat[1].data), dtype=float)

            iir = 0

              
      # need to loop and calculate each result

            Te_arr = numpy.array(Te)
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
                irrate[iir] = rate*N_e
                iir += 1



              
                

              


      # now merge the de-excitation data
            irup = irup[:iir]
            irlo = irlo[:iir]
            irrate = irrate[:iir]
            index = 0.1 

            for num in range(len(upper_levels)):
              upper_level = int(upper_levels[num])
              lower_level = int(lower_levels[num])


              
              index_up = get_indecies(upper_level-1,irup)
              index_lo = get_indecies(lower_level-1,irlo)

              for i in range(len(index_up)):
                if index_up[i] in index_lo:
                  index = index_up[i]
                  break
              
              if errortype =='IR':
                if index == 0.1:   ### checking if transition exists
                  irrate = irrate
                else:
                  irrate[index] = irrate[index] * factor_IR[num]
                  change_num_ir = index
                  
              

              
              
#       for i in range(len(irup)):
#         print "IR %i %i = %e"%(irup[i], irlo[i], irrate[i])
#       if not(has_sum_lv):
            for i in range(len(irup)):
              diagterms[irup[i]] +=irrate[i]
          t2=time.time()
        

        tmp={}
        tmp['up']={}
        tmp['lo']={}
        tmp['rate']={}
        tmp['up']['ec'] = ecup
        tmp['lo']['ec'] = eclo
        tmp['rate']['ec'] = ecrate
         
        

        tmp['up']['pc'] = pcup
        tmp['lo']['pc'] = pclo
        tmp['rate']['pc'] = pcrate
          

        tmp['up']['ir'] = irup
        tmp['lo']['ir'] = irlo
        tmp['rate']['ir'] = irrate

          #print(irup,irlo)
          

        tmp['up']['ai'] = aiup
        tmp['lo']['ai'] = ailo
        tmp['rate']['ai'] = airate
          

        tmp['up']['la'] = laup
        tmp['lo']['la'] = lalo
        tmp['rate']['la'] = larate
          #print(laup,lalo)
          

        tmp['up']['diag'] = numpy.arange(nlev, dtype=int)
        tmp['lo']['diag'] = numpy.arange(nlev, dtype=int)
        tmp['rate']['diag'] = diagterms*-1

#   pickle.dump(tmp, open('tmp_gather_%i_%i.pkl'%(Z,z1),'wb'))



        t1=time.time()
        up_out = numpy.append(laup, numpy.append(aiup, numpy.append(ecup, numpy.append(pcup, irup))))
        lo_out = numpy.append(lalo, numpy.append(ailo, numpy.append(eclo, numpy.append(pclo, irlo))))
        

        if len(ecrate_new) == 0:
          ecrate_new = ecrate
        if errortype == 'EC':
          if len(ecrate) != 0:
            
            ecrate_new[change_num_ec] = ecrate[change_num_ec]
              #ecrate_new[change_num_ec] = factor_EC * ecrate_new[change_num_ec]
            
          
          #rate_out = numpy.append(larate, numpy.append(airate, numpy.append(new_array, numpy.append(pcrate, irrate))))
          #print(change_num)

        if len(larate_new) == 0:
          larate_new = larate
        if errortype == 'LA':
          if len(larate) != 0:
            
            larate_new[change_num_la] = larate[change_num_la]
              #larate_new[change_num_la] = factor_LA * larate_new[change_num_la]
            
          
          #rate_out = numpy.append(new_array, numpy.append(airate, numpy.append(ecrate, numpy.append(pcrate, irrate))))

        
        if len(pcrate_new) == 0:
          pcrate_new = pcrate
        if errortype == 'PC':
          if len(pcrate) != 0:
            pcrate_new[change_num_pc] = pcrate[change_num_pc]
              #pcrate_new[change_num_pc] = factor_PC * pcrate_new[change_num_pc]
            
          

          #rate_out = numpy.append(larate, numpy.append(airate, numpy.append(new_array, numpy.append(pcrate, irrate))))
          #print(change_num)

        
        if len(airate_new) == 0:
          airate_new = airate
        if errortype == 'AI':
          if len(airate) != 0:
            
            airate_new[change_num_ai] = airate[change_num_ai]
              #airate_new[change_num_ai] = factor_AI * airate_new[change_num_ai]


        if len(irrate_new) == 0:
          irrate_new = irrate
        if errortype == 'IR':
          if len(irrate) != 0:
            
            irrate_new[change_num_ir] = irrate[change_num_ir]
              #irrate_new[change_num_ir] = factor_IR * irrate_new[change_num_ir]

            
            
          

          #rate_out = numpy.append(larate, numpy.append(airate, numpy.append(new_array, numpy.append(pcrate, irrate))))
          #print(change_num)


        #rate_out = numpy.append(larate, numpy.append(airate, numpy.append(ecrate, numpy.append(pcrate, irrate))))

        up_out = numpy.append(up_out, numpy.arange(nlev, dtype=int))
        lo_out = numpy.append(lo_out, numpy.arange(nlev, dtype=int))


      #print(new_array)


      '''
      if errortype == 'EC':
        rate_out = numpy.append(larate, numpy.append(airate, numpy.append(new_array, numpy.append(pcrate, irrate))))
        rate_out = numpy.append(rate_out, diagterms*-1)

      if errortype == 'LA':
        rate_out = numpy.append(new_array, numpy.append(airate, numpy.append(ecrate, numpy.append(pcrate, irrate))))
        rate_out = numpy.append(rate_out, diagterms*-1)

      if errortype == 'PC':
        rate_out = numpy.append(larate, numpy.append(airate, numpy.append(ecrate, numpy.append(new_array, irrate))))
        rate_out = numpy.append(rate_out, diagterms*-1)

      if errortype == 'AI':
        rate_out = numpy.append(larate, numpy.append(new_array, numpy.append(ecrate, numpy.append(pcrate, irrate))))
        rate_out = numpy.append(rate_out, diagterms*-1)
      '''

      
      print(ecrate_new)
      #print(larate_new)
      rate_out = numpy.append(larate_new, numpy.append(airate_new, numpy.append(ecrate_new, numpy.append(pcrate_new, irrate_new))))
      rate_out = numpy.append(rate_out, diagterms*-1)
      t2=time.time()

        #if len(new_array) == 0:
        #  new_array = rate_out
        
        #if numpy.array_equal(new_array,rate_out) == False:
        #  for ran in range(len(rate_out)):
        #    if new_array[ran] != rate_out[ran]:
        #new_array[change_num] = rate_out[change_num]

      




      
      return up_out, lo_out, rate_out









    
    
    


    





    





    def calc_elem_ionbal_delta(self, Z, Te, factor, errortype,  tau=False, init_pop='ionizing', teunit='K',\
                    extrap=True, settings=False, datacache=False):
      """
      Calculate the ionization balance for all the elements in Zlist.

      One of init_pop or Te_init should be set. If neither is set, assume
      all elements start from neutral.

c
      Parameters
      ----------
      Z : int
        nuclear charge to include in calculation (e.g. 8 for oxygen)
      Te : float
        electron temperature in keV or K (default K)
      tau : float
        N_e * t for the non-equilibrium ioniziation (default False, i.e. CIE)
      init_pop : string or float
        if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
        if 'recombining': all recombining from ionized (so[...0,0,1])
        if array: actual populations (e.g. [0, 0.1, 0.3, 0.5, 0.1, 0, 0])
        if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
        if single float : the temperature (same units as Te)
      teunit : {'K' , 'keV'}
        units of temperatures (default K)
      extrap : bool
        Extrappolate rates to values outside their given range. (default False)

      Returns
      -------
      final_pop : array
        final population. E.g. [0.1,0.2,0.3,0.2,0.2,0.0,0.0]

      """

      kT = util.convert_temp(Te, teunit, 'keV')
      if tau==False:
        cie = True
        init_pop_calc=False
      else:
        cie = False
      if not cie:
          # if it's not equilibrium, get the initial population
        if isinstance(init_pop, str):
          if init_pop.lower() == 'ionizing':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[0] = 1.0
          elif init_pop.lower() == 'recombining':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[-1] = 1.0
          else:
            raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
             (init_pop))
        elif isinstance(init_pop, float):
      # this is an initial temperature
          kT_init = util.convert_temp(init_pop, teunit, 'keV')

      # rerun this routine in equilibrium mode to find the initial ion pop
          init_pop_calc = variableapec.return_ionbal_delta(self,Z, kT_init, factor, errortype, \
                                    teunit='keV', \
                                    datacache=datacache,fast=False,
                                    settings = settings, extrap=extrap)


        elif isinstance(init_pop, numpy.ndarray) or isinstance(init_pop, list):
          init_pop_calc = init_pop
        elif isinstance(init_pop, dict):
          init_pop_calc = init_pop[Z]
        else:
          raise util.OptionError("Error: invalid type for init_pop")


  # get the end point population

      ionrate = numpy.zeros(Z, dtype=float)
      recrate = numpy.zeros(Z, dtype=float)



      for z1 in range(1,Z+1):
            if errortype=='ION': 
                ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)
                ionrate[z1-1], recrate[z1-1] = factor*ionrate[z1-1], recrate[z1-1]
                #print(ionrate[z1-1])

            elif errortype=='REC': 
                ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)
                ionrate[z1-1], recrate[z1-1]= ionrate[z1-1], factor*recrate[z1-1]
                

            else:
                ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)

      if cie:
        final_pop = setup.solve_ionbal(self,ionrate, recrate)
      else:
        final_pop = setup.solve_ionbal(self, ionrate, recrate, init_pop=init_pop_calc, tau=tau)

      return final_pop




    def return_ionbal_delta(self, Z, Te, factor, errortype, init_pop=False, tau=False,\
                       teunit='K', \
                       filename=False, datacache=False, fast=True,
                       settings= False, debug=False, extrap=True):

      """
      Solve the ionization balance for a element Z.

      Parameters
      ----------
      Z : int
        atomic number of element
      Te : float or array
        electron temperature(s), default in K
      init_pop : float array
        initial population of ions for non-equlibrium calculations. Will be renormalised to 1.
      tau : float or array
        N_e * t for the non-equilibrium ioniziation, in cm^3 s^-1.
      Te_init : float
        initial ionization balance temperature, same units as Te
      teunit : {'K' , 'keV'}
        units of temperatures (default K)
      filename : string
        Can optionally point directly to the file in question, i.e. to look at older data
        look at $HEADAS/../spectral/modelData/eigenELSYMB_v3.0.fits.
        If not set, download from AtomDB FTP site.
      datacache : dict
        Used for caching the data. See description in atomdb.get_data
      fast : bool
        If true, use precalculated eigenvector files to obtain CIE and NEI results

      Returns
      -------
      final_pop : float array
        final populations.

      """

      if fast:

        
        ionbal = variableapec._solve_ionbal_eigen_delta(self, Z, Te, factor, errortype, init_pop=init_pop, tau=tau, \
                       teunit=teunit, \
                       filename=filename, datacache=datacache, debug=debug)
        return ionbal

      else:
        ionbal = variableapec.calc_elem_ionbal_delta(self, Z, Te, factor, errortype, tau=tau, init_pop=init_pop, teunit=teunit,\
                               extrap=extrap, settings=settings, datacache=datacache)
        return ionbal




    def _solve_ionbal_eigen_delta(self, Z, Te, factor, errortype, init_pop=False, tau=False, \
                       teunit='K', \
                       filename=False, datacache=False, debug=False):
      """
      Solve the ionization balance for element Z using the eigenvector
      approach and files as distributed in XSPEC.

      Parameters
      ----------
      Z : int
        atomic number of element
      Te : float or array
        electron temperature(s), default in K
      init_pop : float array
        initial population of ions for non-equlibrium calculations. Will be renormalised to 1.
      tau : float or array
        N_e * t for the non-equilibrium ioniziation, in cm^3 s^-1.

      teunit : {'K' , 'keV'}
        units of temperatures (default K)
      filename : string
        Can optionally point directly to the file in question, i.e. to look at older data
        look at $HEADAS/../spectral/modelData/eigenELSYMB_v3.0.fits.
        If not set, download from AtomDB FTP site.
      datacache : dict
        Used for caching the data. See description in atomdb.get_data

      Returns
      -------
      final_pop : float array
        final populations.

      """
    #
    #  Version 0.1 Initial Release
    #  Adam Foster 16th September 2015
    #


      kT = util.convert_temp(Te, teunit, 'keV')

      if type(tau)==bool:
        if tau==False:
          cie = True
          init_pop_calc=False

        else:
          raise ValueError("Error: tau should be False, a float, or an array of floats. Received "+repr(Tau))
      else:
        cie = False

      if not cie:
      # if it's not equilibrium, get the initial population
        if isinstance(init_pop, str):
          if init_pop.lower() == 'ionizing':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[0] = 1.0
          elif init_pop.lower() == 'recombining':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[-1] = 1.0
          else:
            raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
             (init_pop))
        elif isinstance(init_pop, float):
      # this is an initial temperature
          kT_init = util.convert_temp(init_pop, teunit, 'keV')
          init_pop_calc = variableapec.return_ionbal_delta(self, Z, kT_init, factor, errortype,\
                                            teunit='keV', \
                                            datacache=datacache,fast=True)


        elif isinstance(init_pop, numpy.ndarray) or isinstance(init_pop, list):
          init_pop_calc = init_pop
        elif isinstance(init_pop, dict):
          init_pop_calc = init_pop[Z]
        else:
          raise util.OptionError("Error: invalid type for init_pop")


  # open the eigenvector data file

      if util.keyword_check(filename):
    # we have a filename specified!
        fname = os.path.expandvars(filename)
        if not os.path.isfile(fname):
          print("Specified file %s does not exist. Exiting"%(fname))
          return
        d = pyfits.open(fname)
      else:
        Atomic_data =  ion.datacache(self, Z, z1)
        d = Atomic_data['data'][Z][False]['eigen'] 
        #d = atomdb.get_data(Z, False, 'eigen', datacache=datacache)
      telist = numpy.logspace(4,9,1251)
      kTlist=telist*const.KBOLTZ

  # if we are looking for equilibrium, return the nearest data
      if cie:

        ikTlist = numpy.argsort(numpy.abs(kTlist-kT))
        ite = [min(ikTlist[:2]), max(ikTlist[:2])]
        Tdiff = kTlist[ite[1]] - kTlist[ite[0]]
        if Tdiff > 0.0:
          factorlow = (kTlist[ite[1]]-kT)/Tdiff
          factorhigh = (kT-kTlist[ite[0]])/Tdiff
          equilib = factorlow * d['EIGEN'].data['FEQB'][ite[0]]+\
                factorhigh * d['EIGEN'].data['FEQB'][ite[1]]
        else:
          equilib = d['EIGEN'].data['FEQB'][ite[0]]

    #renormalize
        equilib /= sum(equilib)
        return equilib

  # now do the non-equilibrium data


    # renormalize
    # make Te into a vector
      kT_vec, kT_isvec = util.make_vec(kT)
      tau_vec, tau_isvec = util.make_vec(tau)
      frac_out = numpy.zeros([len(kT_vec),len(tau_vec),Z+1], dtype=float)
      for ikT, kT in enumerate(kT_vec):
        kTindex = numpy.argmin(numpy.abs(kTlist-kT))

        lefteigenvec = numpy.zeros([Z,Z], dtype=float)
        righteigenvec = numpy.zeros([Z,Z], dtype=float)
        if Z==1:
          for i in range(Z):
            for j in range(Z):
              lefteigenvec[i,j] = d['EIGEN'].data['VL'][kTindex]
              righteigenvec[i,j] = d['EIGEN'].data['VR'][kTindex]
        else:
          for i in range(Z):
            for j in range(Z):
              lefteigenvec[i,j] = d['EIGEN'].data['VL'][kTindex][i*Z+j]
              righteigenvec[i,j] = d['EIGEN'].data['VR'][kTindex][i*Z+j]


        work = numpy.array(init_pop_calc[1:] - d['EIGEN'].data['FEQB'][kTindex][1:], dtype=float)

        fspectmp = numpy.matrix(lefteigenvec) * numpy.matrix(work).transpose()

        delt = 1.0

        worktmp = numpy.zeros(Z)

        for itau, ttau in enumerate(tau_vec):
          if Z >1:
            for i in range(Z):
              worktmp[i] = fspectmp[i]*numpy.exp(d['EIGEN'].data['EIG'][kTindex,i]*delt*ttau)

          else:
            worktmp[0] = fspectmp[0]*numpy.exp(d['EIGEN'].data['EIG'][kTindex]*delt*ttau)

          frac = numpy.zeros(Z+1)
          for i in range(Z):
            for j in range(Z):
              frac[i+1] += worktmp[j]*righteigenvec[j][i]
            frac[i+1] += d['EIGEN'].data['FEQB'][kTindex][i+1]

          if debug:
            frac_out[ikT, itau,:] = frac
          frac[frac<0.0] = 0.0

          if sum(frac)> 1.0:
            frac = frac/sum(frac)
          frac[0] = 1-sum(frac[1:])
          if not(debug):
            frac_out[ikT,itau,:]=frac


      if not tau_isvec:
        frac_out = frac_out.sum(1)
      if not kT_isvec:
        frac_out = frac_out.sum(0)


      return frac_out



    def ion_fraction_delta(self, Z, z1, Te, factor, errortype):

        """
        Ion fraction including the effects of atomic data uncertainty
        Parameters
        ----------
        Z: int
          The nuclear charge of the element
        z1 : int
          ion charge
        Te: float
          Temperature in K
        factor: float
          This determines which atomic data (e.g. collision strength) is adjusted by what percentage. factor = 1.1 means 10% increase in the atomic data 
        errortype: 
          type of error, example, LA, EC,PC etc.

        Returns: ionfraction for ion charge z1
        """

        
        
        for z1 in range(1,Z+1):
            ionfrac = element.calc_elem_ionbal_delta(self, Z, Te, factor, errortype)
 
        return ionfrac[z1-2]



    def level_population_delta(self, Z, z1, Te, N_e, factor, errortype, settings, upper_level, lower_level):


        """
        Calculate level populations for all the levels of Z, z1

        Parameters
        ----------
        Z: int
          The nuclear charge of the element
        z1 : int
          ion charge
        Te : float
          Temperature in K
        N_e : float
          Electron Density (cm^-3)
        Te: float
          Temperature in K
        factor: float
          This determines which atomic data (e.g. collision strength) is adjusted by what percentage. factor = 1.1 means 10% increase in the atomic data 
        errortype: 
          type of error, example, LA, EC,PC etc.
        settings : dict
        """


        
        datacache = {}
        
        up, lo, rates = variableapec.gather_rates_delta(self, Z, z1, Te, N_e,factor, errortype, upper_level, lower_level)
        lev_pop = ion.solve_level_pop (self, up,lo,rates,settings)
        #lev_pop *= ion.ion_fraction(self, Z, z1)
        return lev_pop




    def do_lines_delta(self, Z, z1, lev_pop, N_e, factors, errortypes, upper_levels, lower_levels, datacache=False, settings=False, z1_drv_in=-1):
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
      #print("starting do_lines at %s"%(time.asctime()))
      tstart=time.time()
      
      Atomic_data =  ion.datacache(self, Z, z1)
      #ladat = atomdb.get_data(Z,z1,'LA', datacache=datacache, settings=settings)
      ladat = Atomic_data['data'][Z][z1]['LA']
      lvdat = Atomic_data['data'][Z][z1]['LV']

      
      ebins = setup.make_vector_nbins(self, settings['LinearGrid'], \
                              settings['GridMinimum'], \
                              settings['GridMaximum'], \
                              settings['NumGrid'])
      twoph = numpy.zeros(settings['NumGrid'], dtype=float)
      linelist = numpy.zeros(len(ladat[1].data), \
             dtype= setup.generate_datatypes(self,'linetype'))
      goodlines = numpy.ones(len(linelist), dtype=bool)

      if z1_drv_in < 0:
        z1_drv = z1
      else:
        z1_drv = z1_drv_in

  # now do this
      linelist['epsilon'] = ladat[1].data.field('einstein_a') * \
                        lev_pop[ladat[1].data.field('upper_lev')-1]/N_e
      linelist['lambda'] = ladat[1].data.field('wavelen')


      igood = numpy.isfinite(ladat[1].data['wave_obs'])

      if sum(igood) > 0 :
        igood = numpy.where(igood==True)[0]

        igood = igood[numpy.where(ladat[1].data['wave_obs'][igood]>0)[0]]

      if len(igood)  >0:
        linelist['lambda'][igood]=ladat[1].data['wave_obs'][igood]
#  for iline, line in enumerate(ladat[1].data):
#    if numpy.isfinite(line['wave_obs']):
#      if line['wave_obs'] > 0:
#        linelist[iline]['lambda'] = line['wave_obs']
#      else:
#        linelist[iline]['lambda'] = line['wavelen']

#    else:
#      linelist[iline]['lambda'] = line['wavelen']

      linelist['lambda_err'] = numpy.nan
      linelist['epsilon_err'] = numpy.nan
      linelist['element'] = Z
      linelist['elem_drv'] = Z
      linelist['ion'] = z1
      linelist['ion_drv']= z1_drv
      linelist['upperlev'] = ladat[1].data['upper_lev']
      

      linelist['lowerlev'] = ladat[1].data['lower_lev']


     





      for err in range(len(errortypes)):

        errortype = errortypes[err]
        #print(errortype)

        for num in range(len(upper_levels)):
          
          upper_level = int(upper_levels[num])
          lower_level = int(lower_levels[num])
          
          
      

          for i in range(len(linelist['upperlev'])):
            #print(i)

            if errortype =='LA':
              factor_LA = factors[err][:]
              
              #print(linelist['epsilon'][i])
              if linelist['upperlev'][i] == upper_level:
                if linelist['lowerlev'][i] == lower_level:
                  linelist['epsilon'][i]= factor_LA[num]*linelist['epsilon'][i]
          

                  
                  


      

  # I have a linelist. Yay.
      t1=time.time()
      #print("finished making linelist at %s: took %f seconds"%(time.asctime(), t1-tstart))

  # now check for 2 photon transitions
      #print("starting check for 2 photon transitions at %s"%(time.asctime()))
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
      t2=time.time()
      #print("finished checking two photon transitions at %s: took %f seconds"%(time.asctime(), t2-t1))

      linelist = linelist[goodlines]

      tfinish=time.time()
      #print("finished do_lines at %s, took %f seconds"%(time.asctime(), tfinish-tstart))

      return linelist, twoph
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------




    def calc_recomb_popn_delta(self, levpop, Z, z1, z1_drv, Te, N_e, factor, errortype, upper_level, lower_level, drlevrates, rrlevrates,\
                     settings=False, datacache=False, dronly=False,\
                     rronly=False):
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
      Te: electron temperature (K)
      N_e: electron density (cm^-3)
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

      Atomic_data =  ion.datacache(self, Z, z1)
      lvdat =  Atomic_data['data'][Z][z1]['LV']
      
      if not lvdat:
        nlev = 1
        levpop = numpy.zeros(1, dtype=float)
        return levpop
      nlev = len(lvdat[1].data)
      Tarr, dummy = util.make_vec(Te)

      if nlev > const.NLEV_NOSPARSE:
        print("using sparse solver for recomb")

    
        aidat =  Atomic_data['data'][Z][z1]['AI']
        
        if aidat:
          ailev = numpy.array(util.unique(aidat[1].data['level_init']))-1
          nailev = len(ailev)
          isbound = numpy.ones(nlev, dtype=bool)
          isbound[ailev]=False
        else:
          nailev = 0
          isbound = numpy.ones(nlev, dtype=bool)

    

        recombrate = numpy.zeros(nlev, dtype=float)

        #get the recomb data
        irdat =  Atomic_data['data'][Z][z1]['IR']
        

        for iir, ir in enumerate(irdat[1].data):
        #check we have the right data types
          if ir['TR_TYPE'] in ['RR','DR','XR']:
            

            if errortype=='REC':
              recrate = factor*atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)*N_e
            else:
              recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)*N_e

            if not (numpy.isfinite(recrate)):
              print("iir=%i, recrate is not finite!"%(iir))
            else:
              recombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]

        maxlev = numpy.where(recombrate > 0)[0]
        if len(maxlev) == 0:  # so recombrate has all the influxes
          return numpy.zeros(nlev)
        maxlev=maxlev[-1]
        matrixB = recombrate
        #print "Sum recomb rates in level>1:", sum(matrixB[1:])



        matrixA = {}
        matrixA['init'], matrixA['final'], matrixA['rate']=\
          variableapec.gather_rates_delta(self, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)


        matrixB *= -1
        nlev = len(matrixB)



        i = (matrixA['final']>0) & (matrixA['init']>0)
        matrixA['final'] = matrixA['final'][i]
        matrixA['init'] = matrixA['init'][i]
        matrixA['rate'] = matrixA['rate'][i]

    #remove all matrix A from and to high levels
        i = (matrixA['final']<maxlev+1) & (matrixA['init']<maxlev+1)
        matrixA['final'] = matrixA['final'][i]
        matrixA['init'] = matrixA['init'][i]
        matrixA['rate'] = matrixA['rate'][i]

    #subtract 1 from the levels
        matrixA['init']-= 1
        matrixA['final']-= 1

        A  = numpy.zeros([maxlev,maxlev])
        for i in range(len(matrixA['final'])):
          A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]


        levpop_this = numpy.zeros(len(matrixB))

        if sum(matrixB[1:] < 0):
          levpop_this[1:maxlev+1] = numpy.linalg.solve(A, matrixB[1:maxlev+1])



      else:

        #print("using regular solver for recomb")

        rrrecombrate = numpy.zeros(nlev, dtype=float)
        drrecombrate = numpy.zeros(nlev, dtype=float)
        irdat =  Atomic_data['data'][Z][z1]['IR']
        

        havedrrate=False
        haverrrate=False
        for iir, ir in enumerate(irdat[1].data):
      # check we have the right data types
          if ir['TR_TYPE'] in ['RR','XR']:

            recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)

            rrrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*N_e

            if ((ir['TR_TYPE'] in ['RR','XR']) & (ir['level_final']>1)):
              haverrrate=True
          if ir['TR_TYPE'] in ['DR','XD']:
            recrate = atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdat)
            drrecombrate[ir['level_final']-1] += recrate*levpop[ir['level_init']-1]*N_e
            if ((ir['TR_TYPE'] in ['DR','XD']) & (ir['level_final']>1)):
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

        #print("DR: sum from satellite lines: %e, sum from IR file: %e" %\
              #(sumdrlevrates, sum(drrecombrate)))
        #print("RR: sum from PI xsections: %e, sum from IR file: %e" %\
              #(sumrrlevrates, sum(rrrecombrate)))

        matrixB = rrrecombrate+drrecombrate+tmpdrlevrates+tmprrlevrates
        if dronly:
          matrixB = drrecombrate+tmpdrlevrates
        if rronly:
          matrixB = rrrecombrate+tmprrlevrates


        matrixA = numpy.zeros([nlev,nlev],dtype=float)

        ladat =  Atomic_data['data'][Z][z1]['LA']
        

        matrixA_in = {}
        matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
          variableapec.gather_rates_delta(self, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

        #datacache={}


        for i in range(len(matrixA_in['init'])):
          matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]

    # solve unless matrixB ==0
        if sum(matrixB[1:])>0:
          matrixB = -1*matrixB
          levpop_this = setup.calc_cascade_population(self, matrixA, matrixB)
        else:
          levpop_this = numpy.zeros(nlev)


      #print("level population for recombination into Z=%i, z1=%i, z1_drv=%i, levpop=%i"%\
            #(Z, z1, z1_drv,levpop))
      #for i in range(len(levpop_this)):
        #print(i, levpop_this[i])
      return levpop_this






    def calc_ioniz_popn_delta(self, levpop, Z, z1, z1_drv, Te, N_e, factor, errortype, upper_level, lower_level, settings=False, \
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
        factor: float  
        errortype:str
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

        datacache = {}
        #print("Starting calc_ioniz_popn at %s"%(time.asctime()))
        Atomic_data =  ion.datacache(self, Z, z1)
        lvdat =  Atomic_data['data'][Z][z1]['LV']
        
        #lvdat = atomdb.get_data(Z,z1,'LV', settings=settings, datacache=datacache)

  # if we have no lv data, ignore.
        if not util.keyword_check(lvdat):
            nlev = 1
            return numpy.array([0.0])
        nlev = len(lvdat[1].data)

  # get populating rate from previous ion
        Atomic_data_prev =  ion.datacache(self, Z, z1-1)
        aidat = Atomic_data_prev['data'][Z][z1-1]['AI']
        ionizrateai=numpy.zeros(nlev, dtype=float)
        ionizrateir=numpy.zeros(nlev, dtype=float)

        #print("Starting calc_ioniz_popn aidat loop at %s"%(time.asctime()))
        if aidat:
            tmp_pop = levpop[aidat[1].data['level_init']-1]
            for iai in range(len(aidat[1].data)-1):

              if  errortype=='ION':
                ionizrateai[aidat[1].data['level_final'][iai]-1] += \
                    tmp_pop[iai]*aidat[1].data['auto_rate'][iai]*factor

              else:
                ionizrateai[aidat[1].data['level_final'][iai]-1] += \
                    tmp_pop[iai]*aidat[1].data['auto_rate'][iai]
              
    #aidat.close()
        #print("Finished calc_ioniz_popn aidat loop at %s"%(time.asctime()))


        #print("Starting calc_ioniz_popn xidat loop at %s"%(time.asctime()))
        if do_xi:

            irdat = Atomic_data_prev['data'][Z][z1-1]['IR']
            #irdat = atomdb.get_data(Z, z1-1, 'IR', settings=settings, datacache=datacache)
            ionpot = float(irdat[1].header['ionpot'])
            if z1 >1:
                lvdatm1 = Atomic_data_prev['data'][Z][z1-1]['LV']
  # go through each excitation, have fun

            for iir, ir in enumerate(irdat[1].data):
                if ir['TR_TYPE'] in ['XI']:
                    Tarr =  numpy.array([Te])
                    ionrate=atomdb.get_maxwell_rate(Tarr, irdat, iir, lvdatm1, \
                                     lvdatap1=lvdat, ionpot=ionpot)
                    
                    ionizrateir[ir['level_final']-1] += levpop[ir['level_init']-1]*\
                                       ionrate




        ionizrate=ionizrateir+ionizrateai
        matrixB = ionizrate

  # save some time if there is nothing to ionize.

        if sum(matrixB[1:]) ==0:
            levpop_this = numpy.zeros(len(matrixB))
            return levpop_this

        maxlev = numpy.where(matrixB > 1e-40)[0]
        if len(maxlev)==0:
            popn = numpy.zeros(len(matrixB))
            return popn
        maxlev=maxlev[-1]
  
        matrixA_in={}
        matrixA_in['init'], matrixA_in['final'], matrixA_in['rate'] = \
            variableapec.gather_rates_delta(self, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings,\
                        do_la=True, do_ai=True, do_ec=False, do_pc=False,\
                        do_ir=False)

        i = (matrixA_in['init']<=maxlev) & (matrixA_in['final']<=maxlev)
        matrixA_in['init']=matrixA_in['init'][i]
        matrixA_in['final']=matrixA_in['final'][i]
        matrixA_in['rate']=matrixA_in['rate'][i]

  # fix the rates
        for i in range(len(matrixA_in['init'])):
            if matrixA_in['init'][i]==matrixA_in['final'][i]:
                if matrixA_in['rate'][i] >=0.0:
                    matrixA_in['rate'][i] -=1e10
                    #print("CTieing level %i to ground with rate 1e10"%(i))

        if (maxlev <= const.NLEV_NOSPARSE):
    # convert to a regular solver
            #print("regular solver")
            matrixA = numpy.zeros([maxlev+1,maxlev+1], dtype=float)

            for i in range(len(matrixA_in['init'])):
                matrixA[matrixA_in['final'][i], matrixA_in['init'][i]] += matrixA_in['rate'][i]


    # bug-u-fix
            for i in range(1, maxlev):
                if matrixA[i,i] >= 0:
                    matrixA[i,i]=-1e10
                    #print("FIXING matrixA[%i,%i] = -1.0"%(i,i))

            popn = numpy.zeros(nlev)

            matrixB*=-1

            try:
                popn[1:maxlev] = numpy.linalg.solve(matrixA[1:maxlev,1:maxlev], matrixB[1:maxlev])
            except numpy.linalg.linalg.LinAlgError:
                "EEK ERROR!"
                raise

    
        else:
    # add into sparse solver
            #print("Using sparse solver")
            matrixA={}
            matrixB *= -1
            nlev = len(matrixB)

            if sum(matrixB)>=0:
                return numpy.zeros(len(matrixB))

    # remove ground level
            i = (matrixA_in['init']>0) & (matrixA_in['final']>0)

            matrixA['init'] = matrixA_in['init'][i]
            matrixA['final'] = matrixA_in['final'][i]
            matrixA['rate'] = matrixA_in['rate'][i]

            i = (matrixA['init']<=maxlev+1) & (matrixA['final']<=maxlev+1)

            matrixA['init'] = matrixA['init'][i]
            matrixA['final'] = matrixA['final'][i]
            matrixA['rate'] = matrixA['rate'][i]


    # subtract 1 from the levels
            matrixA['init']-=1
            matrixA['final']-=1

#    A = sparse.coo_matrix((matrixA['rate'],\
#                           (matrixA['final'],matrixA['init'])), \
#                           shape=(maxlev,maxlev)).tocsr()
            A = numpy.zeros([maxlev, maxlev])
            for i in range(len(matrixA['final'])):
                A[matrixA['final'][i], matrixA['init'][i]]+=matrixA['rate'][i]

            popn = numpy.zeros(len(matrixB))
#    popn[1:maxlev+1] = spsolve(A, matrixB[1:maxlev+1])
            popn[1:maxlev+1] = numpy.linalg.solve(A, matrixB[1:maxlev+1])

            popn_bak = popn*1.0

            popn[popn<0] = 0.0
    
  
        
        return popn







    def run_apec_ion_delta(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund):
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
        
        settings['WriteIonFname'] = "Z_%i_z1_%i_T_%i_N_%i"%(Z,z1,Te,N_e)
        settings['filemap'] = settings['FileMap']
        settings['atomdbroot'] = os.path.expandvars('$ATOMDB')



  # get the output energy bins for the continuum
        ebins = setup.make_vector_nbins(self, settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])

  ## FIXME CUTOFF FOR MIN IONPOP
        linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
        pseudo = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum = {}
        continuum['brems'] = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum['twophot'] = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum['rrc'] = numpy.zeros(settings['NumGrid'], dtype=float)


  ## FIXME CUTOFF FOR MIN 
        #if ionfrac[z1_drv-1] < const.MIN_IONPOP:
        #    return  linelist, continuum, pseudo
        

  # set up the datacache

        datacache = {}

  # Find the number of levels
        Atomic_data =  ion.datacache(self, Z, z1)
        #print(Atomic_data)
        lvdat = Atomic_data['data'][Z][z1]['LV']
        
        
        linelist_exc = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        linelist_dr = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        linelist_rec = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        if lvdat!= False:
    # find the number of levels
            nlev = len(lvdat[1].data)

    # check if we need to do any of the line-related calculations

            if (settings['EmissionLines'] or settings['TwoPhoton']):


    # gather all the level to level rates
                
                
                up, lo, rates = variableapec.gather_rates_delta(self, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
                lev_pop = variableapec.level_population_delta(self, Z, z1, Te, N_e, factor, errortype, settings, upper_level, lower_level)*Abund[Z]*ionfrac[z1-1]
                
                
                

    # just in case, add zeros to lengthen the lev_pop appropriately
                if len(lev_pop) < nlev:
                    lev_pop = numpy.append(lev_pop, numpy.zeros(nlev-len(lev_pop), dtype=float))

    # fix any sub-zero level populations
                lev_pop[lev_pop<0] = 0.0
                
    
      
                linelist_exc,  continuum['twophot'] = variableapec.do_lines_delta(self, Z, z1, lev_pop, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings, z1_drv_in=z1_drv)
                


            else:
                # skipping the exact level calculation, fill it with zeros, ground state with 1.
                lev_pop = numpy.zeros(nlev, dtype=float)
                lev_pop[0] = 1.0*Abund[Z]*ionfrac[z1-1]

        else:
            lev_pop=numpy.ones(1, dtype=float)*Abund[Z]*ionfrac[z1-1]


        if settings['Bremsstrahlung'] ==True:
            brems = ion.do_brems(self, Z, z1, Te, N_e, settings['BremsType'], ebins)
    # scale for  ion and element abundance.
            continuum['brems']=brems*Abund[Z]*ionfrac[z1-1]

    
        else:
            continuum['brems']=numpy.zeros(len(ebins)-1, dtype=float)
    # now look at the neighbouring ions
        if z1_drv>1:
            z1=z1_drv-1
            if settings['DRSatellite']:
                linelist_dr, drlevrates = ion.calc_satellite(self, Z, z1, Te, datacache=datacache, settings=settings)
                linelist_dr['epsilon']*=Abund[Z]*ionfrac[z1_drv-1]
                drlevrates *=Abund[Z]*ionfrac[z1_drv-1]

            else:
                linelist_dr = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
                drlevrates = 0.0

    # Radiative Recombination
            if settings['RRC']:
      
                rrc, rrlevrates = atomdb.calc_rad_rec_cont(Z, z1, z1_drv, Te, ebins, settings=settings, datacache=datacache)
                continuum['rrc'] = rrc*Abund[Z]*ionfrac[z1_drv-1]
                rrlevrates*=Abund[Z]*ionfrac[z1_drv-1]

            else:
                continuum['rrc'] = numpy.zeros(len(ebins)-1, dtype=float)
                rrlevrates=0.0


    # if there is recombination to process:
            tmpdrlevrates,xxx = util.make_vec(drlevrates)
            tmprrlevrates,xxx = util.make_vec(rrlevrates)

            
            if sum(tmpdrlevrates) + sum(tmprrlevrates)>0:
      

                levpop_recomb= variableapec.calc_recomb_popn_delta(self, lev_pop, Z, z1,\
                                    z1_drv, Te, N_e, factor, errortype, upper_level, lower_level, drlevrates,\
                                     rrlevrates,\
                                     datacache=datacache, settings=settings)

                


                linelist_rec, tmptwophot =variableapec.do_lines_delta(self, Z, z1, levpop_recomb , N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings, z1_drv_in=z1_drv)  
        

                continuum['twophot']+= tmptwophot
            
        # now do the ionizing cases
        linelist_ion = numpy.zeros(0,dtype= setup.generate_datatypes(self,'linetype'))
        if z1_drv < Z:
    #datacache={}
            z1=z1_drv+1
            lev_pop_parent = lev_pop*1.0

            while (sum(lev_pop_parent[1:]) > 1e-40) &\
                (z1 <= Z):

                if z1== z1_drv+1:
                    do_xi = True
                else:
                    do_xi = False

                lev_pop = variableapec.calc_ioniz_popn_delta(self, lev_pop_parent, Z, z1, z1_drv, Te, N_e, factor, errortype, upper_level, lower_level, \
                                settings=settings, datacache=datacache, \
                                do_xi=do_xi)

                lev_pop[lev_pop<const.MIN_LEVPOP] = 0.0
                if sum(lev_pop[1:]) > 0:
                    linelist_ion_tmp, tmptwophot = \
                        variableapec.do_lines_delta(self, Z, z1, lev_pop, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings,  z1_drv_in=z1_drv)

                    linelist_ion = numpy.append(linelist_ion, linelist_ion_tmp)
                    continuum['twophot']+=tmptwophot

                lev_pop_parent = lev_pop
                z1+=1

        linelist = numpy.append(linelist_exc, numpy.append(linelist_dr, numpy.append(linelist_ion, linelist_rec)))
        
        MinEpsilon = settings['MinEpsilon']
        if settings['Ionization']=='CIE':
            MinEpsilon*=0.001

        pseudocont = numpy.zeros(len(ebins)-1, dtype=float)

        if len(linelist) > 0:
            weaklines = linelist[(linelist['epsilon']< MinEpsilon) &\
                         (linelist['lambda']>const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                         (linelist['lambda']<const.HC_IN_KEV_A /settings['GridMinimum'])]

            for line in weaklines:
                e = const.HC_IN_KEV_A /line['lambda']
                ibin = numpy.where(ebins>e)[0][0] - 1
                pseudocont[ibin]+=line['epsilon']

            linelist = linelist[linelist['epsilon'] > MinEpsilon]
        
      
        if settings['WriteIon']==True:
        
            
            ret = {}
            ret['lines'] = linelist
            ret['continuum'] = continuum
            ret['pseudocont'] = pseudocont
            ret['ionfrac'] = ionfrac
            ret['te'] = Te
            ret['dens'] = N_e
            ret['settings'] = settings
            ret['abund'] = Abund[Z]
            fname = settings['WriteIonFname']+'.pkl'
        
            pickle.dump(ret, open(fname, 'wb'))

        
        return linelist





    





    def run_apec_ion_delta_exc(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund):

        
        z1_drv=z1*1
        
        settings['WriteIonFname'] = "Z_%i_z1_%i_T_%i_N_%i"%(Z,z1,Te,N_e)
        settings['filemap'] = settings['FileMap']
        settings['atomdbroot'] = os.path.expandvars('$ATOMDB')



  # get the output energy bins for the continuum
        ebins = setup.make_vector_nbins(self, settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])

  ## FIXME CUTOFF FOR MIN IONPOP
        linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
        pseudo = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum = {}
        continuum['brems'] = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum['twophot'] = numpy.zeros(settings['NumGrid'], dtype=float)
        continuum['rrc'] = numpy.zeros(settings['NumGrid'], dtype=float)


  ## FIXME CUTOFF FOR MIN 
        #if ionfrac[z1_drv-1] < const.MIN_IONPOP:
        #    return  linelist, continuum, pseudo
        

  # set up the datacache

        datacache = {}

  # Find the number of levels
        Atomic_data =  ion.datacache(self, Z, z1)
        #print(Atomic_data)
        lvdat = Atomic_data['data'][Z][z1]['LV']
        
        
        linelist_exc = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        
        if lvdat!= False:
    # find the number of levels
            nlev = len(lvdat[1].data)

    # check if we need to do any of the line-related calculations

            if (settings['EmissionLines'] or settings['TwoPhoton']):


    # gather all the level to level rates
                
                
                up, lo, rates = variableapec.gather_rates_delta(self, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
                lev_pop = variableapec.level_population_delta(self, Z, z1, Te, N_e, factor, errortype, settings, upper_level, lower_level)*Abund[Z]*ionfrac[z1-1]
                
                
                

    # just in case, add zeros to lengthen the lev_pop appropriately
                if len(lev_pop) < nlev:
                    lev_pop = numpy.append(lev_pop, numpy.zeros(nlev-len(lev_pop), dtype=float))

    # fix any sub-zero level populations
                lev_pop[lev_pop<0] = 0.0
                
    
      
                linelist_exc,  continuum['twophot'] = variableapec.do_lines_delta(self, Z, z1, lev_pop, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings, z1_drv_in=z1_drv)
                


            else:
                # skipping the exact level calculation, fill it with zeros, ground state with 1.
                lev_pop = numpy.zeros(nlev, dtype=float)
                lev_pop[0] = 1.0*Abund[Z]*ionfrac[z1-1]

        else:
            lev_pop=numpy.ones(1, dtype=float)*Abund[Z]*ionfrac[z1-1]



        return linelist_exc




    def run_apec_ion_delta_single_exc(self, settings, Z, z1, Te, N_e, uls, lls, factor, errortype, upper_level, lower_level, ionfrac, Abund):
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
        

        
        Atomic_data =  ion.datacache(self, Z, z1)
        ladat = Atomic_data['data'][Z][z1]['LA']
        ULV =  ladat[1].data.field('upper_lev')
        LLV = ladat[1].data.field('lower_lev')
        linelist=variableapec.run_apec_ion_delta_exc(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund)
        #Wavelen = ladat[1].data.field('wavelen')
        #print(Wavelen)
        #print(linelist['lambda'])

        
       
          
        
        Emissivity = []
        for i in range(len(linelist)):
          #if Wavelen[i] == linelist['lambda'][i]:
            if linelist[i][8] in uls and linelist[i][9] in lls:
              Emissivity.append(linelist['epsilon'][i])
              lines = linelist[i]
        return Emissivity

              


        
        #print(LLV,len(LLV))
        #print(ULV)
        

        

  ## FIXME CUTOFF FOR MIN 
        #if ionfrac[z1_drv-1] < const.MIN_IONPOP:
        #    return  linelist, continuum, pseudo
        

        












    def run_apec_ion_delta_rec(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund):

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
        
        

    # gather all the level to level rates
                
                
        datacache={}
        #Atomic_data =  ion.datacache(self, Z, z1)
        #print(Atomic_data)
        #lvdat = Atomic_data['data'][Z][z1]['LV']

        ebins = setup.make_vector_nbins(self, settings['LinearGrid'], settings['GridMinimum'], settings['GridMaximum'], settings['NumGrid'])


        #linelist_rec = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        
  
            
        up, lo, rates = variableapec.gather_rates_delta(self, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
        lev_pop = variableapec.level_population_delta(self, Z, z1, Te, N_e, factor, errortype, settings, upper_level, lower_level)*Abund[Z]*ionfrac[z1-1]


            

        linelist_dr, drlevrates = ion.calc_satellite(self, Z, z1, Te)
        rrc, rrlevrates = atomdb.calc_rad_rec_cont(Z, z1, z1+1, Te,ebins)
        levpop_recomb=variableapec.calc_recomb_popn_delta(self, lev_pop, Z, z1,\
                                    z1+1, Te, N_e, factor, errortype, upper_level, lower_level, drlevrates,\
                                     rrlevrates,\
                                     datacache=datacache, settings=settings)

        linelist_rec, tmptwophot = variableapec.do_lines_delta(self, Z, z1, levpop_recomb, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings, z1_drv_in=z1+1)

                

        
        return linelist_rec







    def run_apec_ion_delta_single_rec(self, settings, Z, z1, Te, N_e, uls, lls, factor, errortype, upper_level, lower_level, ionfrac, Abund):

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
        Returns
        -------
        Emissivity: Returns emissivity when driving ion has the charge z1+1
    
        """

        # get the data.
        
        

    # gather all the level to level rates
                
                
        datacache={}
        #Atomic_data =  ion.datacache(self, Z, z1)
        #print(Atomic_data)
        #lvdat = Atomic_data['data'][Z][z1]['LV']

        ebins = setup.make_vector_nbins(self, settings['LinearGrid'], settings['GridMinimum'], settings['GridMaximum'], settings['NumGrid'])


        #linelist_rec = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        
  
            
        up, lo, rates = variableapec.gather_rates_delta(self, Z, z1+1, Te, N_e, factor, errortype, upper_level, lower_level, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
        lev_pop = variableapec.level_population_delta(self, Z, z1+1, Te, N_e, factor, errortype, settings, upper_level, lower_level)*Abund[Z]*ionfrac[z1-1]


            

        linelist_dr, drlevrates = ion.calc_satellite(self, Z, z1, Te)
        rrc, rrlevrates = atomdb.calc_rad_rec_cont(Z, z1, z1+1, Te,ebins)
        levpop_recomb=variableapec.calc_recomb_popn_delta(self, lev_pop, Z, z1,\
                                    z1+1, Te, N_e, factor, errortype, upper_level, lower_level, drlevrates,\
                                     rrlevrates,\
                                     datacache=datacache, settings=settings)

        linelist_rec, tmptwophot = variableapec.do_lines_delta(self, Z, z1, levpop_recomb, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings, z1_drv_in=z1+1)
        
        Emissivity_rec=[]

        for i in range(len(linelist_rec)):
          #if Wavelen[i] == linelist['lambda'][i]:
            if linelist_rec[i][8] in uls and linelist_rec[i][9] in lls:
              Emissivity_rec.append(linelist_rec['epsilon'][i])
              lines = linelist_rec[i]
              


                

        
        return Emissivity_rec










    def run_apec_ion_delta_single_ion(self, settings, Z, z1, Te, N_e, uls, lls, factor, errortype, upper_level, lower_level, ionfrac, Abund):

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
        
        

    # gather all the level to level rates
                
                
        datacache={}
        z1_drv=z1-1

        #linelist_rec = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        
  
            
        up, lo, rates = variableapec.gather_rates_delta(self, Z, z1-1, Te, N_e, factor, errortype, upper_level, lower_level, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
        lev_pop = variableapec.level_population_delta(self, Z, z1-1, Te, N_e, factor, errortype, settings, upper_level, lower_level)*Abund[Z]*ionfrac[z1-1]


        linelist_ion = numpy.zeros(0,dtype= setup.generate_datatypes(self,'linetype'))
        if z1_drv < Z:
          z1=z1_drv+1
          lev_pop_parent = lev_pop*1.0
          if z1== z1_drv+1:
            do_xi = True
          else:
            do_xi = False



          

          lev_pop = variableapec.calc_ioniz_popn_delta(self, lev_pop_parent, Z, z1, z1_drv, Te, N_e, factor, errortype, upper_level, lower_level, \
                                settings=settings, datacache=datacache, \
                                do_xi=True)



          lev_pop[lev_pop<const.MIN_LEVPOP] = 0.0
          if sum(lev_pop[1:]) > 0:
            linelist_ion_tmp, tmptwophot = \
                        variableapec.do_lines_delta(self, Z, z1, lev_pop, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings,  z1_drv_in=z1_drv)

            linelist_ion = numpy.append(linelist_ion, linelist_ion_tmp)


          Emissivity_ion=[]

          for i in range(len(linelist_ion)):
          #if Wavelen[i] == linelist['lambda'][i]:
            if linelist_ion[i][8] in uls and linelist_ion[i][9] in lls:
              Emissivity_ion.append(linelist_ion['epsilon'][i])
              #lines = linelist_ion[i]
              


                

        
          return Emissivity_ion


        
          
  



    def run_apec_ion_delta_ion(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund):

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
        
        

    # gather all the level to level rates
                
                
        datacache={}
        z1_drv=z1-1

        #linelist_rec = numpy.zeros(0, dtype= setup.generate_datatypes(self,'linetype'))
        
  
            
        up, lo, rates = variableapec.gather_rates_delta(self, Z, z1-1, Te, N_e, factor, errortype, upper_level, lower_level, do_la=True, \
                                                    do_ec=True, do_ir=True, do_pc=True, do_ai=True, datacache=datacache)
      
                
                
        lev_pop = variableapec.level_population_delta(self, Z, z1-1, Te, N_e, factor, errortype, settings, upper_level, lower_level)*Abund[Z]*ionfrac[z1-1]


        linelist_ion = numpy.zeros(0,dtype= setup.generate_datatypes(self,'linetype'))
        if z1_drv < Z:
          z1=z1_drv+1
          lev_pop_parent = lev_pop*1.0
          if z1== z1_drv+1:
            do_xi = True
          else:
            do_xi = False



          

          lev_pop = variableapec.calc_ioniz_popn_delta(self, lev_pop_parent, Z, z1, z1_drv, Te, N_e, factor, errortype, upper_level, lower_level, \
                                settings=settings, datacache=datacache, \
                                do_xi=True)



          lev_pop[lev_pop<const.MIN_LEVPOP] = 0.0
          if sum(lev_pop[1:]) > 0:
            linelist_ion_tmp, tmptwophot = \
                        variableapec.do_lines_delta(self, Z, z1, lev_pop, N_e, factor, errortype, upper_level, lower_level, datacache=datacache, settings=settings,  z1_drv_in=z1_drv)

            linelist_ion = numpy.append(linelist_ion, linelist_ion_tmp)


        
          return linelist_ion
  










    




    




    

    









    def run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, ul, ll, factor, errortype, upper_level, lower_level, ionfrac, Abund):
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
        

        
        Atomic_data =  ion.datacache(self, Z, z1)
        ladat = Atomic_data['data'][Z][z1]['LA']
        ULV =  ladat[1].data.field('upper_lev')
        LLV = ladat[1].data.field('lower_lev')
        linelist=variableapec.run_apec_ion_delta(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund)
        #Wavelen = ladat[1].data.field('wavelen')
        #print(Wavelen)
        #print(linelist['lambda'])

        
       
          
        

        for i in range(len(linelist)):
          #if Wavelen[i] == linelist['lambda'][i]:
            if linelist[i][8] == ul and linelist[i][9] == ll:
              Emissivity = linelist['epsilon'][i]
              lines = linelist[i]
              break


        
        #print(LLV,len(LLV))
        #print(ULV)
        

        

  ## FIXME CUTOFF FOR MIN 
        #if ionfrac[z1_drv-1] < const.MIN_IONPOP:
        #    return  linelist, continuum, pseudo
        

        return Emissivity








    
    def z_by_w(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund):


 
        z = variableapec.run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, 2, 1, factor, errortype, upper_level, lower_level, ionfrac, Abund) 
        

        w = variableapec.run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, 7, 1, factor, errortype, upper_level, lower_level, ionfrac, Abund)
        

        ratio = z/w

        return ratio





    def G_ratio(self, settings, Z, z1, Te, N_e, factor, errortype, upper_level, lower_level, ionfrac, Abund):


 
        z = variableapec.run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, 2, 1, factor, errortype, upper_level, lower_level, ionfrac, Abund) 
        x= variableapec.run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, 6, 1, factor, errortype, upper_level, lower_level, ionfrac, Abund) 
        y= variableapec.run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, 5, 1, factor, errortype, upper_level, lower_level, ionfrac, Abund) 
        w = variableapec.run_apec_ion_delta_single(self, settings, Z, z1, Te, N_e, 7, 1, factor, errortype, upper_level, lower_level, ionfrac, Abund)
        

        G = (x+y+z)/w

        return G




    





    












class element():


    '''
    Z=int(input("Z="))
    Te=float(input("Te="))
    N_e=float(input("N_e="))
  #upper_level=int(input("upper_level="))
  #lower_level=int(input("lower_level="))
    factor=float(input("factor="))
    errortype=str(input("errortype="))
    ind=int(input("index (if Te and N_e are known put ind=0)="))
    '''


    
    
    

    def __init__(self, Z, Te, N_e, ind, factor, errortype, Abund):
        '''
        Z=self.Z
        Te=self.Te
        N_e=self.N_e
        factor=self.factor
        errortype=self.errortype
        ind=self.ind
        '''
        #ionfrac= calc_elem_ionbal_delta(self, Z, Te, factor, errortype)
        
        Abund= atomdb.get_abundance(abundfile=False, \
                abundset='AG89', element=[-1], datacache=False, settings=False, show=False)

        fname = os.path.expandvars('$ATOMDB/apec.par')
        settings = parse_par_file(fname)



        



        






    def datacache(self, Z, datacache=False, settings=False):


        """
        Creating an empty dictionary datacache
    
        Reading atomic data for:
        ‘IR’ - ionization and recombination
        ‘LV’ - energy levels
        ‘LA’ - radiative transition data (lambda and A-values)
        ‘EC’ - electron collision data
        ‘PC’ - proton collision data
        ‘DR’ - dielectronic recombination satellite line data
        ‘PI’ - XSTAR photoionization data
        ‘AI’ - autoionization data

        Parameters
        ----------
        Z: int
          The nuclear charge of the element
        z1 : int
          ion charge
        datacache : dict
          Used for caching the data. See description in atomdb.get_data
        settings : dict
          See description in atomdb.get_data

        Returns: datacache
        Storing in datacache for later use
        """
       
        datacache={}
        
        for z1 in range(1, Z+1):
          Alldat = atomdb.get_data(Z, z1, 'ALL', datacache=datacache, \
                            settings = settings)
        
        return datacache



    def extract_gauntff(self, Z, gamma2, gaunt_U, gaunt_Z, gaunt_Ng, gaunt_g2, gaunt_gf):
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
      for iii in ii:
        ng = gaunt_Ng[i[iii]]

        GauntFFvec[iii] = 10**(numpy.interp(gamma2,gaunt_g2[i[iii]][:ng], \
                                        numpy.log10(gaunt_gf[i[iii]][:ng])))
      ii = numpy.where(gamma2<gaunt_g2[i,0])[0]
      GauntFFvec[ii]=gaunt_gf[i[ii],0]

      ii = numpy.where(gamma2>gaunt_g2[i,gaunt_Ng[i]-1])[0]
      GauntFFvec[ii]=gaunt_gf[i[ii],gaunt_Ng[i[ii]]-1]
      return Uvec, GauntFFvec



#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


    def generate_nei_outputs(self, settings, Z, linelist, contlist, pseudolist, ionfrac_nei):
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
      print("initially we have %i lines for Z =%i"%(len(linelist), Z))
      for z1 in range(1, Z+2):
        ionfrac = ionfrac_nei[z1-1]
        mineps = settings['NEIMinEpsilon'][0]
        pseudotmp = numpy.zeros(len(ebins)-1, dtype=float)

        pseudotmp += pseudolist[z1]

        for i in range(len(settings['NEIMinFrac'])):
          if ionfrac < settings['NEIMinFrac'][i]:
            mineps = settings['NEIMinEpsilon'][i+1]
        print("z1 = %i. Ionfrac = %e. mineps = %e"%(z1,ionfrac, mineps))


        weaklines = linelist[(linelist['element']==Z) &\
                (linelist['ion_drv']==z1) &\
                (linelist['epsilon']<mineps) &\
                (linelist['lambda']>const.HC_IN_KEV_A /settings['GridMaximum']) &\
                (linelist['lambda']<const.HC_IN_KEV_A /settings['GridMinimum'])]
        print("identified %i weak lines"%(len(weaklines)))

        for line in weaklines:
          e = const.HC_IN_KEV_A /line['lambda']
          ibin = numpy.where(ebins>e)[0][0] - 1
          pseudotmp[ibin]+=line['epsilon']

        igood[(linelist['element']==Z) &\
              (linelist['ion_drv']==z1) &\
              (linelist['epsilon']<mineps)] = False
        print("Filtered by mineps %e: from %i to %i lines"%(mineps, len(igood), sum(igood)))
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

      print("Filtering lines on wavelength: keeping %i of %i lines"%\
            (sum(igood), len(igood)))
      linelist = linelist[igood]

      ret={}
      ret['lines']=linelist
      maxnpseudo = 0
      maxncont = 0

      for i in list(cont.keys()):
        if len(cont[i]['E_Cont'])>maxncont:
          maxncont= len(cont[i]['E_Cont'])
        if len(cont[i]['E_Pseudo'])>maxnpseudo:
          maxnpseudo= len(cont[i]['E_Pseudo'])

      ret['cont'] = numpy.zeros(len(list(cont.keys())), dtype=generate_datatypes('continuum', npseudo=maxnpseudo, ncontinuum=maxncont))

      for iz1, z1 in enumerate(cont.keys()):
        ret['cont']['Z'][iz1] = Z
        ret['cont']['rmJ'][iz1] = z1
        ret['cont']['N_Cont'][iz1] = len(cont[z1]['E_Cont'])
        ret['cont']['E_Cont'][iz1][:ret['cont']['N_Cont'][iz1]] = cont[z1]['E_Cont']
        ret['cont']['Continuum'][iz1][:ret['cont']['N_Cont'][iz1]] = cont[z1]['Cont']
        ret['cont']['N_Pseudo'][iz1] = len(cont[z1]['E_Pseudo'])
        ret['cont']['E_Pseudo'][iz1][:ret['cont']['N_Pseudo'][iz1]] = cont[z1]['E_Pseudo']
        ret['cont']['Pseudo'][iz1][:ret['cont']['N_Pseudo'][iz1]] = cont[z1]['Pseudo']

      print("returning ret['lines'] with length %i"%(len(ret['lines'])))
      return ret




    def generate_cie_outputs(self, settings, Z, linelist, contlist, pseudolist):
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

#  pickle.dump(linelist, open('whole_linedebugZ%i.'%(Z), 'wb'))

      for i in range(1,len(linelist)):
        if ((linelist['element'][i]==linelist['element'][i-1]) &\
            (linelist['ion'][i]==linelist['ion'][i-1]) &\
            (linelist['upperlev'][i]==linelist['upperlev'][i-1]) &\
            (linelist['lowerlev'][i]==linelist['lowerlev'][i-1])):
          linelist['epsilon'][i] += linelist['epsilon'][i-1]
          igood[i-1] = False

      linelist= linelist[igood]

      linelist_equ = numpy.zeros(len(linelist), dtype=setup.generate_datatypes(self, 'linelist_cie'))

      for i in linelist_equ.dtype.names:
        linelist_equ[i] = linelist[i]


      ebins =  setup.make_vector_nbins(self,settings['LinearGrid'], \
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

      for z1 in  list(contlist.keys()):
        cont_rrc += contlist[z1]['rrc']
        cont_2ph += contlist[z1]['twophot']
        cont_bre += contlist[z1]['brems']

      for z1 in  list(pseudolist.keys()):
        pseudo += pseudolist[z1]

  

  # create the real outputs

      econt, cont = setup.compress_continuum(self,ebins, cont_rrc+cont_2ph+cont_bre, const.TOLERANCE, minval=1e-38)
      epseudo, pseudo = setup.compress_continuum(self, ebins, pseudo, const.TOLERANCE, minval=1e-38)

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

      ret['cont'] = numpy.zeros(1, dtype=setup.generate_datatypes(self, 'continuum', npseudo=maxnpseudo, ncontinuum=maxncont))

      ret['cont']['Z'][0] = Z
      ret['cont']['rmJ'][0] = 0
      ret['cont']['N_Cont'][0] = maxncont
      ret['cont']['E_Cont'][0][:maxncont] = econt
      ret['cont']['Continuum'][0][:maxncont] = cont
      ret['cont']['N_Pseudo'][0] = maxnpseudo
      ret['cont']['E_Pseudo'][0][:maxnpseudo] = epseudo
      ret['cont']['Pseudo'][0][:maxnpseudo] = pseudo

      return ret




    def return_ionbal(self, Z, Te, init_pop=False, tau=False,\
                       teunit='K', \
                       filename=False, datacache=False, fast=True,
                       settings= False, debug=False, extrap=True):

      """
      Solve the ionization balance for a element Z.

      Parameters
      ----------
      Z : int
        atomic number of element
      Te : float or array
        electron temperature(s), default in K
      init_pop : float array
        initial population of ions for non-equlibrium calculations. Will be renormalised to 1.
      tau : float or array
        N_e * t for the non-equilibrium ioniziation, in cm^3 s^-1.
      Te_init : float
        initial ionization balance temperature, same units as Te
      teunit : {'K' , 'keV'}
        units of temperatures (default K)
      filename : string
        Can optionally point directly to the file in question, i.e. to look at older data
        look at $HEADAS/../spectral/modelData/eigenELSYMB_v3.0.fits.
        If not set, download from AtomDB FTP site.
      datacache : dict
        Used for caching the data. See description in atomdb.get_data
      fast : bool
        If true, use precalculated eigenvector files to obtain CIE and NEI results

      Returns
      -------
      final_pop : float array
        final populations.

      """

      if fast:
        ionbal = element._solve_ionbal_eigen(self, Z, Te, init_pop=init_pop, tau=tau, \
                       teunit=teunit, \
                       filename=filename, datacache=datacache, debug=debug)
        return ionbal

      else:
        ionbal = element.calc_elem_ionbal(self, Z, Te, tau=tau, init_pop=init_pop, teunit=teunit,\
                               extrap=extrap, settings=settings, datacache=datacache)
        return ionbal



    def _solve_ionbal_eigen(self, Z, Te, init_pop=False, tau=False, \
                       teunit='K', \
                       filename=False, datacache=False, debug=False):
      """
      Solve the ionization balance for element Z using the eigenvector
      approach and files as distributed in XSPEC.

      Parameters
      ----------
      Z : int
        atomic number of element
      Te : float or array
        electron temperature(s), default in K
      init_pop : float array
        initial population of ions for non-equlibrium calculations. Will be renormalised to 1.
      tau : float or array
        N_e * t for the non-equilibrium ioniziation, in cm^3 s^-1.

      teunit : {'K' , 'keV'}
        units of temperatures (default K)
      filename : string
        Can optionally point directly to the file in question, i.e. to look at older data
        look at $HEADAS/../spectral/modelData/eigenELSYMB_v3.0.fits.
        If not set, download from AtomDB FTP site.
      datacache : dict
        Used for caching the data. See description in atomdb.get_data

      Returns
      -------
      final_pop : float array
        final populations.

      """
    #
    #  Version 0.1 Initial Release
    #  Adam Foster 16th September 2015
    #


      kT = util.convert_temp(Te, teunit, 'keV')

      if type(tau)==bool:
        if tau==False:
          cie = True
          init_pop_calc=False

        else:
          raise ValueError("Error: tau should be False, a float, or an array of floats. Received "+repr(Tau))
      else:
        cie = False

      if not cie:
      # if it's not equilibrium, get the initial population
        if isinstance(init_pop, str):
          if init_pop.lower() == 'ionizing':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[0] = 1.0
          elif init_pop.lower() == 'recombining':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[-1] = 1.0
          else:
            raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
             (init_pop))
        elif isinstance(init_pop, float):
      # this is an initial temperature
          kT_init = util.convert_temp(init_pop, teunit, 'keV')
          init_pop_calc = element.return_ionbal(self,Z, kT_init, \
                                            teunit='keV', \
                                            datacache=datacache,fast=True)


        elif isinstance(init_pop, numpy.ndarray) or isinstance(init_pop, list):
          init_pop_calc = init_pop
        elif isinstance(init_pop, dict):
          init_pop_calc = init_pop[Z]
        else:
          raise util.OptionError("Error: invalid type for init_pop")


  # open the eigenvector data file

      if util.keyword_check(filename):
    # we have a filename specified!
        fname = os.path.expandvars(filename)
        if not os.path.isfile(fname):
          print("Specified file %s does not exist. Exiting"%(fname))
          return
        d = pyfits.open(fname)
      else:
        Atomic_data =  ion.datacache(self, Z, z1)
        d = Atomic_data['data'][Z][False]['eigen'] 
        #d = atomdb.get_data(Z, False, 'eigen', datacache=datacache)
      telist = numpy.logspace(4,9,1251)
      kTlist=telist*const.KBOLTZ

  # if we are looking for equilibrium, return the nearest data
      if cie:

        ikTlist = numpy.argsort(numpy.abs(kTlist-kT))
        ite = [min(ikTlist[:2]), max(ikTlist[:2])]
        Tdiff = kTlist[ite[1]] - kTlist[ite[0]]
        if Tdiff > 0.0:
          factorlow = (kTlist[ite[1]]-kT)/Tdiff
          factorhigh = (kT-kTlist[ite[0]])/Tdiff
          equilib = factorlow * d['EIGEN'].data['FEQB'][ite[0]]+\
                factorhigh * d['EIGEN'].data['FEQB'][ite[1]]
        else:
          equilib = d['EIGEN'].data['FEQB'][ite[0]]

    #renormalize
        equilib /= sum(equilib)
        return equilib

  # now do the non-equilibrium data


    # renormalize
    # make Te into a vector
      kT_vec, kT_isvec = util.make_vec(kT)
      tau_vec, tau_isvec = util.make_vec(tau)
      frac_out = numpy.zeros([len(kT_vec),len(tau_vec),Z+1], dtype=float)
      for ikT, kT in enumerate(kT_vec):
        kTindex = numpy.argmin(numpy.abs(kTlist-kT))

        lefteigenvec = numpy.zeros([Z,Z], dtype=float)
        righteigenvec = numpy.zeros([Z,Z], dtype=float)
        if Z==1:
          for i in range(Z):
            for j in range(Z):
              lefteigenvec[i,j] = d['EIGEN'].data['VL'][kTindex]
              righteigenvec[i,j] = d['EIGEN'].data['VR'][kTindex]
        else:
          for i in range(Z):
            for j in range(Z):
              lefteigenvec[i,j] = d['EIGEN'].data['VL'][kTindex][i*Z+j]
              righteigenvec[i,j] = d['EIGEN'].data['VR'][kTindex][i*Z+j]


        work = numpy.array(init_pop_calc[1:] - d['EIGEN'].data['FEQB'][kTindex][1:], dtype=float)

        fspectmp = numpy.matrix(lefteigenvec) * numpy.matrix(work).transpose()

        delt = 1.0

        worktmp = numpy.zeros(Z)

        for itau, ttau in enumerate(tau_vec):
          if Z >1:
            for i in range(Z):
              worktmp[i] = fspectmp[i]*numpy.exp(d['EIGEN'].data['EIG'][kTindex,i]*delt*ttau)

          else:
            worktmp[0] = fspectmp[0]*numpy.exp(d['EIGEN'].data['EIG'][kTindex]*delt*ttau)

          frac = numpy.zeros(Z+1)
          for i in range(Z):
            for j in range(Z):
              frac[i+1] += worktmp[j]*righteigenvec[j][i]
            frac[i+1] += d['EIGEN'].data['FEQB'][kTindex][i+1]

          if debug:
            frac_out[ikT, itau,:] = frac
          frac[frac<0.0] = 0.0

          if sum(frac)> 1.0:
            frac = frac/sum(frac)
          frac[0] = 1-sum(frac[1:])
          if not(debug):
            frac_out[ikT,itau,:]=frac


      if not tau_isvec:
        frac_out = frac_out.sum(1)
      if not kT_isvec:
        frac_out = frac_out.sum(0)


      return frac_out




    def gather_rates(self, Z, Te, N_e, datacache=False, settings=False,\
                 do_la=True, do_ai=True, do_ec=True, do_pc=True,\
                 do_ir=True):


        
        gather_rates={}
        for z1 in range(1,Z+1):
          gather_rates[z1]=ion.gather_rates(self, Z, z1, Te, N_e)

          return gather_rates



    def level_population(self, Z, Te, N_e, settings):
        level_pop={}
        for z1 in range(1,Z+1):
          level_pop[z1] = ion.level_population(self, Z, z1, Te, N_e, settings)

        return level_pop


  







    def calc_elem_ionbal(self, Z, Te, tau=False, init_pop='ionizing', teunit='K',\
                    extrap=True, settings=False, datacache=False):
      """
      Calculate the ionization balance for all the elements in Zlist.

      One of init_pop or Te_init should be set. If neither is set, assume
      all elements start from neutral.

c
      Parameters
      ----------
      Z : int
        nuclear charge to include in calculation (e.g. 8 for oxygen)
      Te : float
        electron temperature in keV or K (default K)
      tau : float
        N_e * t for the non-equilibrium ioniziation (default False, i.e. CIE)
      init_pop : string or float
        if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
        if 'recombining': all recombining from ionized (so[...0,0,1])
        if array: actual populations (e.g. [0, 0.1, 0.3, 0.5, 0.1, 0, 0])
        if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
        if single float : the temperature (same units as Te)
      teunit : {'K' , 'keV'}
        units of temperatures (default K)
      extrap : bool
        Extrappolate rates to values outside their given range. (default False)

      Returns
      -------
      final_pop : array
        final population. E.g. [0.1,0.2,0.3,0.2,0.2,0.0,0.0]

      """

      kT = util.convert_temp(Te, teunit, 'keV')
      if tau==False:
        cie = True
        init_pop_calc=False
      else:
        cie = False
      if not cie:
          # if it's not equilibrium, get the initial population
        if isinstance(init_pop, str):
          if init_pop.lower() == 'ionizing':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[0] = 1.0
          elif init_pop.lower() == 'recombining':
            init_pop_calc = numpy.zeros(Z+1)
            init_pop_calc[-1] = 1.0
          else:
            raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
             (init_pop))
        elif isinstance(init_pop, float):
      # this is an initial temperature
          kT_init = util.convert_temp(init_pop, teunit, 'keV')

      # rerun this routine in equilibrium mode to find the initial ion pop
          init_pop_calc = element.return_ionbal(self,Z, kT_init, \
                                    teunit='keV', \
                                    datacache=datacache,fast=False,
                                    settings = settings, extrap=extrap)


        elif isinstance(init_pop, numpy.ndarray) or isinstance(init_pop, list):
          init_pop_calc = init_pop
        elif isinstance(init_pop, dict):
          init_pop_calc = init_pop[Z]
        else:
          raise util.OptionError("Error: invalid type for init_pop")


  # get the end point population

      ionrate = numpy.zeros(Z, dtype=float)
      recrate = numpy.zeros(Z, dtype=float)
      for z1 in range(1,Z+1):
        ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)

      if cie:
        final_pop = setup.solve_ionbal(self,ionrate, recrate)
      else:
        final_pop = setup.solve_ionbal(self, ionrate, recrate, init_pop=init_pop_calc, tau=tau)

      return final_pop







    def calc_elem_ionbal_delta(self, Z, Te, factor, errortype, tau=False, init_pop='ionizing', teunit='K',\
                    extrap=True, settings=False, datacache=False):
    
        

        kT = util.convert_temp(Te, teunit, 'keV')
        if tau==False:
            cie = True
            init_pop_calc=False
        else:
            cie = False
        if not cie:
            if isinstance(init_pop, str):
                if init_pop.lower() == 'ionizing':
                    init_pop_calc = numpy.zeros(Z+1)
                    init_pop_calc[0] = 1.0
                elif init_pop.lower() == 'recombining':
                    init_pop_calc = numpy.zeros(Z+1)
                    init_pop_calc[-1] = 1.0
                else:
                    raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
                        (init_pop))
            elif isinstance(init_pop, float):
                kT_init = util.convert_temp(init_pop, teunit, 'keV')

      # rerun this routine in equilibrium mode to find the initial ion pop
                init_pop_calc = element.return_ionbal(self, Z, kT_init, \
                                    teunit='keV', \
                                    datacache=datacache,fast=False,
                                    settings = settings, extrap=extrap)


            elif isinstance(init_pop, numpy.ndarray) or isinstance(init_pop, list):
                init_pop_calc = init_pop
            elif isinstance(init_pop, dict):
                init_pop_calc = init_pop[Z]
            else:
                raise util.OptionError("Error: invalid type for init_pop")


  # get the end point population
        
        ionrate = numpy.zeros(Z, dtype=float)
        recrate = numpy.zeros(Z, dtype=float)
        
           
        for z1 in range(1,Z+1):
            if errortype=='ION': 
                ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)
                ionrate[z1-1], recrate[z1-1] = factor*ionrate[z1-1], recrate[z1-1]

            elif errortype=='REC': 
                ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)
                ionrate[z1-1], recrate[z1-1]= ionrate[z1-1], factor*recrate[z1-1]
                

            else:
                ionrate[z1-1], recrate[z1-1] = atomdb.get_ionrec_rate(kT, False, Te_unit='keV', \
                     Z=Z, z1=z1, datacache=datacache, extrap=extrap,\
                     settings=settings)



        
        if cie:
            final_pop = setup.solve_ionbal(self,ionrate, recrate)
        else:
            final_pop = setup.solve_ionbal(self, ionrate, recrate, init_pop=init_pop_calc, tau=tau)

        
        

        return final_pop



    


    def run_apec_element(self, Z, Te, N_e, settings):
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
        z1list = list(range(1, Z+2))
        ionfrac = numpy.ones(len(z1list), dtype=float)

      elif settings['Ionization']=='CIE':
        z1list = list(range(1, Z+2))

    # time to calculate the ionization balance
        if settings['UseIonBalanceTable']:
      # read the ionization balance table
          ionfrac = atomdb._get_precalc_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, Te)
          ionfrac=ionfrac[Z]
        else:
      # calculate the ionization balance
          ionftmp = allelement.calc_full_ionbal(self, Te, 1e14, Te_init=Te, Zlist=[Z], extrap=True)
          ionfrac = ionftmp[Z]

      else:
        print("ERROR: settings['Ionization'] must be CIE or NEI, not %s"%(settings['Ionization']))


      abundfile = atomdb.get_filemap_file('abund',\
                                      Z,\
                                      False,\
                                      fmapfile=settings['FileMap'],\
                                      atomdbroot=os.path.expandvars('$ATOMDB'),\
                                      misc=True)

      abundances = atomdb.get_abundance(abundfile = abundfile, abundset=settings['Abundances'])
      Abund=abundances

  # create placeholders for all the data

      linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
      contlist = {}
      pseudolist = {}

  # now go through each ion and assemble the data
      ebins =  setup.make_vector_nbins(self, settings['LinearGrid'], \
                             settings['GridMinimum'], \
                             settings['GridMaximum'], \
                             settings['NumGrid'])
      ecent = (ebins[1:]+ebins[:-1])/2

      z1_drv_list = numpy.arange(1,Z+2, dtype=int)
      for z1_drv in range(1, Z+2):

        tmplinelist, tmpcontinuum, tmppseudocont = ion.run_apec_ion(self, settings, Z, z1_drv, Te, N_e, ionfrac, Abund)
        linelist = numpy.append(linelist, tmplinelist)
        contlist[z1_drv] = tmpcontinuum
        pseudolist[z1_drv] = tmppseudocont


  # now merge these together.
      if settings['Ionization']=='CIE':

        cieout = element.generate_cie_outputs(self, settings, Z, linelist, contlist, pseudolist)
        return cieout
      
      elif settings['Ionization']=='NEI':
        ionftmp= allelement.calc_full_ionbal(self,Te, 1e14, Te_init=Te, Zlist=[Z], extrap=True)
        ionfrac_nei = ionftmp[Z]
        neiout = element.generate_nei_outputs(self,settings, Z, linelist, contlist, pseudolist, ionfrac_nei)
        return neiout



    def wrap_element_directly(self, ind, Te, N_e, Z):


      wrap_element_directly={}
      for z1 in range(1,Z+2):
        wrap_element_directly[z1]=ion.wrap_ion_directly(self,  ind, Te, N_e, Z, z1)
        print(z1)
      return 

      





    def run_apec_element_delta(self, settings, Te, N_e, Z, factor, errortype):
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
        z1list = list(range(1, Z+2))
        ionfrac = numpy.ones(len(z1list), dtype=float)

      elif settings['Ionization']=='CIE':
        z1list = list(range(1, Z+2))

    # time to calculate the ionization balance
        if settings['UseIonBalanceTable']:
      # read the ionization balance table
          ionfrac = atomdb._get_precalc_ionfrac(os.path.expandvars(settings['IonBalanceTable']), Z, Te)
          ionfrac=ionfrac[Z]
        else:
      # calculate the ionization balance
          ionftmp = allelement.calc_full_ionbal(self, Te, 1e14, Te_init=Te, Zlist=[Z], extrap=True)
          ionfrac = ionftmp[Z]

      else:
        print("ERROR: settings['Ionization'] must be CIE or NEI, not %s"%(settings['Ionization']))


      abundfile = atomdb.get_filemap_file('abund',\
                                      Z,\
                                      False,\
                                      fmapfile=settings['FileMap'],\
                                      atomdbroot=os.path.expandvars('$ATOMDB'),\
                                      misc=True)

      abundances = atomdb.get_abundance(abundfile = abundfile, abundset=settings['Abundances'])
      Abund=abundances

  # create placeholders for all the data

      linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
      contlist = {}
      pseudolist = {}

  # now go through each ion and assemble the data
      ebins =  setup.make_vector_nbins(self, settings['LinearGrid'], \
                             settings['GridMinimum'], \
                             settings['GridMaximum'], \
                             settings['NumGrid'])
      ecent = (ebins[1:]+ebins[:-1])/2

      z1_drv_list = numpy.arange(1,Z+2, dtype=int)
      for z1_drv in range(1, Z+2):

        tmplinelist, tmpcontinuum, tmppseudocont = ion.run_apec_ion_delta(self, settings, Z, z1_drv, Te, N_e, factor, errortype, ionfrac, Abund)
        linelist = numpy.append(linelist, tmplinelist)
        contlist[z1_drv] = tmpcontinuum
        pseudolist[z1_drv] = tmppseudocont


  # now merge these together.
      if settings['Ionization']=='CIE':

        cieout = element.generate_cie_outputs(self, settings, Z, linelist, contlist, pseudolist)
        return cieout
      elif settings['Ionization']=='NEI':
        ionftmp= allelement.calc_full_ionbal(self,Te, 1e14, Te_init=Te, Zlist=[Z], extrap=True)
        ionfrac_nei = ionftmp[Z]
        neiout = element.generate_nei_outputs(self,settings, Z, linelist, contlist, pseudolist, ionfrac_nei)
        return neiout














    def wrap_run_apec_element(self, settings, Te, N_e, ind, Z, writepickle=True, readpickle=True):
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
      
    
      if ind==0:

        if readpickle:
          setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,Te,N_e)

        print("Trying to read in the whole element pickle file %s"%(setpicklefname))
        if os.path.exists(setpicklefname):
          ret = pickle.load(open(setpicklefname,'rb'))
          print("read %s"%(setpicklefname))
          return ret
        else:
          print("Not found. Going to assemble it.")


        linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
        contlist = {}
        pseudolist = {}


        ebins =  setup.make_vector_nbins(self, settings['LinearGrid'], \
                             settings['GridMinimum'], \
                             settings['GridMaximum'], \
                             settings['NumGrid'])

        ecent = (ebins[1:]+ebins[:-1])/2



        for z1_drv in range(1,Z+2):
          
          setpicklefname = "%s_Z_%i_z1_%i_iT_%iiN_%i.pkl"%(settings['OutputFileStem'], Z, z1_drv, Te, N_e)
          
          print("loading %s"%(setpicklefname))
          if not os.path.exists(setpicklefname):
            print("Warning: no such file: %s"%(setpicklefname))
            contlist[z1_drv]={}
            contlist[z1_drv]['rrc']=numpy.zeros(settings['NumGrid'], dtype=float)
            contlist[z1_drv]['twophot']=numpy.zeros(settings['NumGrid'], dtype=float)
            contlist[z1_drv]['brems']=numpy.zeros(settings['NumGrid'], dtype=float)

            pseudolist[z1_drv] = numpy.zeros(settings['NumGrid'], dtype=float)
          else:
            dat= pickle.load(open(setpicklefname,'rb'))
            tmplinelist, tmpcontinuum, tmppseudocont = dat['data']


            nlines = len(tmplinelist)
            print(nlines)

            ngoodlines = sum(numpy.isfinite(tmplinelist['epsilon']))
            if nlines != ngoodlines:
              print("Bad lines found in %s"%(setpicklefname))
            linelist = numpy.append(linelist, tmplinelist)
            for key in list(tmpcontinuum.keys()):
              tmpncont = len(tmpcontinuum[key])
              if tmpncont != sum(numpy.isfinite(tmpcontinuum[key])):
                print("Bad continuum found in %s %s"%(key, setpicklefname), end=' ')
#          if key=='rrc':
#            tmpcontinuum['rrc'][numpy.isnan(tmpcontinuum['rrc'])]=0.0
                print("")

            contlist[z1_drv] = tmpcontinuum

            tmpncont = len(tmppseudocont)
            if tmpncont != sum(numpy.isfinite(tmppseudocont)):
              print("Bad pseudocont found in %s"%( setpicklefname))
            pseudolist[z1_drv] = tmppseudocont

  # now merge these together.
        if settings['Ionization']=='CIE':

          cieout = element.generate_cie_outputs(self,settings, Z, linelist, contlist, pseudolist)
          if writepickle:
            setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,Te,N_e)
            pickle.dump(cieout, open(setpicklefname,'wb'))
            print("wrote %s"%(setpicklefname))
          return cieout
        elif settings['Ionization']=='NEI':
          ionftmp= allelement.calc_full_ionbal(self,te, extrap=True, cie=True, settings=settings, Zlist=[Z])

          ionfrac_nei = ionftmp[Z]
          neiout = element.generate_nei_outputs(self,settings, Z, linelist, contlist, pseudolist, ionfrac_nei)
          if writepickle:
            setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)
            pickle.dump(neiout, open(setpicklefname,'wb'))
            print("wrote %s"%(setpicklefname))
          return neiout


      else:


        if readpickle:
          setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)

          print("Trying to read in the whole element pickle file %s"%(setpicklefname))
          if os.path.exists(setpicklefname):
            ret = pickle.load(open(setpicklefname,'rb'))
            print("read %s"%(setpicklefname))
            return ret
          else:
            print("Not found. Going to assemble it.")


        linelist = numpy.zeros(0, dtype=setup.generate_datatypes(self,'linetype'))
        contlist = {}
        pseudolist = {}
 




  
        ebins =  setup.make_vector_nbins(self, settings['LinearGrid'], \
                             settings['GridMinimum'], \
                             settings['GridMaximum'], \
                             settings['NumGrid'])

        ecent = (ebins[1:]+ebins[:-1])/2
        dens = setup.make_vector(self,settings['LinearDens'], \
                     settings['DensStart'], \
                     settings['DensStep'], \
                     settings['NumDens'])

        ite = ind //len(dens)
        idens = ind%len(dens)

        for z1_drv in range(1,Z+2):
          setpicklefname = "%s_Z_%i_z1_%i_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,z1_drv,ite,idens)
          print("loading %s"%(setpicklefname))
          if not os.path.exists(setpicklefname):
            print("Warning: no such file: %s"%(setpicklefname))
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
            print(nlines)

            ngoodlines = sum(numpy.isfinite(tmplinelist['epsilon']))
            if nlines != ngoodlines:
              print("Bad lines found in %s"%(setpicklefname))
            linelist = numpy.append(linelist, tmplinelist)
            for key in list(tmpcontinuum.keys()):
              tmpncont = len(tmpcontinuum[key])
              if tmpncont != sum(numpy.isfinite(tmpcontinuum[key])):
                print("Bad continuum found in %s %s"%(key, setpicklefname), end=' ')
#          if key=='rrc':
#            tmpcontinuum['rrc'][numpy.isnan(tmpcontinuum['rrc'])]=0.0
                print("")

            contlist[z1_drv] = tmpcontinuum

            tmpncont = len(tmppseudocont)
            if tmpncont != sum(numpy.isfinite(tmppseudocont)):
              print("Bad pseudocont found in %s"%( setpicklefname))
            pseudolist[z1_drv] = tmppseudocont

  # now merge these together.
        if settings['Ionization']=='CIE':

          cieout = element.generate_cie_outputs(self,settings, Z, linelist, contlist, pseudolist)
          if writepickle:
            setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)
            pickle.dump(cieout, open(setpicklefname,'wb'))
            print("wrote %s"%(setpicklefname))
          return cieout
        elif settings['Ionization']=='NEI':
          ionftmp= allelement.calc_full_ionbal(self,te, extrap=True, cie=True, settings=settings, Zlist=[Z])

          ionfrac_nei = ionftmp[Z]
          neiout = element.generate_nei_outputs(self, settings, Z, linelist, contlist, pseudolist, ionfrac_nei)
          if writepickle:
            setpicklefname = "%s_Z_%i_elem_iT_%iiN_%i.pkl"%(settings['OutputFileStem'],Z,ite,idens)
            pickle.dump(neiout, open(setpicklefname,'wb'))
            print("wrote %s"%(setpicklefname))
          return neiout



    
  

        



class allelement():



    
    '''
    Te=float(input("Te="))
    N_e=float(input("N_e="))
  #upper_level=int(input("upper_level="))
  #lower_level=int(input("lower_level="))
    ind=int(input("index (if Te and N_e are known put ind=0)="))
    '''




    

    def __init__(self, Te, N_e, ind):
        
        '''
        Te=self.Te
        Ne=self.N_e
        ind=self.ind
        '''
        

    
    def calc_full_ionbal(self, Te, tau=False, init_pop='ionizing', Te_init=False, Zlist=False, teunit='K',\
                    extrap=True, cie=True, settings=False, datacache=False):
      """
      Calculate the ionization balance for all the elements in Zlist.

      One of init_pop or Te_init should be set. If neither is set, assume
      all elements start from neutral.


      Parameters
      ----------
      Te : float
        electron temperature in keV or K (default K)
      tau : float
        N_e * t for the non-equilibrium ioniziation (default False, i.e. off)
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
        print("*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunits))

  # input checking
      if util.keyword_check(Te_init):
        if teunit.lower() == 'kev':
          kT_init = Te_init*1.0
        elif teunit.lower() == 'ev':
          kT_init = Te_init/1000.0
        elif teunit.lower() == 'k':
          kT_init = Te_init*const.KBOLTZ
        else:
          print("*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
            (teunits))

      if not Zlist:
        Zlist = list(range(1,29))

      if (not util.keyword_check(Te_init)) & (not util.keyword_check(init_pop)) &\
          (not util.keyword_check(cie)):
        print("Warning: you have not specified an initial temperature or "+\
              "ion population: assuming everything is neutral")
        init_pop={}
        for Z in Zlist:
          init_pop[Z] = numpy.zeros(Z+1, dtype=float)
          init_pop[Z][0] = 1.0


      if (util.keyword_check(Te_init)!=False) & (util.keyword_check(init_pop)!=False):
        print("Warning: you have specified both an initial temperature and "+\
          "ion population: using ion population")




      if cie:
    # do some error checking
        if util.keyword_check(Te_init):
          print("Warning: you have specified both Te and Te_init for a CIE calculation. "+\
            "Using Te for ionization balance calculation")


      if (not init_pop) | (cie):
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

          init_pop[Z] = setup.solve_ionbal(self,ionrate, recrate)
      if cie:
    # Exit here if the inital

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


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
    


    




    def wrap_allelement_directly(self, ind, Te, N_e):

      Zlist = list(range(6,9))
      wrap_allelement_directly={}


      for Z in Zlist:
        for z1 in range(1,Z+2):
          wrap_allelement_directly[z1]=ion.wrap_ion_directly(self,  ind, Te, N_e, Z, z1)

      return




    def create_lhdu_cie(self,linedata):

      # sort the data

      tmp = numpy.zeros(len(linedata), dtype=numpy.dtype({'names':['negLambda','Element','Ion'],\
                                                       'formats':[float, float, float]}))
      tmp['Element']= linedata['element']
      tmp['Ion']= linedata['ion']
      tmp['negLambda']= linedata['lambda']*-1

      srt = numpy.argsort(tmp, order=['Element','Ion','negLambda'])
      linedata = linedata[srt]



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
      #print("Created a linelist HDU with %i lines"%(len(tbhdu.data)))
      return tbhdu




    def create_lhdu_nei(self, linedata):

      # sort the data, as required by XSPEC. Cannot sort in reverse order, so hackery ensues.

      tmp = numpy.zeros(len(linedata), dtype=numpy.dtype({'names':['negLambda','Element','Ion'],\
                                                       'formats':[float, float, float]}))
      tmp['Element']= linedata['element']
      tmp['Ion']= linedata['ion']
      tmp['negLambda']= linedata['lambda']*-1

      srt = numpy.argsort(tmp, order=['Element','Ion','negLambda'])
      linedata = linedata[srt]

      cols = []
      cols.append(pyfits.Column(name='Lambda', format='1E', unit="A", array=linedata['lambda']))
      cols.append(pyfits.Column(name='Lambda_Err', format='1E', unit="A", array=linedata['lambda_err']))
      cols.append(pyfits.Column(name='Epsilon', format='1E', unit="photons cm^3 s^-1", array=linedata['epsilon']))
      cols.append(pyfits.Column(name='Epsilon_Err', format='1E', unit="photons cm^3 s^-1", array=linedata['epsilon_err']))
      cols.append(pyfits.Column(name='Element', format='1J',  array=linedata['element']))
      cols.append(pyfits.Column(name='Elem_drv', format='1J',  array=linedata['element']))
      cols.append(pyfits.Column(name='Ion', format='1J',  array=linedata['ion']))
      cols.append(pyfits.Column(name='Ion_drv', format='1J',  array=linedata['ion_drv']))
      cols.append(pyfits.Column(name='UpperLev', format='1J', array=linedata['upperlev']))
      cols.append(pyfits.Column(name='LowerLev', format='1J',  array=linedata['lowerlev']))

      coldefs = pyfits.ColDefs(cols)
      tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
      return tbhdu

    


    def wrap_run_apec(self, settings, readpickle=False, writepickle=False):
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
      
      
      # need to parse the IncAtoms parameter - this defines which atoms
      # we will include
      #
      # we will transfer this to a "Zlist" parameter, in the form of
      # nuclear charges

    #  if len(Zlist)==0:
      #Zlist = settings['Zlist']
      Zlist=[6,7,8]

      print("I will be running Z=", Zlist)
  # run for each element, temperature, density

      lhdulist = []
      chdulist = []

      seclHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=setup.generate_datatypes(self,'lineparams'))
      seccHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=setup.generate_datatypes(self,'cocoparams'))


      for iTe in range(settings['NumTemp']):

        te = setup.make_vector(self,settings['LinearTemp'], \
                     settings['TempStart'], \
                     settings['TempStep'], \
                     settings['NumTemp'])[iTe]

        if settings['TempUnits']=='keV':
          te /= const.KBOLTZ


        for iDens in range(settings['NumDens']):
          dens = setup.make_vector(self, settings['LinearDens'], \
                         settings['DensStart'], \
                         settings['DensStep'], \
                         settings['NumDens'])[iDens]

      # AT THIS POINT, GENERATE SHELL WHICH WILL GO IN THE HDU OF CHOICE

          if settings['Ionization']=='CIE':
            linedata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'linelist_cie'))
            cocodata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'continuum', ncontinuum=0, npseudo=0))
          if settings['Ionization']=='NEI':
            linedata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'linetype'))
            cocodata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'continuum', ncontinuum=0, npseudo=0))

          for Z in Zlist:
            print("Calling run_apec_element for Z=%i Te=%e dens=%e at %s"%(Z, te, dens, time.asctime()))
            dat = element.wrap_run_apec_element(self, settings, Te, N_e, ind, Z,readpickle=readpickle)
        # append this data to the output
        #pickle.dump(dat, open('dump_%i.pkl'%(Z),'wb'))
            linedata = numpy.append(linedata, dat['lines'])
            cocodata = setup.continuum_append(self, cocodata, dat['cont'])
            print("Z=%i, nlines=%i"%(Z, len(dat['lines'])))


      # now make an HDU for all of this
          if settings['Ionization']=='CIE':
            LHDUdat = allelement.create_lhdu_cie(self,linedata)
          elif settings['Ionization']=='NEI':
#        pickle.dump(linedata,open('argh.pkl','wb'))
            LHDUdat = allelement.create_lhdu_nei(self, linedata)

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
          CHDUdat = setup.create_chdu_cie(self,cocodata)

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
          tot_emiss = setup.calc_total_coco(self,cocodata, settings)

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
          seclHDU = setup.create_lparamhdu_cie(self,seclHDUdat)
          seclHDU.header['EXTNAME']=('PARAMETERS','name of this binary table extension')
          seclHDU.header['HDUCLASS']=('Proposed OGIP','Proposed OGIP standard')
          seclHDU.header['HDUCLAS1']=('LINE MODEL','line emission spectra model')
          seclHDU.header['HDUCLAS2']=('Parameters','extension containing parameter info')
          seclHDU.header['HDUVERS1']=('1.0.0','version of format')


          seccHDU = setup.create_cparamhdu_cie(self,seccHDUdat)
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
          setup.generate_apec_headerblurb(self,settings, tmplhdulist, tmpchdulist)

  # write out results
          tmplhdulist.writeto('%s_line.fits'%(fileroot), clobber=True, checksum=True)

          if settings['Ionization']=='CIE':
            tmpchdulist.writeto('%s_coco.fits'%(fileroot), clobber=True, checksum=True)
          elif settings['Ionization']=='NEI':
            tmpchdulist.writeto('%s_comp.fits'%(fileroot), clobber=True, checksum=True)





    def run_apec(self, settings):
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
      #settings = parse_par_file(fname)

      # need to parse the IncAtoms parameter - this defines which atoms
      # we will include
      #
      # we will transfer this to a "Zlist" parameter, in the form of
      # nuclear charges

      #Zlist = settings['Zlist']
      Zlist=[6,7,8]
      print("I will be running Z=", Zlist)
      # run for each element, temperature, density

      lhdulist = []
      chdulist = []

      seclHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=setup.generate_datatypes(self,'lineparams'))
      seccHDUdat = numpy.zeros(settings['NumTemp']*settings['NumDens'], \
                          dtype=setup.generate_datatypes(self,'cocoparams'))


      for iTe in range(settings['NumTemp']):

        te = setup.make_vector(self,settings['LinearTemp'], \
                     settings['TempStart'], \
                     settings['TempStep'], \
                     settings['NumTemp'])[iTe]

        if settings['TempUnits']=='keV':
          te /= const.KBOLTZ


        for iDens in range(settings['NumDens']):
          dens = setup.make_vector(self,settings['LinearDens'], \
                         settings['DensStart'], \
                         settings['DensStep'], \
                         settings['NumDens'])[iDens]

      # AT THIS POINT, GENERATE SHELL WHICH WILL GO IN THE HDU OF CHOICE

          if settings['Ionization']=='CIE':
            linedata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'linelist_cie'))
            cocodata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'continuum', ncontinuum=0, npseudo=0))
          if settings['Ionization']=='NEI':
            linedata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'linetype'))
            cocodata = numpy.zeros(0,dtype=setup.generate_datatypes(self,'continuum', ncontinuum=0, npseudo=0))

          for Z in Zlist:
            print("Calling run_apec_element for Z=%i Te=%e dens=%e at %s"%(Z, te, dens, time.asctime()))
            dat = element.run_apec_element(self, settings, te, dens, Z)
        # append this data to the output
            linedata = numpy.append(linedata, dat['lines'])
            cocodata = setup.continuum_append(self,cocodata, dat['cont'])



      # now make an HDU for all of this
          if settings['Ionization']=='CIE':
            LHDUdat = allelement.create_lhdu_cie(self,linedata)
          elif settings['Ionization']=='NEI':
            LHDUdat = allelement.create_lhdu_nei(self,linedata)
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
          CHDUdat = setup.create_chdu_cie(self,cocodata)

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
          tot_emiss = setup.calc_total_coco(self,cocodata, settings)

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
      seclHDU = setup.create_lparamhdu_cie(self,seclHDUdat)
      seclHDU.header['EXTNAME']=('PARAMETERS','name of this binary table extension')
      seclHDU.header['HDUCLASS']=('Proposed OGIP','Proposed OGIP standard')
      seclHDU.header['HDUCLAS1']=('LINE MODEL','line emission spectra model')
      seclHDU.header['HDUCLAS2']=('Parameters','extension containing parameter info')
      seclHDU.header['HDUVERS1']=('1.0.0','version of format')


      seccHDU = setup.create_cparamhdu_cie(self,seccHDUdat)
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
        tmpchdulist.writeto('%s_coco_1.fits'%(fileroot), clobber=True, checksum=True)
      elif settings['Ionization']=='NEI':
        tmpchdulist.writeto('%s_comp_1.fits'%(fileroot), clobber=True, checksum=True)





class setup:

    






    def generate_datatypes(self, dtype, npseudo=0, ncontinuum=0):
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
                                 'elem_drv',\
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
      elif dtype =='linelist_cie_spectrum':
        ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Ion', \
                                 'UpperLev',\
                                 'LowerLev'],\
                        'formats':[numpy.float32,\
                                   numpy.float32,\
                                   numpy.float32,\
                                   numpy.float32,\
                                   numpy.int32,\
                                   numpy.int32,\
                                   numpy.int32,\
                                   numpy.int32]})
      elif dtype =='linelist_nei_spectrum':
        ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Elem_drv',\
                                 'Ion', \
                                 'Ion_drv', \
                                 'UpperLev',\
                                 'LowerLev'],\
                        'formats':[numpy.float32,\
                                   numpy.float32,\
                                   numpy.float32,\
                                   numpy.float32,\
                                   numpy.int32,\
                                   numpy.int32,\
                                   numpy.int32,\
                                   numpy.int32,\
                                   numpy.int32,\
                                   numpy.int32]})
      elif dtype == 'linetype_cap':
        ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Ion', \
                                 'Elem_Drv',\
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
      elif dtype == 'ecdat':
    # Electron collisional data
        ret = numpy.dtype({'names':['lower_lev','upper_lev','coeff_type',\
                                'min_temp','max_temp', 'temperature', \
                                'effcollstrpar','inf_limit','reference'],\
                       'formats':[int, int, int, \
                                  float, float, (float,20), \
                                  (float,20), float, '|S20']})


      else:
        print("Unknown dtype %s in generate_datatypes"%(dtype))
      return ret

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------



    def calc_cascade_population(self, matrixA, matrixB):

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
            print("failed again")

        # look for levels with no way to ground. Put in a -1rate
            for i in range(ma.shape[0]):
              if ma[i,i]>= 0.0:
                ma[i,i]=-1.0

            try:
              popn = numpy.linalg.solve(ma,mb)
            except:
              print('triple fail')
              print(ma)
              print(mb)
              raise

  #check
      soln = numpy.allclose(numpy.dot(ma, popn), mb)

      if soln==False:
        print("ERROR Solving population matrix!")
      popn=numpy.append(numpy.array([0.0]), popn)
      return popn




    def make_vector(self, linear, minval, step, nstep):
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





#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
    def make_vector_nbins(self, linear, minval, maxval, nstep):
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




    def compress_continuum(self, xin, yin, tolerance, minval = 0.0):
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

      try:
        from pyatomdb import liblinapprox
      except OSError:
        print("Unable to open liblinapprox. Hopefully you are on Readthedocs")
        return
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



    def continuum_append(self,a,b):
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

      c = numpy.zeros(nlines, dtype=setup.generate_datatypes(self,'continuum',npseudo=npseudomax, ncontinuum=ncontmax))
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



    def create_chdu_cie(self,cocodata):

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


    def create_lparamhdu_cie(self,linedata):

      cols = []
      cols.append(pyfits.Column(name='kT', format='1E', unit="keV", array=linedata['kT']))
      cols.append(pyfits.Column(name='EDensity', format='1E', unit="cm**-3", array=linedata['EDensity']))
      cols.append(pyfits.Column(name='Time', format='1E', unit="s", array=linedata['Time']))
      cols.append(pyfits.Column(name='Nelement', format='1J', array=linedata['Nelement']))
      cols.append(pyfits.Column(name='Nline', format='1J',  array=linedata['Nline']))

      coldefs = pyfits.ColDefs(cols)
      tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
      return tbhdu

    def create_cparamhdu_cie(self,cocodata):
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

    def calc_total_coco(self, cocodata, settings):
      """
      Calculate the total emission in erg cm^3 s^-1

      """
      emiss = 0.0

      ebins = setup.make_vector_nbins(self,settings['LinearGrid'], \
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


    def generate_apec_headerblurb(self, settings, linehdulist, cocohdulist):
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
              Atomic_data =  ion.datacache(self, Z, z1)
              a = Atomic_data['data'][Z][z1][ftype] 
              #a=atomdb.get_data(Z,z1,ftype)
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


    def solve_ionbal(self, ionrate, recrate, init_pop=False, tau=False):
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
      b = numpy.zeros(Z+1, dtype=float)
      a = numpy.zeros((Z+1,Z+1), dtype=float)


      for iZ in range(0,Z):
        a[iZ,iZ] -= (ionrate[iZ])
        a[iZ+1,iZ] += (ionrate[iZ])

        a[iZ,iZ+1] += (recrate[iZ])
        a[iZ+1,iZ+1] -= (recrate[iZ])

  # conservation of population
      for iZ in range(0,Z+1):
        a[0,iZ]=1.0
      b[0]=1.0

      eqpop=numpy.linalg.solve(a,b)

      eqpop[eqpop<0] = 0.0

      eqpop[0] = max([1.0-sum(eqpop[1:]), 0])

      if do_equilib == True:
        return eqpop

  # now the NEI part
  #remake matrix a
      Z=len(ionrate)+1
      ndim=Z
      AA = numpy.zeros((ndim-1,ndim-1), dtype=float)
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


      w,v = numpy.linalg.eig(AA)


  # now copy VR to AA
      AA=v

  # b_in is difference betwen initial and final populations
      b_in=init_pop[1:]-eqpop[1:]

  # solve for initial conditions
      b_out=numpy.linalg.solve(AA,b_in)

  # include exponential decay term
      C = b_out*numpy.exp( w*tau )

  # get change
      G = numpy.dot(v,C)

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
              print("Error in paramater file %s. Item %s: %s is not in allowed parameter list %s"%\
                (fname, name, value, minval))
              print(options)
            else:
              data[name] = value
          else:
            data[name] = value
        elif dtype =='b':
      #boolean, yes or no.
          if not value in ['yes','no']:
            print("Error in paramater file %s. Item %s: should be boolean yes or no, supplied value is %s"%\
              (fname, name, value))
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
              print("Error in paramater file %s. Item %s: should be > %e, supplied value is %e"%\
              (fname, name, minval, value))

          if len(maxval) > 0:
            maxval = float(maxval)
            if value > maxval:
              print("Error in paramater file %s. Item %s: should be < %e, supplied value is %e"%\
              (fname, name, maxval, value))

          data[name] = value
        elif dtype=='i':
      # integer
          value = int(value)
          if len(minval) > 0:
            minval = int(minval)
            if value < minval:
              print("Error in paramater file %s. Item %s: should be > %i, supplied value is %i"%\
              (fname, name, minval, value))

          if len(maxval) > 0:
            maxval = int(maxval)
            if value > maxval:
              print("Error in paramater file %s. Item %s: should be < %i, supplied value is %i"%\
              (fname, name, maxval, value))
          data[name] = value
        else:
          print("Error in paramater file %s. Item %s: unknown data type %s"%\
          (fname, name, dtype))

    # now some massaging of the parameters
      if 'BremsType' in list(data.keys()):
        if data['BremsType']=='Hummer':
          data['BremsType']=const.HUMMER
        elif data['BremsType']=='Kellogg':
          data['BremsType']=const.KELLOGG
        elif data['BremsType']=='Relativistic':
          data['BremsType']=const.RELATIVISTIC
        else:
          print("UNKNOWN BREMS TYPE: %s" %(data['BremsType']))

      if 'NEIMinEpsilon' in list(data.keys()):
          tmp = data['NEIMinEpsilon'].split(',')
          tmp2 = []
          for i in tmp:
            tmp2.append(float(i))
          data['NEIMinEpsilon']=tmp2

      if 'NEIMinFrac' in list(data.keys()):
          tmp = data['NEIMinFrac'].split(',')
          tmp2 = []
          for i in tmp:
            tmp2.append(float(i))
          data['NEIMinFrac']=tmp2

      if 'IncAtoms' in list(data.keys()):
        elsymblist = data['IncAtoms'].split(',')
        Zlist=[]
        for iel in elsymblist:
          Zlist.append(atomic.elsymb_to_Z(iel))
        data['Zlist']=Zlist

      if not ('WriteIon' in list(data.keys())):
        data['WriteIon'] = False


      return data






















        

