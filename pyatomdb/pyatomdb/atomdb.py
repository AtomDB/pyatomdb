import os, datetime, numpy, re, time
import util, atomic, spectrum
import astropy.io.fits as pyfits
from scipy import stats, integrate

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def write_filemap(d, filemap, atomdbroot=''):
#
# filemap - name of filemap file to writeto
#


  # open file for writing
  a=open(filemap, 'w')

  # cycle through all the "misc" files
  for i in range(0, len(d['misc'])):
    # search for the atomdbroot substring, if it's present, replace with
    # $ATOMDB. If atomdbroot not set, use ATOMDB environment variable instead
    # If that's not set, do nothing.
    
    if len(atomdbroot) > 0:
      d['misc'][i] = re.sub(atomdbroot,'$ATOMDB',d['misc'][i])
    else:
      # see if the ATOMDB environment variable is set
      if 'ATOMDB' in os.environ.keys():
        d['misc'][i] = re.sub(os.environ['ATOMDB'],'$ATOMDB',d['misc'][i])
    a.write('%2i %2i %2i ' % (d['misc_type'][i], 0, -1)+\
                              d['misc'][i]+'\n')

  tlist = ['ir','lv','la','ec','pc','dr','em','pi','ai','ci']
  for i in range(0, len(d['z0'])):

    for t in range(1, len(tlist)+1):
      tt = tlist[t-1]
      if d[tt][i] != '':

        # search for the atomdbroot substring, if it's present, replace with
        # $ATOMDB. If atomdbroot not set, use ATOMDB environment variable instead
        # If that's not set, do nothing.

        if len(atomdbroot) > 0:
          d[tt][i] = re.sub(atomdbroot,'$ATOMDB',d[tt][i])
        else:
          # see if the ATOMDB environment variable is set
          if 'ATOMDB' in os.environ.keys():
            d[tt][i] = re.sub(os.environ['ATOMDB'],'$ATOMDB',d[tt][i])


    

        y = t
        if t>=10: y+=10  
        a.write('%2i %2i %2i ' % (y, d['z0'][i], d['z1'][i])+\
                              d[tt][i]+'\n')

  a.close()
  return


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def extract_n(conf_str):
  # to get the maximum n from the configuration string, e.g. "1s2 3p1", n=3

  #split by ' '
  tmp = re.split(' ',conf_str)
  n_out = -1
  for i in tmp:
    n_tmp = int( re.search('^[0-9]+', i).group(0))
    if n_tmp > n_out:
      n_out = n_tmp
  return n_out
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def read_filemap(filemap=False, atomdbroot=False):
  """
  filemap - name of filemap file to read. If zero length, use "$ATOMDB/filemap"
  atomdbroot - replace any $ATOMDB in the file names with this. If not provided
              will use "ATOMDB" environment variable instead
  """              

  if not(atomdbroot):
    try:
      atomdbroot = os.environ['ATOMDB']
    except KeyError:
      print '*** ERROR: atomdbroot not set, use environment variable ATOMDB or '+\
          'pass atomdbroot to read_filemap'
      return False
  
  if not(filemap):
    filemap = atomdbroot+'/filemap'
  if not os.path.exists(filemap):
    print "*** ERROR: Filemap %s does not exist, cannot fetch file" %(filemap)
    return False

  f = open(filemap,'r')

  filedtype = numpy.dtype({'names':['z0','z1','ec','lv','ir','pc','dr',\
                                    'la','em','pi','ai'],\
                           'formats':[int, int, '|S160','|S160','|S160',\
                                      '|S160','|S160','|S160','|S160',\
                                      '|S160','|S160','|S160']})

  miscdtype = numpy.dtype({'names':['misc_type','file'],\
                           'formats':[int, '|S160']})

  ret = {}
  ret['ionfiles'] = numpy.zeros(0, dtype=filedtype)
  ret['miscfiles'] = numpy.zeros(0, dtype=miscdtype)
  
  
  for i in f:
    splt = i.split()
    fname = re.sub('\$ATOMDB',atomdbroot,splt[3])
    z0_tmp = int(splt[1])
    z1_tmp = int(splt[2])
    
    if z0_tmp < 1:
      # in this case, we have a "misc" datatype, not corresponding to a particular ion
      misc_type = int(splt[0])
      misc_file = fname
      ret['miscfiles'] = numpy.append(ret['miscfiles'], numpy.array((misc_type, misc_file),\
                                                                    dtype=miscdtype))
      
    else:
      j = numpy.where((ret['ionfiles']['z0']==z0_tmp) & \
                      (ret['ionfiles']['z1']==z1_tmp))[0]
      if len(j)==0:
        j = len(ret['ionfiles'])

        ret['ionfiles'] = numpy.append(ret['ionfiles'], numpy.zeros(1,\
                                                                    dtype=filedtype))
        ret['ionfiles']['z0'][j]=z0_tmp
        ret['ionfiles']['z1'][j]=z1_tmp
      else:
        j = j[0]
        
      if int(splt[0]) == 1:
        ret['ionfiles']['ir'][j] = fname
      if int(splt[0]) == 2:
        ret['ionfiles']['lv'][j] = fname
      if int(splt[0]) == 3:
        ret['ionfiles']['la'][j] = fname
      if int(splt[0]) == 4:
        ret['ionfiles']['ec'][j] = fname
      if int(splt[0]) == 5:
        ret['ionfiles']['pc'][j] = fname
      if int(splt[0]) == 6:
        ret['ionfiles']['dr'][j] = fname
      if int(splt[0]) == 7:
        ret['ionfiles']['em'][j] = fname
      if int(splt[0]) == 8:
        ret['ionfiles']['pi'][j] = fname
      if int(splt[0]) == 9:
        ret['ionfiles']['ai'][j] = fname
#      if int(splt[0]) == 10:
#        cilist[j] = fname


#  ret={}
#  ret['z0'] = numpy.array(z0)
#  ret['z1'] = numpy.array(z1)
#  ret['ec'] = numpy.array(eclist, dtype='|S160')
#  ret['lv'] = numpy.array(lvlist, dtype='|S160')
#  ret['ir'] = numpy.array(irlist, dtype='|S160')
#  ret['pc'] = numpy.array(pclist, dtype='|S160')
#  ret['dr'] = numpy.array(drlist, dtype='|S160')
#  ret['la'] = numpy.array(lalist, dtype='|S160')
#  ret['em'] = numpy.array(emlist, dtype='|S160')
#  ret['pi'] = numpy.array(pilist, dtype='|S160')
#  ret['ai'] = numpy.array(ailist, dtype='|S160')
#  ret['ci'] = numpy.array(cilist, dtype='|S160')
#  ret['misc'] = numpy.array(misc, dtype='|S160')
#  ret['misc_type'] = numpy.array(misc_type)
#
  return ret

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def ci_younger(Te, c):
  """Te is an array in Kelvin
     c is the ionrec_par
     returns ionization rate in cm^3 s^-1"""

  KBOLTZ = 8.617385e-8  #keV/K

  T_eV = 1.e3*KBOLTZ*Te
  x = c[0]/T_eV
  ci = numpy.zeros(len(Te), dtype=float)
  for ite in range(len(Te)):
  
    if (x[ite] <= 30.0):
      f1_val = f1_fcn(x[ite])
      ci[ite] = (numpy.exp(-x[ite])/x[ite])* \
                ( c[1]*( 1 - x[ite]*f1_val ) +\
                  c[2]*( 1 + x[ite] - x[ite]*(x[ite]+2)*f1_val )+\
                  c[3]*f1_val + \
                  c[4]*x[ite]*f2_fcn(x[ite]) )

  ci *= 6.69e-07/(T_eV*numpy.sqrt(T_eV))
  return ci
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def ea_mazzotta(Te, c, par_type):
  """Te is an array in Kelvin
     c is the ionrec_par
     par_type is the number denoting the type of the parameter
     returns excitation-autoionization rate in cm^3 s^-1"""

  KBOLTZ = 8.617385e-8  #keV/K

  T_eV = 1.e3*KBOLTZ*Te
  ea = numpy.zeros(len(Te), dtype=float)
  if par_type == 61:
    # 1985A&AS...60..425A (Arnaud & Rothenflug)

    # Note: c[1] = b/Zeff^2 in Arnaud & Rothenflug
    y = c[0]/T_eV
    for iy in range(len(y)):
      if (y[iy] < 50.0): # otherwise, rate is effectively zero
        yf1 = y[iy]*f1_fcn(y[iy])
        G = 2.22*f1_fcn(y[iy]) + 0.67*(1 - yf1) + 0.49*yf1 + 1.2*y[iy]*(1 - yf1)
        ea[iy] = c[2] * c[1]  * 1.92e-07 * numpy.exp(-y[iy]) * G/numpy.sqrt(T_eV[iy])

  
  elif par_type == 62:
  # 1985A&AS...60..425A (Arnaud & Rothenflug)
    I_ea = c[0]
    a = c[1]
    b = c[2]
    y =  I_ea/T_eV
    for iy in range(len(y)):
      if (y[iy] < 50.0):
        f1_val = f1_fcn(y[iy])
        ea[iy] = 6.69e7 * a * (I_ea/sqrt(T_eV[iy])) * exp(-y[iy]) * \
             (1.0 + b*f1_val - c[3]*y[iy]*f1_val - \
             c[4]*0.5*(y[iy] - y[iy]*y[iy] + y[iy]*y[iy]*y[iy]*f1_val))

  elif par_type == 63:
  # 1998A&AS..133..403M (Mazzotta etal 1998) */
    y=c[0]/T_eV
    for iy in range(len(y)):
      if (y[iy] < 50.0):
        f1_val = f1_fcn(y[iy]);
        ea[iy] = (6.69e7/numpy.sqrt(T_eV[iy]))*numpy.exp(-y[iy])*1.e-16*\
             (c[1]+\
              c[2]*(1-y[iy]*f1_val)+\
              c[3]*(1-y[iy]*(1-y[iy]*f1_val))+\
              c[4]*(1-0.5*(y[iy] - y[iy]*y[iy] + y[iy]*y[iy]*y[iy]*f1_val))+\
              c[5]*f1_val)

  return ea
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def rr_verner(Te, c):

  tt0=numpy.sqrt(Te/c[2])
  tt1=numpy.sqrt(Te/c[3])
  rr = c[0]/( (tt0) * ((tt0+1)**(1-c[1])) * ((1+tt1)** (1+c[1])))

  return rr

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def ea_mazzotta_iron(T_eV, c):
  ret = numpy.zeros(len(T_eV), dtype = float)

  for i in range(len(T_eV)):
    y = c[0]/T_eV[i]
    if (y < 50.0):
      f1_val = f1_fcn(y)
      ea = (6.69e7/numpy.sqrt(T_eV[i]))* numpy.exp(-y) * 1.0e-16 * \
           (c[1]+\
            c[2]*(1-y*f1_val)+\
            c[3]*(1-y*(1-y*f1_val)) +\
            c[4]*(1-0.5*(y - y*y + y*y*y*f1_val)) +\
            c[5]*f1_val)
    else:
      ea = 0.0
    ret[i] = ea
  return ret 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def rr_shull(Te, c):
  KBOLTZ = 8.617385e-8  #keV/K

  rr = c[0]* (Te/1e4)**(-c[1])
  return rr
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def f1_fcn(x):

#  /* Analytic approx.  to f_1(x) = exp(x) \int_1^{\infty} (exp(-xt)/t) dt */
#  /* Commented-out version from                                           */
#  /*   Arnaud \& Rothenflug, 1985 A&A, Appendix B, p. 436                 */
#  /* Version used from Mazzotta's original code.                          */

  result=0.0

  a=numpy.array( [-0.57721566,0.99999193,-0.24991055,
     0.05519968,-0.00976004,0.00107857])
  b=numpy.array( [8.5733287401,18.0590169730,8.6347608925,0.2677737343,
     9.5733223454,25.6329561486,21.0996530827,3.9584969228])

  if (x < 0.0 ) :
    errmess("f1_fcn","Negative values of x not allowed")

  if (x < 1.0) :
    result=(a[0]+a[1]*x+a[2]*x*x + a[3]* x**3 + a[4]*x**4 +
      a[5]*x**5-numpy.log(x))*numpy.exp(x)
  else:
    if (x < 1.0e9) :
      result=(x**4+b[0]*x**3+b[1]*x**2+b[2]*x+b[3])/ \
  (x**4+b[4]*x**3+b[5]*x**2+b[6]*x+b[7])/x
    else:
      result=1.0/x

  return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def f2_fcn(x):

  # Analytic approximation to                                     */
  # f_2(x) = exp(x) \int_1^{\infty} (log(t) exp(-xt)/t) dt      */
  # Taken from Arnaud \& Rothenflug, 1985 A&A, Appendix B, p. 436 */


  pc= numpy.array([1.0000e0 ,2.1658e2 ,2.0336e4 ,1.0911e6 ,3.7114e7 ,
       8.3963e8 ,1.2889e10,1.3449e11,9.4002e11,4.2571e12,
       1.1743e13,1.7549e13,1.0806e13,4.9776e11,0.0000e0])
  qc= numpy.array([1.0000e0 ,2.1958e2 ,2.0984e4 ,1.1517e6 ,4.0349e7 ,
       9.4900e8 ,1.5345e10,1.7182e11,1.3249e12,6.9071e12,
       2.3531e13,4.9432e13,5.7760e13,3.0225e13,3.3641e12])

  if (x > 0.27):
    p=0.0
    q=0.0
    xp = 1
    for i in range(0,15):
      p = p+pc[i]/xp
      q = q+qc[i]/xp
      xp = x*xp

    result=((p/q)/x)/x
  else:
    result=(((numpy.log(x)*numpy.log(x))/2.0)+(0.57722*numpy.log(x)))+1.0

  return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def dr_mazzotta(Te, c):
  KBOLTZ = 8.617385e-8
  T_eV = 1.e3*KBOLTZ*Te
  dr = (c[0]/(T_eV**1.5))*(c[5]*numpy.exp(-c[1]/T_eV) +\
                               c[6]*numpy.exp(-c[2]/T_eV) +\
                               c[7]*numpy.exp(-c[3]/T_eV) +\
                               c[8]*numpy.exp(-c[4]/T_eV))
  return dr


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def addline(xbins, yvals, wv, amp, dx):
  
  # create the cdf
  cdf = stats.norm(wv,dx).cdf(xbins)
  vec = (cdf[1:]-cdf[:-1])*amp
#  tmpsumyin = sum(yvals)
#  print "sum(yin)", tmpsumyin
#  print vec
  return vec
  yvals = yvals + vec
#  print sum(vec), amp, sum(yvals), sum(yvals)-sum(vec)
#  zzz=raw_input()
  return yvals


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def addline2(xbins, wv, amp, dx):
  
  # create the cdf
  cdf = stats.norm(wv,dx).cdf(xbins)
  vec = (cdf[1:]-cdf[:-1])*amp
#  yvals = yvals + vec
  return vec
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_ionfrac(ionbalfile, element, te, ion=-1):
  """ Gets the ionization fraction of a given ion at a given Te
      Assuming ionization equilibrium
      ionbalfile : string, location of ionization balance file
      element    : int, z0 of element (e.g. 6 for carbon)
      ion        : if provided, z+1 of ion (e.g. 5 for O V)
                   if omitted, returns ionization fraction for all ions
                   of element
      returns
      ionization fraction, or if ion not specified, list of values
  """
  
  dat = get_ionbal(ionbalfile, element)
  it = numpy.argmin(abs(dat['te']-te))
  print it, dat['te'][it]
  if ion==-1:
    ion = range(0,element+1)
  frac = dat['frac'][it,ion]
  
  return frac
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_ionbal(ionbalfile, element, ion=-1):

  # ionbalfile : string, location of ionization balance file
  # element    : int, z0 of element (e.g. 6 for carbon)
  # ion        : if provided, z+1 of ion (e.g. 5 for O V, 3 for Ne III)
  #              if omitted, returns ionization balance for all ions of element
  #
  # returns dictionary with keys:
  # ['te']     : float array, electron temperature
  # ['dens']   : float array, electron densities
  # ['time']   : float array, elapsed times, if in original file
  # ['frac']   : float array, fractional abundance of ion
  #            : if ion not set, returns ['frac'][te,z]
  

  # check ion makes sense for the given element
  if ion ==0:
    print 'ion should be in range 1 to z0, in this case '+repr(z0)
    print 'return whole element ionization balance in stead'
  if ion > element + 1:
    print 'ion should be in range 1 to z0, in this case '+repr(z0)
    print 'ERROR, returning -1'
    return -1
  


  a=pyfits.open(ionbalfile)
  
#  find data type
  if a[1].header['HDUCLASS']=='ION_BAL':
  
    # find index for correct ion
    
    i = 0
    
    z_element = a[1].data.field('z_element')[0]
    
    # find the correct column index for the element
    for j in range(1, element):
      if j in z_element: i = i + j+1
  
    # if ion specified, find the column for the ion
    if ion <= 0:
      i = i+numpy.arange(element+1)
    else:
      i = i + ion-1
  
    out = a[1].data.field('x_ionpop')[:,i]
    
    # prepare return dictionary
    ret={}
    ret['frac'] = out
    ret['te'] = a[1].data.field('temperature')
    ret['dens'] = a[1].data.field('density')

  else:
    
    i = 0
    
    z_element = a[1].data.field('element')[0]
    z_element = numpy.array(z_element)
    if a[1].data.field('NZ')[0]==1:
      j = a[1].data.field('indx')[0]
    else:
      j = numpy.where(z_element==element)[0][0]
      j = a[1].data.field('indx')[0][j]
      
    
#    for j in range(1, element):
#      if j in z_element: i = i + j+1
    
    if ion <= 0:
      i = j+numpy.arange(element+1)
    else:
      i = j + ion-1

    out = a[1].data.field('Ionbal')[:,i]
    ret={}
    ret['frac'] = out
    ret['te'] = a[1].data.field('temperature')
    ret['dens'] = a[1].data.field('density')
    if "Time" in a[1].data.names:
      ret['time'] = a[1].data.field('time')
    
  return ret
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def get_filemap_file(ftype, z0, z1, fmapfile=False,atomdbroot=False, quiet=False, misc=False):
  """
  get_filemap_file: find the correct file from the database for atomic data of type ftype
  for ion with nuclear charge z0 and ioncharge+1 = z1
  
  inputs:
  ftype - file type of interest:
    'ir' = ionization & recombination data
    'lv' = energy levels
    'la' = wavelength and transition probabilities (lambda & a-values)
    'ec' = electron collision rates
    'pc' = proton collision rates
    'dr' = dielectronic recombination satellite line information
    'ai' = autoionization rate data
    'pi' = XSTAR photoionization data
    'em' = emission feature data (currently unused)
    
  z0  - element atomic number (=6 for C+4)
  z1  - ion charge +1 (=5 for C+4)
  
  kwargs:
  fmapfile = specific filemap to use. Otherwise defaults to atomdbroot+'/filemap'
  atomdbroot = location of ATOMDB database. Defaults to ATOMDB environment variable.
               all $ATOMDB in the filemap will be expanded to this value
  quiet = if true, suppress warnings about files not being present for certain ions
  misc = if requesting "misc" data, i.e. the Bremsstrahlung inputs, use this. This is
         for non ion-specific data.

  returns: the filename for the relevant file, with all $ATOMDB expanded. If 
           no file exists, returns zero length string.
  """
  
  fmap=read_filemap(filemap=fmapfile, atomdbroot=atomdbroot)
  
  
  if misc:
    # looking for type 10,11 or 12 data
    if ftype in [10,11,12,13]:
      i = numpy.where(fmap['miscfiles']['misc_type']==ftype)[0]
      if len(i)==0:
        print "Error: file type: %i not found in filemap %s" %(ftype,fmapfile)
        ret = ''        
      else:
        ret = fmap['miscfiles']['file'][i[0]]
  else:
    i = numpy.where((fmap['ionfiles']['z0']==z0)&(fmap['ionfiles']['z1']==z1))[0]
    ret=''
    if len(i)==0:
      if not quiet :
        print "WARNING: there is no data for the ion "+\
               adbatomic.spectroscopic_name(z0,z1-1)
      ret=''
           
    if len(i)>1:
      print "ERROR: there are multiple entries for the ion "+\
             adbatomic.spectroscopic_name(z0,z1-1)
      ret=''
    
    if len(i)==1:
      i=i[0]
      ftypel = ftype.lower()
    
      if not ftypel in fmap['ionfiles'].dtype.names:
      
        print "Error: invalid file type: "+ftype
        ret = ''
      
      else:
        ret = fmap['ionfiles'][ftypel][i]
        if len(ret)==0:
          if not quiet :
            print "WARNING: no data of type "+ftype+" exists for ion "+\
                adbatomic.spectroscopic_name(z0,z1-1)

  return ret
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def make_level_descriptor(lv):
  cfgstr = lv['elec_config']
  #check if L & S are not -1, to determine if a term symbol can be made
  
  lsymb = 'spdfghiklmnopqrstuvwxyz'

  if lv['lev_deg'] % 2 != 0:
    Jsym = "%i" % ((lv['lev_deg']-1.0)/2.0)
  else:
    Jsym = "%i/%1i" % (lv['lev_deg']-1,2)

  if ((lv['s_quan']==-1) | (lv['l_quan']==-1)):
    # jj coupling
    LSsym=''
  else:
    LSsym= "^{%i}%s" % ((2*lv['s_quan'])+1,lsymb[lv['l_quan']].upper())
  
  ret = cfgstr+' '+LSsym+"_{"+Jsym+"}"
  return ret
    
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_ionpot(z0, z1, filemap, atomdbroot=''):
  
  irfile = get_filemap_file(filemap,'IR', z0, z1, atomdbroot=atomdbroot)
  
  irdat = pyfits.open(irfile)
  
  ionpot = irdat[1].header['IONPOT']
  
  return ionpot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def get_level_details(level, z0=-1, z1=-1, filename='', \
                      filemap='', atomdbroot=''):
  """Function returns the details in the level file for the specified 
  level. LV file can be specified by filename, or by filemap, z0, z1"""
  
  if filename=='':
    # get the filename from the other variables
    if z0 < 0:
      print "Error in get_level_details: must specify filename or "+\
            "z0, z1 and filemap"
      return -1

    if z1 < 0:
      print "Error in get_level_details: must specify filename or "+\
            "z0, z1 and filemap"
      return -1

    if filemap =='':
      print "Error in get_level_details: must specify filename or "+\
            "z0, z1 and filemap"
      return -1
    filename = get_filemap_file(filemap, 'LV', z0, z1,\
                                atomdbroot=atomdbroot)

  # open the file
  lvdat = pyfits.open(filename)
  
  ret = {}
  ret['elec_config'] = lvdat[1].data.field('elec_config')[level-1]
  ret['energy']      = lvdat[1].data.field('energy')[level-1]
  ret['e_error']     = lvdat[1].data.field('e_error')[level-1]
  ret['n_quan']      = lvdat[1].data.field('n_quan')[level-1]
  ret['l_quan']      = lvdat[1].data.field('l_quan')[level-1]
  ret['s_quan']      = lvdat[1].data.field('s_quan')[level-1]
  ret['lev_deg']     = lvdat[1].data.field('lev_deg')[level-1]
  ret['phot_type']   = lvdat[1].data.field('phot_type')[level-1]
  ret['phot_par']    = lvdat[1].data.field('phot_par')[level-1]
  ret['energy_ref']  = lvdat[1].data.field('energy_ref')[level-1]
  ret['phot_ref']    = lvdat[1].data.field('phot_ref')[level-1]
  ret['string']      = "%s, E=%.4f, S=%.1f, L=%i, lev_deg=%i, ref=%s"%\
                      (lvdat[1].data.field('elec_config')[level-1],\
                       lvdat[1].data.field('energy')[level-1],\
                       lvdat[1].data.field('s_quan')[level-1],\
                       lvdat[1].data.field('l_quan')[level-1],\
                       lvdat[1].data.field('lev_deg')[level-1],\
                       lvdat[1].data.field('energy_ref')[level-1])
                       
                         
  
  return ret 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_abundance(abundfile, abundset, element=[-1]):

  a = pyfits.open(abundfile)
  if element[0]==-1:
    element = range(1,29)
    
  ind = numpy.where(a[1].data.field('Source')==abundset)[0]
  if len(ind)==0:
    print "Invalid Abundance Set chosen: select from ", \
          a[1].data.field('Source')
    return -1
  ret = {}
  for z0 in element:
    elsymb = adbatomic.z0toelsymb(z0)
    ret[z0]=10**(a[1].data.field(elsymb)[ind[0]])/1e12
  return ret
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_emissivity(linefile, elem, ion, upper, lower, kT=[-1], \
                   hdu = [-1], kTunits='keV'):
                     
  a = pyfits.open(linefile)
  # find the blocks of interest
  if ((kT[0] != -1) & (hdu[0] != -1)):
    print "Error in get_emissivity: provide either the hdu numbers "+\
          "(starting from 2) or the temperatures"
    return -1
  interpKt = False
  if (kT[-1] != -1):
    interpKt = True
    hdulist = []
    if kTunits.lower() == "k":
      kT_keV = numpy.array(kT)/(1e3*11604.5)
    elif kTunits.lower()== "kev":
      kT_keV = numpy.array(kT)
    else:
      print "Error in get_emissivity: kTunits must be K or keV"
      return -1
        
    ikT_min = numpy.where(a[1].data.field('kT') >= min(kT_keV))[0]
    ikT_max = numpy.where(a[1].data.field('kT') <= max(kT_keV))[0]
    if ((len(ikT_min)==0) | (ikT_min[0] == 0)):
      print "Error in get_emissivity: kT=%e out of range %e:%e keV" %\
             (min(kT_keV), a[1].data.field('kT')[0], \
              a[1].data.field('kT')[-1])
      return -1

    if ((len(ikT_max)==0) | (ikT_max[-1] == len(a[1].data.field('kT'))-1)):
      print "Error in get_emissivity: kT=%e out of range %e:%e keV" %\
             (max(kT_kev), a[1].data.field('kT')[0], \
              a[1].data.field('kT')[-1])
      return -1
    ikT = numpy.arange(ikT_min[0]-1, ikT_max[-1]+2)
  elif (hdu[0] != -1):
    ikT = numpy.array(hdu)-2
    kT_keV = a[1].data.field('kT')[ikT]
  else:
    ikT = numpy.arange(len(a[1].data.field('kT')), dtype=int)
    kT_keV = a[1].data.field('kT')[ikT]
  
  # ok. Get the numbers
  emiss_grid = numpy.zeros(len(ikT), dtype=float)
  for iemiss, i in enumerate(ikT):
    ii = i+2
    j = numpy.where((a[ii].data.field('element')==elem) &\
                    (a[ii].data.field('ion')==ion) &\
                    (a[ii].data.field('upperlev')==upper) &\
                    (a[ii].data.field('lowerlev')==lower))[0]

    if len(j) > 0:
      emiss_grid[iemiss] = a[ii].data.field('epsilon')[j[0]]

  
  # ok, now get the numbers on a nice grid
  if (interpKt):
    emiss = numpy.interp(kT_keV, a[1].data.field('kT')[ikT], emiss_grid)
  else:
    emiss = emiss_grid
  return kT_keV, emiss

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def interp_rate(Te, npar, Te_grid,ionrec_par):
  from scipy.interpolate import interp1d

  Te_log = numpy.log(Te)
  Te_grid_log = numpy.log(Te_grid[:npar])
  ionrec_par_log = numpy.log(ionrec_par[:npar])
  #print Te_grid_log
  #print ionrec_par_log
  tmp = interp1d(Te_grid_log, ionrec_par_log, kind='cubic')
  #print Te_log
  ret = numpy.exp(tmp(Te_log))
#  print Te_grid, ionrec_par
  #print ret
#  zzz=raw_input()
  return ret

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_ion_lines(linefile, z0, z1, fullinfo=False):
  a = pyfits.open(linefile)
  kT = a[1].data.field('kT')
  nlines = numpy.zeros(len(kT), dtype=int)
  for ikT in range(len(kT)):
    iikT = ikT + 2

    j = numpy.where((a[iikT].data.field("element") == z0) &\
                    (a[iikT].data.field("ion") == z1))[0]
    nlines[ikT] = len(j)
    if (fullinfo):
      print "kT = %.2e, nlines = %i" % (kT[ikT], nlines[ikT])
      if nlines[ikT] > 0:
        for jj in j:
          print "%5i: %10f %.4e" % (jj, \
                                    a[iikT].data.field('lambda')[jj], \
                                    a[iikT].data.field('epsilon')[jj])
                                    
  return nlines

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_line_emissivity(linefile, z0, z1, upind, loind,ion_drv=False, elem_drv=False):
  a = pyfits.open(linefile)
  kT = a[1].data.field('kT')
  dens = a[1].data.field('eDensity')
  time = a[1].data.field('time')
  
  epsilon = numpy.zeros(len(kT), dtype=float)
  for ikT in range(len(kT)):
    iikT = ikT + 2
    if ion_drv:
      j = numpy.where((a[iikT].data.field("element") == z0) &\
                      (a[iikT].data.field("ion") == z1) &\
                      (a[iikT].data.field("ion_drv") == ion_drv) &\
                      (a[iikT].data.field("UpperLev") == upind) &
                      (a[iikT].data.field("LowerLev") == loind))[0]
    else:
      j = numpy.where((a[iikT].data.field("element") == z0) &\
                      (a[iikT].data.field("ion") == z1) &\
                      (a[iikT].data.field("UpperLev") == upind) &
                      (a[iikT].data.field("LowerLev") == loind))[0]
    if len(j) > 0:
      epsilon[ikT] += sum(a[iikT].data.field("Epsilon")[j])
  ret = {}
  ret['kT'] = kT
  ret['dens'] = dens
  ret['time'] = time
  ret['epsilon'] = epsilon
  
  return ret
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_burgess_tully_transition_type(lolev, uplev, Aval):

  ll=lolev['l_quan']
  lu=uplev['l_quan']
  sl=lolev['s_quan']
  su=uplev['s_quan']
  degl=lolev['lev_deg']
  degu=uplev['lev_deg']
  dE = (uplev['energy']-lolev['energy'])/13.6058
  if Aval==0:
    S=0.0
  else:
    if dE<1e-10: dE=1e-10
    S = 3.73491e-10*degu*Aval/(dE**3)
  FIN = dE*S/(3.0*degl)
  FBIG = 0.01
  FZERO = 1e-4
  if su==sl:
    if (abs(lu-ll)<=1) & (FIN>=FBIG):
      ret = 1
    else:
      ret = 2
  else:
    if ((FIN>FZERO) & (FIN < FBIG)):
      ret = 4
    else:
      ret = 3
  return ret

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_burgess_tully_extrap(bttype, lolev, uplev, Aval, Tarr, om, TTarg):
  from scipy import interpolate
  dE = (uplev['energy']-lolev['energy'])/13.6058
  Tarr_ryd = (Tarr/11604.5)/13.6058
  TTarg_ryd = (TTarg/11604.5)/13.6058
  if dE<1e-10: dE=1e-10
  
  C = 1.5
  e = 1.0
  if bttype ==1:
    x = 1 - numpy.log(C)/(numpy.log((Tarr_ryd/dE) +C))
    x = numpy.append(0,x)
    x = numpy.append(x,1)
    y = om/(numpy.log((Tarr_ryd/dE)+1))
    y = numpy.append(0,y)
    S = 3.73491e-10*uplev['lev_deg']*Aval/(dE**3)
    FIN = dE*S/(3.0*lolev['lev_deg'])
    y = numpy.append(y, 4*lolev['lev_deg']*FIN/dE)
    
    xtarg = 1 - numpy.log(C)/(numpy.log((TTarg_ryd/dE) +C))
  elif bttype ==2:
    x = (Tarr_ryd/dE) /((Tarr_ryd/dE)+C)
    x = numpy.append(0,x)
    x = numpy.append(x,1)
    y = om
    y0 = 0
    y1 = 0
    y = numpy.append(y0,y)
    y = numpy.append(y,y1)

    xtarg = (TTarg_ryd/dE) /((TTarg_ryd/dE)+C)

  elif bttype ==3:
    x = (Tarr_ryd/dE) /((Tarr_ryd/dE)+C)
    x = numpy.append(0,x)
    x = numpy.append(x,1)
    
    y = ((Tarr_ryd/dE)+1)*om
    y = numpy.append(y[0],y)
    y = numpy.append(y,0)

    xtarg = (TTarg_ryd/dE) /((TTarg_ryd/dE)+C)
  elif bttype ==4:
    x = 1 - numpy.log(C)/(numpy.log((Tarr_ryd/dE) +C))
    x = numpy.append(0,x)
    x = numpy.append(x,1)
    
    y = om/(numpy.log( (Tarr_ryd/dE)+C))
    
    y = numpy.append(0,y)
    S = 3.73491e-10*uplev['lev_deg']*Aval/(dE**3)
    FIN = dE*S/(3.0*lolev['lev_deg'])
    y = numpy.append(y, 4*lolev['lev_deg']*FIN/dE)
    
    xtarg = 1 - numpy.log(C)/(numpy.log((TTarg_ryd/dE) +C))
  # now    
  
  ytarg = numpy.interp(xtarg,x,y)
  
  
  
  if bttype ==1:
    upsout = ytarg*(numpy.log((TTarg_ryd/dE)+C))
  elif bttype ==2:
    upsout = ytarg
  elif bttype ==3:
    upsout = ytarg/(TTarg_ryd/dE+1)
  elif bttype ==4:
    upsout = ytarg * numpy.log((TTarg_ryd/dE)+C)
  
  
  if numpy.isscalar(upsout):
    if upsout < 0:
      upsout = 0.0
  else:
    upsout[upsout < 0] = 0.0
  return upsout

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def get_bt_approx(om, Tin, Tout, uplev, lolev, levdat,ladat):
  # Tin, Tout in K
  # om is effective collision strength
#  print " In BT"
  i = numpy.where((ladat[1].data['upper_lev']==uplev) & \
                    (ladat[1].data['lower_lev']==lolev))[0]
  if len(i)==0:
    Aval = 0.0
  else:
    Aval = ladat[1].data['einstein_a'][i[0]]
                    
  # find the type:
  bttype = get_burgess_tully_transition_type(levdat[1].data[lolev],\
                                             levdat[1].data[uplev],\
                                             Aval)

  
    # do the extrappolation etc
  btval = get_burgess_tully_extrap(bttype, \
                                   levdat[1].data[lolev], \
                                   levdat[1].data[uplev], \
                                   Aval, \
                                   Tin, \
                                   om, \
                                   Tout)

#  print bttype, btval
#  print Tin, om
  return btval

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def calc_maxwell_rates(coll_type, min_T, max_T, Tarr, \
                       om, dE, T, Z, degl, degu, quiet=False, \
                       levdat=False, ladat=False, \
                       lolev=False, uplev=False, \
                       force_extrap=False, did_extrap=False):

  import calc_maxwell_rates_constants as const
  from scipy.special import expn
  from scipy import interpolate
  
  xs = numpy.array([0.00, 0.25, 0.50, 0.75, 1.00])
  xs9 = numpy.array([0.00, 0.125, 0.25, 0.375, 0.50, 0.675, 0.80, 0.925, 1.00])
#  print Tarr, om, coll_type
  calc_type = -1
  chi  = dE / (const.KBOLTZ*T)
#  print "chi",chi
  if dE != 0.0:
    chi_inv = (const.KBOLTZ*T) / dE
  else:
    chi_inv = 1e6

  logT = numpy.log10(T)
  st = 0

  if ((ladat!=False) & (force_extrap==True)):
    fill_value=numpy.nan
  else:
    fill_value=0.0

#  if ((min(T) < min_T) | (max(T) > max_T)):
#    if not quiet:
#      print "WARNING: one of T is outside of valid range %e to %e K" %(min_T, max_T)
#      print " T=", T, "coll_type = ", coll_type
#    return False;


  #print "Coll_type=%i"%(coll_type)


  if (coll_type == const.BURGESSTULLY):
    expint_2 = expn(1,chi)*numpy.exp(chi)
    upsilon = om[0] +\
              om[1]*chi*expint_2 +\
              om[2]*chi*(1-chi*expint_2) +\
              om[3]*(chi/2.0)*(1-chi*(1-chi*expint_2)) +\
              om[4]*expint_2
    calc_type = E_UPSILON

  elif coll_type in [const.CHIANTI_1,\
                     const.CHIANTI_2,\
                     const.CHIANTI_3,\
                     const.CHIANTI_4,\
                     const.CHIANTI_5,\
                     const.CHIANTI_6]:
    c = om[0]
#    print "CHIANTI6", om, c, chi
    
    
    if (coll_type == const.CHIANTI_1):
      st = 1 - (numpy.log(c)/numpy.log(chi_inv+c))
    elif (coll_type == const.CHIANTI_2):
      st = chi_inv/(chi_inv+c)
    elif (coll_type == const.CHIANTI_3):
      st = chi_inv/(chi_inv+c)
    elif (coll_type == const.CHIANTI_4):
      st = 1 - (numpy.log(c)/numpy.log(chi_inv+c))
    elif (coll_type == const.CHIANTI_5):
      st = chi_inv/(chi_inv+c)
    elif (coll_type == const.CHIANTI_6):
      st = chi_inv/(chi_inv+c)
                     
    
#    upsilon = interpolate.interp1d(xs, om[1:6], kind='linear', \
#          bounds_error=False, fill_value=0.0)(st)

    y2 = prep_spline_atomdb(xs, om[1:6], 5)
    upsilon = calc_spline_atomdb(xs, om[1:6], y2, 5, st)

#    upsilon = interpolate.interp1d(xs, om[1:6], kind='linear', \
#          bounds_error=False, fill_value=0.0)(st)

    if (coll_type == const.CHIANTI_1):
      upsilon *= numpy.log(chi_inv + const.M_E)
    elif (coll_type == const.CHIANTI_2):
      pass
    elif (coll_type == const.CHIANTI_3):
      upsilon /= (chi_inv+1)
    elif (coll_type == const.CHIANTI_4):
      upsilon *= numpy.log(chi_inv + c)
    elif (coll_type == const.CHIANTI_5):
      upsilon /= chi_inv
    elif (coll_type == const.CHIANTI_6):
      upsilon = 10.0**upsilon
      
    if (coll_type == const.CHIANTI_6):
      calc_type = const.P_UPSILON
    else:
      calc_type = const.E_UPSILON



  elif coll_type in [const.CHIANTI4_1,\
                     const.CHIANTI4_2,\
                     const.CHIANTI4_3,\
                     const.CHIANTI4_4,\
                     const.CHIANTI4_5,\
                     const.CHIANTI4_6]:
    c = om[1]

    if (coll_type == const.CHIANTI4_1):
      st = 1 - (numpy.log(c)/numpy.log(chi_inv+c))
    elif (coll_type == const.CHIANTI4_2):
      st = chi_inv/(chi_inv+c)
    elif (coll_type == const.CHIANTI4_3):
      st = chi_inv/(chi_inv+c)
    elif (coll_type == const.CHIANTI4_4):
      st = 1 - (numpy.log(c)/numpy.log(chi_inv+c))
    elif (coll_type == const.CHIANTI4_5):
      st = chi_inv/(chi_inv+c)
    elif (coll_type == const.CHIANTI4_6):
      st = chi_inv/(chi_inv+c)
                     
#    print xs
#    print om[2:2+om[0]]
#    print st
#    print "COLL_TYPE: %i"%(coll_type)
#    print "om[0]= ", om[0]
    if int(om[0]) == 5:
#      upsilon = interpolate.interp1d(xs, om[2:2+om[0]], kind='cubic',\
#       bounds_error=False, fill_value=0.0)(st)
      y2 = prep_spline_atomdb(xs, om[2:2+om[0]], 5)
      upsilon = calc_spline_atomdb(xs, om[2:2+om[0]], y2, 5, st)
    elif int(om[0])== 9:
      upsilon = interpolate.interp1d(xs9, om[2:2+om[0]], kind='cubic',\
       bounds_error=False, fill_value=0.0)(st)
      y2 = prep_spline_atomdb(xs9, om[2:2+om[0]], 9)
      upsilon = calc_spline_atomdb(xs9, om[2:2+om[0]], y2, 9, st)

    if (coll_type == const.CHIANTI4_1):
      upsilon *= numpy.log(chi_inv + const.M_E)
    elif (coll_type == const.CHIANTI4_2):
      pass
    elif (coll_type == const.CHIANTI4_3):
      upsilon /= (chi_inv+1)
    elif (coll_type == const.CHIANTI4_4):
      upsilon *= numpy.log(chi_inv + c)
    elif (coll_type == const.CHIANTI4_5):
      upsilon /= chi_inv
    elif (coll_type == const.CHIANTI4_6):
      upsilon = 10.0**upsilon
      
    if (coll_type == const.CHIANTI4_6):
      calc_type = const.P_UPSILON
    else:
      calc_type = const.E_UPSILON

  elif (coll_type == const.SGC_1): # S-type transitions, He-like
    upsilon = calc_sampson_s(om, Z, T)
    calc_type = const.E_UPSILON
  elif (coll_type == const.SGC_2): # P-type transitions, He-like
    upsilon = calc_sampson_p(om, Z, T)
    calc_type = const.E_UPSILON
  elif (coll_type == const.SGC_3): # S-type transitions, H-like
    upsilon = calc_sampson_h(om, Z, T)
    calc_type = const.E_UPSILON

  elif (coll_type == const.KATO_NAKAZAKI_1):
    upsilon = calc_kato(1, om, Z, T)
    calc_type = const.E_UPSILON

  elif (coll_type == const.KATO_NAKAZAKI_2):
    upsilon = calc_kato(2, om, Z, T)
    calc_type = const.E_UPSILON

  elif (coll_type == const.PROTON_BT):
    a = om[0]
    b = om[1]
    c = om[2]
#    /* 0.8 for the np/ne ratio, */
    if ((logT > om[3]) & ( logT < om[4])):
      rate_coeff=0.8*(10**(a + b*(logT) + c*logT*logT) )

    calc_type = const.P_RATE_COEFF

  
  elif ((coll_type >= const.INTERP_E_UPSILON) & \
      (coll_type <= const.INTERP_E_UPSILON + const.MAX_UPS)):
    
    N_interp = coll_type - const.INTERP_E_UPSILON
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
                     
    upsilon = numpy.zeros(len(T),dtype=float)
    upsilon[:]=numpy.nan

    if len(it) > 0:
      if sum(om<0)>0:
        om[om<0]=0.0
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]+1e-30), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
      
      upsilon[it] = numpy.exp(tmp)-1e-30
      # fix the nans
    inan = numpy.isnan(upsilon)
    if sum(inan)>0:
      if force_extrap:
        upsilon[inan]=get_bt_approx(om[:N_interp], Tarr[:N_interp], T[inan], uplev, lolev, levdat,ladat)
        did_extrap=True
      else:
        upsilon[inan]= 0.0
    
    
    calc_type = const.E_UPSILON

  elif ((coll_type >= const.INTERP_I_UPSILON) & \
      (coll_type <= const.INTERP_I_UPSILON + const.MAX_UPS)):
    
    N_interp = coll_type - const.INTERP_I_UPSILON
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      if sum(om<0)>0:
        om[om<0]=0.0
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]+1e-30), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)-1e-30
    calc_type = const.EI_UPSILON


  elif ((coll_type >= const.INTERP_E_UPS_OPEN) & \
      (coll_type <= const.INTERP_E_UPS_OPEN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_UPS_OPEN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    upsilon[:]=numpy.nan
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)

    # fix the nans
    inan = numpy.isnan(upsilon)
    if sum(inan)>0:
      if force_extrap:
        upsilon[inan]=get_bt_approx(om[:N_interp], Tarr[:N_interp], T[inan], uplev, lolev, levdat,ladat)
        did_extrap=True
      else:
        upsilon[inan]= 0.0
    calc_type = const.E_UPSILON

  elif ((coll_type >= const.INTERP_E_UPS_INC_MIN) & \
      (coll_type <= const.INTERP_E_UPS_INC_MIN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_UPS_INC_MIN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    upsilon[:]=numpy.nan
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)
      # fix the nans
    inan = numpy.isnan(upsilon)
    if sum(inan)>0:
      # do the BT extrappolations
      if force_extrap:
        upsilon[inan]=get_bt_approx(om[:N_interp], Tarr[:N_interp], T[inan], uplev, lolev, levdat,ladat)
        did_extrap=True
      else:
        upsilon[inan]= 0.0
    calc_type = const.E_UPSILON

  elif ((coll_type >= const.INTERP_E_UPS_INC_MAX) & \
      (coll_type <= const.INTERP_E_UPS_INC_MAX + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_UPS_INC_MAX
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    upsilon[:]=numpy.nan
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)
      # fix the nans
    inan = numpy.isnan(upsilon)
    if sum(inan)>0:
      # do the BT extrappolations
      if force_extrap:
        upsilon[inan]=get_bt_approx(om[:N_interp], Tarr[:N_interp], T[inan], uplev, lolev, levdat,ladat)
        did_extrap=True
      else:
        upsilon[inan]= 0.0
    calc_type = const.E_UPSILON

  elif ((coll_type >= const.INTERP_P_UPSILON) & \
      (coll_type <= const.INTERP_P_UPSILON + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_UPSILON
    upsilon = numpy.zeros(len(T),dtype=float)
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)
    calc_type = const.P_UPSILON

  elif ((coll_type >= const.INTERP_P_UPS_OPEN) & \
      (coll_type <= const.INTERP_P_UPS_OPEN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_UPS_OPEN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)
    calc_type = const.P_UPSILON

  elif ((coll_type >= const.INTERP_P_UPS_INC_MIN) & \
      (coll_type <= const.INTERP_P_UPS_INC_MIN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_UPS_INC_MIN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)
    calc_type = const.P_UPSILON

  elif ((coll_type >= const.INTERP_P_UPS_INC_MAX) & \
      (coll_type <= const.INTERP_P_UPS_INC_MAX + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_UPS_INC_MAX
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    upsilon = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      upsilon[it] = numpy.exp(tmp)
    calc_type = const.P_UPSILON

  elif ((coll_type >= const.INTERP_E_RATE_COEFF) & \
      (coll_type <= const.INTERP_E_RATE_COEFF + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_RATE_COEFF
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = const.E_RATE_COEFF

    rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
    calc_type = const.E_RATE_COEFF

  elif ((coll_type >= const.INTERP_E_RATE_OPEN) & \
      (coll_type <= const.INTERP_E_RATE_OPEN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_RATE_OPEN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = const.E_RATE_COEFF

  elif ((coll_type >= const.INTERP_E_RATE_INC_MIN) & \
      (coll_type <= const.INTERP_E_RATE_INC_MIN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_RATE_INC_MIN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = const.E_RATE_COEFF

  elif ((coll_type >= const.INTERP_E_RATE_INC_MAX) & \
      (coll_type <= const.INTERP_E_RATE_INC_MAX + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_E_RATE_INC_MAX
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = E_RATE_COEFF

  elif ((coll_type >= const.INTERP_P_RATE_COEFF) & \
      (coll_type <= const.INTERP_P_RATE_COEFF + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_RATE_COEFF
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmpfn = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]+1e-30), \
                                 bounds_error=False, \
                                 fill_value=fill_value)

      for iit in it:
        if T[iit] in Tarr[:N_interp]:
          rate_coeff[iit]=om[Tarr[:N_interp]==T[iit]]
        else:
          rate_coeff[iit]=numpy.exp(tmpfn(numpy.log(T[iit])))-1e-30
      
    calc_type = const.P_RATE_COEFF

  elif ((coll_type >= const.INTERP_P_RATE_OPEN) & \
      (coll_type <= const.INTERP_P_RATE_OPEN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_RATE_OPEN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = const.P_RATE_COEFF

  elif ((coll_type >= const.INTERP_P_RATE_INC_MIN) & \
      (coll_type <= const.INTERP_P_RATE_INC_MIN + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_RATE_INC_MIN
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = const.P_RATE_COEFF

  elif ((coll_type >= const.INTERP_P_RATE_INC_MAX) & \
      (coll_type <= const.INTERP_P_RATE_INC_MAX + const.MAX_UPS)):
    N_interp = coll_type - const.INTERP_P_RATE_INC_MAX
    it = numpy.where((T>=min(Tarr[:N_interp])) &\
                     (T<=max(Tarr[:N_interp])))[0]
    rate_coeff = numpy.zeros(len(T),dtype=float)
    if (T > min_T) :
      rate_coeff = interpol_huntd(N_interp, Tarr, om, T)
    if len(it) > 0:
      tmp = interpolate.interp1d(numpy.log(Tarr[:N_interp]), \
                                 numpy.log(om[:N_interp]), \
                                 bounds_error=False, \
                                 fill_value=fill_value)(numpy.log(T[it]))
                                 
      rate_coeff[it] = numpy.exp(tmp)
    calc_type = const.P_RATE_COEFF

  if (calc_type == -1):
    print "ERROR: undefined collision type %i" %(coll_type)
    return False

  #print "upsilon:", upsilon
#------------------------------------

  if (calc_type == const.E_UPSILON):
    
    #negative upsilon is unphysical
    upsilon[upsilon < 0] = 0.0
    
    exc_rate = numpy.zeros(len(chi), dtype=float)
    dex_rate = numpy.zeros(len(chi), dtype=float)
#    print numpy.exp(-chi) 
#    print numpy.sqrt(T)
    for i in range(len(chi)):
      if (chi[i] < const.MAX_CHI):
        exc_rate[i] = const.UPSILON_COLL_COEFF * upsilon[i] * \
                   numpy.exp(-chi[i]) / (numpy.sqrt(T[i])*degl)
        dex_rate[i] = const.UPSILON_COLL_COEFF * upsilon[i] / (numpy.sqrt(T[i])*degu)
  

  if (calc_type == const.EI_UPSILON):
    
    #negative upsilon is unphysical
    upsilon[upsilon < 0] = 0.0
    
    exc_rate = numpy.zeros(len(chi), dtype=float)
    dex_rate = numpy.zeros(len(chi), dtype=float)
#    print numpy.exp(-chi) 
#    print numpy.sqrt(T)
    # this is the same as the electron-impact version, but with extra factor of pi
    for i in range(len(chi)):
      if (chi[i] < const.MAX_CHI):
        exc_rate[i] = const.UPSILON_COLL_COEFF * upsilon[i] * \
                   numpy.exp(-chi[i]) / (numpy.pi* numpy.sqrt(T[i])*degl)
   

  if (calc_type == const.P_UPSILON):
    upsilon[upsilon < 0] = 0.0
    
    exc_rate = numpy.zeros(len(upsilon), dtype=float)
    dex_rate = numpy.zeros(len(upsilon), dtype=float)
    print "Can't calculate collision strength for protons."

  
  if ((calc_type == const.P_RATE_COEFF)|(calc_type == const.E_RATE_COEFF)):
    #rate_coeff=numpy.zeros(len(chi), dtype=float)
    rate_coeff[rate_coeff < 0] = 0.0
    exc_rate = rate_coeff
    dex_rate = rate_coeff*numpy.exp(chi)*(degl*1.0/degu)
#  print "exc_rate:", exc_rate
  
#  zzz=raw_input()
  
  if sum(numpy.isnan(numpy.asarray(exc_rate)))>0:
    print "ERROR: calculated excitation rate an NaN"
    print coll_type, min_T, max_T, Tarr, om, dE, T, Z, degl, degu
  if sum(numpy.isnan(numpy.asarray(dex_rate)))>0:
    print "ERROR: calculated excitation rate an NaN"
    print coll_type, min_T, max_T, Tarr, om, dE, T, Z, degl, degu

  if force_extrap:
    return exc_rate, dex_rate, did_extrap
  
  else:  
    return exc_rate, dex_rate


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_sampson_s(om, Z,Te):
  from scipy.special import expn
  import calc_maxwell_rates_constants as const
#  /* These routines come from Sampson, Goett, & Clark, ADNDT 29, 467 */
  
  result=0.0
  
  
  Z2gamma=0.0
  Z2gamma_e=0.0

  dE     = om[0]
  a1_gam = om[1]
  a1e_gam= om[2]
  z2s_h  = om[3]
  a2     = om[4]
  c0     = om[5]
  c1     = om[6]
  c2     = om[7]
  a2e    = om[8]
  cere   = om[9]
  cere1  = om[10]
  re     = int(om[11])
  s      = om[12]
  se     = om[13]
  
  kT = const.KBOLTZ*Te
  y = dE/kT

  a1 = a2+1.0
  a1y = a1*y
  E1y = expn(1,y)
  Ery = expn(1,a1y)
  Er1y= expn(2,a1y)
  term = c1*Ery + c2*Er1y/a1
  if ((a1_gam != 0.0) & (term > 0)):
    Z2gamma = c0 + 1.333*z2s_h*E1y*numpy.exp(y) + y*numpy.exp(a1y)*term
  else:
    Z2gamma = c0 + 1.333*z2s_h*E1y*numpy.exp(y)


  a1 = a2e+1
  a1y = a1*y
  Ery = expn(re,a1y)
  Er1y = expn(re+1,a1y)

  term = cere*Ery/(a1**(re-1)) + cere1*Er1y/(a1**re)
  if ((a1e_gam != 0.0) & (term > 0)):
    Z2gamma_e = y*numpy.exp(a1y)*term
  else:
    Z2gamma_e = 0.0

  Zeff = Z - s
  Zeff_e = Z - se
  
  result = a1_gam*Z2gamma/(Zeff*Zeff) + a1e_gam*Z2gamma_e/(Zeff_e*Zeff_e)
  return result


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def calc_sampson_p(om, Z, Te):
  from scipy.special import expn
  import calc_maxwell_rates_constants as const

  dE  = om[0]
  a   = om[1]
  z2s = om[2]
  c0  = om[3]
  cr  = om[4]
  cr1 = om[5]
  r   = int(om[6])
  s   = om[7]

  kT = const.KBOLTZ*Te
  y  = dE/kT

  a1 = a+1
  a1y = a1*y

  E1y = expn(1,y)
  Ery = expn(r,a1y)
  Er1y= expn(r+1,a1y)

  term = (cr*Ery/(a1**(r-1)) + cr1*Er1y/(a1**r))
  if (term > 0):
    Z2gamma = c0 + 1.333*z2s*E1y*numpy.exp(y) + y*numpy.exp(a1y)*term
  else:
    Z2gamma = c0 + 1.333*z2s*E1y*numpy.exp(y)
  
  Zeff = Z - s
    
  result = Z2gamma / (Zeff*Zeff)

  return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_sampson_h(om, Z, Te):
  from scipy.special import expn
  import calc_maxwell_rates_constants as const

  dE   = om[0]
  z2s  = om[1]
  a    = om[2]
  c0   = om[3]
  c1   = om[4]
  c2   = om[5]
  c_sw = om[6]
  Zeff2 = float(Z*Z) 

  kT  = const.KBOLTZ*Te
  y   = dE/kT
  a1  = a+1
  a1y = a1*y
  E1y = expn(1,y)
  Ery = expn(1,a1y)
  Er1y = expn(2,a1y)

  term = (c1*Ery + c2*Er1y/a1)
  if (term > 0):
    result = c0 + 1.333*z2s*E1y*numpy.exp(y) + y*numpy.exp(a1y)*term
  else:
    result = c0 + 1.333*z2s*E1y*numpy.exp(y)

  result *= 2*c_sw/Zeff2
  
  return result
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def calc_kato(coll_type, par, Z, Te):
  from scipy.special import expn
  import calc_maxwell_rates_constants as const

  # This fit comes from Kato & Nakazaki, 1989, adbatomic Data and Nuclear 
  #   Data Tables 42, 2, 313 
  
  dE = par[0]  # in keV
  A  = par[1]
  B  = par[2]
  C  = par[3]
  D  = par[4]
  E  = par[5]
  P  = par[6]
  Q  = par[7]
  X1 = par[8]
  
  kT = const.KBOLTZ*Te;
  y = dE/kT;

  if (coll_type==1):
    # Simple case (eq 6, in above reference)
    E1y = expn(1,y)
    term1 = A/y + C + (D/2)*(1-y)
    term2 = B - C*y + (D/2)*y*y + E/y
    result = y*(term1 + numpy.exp(y)*E1y*term2)
  elif (coll_type==2):
    # More complex case (eq 10-12 in above reference)
    E1Xy = expn(1,X1*y)
    term3 = numpy.exp(y*(1-X1))
    
    term1 = A/y + C/X1 + (D/2)*(1/(X1*X1) - y/X1) + (E/y)*numpy.log(X1)
    term2 = numpy.exp(y*X1)*E1Xy*(B - C*y + D*y*y/2 + E/y)

    ups_nr = y*term3*( term1 + term2 )
    ups_r = P*((1 + 1/y) - term3*(X1 + 1/y)) + Q*(1-term3)

    result = ups_nr + ups_r

  return result



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def interpolate_ionrec_rate(cidat,Te):
  import calc_maxwell_rates_constants as const
  from scipy import interpolate

  if ((Te > cidat['max_temp']) |\
      (Te < cidat['min_temp'])):
    print "Te outside of CI data range: Te = %e, Te_min=%e, Te_max=%e" %\
          (Te, cidat['min_temp'], cidat['max_temp'])
    return 0.0
  
  #real numbers.     
  if ((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
      (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC)):
    N_interp = cidat['par_type'] - const.INTERP_IONREC_RATE_COEFF
  elif ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
      (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC)):
    N_interp = cidat['par_type'] - const.INTERP_IONREC_RATE_OPEN
  elif ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
      (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC)):
    N_interp = cidat['par_type'] - const.INTERP_IONREC_RATE_OPEN

  elif ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
      (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC)):
    N_interp = cidat['par_type'] - const.INTERP_IONREC_RATE_OPEN

  else:
    print "interpolate_rate failed -- very odd."
    return 0.0
  
  # do interpolation

  tmpci = cidat['ionrec_par'][:N_interp]
  tmpci[tmpci<0] = 0.0

  try:
    tmp = numpy.exp(interpolate.interp1d(numpy.log(cidat['temperature'][:N_interp]), \
                             numpy.log(tmpci+1e-30), \
                             kind='cubic')(numpy.log(Te)))-1e-30
  except ValueError:
#    print "Raising a Value Error"
    if (numpy.log(Te) < numpy.log(cidat['temperature'][0])):
      if (numpy.log(Te) > numpy.log(cidat['temperature'][0])*0.99):
        tmp= tmpci[0]
#        zzz=raw_input('argh1')
      else:
        print cidat['temperature']
        print Te, Tetmp
        zzz=raw_input('ugh1')
        
    elif (numpy.log(Te) > numpy.log(cidat['temperature'][N_interp-1])):
      if (numpy.log(Te) < numpy.log(cidat['temperature'][N_interp-1])*1.01):
        tmp = tmpci[N_interp-1]
#        zzz=raw_input('argh2')
      else:
        print cidat['temperature']
        print Te,Tetmp
        zzz=raw_input('ugh2')
    else:
      print Te, cidat['temperature']
      print numpy.log(Te), numpy.log(cidat['temperature'][:N_interp])
#  print tmp
#  if numpy.isnan(tmp):
#    print '---'
#    print numpy.log(cidat['temperature'][:N_interp])
#    print numpy.log(cidat['ionrec_par'][:N_interp]+1e-30)
#    print cidat['ionrec_par'][:N_interp]
#    print numpy.log(cidat['ionrec_par'][:N_interp])
#    print Te, numpy.log(Te)
#    print '---'
  return tmp



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_ionrec_ci(cidat, Te):
  import calc_maxwell_rates_constants as const
  from scipy import interpolate

  if ((Te < cidat['min_temp'])|\
      (Te > cidat['max_temp'])):
    print " CI data is invalid at this temperature, Te=%e, Em_min = %e, Te_max=%e" %\
         (Te, cidat['min_temp'],cidat['max_temp'])
         
    ci = 0.0
    return ci

  
  # See Arnaud \& Rothenflug, 1985 A&ASS 60, 425 
  if (cidat['par_type'] == const.CI_YOUNGER):
    T_eV = 1.e3*const.KBOLTZ*Te
    x = cidat['ionrec_par'][0]/T_eV
    if (x <= 30.0):
      f1_val = f1_fcn(x)
      ci = (numpy.exp(-x)/x)*\
            ( cidat['ionrec_par'][1]*( 1 - x*f1_val ) +
              cidat['ionrec_par'][2]*( 1 + x - x*(x+2)*f1_val )+
              cidat['ionrec_par'][3]*f1_val + 
              cidat['ionrec_par'][4]*x*f2_fcn(x) )
    
    ci *= 6.69e-07/(T_eV*numpy.sqrt(T_eV))
    return ci
  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|
        ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|
        ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|
        ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
         
    ci = interpolate_ionrec_rate(cidat,Te)
    return ci
  else:
    print "Unknown CI type: %i"%(cidat['par_type'])
    return 0.0


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def calc_ionrec_rr(cidat, Te):
  import calc_maxwell_rates_constants as const
  from scipy import interpolate
  
  if ((Te < cidat['min_temp']) |(Te > cidat['max_temp'])):
    print " RR data is invalid at this temperature, Te=%e, Em_min = %e, Te_max=%e" %\
         (Te, cidat['min_temp'],cidat['max_temp'])
    rr = 0.0
    return rr

  
  
  # Shull & Van Steenberg 1982; 1982ApJS...48...95S  
  if (cidat['par_type'] == const.RR_SHULL):
    rr = cidat['ionrec_par'][0]* (Te/1.e4)**(-1*cidat['ionrec_par'][1])
    return rr


  # Verner & Ferland 1996; 1996ApJS..103..467V */
  elif (cidat['par_type'] == const.RR_VERNER):
    tt0=numpy.sqrt(Te/cidat['ionrec_par'][2])
    tt1=numpy.sqrt(Te/cidat['ionrec_par'][3])
    rr = cidat['ionrec_par'][0]/\
         ( tt0* (tt0+1)**( 1-cidat['ionrec_par'][1]) * \
                (1+tt1)**(1+cidat['ionrec_par'][1]) )
    
    return rr


  # Arnaud & Raymond 1992; 1992ApJ...398..394A  */
  elif (cidat['par_type'] == const.RR_ARNRAY):
    rr = cidat['ionrec_par'][0]*\
          (Te/1.e4)**(-1*(cidat['ionrec_par'][1]+\
                          cidat['ionrec_par'][2]*numpy.log10(Te/1.e4)))
    return rr

  # Now handle 4 different interpolation cases. */
  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
    rr = interpolate_ionrec_rate(cidat,Te)
    return rr  
  
  else:
    print "calc_ionrec_rate: RR Recombination type %i not recognized" %(cidat['par_type'])
    return 0.0

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def calc_ionrec_dr(cidat, Te):
  import calc_maxwell_rates_constants as const
  from scipy import interpolate

  if ((Te < cidat['min_temp']) |(Te > cidat['max_temp'])):
    print " DR data is invalid at this temperature, Te=%e, Em_min = %e, Te_max=%e" %\
         (Te, cidat['min_temp'],cidat['max_temp'])
    dr = 0.0
    return dr

  
  


  # Mazzotta
  if (cidat['par_type'] == const.DR_MAZZOTTA):
    T_eV = 1.e3*const.KBOLTZ*Te
    dr = (cidat['par_type'][0]/(T_eV**1.5))  * \
         (cidat['par_type'][5]*numpy.exp(-cidat['par_type'][1]/T_eV) +\
          cidat['par_type'][6]*numpy.exp(-cidat['par_type'][2]/T_eV) +\
          cidat['par_type'][7]*numpy.exp(-cidat['par_type'][3]/T_eV) +\
          cidat['par_type'][8]*numpy.exp(-cidat['par_type'][4]/T_eV))
    return dr

  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
    dr = interpolate_ionrec_rate(cidat,Te)
    return dr  
  
  else:
    print "calc_ionrec_rate: DR Recombination type %i not recognized" %(cidat['par_type'])
    return 0.0


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_ionrec_ea(cidat, Te):
  import calc_maxwell_rates_constants as const
  from scipy import interpolate


  if ((Te < cidat['min_temp']) |(Te > cidat['max_temp'])):
    print " EA data is invalid at this temperature, Te=%e, Te_min = %e, Te_max=%e" %\
         (Te, cidat['min_temp'],cidat['max_temp'])
    ea = 0.0
    return ea
    
      
  T_eV = 1.e3*const.KBOLTZ*Te

  # 1985A&AS...60..425A (Arnaud & Rothenflug) 
  if (cidat['par_type'] == const.EA_ARNROTH_LITHIUM):
    # Note: c[1] = b/Zeff^2 in Arnaud & Rothenflug 
    y = cidat['ionrec_par'][0]/T_eV
    if (y < 50.0):# otherwise, rate is effectively zero 
      yf1 = y*f1_fcn(y)
      G = 2.22*f1_fcn(y) + 0.67*(1 - yf1) + 0.49*yf1 + 1.2*y*(1 - yf1)
      ea = cidat['ionrec_par'][2] * cidat['ionrec_par'][1]  * \
           1.92e-07 * numpy.exp(-y) * G/numpy.sqrt(T_eV)
    return ea
  
  # 1985A&AS...60..425A (Arnaud & Rothenflug) s
  elif (cidat['par_type'] == const.EA_ARNROTH):
    I_ea = cidat['ionrec_par'][0] # in eV 
    a = cidat['ionrec_par'][1]
    b = cidat['ionrec_par'][2]
    y =  I_ea/T_eV
    if (y < 50.0):
      f1_val = f1_fcn(y)
      ea = 6.69e7 * a * (I_ea/numpy.sqrt(T_eV)) * numpy.exp(-y) *\
  (1.0 + b*f1_val - c[3]*y*f1_val -\
   c[4]*0.5*(y - y*y + y*y*y*f1_val))
    return ea


  # 1998A&AS..133..403M (Mazzotta etal 1998) 
  if (cidat['par_type'] == const.EA_MAZZOTTA_IRON):
    y=cidat['ionrec_par'][0]/T_eV
    if (y < 50.0):
      f1_val = f1_fcn(y)
      ea = (6.69e7/numpy.sqrt(T_eV))*numpy.exp(-y)*1.e-16*\
           (cidat['ionrec_par'][1]+\
            cidat['ionrec_par'][2]*(1-y*f1_val)+\
            cidat['ionrec_par'][3]*(1-y*(1-y*f1_val))+\
            cidat['ionrec_par'][4]*(1-0.5*(y - y*y + y*y*y*f1_val))+\
            cidat['ionrec_par'][5]*f1_val)

    return ea

  # Now handle 4 different interpolation cases. */
  if (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|\
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|\
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|\
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
    ea = interpolate_ionrec_rate(cidat, Te)
    return ea;  


  else:
    print "calc_ionrec_rate: EA type %i not recognized" %(cidat['par_type'])
    return 0.0



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_maxwell_rate(Te, colldata, index, lvdata, Te_unit='K', lvdatap1=False, ionpot = False):
  """
  Get the maxwellian rate for a transition from a file, typically for ionization,
  recombination or excitation.
  inputs:
    Te: electron temperature(s), in K by default
    colldata: the hdulist for the collision file (as returned by pyfits.open('file'))
    index: the line in the data to do the calculation for. Indexed from 0.
    lvdata: the hdulist for the energy level file (as returned by pyfits.open('file'))
    Te_unit: can be K or eV, K by default. 
  
  Returned values are in units of cm^3 s^-1  
  """

  Te_arr = numpy.array(Te)
  if Te_unit.lower()=='ev':
    Te_arr = Te_arr*11604.505
  elif Te_unit.lower() != 'k':
    print 'ERROR: Unknown Te_unit "%s": should be "K" or "keV"' % (Te_unit)
    return False
  else:
    pass
  exconly=False
  
  
  # get the data type 
  if colldata[1].header['HDUCLAS1'] == 'COLL_STR':
    dtype = 'EC'
  elif colldata[1].header['HDUCLAS1'] == 'IONREC':
    # get the type:
    dtype = colldata[1].data.field('TR_TYPE')[index]
    if not dtype in ['RR','DR','CI','EA', 'XI','XR']:
      print "Unknown transition type %s, should be one of RR, DR, CI, EA, XI, XR. "+\
            "Returning." % (dtype)
      return False
    # do things with DR
#    if colldata[1].data.field('PAR_TYPE')[index]>300:
#      d=interp_rate(Te_arr, colldata[1].data.field('PAR_TYPE')[index]-300,\
#                  colldata[1].data.field('TEMPERATURE')[index],\
#                  colldata[1].data.field('IONREC_PAR')[index])
#      return d
    
  else:
    print "ERROR: supplied data is not a collision strength (EC) or "+\
          "ionization/recombination (IR) set. Returning"
    return False
    

  # convert the data.
  
  if dtype=='EC':
    # go through all the different possibilities for collisional excitation data.
    ecdat = colldata[1].data[index]
    upind = ecdat['upper_lev']
    loind = ecdat['lower_lev']
    uplev = lvdata[1].data[upind]
    lolev = lvdata[1].data[loind]
    
    delta_E = uplev['energy']-lolev['energy']
    Z = lvdata[1].header['ELEMENT']
    degu = uplev['lev_deg']
    degl = lolev['lev_deg']

    if delta_E < 0:
      print "WARNING: delta_E < 0, upind = %i, loind = %i, eup =%e, elo=%e"%\
               (upind, loind, uplev['energy'],lolev['energy'])
      delta_E *= -1.0
      degl =uplev['lev_deg']
      degu = lolev['lev_deg']

    exc,dex = calc_maxwell_rates(ecdat['coeff_type'],\
                                 ecdat['min_temp'],\
                                 ecdat['max_temp'],\
                                 ecdat['temperature'],\
                                 ecdat['effcollstrpar'],\
                                 delta_E/1e3, Te_arr, Z, degl, degu)
    
    if numpy.isscalar(exc):
      if dex < 0:
        dex = 0.0
    else:
      dex[dex<0]=0.0
    if exconly:
      return exc
    else:
      return exc, dex

  elif dtype=='CI':
    cidat = colldata[1].data[index]
    
    if ((cidat['par_type']>900) & \
        (cidat['par_type']<=920)):

      
      upind = cidat['upper_lev']
      loind = cidat['lower_lev']
      uplev = lvdatap1[1].data[upind]
      lolev = lvdata[1].data[loind]
      
      delta_E = ionpot + uplev['energy']-lolev['energy']
      Z = lvdata[1].header['ELEMENT']
      degl = lolev['lev_deg']
      degu = uplev['lev_deg']
      
      ci, dex = calc_maxwell_rates(cidat['par_type'],\
                         cidat['min_temp'],\
                         cidat['max_temp'],\
                                 cidat['temperature'],\
                                 cidat['ionrec_par'],\
                                 delta_E/1e3, Te_arr, Z, degl, degu)
      
    
    else:
      ci = calc_ionrec_ci(cidat,Te)
      if (ci < 0.0):
        print "calc_ionrec_rate: CI(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
                   Te,ci)
      
    return ci

  elif dtype=='EA':
    cidat = colldata[1].data[index]
    ea = calc_ionrec_ea(cidat,Te)
    if (ea < 0.0):
      print "calc_ionrec_rate: EA(%10s -> %10s,T=%9.3e) = %8g"%\
                (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
                 adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
                 Te,ea)
      
    return ea

  elif dtype=='DR':
    cidat = colldata[1].data[index]
    dr = calc_ionrec_dr(cidat,Te)
    if (dr < 0.0):
      print "calc_ionrec_rate: DR(%10s -> %10s,T=%9.3e) = %8g"%\
                (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
                 adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
                 Te,dr)
      

    return dr


  elif dtype=='RR':
    cidat = colldata[1].data[index]
    rr = calc_ionrec_rr(cidat,Te)
    if (rr < 0.0):
      print "calc_ionrec_rate: rr(%10s -> %10s,T=%9.3e) = %8g"%\
                (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
                 adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
                 Te,rr)
      
    return rr



  elif dtype=='XR':
    cidat = colldata[1].data[index]
    xr = calc_ionrec_rr(cidat,Te)
    if (xr < 0.0):
      print "calc_ionrec_rate: xr(%10s -> %10s,T=%9.3e) = %8g"%\
                (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
                 adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
                 Te,xr)
      
    return xr


  elif dtype=='XI':
    cidat = colldata[1].data[index]
    
    if ((cidat['par_type']>900) & \
        (cidat['par_type']<=920)):

      # sanity check
      if (lvdata[1].header['ion_stat']+1 != cidat['ion_init']):
        print "ERROR: lvdata and cidat not for matching ions!"
      if (lvdatap1[1].header['ion_stat']+1 != cidat['ion_final']):
        print "ERROR: lvdatap1 and cidat not for matching ions!"
        
      upind = cidat['level_final']
      loind = cidat['level_init']
      uplev = lvdatap1[1].data[upind]
      lolev = lvdata[1].data[loind]
      
      try:
        delta_E = ionpot + uplev['energy']-lolev['energy']
      except IndexError:
        print upind
        print loind
        print ionpot
        print "uplev", uplev
        print uplev['energy']
        print lolev['energy']
        zzz=raw_input("HAD AN INDEXERROR")
      Z = lvdata[1].header['ELEMENT']
      degl = lolev['lev_deg']
      degu = uplev['lev_deg']

      
      xi, dex = calc_maxwell_rates(cidat['par_type'],\
                         cidat['min_temp'],\
                         cidat['max_temp'],\
                         cidat['temperature'],\
                         cidat['ionrec_par'],\
                         delta_E/1e3, Te_arr, Z, degl, degu)
      
#      xitmp = 2.1716e-8 * (13.6068/(Te_arr/11604.5))**0.5 * numpy.exp(delta_E/(Te_arr/11604.5)) /(numpy.pi*degl)
      
    else:
      xi = calc_ionrec_ci(cidat,Te)
      if (xi < 0.0):
        print "calc_ionrec_rate: CI(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
                   Te,xi)
        xi=0.0
      
    return xi


#
#  elif dtype=='XI':
#    cidat = colldata[1].data[index]
#    xi = calc_ionrec_ci(cidat,Te)
#    if (xi < 0.0):
#      print "calc_ionrec_rate: xi(%10s -> %10s,T=%9.3e, %i -> %i) = %8g"%\
#                (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']-1),\
#                 adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']-1),\
#                 Te,\
#                 cidat['level_init'],\
#                 cidat['level_final'],\
#                 xi)
#      xi=0.0
#      
#    return xi




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def sigma_hydrogenic(N,L, Z, Ein):
  """ Ein in keV """
  n = N
  l = L
  Z = Z
  coeff=1.09768e-19
  RYDBERG = 0.013605804
  E = numpy.array(Ein)

  chi = Z*Z*RYDBERG/(n*n)
  Eelec = (E - chi)/RYDBERG
  sigma = numpy.zeros(len(E), dtype=float)

  for iE in range(len(E)):
    if (Eelec[iE] <= 0): continue

    eta = Z/numpy.sqrt(Eelec[iE])
    rho = eta / n
    lp1 = l+1.0

    # Exponential term 
    expterm = numpy.exp(-4*eta*numpy.arctan(1./rho)) / (1-numpy.exp(-2*numpy.pi*eta))

    # Common Factorial
    commfact1 = 1.
    commfact2 = 1.
    if (n+l >= 2*l+2):
      for iFact in range(2*l+3, n+2):
        commfact1 *= iFact
    else:
      for iFact in range(n+l+1, 2*l+3):
        commfact2 *= iFact

    for iFact in range(2,2*l+3):
      commfact2 *= iFact
    for iFact in range(2,n-l):
      commfact2 *= iFact

  # Common Product term
    prodterm = 1.
    for iProd in range(1,l+1):
      prodterm *= (iProd*iProd + eta*eta)

  # Rho term
    rhoterm = (rho**(2*l+4.))/ ((1+rho*rho)**(2.*n))

  # G terms
    Gterm_p1 = (l+1-n)*G_hyd(l+1, n-l-1, eta, rho) +\
    ((l+1+n)/(1+rho*rho))*G_hyd(l+1, n-l, eta, rho)

    if (l!=0):
      Gterm_m1 = G_hyd(l, n-l-1, eta, rho) - \
        G_hyd(l, n-l+1, eta, rho)/((1.+rho*rho)**2)
    else:
      Gterm_m1 = 0.


  
    sigma_p1 = (coeff/E[iE])*((2)**(4*l+6.)/3.)*(lp1*lp1/(2*l+1.))*\
               (commfact1/commfact2)*(prodterm/(lp1*lp1+eta*eta))*\
               expterm*rhoterm*eta*eta*Gterm_p1*Gterm_p1

    sigma_m1 = 0.0

    if (l != 0):
      sigma_m1 = (l*l/(64.*lp1*lp1))*(2*l+1)*(2*l+2)*(2*l+1)*(2*l)*\
        (((1+rho*rho)*(1+rho*rho))/(eta*eta*rho*rho))*\
        ((lp1*lp1+eta*eta)/(l*l+eta*eta))*\
        (Gterm_m1*Gterm_m1/(Gterm_p1*Gterm_p1))*sigma_p1

    sigma[iE] = sigma_m1 + sigma_p1
  return sigma

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def G_hyd(l,m, eta, rho):


  result=0
  
  B_s = B_hyd(2*m, l, m, eta)
  result = B_s

  for s in range(2*m-1,0,-1):
    B_s = B_hyd(s, l, m, eta)
    result = result*rho + B_s

  return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def B_hyd(s, l, m, eta):
  MAX_BHYD = 10000
  d = numpy.zeros(MAX_BHYD, dtype=float)
  if ((s+1)*(s+1) > MAX_BHYD):
    print("B_hyd: Increase MAX_BHYD to at least %d",(s+1)*(s+1))


  d[0*(s+1) + 0] = 1.0
  for iS in range(1, s+1):
    for t in range(0,iS+1):
      term1 = 0.
      if ((iS > 0)&(t>0)): term1 = d[(iS-1)*(s+1)+t-1]
      term2 = 0.
      if ((iS > 1)&((iS-2) >= t)): term2 = d[(iS-2)*(s+1)+t]
      d[iS*(s+1) + t] = - (1./(iS*(iS+2.*l-1.)))*\
      ((4*(iS-1-m)*term1) + (2*m+2-iS)*(2*m+2*l+1-iS)*term2)

  result = d[s*(s+1)+s]
  for t in range(s-1,-1,-1):
    result = result*eta + d[s*(s+1) + t]

  return result


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def calc_rrc(levdat, plevdat, Te, ionpot, Z):
  #Te in K
  #ionpot in keV
  RRC_COEFF = 1.31e8
  ERG_KEV = 1.60219e-9
  
  lev_deg = levdat['lev_deg']
  if levdat['PHOT_TYPE'] == 0:
    parent_lev_deg = 1.0
  else:
    parent_lev_deg = plevdat['lev_deg']
  
  g_ratio = lev_deg*1.0/parent_lev_deg

  const  = g_ratio * RRC_COEFF/ERG_KEV

  rrc_ph_factor = const /(Te**1.5)
  
  
  if levdat['PHOT_TYPE'] == 0:
    sigparam  = {}
    sigparam['Zel'] = Z
    sigparam['nq'] = levdat['N_QUAN']
    sigparam['lq'] = levdat['L_QUAN']
    sigparam['phot_type'] = 0

    
#    sigma = sigma_hydrogenic(levdat['N_QUAN'],\
#                             levdat['L_QUAN'],\
#                             Z, energy)
  args = rrc_ph_factor, ionpot,Te, sigparam
  rr_lev_pop = integrate.quad(rrc_ph_value, ionpot+0.001, numpy.inf,\
                              args=args)
  
  return rr_lev_pop

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def rrc_ph_value (energy, rrc_ph_factor, ionpot, kT, sigparam):
  result = rrc_ph_factor * sigma_photoion(energy, sigparam) * energy**2 * \
           numpy.exp( -(energy-ionpot)/kT)
  return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def sigma_photoion(energy, sigparam):
  if sigparam['phot_type']==0:
    result = sigma_hydrogenic( sigparam['nq'], sigparam['lq'],sigparam['Zel'],\
                              energy)
    return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def calc_ionrec_ea(Z, Te, Ne, \
          par_type, ionrec_par, min_temp, max_temp, temperatures):
  # Te is in K          
            
  import calc_maxwell_rates_constants as const
  
  T_eV= Te*const.KBOLTZ*1e3
  ret = numpy.zeros(len(T_eV), dtype=float)
  if (par_type == const.EA_MAZZOTTA_IRON):
    for i in range(len(T_eV)):
      y=ionrec_par[0]/T_eV[i]
      if (y < 50.0):
         f1_val = f1_fcn(y)
         ea = (6.69e7/numpy.sqrt(T_eV[i]))*numpy.exp(-y)*1.e-16*\
             (ionrec_par[1]+\
              ionrec_par[2]*(1-y*f1_val)+\
              ionrec_par[3]*(1-y*(1-y*f1_val))+\
              ionrec_par[4]*(1-0.5*(y - y*y + y*y*y*f1_val))+\
              ionrec_par[5]*f1_val)
         ret[i] = ea
    
    return ret
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def prep_spline_atomdb(x, y, n):

  y2 = numpy.zeros(n,dtype=float)  
  u = numpy.zeros(n,dtype=float)  

  
  y2[0]=0.0
  u[0]=0.0

  for i in range(1,n-1):
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1])
    p=sig*y2[i-1]+2.0
    y2[i]=(sig-1.0)/p
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
  

  qn=0.0
  un=0.0
  
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0)

  for k in range(n-2,0,-1):
    y2[k]=y2[k]*y2[k+1]+u[k]

  
  return y2

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_spline_atomdb(xa, ya, y2a, n, x):


  klo=0
  khi=n-1

  while (khi-klo > 1):
    k=(khi+klo) >> 1
    if (xa[k] > x):
      khi=k
    else:
      klo=k

  h=xa[khi]-xa[klo]
  if (h == 0.0):
    print("calc_spline: Bad x array input\n")


  a=(xa[khi]-x)/h
  b=(x-xa[klo])/h
  result=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0

  return result

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def interpol_huntd(x, y, z):
   #x=xin
   #y=yin
   #z=xout
  inc = x[-1] > x[0]
  if (( inc & ((z < x[0]) | (z > x[-1]))) |\
      (-inc & ((z > x[0]) | (z < x[-1])))):
    print "interpol_huntd: Asking for %e, min is %e, max is %e"%\
          (z,x[0],x[-1])
    print "interpol_huntd: Cannot extrapolate"
    return numpy.nan
  
  jl = 0
  ju = len(x)
  while (ju - jl > 1):
    
    jm = (ju + jl) / 2
    if ((z > x[jm]) == inc):
      jl = jm
    else:
      ju = jm
    
  
  #/* ------ Z is sandwiched between JL and JU ------ */
  if (((x[jl] > 0.) & (x[ju] > 0.)) & ((y[jl] > 0.) & (y[ju] > 0.))):
    grad = (numpy.log10(y[ju]) - numpy.log10(y[jl])) /\
      (numpy.log10(x[ju]) - numpy.log10(x[jl]))
    df = grad * (numpy.log10(z) - numpy.log10(x[jl]))
    d1 = numpy.log10(y[jl]) + df
    f = 10.0**d1
  else:
    f = y[jl]+(y[ju]-y[jl])*((z-x[jl])/(x[ju]-x[jl]))
  
  return f



