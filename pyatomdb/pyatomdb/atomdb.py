"""
The atomdb module contains several routines for interfacing with the AtomDB
database to extract useful physical quantities, line lists, write new fits
files and more. It is currently a dump of everything I've done with AtomDB.
This should all be considered unstable and possibly susceptible to being
wrong. It will be fixed, including moving many routines out of this library,
as time goes on.

Version 0.1 - initial release
Adam Foster July 17th 2015


Version 0.2 - added PI reading routines and get_data online enhancements.
Adam Foster August 17th 2015

Version 0.3 - added RRC generation routines
Adam Foster August 28th 2015
"""

import os, datetime, numpy, re, time
import util, atomic, spectrum, const, urllib2
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
  for i in range(0, len(d['Z'])):

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
        a.write('%2i %2i %2i ' % (y, d['Z'][i], d['z1'][i])+\
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
  Reads AtomDB filemap data into dictionary
  
  Parameters
  ----------
  filemap : str
    Name of filemap file to read. If zero length, use "$ATOMDB/filemap"
  atomdbroot : str
    Replace any $ATOMDB in the file names with this. If not provided,
    use "ATOMDB" environment variable instead
    
  Returns
  -------
    dict
      dictionary listing all the filenames for each different datatype.
      Files separated into "misctypes" for abundance and bremsstrahlung files,
      and 'ionfiles', which are then split by element, ion, and datatype, e.g.
      ret['ionfiles'][6][5]['lv'] contains the energy level (lv) data for 
      carbon V.
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

  filedtype = numpy.dtype({'names':['Z','z1','ec','lv','ir','pc','dr',\
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
    Z_tmp = int(splt[1])
    z1_tmp = int(splt[2])
    
    if Z_tmp < 1:
      # in this case, we have a "misc" datatype, not corresponding to a particular ion
      misc_type = int(splt[0])
      misc_file = fname
      ret['miscfiles'] = numpy.append(ret['miscfiles'], numpy.array((misc_type, misc_file),\
                                                                    dtype=miscdtype))
      
    else:
      j = numpy.where((ret['ionfiles']['Z']==Z_tmp) & \
                      (ret['ionfiles']['z1']==z1_tmp))[0]
      if len(j)==0:
        j = len(ret['ionfiles'])

        ret['ionfiles'] = numpy.append(ret['ionfiles'], numpy.zeros(1,\
                                                                    dtype=filedtype))
        ret['ionfiles']['Z'][j]=Z_tmp
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
#  ret['Z'] = numpy.array(Z)
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
  """
  Calculates Collisional Ionization Rates from Younger formula
  
  Parameters
  ----------
  Te : array(float)
    Temperatures in Kelvin
  c : the ionrec_par from the transition in the AtomDB IR file
  
  Returns
  -------
  array(float) 
    returns ionization rate in cm^3 s^-1
  """
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

  result=numpy.zeros(len(x), dtype=float)

  a=numpy.array( [-0.57721566,0.99999193,-0.24991055,
     0.05519968,-0.00976004,0.00107857])
  b=numpy.array( [8.5733287401,18.0590169730,8.6347608925,0.2677737343,
     9.5733223454,25.6329561486,21.0996530827,3.9584969228])

  if sum(x < 0.0 )>0 :
    errmess("f1_fcn","Negative values of x not allowed")

  i = numpy.where(x < 1.0)[0]
  if len(i) > 0:
    result[i]=(a[0]+a[1]*x[i]+a[2]*x[i]*x[i] + a[3]* x[i]**3 + a[4]*x[i]**4 +
      a[5]*x[i]**5-numpy.log(x[i]))*numpy.exp(x[i])
  i = numpy.where((1.0 < x) & (x<1e9))[0]
  if len(i) > 0:
    result[i]=(x[i]**4+b[0]*x[i]**3+b[1]*x[i]**2+b[2]*x[i]+b[3])/ \
  (x[i]**4+b[4]*x[i]**3+b[5]*x[i]**2+b[6]*x[i]+b[7])/x[i]
  i = numpy.where((x>=1e9))[0]
  result[i]=1.0/x[i]

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

 
def get_ionfrac(ionbalfile, Z, te, z1=-1):
  """
  Reads the ionization fraction of a given ion at a given Te from an ionbalfile
  Assumes ionization equilibrium
  
  Parameters
  ----------
  
  ionbalfile : str
    location of ionization balance file
  Z : int
    atomic number of element (e.g. 6 for carbon)
  te : float
    electron temperature (in K)
  z1 : int
    if provided, z+1 of ion (e.g. 5 for O V). If omitted, returns ionization 
    fraction for all ions of element
  
  Returns
  -------
  
  ionization fraction of ion or, if not specified, of all ions at Te
  
  """
  
  dat = get_ionbal(ionbalfile, Z)
  it = numpy.argmin(abs(dat['te']-te))
  if z1==-1:
    z1 = range(0,Z+1)
  else:
    z1=z1-1  
  frac = dat['frac'][it,z1]
  return frac
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

 
def get_ionbal(ionbalfile, element, ion=-1):

  # ionbalfile : string, location of ionization balance file
  # element    : int, Z of element (e.g. 6 for carbon)
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
    print 'ion should be in range 1 to Z, in this case '+repr(Z)
    print 'return whole element ionization balance in stead'
  if ion > element + 1:
    print 'ion should be in range 1 to Z, in this case '+repr(Z)
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

def get_filemap_file(ftype, Z, z1, fmapfile=False,atomdbroot=False, quiet=False, misc=False):
  """
  Find the correct file from the database for atomic data of type ftype
  for ion with nuclear charge Z and ioncharge+1 = z1
  
  Parameters
  ----------
  ftype : str
   *  'ir' = ionization & recombination data
   *  'lv' = energy levels
   *  'la' = wavelength and transition probabilities (lambda & a-values)
   *  'ec' = electron collision rates
   *  'pc' = proton collision rates
   *  'dr' = dielectronic recombination satellite line information
   *  'ai' = autoionization rate data
   *  'pi' = XSTAR photoionization data
   *  'em' = emission feature data (currently unused)
  Z : int
    Element atomic number (=6 for C+4)
  z1 : int
    Ion charge +1 (=5 for C+4)
  fmapfile : str
    Specific filemap to use. Otherwise defaults to atomdbroot+'/filemap'
  atomdbroot : str
    Location of ATOMDB database. Defaults to ATOMDB environment variable.
    all $ATOMDB in the filemap will be expanded to this value
  quiet : bool
    If true, suppress warnings about files not being present for certain ions
  misc : bool
    If requesting "misc" data, i.e. the Bremsstrahlung inputs, use this. This is
    for non ion-specific data.

  Returns
  -------
  str
    The filename for the relevant file, with all $ATOMDB expanded. If 
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
    i = numpy.where((fmap['Z']==Z)&(fmap['z1']==z1))[0]
    ret=''
    if len(i)==0:
      if not quiet :
        print "WARNING: there is no data for the ion "+\
               adbatomic.spectroscopic_name(Z,z1)
      ret=''
           
    if len(i)>1:
      print "ERROR: there are multiple entries for the ion "+\
             adbatomic.spectroscopic_name(Z,z1)
      ret=''
    
    if len(i)==1:
      i=i[0]
      ftypel = ftype.lower()
    
      if not ftypel in fmap.keys():
      
        print "Error: invalid file type: "+ftype
        ret = ''
      
      else:
        ret = fmap[ftypel][i]
        if len(ret)==0:
          if not quiet :
            print "WARNING: no data of type "+ftype+" exists for ion "+\
                adbatomic.spectroscopic_name(Z,z1)

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
def get_ionpot(Z, z1, settings=False, 
               datacache=False):
  """
  Get the ionization potential of an ion in eV
  
  Parameters
  ----------
  Z : int
    The atomic number of the element
  z1 : int
    The ion charge + 1 of the ion
  settings : dict
    See description in get_data
  datacache : dict
    Used for caching the data. See description in get_data
  
  Returns
  -------
  float
    The ionization potential of the ion in eV.
    
  """
#
# Version 0.1 Initial Release
# Adam Foster 25 Sep 2015
#
  irdat = get_data(Z,z1,'IR', settings=settings, datacache=datacache)
  ionpot = irdat[1].header['IONPOT']
  
  return ionpot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def get_level_details(level, Z=-1, z1=-1, filename='', \
                      filemap='', atomdbroot=''):
  """Function returns the details in the level file for the specified 
  level. LV file can be specified by filename, or by filemap, Z, z1"""
  
  if filename=='':
    # get the filename from the other variables
    if Z < 0:
      print "Error in get_level_details: must specify filename or "+\
            "Z, z1 and filemap"
      return -1

    if z1 < 0:
      print "Error in get_level_details: must specify filename or "+\
            "Z, z1 and filemap"
      return -1

    if filemap =='':
      print "Error in get_level_details: must specify filename or "+\
            "Z, z1 and filemap"
      return -1
    filename = get_filemap_file(filemap, 'LV', Z, z1,\
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
  for Z in element:
    elsymb = adbatomic.Ztoelsymb(Z)
    ret[Z]=10**(a[1].data.field(elsymb)[ind[0]])/1e12
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

 
def get_ion_lines(linefile, Z, z1, fullinfo=False):
  a = pyfits.open(linefile)
  kT = a[1].data.field('kT')
  nlines = numpy.zeros(len(kT), dtype=int)
  for ikT in range(len(kT)):
    iikT = ikT + 2

    j = numpy.where((a[iikT].data.field("element") == Z) &\
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

 
def get_line_emissivity( Z, z1, upind, loind, \
                         linefile="$ATOMDB/apec_line.fits",\
                         ion_drv=False, elem_drv=False):
                         
  """
  Get the emissivity of a line as fn of temperature from APEC line file
  
  Parameters
  ----------
  Z : int
    Atomic number of element of line
  z1 : int
    Ion charge +1 of ion
  upind : int
    Upper level of transition
  loind : int
    Lower level of transition
  linefile : str
    line emissivity file. defaults to $ATOMDB/apec_line.fits
  ion_drv : int
    if set, return only the contribution from driving ion ion_drv. This is
    useful for non-equilibrium plasma calculations, and requires an 
    nei_line file to be specified in linefile
  elem_drv : int
    same as ion_drv, but specified driving element. Currently this setting
    is pointless, as all transitions have the same driving element as element.

  Returns
  -------
  dict
    dictionary with the following data in it:
  ['kT'] : array(float)
      the electron temperatures, in keV
  ['dens'] : array(float)
      the electron densities, in cm^-3
  ['time'] : array(float)
      the time (for old-style NEI files only, typically all zeros in 
      current files)
  ['epsilon'] : array(float)
      the emissivity in ph cm^3 s^-1
  
  """
  #
  # Version 0.1 - Initial Release
  # Adam Foster 25 Sep 2015
  #
    
  a = pyfits.open(os.path.expandvars(linefile))
  kT = a[1].data.field('kT')
  dens = a[1].data.field('eDensity')
  time = a[1].data.field('time')
  
  epsilon = numpy.zeros(len(kT), dtype=float)
  for ikT in range(len(kT)):
    iikT = ikT + 2
    if ion_drv:
      j = numpy.where((a[iikT].data.field("element") == Z) &\
                      (a[iikT].data.field("ion") == z1) &\
                      (a[iikT].data.field("ion_drv") == ion_drv) &\
                      (a[iikT].data.field("UpperLev") == upind) &
                      (a[iikT].data.field("LowerLev") == loind))[0]
    else:
      j = numpy.where((a[iikT].data.field("element") == Z) &\
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


#  print "Coll_type=%i"%(coll_type)
#  zzz=raw_input("URGH")


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
def interpolate_ionrec_rate(cidat,Te, force_extrap=False):
  from scipy import interpolate
  ret = numpy.zeros(len(Te), dtype=float)
#  if ((Te > cidat['max_temp']) |\
#      (Te < cidat['min_temp'])):
#    print "Te outside of CI data range: Te = %e, Te_min=%e, Te_max=%e" %\
#          (Te, cidat['min_temp'], cidat['max_temp'])
#    return 0.0
  
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
  
#  try:
  tmp = numpy.exp(interpolate.interp1d(numpy.log(cidat['temperature'][:N_interp]), \
                             numpy.log(tmpci+1e-30), \
                             kind=1, bounds_error=False,\
                             fill_value=numpy.nan)(numpy.log(Te)))-1e-30
  
  
#  print Te, tmp, cidat['min_temp'], cidat['max_temp']

#  for i in range(len(Te)):
#    print Te[i],tmp[i]
#  print 'input'
#  for i in range(len(cidat['temperature'])):
#    print cidat['temperature'][i], cidat['ionrec_par'][i]
  
  # let's deal with the near misses
  nantest =  numpy.isnan(tmp)
  for i in range(len(nantest)):
    if nantest[i]:
      if ((Te[i] > 0.99*(cidat['min_temp']))&\
          (Te[i] <= (cidat['min_temp']))):
        tmp[i] = cidat['ionrec_par'][0]
    
      if ((Te[i] < 1.01*(cidat['max_temp']))&\
          (Te[i] >= (cidat['max_temp']))):
        tmp[i] = cidat['ionrec_par'][N_interp-1]
        
  if force_extrap:
  
    nantest =  numpy.isnan(tmp)
    for i in range(len(nantest)):
      if (Te[i] > cidat['max_temp']):
        tmp[i]= tmp[N_interp-1] * (Te[N_interp-1]/Te[i])**1.5

  # set other "nan" to zero
  tmp[numpy.isnan(tmp)] = 0.0
  return tmp



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_ionrec_ci(cidat, Te, extrap=False):
  from scipy import interpolate
  ci = numpy.zeros(len(Te), dtype=float)
  ci[(Te< cidat['min_temp']) | \
     (Te> cidat['max_temp'])]= numpy.nan

  ici = numpy.where((Te<= cidat['max_temp']) & \
                    (Te>= cidat['min_temp']))[0]
#  if ((Te < cidat['min_temp'])|\
#      (Te > cidat['max_temp'])):
#    print " CI data is invalid at this temperature, Te=%e, Em_min = %e, Te_max=%e" %\
#         (Te, cidat['min_temp'],cidat['max_temp'])
#         
#    ci = 0.0
#    return ci

  # See Arnaud \& Rothenflug, 1985 A&ASS 60, 425 
  
  if (cidat['par_type'] == const.CI_YOUNGER):
    T_eV = 1.e3*const.KBOLTZ*Te[ici]
    x = cidat['ionrec_par'][0][ici]/T_eV
    i = numpy.where(x <=30.0)[0]
    if len(i) > 0:
      f1_val = f1_fcn(x[i])
      ci[ici[i]] = (numpy.exp(-x[i])/x[i])*\
            ( cidat['ionrec_par'][1]*( 1 - x[i]*f1_val ) +
              cidat['ionrec_par'][2]*( 1 + x[i] - x[i]*(x[i]+2)*f1_val )+
              cidat['ionrec_par'][3]*f1_val + 
              cidat['ionrec_par'][4]*x[i]*f2_fcn(x[i]) )

    ci[ici] *= 6.69e-07/(T_eV*numpy.sqrt(T_eV))

  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|
        ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|
        ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|
        ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
         (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
         
    ci[ici] = interpolate_ionrec_rate(cidat,Te[ici])
    
  else:
    print "Unknown CI type: %i"%(cidat['par_type'])


  # now see if extrappolation is required:
  if extrap:
    #  get Te, coeff for the 2 lowest points.
    ilow = numpy.where(Te<cidat['min_temp'])[0]
    
    if len(ilow) > 0:
      # calculate the value at Te_min
      Te_min = numpy.array([cidat['min_temp']])
      cimin = calc_ionrec_ci(cidat, Te_min)
      # if log of this is < 46.0, this is a small number, just repeat this
      if numpy.log(cimin[0])<46.0:
        ci[ilow]=cimin
      
      # otherwise, calculate the value at a range of near-minimum temperatures,
      # and use to construct a good second derivative.
      else:
        tetmp = numpy.logspace(numpy.log10(Te_min), numpy.log10(Te_min)+1,4)
        citmp=calc_ionrec_ci(cidat, tetmp)
      
        tetmpl = numpy.log(tetmp)
        citmpl = numpy.log(citmp)
        
        dci = (citmpl[1:]-citmpl[:-1])/(tetmpl[1:]-tetmpl[:-1])
        ddci = (dci[1:]-dci[:-1])/(tetmpl[1:-1]-tetmpl[:-2])
        dddci = (ddci[1:]-ddci[:-1])/(tetmpl[2:-1]-tetmpl[1:-2])

        a=dddci[0]/6.0
        b=-3.0*a*tetmpl[0]
        c=dci[0]-(3*a*tetmpl[0]**2+2*b*tetmpl[0])
        d = citmp[0]-(a*tetmpl[0]**3+b*tetmpl[0]**2+c*tetmpl[0])
      
        cinew= a * numpy.log(Te[ilow])**3 +\
               b * numpy.log(Te[ilow])**2 +\
               c * numpy.log(Te[ilow]) +\
               d
        ci[ilow]=cinew
        
#        Te_min = numpy.array([cidat['min_temp']])
#        cimin = calc_ionrec_ci(cidat, Te_min)
#        ci[ilow]=cimin[0]*(Te_min[0]/Te[ilow])**-4.5

    
    # calculate the value at "max_temp"
    ihigh = numpy.where(Te>cidat['max_temp'])[0]
    if len(ihigh) > 0:
      Te_max = numpy.array([cidat['max_temp']])
      cimax = calc_ionrec_ci(cidat, Te_max)
      ci[ihigh]=cimax[0]*(Te_max[0]/Te[ihigh])**0.5



  return ci


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def calc_ionrec_rr(cidat, Te, extrap=False):
  from scipy import interpolate

#  if ((Te < cidat['min_temp']) |(Te > cidat['max_temp'])):
#    print " RR data is invalid at this temperature, Te=%e, Em_min = %e, Te_max=%e" %\
#         (Te, cidat['min_temp'],cidat['max_temp'])
#    rr = 0.0
#    return rr
  rr = numpy.zeros(len(Te), dtype=float)
  
  irr = numpy.where((Te >= cidat['min_temp']) & \
                    (Te <= cidat['max_temp']))[0]
                    
  
  # Shull & Van Steenberg 1982; 1982ApJS...48...95S  
  if (cidat['par_type'] == const.RR_SHULL):
    rr[irr] = cidat['ionrec_par'][0]* (Te[irr]/1.e4)**(-1*cidat['ionrec_par'][1])


  # Verner & Ferland 1996; 1996ApJS..103..467V */
  elif (cidat['par_type'] == const.RR_VERNER):
    tt0=numpy.sqrt(Te[irr]/cidat['ionrec_par'][2])
    tt1=numpy.sqrt(Te[irr]/cidat['ionrec_par'][3])
    rr[irr] = cidat['ionrec_par'][0]/\
         ( tt0* (tt0+1)**( 1-cidat['ionrec_par'][1]) * \
                (1+tt1)**(1+cidat['ionrec_par'][1]) )
  # Arnaud & Raymond 1992; 1992ApJ...398..394A  */
  elif (cidat['par_type'] == const.RR_ARNRAY):
    rr[irr] = cidat['ionrec_par'][0]*\
          (Te[irr]/1.e4)**(-1*(cidat['ionrec_par'][1]+\
                          cidat['ionrec_par'][2]*numpy.log10(Te[irr]/1.e4)))


  # Now handle 4 different interpolation cases. */
  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
    rr[irr] = interpolate_ionrec_rate(cidat,Te[irr])
  
  else:
    print "calc_ionrec_rate: RR Recombination type %i not recognized" %(cidat['par_type'])
  
  if extrap:
    ilow = numpy.where(Te<cidat['min_temp'])[0]

    if len(ilow) > 0:
      # calculate the value at Te_min
      Te_min = numpy.array([cidat['min_temp']])
      cimin = calc_ionrec_rr(cidat, Te_min)
      # if log of this is < 46.0, this is a small number, just repeat this
      if numpy.log(cimin[0])<46.0:
        rr[ilow]=cimin
      
      # otherwise, calculate the value at a range of near-minimum temperatures,
      # and use to construct a good second derivative.
      else:
        tetmp = numpy.logspace(numpy.log10(Te_min), numpy.log10(Te_min)+1,4)
        citmp=calc_ionrec_rr(cidat, tetmp)
      
        tetmpl = numpy.log(tetmp)
        citmpl = numpy.log(citmp)
        
        dci = (citmpl[1:]-citmpl[:-1])/(tetmpl[1:]-tetmpl[:-1])
        ddci = (dci[1:]-dci[:-1])/(tetmpl[1:-1]-tetmpl[:-2])
        dddci = (ddci[1:]-ddci[:-1])/(tetmpl[2:-1]-tetmpl[1:-2])

        a=dddci[0]/6.0
        b=-3.0*a*tetmpl[0]
        c=dci[0]-(3*a*tetmpl[0]**2+2*b*tetmpl[0])
        d = citmp[0]-(a*tetmpl[0]**3+b*tetmpl[0]**2+c*tetmpl[0])
      
        cinew= a * numpy.log(Te[ilow])**3 +\
               b * numpy.log(Te[ilow])**2 +\
               c * numpy.log(Te[ilow]) +\
               d
        rr[ilow]=cinew



#    if len(ilow) > 0:
#      citmp=numpy.log(cidat['ionrec_par'][:4])
#      tetmp=numpy.log(cidat['temperature'][:4])
#      if citmp[0]< 46.0:
#        rr[ilow]=cidat['ionrec_par'][0]
#      else:
#        dci = (citmp[1:]-citmp[:-1])/(tetmp[1:]-tetmp[:-1])
#        ddci = (dci[1:]-dci[:-1])/(tetmp[1:-1]-tetmp[:-2])
#        dddci = (ddci[1:]-ddci[:-1])/(tetmp[2:-1]-tetmp[1:-2])
#
#        a=dddci[0]/6.0
#        b=-3.0*a*tetmp[0]
#        c=dci[0]-(3*a*tetmp[0]**2+2*b*tetmp[0])
#        d = citmp[0]-(a*tetmp[0]**3+b*tetmp[0]**2+c*tetmp[0])
#      
#        cinew= a * numpy.log(Te[ilow])**3 +\
#               b * numpy.log(Te[ilow])**2 +\
#               c * numpy.log(Te[ilow]) +\
#               d
#        rr[ilow] = numpy.exp(cinew)
#

    
    # calculate the value at "max_temp"
    ihigh = numpy.where(Te>cidat['max_temp'])[0]

    if len(ihigh) > 0:
      Te_max = numpy.array([cidat['max_temp']])
      rrmax = calc_ionrec_rr(cidat, Te_max)
      rr[ihigh]=rrmax[0]*(Te_max[0]/Te[ihigh])**1.5
  return rr
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def calc_ionrec_dr(cidat, Te, extrap=False):
  from scipy import interpolate

#  if ((Te < cidat['min_temp']) |(Te > cidat['max_temp'])):
#    print " DR data is invalid at this temperature, Te=%e, Em_min = %e, Te_max=%e" %\
#         (Te, cidat['min_temp'],cidat['max_temp'])
#    dr = 0.0
#    return dr
  dr = numpy.zeros(len(Te))
  # set values outside range to NAN
  idr = numpy.where((Te >= cidat['min_temp']) & (Te <= cidat['max_temp']))[0]
  dr[:] = numpy.nan
  dr[idr] = 0.0

  
  # Mazzotta
  if (cidat['par_type'] == const.DR_MAZZOTTA):
    T_eV = 1.e3*const.KBOLTZ*Te[idr]
    dr[idr] = (cidat['par_type'][0]/(T_eV**1.5))  * \
         (cidat['par_type'][5]*numpy.exp(-cidat['par_type'][1]/T_eV) +\
          cidat['par_type'][6]*numpy.exp(-cidat['par_type'][2]/T_eV) +\
          cidat['par_type'][7]*numpy.exp(-cidat['par_type'][3]/T_eV) +\
          cidat['par_type'][8]*numpy.exp(-cidat['par_type'][4]/T_eV))

  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & 
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
    dr[idr] = interpolate_ionrec_rate(cidat,Te[idr])
    
  else:
    print "calc_ionrec_rate: DR Recombination type %i not recognized" %(cidat['par_type'])

  # now extrappolate if required
  if extrap:
    ilow = numpy.where(Te<cidat['min_temp'])[0]
    if len(ilow) > 0:
      # calculate the value at Te_min
      Te_min = numpy.array([cidat['min_temp']])
      cimin = calc_ionrec_dr(cidat, Te_min)
      # if log of this is < 46.0, this is a small number, just repeat this
      if numpy.log(cimin[0])<46.0:
        dr[ilow]=cimin
      
      # otherwise, calculate the value at a range of near-minimum temperatures,
      # and use to construct a good second derivative.
      else:
        tetmp = numpy.logspace(numpy.log10(Te_min), numpy.log10(Te_min)+1,4)
        citmp=calc_ionrec_dr(cidat, tetmp)
      
        tetmpl = numpy.log(tetmp)
        citmpl = numpy.log(citmp)
        
        dci = (citmpl[1:]-citmpl[:-1])/(tetmpl[1:]-tetmpl[:-1])
        ddci = (dci[1:]-dci[:-1])/(tetmpl[1:-1]-tetmpl[:-2])
        dddci = (ddci[1:]-ddci[:-1])/(tetmpl[2:-1]-tetmpl[1:-2])

        a=dddci[0]/6.0
        b=-3.0*a*tetmpl[0]
        c=dci[0]-(3*a*tetmpl[0]**2+2*b*tetmpl[0])
        d = citmp[0]-(a*tetmpl[0]**3+b*tetmpl[0]**2+c*tetmpl[0])
      
        cinew= a * numpy.log(Te[ilow])**3 +\
               b * numpy.log(Te[ilow])**2 +\
               c * numpy.log(Te[ilow]) +\
               d
        dr[ilow]=cinew
        
#        Te_min = numpy.array([cidat['min_temp']])
#        cimin = calc_ionrec_ci(cidat, Te_min)
#        ci[ilow]=cimin[0]*(Te_min[0]/Te[ilow])**-4.5

    
    # calculate the value at "max_temp"
    ihigh = numpy.where(Te>cidat['max_temp'])[0]
    if len(ihigh) > 0:
      Te_max = numpy.array([cidat['max_temp']])
      cimax = calc_ionrec_dr(cidat, Te_max)
      dr[ihigh]=cimax[0]*(Te_max[0]/Te[ihigh])**1.5

  return dr

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_ionrec_ea(cidat, Te, extrap=False):
  from scipy import interpolate
  ea = numpy.zeros(len(Te))
  # set values outside range to NAN
  iea = numpy.where((Te >= cidat['min_temp']) & (Te <= cidat['max_temp']))[0]
  ea[:] = numpy.nan
  ea[iea] = 0.0


#  if ((Te < cidat['min_temp']) |(Te > cidat['max_temp'])):
#    print " EA data is invalid at this temperature, Te=%e, Te_min = %e, Te_max=%e" %\
#         (Te, cidat['min_temp'],cidat['max_temp'])
#    ea = 0.0
#    return ea
    
      
  T_eV = 1.e3*const.KBOLTZ*Te
  # 1985A&AS...60..425A (Arnaud & Rothenflug) 
  if (cidat['par_type'] == const.EA_ARNROTH_LITHIUM):
    # Note: c[1] = b/Zeff^2 in Arnaud & Rothenflug 
    y = cidat['ionrec_par'][0]/T_eV[iea]
    i = numpy.where(y<50.0)[0]
    if len(i)>0:# otherwise, rate is effectively zero 
      yf1 = y[i]*f1_fcn(y[i])
      G = 2.22*f1_fcn(y[i]) + 0.67*(1 - yf1) + 0.49*yf1 + 1.2*y[i]*(1 - yf1)
      ea[i] = cidat['ionrec_par'][2] * cidat['ionrec_par'][1]  * \
           1.92e-07 * numpy.exp(-y[i]) * G/numpy.sqrt(T_eV[iea[i]])
  
  # 1985A&AS...60..425A (Arnaud & Rothenflug) s
  elif (cidat['par_type'] == const.EA_ARNROTH):
    I_ea = cidat['ionrec_par'][0] # in eV 
    a = cidat['ionrec_par'][1]
    b = cidat['ionrec_par'][2]
    y =  I_ea/T_eV[iea]
    i = numpy.where(y < 50.0)[0]
    if len(y) > 0:
      f1_val = f1_fcn(y[i])
      ea[iea[i]] = 6.69e7 * a * (I_ea/numpy.sqrt(T_eV[iea[i]])) * numpy.exp(-y[i]) *\
  (1.0 + b*f1_val - cidat['ionrec_par'][3]*y[i]*f1_val -\
   cidat['ionrec_par'][4]*0.5*(y[i] - y[i]*y[i] + y[i]*y[i]*y[i]*f1_val))



  # 1998A&AS..133..403M (Mazzotta etal 1998) 
  elif (cidat['par_type'] == const.EA_MAZZOTTA_IRON):
    y=cidat['ionrec_par'][0]/T_eV[iea]
    i = numpy.where(y < 50.0)[0]
    if len(i)>0:
      f1_val = f1_fcn(y[i])
      ea[iea[i]] = (6.69e7/numpy.sqrt(T_eV[iea[i]]))*numpy.exp(-y[i])*1.e-16*\
           (cidat['ionrec_par'][1]+\
            cidat['ionrec_par'][2]*(1-y[i]*f1_val)+\
            cidat['ionrec_par'][3]*(1-y[i]*(1-y[i]*f1_val))+\
            cidat['ionrec_par'][4]*(1-0.5*(y[i] - y[i]*y[i] + y[i]*y[i]*y[i]*f1_val))+\
            cidat['ionrec_par'][5]*f1_val)



  # Now handle 4 different interpolation cases. */
  elif (((cidat['par_type'] >= const.INTERP_IONREC_RATE_COEFF) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_COEFF + const.MAX_IONREC))|\
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_OPEN) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_OPEN + const.MAX_IONREC))|\
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MIN) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MIN + const.MAX_IONREC))|\
      ((cidat['par_type'] >= const.INTERP_IONREC_RATE_INC_MAX) & \
       (cidat['par_type'] <= const.INTERP_IONREC_RATE_INC_MAX + const.MAX_IONREC))):
    ea[iea] = interpolate_ionrec_rate(cidat, Te[iea])



  else:
    print "calc_ionrec_rate: EA type %i not recognized" %(cidat['par_type'])
  # now extrappolate if required
  
  if extrap:
    ilow = numpy.where(Te<cidat['min_temp'])[0]
    ilow = numpy.where(Te<cidat['min_temp'])[0]
    if len(ilow) > 0:
      # calculate the value at Te_min
      Te_min = numpy.array([cidat['min_temp']])
      cimin = calc_ionrec_ea(cidat, Te_min)
      # if log of this is < 46.0, this is a small number, just repeat this
      if numpy.log(cimin[0])<46.0:
        ea[ilow]=cimin
      
      # otherwise, calculate the value at a range of near-minimum temperatures,
      # and use to construct a good second derivative.
      else:
        tetmp = numpy.logspace(numpy.log10(Te_min), numpy.log10(Te_min)+1,4)
        citmp=calc_ionrec_ea(cidat, tetmp)
      
        tetmpl = numpy.log(tetmp)
        citmpl = numpy.log(citmp)
        
        dci = (citmpl[1:]-citmpl[:-1])/(tetmpl[1:]-tetmpl[:-1])
        ddci = (dci[1:]-dci[:-1])/(tetmpl[1:-1]-tetmpl[:-2])
        dddci = (ddci[1:]-ddci[:-1])/(tetmpl[2:-1]-tetmpl[1:-2])

        a=dddci[0]/6.0
        b=-3.0*a*tetmpl[0]
        c=dci[0]-(3*a*tetmpl[0]**2+2*b*tetmpl[0])
        d = citmp[0]-(a*tetmpl[0]**3+b*tetmpl[0]**2+c*tetmpl[0])
      
        cinew= a * numpy.log(Te[ilow])**3 +\
               b * numpy.log(Te[ilow])**2 +\
               c * numpy.log(Te[ilow]) +\
               d
        ea[ilow]=cinew
        
#        Te_min = numpy.array([cidat['min_temp']])
#        cimin = calc_ionrec_ci(cidat, Te_min)
#        ci[ilow]=cimin[0]*(Te_min[0]/Te[ilow])**-4.5

    
    # calculate the value at "max_temp"
    ihigh = numpy.where(Te>cidat['max_temp'])[0]
    if len(ihigh) > 0:
      Te_max = numpy.array([cidat['max_temp']])
      cimax = calc_ionrec_ea(cidat, Te_max)
      ea[ihigh]=cimax[0]*(Te_max[0]/Te[ihigh])**0.5
  return ea


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_ionrec_rate(Te_in, irdat_in, lvdat_in=False, Te_unit='K', \
                     lvdatp1_in=False, ionpot=False, separate=False,\
                     Z=-1, z1=-1, settings=False, datacache=False,\
                     extrap=False):
  """
  Get the ionization and recombination rates at temperture(s) Te from 
  ionization and recombination rate data file irdat.
  
  Parameters
  ----------
  Te_in : float or arr(float)
    electron temperature in K (default), eV, or keV
  irdat_in : HDUList
    ionization and recombination rate data
  lvdat_in : HDUList
    level data for ion with lower charge (i.e. ionizing ion or recombined ion)
  Te_unit : {'K' , 'keV' , 'eV'} 
    temperature unit
  lvdatp1_in : HDUList
    level data for the ion with higher charge (i.e ionized or recombining ion)
  ionpot : float
    ionization potential of ion (eV).
  separate : bool
    if set, return DR, RR, EA and CI rates seperately. (DR = dielectronic recombination, RR = radiative recombination, EA = excitaiton autoionization, CI = collisional ionization)
    Note that EA & CI are not stored separately in all cases, so 
    may return zeros for EA as the data is incorporated into CI
    rates.
  Z : int
    Element charge to get rates for (ignores "irdat_in")
  z1  : int
    Ion charge +1 to get rates for (ignores "irdat_in")
    e.g. Z=6,z1=4 for C IV (C 3+)
  settings : dict
    See description in read_data
  datacache : dict
    See description in read_data
  extrap : bool
    Extrappolate rates to Te ranges which are off the provided 
    scale
  
  Returns
  -------
  float, float:
    (ionization rate coeff., recombination rate coeff.) in cm^3 s^-1
    *unless* separate is set, in which case:
  float, float, float, float:
    (CI, EA, RR, DR rate coeffs) in cm^3 s^-1
    Note that these assume low density & to get the real rates you need to
    multiply by N_e N_ion. 
  """
#  Version 0.1 Initial Release
#  Adam Foster 31st August 2015
#
  
  # input checking
  
  # check if Te is iterable
  Te = Te_in * 1.0

  isiter=True
  try:
    _ = (e for e in Te)
  except TypeError:
    isiter=False
    Te = numpy.array([Te])

  # check units of Te
  if Te_unit.lower()=='k':
    pass
  elif Te_unit.lower()=='ev':
    Te /=KBOLTZ
    Te *=1000.0
  elif Te_unit.lower()=='kev':
    Te /=KBOLTZ
  else:
    print "ERROR: units should be k, eV or keV"
    return -1

   # check IR data
  # see the type of the IR data

  if (z1>=0) &( Z>=0):
    irdat = get_data(Z,z1, 'IR',settings=settings, datacache=datacache)
  elif isinstance(irdat_in, basestring):
    #string (assumed filename)
    irdat = pyfits.open(irdat_in)
  elif type(irdat_in)==pyfits.hdu.hdulist.HDUList:
    # already opened IR file
    irdat = irdat_in
  else:
    print "ERROR: unable to process irdat as supplied."
    return -1

  if (z1>=0) &( Z>=0):
    lvdat = get_data(Z,z1, 'LV',settings=settings, datacache=datacache)
  elif isinstance(lvdat_in, basestring):
    #string (assumed filename)
    lvdat = pyfits.open(lvdat_in)
  elif type(lvdat_in)==pyfits.hdu.hdulist.HDUList:
    # already opened IR file
    lvdat = lvdat_in
  else:
    print "ERROR: unable to process lvdat as supplied."
    return -1

  if (z1>=0) &( Z>=0):
    lvdatp1 = get_data(Z,z1+1, 'LV',settings=settings, datacache=datacache)
  elif isinstance(lvdatp1_in, basestring):
    #string (assumed filename)
    lvdatp1 = pyfits.open(lvdatp1_in)
  elif type(lvdatp1_in)==pyfits.hdu.hdulist.HDUList:
    # already opened IR file
    lvdatp1 = lvdatp1_in
  else:
    lvdatp1 = False

  # -- OK, so now inputs should be *largely* settled -- *
  
  ciret = numpy.zeros(len(Te), dtype=float)
  earet = numpy.zeros(len(Te), dtype=float)
  rrret = numpy.zeros(len(Te), dtype=float)
  drret = numpy.zeros(len(Te), dtype=float)
  
  ionpot = irdat[1].header['IONPOT']
  # Start with the CI data
  ici = numpy.where(irdat[1].data['TR_TYPE']=='CI')[0]

  for i in ici:
    tmp=get_maxwell_rate(Te, irdat, i, lvdat, Te_unit='K', \
                         lvdatap1=lvdatp1, ionpot = False,\
                         force_extrap=extrap)  
    ciret += tmp

  iea = numpy.where(irdat[1].data['TR_TYPE']=='EA')[0]
  for i in iea:
    tmp=get_maxwell_rate(Te, irdat, i, lvdat, Te_unit='K', \
                         lvdatap1=lvdatp1, ionpot = False,\
                         force_extrap=extrap)  
    earet += tmp

  irr = numpy.where(irdat[1].data['TR_TYPE']=='RR')[0]

  for i in irr:
    tmp=get_maxwell_rate(Te, irdat, i, lvdat, Te_unit='K', \
                         lvdatap1=lvdatp1, ionpot = False,\
                         force_extrap=extrap)  
    rrret += tmp

  idr = numpy.where(irdat[1].data['TR_TYPE']=='DR')[0]
  for i in idr:

    tmp=get_maxwell_rate(Te, irdat, i, lvdat, Te_unit='K', \
                         lvdatap1=lvdatp1, ionpot = False,\
                         force_extrap=extrap)  

    drret += tmp

  # convert back to scalars if required
  if not isiter:
    ciret = ciret[0]
    earet = earet[0]
    rrret = rrret[0]
    drret = drret[0]

  # return the requested data
  if separate:
    return ciret, earet, rrret, drret
  else:  
    return ciret+ earet, rrret+ drret


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_maxwell_rate(Te, colldata, index, lvdata, Te_unit='K', \
                     lvdatap1=False, ionpot = False, \
                     force_extrap=False, silent=True):
  """
  Get the maxwellian rate for a transition from a file, typically for ionization,
  recombination or excitation.
  
  Parameters
  ----------
  
  Te : float
    electron temperature(s), in K by default
  colldata : HDUList
    The collisional data of interest
  index : int
    The line in the data to do the calculation for. Indexed from 0.
  lvdata : HDUList
    the hdulist for the energy level file (as returned by pyfits.open('file'))
  lvdatap1 : HDUList
    The level data for the recombining or ionized data.
  Te_unit : {'K','eV','keV'}
    Units of temperature grid.

  Returns
  -------
  float or array(float)  
    Maxwellian rate coefficient, in units of cm^3 s^-1  
  """
  isiter=True
  try:
    _ = (e for e in Te)
  except TypeError:
    isiter=False
    Te_arr = numpy.array([Te])
  else:
    Te_arr = Te
    
#  Te_arr = numpy.array(Te)
  if Te_unit.lower()=='ev':
    Te_arr = Te_arr*11604.505
  elif Te_unit.lower() != 'kev':
    Te_arr = Te_arr*11604.505*1000.0
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
      ci = calc_ionrec_ci(cidat,Te_arr, extrap=force_extrap)
      if sum(numpy.isnan(ci))>0:
        if not silent:
          print "calc_ionrec_rate: CI(%10s -> %10s): Te out of range min->max=%e->%e:"%\
                    (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                     atomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                     cidat['min_temp'], cidat['max_temp']),\
                     Te_arr[numpy.isnan(ci)]
      if sum(ci[numpy.isfinite(ci)] < 0)>0:
        if not silent:
          s= "calc_ionrec_rate: CI(%10s -> %10s): negative CI found: =%e->%e:"%\
                    (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                     atomic.spectroscopic_name(cidat['element'],cidat['ion_final']))
          for i in  numpy.where(ci[numpy.isfinite(ci)] < 0)[0]:
            s += " %e:%e, " % (Te_arr[i],ci[i])
          print s
      
      
    return ci

  elif dtype=='EA':
    cidat = colldata[1].data[index]
    ea = calc_ionrec_ea(cidat,Te, extrap=force_extrap)
    if sum(numpy.isnan(ea))>0:
      if not silent:
        print "calc_ionrec_rate: EA(%10s -> %10s): Te out of range min->max=%e->%e:"%\
                 (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                 atomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                 cidat['min_temp'], cidat['max_temp']),\
                 Te_arr[numpy.isnan(ea)]
    if sum(ea[numpy.isfinite(ea)] < 0)>0:
      if not silent:
        s= "calc_ionrec_rate: EA(%10s -> %10s): negative EA found: =%e->%e:"%\
                (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                 atomic.spectroscopic_name(cidat['element'],cidat['ion_final']))
        for i in  numpy.where(ea[numpy.isfinite(ea)] < 0)[0]:
          s += " %e:%e, " % (Te_arr[i],ea[i])
        print s
    return ea
  elif dtype=='DR':
    cidat = colldata[1].data[index]
    dr = calc_ionrec_dr(cidat,Te, extrap=force_extrap)
    if sum(numpy.isnan(dr))>0:
      if not silent:

        print "calc_ionrec_rate: DR(%10s -> %10s): Te out of range min->max=%e->%e:"%\
                (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                 atomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                 cidat['min_temp'], cidat['max_temp']),\
                 Te_arr[numpy.isnan(dr)]
    if sum(dr[numpy.isfinite(dr)] < 0)>0:
      if not silent:

        s= "calc_ionrec_rate: DR(%10s -> %10s): negative DR found: =%e->%e:"%\
                  (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   atomic.spectroscopic_name(cidat['element'],cidat['ion_final']))
        for i in  numpy.where(dr[numpy.isfinite(dr)] < 0)[0]:
          s += " %e:%e, " % (Te_arr[i],dr[i])
        print s
    return dr


  elif dtype=='RR':
    cidat = colldata[1].data[index]
    rr = calc_ionrec_rr(cidat,Te, extrap=force_extrap)
    if sum(numpy.isnan(rr))>0:
      if not silent:
        print "calc_ionrec_rate: RR(%10s -> %10s): Te out of range min->max=%e->%e:"%\
                (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                 atomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                 cidat['min_temp'], cidat['max_temp']),\
                 Te_arr[numpy.isnan(r)]
    if sum(rr[numpy.isfinite(rr)] < 0)>0:
      if not silent:
        s= "calc_ionrec_rate: RR(%10s -> %10s): negative RR found: =%e->%e:"%\
                  (atomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   atomic.spectroscopic_name(cidat['element'],cidat['ion_final']))
        for i in  numpy.where(rr[numpy.isfinite(rr)] < 0)[0]:
          s += " %e:%e, " % (Te_arr[i],rr[i])
        print s
    return rr



  elif dtype=='XR':
    cidat = colldata[1].data[index]
    xr = calc_ionrec_rr(cidat,Te, extrap=force_extrap)
    if sum(numpy.isnan(xr))>0:
      if not silent:

        print "calc_ionrec_rate: xr(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
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
      xi = calc_ionrec_ci(cidat,Te, extrap=force_extrap)
      if (xi < 0.0):
        print "calc_ionrec_rate: CI(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                   Te,xi)
        xi=0.0
      
    return xi


#
#  elif dtype=='XI':
#    cidat = colldata[1].data[index]
#    xi = calc_ionrec_ci(cidat,Te)
#    if (xi < 0.0):
#      print "calc_ionrec_rate: xi(%10s -> %10s,T=%9.3e, %i -> %i) = %8g"%\
#                (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
#                 adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
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
  """
  Calculate the PI cross sections of type hydrogenic.

  Parameters
  ----------
  N : int
    n shell
  L : int
    l quantum number
  Z : int
    nuclear charge
  Ein : array(float)
    energy grid for PI cross sections (in keV)

  Returns
  -------
  array(float)
    Photoionization cross section (in cm^2)

  """
#
# Version 0.1 - initial release
# Adam Foster August 28th 2015
#
  
  n = N
  l = L
  z = Z
  coeff=1.09768e-19
  RYDBERG = 0.013605804
  E = numpy.array(Ein)

  chi = Z*Z*1.0*RYDBERG/(n*n)
  print chi, Z, n
  zzz=raw_input()
  Eelec = (E - chi)/RYDBERG
  sigma = numpy.zeros(len(E), dtype=float)

  iE = numpy.where(Eelec > 0)[0]
  if len(iE > 0):

    eta = Z/numpy.sqrt(Eelec[iE])
    rho = eta / n
    lp1 = l+1.0
    print eta
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
def calc_rrc(Z, z1, eedges, Te, lev, xstardat=False, \
             xstarlevfinal=1,  settings=False, datacache=False):

  """
  Calculate the radiative recombination continuum for a given ion
  
  Parameters
  ----------
  Z : int
    Atomic number
  z1 : int
    recombined ion charge
  eedges : array(float)
    the bin edges for the spectrum to be calculated on (keV)
  Te : float
    The electron temperature (keV)
  lev : int
    The level of the ion for the recombination to be calculated into
  xstardat : dict or HDUList
    The xstar PI data. This can be an already sorted dictionary, as returned by
    sort_xstar_data, or the raw results of opening the PI file
  xstarlevfinal : int
    If you need to identify the recombining level, you can do so here. Should normally
    be 1.
  settings : dict
    See description in read_data
  datacache : dict
    See description in read_data  
  """
 
  ### WORKING HERE $$$$
  
  # calculate the RRC on bins of energy eedges (in keV)
  # Te (K)
  # I_e ionization potential (keV)
  # lev is the level in question
  # xstardat is the xstardata, if required
  
  #OK, let's get the line data.
  
  rrc = numpy.zeros(len(eedges)-1)
  
  lvdat = get_data(Z,z1,'LV',settings=settings, datacache=datacache)
  ldat = lvdat[1].data[lev-1]
#rrc_ph_value(E, Z, z1, rrc_ph_factor, IonE, kT, levdat, \
#                 xstardata=False, xstarfinallev=False):
  #Find the I_e for each case, and the gratio (recombined/recombining)
  I_e=0.0
  if ldat['PHOT_TYPE']==const.HYDROGENIC:
    I_e = const.RYDBERG*Z**2*1.0/(ldat['n_quan']**2)
    gfin = 1.0 * ldat['lev_deg']
    ginit = 1.0
    if Z-z1>0:
      plvdat = get_data(Z, z1+1, 'LV', settings==settings, datacache=datacache)
      if plvdat:
        ginit = 1.0* plvdat[1].data['lev_deg'][0]
    gratio = gfin/ginit
  elif ldat['PHOT_TYPE']==const.CLARK:
    I_e = ldat['phot_par'][1]
    gfin = 1.0 * ldat['lev_deg']
    ginit = 1.0
    if Z-z1>0:
      plvdat = get_data(Z, z1+1, 'LV', settings==settings, datacache=datacache)
      if plvdat:
        ginit = 1.0* plvdat[1].data['lev_deg'][0]
    gratio = gfin/ginit
  elif ldat['PHOT_TYPE']==const.VERNER:
    I_e = ldat['phot_par'][0]
    gfin = 1.0 * ldat['lev_deg']
    ginit = 1.0
    if Z-z1>0:
      plvdat = get_data(Z, z1+1, 'LV', settings==settings, datacache=datacache)
      if plvdat:
        ginit = 1.0* plvdat[1].data['lev_deg'][0]
    gratio = gfin/ginit

  elif ldat['PHOT_TYPE']==const.XSTAR:
    # get the xstar data
    if not xstardat:
      xstardat = get_data(Z,z1,'PI', settings=settings, datacache=datacache)
    if not isinstance(xstardat,dict):
      xstardat = sort_pi_data(xstardat, int(ldat['PHOT_PAR'][0]),\
                               xstarlevfinal)
    gratio = xstardat['g_ratio']
    I_e = (get_ionpot(Z, z1) - ldat['energy'])/1000.0
    
  else:
    #no photoionization data
    rrc = numpy.zeros(len(eedges)-1,dtype=float)
    return rrc  
  
  rrc_ph_factor = (const.RRC_COEFF/const.ERG_KEV)*gratio/(Te**1.5)
  rrc_erg_factor = const.RRC_COEFF*gratio/(Te**1.5) 
    
  # TOTAL INTEGRAL
#  rr_lev_pop,err = integrate.quad(rrc_ph_value, I_e+1.e-4, numpy.inf,\
#                                    epsabs=const.TOT_ABSACC,\
#                                    epsrel=const.TOT_RELACC, \
#                            args=(Z, z1, rrc_ph_factor, I_e, Te, ldat,
#                                  xstardat, xstarlevfinal), limit=100000)

  # VALUE AT EACH BIN EDGE
  emission_edges = rrc_ph_value(eedges, Z, z1,  rrc_ph_factor, I_e, \
                                Te, ldat, xstardat, xstarlevfinal)

  edgebin = numpy.argmin(eedges>I_e)-1

  if edgebin > -1:
    # we have an edge to deal with
    rrc[edgebin] = integrate.quad(rrc_ph_value, I_e+1.e-4, eedges[edgebin+1],\
                                    epsabs=const.TOT_ABSACC,\
                                    epsrel=const.TOT_RELACC, \
                            args=(Z, z1, rrc_ph_factor, I_e, Te, ldat,
                                  xstardat, xstarlevfinal), limit=100000)
  rrc[edgebin+1:] = (emission_edges[edgebin+1:-1]+emission_edges[edgebin+2:])/2.0
  

  return rrc

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calc_rad_rec_cont(Z, z1, z1_drv, T, ebins, abund=1.0, ion_pop=1.0, \
                      settings=False, datacache=False):
  """
  Calculate the radiative recombination continuum for an ion at temperature T
  
  Parameters
  ----------
  Z : int
    nuclear charge
  z1 : int
    recombined ion charge+1
  z1_drv : int
    recombining ion charge+1
  T : float
    temperautre (keV)
  ebins : array(float)
    energy bins (in keV) on which to caclulate the spectrum
  abund : float
    elemental abundance, relative to hydrogen
  ion_pop : float
    the ion's population fraction of that element 
    (i.e. sum of all ion_pop for an element = 1)

  Returns
  -------
  array(float)
    RRC in photons cm^3 s^-1 bin^-1, in an array of length(ebins)-1
  """
#  
#  Version 0.1 Initial Release
#  Adam Foster 25 Sep 2015
# 
  initlevels = get_data(Z,z1_drv,'LV', settings=settings, datacache = datacache)
  if (initlevels):
    hasparent = True
  else:
    hasparent = False
  finallevels = get_data(Z,z1,'LV', settings=settings, datacache=datacache)
  rrc = numpy.zeros(len(ebins)-1, dtype=float)
  if not finallevels:
    return rrc, 0.0
  nlev= len(finallevels[1].data)
  tot_rec_rate = 0.0
  rrc = numpy.zeros(len(ebins)-1, dtype=float)
  LevelRecombRate = numpy.zeros(nlev, dtype=float)
  
  for iLev in range(nlev):
      
    finlev = finallevels[1].data[iLev]
    if finlev['phot_type'] >=0:
      print "iLev=%i, phot_type=%i"%(iLev, finlev['phot_type'])
    if finlev['phot_type']==const.NO_PHOT:
      continue
    elif finlev['phot_type']==const.HYDROGENIC:
      if (Z-z1==0):
        I_e = const.RYDBERG*Z**2*1.0/(finlev['n_quan']**2)
        lev_deg = finlev['lev_deg']
        parent_lev_deg = 1.0
        sig_type = const.HYDROGENIC
        sigma_coeff=numpy.zeros(2, dtype=float)
        sigma_coeff[0] = finlev['n_quan']
        sigma_coeff[1] = finlev['l_quan']

      else:
        continue
    elif finlev['phot_type']==const.CLARK:
      I_e = finlev['phot_par'][1]
      lev_deg = finlev['lev_deg']
      if ((hasparent) & (len(initlevels[1].data)>0)):
        parent_lev_deg = initlevels[1].data['lev_deg'][0]*1.0
      else:
        parent_lev_deg=1.0
      sig_type = const.CLARK
      sigma_coeff = finlev['phot_par']

    elif finlev['phot_type']==const.VERNER:
      I_e = finlev['phot_par'][0]
      lev_deg = finlev['lev_deg']
      if ((hasparent) & (len(initlevels[1].data)>0)):
        parent_lev_deg = initlevels[1].data['lev_deg'][0]*1.0
      else:
        parent_lev_deg=1.0
      sig_type = const.VERNER
      sigma_coeff = finlev['phot_par']

    elif finlev['phot_type']==const.XSTAR:
      I_e = (get_ionpot(Z, z1, settings=settings, \
                               datacache=datacache) -\
             finlev['energy'])/1000.0
      lev_deg = finlev['lev_deg']*1.0
      if (hasparent):
        parent_lev_deg = initlevels[1].data['lev_deg'][0]*1.0
      else:
        parent_lev_deg = 1.0

      #### FIXING HERE
      xstarlevinit = int(finlev['phot_par'][0])
      pidat = get_data(Z,z1,'PI', settings=settings, \
                               datacache=datacache)
      print pidat[1].data['lev_final']
      print pidat[1].data['lev_init']==xstarlevinit
      print pidat[1].data['lev_final'][pidat[1].data['lev_init']==xstarlevinit]
      # find the possible places to photoionize to:
      finlev_possible = util.unique(\
        pidat[1].data['lev_final'][pidat[1].data['lev_init']==xstarlevinit])

      finlev_possible = numpy.array([min(finlev_possible)])
      print 'HI!', finlev_possible
      for xstarlevfinal in finlev_possible:
        # get the numbers
        sigma_coeff = sort_pi_data(pidat, xstarlevinit, xstarlevfinal)
        print sigma_coeff
        if not sigma_coeff:
          print 'oops'
          continue
        g_ratio = sigma_coeff['g_ratio']
        
        tmprrc = calc_rrc(Z, z1, ebins, T, iLev+1, xstardat=sigma_coeff, \
              settings=settings, datacache=datacache)
        print 'tmprrc', tmprrc
        tmprrc *= abund*ion_pop
        rr_lev_rate = sum(tmprrc)
        LevelRecombRate[iLev] = rr_lev_rate
        tot_rec_rate += rr_lev_rate
        rrc += tmprrc
    else:
      print "Uh..."
      continue
    rr_lev_rate = 0.0
    g_ratio = lev_deg/parent_lev_deg

    if ((finlev['phot_type'] != const.NO_PHOT) &\
        (finlev['phot_type'] != const.XSTAR)):

      tmprrc = calc_rrc(Z, z1, ebins, T, lev, \
              settings=settings, datacache=datacache)
      tmprrc *= abund*ion_pop
      rr_lev_rate = sum(tmprrc)
        
      LevelRecombRate[iLev] += rr_lev_rate
      tot_rec_rate += rr_lev_rate
      rrc += tmprrc
#      for i in xrange(len(tmprrc)):
        #print "%e %e %e" %(ebins[i], ebins[i+1], tmprrc[i])
      #print "phot_type= %i"%(finlev['phot_type'])
      #print finlev
  return rrc, LevelRecombRate

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def read_filemap(filemap="$ATOMDB/filemap", atomdbroot="$ATOMDB"):
  """
  Reads the AtomDB filemap file in to memory. By default, tries to read
  $ATOMDB/filemap, replacing all instances of $ATOMDB in the filemap
  file with the value of the environment variable $ATOMDB
  
  Parameters
  ----------
  
  filemap:  str
    the filemap file to read
  atomdbroot: str
    location of files, if not $ATOMDB.
  """
#  
#  Version 0.1 - initial release
#  Adam Foster August 15th 2015
#  
  
  # parse the options here.
  
  fmapfile = os.path.expandvars(filemap)
  
  f = open(fmapfile,'r')

  Z=[]
  z1=[]
  eclist=[]
  lvlist=[]
  pclist=[]
  lalist=[]
  drlist=[]
  irlist=[]
  emlist=[]
  pilist=[]
  ailist=[]
  cilist=[]
  misc=[]
  misc_type=[]
  
  # temporarily store ATOMDB environment variable
  os.environ['OLDATOMDB'] = os.environ['ATOMDB']  
  # set ATOMDB to new value
  os.environ['ATOMDB'] = os.path.expandvars(atomdbroot)
  
  for i in f:
#    print i
    splt = i.split()
#    print splt
    fname = os.path.expandvars(splt[3])
    Z_tmp = int(splt[1])
    z1_tmp = int(splt[2])

    if Z_tmp < 1:
      misc.append(fname)
      misc_type.append(int(splt[0]))
    else:
      j = numpy.where((numpy.array(Z)==Z_tmp) & \
                      (numpy.array(z1)==z1_tmp))[0]
      if len(j)==0:
        Z.append(Z_tmp)
        z1.append(z1_tmp)
        eclist.append('')
        lvlist.append('')
        lalist.append('')
        pclist.append('')
        drlist.append('')
        irlist.append('')
        emlist.append('')
        pilist.append('')
        ailist.append('')
        cilist.append('')
        j = numpy.where((numpy.array(Z)==Z_tmp) & \
                        (numpy.array(z1)==z1_tmp))[0]

      if int(splt[0]) == 1:
        irlist[j] = fname
      if int(splt[0]) == 2:
        lvlist[j] = fname
      if int(splt[0]) == 3:
        lalist[j] = fname
      if int(splt[0]) == 4:
        eclist[j] = fname
      if int(splt[0]) == 5:
        pclist[j] = fname
      if int(splt[0]) == 6:
        drlist[j] = fname
      if int(splt[0]) == 7:
        emlist[j] = fname
      if int(splt[0]) == 8:
        pilist[j] = fname
      if int(splt[0]) == 9:
        ailist[j] = fname
      if int(splt[0]) == 10:
        cilist[j] = fname


  ret={}
  ret['Z'] = numpy.array(Z)
  ret['z1'] = numpy.array(z1)
  ret['ec'] = numpy.array(eclist, dtype='|S160')
  ret['lv'] = numpy.array(lvlist, dtype='|S160')
  ret['ir'] = numpy.array(irlist, dtype='|S160')
  ret['pc'] = numpy.array(pclist, dtype='|S160')
  ret['dr'] = numpy.array(drlist, dtype='|S160')
  ret['la'] = numpy.array(lalist, dtype='|S160')
  ret['em'] = numpy.array(emlist, dtype='|S160')
  ret['pi'] = numpy.array(pilist, dtype='|S160')
  ret['ai'] = numpy.array(ailist, dtype='|S160')
  ret['ci'] = numpy.array(cilist, dtype='|S160')
  ret['misc'] = numpy.array(misc, dtype='|S160')
  ret['misc_type'] = numpy.array(misc_type)


  # restore the ATOMDB variable
  os.environ['ATOMDB']=os.environ['OLDATOMDB']
  x=os.environ.pop('OLDATOMDB')
  
  return ret



#-------------------------------------------------------------------------------
#--#-------------------------------------------------------------------------------
#--#-------------------------------------------------------------------------------
#--
#--def get_filemap_file(ftype, Z, z1, filemapfile=False, atomdbroot=False,\
#--                     quiet=False, misc=False):
#--  """
#--  Gets the filename for the file of type ftype for the ion Z,z1 from the 
#--  given filemap.
#--  
#--  INPUTS
#--  ftype: string: type of data to read. Currently available:
#--                'IR' - ionization and recombination
#--                'LV' - energy levels
#--                'LA' - radiative transition data (lambda and A-values)
#--                'EC' - electron collision data
#--                'PC' - proton collision data
#--                'DR' - dielectronic recombination satellite line data
#--                'PI' - XSTAR photoionization data
#--                'AI' - autoionization data
#--  Z: int : nuclear charge
#--  z1: int : ion charge +1 e.g. 5 for C+4, a.k.a. C V
#--  
#--  KWARGS
#--  filemapfile: string: if a particular filemap should be read. Otherwise
#--                       defaults to $ATOMDB/filemap
#--  atomdbroot: string : string to replace $ATOMDB by. If not set, use enviroment
#--                       variable $ATOMDB instead.
#--  quiet : If not set, code will provide warnings where data types not found
#--  misc  : If data requested is not ion specific but generic data, e.g. 
#--          bremsstrahlung data sets, set this to true and ftype becomes integer
#--          10 =  Abundances
#--          11 =  Bremsstrahlung: Hummer
#--          12 =  Bremsstrahlung: Kellogg
#--          13 =  Bremsstrahlung: Relativistic
#--
#--  RETURNS
#--     string filename if file exists.
#--     zero length string ('') if file does not exist in filemap.
#--
#--  Version 0.1 - initial release
#--  Adam Foster August 15th 2015
#--  """
#--  fmap=read_filemap(fmapfile, atomdbroot=atomdbroot)
#--  
#--  
#--  if misc:
#--    # looking for type 10,11 or 12 data
#--    if ftype in [10,11,12,13]:
#--      i = numpy.where(fmap['misc_type']==ftype)[0]
#--      if len(i)==0:
#--        print "Error: file type: %i not found in filemap %s" %(ftype,fmapfile)
#--        ret = ''        
#--      else:
#--        ret = fmap['misc'][i[0]]
#--  else:
#--    i = numpy.where((fmap['Z']==Z)&(fmap['z1']==z1))[0]
#--    ret=''
#--    if len(i)==0:
#--      if not quiet :
#--        print "WARNING: there is no data for the ion "+\
#--               atomic.spectroscopic_name(Z,z1)
#--      ret=''
#--           
#--    if len(i)>1:
#--      print "ERROR: there are multiple entries for the ion "+\
#--             atomic.spectroscopic_name(Z,z1)
#--      ret=''
#--    
#--    if len(i)==1:
#--      i=i[0]
#--      ftypel = ftype.lower()
#--    
#--      if not ftypel in fmap.keys():
#--      
#--        print "Error: invalid file type: "+ftype
#--        ret = ''
#--      
#--      else:
#--        ret = fmap[ftypel][i]
#--        if len(ret)==0:
#--          if not quiet :
#--            print "WARNING: no data of type "+ftype+" exists for ion "+\
#--                atomic.spectroscopic_name(Z,z1)
#--
#--  return ret
#--
#--#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_data(Z, z1, ftype, datacache=False, \
             settings=False, indexzero=False, offline=False):
  """
  Read AtomDB data of type ftype for ion rmJ of element Z.
  
  If settings are set, the filemap can be overwritten (see below), otherwise
  $ATOMDB/filemap will be used to locate the file.
  If indexzero is set, all levels will have 1 subtracted from them (AtomDB
  indexes lines from 1, but python and C index arrays from 0, so this can be
  useful)

  Parameters
  ----------
  
  Z : int
    Element nuclear charge
  rmJ : int
    Ion charge +1 (e.g. 5 for C^{4+}, a.k.a. C V)
  ftype : string
    type of data to read. Currently available
       *           'IR' - ionization and recombination
       *           'LV' - energy levels
       *           'LA' - radiative transition data (lambda and A-values)
       *           'EC' - electron collision data
       *           'PC' - proton collision data
       *           'DR' - dielectronic recombination satellite line data
       *           'PI' - XSTAR photoionization data
       *           'AI' - autoionization data
  filemap : string
    The filemap to use, if you do not want to use the default one.
    
  settings : dict
    This will let you override some standard inputs for get_data:
    
    * settings['filemap']: the filemap to use if you do not want to use
      the default $ATOMDB/filemap
      
    * settings['atomdbroot']: If you have files in non-standard locations
      you can replace $ATOMDB with this value
    
  datacache : dict
    This variable will hold the results of the read in a dictionary.
    It will also be checked to see if the requested data has already been
    cached here before re-reading from the disk. If you have not yet read
    in any data but want to start caching, provide it as an empty dictionary
    i.e. mydatacache={}
    
    2 parts of the data ares stored here:
     
    * Settings['data'] will store a copy of the data
      you read in. This means that if your code ends up calling for the same
      file multiple times, rather than re-reading from the disk, it will just
      point to this data already in memory. To clear the read files, just reset
      the data dictionary (e.g. settings['data'] ={})
  
    * settings['datasums'] stores the datasum when read in. Can be used later 
      to check files are the same.

    Both data and datasums store the data in identical trees, e.g.:
    settings['data'][Z][z1][ftype] will have the data.
    
  indexzero: bool
    If True, subtract 1 from all level indexes as python indexes from 0, 
    while AtomDB indexes from 1.
    
  offline: bool
    If True, do not search online to download data files - just return
    as if data does not exist
  
  Returns
  -------
  HDUlist  
    the opened pyfits hdulist if succesful. False if file doesn't exist
  """
#  
#  
#  Version 0.1 - initial release
#  Adam Foster August 15th 2015

#  Version 0.2 - separated "settings" and "datacache"
#  Adam Foster September 24th 2015

  d = False
  didurl=False

  if datacache != False:
    # make sure that the relevant dictionaries are ready to receive the data
    if not 'data' in datacache.keys():
      datacache['data']={}
    if not Z in datacache['data'].keys():
      datacache['data'][Z]={}
    if not z1 in datacache['data'][Z].keys():
      datacache['data'][Z][z1]={}

    if not 'datasums' in datacache.keys():
      datacache['datasums']={}
    if not Z in datacache['datasums'].keys():
      datacache['datasums'][Z]={}
    if not z1 in datacache['datasums'][Z].keys():
      datacache['datasums'][Z][z1]={}


    if ftype.upper() in datacache['data'][Z][z1].keys():
      # this means we have the data cached, no need to fetch it
      pass
    else:
      # check for file location overrides
      fmapfile = False
      atomdbroot=False
      if settings:
        if settings['filemap']:
          fmapfile = settings['filemap']
        if settings['atomdbroot']:
          atomdbroot = settings['atomdbroot']
          
      fname = atomdb.get_filemap_file(ftype, Z, z1, fmapfile=fmapfile,\
                             atomdbroot=atomdbroot, quiet=True)

      if fname=='':
        # no data exists
        datacache['data'][Z][z1][ftype.upper()] = False

      else:
        # Try and open the files in the following order:
        # (1) filename
        # (2) filename+'.gz'
        # (3) atomdburl/filename+'.gz'
        
        try:
          d = pyfits.open(fname)
        except IOError:
          try:
            d = pyfits.open(fname+'.gz')
          except IOError:
            if offline:
              d = False
            else:
              url = re.sub(atomdbroot,\
                           'ftp://sao-ftp.harvard.edu/AtomDB',fname)+'.gz'
              try:
                d = pyfits.open(url)
                didurl=True
                util.record_upload(re.sub(atomdbroot,'',fname))
              except urllib2.URLError:
                print "Error trying to open file %s. Not found locally or on"%(fname)+\
                    " server." 
                d=False

  else:
    fmapfile = False
    atomdbroot=os.environ['ATOMDB']
    if settings:
      if settings['filemap']:
        fmapfile = settings['filemap']
      if settings['atomdbroot']:
        atomdbroot = settings['atomdbroot']

    fname = get_filemap_file(ftype, Z, z1, \
                             quiet=True, fmapfile=fmapfile,\
                             atomdbroot=atomdbroot)

    if fname=='':
      pass
    else:
        # Try and open the files in the following order:
        # (1) filename
        # (2) filename+'.gz'
        # (3) atomdburl/filename+'.gz'
      try:
        d = pyfits.open(fname)
      except IOError:
        try:
          d = pyfits.open(fname+'.gz')
        except IOError:
          if offline:
            d = False
          else:
            url = re.sub(atomdbroot,\
                         'ftp://sao-ftp.harvard.edu/AtomDB',fname)+'.gz'
            try:
              d = pyfits.open(url)
              didurl=True
              util.record_upload(re.sub(atomdbroot,'',fname))
            except urllib2.URLError:
              print "Error trying to open file %s. Not found locally or on"%(fname)+\
                    " server." 
              d=False
            

  if didurl:
    # cache file locally
    
    # make directory if required
    util.mkdir_p(fname.rsplit('/',1)[0])
    # save file
    d.writeto(fname, output_verify='warn')
    print "wrote file locally to %s"%(fname)

  if d:  
    if indexzero:
      # python indexes from zero. Easiest to just subtract here.
      if ftype.upper()=='LA':
        d[1].data.field('lower_lev')[:] -= 1
        d[1].data.field('upper_lev')[:] -= 1
      elif ftype.upper()=='EC':
        d[1].data.field('lower_lev')[:] -= 1
        d[1].data.field('upper_lev')[:] -= 1
      elif ftype.upper()=='PC':
        d[1].data.field('lower_lev')[:] -= 1
        d[1].data.field('upper_lev')[:] -= 1
      elif ftype.upper()=='DR':
        d[1].data.field('lower_lev')[:] -= 1
        d[1].data.field('upper_lev')[:] -= 1
      elif ftype.upper()=='AI':
        d[1].data.field('level_init')[:] -= 1
        d[1].data.field('level_final')[:] -= 1
      elif ftype.upper()=='IR':
        d[1].data.field('level_init')[:] -= 1
        d[1].data.field('level_final')[:] -= 1
      elif ftype.upper()=='LV':
        pass
      elif ftype.upper()=='PI':
        d[1].data.field('lev_init')[:] -= 1
        d[1].data.field('lev_final')[:] -= 1
      else:
        print "Unknown filetype: %s"%(ftype)

    # rename columns in older versions of EC & PC files
    if ftype.upper() in ['PC','EC']:
      if d[1].header['HDUVERS1']=='1.0.0':
        d[1].columns.change_name('COEFF_OM','EFFCOLLSTRPAR')
        d[1].columns[d[1].data.names.index('COEFF_OM')].name='EFFCOLLSTRPAR'
    if datacache:
      datacache['data'][Z][z1][ftype.upper()] = d
      datacache['datasums'][Z][z1][ftype.upper()] = d[1].header['DATASUM']
   
      return datacache['data'][Z][z1][ftype.upper()]
    else:
      return d
  else:
    return False
    
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def sort_pi_data(pidat, lev_init, lev_final):
  """
  Given the pidat (returned by opening the PI data file, i.e.
  pyfits.open('XX_YY_PI.fits'), and the initial and final levels, return
  the PI cross section data.

  Parameters
  ----------
  pidat : hdulist
    The photoionization data for the ion
  lev_init : int
    The initial level
  lev_final : int
    The final level
  
  Returns
  --------
  dict:
    which contains the following information:
    pi['ion_init'] - the initial ion charge +1
    pi['lev_init'] - the initial level
    pi['ion_final'] - the final ion charge+1 (should be ion_init+1)
    pi['lev_final'] - the final level
    pi['pi_type'] - the type. (best to ignore)
    pi['g_ratio'] - the ratio of the statistical weight of the intitial and final levels
    pi['energy'] - the array of energies (keV)
    pi['pi_param'] - the array of pi cross sections in Mbarn.
"""
#
#  Version 0.1 - initial release
#  Adam Foster August 15th 2015
#
  print 'lev_init=', lev_init, ', lev_final=', lev_final
  i = numpy.where((pidat[1].data['lev_init'] ==lev_init) &\
                  (pidat[1].data['lev_final'] ==lev_final))[0]
  print 'i=', i
  energy = numpy.zeros(0,dtype=float)
  pi_param = numpy.zeros(0,dtype=float)
  if len(i)==0:
    return False
  for ii in i:
    energy = numpy.append(energy, pidat[1].data['energy'][ii])
    pi_param = numpy.append(pi_param, pidat[1].data['pi_param'][ii])
  n_expected = sum(numpy.array(util.unique(pidat[1].data['pi_type'][i]))%10000)
  
  n_found = sum(energy>0)
  if n_found != n_expected:
    print "WARNING: we do not have the same length of expected and found parameters"
    print "Expected: %i, found %i"%(n_expected, n_found)
  pi_param = pi_param[energy>0]
  energy = energy[energy>0]/1000.0 # convert to keV
  
  i  = numpy.argsort(energy)
  pi_param = pi_param[i]
  energy = energy[i]
  pi = {}
  pi['ion_init'] = pidat[1].data['ion_init'][ii]
  pi['lev_init'] = pidat[1].data['lev_init'][ii]
  pi['ion_final'] = pidat[1].data['ion_final'][ii]
  pi['lev_final'] = pidat[1].data['lev_final'][ii]
  pi['pi_type'] = pidat[1].data['pi_type'][ii]
  pi['g_ratio'] = pidat[1].data['g_ratio'][ii]
  pi['energy'] = energy
  pi['pi_param'] = pi_param
  

  return pi
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def rrc_ph_value(E, Z, z1, rrc_ph_factor, IonE, kT, levdat, \
                 xstardata=False, xstarfinallev=False):
  isiter=True
  try:
    _ = (e for e in E)
  except TypeError:
    isiter=False
    E = numpy.array([E])

  igood = numpy.where((-(E - IonE)/kT)>const.MIN_RRC_EXPONENT)[0]

  result = numpy.zeros(len(E))
  result[igood] = rrc_ph_factor*sigma_photoion(E[igood], 
                                               Z,\
                                               z1,\
                                               levdat['PHOT_TYPE'],\
                                               levdat['PHOT_PAR'],\
                                               xstardata=xstardata,\
                                               xstarfinallev=xstarfinallev)*\
                  E[igood]**2*numpy.exp(-(E[igood] - IonE)/kT)
  if not isiter:
    result=result[0]
  return result
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


def sigma_photoion(E, Z, z1, pi_type, pi_coeffts, xstardata=False, xstarfinallev=1):
  """
  Returns the photoionization cross section at E, given an input of sig_coeffts.
   
  Parameters
  ----------
  E: float or array of floats
    Energy/ies to find PI cross section at (keV)
  Z: int
    Atomic number of element (i.e. 8 for Oxygen)
  pi_type : int
    the "PI_TYPE" from the energy level file for this level, can be:
   
     -1. No PI data
     0.  Hydrogenic
     1.  Clark
     2.  Verner
     3.  XSTAR
     
  pi_coeffts : array(float)
    the "PI_PARAM" array for this level from the LV file
  xstardata : dict, str or HDUList
    if the data is XSTAR data (pi_type=3), supply the xstardata. This can be 
    a dictionary with 2 arrays, one "Energy", one "sigma", the file name, or
    the entire PI file (already loaded)::
    
      # load level data
      lvdata = atomdb.get_data(26, 24, 'LV', settings)
    
      # load XSTAR PI data if it exists
      pidata = atomdb.get_data(26, 24, 'PI', settings)
    
      # get pi xsection at energy E for the ground state to ground state
      sigma_photoion(E,
                     lvdata[1].data['pi_type'][0],
                     lvdata[1].data['pi_param'][0],
                     xstardata=pidata,
                     xstarfinallev=1)
  xstarfinallev: the level to ionize in to. Defaults to 1.
                    
  Returns
  -------
  array(float)
    pi cross section in cm^2 at energy E. 
     
  """ 
#   
#  Version 0.1 - initial release
#  Adam Foster August 15th 2015
#
#  Version 0.2 - 
#  now accepts vector energy input.
#  Will attempt to download the PI cross section files if not supplied.
#  Adam Foster August 28th 2015

  # determine whether input is scalar or vector
  isvec = False
  try:
    _ = (e for e in E)
  except TypeError:
    isvec = False
  Evec = numpy.array(E)
  result = numpy.zeros(len(Evec), dtype=float)
  
  # set up the sigma coefficients
  if pi_type==const.NO_PHOT:
    pass
  elif pi_type==const.VERNER:
    E_th = pi_coeffts[0]
    E_0 = pi_coeffts[1]
    sigma0 = pi_coeffts[2]
    ya = pi_coeffts[3]
    P = pi_coeffts[4]
    yw = pi_coeffts[5]
    l1 = pi_coeffts[6]
    Q = 5.5 + l1 - 0.5*P
    sig_coeffts={}
    sig_coeffts['E_th']   = E_th
    sig_coeffts['E_0']    = E_0
    sig_coeffts['sigma0'] = sigma0
    sig_coeffts['ya']     = ya
    sig_coeffts['P']      = P
    sig_coeffts['yw']     = yw
    sig_coeffts['l1']     = l1
    sig_coeffts['Q']      = Q
    
    
  elif pi_type==const.HYDROGENIC:
    nq = int(round(pi_coeffts[0]))
    lq = int(round(pi_coeffts[1]))
    Zel = Z
    sig_coeffts={}
    sig_coeffts['nq']=nq
    sig_coeffts['lq']=lq
    sig_coeffts['Zel']=Zel
    
  elif pi_type==const.CLARK:
    Zel = Z
    frac = pi_coeffts[0]
    IonE = pi_coeffts[1]
    ip = pi_coeffts[2]
    n1 = pi_coeffts[3]
    b1 = pi_coeffts[4]
    d1 = pi_coeffts[5]
    b2 = pi_coeffts[6]
    d2 = pi_coeffts[7]
    c = numpy.array(pi_coeffts[8:12])

    sig_coeffts={}
    sig_coeffts['Zel']  = Zel
    sig_coeffts['frac'] = frac
    sig_coeffts['ip']   = ip
    sig_coeffts['n1']   = n1
    sig_coeffts['b1']   = b1
    sig_coeffts['d1']   = d1
    sig_coeffts['b2']   = b2
    sig_coeffts['d2']   = d2
    sig_coeffts['c']    = c
    sig_coeffts['I_e'] = IonE

  elif pi_type==const.XSTAR:
    # now look into type of xstardata
    
    if not xstardata:
      # get the data
      pidat = get_data(Z, z1, 'PI')
      if not pidat:
        print "ERROR: cannot find photoionization data requested"
        return False
    # just filename

    if isinstance(xstardata,dict):
      sig_coeffts = xstardata

    else:
      if isinstance(xstardata, basestring):
      #yay.
        pidat = pyfits.open(xstardata)
        initlevel = int(pi_coeffts[0])

      sigma_coeff = sort_pi_data(pidat, initlevel, xstarfinallev)
      sig_coeffts = sigma_coeff


  # now calculate the sigma.

  if pi_type==const.NO_PHOT:
    result = 0.0

  elif pi_type == const.VERNER:
    iE = numpy.where(E > sig_coeffts['E_th'])[0]
    if len(iE) > 0:
      y = E[iE]/sig_coeffts['E_0']
      F_y = ((y-1)**2. + sig_coeffts['yw']**2.) * y**(-sig_coeffts['Q']) *\
            (1.+numpy.sqrt(y/sig_coeffts['ya']))**(-sig_coeffts['P'])
      result[iE] = sig_coeffts['sigma0'] * F_y * 1e-18 # conversion to cm2 

  elif pi_type == const.HYDROGENIC:
    if (sig_coeffts['nq']<6) & (sig_coeffts['nq']!=0):
      result = sigma_hydrogenic(sig_coeffts['Zel'],\
                                sig_coeffts['nq'],\
                                sig_coeffts['lq'],\
                                E)


  elif pi_type == const.CLARK:
    iE = numpy.where(E > sig_coeffts['I_e'])[0]
    if len(iE) > 0:
      x = E[iE]/sig_coeffts['I_e']
      term1 = 1/((sig_coeffts['Zel'] + sig_coeffts['b1'] + \
                  sig_coeffts['d1']/sig_coeffts['Zel']) * \
                  (sig_coeffts['Zel'] + sig_coeffts['b1'] +\
                   sig_coeffts['d1']/sig_coeffts['Zel']))

      sum1=0.0
      for ii in numpy.arange(sig_coeffts['n1'], dtype=int):
          sum1 += sig_coeffts['c'][ii]* \
                  x**(sig_coeffts['ip']+ii)
 
      term2 = 1/((sig_coeffts['Zel'] + sig_coeffts['b2'] + \
                  sig_coeffts['d2']/sig_coeffts['Zel']) * \
                  (sig_coeffts['Zel'] + sig_coeffts['b2'] +\
                   sig_coeffts['d2']/sig_coeffts['Zel']))
      sum2=0.0
      for ii in numpy.arange(sig_coeffts['n1'],4, dtype=int):
          sum2 += sig_coeffts['c'][ii]* \
                  x**(sig_coeffts['ip']+ii)

      result[iE] = term1*sum1 + term2*sum2

  elif pi_type == const.XSTAR:

    tmp2=numpy.interp(numpy.log(E), \
                      numpy.log(sig_coeffts['energy']),\
                      numpy.log(sig_coeffts['pi_param']),\
                      left=numpy.nan,right=numpy.inf)

    # deal with the edge cases: zero below ionization potential,
    # extrappolate as E^-3 at high energy
    inan = numpy.isnan(tmp2)
    iinf = numpy.isinf(tmp2)
    ifin = numpy.isfinite(tmp2)
    result[inan][:] = 0.0
    result[iinf] = sig_coeffts['pi_param'][-1]* \
                              ((E[iinf]/\
                                sig_coeffts['energy'][-1])**-3.0)*1e-18
    result[ifin] = 1e-18  * numpy.exp(tmp2[ifin])

  else:
    print "Error"

  return result



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#def calc_ionrec_ea(Z, Te, Ne, \
          #par_type, ionrec_par, min_temp, max_temp, temperatures):
  ## Te is in K          
            
  
  #T_eV= Te*const.KBOLTZ*1e3
  #ret = numpy.zeros(len(T_eV), dtype=float)
  #if (par_type == const.EA_MAZZOTTA_IRON):
    #for i in range(len(T_eV)):
      #y=ionrec_par[0]/T_eV[i]
      #if (y < 50.0):
         #f1_val = f1_fcn(y)
         #ea = (6.69e7/numpy.sqrt(T_eV[i]))*numpy.exp(-y)*1.e-16*\
             #(ionrec_par[1]+\
              #ionrec_par[2]*(1-y*f1_val)+\
              #ionrec_par[3]*(1-y*(1-y*f1_val))+\
              #ionrec_par[4]*(1-0.5*(y - y*y + y*y*y*f1_val))+\
              #ionrec_par[5]*f1_val)
         #ret[i] = ea
    
    #return ret
  
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



