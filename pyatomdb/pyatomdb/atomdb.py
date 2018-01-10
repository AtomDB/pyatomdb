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

import os, datetime, numpy, re, time, getpass
import util, atomic, spectrum, const, urllib2, apec
import astropy.io.fits as pyfits
from scipy import stats, integrate

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def write_filemap(d, filemap, atomdbroot=''):
  """
  Write filemap to file

  Parameters
  ----------
  d : dict
    Dictionary with filemap data in it. Structure defined as return value from
    read_filemap.
  filemap : str
    Name of filemap file to read. If zero length, use "$ATOMDB/filemap"
  atomdbroot : str
    Replace any $ATOMDB in the file names with this. If not provided,
    use "ATOMDB" environment variable instead

  Returns
  -------
  none
  """

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

#def read_filemap(filemap=False, atomdbroot=False):
  #"""
  #Reads AtomDB filemap data into dictionary

  #Parameters
  #----------
  #filemap : str
    #Name of filemap file to read. If zero length, use "$ATOMDB/filemap"
  #atomdbroot : str
    #Replace any $ATOMDB in the file names with this. If not provided,
    #use "ATOMDB" environment variable instead

  #Returns
  #-------
    #dict
      #dictionary listing all the filenames for each different datatype.
      #Files separated into "misctypes" for abundance and bremsstrahlung files,
      #and 'ionfiles', which are then split by element, ion, and datatype, e.g.
      #ret['ionfiles'][6][5]['lv'] contains the energy level (lv) data for
      #carbon V.
  #"""

  #if not(atomdbroot):
    #try:
      #atomdbroot = os.environ['ATOMDB']
    #except KeyError:
      #print '*** ERROR: atomdbroot not set, use environment variable ATOMDB or '+\
          #'pass atomdbroot to read_filemap'
      #return False

  #if not(filemap):
    #filemap = atomdbroot+'/filemap'
  #if not os.path.exists(filemap):
    #print "*** ERROR: Filemap %s does not exist, cannot fetch file" %(filemap)
    #return False

  #f = open(filemap,'r')

  #filedtype = numpy.dtype({'names':['Z','z1','ec','lv','ir','pc','dr',\
                                    #'la','em','pi','ai'],\
                           #'formats':[int, int, '|S160','|S160','|S160',\
                                      #'|S160','|S160','|S160','|S160',\
                                      #'|S160','|S160','|S160']})

  #miscdtype = numpy.dtype({'names':['misc_type','file'],\
                           #'formats':[int, '|S160']})

  #ret = {}
  #ret['ionfiles'] = numpy.zeros(0, dtype=filedtype)
  #ret['miscfiles'] = numpy.zeros(0, dtype=miscdtype)


  #for i in f:
    #splt = i.split()
    #fname = re.sub('\$ATOMDB',atomdbroot,splt[3])
    #Z_tmp = int(splt[1])
    #z1_tmp = int(splt[2])

    #if Z_tmp < 1:
      ## in this case, we have a "misc" datatype, not corresponding to a particular ion
      #misc_type = int(splt[0])
      #misc_file = fname
      #ret['miscfiles'] = numpy.append(ret['miscfiles'], numpy.array((misc_type, misc_file),\
                                                                    #dtype=miscdtype))

    #else:
      #j = numpy.where((ret['ionfiles']['Z']==Z_tmp) & \
                      #(ret['ionfiles']['z1']==z1_tmp))[0]
      #if len(j)==0:
        #j = len(ret['ionfiles'])

        #ret['ionfiles'] = numpy.append(ret['ionfiles'], numpy.zeros(1,\
                                                                    #dtype=filedtype))
        #ret['ionfiles']['Z'][j]=Z_tmp
        #ret['ionfiles']['z1'][j]=z1_tmp
      #else:
        #j = j[0]

      #if int(splt[0]) == 1:
        #ret['ionfiles']['ir'][j] = fname
      #if int(splt[0]) == 2:
        #ret['ionfiles']['lv'][j] = fname
      #if int(splt[0]) == 3:
        #ret['ionfiles']['la'][j] = fname
      #if int(splt[0]) == 4:
        #ret['ionfiles']['ec'][j] = fname
      #if int(splt[0]) == 5:
        #ret['ionfiles']['pc'][j] = fname
      #if int(splt[0]) == 6:
        #ret['ionfiles']['dr'][j] = fname
      #if int(splt[0]) == 7:
        #ret['ionfiles']['em'][j] = fname
      #if int(splt[0]) == 8:
        #ret['ionfiles']['pi'][j] = fname
      #if int(splt[0]) == 9:
        #ret['ionfiles']['ai'][j] = fname
##      if int(splt[0]) == 10:
##        cilist[j] = fname


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
  T_eV = 1.e3*const.KBOLTZ*Te
  dr = (c[0]/(T_eV**1.5))*(c[5]*numpy.exp(-c[1]/T_eV) +\
                               c[6]*numpy.exp(-c[2]/T_eV) +\
                               c[7]*numpy.exp(-c[3]/T_eV) +\
                               c[8]*numpy.exp(-c[4]/T_eV))
  return dr
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



def dr_badnell(Te, c):
  """
  Convert data from Badnell constants into a DR Rate

  Parameters
  ----------
  Te : float or array(float)
    Electron temperature[s] in K

  c : array
    Constants from DR rates. Stored as alternating pairs in AtomDB, so
    c1,e1,c2,e2,c3,e3 etc in the IONREC_PAR column

  Returns
  -------
  float
    DR rate in cm^3 s-1

  References
  ----------
  See http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/
  """

  Te_in, wasvec = util.make_vec(Te)

  ret = numpy.zeros(len(Te_in))
  for i in range(len(c)/2):
    if c[2*i] != 0.0:
      ret += c[2*i] * numpy.exp(-c[2*i+1]/Te_in)
  ret *= Te_in**-1.5


  if wasvec==False:
    ret=ret[0]


  return ret

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def rr_badnell(Te, c):
  """
  Convert data from Badnell constants into a RR Rate

  Parameters
  ----------
  Te : float or array(float)
    Electron temperature[s] in K

  c : array
    Constants from DR rates. Stored as alternating pairs in AtomDB, so
    c1,e1,c2,e2,c3,e3 etc in the IONREC_PAR column

  Returns
  -------
  float
    RR rate in cm^3 s-1

  References
  ----------
  See http://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/
  """

  Te_in, wasvec = util.make_vec(Te)

  ret = numpy.zeros(len(Te_in))

  A = c[0]
  B = c[1]
  T0 = c[2]
  T1=c[3]
  C = c[4]
  T2 = c[5]

  if C != 0:
    B += C*numpy.exp(-T2/Te_in)

  ret = A/ ((Te_in/T0)**0.5*\
            (1+(Te_in/T0)**0.5)**(1-B) *\
            (1+(Te_in/T1)**0.5)**(1+B))

  if wasvec==False:
    ret=ret[0]

  return ret

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



  a=get_data(1,1,'IONBAL')
#  find data type
  for i in a[1].header:
    print i
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

def get_filemap_file(ftype, Z, z1, fmapfile="$ATOMDB/filemap",\
                     atomdbroot="$ATOMDB", quiet=False, misc=False):
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
    for non ion-specific data, therefore Z,z1 are ignored.
    types are:
    10 or 'abund': elemental abundances
    11 or 'hbrems': Hummer bremstrahlung gaunt factor coefficients
    13 or 'rbrems': Relativistic bremstrahlung gaunt factor coefficients


  Returns
  -------
  str
    The filename for the relevant file, with all $ATOMDB expanded. If
    no file exists, returns zero length string.
  """

  fmap=read_filemap(filemap=fmapfile, atomdbroot=atomdbroot)


  if misc:
    # looking for type 10,11 or 13 data
    if ftype.lower() in ['10','11','13','abund','hbrems','rbrems']:
      if ftype.lower() in ['abund','hbrems','rbrems']:
        if ftype.lower() == 'abund':
          ftype_misc = 10
        if ftype.lower() == 'hbrems':
          ftype_misc = 11
        if ftype.lower() == 'rbrems':
          ftype_misc = 13
      else:
        ftype_misc = int(ftype)

      i = numpy.where(fmap['misc_type']==ftype_misc)[0]
      if len(i)==0:
        print "Error: file type: %i not found in filemap %s" %(ftype,fmapfile)
        ret = ''
      else:
        ret = fmap['misc'][i[0]]
    else:
      if not ftype.lower() in ['eigen','ionbal']:
        print "Error: unknown file type: %s not recognized" %(ftype)
      ret = ''
  else:
    i = numpy.where((fmap['Z']==Z)&(fmap['z1']==z1))[0]
    ret=''
    if len(i)==0:
      if not quiet :
        print "WARNING: there is no data for the ion "+\
               atomic.spectroscopic_name(Z,z1)
      ret=''

    if len(i)>1:
      print "ERROR: there are multiple entries for the ion "+\
             atomic.spectroscopic_name(Z,z1)
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
                atomic.spectroscopic_name(Z,z1)

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


def get_abundance(abundfile=False, abundset='AG89', element=[-1],\
                  datacache=False, settings = False):
  """
  Get the elemental abundances, relative to H (H=1.0)

  Parameters
  ----------
  abundfile : string
    special abundance file, if not using the default from filemap
  abundset : string
    Abundance set. Available:

    * Allen: Allen, C. W.  Astrophysical Quantities, 3rd Ed.,  1973 (London: Athlone Press)

    * AG89: Anders, E. and Grevesse, N. 1989, Geochimica et Cosmochimica Acta, 53, 197

    * GA88: Grevesse, N, and Anders, E.1988, Cosmic abundances of matter, ed. C. J. Waddington, AIP Conference, Minneapolis, MN

    * Feldman: Feldman, U., Mandelbaum, P., Seely, J.L., Doschek, G.A.,Gursky H., 1992, ApJSS, 81,387

    Default is AG89
  element : list of int
    Elements to find abundance for. If not specified, return all.
  datacache : dict
    See get_data
  datacache : settings
    See get_data

  Returns
  -------
  dict
    abundances in dictionary, i.e :

    {1: 1.0,\n
     2: 0.097723722095581111,\n
     3: 1.4454397707459272e-11,\n
     4: 1.4125375446227541e-11,\n
     5: 3.9810717055349735e-10,\n
     6: 0.00036307805477010178,...\n

  """

  if not abundfile:
    abunddata = get_data(False, False, 'abund', \
                                datacache=datacache,\
                                settings = settings)
  else:
    abunddata = pyfits.open(abundfile)

  if element[0]==-1:
    element = range(1,31)

  ind = numpy.where(abunddata[1].data.field('Source')==abundset)[0]
  if len(ind)==0:
    print "Invalid Abundance Set chosen: select from ", \
          abunddata[1].data.field('Source')
    return -1
  ret = {}
  for Z in element:
    elsymb = atomic.Ztoelsymb(Z)
    ret[Z]=10**(abunddata[1].data.field(elsymb)[ind[0]])/1e12
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
                         ion_drv=False, elem_drv=False, use_nei=False):

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
  use_nei : bool
    This can be useful when trying to get line emissivities which fall
    below the 1e-20 cut off. Applying this flag, the NEI file will be used
    by default and an ionization balance applied. This should give the
    same results as normal for strong emissivities, but go to a lower
    emissivity before being set to zero. Use with caution...

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
  # Version 0.2 - Added use_nei
  # Adam Foster 28 Apr 2017

  #


  if use_nei:
    if linefile=="$ATOMDB/apec_line.fits":
      linefile = "$ATOMDB/apec_nei_line.fits"

  a = pyfits.open(os.path.expandvars(linefile))

  kT = a[1].data.field('kT')
  dens = a[1].data.field('eDensity')
  time = a[1].data.field('time')

  epsilon = numpy.zeros(len(kT), dtype=float)

  datacache={}

  for ikT in range(len(kT)):

    iikT = ikT + 2
    if use_nei:
      j = numpy.where((a[iikT].data.field("element") == Z) &\
                      (a[iikT].data.field("ion") == z1) &\
                      (a[iikT].data.field("UpperLev") == upind) &
                      (a[iikT].data.field("LowerLev") == loind))[0]
      if len(j) == 0: continue
      ionbal = apec.solve_ionbal_eigen(Z,kT[ikT],teunit='keV', datacache=datacache)
      for jj in j:
        print ikT, a[iikT].data['Ion_drv'][jj], a[iikT].data['Epsilon'][jj], ionbal[a[iikT].data['Ion_drv'][jj]-1]

        epsilon[ikT]+= a[iikT].data['Epsilon'][jj] * ionbal[a[iikT].data['Ion_drv'][jj]-1]
    else:
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
  bttype = get_burgess_tully_transition_type(levdat[1].data[lolev-1],\
                                             levdat[1].data[uplev-1],\
                                             Aval)


    # do the extrappolation etc
  btval = get_burgess_tully_extrap(bttype, \
                                   levdat[1].data[lolev-1], \
                                   levdat[1].data[uplev-1], \
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
                       force_extrap=False, did_extrap=False, \
                       datacache=False):

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
      y2 = prep_spline_atomdb(xs, om[2:2+int(om[0])], 5)
      upsilon = calc_spline_atomdb(xs, om[2:2+int(om[0])], y2, 5, st)
    elif int(om[0])== 9:
      upsilon = interpolate.interp1d(xs9, om[2:2+int(om[0])], kind='cubic',\
       bounds_error=False, fill_value=0.0)(st)
      y2 = prep_spline_atomdb(xs9, om[2:2+int(om[0])], 9)
      stvec, isstvec = util.make_vec(st)
      upsilon = numpy.zeros(len(stvec))
      for ist, st in enumerate(stvec):
        upsilon[ist] = calc_spline_atomdb(xs9, om[2:2+int(om[0])], y2, 9, st)
      if isstvec==False:
        upsilon = upsilon[0]

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

      # this is for values outside the range
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
  ci_in = numpy.double(cidat['temperature'][:N_interp])
  tmpci = numpy.double(tmpci)

  tmp = numpy.exp(interpolate.interp1d(numpy.log(ci_in), \
                             numpy.log(tmpci+1e-30), \
                             kind=1, bounds_error=False,\
                             fill_value=numpy.nan)(numpy.log(Te)))-1e-30

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

def calc_ionrec_ci(cidat, Te, extrap=False, ionpot=False):
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
  if len(ici) > 0:
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

    elif ((cidat['par_type']> const.CI_DERE) &\
          (cidat['par_type']<= const.CI_DERE+20)):
      npts = cidat['par_type']-const.CI_DERE
      ci[ici] = calc_ci_dere(Te[ici], ionpot, cidat['Temperature'][:npts], \
                             cidat['ionrec_par'][:npts])


    else:
      print "Unknown CI type: %i"%(cidat['par_type'])


  # now see if extrappolation is required:
  if extrap:
    #  get Te, coeff for the 2 lowest points.
    ilow = numpy.where(Te<cidat['min_temp'])[0]

    if len(ilow) > 0:
      # calculate the value at Te_min
      Te_min = numpy.array([cidat['min_temp']])

      cimin = calc_ionrec_ci(cidat, Te_min, ionpot=ionpot)
      # if log of this is < 46.0, this is a small number, just repeat this
      if numpy.log(cimin[0])<46.0:
        ci[ilow]=cimin

      # otherwise, calculate the value at a range of near-minimum temperatures,
      # and use to construct a good second derivative.
      else:
        tetmp = numpy.logspace(numpy.log10(Te_min), numpy.log10(Te_min)+1,4)

        citmp=calc_ionrec_ci(cidat, tetmp, ionpot=ionpot)

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
      cimax = calc_ionrec_ci(cidat, Te_max, ionpot=ionpot)
      ci[ihigh]=cimax[0]*(Te_max[0]/Te[ihigh])**0.5



  return ci
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def calc_ci_dere(Te, ionpot, Tscal, Upsscal):

  """
  Calculate the collisional ionization rates using the Dere 2007 method

  Parameters
  ----------
  Te : float or array(float)
    Electron temperature (K)
  ionpot : float
    Ionization potential (eV)
  Tscal : array(float)
    scaled temperatures
  Upsscal : array(float)
    scaled upsilons

  Returns
  -------
  float or array(float)
    Ionization rate in cm^3 s^-1

  References
  ----------
  2007A&A...466..771D
  """
  from scipy import interpolate
  from scipy.special import exp1

  tin, wasvec = util.make_vec(Te)
  tinscal = tin*const.KBOLTZ*1000/ionpot

  f = 2.0

  xin = (1 - numpy.log10(f) / numpy.log10(tinscal+f))
  xdat = Tscal
  ydat = Upsscal

  tck = interpolate.splrep(xdat,ydat)
  yout = interpolate.splev(xin,tck)
  R = tinscal**(-0.5) * ionpot**(-1.5) * exp1(1/tinscal)*yout

  if not wasvec:
    R=R[0]
  return R


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
  # Badnell
  elif (cidat['par_type'] == const.RR_BADNELL):
    rr[irr] = rr_badnell(Te[irr],cidat['ionrec_par'])


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

  """
  Calculate the DR rate at temperature Te for a transition cidat

  Parameters
  ----------

  cidat: array
    A line from the IR file with the data on the transition to be calcualted
  Te: array(float)
    Electron temperature of interest (in K)
  extrap : bool
    Whether to perform extrappolation or not


  Returns
  -------
  float
    The DR rate in cm^3 s^-1
  """

  dr = numpy.zeros(len(Te))
  # set values outside range to NAN
  idr = numpy.where((Te >= cidat['min_temp']) & (Te <= cidat['max_temp']))[0]
  dr[:] = numpy.nan
  dr[idr] = 0.0

  if len(idr) > 0:
    # Mazzotta
    if (cidat['par_type'] == const.DR_MAZZOTTA):
      T_eV = 1.e3*const.KBOLTZ*Te[idr]
      dr[idr] = (cidat['par_type'][0]/(T_eV**1.5))  * \
           (cidat['par_type'][5]*numpy.exp(-cidat['par_type'][1]/T_eV) +\
            cidat['par_type'][6]*numpy.exp(-cidat['par_type'][2]/T_eV) +\
            cidat['par_type'][7]*numpy.exp(-cidat['par_type'][3]/T_eV) +\
            cidat['par_type'][8]*numpy.exp(-cidat['par_type'][4]/T_eV))

    elif (cidat['par_type'] == const.DR_BADNELL):
      dr[idr] = dr_badnell(Te[idr], cidat['ionrec_par'])

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
                     extrap=True):
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
    Te /=const.KBOLTZ
    Te *=1000.0
  elif Te_unit.lower()=='kev':
    Te /=const.KBOLTZ
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
def get_maxwell_rate(Te, colldata=False, index=-1, lvdata=False, Te_unit='K', \
                     lvdatap1=False, ionpot = False, \
                     force_extrap=False, silent=True,\
                     finallev=False, initlev=False,\
                     Z=-1, z1=-1, dtype=False, exconly=False,\
                     datacache=False, settings=False):
  """
  Get the maxwellian rate for a transition from a file, typically for ionization,
  recombination or excitation.

  Parameters
  ----------

  Te : float
    electron temperature(s), in K by default
  colldata : HDUList
    If provided, the HDUList for the collisional data
  index : int
    The line in the HDUList to do the calculation for. Indexed from 0.
  lvdata : HDUList
    the hdulist for the energy level file (as returned by pyfits.open('file'))
  Te_unit : {'K' , 'eV' , 'keV'}
    Units of temperature grid.
  lvdatap1 : HDUList
    The level data for the recombining or ionized data.
  ionpot : float
    The ionization potential in eV (required for some calculations, if
    not provided, it will be looked up)
  force_extrap : bool
    Force extrappolation to occur for rates outside the nominal range
    of the input data
  silent : bool
    Turn off notifications
  finallev : int
    Instead of specifying the index, can use upperlev, lowerlev instead.
  initlev : int
    Instead of specifying the index, can use upperlev, lowerlev instead
  Z : int
    Instead of providing colldata, can provide Z & z1. Z is the atomic
    number of the element.
  z1 : int
    Instead of providing colldata, can provide Z & z1. z1 is the ion
    charge +1 for the initial ion
  dtype : str
    data type. One of:\n
    'EC' : electron impact excitation\n
    'PC' : proton impact excitation\n
    'CI' : collisional ionization\n
    'EA' : excitation-autoionization\n
    'XI' : excluded ionization\n
    'XR' : excluded recombination\n
    'RR' : radiative recombination\n
    'DR' : dielectronic recombination
  exconly : bool
    For collisional excitation, return only the excitation rate, not
    the de-excitation rate.
  settings : dict
    See description in read_data
  datacache : dict
    See description in read_data

  Returns
  -------
  float or array(float)
    Maxwellian rate coefficient, in units of cm^3 s^-1
    For collisional excitation (proton or electron) returns
    excitation, dexcitation rates

  Examples
  --------
  Te = numpy.logspace(4,9,20)

  (1) Get excitation rates for row 12 of an Fe XVII file
  colldata = pyatomdb.atomdb.get_data(26,17,'EC')
  exc, dex = get_maxwell_rate(Te, colldata=colldata, index=12)

  (2) Get excitation rates for row 12 of an Fe XVII file
  exc, dex = get_maxwell_rate(Te, Z=26,z1=17, index=12)

  (3) Get excitation rates for transitions from level 1 to 15 of FE XVII
  exc, dex = get_maxwell_rate(Te, Z=26, z1=17, dtype='EC', finallev=15, initlev=1)

  """
# Note interface update 03-Apr-2016

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
  elif Te_unit.lower() == 'kev':
    Te_arr = Te_arr*11604.505*1000.0
  elif Te_unit.lower() == 'k':
    Te_arr = 1.0* Te_arr
  else:
    print 'ERROR: Unknown Te_unit "%s": should be "K" or "keV"' % (Te_unit)
    return False

  # CHECK THE INPUTS
  # 1: the collional excitation data
  if colldata:
    if Z>=0 | z1 >=0:
      print 'ERROR: specified colldata, Z and z1. Either specify colldata or Z, z1 & dtype'
      return False
    # set the dtype
    if dtype==False:
      if colldata[1].header['HDUCLAS1']=='COLL_STR':
        dtype='EC'
    # get the Z, z1
    if dtype=='EC':
      Z=colldata[1].header['ELEMENT']
      z1=colldata[1].header['ION_STAT']+1



  else:
    # Check we have Z & z1 & dtype specified
    if ((dtype != False) & (Z>0) & (z1>0)):
      # we should have data...
      if dtype in ['CI','EA','XI']:
        colldata = get_data(Z,z1,'IR', settings=settings, datacache=datacache)
      elif dtype in ['RR','DR','XR','XD']:
        colldata = get_data(Z,z1-1,'IR', settings=settings, datacache=datacache)
      elif dtype == 'EC':
        colldata = get_data(Z,z1,'EC', settings=settings, datacache=datacache)
      elif dtype == 'PC':
        colldata = get_data(Z,z1,'PC', settings=settings, datacache=datacache)
      else:
        print "Error: unknown dtype %s"%(dtype)
        return False


    else:
      print "Error: must specify dtype, Z and z1 to load file"
      return False
  # OK, at this point we have the colldata. Yay.
  # check colldata is real
  if colldata==False:
    print "No collisional data found. Returning zeros"
    ret = numpy.zeros(len(Te_arr), dtype=float)
    if dtype in ['EC','PC']:
      if not(isiter):
        ret = ret[0]
      if exconly:
        return ret
      else:
        return ret, ret
    else:
      return ret

  # now to sort out Z, z1
  if Z<0:
    Z = colldata[1].header['ELEMENT']
  if z1<0:
    if dtype in ['EC','PC','EA','CI','XI']:
      z1 = colldata[1].header['ION_STAT']+1
    elif dtype in ['XR','DR','RR','XD']:
      z1 = colldata[1].header['ION_STAT']+2

  # Now get the correct transition
  if index>=0:
    if (finallev | initlev):
      print "Error: specify index or upperlev and lowerlev"
      return False
    # If I have the index, set the dtype
    if colldata[1].header['HDUCLAS1'] == 'COLL_STR':
      dtype = 'EC'
    elif colldata[1].header['HDUCLAS1'] == 'IONREC':
      # get the type:
      dtype = colldata[1].data.field('TR_TYPE')[index]
    else:
      print "ERROR: supplied data is not a collision strength (EC) or "+\
            "ionization/recombination (IR) set. Returning"
      return False

    if z1<0:
      if dtype in ['EC','PC','EA','CI','XI']:
        z1 = colldata[1].header['ION_STAT']+1
      elif dtype in ['XR','DR','RR','XD']:
        z1 = colldata[1].header['ION_STAT']+2


  else:

    if (finallev==False | initlev == False):
      print "Error: if not specifying index, must specify finallev and ",
      print "initlev"
      return False
    else:

      if dtype=='EC':
        index = numpy.where((colldata[1].data['lower_lev']==initlev) &\
                            (colldata[1].data['upper_lev']==finallev))[0]
        if len(index)==0:
          print "Warning: no data found for electron excitation from "+\
                "level %i to %i"%(initlev, finallev)
          ret = numpy.zeros(len(Te_arr), dtype=float)
          if not isiter:
            ret = ret[0]
          if exconly:
            return ret
          else:
            return ret, ret
        else:
          index = index[0]
      elif dtype =='PC':
        index = numpy.where((colldata[1].data['lowerlev']==initlev) &\
                            (colldata[1].data['upperlev']==finallev))[0]
        if len(index)==0:
          print "Warning: no data found for proton excitation from "+\
                "level %i to %i"%(initlev, finallev)
          ret = numpy.zeros(len(Te_arr), dtype=float)
          if not isiter:
            ret = ret[0]
          if exconly:
            return ret
          else:
            return ret, ret
        else:
          index = index[0]
      elif dtype in ['CI','XI','EA']:
        index = numpy.where((colldata[1].data['level_init']==initlev) &\
                            (colldata[1].data['level_final']==finallev) &\
                            (colldata[1].data['tr_type']==dtype) &\
                            (colldata[1].data['ion_init']==z1) &\
                            (colldata[1].data['ion_final']==z1+1))[0]
        if len(index)==0:
          print "Warning: no data found for collisional ionization from "+\
                "ion %i, level %i to ion %i, level %i"%\
                 (z1, initlev, z1+1, finallev)
          ret = numpy.zeros(len(Te_arr), dtype=float)
          if not isiter:
            ret = ret[0]
          return ret
        else:
          index = index[0]

      elif dtype in ['RR','XR','DR','XD']:
        index = numpy.where((colldata[1].data['level_init']==initlev) &\
                            (colldata[1].data['level_final']==finallev) &\
                            (colldata[1].data['tr_type']==dtype) &\
                            (colldata[1].data['ion_init']==z1) &\
                            (colldata[1].data['ion_final']==z1-1))[0]
        if len(index)==0:
          print "Warning: no data found for collisional ionization from "+\
                "ion %i, level %i to ion %i, level %i"%\
                 (z1, initlev, z1-1, finallev)
          ret = numpy.zeros(len(Te_arr), dtype=float)
          if not isiter:
            ret = ret[0]
          return ret
        else:
          index = index[0]

  # convert the data.

  if dtype=='EC':
    # go through all the different possibilities for collisional excitation data.
    ecdat = colldata[1].data[index]
    upind = ecdat['upper_lev']
    loind = ecdat['lower_lev']

    if not(lvdata):
      lvdata = get_data(Z,z1,'LV', settings=settings, datacache=datacache)

    uplev = lvdata[1].data[upind-1]
    lolev = lvdata[1].data[loind-1]

    delta_E = uplev['energy']-lolev['energy']
    Ztmp = lvdata[1].header['ELEMENT']

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
                                 delta_E/1e3, Te_arr, Ztmp, degl, degu)

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

    if ((cidat['par_type']>const.INTERP_I_UPSILON) & \
        (cidat['par_type']<=const.INTERP_I_UPSILON+20)):



      upind = cidat['upper_lev']
      loind = cidat['lower_lev']
      if not(lvdata):
        lvdata = get_data(z,z1,'LV', settings=settings, datacache=datacache)
      if not(lvdatap1):
        lvdatap1 = get_data(z,z1+1,'LV', settings=settings, datacache=datacache)

      uplev = lvdatap1[1].data[upind-1]
      lolev = lvdata[1].data[loind-1]

      # get the ionization potential
      if not(ionpot):
#        print "fixing ionpot"
        ionpot = colldata[1].header['IONPOT']
#      else:
#        print "I don't need to fix ionpot, it is ", ionpot

      delta_E = ionpot + uplev['energy']-lolev['energy']
      Ztmp = lvdata[1].header['ELEMENT']
      degl = lolev['lev_deg']
      degu = uplev['lev_deg']


      ci, dex = calc_maxwell_rates(cidat['par_type'],\
                         cidat['min_temp'],\
                         cidat['max_temp'],\
                                 cidat['temperature'],\
                                 cidat['ionrec_par'],\
                                 delta_E/1e3, Te_arr, Ztmp, degl, degu)

    elif ((cidat['par_type']>const.CI_DERE) &\
          (cidat['par_type']<=const.CI_DERE+20)):
      ionpot = float(colldata[1].header['IP_DERE'])
      ci = calc_ionrec_ci(cidat, Te_arr, extrap=force_extrap, ionpot=ionpot)
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
    ea = calc_ionrec_ea(cidat,Te_arr, extrap=force_extrap)
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
    dr = calc_ionrec_dr(cidat,Te_arr, extrap=force_extrap)
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
    rr = calc_ionrec_rr(cidat,Te_arr, extrap=force_extrap)
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
    xr = calc_ionrec_rr(cidat,Te_arr, extrap=force_extrap)
    if sum(numpy.isnan(xr))>0:
      if not silent:

        print "calc_ionrec_rate: xr(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                   Te,xr)

    return xr

  elif dtype=='XD':
    cidat = colldata[1].data[index]
    xr = calc_ionrec_dr(cidat,Te_arr, extrap=force_extrap)
    if sum(numpy.isnan(xr))>0:
      if not silent:

        print "calc_ionrec_rate: xd(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                   Te,xr)

    return xr


  elif dtype=='XI':
    cidat = colldata[1].data[index]

    if ((cidat['par_type']>900) & \
        (cidat['par_type']<=920)):

      # sanity check
      if not(lvdata):
#        print "calling get lvdata"
#        print lvdata, Z, z1
        lvdata = get_data(Z,z1,'LV', settings=settings, datacache=datacache)
#        print lvdata

      if not(lvdatap1):
        lvdatap1 = get_data(Z,z1+1,'LV', settings=settings, datacache=datacache)

      if (lvdata[1].header['ion_stat']+1 != cidat['ion_init']):
        print "ERROR: lvdata and cidat not for matching ions!"
      if (lvdatap1[1].header['ion_stat']+1 != cidat['ion_final']):
        print "ERROR: lvdatap1 and cidat not for matching ions!"

      upind = cidat['level_final']
      loind = cidat['level_init']


      uplev = lvdatap1[1].data[upind-1]
      lolev = lvdata[1].data[loind-1]

      # get the ionization potential
      if not(ionpot):
        ionpot = colldata[1].header['IONPOT']

      delta_E = ionpot + uplev['energy']-lolev['energy']

      Ztmp = lvdata[1].header['ELEMENT']
      degl = lolev['lev_deg']
      degu = uplev['lev_deg']


      xi, dex = calc_maxwell_rates(cidat['par_type'],\
                         cidat['min_temp'],\
                         cidat['max_temp'],\
                         cidat['temperature'],\
                         cidat['ionrec_par'],\
                         delta_E/1e3, Te_arr, Ztmp, degl, degu)

    else:
      xi = calc_ionrec_ci(cidat,Te_arr, extrap=force_extrap)
      if (xi < 0.0):
        print "calc_ionrec_rate: CI(%10s -> %10s,T=%9.3e) = %8g"%\
                  (adbatomic.spectroscopic_name(cidat['element'],cidat['ion_init']),\
                   adbatomic.spectroscopic_name(cidat['element'],cidat['ion_final']),\
                   Te,xi)
        xi=0.0
    return xi


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def sigma_hydrogenic(Z, N,L,  Ein):
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
#  print chi, Z, n
#  zzz=raw_input()
  Eelec = (E - chi)/RYDBERG
  sigma = numpy.zeros(len(E), dtype=float)
  if N>5:
    return sigma
  iE = numpy.where(Eelec > 0)[0]

  if len(iE > 0):

    eta = Z/numpy.sqrt(Eelec[iE])
    rho = eta / n
    lp1 = l+1.0
#    print 'eta', eta
    # Exponential term
    expterm = numpy.exp(-4*eta*numpy.arctan(1./rho)) / (1-numpy.exp(-2*numpy.pi*eta))
    # Common Factorial
    commfact1 = 1.
    commfact2 = 1.

    if (n+l >= 2*l+2):
      for iFact in range(2*l+3, n+l+1):
        commfact1 *= iFact
    else:
      for iFact in range(n+l+1, 2*l+3):
        commfact2 *= iFact

    for iFact in range(2,2*l+2):
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
#    for i in range(len(iE)):
#      print "%6i %e %e %e %e"%(iE,Ein[iE[i]], sigma[iE[i]], sigma_m1, sigma_p1)
  return sigma

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def G_hyd(l,m, eta, rho):


  result=0

  B_s = B_hyd(2*m, l, m, eta)
  result = B_s

  for s in range(2*m-1,-1,-1):
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
             xstarlevfinal=1,  settings=False, datacache=False,\
             returntotal = False):

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
    The electron temperature (K)
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
  returntotal : bool
    If true, return the total recombination rate as well

  Returns
  -------
  array(float)
    The rrc in photons cm^3 s^-1 keV^-1
  optional float
    If returntotal is set, also return total RRC calculated by
    separate integral from the ionization edge to infinity.
  """
  # Version 0.1 Initial Release
  #
  # Version 0.2
  # Fixed unit issues
  # Adam Foster 2015-Oct-23

  #OK, let's get the line data.

  rrc = numpy.zeros(len(eedges)-1)

  # get temperature in keV
  kT = Te * const.KBOLTZ
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
      plvdat = get_data(Z, z1+1, 'LV', settings=settings, datacache=datacache)
      if plvdat:
        ginit = 1.0* plvdat[1].data['lev_deg'][0]
    gratio = gfin/ginit
    if ldat['n_quan']>5:
      if returntotal:
        return rrc, 0.0
      else:
        return rrc
  elif ldat['PHOT_TYPE']==const.CLARK:
    I_e = ldat['phot_par'][1]
    gfin = 1.0 * ldat['lev_deg']
    ginit = 1.0
    if Z-z1>0:
      plvdat = get_data(Z, z1+1, 'LV', settings=settings, datacache=datacache)
      if plvdat:
        ginit = 1.0* plvdat[1].data['lev_deg'][0]
    gratio = gfin/ginit
  elif ldat['PHOT_TYPE']==const.VERNER:
    I_e = ldat['phot_par'][0]
    gfin = 1.0 * ldat['lev_deg']
    ginit = 1.0
    if Z-z1>0:
      plvdat = get_data(Z, z1+1, 'LV', settings=settings, datacache=datacache)
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
    if xstardat==False:
      #no photoionization data
      gratio = 0
      I_e = 0
      rrc = numpy.zeros(len(eedges)-1,dtype=float)
      if returntotal:
        return rrc , 0.0
      else:
        return rrc

    gratio = xstardat['g_ratio']
    I_e = (get_ionpot(Z, z1) - ldat['energy'])/1000.0

  else:
    #no photoionization data
    rrc = numpy.zeros(len(eedges)-1,dtype=float)
    if returntotal:
      return rrc , 0.0
    else:
      return rrc

  rrc_ph_factor = (const.RRC_COEFF/const.ERG_KEV)*gratio/(Te**1.5)
  rrc_erg_factor = const.RRC_COEFF*gratio/(Te**1.5)

  # TOTAL INTEGRAL
  if returntotal:
    rr_lev_pop,err = integrate.quad(rrc_ph_value, I_e+1.e-4, numpy.inf,\
                                      epsabs=const.TOT_ABSACC,\
                                      epsrel=const.TOT_RELACC, \
                              args=(Z, z1, rrc_ph_factor, I_e, kT, ldat,
                                    xstardat, xstarlevfinal), limit=100000)
  # VALUE AT EACH BIN EDGE
  emission_edges = rrc_ph_value(eedges, Z, z1,  rrc_ph_factor, I_e, \
                                kT, ldat, xstardat, xstarlevfinal)


  edgebin = numpy.argmin(eedges>I_e)-1

  if edgebin > -1:
    # we have an edge to deal with
    rrc[edgebin] = integrate.quad(rrc_ph_value, I_e+1.e-4, eedges[edgebin+1],\
                                    epsabs=const.TOT_ABSACC,\
                                    epsrel=const.TOT_RELACC, \
                            args=(Z, z1, rrc_ph_factor, I_e, kT, ldat,
                                  xstardat, xstarlevfinal), limit=100000)
  rrc[edgebin+1:] = (emission_edges[edgebin+1:-1]+emission_edges[edgebin+2:])/2.0



  if returntotal:
    return rrc, rr_lev_pop
  else:
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
    temperautre (K)
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
  array(float)
    Recombination rates into the excited levels, in s-1
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
  LevelRecombRate2 = numpy.zeros(nlev, dtype=float)
  binwidth = ebins[1:]-ebins[:-1]

  for iLev in range(nlev):

    finlev = finallevels[1].data[iLev]
    if finlev['phot_type'] >=0:
      pass
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
      hasinitlevels =  util.keyword_check(initlevels)
      parent_lev_deg=1.0
      if hasinitlevels:
        if ((hasparent) & (len(initlevels[1].data)>0)):
          parent_lev_deg = initlevels[1].data['lev_deg'][0]*1.0
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

      # find the possible places to photoionize to:
      finlev_possible = util.unique(\
        pidat[1].data['lev_final'][pidat[1].data['lev_init']==xstarlevinit])

      finlev_possible = numpy.array([min(finlev_possible)])

      for xstarlevfinal in finlev_possible:
        # get the numbers
        sigma_coeff = sort_pi_data(pidat, xstarlevinit, xstarlevfinal)

        if not sigma_coeff:
          continue
        g_ratio = sigma_coeff['g_ratio']

        tmprrc,total = calc_rrc(Z, z1, ebins, T, iLev+1, xstardat=sigma_coeff, \
              settings=settings, datacache=datacache, returntotal=True)

        tmprrc *= abund*ion_pop*binwidth
#        print "Calculated RRC for Z %i  z1 %i  iLev %i"%\
#              (Z, z1, iLev)
#        for i in range(len(tmprrc)):
#          print ebins[i],ebins[i+1], tmprrc[i]
#        print "ENDCalculated RRC for Z %i  z1 %i  iLev %i"%\
#              (Z, z1, iLev)
        rr_lev_rate = sum(tmprrc)
        LevelRecombRate[iLev] += total*abund*ion_pop
        tot_rec_rate += total*abund*ion_pop
        rrc += tmprrc
    else:

      continue
    rr_lev_rate = 0.0
    g_ratio = lev_deg/parent_lev_deg

    if ((finlev['phot_type'] != const.NO_PHOT) &\
        (finlev['phot_type'] != const.XSTAR)):
      tmprrc, total = calc_rrc(Z, z1, ebins, T, iLev+1, \
              settings=settings, datacache=datacache, returntotal=True)
      tmprrc *= abund*ion_pop*binwidth

      rr_lev_rate = sum(tmprrc)
#      print "Calculated RRC for Z %i  z1 %i  iLev %i"%\
#            (Z, z1, iLev)
#      for i in range(len(tmprrc)):
#        print ebins[i],ebins[i+1], tmprrc[i]
#      print "ENDCalculated RRC for Z %i  z1 %i  iLev %i"%\
#            (Z, z1, iLev)

      LevelRecombRate[iLev] += total*abund*ion_pop

      tot_rec_rate += total*abund*ion_pop

      rrc += tmprrc
#      for i in xrange(len(tmprrc)):
        #print "%e %e %e" %(ebins[i], ebins[i+1], tmprrc[i])
      #print "phot_type= %i"%(finlev['phot_type'])
      #print finlev
#  for i in range(len(LevelRecombRate)):
#      print "%i %e %e %e"%(i+1, LevelRecombRate[i], LevelRecombRate2[i], LevelRecombRate[i]/ LevelRecombRate2[i])
#  zzz=raw_input('argh')
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
#
#  Version 0.2 - bugfix
#  added check for "filemap==False"
#  Adam Foster October 9th 2015
#

  # parse the options here.

  if filemap==False:
    filemap = "$ATOMDB/filemap"

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

      j=j[0]

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

       Or, for non-ion-specific data (abundances and bremstrahlung coeffts)
       *           'ABUND' - abundance tables
       *           'HBREMS' - Hummer bremstrahlung coefficients
       *           'RBREMS' - relativistic bremstrahlung coefficitients
       *           'IONBAL' - ionization balance tables
       *           'EIGEN'  - eigenvalue files

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

#  Version 0.3 - Fixed fmapfile and atomdbroot use to avoid returning
#  "False" and triggering errors
#  Adam Foster October 27th 2015

  d = False
  didurl=False
  fmapfile = "$ATOMDB/filemap"
  atomdbroot="$ATOMDB"

  # check if data type requested in "miscellanous", i.e. a bremstrahlung
  # or abundance related file, not an ion-specific file.

  ismisc = False
  if ftype.lower() in ['abund','hbrems','rbrems','ionbal','eigen']:
    ismisc = True


  if datacache != False:
    # make sure that the relevant dictionaries are ready to receive the data
    if not 'data' in datacache.keys():
      datacache['data']={}
    if not ismisc:
      if not Z in datacache['data'].keys():
        datacache['data'][Z]={}
      if not z1 in datacache['data'][Z].keys():
        datacache['data'][Z][z1]={}
    else:
      if not 'misc' in datacache['data'].keys():
        datacache['data']['misc'] = {}

    if not 'datasums' in datacache.keys():
      datacache['datasums']={}
    if not ismisc:
      if not Z in datacache['datasums'].keys():
        datacache['datasums'][Z]={}
      if not z1 in datacache['datasums'][Z].keys():
        datacache['datasums'][Z][z1]={}
    else:
      if not 'misc' in datacache['datasums'].keys():
        datacache['datasums']['misc'] = {}

    if ismisc:
      havedata = False

#      if not 'misc' in datacache['data'].keys():
#        print "ping1"
#        datacache['data']['misc']={}


      if ftype.upper() =='EIGEN':
        if not 'EIGEN' in datacache['data']['misc'].keys():
          datacache['data']['misc']['EIGEN']={}
          datacache['datasums']['misc']['EIGEN']={}
          havedata = False

      if ftype.upper() =='EIGEN':
        if Z in datacache['data']['misc']['EIGEN'].keys():
          havedata=True


      elif ((ftype.upper() in datacache['data']['misc'].keys()) &\
            (ftype.upper() != 'EIGEN')):
      # this means we have the data cached, no need to fetch it
        havedata=True

      if not(havedata):
        if settings:
          if settings['filemap']:
            fmapfile = settings['filemap']
          if settings['atomdbroot']:
            atomdbroot = settings['atomdbroot']

        fname = get_filemap_file(ftype, False, False, fmapfile=fmapfile,\
                               atomdbroot=atomdbroot, quiet=True,\
                               misc = True)
        if fname=='':
        # no data exists
          # This is expected if it's an ionbal file

          if ftype.lower()=='ionbal':
            # conversion here:
            curversion = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]

            if curversion in ['2.0.0', '2.0.1', '2.0.2','3.0.0','3.0.1','3.0.2','3.0.3']:
              fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v2.0.2_ionbal.fits'
            elif curversion in ['3.0.4','3.0.5','3.0.6']:
              fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v3.0.4_ionbal.fits'
            else:
              fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v3.0.7_ionbal.fits'
          elif ftype.lower()=='eigen':
            # conversion here:
            curversion = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]

            if curversion in ['2.0.0', '2.0.1', '2.0.2','3.0.0','3.0.1','3.0.2','3.0.3']:
              fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/eigen/eigen%s_v3.0.fits'%(atomic.Ztoelsymb(Z).lower())
            elif curversion in ['3.0.4','3.0.5','3.0.6']:
              fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/eigen/eigen%s_v3.0.4.fits'%(atomic.Ztoelsymb(Z).lower())
            else:
              fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/eigen/eigen%s_v3.0.7.fits'%(atomic.Ztoelsymb(Z).lower())
            if not 'EIGEN' in datacache['data']['misc'].keys():
              datacache['data']['misc'][ftype.upper()]={}
              datacache['datasums']['misc'][ftype.upper()]={}
          else:
            datacache['data']['misc'][ftype.upper()] = False

        if fname!='':
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
                url = re.sub(os.path.expandvars(atomdbroot),\
                             'ftp://sao-ftp.harvard.edu/AtomDB',fname)+'.gz'
                try:
                  d = pyfits.open(url, cache=False)
                  didurl=True
                  util.record_upload(re.sub(os.path.expandvars(atomdbroot),'',fname))
                except urllib2.URLError:
                  print "Error trying to open file %s. Not found locally or on"%(fname)+\
                      " server."
                  d=False
#          datacache['data']['misc'][ftype.upper()] = False



    else:
      if ftype.upper() in datacache['data'][Z][z1].keys():
      # this means we have the data cached, no need to fetch it
        pass
      else:
      # check for file location overrides
        if settings:
          if 'filemap' in settings.keys():
            if settings['filemap']:
              fmapfile = settings['filemap']
          if 'atomdbroot' in settings.keys():
            if settings['atomdbroot']:
              atomdbroot = settings['atomdbroot']

        fname = get_filemap_file(ftype, Z, z1, fmapfile=fmapfile,\
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
                url = re.sub(os.path.expandvars(atomdbroot),\
                             'ftp://sao-ftp.harvard.edu/AtomDB',fname)+'.gz'
                try:
                  d = pyfits.open(url, cache=False)
                  didurl=True
                  util.record_upload(re.sub(os.path.expandvars(atomdbroot),'',fname))
                except urllib2.URLError:
                  print "Error trying to open file %s. Not found locally or on"%(fname)+\
                      " server."
                  d=False



  else:
    if settings:
      if settings['filemap']:
        fmapfile = settings['filemap']
      if 'atomdbroot' in settings.keys():
        if settings['atomdbroot']:
          atomdbroot = settings['atomdbroot']

    if ismisc:

      fname = get_filemap_file(ftype, False, False, \
                             quiet=True, fmapfile=fmapfile,\
                             atomdbroot=atomdbroot, misc = True)
    else:
      fname = get_filemap_file(ftype, Z, z1, \
                             quiet=True, fmapfile=fmapfile,\
                             atomdbroot=atomdbroot)







#    if ftype.lower()=='ionbal':
#      fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v2.0.2_ionbal.fits'

    if ftype.lower()=='ionbal':
      # conversion here:
      curversion = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]

      if curversion in ['2.0.0', '2.0.1', '2.0.2','3.0.0','3.0.1','3.0.2','3.0.3']:
        fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v2.0.2_ionbal.fits'
      elif curversion in ['3.0.4','3.0.5','3.0.6']:
        fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v3.0.4_ionbal.fits'
      else:
        fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/v3.0.7_ionbal.fits'

    elif ftype.lower()=='eigen':
      # conversion here:
      curversion = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]

      if curversion in ['2.0.0', '2.0.1', '2.0.2','3.0.0','3.0.1','3.0.2','3.0.3']:
        fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/eigen/eigen%s_v3.0.fits'%(atomic.Ztoelsymb(Z).lower())
      elif curversion in ['3.0.4','3.0.5','3.0.6']:
        fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/eigen/eigen%s_v3.0.4.fits'%(atomic.Ztoelsymb(Z).lower())
      else:
        fname = os.path.expandvars(atomdbroot)+'/APED/ionbal/eigen/eigen%s_v3.0.7.fits'%(atomic.Ztoelsymb(Z).lower())


    if fname=='':
          # This is expected if it's an ionbal file
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
            url = re.sub(os.path.expandvars(atomdbroot),\
                         'ftp://sao-ftp.harvard.edu/AtomDB',fname)+'.gz'
            print "trying URL %s"%(url)
            try:
              d = pyfits.open(url, cache=False)
              didurl=True
              util.record_upload(re.sub(os.path.expandvars(atomdbroot),'',fname))
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
        try:
          d[1].columns.change_name('COEFF_OM','EFFCOLLSTRPAR')
        except:
          d[1].columns[d[1].data.names.index('COEFF_OM')].name='EFFCOLLSTRPAR'

    if datacache:
      if ismisc:
        if ftype.upper()=='EIGEN':
          if not Z in datacache['data']['misc']['EIGEN'].keys():
            datacache['data']['misc']['EIGEN'][Z]=d
            datacache['datasums']['misc']['EIGEN'][Z] = d[1].header['DATASUM']
            return datacache['data']['misc'][ftype.upper()][Z]

        else:
          datacache['data']['misc'][ftype.upper()] = d
          datacache['datasums']['misc'][ftype.upper()] = d[1].header['DATASUM']
          return datacache['data']['misc'][ftype.upper()]

      else:
        datacache['data'][Z][z1][ftype.upper()] = d
        datacache['datasums'][Z][z1][ftype.upper()] = d[1].header['DATASUM']

        return datacache['data'][Z][z1][ftype.upper()]
    else:
      return d
  else:
    if datacache:
      if ismisc:
        if ftype.upper()=='EIGEN':
          return datacache['data']['misc'][ftype.upper()][Z]
        else:
          return datacache['data']['misc'][ftype.upper()]
      else:
        return datacache['data'][Z][z1][ftype.upper()]


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
  i = numpy.where((pidat[1].data['lev_init'] ==lev_init) &\
                  (pidat[1].data['lev_final'] ==lev_final))[0]
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
  """
  Returns RRC in photons cm3 s-1 keV-1

  Parameters
  ----------
  E:
  Z: int
    Atomic number of element (i.e. 8 for Oxygen)
  z1: int
    Ion charge +1 e.g. 5 for C+4, a.k.a. C V
  rrc_ph_factor: float
    Conversion factor for RRC.
  IonE: float
    Ionization potential of ion
  kT: float
    Temperature (keV)
  levdat: lvdat line
    Line from the lvdat file
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
  float
    The RRC in photons cm3 s-1 keV-1 at energy(ies) E.

  """



  isiter=True
  try:
    _ = (e for e in E)
  except TypeError:
    isiter=False
    E = numpy.array([E])

  igood = numpy.where(((-(E - IonE)/kT)>const.MIN_RRC_EXPONENT)&\
                       ((-(E - IonE)/kT)<const.MAX_RRC_EXPONENT))[0]

  ifinite = numpy.isfinite(E[igood]**2*numpy.exp(-(E[igood] - IonE)/kT))
  #print kT
  if sum(ifinite)!=len(igood):
    print "found %i not finite indices!"%(len(igood)-sum(ifinite))
  if levdat['PHOT_TYPE']==const.HYDROGENIC:
    levdat['PHOT_PAR'][0] = levdat['n_quan']
    levdat['PHOT_PAR'][1] = levdat['l_quan']
    if levdat['PHOT_PAR'][0]>5:
      return  numpy.zeros(len(E))
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

  Evec,isvec=util.make_vec(E)
#  isvec = True
#  try:
#    _ = (e for e in E)
#  except TypeError:
#    isvec = False
#  Evec = numpy.array([E])
#
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
      pidat = get_data(Z, z1, 'PI', settings=settings, datacache=datacache)
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
      else:
        pidat = xstardata
        initlevel = int(pi_coeffts[0])
      sigma_coeff = sort_pi_data(pidat, initlevel, xstarfinallev)
      sig_coeffts = sigma_coeff


  # now calculate the sigma.

  if pi_type==const.NO_PHOT:
    result = 0.0

  elif pi_type == const.VERNER:
    iE = numpy.where(Evec > sig_coeffts['E_th'])[0]
    if len(iE) > 0:
      y = Evec[iE]/sig_coeffts['E_0']
      F_y = ((y-1)**2. + sig_coeffts['yw']**2.) * y**(-sig_coeffts['Q']) *\
            (1.+numpy.sqrt(y/sig_coeffts['ya']))**(-sig_coeffts['P'])
      result[iE] = sig_coeffts['sigma0'] * F_y * 1e-18 # conversion to cm2

  elif pi_type == const.HYDROGENIC:
    if (sig_coeffts['nq']<6) & (sig_coeffts['nq']!=0):
      result = sigma_hydrogenic(sig_coeffts['Zel'],\
                                sig_coeffts['nq'],\
                                sig_coeffts['lq'],\
                                Evec)
#      print "Zel=%i, n=%i, l=%i"%(sig_coeffts['Zel'],\
#                                sig_coeffts['nq'],\
#                                sig_coeffts['lq'])
#      for i in range(len(Evec)):
#        print "%6i, %e %e"%(i, Evec[i], result[i])

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

    tmp2=numpy.interp(numpy.log(Evec), \
                      numpy.log(sig_coeffts['energy']),\
                      numpy.log(sig_coeffts['pi_param']),\
                      left=numpy.nan,right=numpy.inf)

    # deal with the edge cases: zero below ionization potential,
    # extrappolate as E^-3 at high energy
    inan = numpy.isnan(tmp2)
    iinf = numpy.isinf(tmp2)
    ifin = numpy.isfinite(tmp2)
    if sum(inan) > 0:
      result[inan] = 0.0
    if sum(iinf) > 0:
      result[iinf] = sig_coeffts['pi_param'][-1]* \
                                ((Evec[iinf]/\
                                  sig_coeffts['energy'][-1])**-3.0)*1e-18
    if sum(ifin) > 0:
#      print ifin
      result[ifin] = 1e-18  * numpy.exp(tmp2[ifin])

  else:
    print "Error"
  if not isvec:
    result=result[0]
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
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
def A_twoph(A,E0,E):
  """
  Convert the A value into energy distribution for 2-photon transitions

  Parameters
  ----------
  A : float
    Einstein A for transition
  E0 : float
    Energy in keV of transition
  E : array(float)
    Energies of each bin to output continuum at (keV)
  Returns
  -------
  array(float)
    Distribution of transtion rate amongst bins E (s^-1)

  References
  ----------
  From Nussbaumer & Schmutz, 1984, A+A, 138,495
  Z is the element, and E is the energy of the bin, in keV
  y is unitless, and is equal to nu/nu0 = lambda0/lambda, where
  lambda0 = 1215.7 A for hydrogen--the base wavelength of the 2s->1s
  transition.  This fit is accurate to better than 0.6% for
  0.01 < y < 0.99

  The A_norm is the A value for neutral hydrogen for this transition.
  For other transitions, we renormalize to the appropriate A value.

  This routine is used for BOTH hydrogenic and He-like two-photon
  distributions.  This is justified using the result of
  Derevianko & Johnson, 1997, Phys Rev A, 56, 1288 who show in
  Figures 5 and 2 of that paper that the difference is everywhere
  less than 10% between these two for Z=6-28 -- it is about 5% or so.
  """

  C     = 202.0 #s^-1
  alpha = 0.88
  beta  = 1.53
  gamma = 0.80
  A_norm= 8.2249 #s^-1


  y = E/E0
  i = numpy.where((y>0) & (y<1))[0]

  result = numpy.zeros(len(E), dtype=float)

  x = y[i]*(1-y[i])

  z = (4*x)**gamma

  result[i] = C*(x*(1-z) + alpha* (x**beta)*z)

  result *= (A/A_norm)  # Also need R_Z/R_H, but even for Z=26, this is
                         # only 1.0005, so we'll ignore it.

  return result   # in s^-1

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def calc_two_phot(wavelength, einstein_a, lev_pop, ebins):
  """
  Calculate two photon spectrum

  Parameters
  ----------
  wavelength : float
    Wavelength of the transition (Angstroms)
  einstein_a : float
    The Einstein_A paramater for the transition
  lev_pop : float
    The level population for the upper level
  ebins : array(float)
    The bin edges for the spectrum (in keV)

  Returns
  -------
  array(float)
    The flux in photons cm-3 s-1 bin-1
    array is one element shorter than ebins.
  """
  emission = numpy.zeros(len(ebins)-1, dtype=float)
  E = (ebins[1:]+ebins[:-1])/2.0
  E0 = const.HC_IN_KEV_A/wavelength
  dE= ebins[1:]-ebins[:-1]
  A_E=A_twoph(einstein_a,E0,E)

#    emission = ldat['lev_pop']*(E/E0)*A_E*dE*const.ERG_KEV
  emission = lev_pop*(E/E0)*A_E*dE*const.ERG_KEV
  emission[E>=E0]=0.0
    # at this point, emission is in erg cm^3/s/bin

    # convert to photons cm^3/s/bin
  emission/= (const.ERG_KEV*E)

#  print "A=%e, E0=%e, lev_pop=%e, wavelength=%e"%(einstein_a, E0, lev_pop, wavelength)

  i = numpy.where(E<E0)[0]
#  for ii in i:
#    print "iBin: %i E: %e A_E: %e Emiss: %e"%(ii, E[ii], A_E[ii], emission[ii])
#  print "ENDtwo_photon"
  return emission
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

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def get_oscillator_strength(Z, z1, upperlev, lowerlev, datacache=False):


  """
  Get the oscillator strength f_{ij} of a transition

  Parameters
  ----------
  Z : int
    The atomic number of the element
  z1 : int
    The ion charge + 1 of the ion
  upperlev : int
    The upper level, indexed from 1
  lowerlev : int
    The lower level, indexed from 1
  datacache : dict
    Used for caching the data. See description in get_data

  Returns
  -------
  float
    The oscillator strength. Returns 0 if transition not found.
    If transition is not found but the inverse transition is present the
    oscillator strength is calculated for this instead.

  """
#
# Version 0.1 Initial Release
# Adam Foster 28 Jul 2016
#


  lvdat = get_data(Z, z1, 'LV', datacache=datacache)
  ladat = get_data(Z, z1, 'LA', datacache=datacache)

  # get the various things

  i = numpy.where((ladat[1].data['upper_lev']==upperlev) &\
                  (ladat[1].data['lower_lev']==lowerlev))[0]

  if len(i) ==0:
    i = numpy.where((ladat[1].data['upper_lev']==lowerlev) &\
                    (ladat[1].data['lower_lev']==upperlev))[0]

    if len(i) > 0:
      print "WARNING: no transition information found for transition %i->%i"%\
       (upperlev, lowerlev),
      print " but found data for reverse transition. Using this instead."
      up = lowerlev
      lo = upperlev
    else:
      print "WARNING: no transition information found for transition %i->%i"%\
       (upperlev, lowerlev),
      return 0.0
  else:
    up = upperlev
    lo = lowerlev

  i = i[0]
  Aji = ladat[1].data['einstein_a'][i]
  g_j = lvdat[1].data['lev_deg'][up-1]
  g_i = lvdat[1].data['lev_deg'][lo-1]
  lam = ladat[1].data['wavelen'][i]


  f_ij = Aji * (g_j*1.0/g_i) * (lam**2/6.6702e15)

  return f_ij


def make_lorentz(version = False, do_all=True, cie=False, power=False,\
                 stronglines=False, neicsd=False, neilines=False,\
                 neicont=False, levpop=False):
  """
  This makes all the Lorentz data comparison files from the Astrophysical
  Collisional Plasma Test Suite, version 0.4.0

  Parameters
  ----------
  version : string (optional)
    e.g. "3.0.7" to run the suite for v3.0.7. Otherwise uses latest version.

  Returns
  -------
  none
  """

  # set the version
  if util.keyword_check(version):
    util.switch_version(version)
  else:
    version = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]

  # if one (or more) in particular is specified, turn off "do_all"
  if sum([cie, power, stronglines, neicsd,neilines, neicont, levpop])>0:
    do_all=False

  # if do_all, turn them all on
  if do_all:
    cie=True
    power=True
    stronglines=True
    neicsd=True
    neilines=True
    neicont = True
    levpop = True

  # run the data
  if cie:
    lorentz_cie(version)
  if power:
    lorentz_power(version)
  if stronglines:
    lorentz_stronglines(version)
  if neicsd:
    lorentz_neicsd(version)
  if neilines:
    lorentz_neilines(version)
  if neicont:
    lorentz_neicont(version)
  if levpop:
    lorentz_levpop(version)

def lorentz_cie(version):
  """
  Calculate the CSD of equilibrium plasmas at 1e6, 6e6K and 4keV.

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """

  # open the output files
  f = open('CIE-CSD_atomdb_%s.dat'%(version),'w')
  f.write('Z   Ion   1e6K   6e6K   4.642e7K\n')
  Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]
  ionbal_out = []
  for kT in [1e6, 6e6, 4/const.KBOLTZ]:
    ionbal = apec.calc_full_ionbal(kT, extrap=True, cie=True, \
                      Zlist=Zlist)
    ionbal_out.append(ionbal)

  for Z in Zlist:
    for z1 in range(Z+1):
      f.write("%2i %2i %10.6f %10.6f %10.6f\n"%\
              (Z, z1, \
               max(-20, numpy.log10(ionbal_out[0][Z][z1])),\
               max(-20, numpy.log10(ionbal_out[1][Z][z1])),\
               max(-20, numpy.log10(ionbal_out[2][Z][z1]))))
  f.close()



def lorentz_power(version):
  """
  Calculate the power emitted from 13.6eV to 13.6keV in a 1m^3 slab of
  plasma with n_e=1e6m^-3.

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """
  Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]
  ag89=get_abundance(abundset='AG89')
  lodd=get_abundance(abundset='Lodd09')

  # open the output files
  f = open('Power_atomdb_%s.dat'%(version),'w')
  ebins = numpy.linspace(0.0136, 13.6, 10001)
  energy = const.ERG_KEV*(ebins[1:]+ebins[:-1])/2
  s = "Z "
  for iT in range(51):
    flt = 4.0+(0.1*iT)
    s += "        %.1f"%(flt)

  f.write("%s\n"%(s))

  ses = spectrum.Session(elements=Zlist)
  ses.set_specbins(ebins, specunits='keV')
  for iT in range(51):
    kT = 4.0+(iT*0.1)
    kT = 10**kT
    if iT == 0:
      kT+=1
    ses.return_spectra(kT, teunit='K', nearest=True)
    print iT

  for Z in Zlist:
    tot_e = numpy.zeros(51)
    s = "%2i"%(Z)
    print Z
    for iT in range(2,53):
      spec = ses.spectra[iT].spectrum_by_Z[Z]
      e = spec*energy
      # add corrections for NH != 1, to 1m3 volume, and abundance set
      # from AG89 to Lodders 2009
      tot_e[iT-2] = sum(e)*0.8365*1e6*lodd[Z]/ag89[Z]
      s+= " %10.6f"%(numpy.log10(tot_e[iT-2]))
    f.write("%s\n"%(s))
  f.close()


def lorentz_stronglines(version):
  """
  Calculate the 100 strongest lines below 1000A from a 1m3 slab of plasma
  with n_e = 1e6m-3, at 3 different temperatures: 10^6K, 6e6K, 4.642e7K

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """

  ldat = pyfits.open(os.path.expandvars("$ATOMDB/apec_v%s_line.fits"%(version)))

  # get the strongest 500 lines
  llist_type = numpy.dtype({'names':['Lambda','Z','Ion','Flux', 'UpperLev', 'LowerLev'],\
                              'formats':[float, int, int, float, int, int]})
  lodd=get_abundance(abundset='Lodd09')
  ag89=get_abundance(abundset='AG89')
  for ikT, kT in enumerate([1e6, 6e6, 4/const.KBOLTZ]):
    upind = numpy.where(ldat[1].data['kT']>kT*const.KBOLTZ)[0][0]+2
    loind = upind-1

    llist_lo = ldat[loind].data
    llist_up = ldat[upind].data
    # only include lines below 1000A
    llist_lo = llist_lo[llist_lo['Lambda']<1000]
    llist_up = llist_up[llist_up['Lambda']<1000]
    fname = 'StrongLines%i_atomdb_%s.dat'%(ikT+1, version)
    f = open(fname, 'w')
    # now interpolate
    t1 = numpy.log(ldat[1].data['kT'][loind-2])
    t2 = numpy.log(ldat[1].data['kT'][upind-2])



    r1 = 1- (numpy.log(kT*const.KBOLTZ)-t1)/(t2-t1)
    r2 = 1- r1


    llist_lo.sort(order='Epsilon')
    llist_lo=llist_lo[::-1]

    llist_up.sort(order='Epsilon')
    llist_up=llist_up[::-1]

    llist_out = numpy.zeros(500, dtype = llist_type)
    llist_out['Lambda'] = llist_lo['Lambda'][:500]
    llist_out['Z'] = llist_lo['Element'][:500]
    llist_out['Ion'] = llist_lo['Ion'][:500]-1
    llist_out['Flux'] = numpy.log(llist_lo['Epsilon'][:500])*r1
    llist_out['UpperLev'] = llist_lo['UpperLev'][:500]
    llist_out['LowerLev'] = llist_lo['LowerLev'][:500]

    for i in range(len(llist_out)):
      ii = numpy.where((llist_up['Element']==llist_out['Z'][i]) &\
                       (llist_up['Ion'] == llist_out['Ion'][i]+1)&\
                       (llist_up['UpperLev'] == llist_out['UpperLev'][i])&\
                       (llist_up['LowerLev'] == llist_out['LowerLev'][i]))[0]
      if len(ii) > 0:
        llist_out['Flux'][i] += numpy.log(llist_up['Epsilon'][ii[0]])*r2

    # do abundance
    abvec = numpy.zeros(max(lodd.keys())+1)
    for Z in range(1,len(abvec)):
      abvec[Z] = lodd[Z]/ag89[Z]
    ab = abvec[llist_out['Z']]

    llist_out['Flux']=numpy.exp(llist_out['Flux'])
    llist_out['Flux']*=ab
    llist_out.sort(order='Flux')
    llist_out = llist_out[::-1]

    now = datetime.datetime.now()

    f.write('# Generated %s by %s\n'%(now.strftime("%c"), os.getenv('USER')))
    f.write('# 100 strongest lines in collisional plasma with T=%eK\n'%(kT))
    f.write('Indx Lambda Z Ion Flux\n')

    for i in range(100):
      print "%3i %14f %2i %2i %10.6f"%\
            (i+1,\
             llist_out['Lambda'][i],\
             llist_out['Z'][i],\
             llist_out['Ion'][i],\
             numpy.log10(llist_out['Flux'][i]*0.8365*1e6))

      f.write("%3i %14f %2i %2i %10.6f\n"%\
            (i+1,\
             llist_out['Lambda'][i],\
             llist_out['Z'][i],\
             llist_out['Ion'][i],\
             numpy.log10(llist_out['Flux'][i]*0.8365*1e6)))
    f.close()

def lorentz_neicsd(version):
  """
  Charge state distribution of a gas ionizing from 1e4K to 2.321e7K (=2keV)
  at a fluence ($n_e$ * t, or $\tau$) of $10^{10}$ cm$^-3$ s

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """
  Te_init = 1e4
  Te_final = 2.321e7
  tau = 1e10
  now = datetime.datetime.now()
  util.switch_version(version)
  f = open('NEI-CSD_atomdb_v%s.dat'%(version),'w')
  f.write('# Generated %s by %s\n'%(now.strftime("%c"), os.getenv('USER')))
  f.write('# CSD of plasma ionizing from 1e4K to 2.321e7K for ne*t = 1e10 cm^-3 s\n')
  f.write('Z Ion Pop\n')
  for Z in [1,2,6,7,8,10,12,14,16,18,20,26,28]:
    print Z
    ionbal = apec.solve_ionbal_eigen(Z, Te_final,  tau=tau, Te_init=Te_init, \
                         teunit='K')

    for i in range(len(ionbal)):
      f.write('%2i %2i %e\n'%(Z,i,ionbal[i]))
  f.close()


def lorentz_neilines(version):
  """
  100 strongest lines with wavelength < 1000A for a 1cm^3 plasma
  (1) starting at 1e4K, going to 2.321e7K
  at a fluence ($n_e$ * t, or $\tau$) of $10^{10}$ cm$^-3$ s
  (2) starting at 3.5keV, going to 1.5keV
  at a fluence ($n_e$ * t, or $\tau$) of $10^{10}$ cm$^-3$ s

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """
  import pickle
  lodd=get_abundance(abundset='Lodd09')
  ag89=get_abundance(abundset='AG89')

  abund = numpy.zeros(len(ag89)+1)
  for i in lodd.keys():
    abund[i]=lodd[i]/ag89[i]
  Te_init_list = [1e4, 3.5/const.KBOLTZ]
  Te_final_list = [2.321e7, 1.5/const.KBOLTZ]
  tau_list = [1e10, 1e10]
  now = datetime.datetime.now()
  util.switch_version(version)



  for irun in range(len(Te_init_list)):


    Te_init = Te_init_list[irun]
    Te_final = Te_final_list[irun]
    tau = tau_list[irun]


    f = open('NEI-Lines%i_atomdb_v%s.dat'%(irun+1,version),'w')
    f.write('# Generated %s by %s\n'%(now.strftime("%c"), os.getenv('USER')))
    f.write('# Stong Lines of 1m3 slab of plasma ionizing from %.3eK to %.3eK for ne*t = %.0e cm^-3 s\n'%\
             (Te_init, Te_final, tau))
    f.write('# Lambda is in Angstroms\n')
    f.write('Lambda Z Ion Flux\n')

    ldat = pyfits.open(os.path.expandvars("$ATOMDB/apec_v%s_nei_line.fits"%(version)))

    upind = numpy.where(ldat[1].data['kT']>Te_final*const.KBOLTZ)[0][0]+2
    loind = upind-1
    print "upind = %i, %eK" %(upind, ldat[upind].header['TEMPERATURE'])
    print "loind = %i, %eK"%(loind, ldat[loind].header['TEMPERATURE'])
    print "Te_final = %eK"%(Te_final)


    t1 = numpy.log(ldat[1].data['kT'][loind-2])
    t2 = numpy.log(ldat[1].data['kT'][upind-2])
    print "t1 = ", t1
    print "t2 = ", t2
    print "log(tefinal)", numpy.log(Te_final*const.KBOLTZ)
    r1 = 1- (numpy.log(Te_final*const.KBOLTZ)-t1)/(t2-t1)
    r2 = 1- r1
    print "r1= ",r1, "r2  ", r2



    llist_lo = ldat[loind].data
    llist_up = ldat[upind].data
    # only include lines below 1000A
    llist_lo = llist_lo[llist_lo['Lambda']<1000]
    llist_up = llist_up[llist_up['Lambda']<1000]

    ionbal = {}
    for Z in range(1,31):
      ionbal[Z] = numpy.zeros(Z+1)

    for Z in [1,2,6,7,8,10,12,14,16,18,20,26,28]:

      ionbal[Z] = apec.solve_ionbal_eigen(Z, Te_final,  tau=tau, Te_init=Te_init, \
                           teunit='K')

    ionbal_square = numpy.zeros([31,31])
    for Z in [1,2,6,7,8,10,12,14,16,18,20,26,28]:
      ionbal_square[Z,:Z+1] = ionbal[Z]

    scale = abund[llist_lo['ELEMENT']] * ionbal_square[llist_lo['ELEMENT'],llist_lo['ION_DRV']-1]
    llist_lo['EPSILON'] *= scale

    print "max scale:", max(scale)
    print "max eps:", max(llist_lo['EPSILON'])

    scale = abund[llist_up['ELEMENT']] * ionbal_square[llist_up['ELEMENT'],llist_up['ION_DRV']-1]
    llist_up['EPSILON'] *= scale
    print "max scale:", max(scale)
    print "max eps:", max(llist_up['EPSILON'])
      # remove weak lines
    llist_lo = llist_lo[llist_lo['EPSILON']>1e-30]
    llist_up = llist_up[llist_up['EPSILON']>1e-30]
    print "len lo = ", len(llist_lo)
    print "len up = ", len(llist_up)


#    llist_lo = numpy.array(llist_lo)

    llist_lo.sort(order = ['Element','Ion','UpperLev','LowerLev', 'Ion_drv'])
    llist_out_dtype=numpy.dtype({'names':['Lambda','Element','Ion','UpperLev','LowerLev','Epsilon'],\
                                 'formats':[float, int, int, int, int, float]})
    llist_lo_out = numpy.zeros(len(llist_lo), dtype=llist_out_dtype)
    iline = -1
    print "summing NEI lines for lo ind"
    for i in range(len(llist_lo)):

      if ( (llist_lo['Element'][i] == llist_lo_out['Element'][iline]) &\
           (llist_lo['UpperLev'][i] == llist_lo_out['UpperLev'][iline]) &\
           (llist_lo['LowerLev'][i] == llist_lo_out['LowerLev'][iline]) &\
           (llist_lo['Ion'][i] == llist_lo_out['Ion'][iline])):
        llist_lo_out['Epsilon'][iline]+=llist_lo['Epsilon'][i]

      else:
        iline+=1
        llist_lo_out['Element'][iline] = llist_lo['Element'][i]
        llist_lo_out['Ion'][iline] = llist_lo['Ion'][i]
        llist_lo_out['UpperLev'][iline] = llist_lo['UpperLev'][i]
        llist_lo_out['LowerLev'][iline] = llist_lo['LowerLev'][i]
        llist_lo_out['Epsilon'][iline] = llist_lo['Epsilon'][i]
        llist_lo_out['Lambda'][iline] = llist_lo['Lambda'][i]

    llist_lo_out=llist_lo_out[:iline+1]
    llist_lo_out=llist_lo_out[llist_lo_out['Epsilon']>1e-20]
    llist_lo_out.sort(order='Epsilon')

    print "summing NEI lines for up ind"

    llist_up.sort(order = ['Element','Ion','UpperLev','LowerLev', 'Ion_drv'])
    llist_up_out = numpy.zeros(len(llist_up), dtype=llist_out_dtype)
    iline = -1

    for i in range(len(llist_up)):

      if ( (llist_up['Element'][i] == llist_up_out['Element'][iline]) &\
           (llist_up['UpperLev'][i] == llist_up_out['UpperLev'][iline]) &\
           (llist_up['LowerLev'][i] == llist_up_out['LowerLev'][iline]) &\
           (llist_up['Ion'][i] == llist_up_out['Ion'][iline])):
        llist_up_out['Epsilon'][iline]+=llist_up['Epsilon'][i]

      else:
        iline+=1
        llist_up_out['Element'][iline] = llist_up['Element'][i]
        llist_up_out['Ion'][iline] = llist_up['Ion'][i]
        llist_up_out['UpperLev'][iline] = llist_up['UpperLev'][i]
        llist_up_out['LowerLev'][iline] = llist_up['LowerLev'][i]
        llist_up_out['Epsilon'][iline] = llist_up['Epsilon'][i]
        llist_up_out['Lambda'][iline] = llist_up['Lambda'][i]

    llist_up_out=llist_up_out[:iline+1]
    llist_up_out=llist_up_out[llist_up_out['Epsilon']>1e-20]
    llist_up_out.sort(order='Epsilon')

    print "Combining emissivities"

    # trim to 1000 strongest lines
    llist_up_out = llist_up_out[-1000:]
    llist_lo_out = llist_lo_out[-1000:]

    print 'saving llist_up_out to llist_up_out.pkl'
    pickle.dump(llist_up_out, open('llist_up_out_%i.pkl'%(irun),'w'))
    print 'saving llist_lo_out to llist_lo_out.pkl'
    pickle.dump(llist_lo_out, open('llist_lo_out_%i.pkl'%(irun),'w'))



    llist_out = numpy.zeros(2000,dtype=llist_out_dtype)
    llist_out[:1000] = llist_lo_out
    llist_out['Epsilon']*=r1

    iline = 1000
    for i in range(len(llist_up_out)):
      j = numpy.where( (llist_out['Element'] == llist_up_out['Element'][i]) &\
                       (llist_out['Ion'] == llist_up_out['Ion'][i]) &\
                       (llist_out['UpperLev'] == llist_up_out['UpperLev'][i]) &\
                       (llist_out['LowerLev'] == llist_up_out['LowerLev'][i]))[0]
      if len(j) == 0:
        llist_out['Element'][iline] = llist_up_out['Element'][i]
        llist_out['Lambda'][iline] = llist_up_out['Lambda'][i]
        llist_out['Ion'][iline] = llist_up_out['Ion'][i]
        llist_out['UpperLev'][iline] = llist_up_out['UpperLev'][i]
        llist_out['LowerLev'][iline] = llist_up_out['LowerLev'][i]
        llist_out['Epsilon'][iline] = llist_up_out['Epsilon'][i]*r2
        iline+=1
      else:
        print "adding epsilon=%e to %e line"%(llist_up_out['Epsilon'][i]*r2,\
                llist_out['Epsilon'][j[0]])
        llist_out['Epsilon'][j[0]] += llist_up_out['Epsilon'][i]*r2
        print llist_out[j[0]]
    print "Final filtering, keep %i lines"%(iline)
    llist_out = llist_out[:iline]

    llist_out.sort(order='Epsilon')

    llist_out=llist_out[::-1]

    llist_out = llist_out[:100]

    # scale for NH < 1
    llist_out['Epsilon']*=0.8365*1e6

    for i in range(len(llist_out)):
      print llist_out[i]
      f.write('%9.5f %2i %2i %e\n'%\
              (llist_out['Lambda'][i],\
               llist_out['Element'][i],\
               llist_out['Ion'][i]-1,\
               llist_out['Epsilon'][i]))

    f.close()




def lorentz_neicont(version):
  """
  Full spectrum of a gas ionizing from 1e4K to 2.321e7K (=2keV)
  at a fluence ($n_e$ * t, or $\tau$) of $10^{10}$ cm$^-3$ s

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """
  Te_init = 3.5/(const.KBOLTZ)
  Te_final = 1.5/(const.KBOLTZ)
  tau = 1e10
  # set up 1 ev bins from 10eV to 10keV
  ebins = numpy.arange(0.01, 10.001, 0.001)
  now = datetime.datetime.now()
  util.switch_version(version)

  # make the spectrum.
  speclo = numpy.zeros(len(ebins)-1)
  specup = numpy.zeros(len(ebins)-1)
  ag89 = get_abundance(abundset='AG89')
  lodd = get_abundance(abundset='Lodd09')
  ldat = pyfits.open(os.path.expandvars("$ATOMDB/apec_nei_line.fits"))
  cdat = pyfits.open(os.path.expandvars("$ATOMDB/apec_nei_comp.fits"))


  upind = numpy.where(ldat[1].data['kT']>Te_final*const.KBOLTZ)[0][0]+2
  loind = upind-1

  for Z in range(1,31):
    print "starting element %s"%(atomic.Ztoelname(Z))
    ionbal = apec.solve_ionbal_eigen(Z, Te_final,  tau=tau, Te_init=Te_init, \
                           teunit='K')

    abund = lodd[Z]/ag89[Z]
    for z in range(len(ionbal)):
      z1 = z+1
      if ionbal[z] > 1e-10:
        tmp = spectrum.make_ion_spectrum(ebins, loind, Z, z1, linefile=ldat,\
                                         cocofile=cdat)
        speclo+=tmp*ionbal[z]*abund

        tmp = spectrum.make_ion_spectrum(ebins, upind, Z, z1, linefile=ldat,\
                                         cocofile=cdat)
        specup+=tmp*ionbal[z]*abund

  # now interpolate
  t1 = numpy.log(ldat[1].data['kT'][loind-2])
  t2 = numpy.log(ldat[1].data['kT'][upind-2])
  print "t1 = ", t1
  print "t2 = ", t2
  print "log(tefinal)", numpy.log(Te_final*const.KBOLTZ)
  r1 = 1- (numpy.log(Te_final*const.KBOLTZ)-t1)/(t2-t1)
  r2 = 1- r1
  print "r1= ",r1, "r2  ", r2
  spec = speclo*r1+specup*r2

  # now scale spectrum by NH to get correct norm, and 1e6 to get to 1m^3
  spec *= 0.8365*1e6

  f = open('NEI-Cont_atomdb_%s.dat'%(version),'w')
  f.write('# Generated %s by %s\n'%(now.strftime("%c"), os.getenv('USER')))
  f.write('# Full spectrum of 1m3 plasma slab with n_e=1e6m^-3\n')
  f.write('# recombining from %.3eK to %.3eK for ne*t = %.0e cm^-3 s\n'%\
             (Te_init, Te_final, tau))
  f.write('# listed on 1eV bins, no broadening applied\n')
  f.write('# flux in ph cm^3 s^-1 bin^-1\n')
  f.write('# BinLo BinHi Flux\n')
  for i in range(len(spec)):
    f.write("%6.3f %6.3f %e\n"%(ebins[i], ebins[i+1], spec[i]))
  f.close()


def get_lorentz_levpop(Z,z1,up,lo, Te, Ne, version, linelabel):
  """
  calculate the level population for a particular ion

  """
  abund = get_abundance(abundset='Lodd09')[Z]
  util.switch_version(version)
  # first, get the ionization balance
  datacache={}
  lvdat = get_data(Z,z1,'LV', datacache=datacache)
  ionbal = apec.solve_ionbal_eigen(Z,Te, datacache=datacache)
  settings = apec.parse_par_file(os.path.expandvars('$ATOMDB/apec_v%s.par'%\
                                (version)))

  # now find the total rate out of the level
  la_rates_up, la_rates_lo, la_rates_rates =\
                          apec.gather_rates(Z, z1, Te, Ne, datacache=datacache, \
                          settings=settings, do_la=True, do_ai=False, \
                          do_ec=False, do_pc=False,  do_ir=False)

  ec_rates_up,ec_rates_lo,ec_rates_rates = \
                          apec.gather_rates(Z, z1, Te, Ne, datacache=datacache, \
                          settings=settings, do_la=False, do_ai=False, \
                          do_ec=True, do_pc=False,  do_ir=False)

  i = numpy.where((la_rates_up==up-1) &\
                  (la_rates_lo != up-1))[0]
  print la_rates_rates[i]
  out = sum(la_rates_rates[i])
  print "calc out from LA file %e " %(out)
  if 'ARAD_TOT' in lvdat[1].data.names:
    out = lvdat[1].data['ARAD_TOT'][up-1]
    print "calc out from LV file %e " %(out)
  i = numpy.where((ec_rates_up==up-1) &\
                  (ec_rates_lo != up-1))[0]
  print "calc out from EC file %e " %(sum(ec_rates_rates[i]))
  out += sum(ec_rates_rates[i])


  # ionization
  if (z1 > 1):
    ioniz_ionpop = ionbal[z1-2]

    ionup, ionlo, ionrates = apec.gather_rates(Z, z1-1, Te, Ne, \
                           datacache=datacache, settings=settings)
    lev_pop = apec.solve_level_pop(ionup,ionlo,ionrates, settings)

    # scale by ion pop, abundance
    lev_pop *= abund*ioniz_ionpop

    for i in range(len(lev_pop)):
      print i, lev_pop[i]

    lvdat = get_data(Z,z1-1,'LV', datacache=datacache)
    if 'AAUT_TOT' in lvdat[1].data.names:
      iaut = numpy.where(lvdat[1].data['AAUT_TOT']>0)[0]
      print iaut
    aidat = get_data(Z,z1-1,'AI', datacache=datacache)
    lvdat2 = get_data(Z,z1,'LV', datacache=datacache)

#    for i in range(len(lvdat2[1].data)):
#      k = numpy.where(aidat[1].data['LEVEL_FINAL']==i+1)[0]
#      if len(k) > 0:
#        print i+1, sum(aidat[1].data['AUTO_RATE'][k]*
#                       lev_pop[aidat[1].data['LEVEL_INIT'][k]-1])
#
#    zzz=raw_input()

    # now find the ionization rates
    do_xi = True
    lev_pop_xi=apec.calc_ioniz_popn(lev_pop, Z, z1, z1-1, Te, Ne, \
                            settings=settings, datacache=datacache, \
                             do_xi=True)
    lev_pop_noxi=apec.calc_ioniz_popn(lev_pop, Z, z1, z1-1, Te, Ne, \
                            settings=settings, datacache=datacache, \
                             do_xi=False)
    lev_pop_ci = lev_pop_xi-lev_pop_noxi

    excitauto = lev_pop_noxi[up-1]*out
    direction = lev_pop_ci[up-1]*out

  else:
    ionlevpop = numpy.zeros(len(lvdat[1].data))
    excitauto = 0.0
    direction = 0.0


# excitation
  if True:
    ionpop = ionbal[z1-1]
    lvdat=get_data(Z,z1,'LV', datacache=datacache)

    excup, exclo, excrates = apec.gather_rates(Z, z1, Te, Ne, datacache=datacache,\
                                      settings=settings)
    lev_pop = apec.solve_level_pop(excup,exclo,excrates, settings)

    # scale by ion pop, abundance
    lev_pop *= abund*ionpop


    # now gather rates by individual process:
    ecup, eclo, ecrates = apec.gather_rates(Z, z1, Te, Ne, datacache=datacache,\
                                      settings=settings, do_la=False, \
                                      do_ai=False, \
                          do_ec=True, do_pc=False,  do_ir=False)

    k = numpy.where((ecup==0) & (eclo==up-1))[0]
    print ecrates[k]

    laup, lalo, larates = apec.gather_rates(Z, z1, Te, Ne, datacache=datacache,\
                                      settings=settings, do_la=True, \
                                      do_ai=False, \
                          do_ec=False, do_pc=False,  do_ir=False)
    pcup, pclo, pcrates = apec.gather_rates(Z, z1, Te, Ne, datacache=datacache,\
                                      settings=settings, do_la=False, \
                                      do_ai=False, \
                          do_ec=False, do_pc=True,  do_ir=False)
    # now multiply these rates by their driving populations to get the
    # contribution

    # electron excitation and de-excitation
    iec= numpy.where(ecup == up-1)[0]
    eexcout = 0.0
    edexout = 0.0

    eexcin = 0.0
    edexin = 0.0
    print ecup[iec]
    print eclo[iec]
    print len(lvdat[1].data)
    print "HMM"
    for ii in iec:
      if eclo[ii]==ecup[ii]: continue
      if lvdat[1].data['ENERGY'][ecup[ii]] < lvdat[1].data['ENERGY'][eclo[ii]]:
        eexcout += ecrates[ii]*lev_pop[ecup[ii]]
      else:
        edexout += ecrates[ii]*lev_pop[ecup[ii]]

    iec= numpy.where(eclo == up-1)[0]
    for ii in iec:
      if eclo[ii]==ecup[ii]: continue
      if lvdat[1].data['ENERGY'][ecup[ii]] < lvdat[1].data['ENERGY'][eclo[ii]]:
        eexcin += ecrates[ii]*lev_pop[ecup[ii]]
      else:
        edexin += ecrates[ii]*lev_pop[ecup[ii]]


    # same for proton collisions
    ipc= numpy.where(pcup == up-1)[0]
    pexcin = 0.0
    pdexin = 0.0
    pexcout = 0.0
    pdexout = 0.0

    for ii in ipc:
      if pclo[ii]==pcup[ii]: continue
      if lvdat[1].data['ENERGY'][pcup[ii]] < lvdat[1].data['ENERGY'][pclo[ii]]:
        pexcout += pcrates[ii]*lev_pop[pcup[ii]]
      else:
        pdexout += pcrates[ii]*lev_pop[pcup[ii]]

    ipc= numpy.where(pclo == up-1)[0]
    for ii in ipc:
      if pclo[ii]==pcup[ii]: continue
      if lvdat[1].data['ENERGY'][pcup[ii]] < lvdat[1].data['ENERGY'][pclo[ii]]:
        pexcin += pcrates[ii]*lev_pop[pcup[ii]]
      else:
        pdexin += pcrates[ii]*lev_pop[pcup[ii]]

    # something something radiative transitions
    ila= numpy.where((lalo == up-1)& (lalo!=laup))[0]

    lain = sum(larates[ila] * lev_pop[laup[ila]])

    for iila in ila:
      if larates[iila] * lev_pop[laup[iila]] > 1e-17:
        print larates[iila], laup[iila], lalo[iila], lev_pop[laup[iila]], larates[iila] * lev_pop[laup[iila]]

  else:
    lain = 0.0
    eexcin = 0.0
    eexcout = 0.0
    edexin = 0.0
    edexout = 0.0
    pexcin = 0.0
    pexcout = 0.0
    pdexin = 0.0
    pdexout = 0.0

  # now do recombination!
  if z1 < Z+1:
    # recombination is always from the groudn state so we are ok
    # no need to caculate level pop
    ionpop = ionbal[z1]
    levpop =numpy.ones(1)
    levpop *= ionpop*abund


    ebins = apec.make_vector_nbins(settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])
    # do the DR satellite lines
    if settings['DRSatellite']:
      print "Start calc_satellte run_apec_ion at %s"%(time.asctime())

      linelist_dr, drlevrates = apec.calc_satellite(Z, z1, Te, \
                                datacache=datacache, settings=settings)
      drlevrates *=ionpop*abund

    else:
      linelist_dr = numpy.zeros(0, dtype= generate_datatypes(linetype))
      drlevrates = 0.0
    print "drlevrates"
    print drlevrates
    # Radiative Recombination
    if settings['RRC']:
      rrc, rrlevrates = calc_rad_rec_cont(Z, z1, z1+1, Te, \
                        ebins, settings=settings, datacache=datacache)
      rrlevrates*=ionpop*abund

    else:
      rrlevrates=0.0
    print "rrlevrates"
    print rrlevrates

    # if there is recombination to process:
    tmpdrlevrates,xxx = util.make_vec(drlevrates)
    tmprrlevrates,xxx = util.make_vec(rrlevrates)

    sum_rr_out=0.0
    sum_dr_out = 0.0
    if sum(tmpdrlevrates) + sum(tmprrlevrates)>0:
      print "Start calc_recomb_popn at %s"%(time.asctime())

      levpop_dr=apec.calc_recomb_popn(levpop, Z, z1,\
                                      z1+1, Te, Ne, drlevrates,\
                                      rrlevrates,\
                                      datacache=datacache, \
                                      settings=settings,
                                      dronly=True)

      levpop_rr=apec.calc_recomb_popn(levpop, Z, z1,\
                                      z1+1, Te, Ne, drlevrates,\
                                      rrlevrates,\
                                      datacache=datacache, \
                                      settings=settings,
                                      rronly=True)

      # make into lines
      ladat = get_data(Z,z1,'LA', datacache=datacache)

      j = numpy.where(ladat[1].data['LOWER_LEV']==lo)[0]
      print "levpop_rr, levopo_dr"
      for i in range(len(levpop_rr)):
        print i, levpop_rr[i], levpop_dr[i]

      print "levpop", levpop
      sum_rr_in = sum(ladat[1].data['EINSTEIN_A'][j] *\
                      levpop_rr[ladat[1].data['UPPER_LEV'][j]-1])

      sum_dr_in = sum(ladat[1].data['EINSTEIN_A'][j] *\
                      levpop_dr[ladat[1].data['UPPER_LEV'][j]-1])

      j = numpy.where(ladat[1].data['UPPER_LEV']==up)[0]
      print "up=", up
      print ladat[1].data['EINSTEIN_A'][j]
      sum_rr_out = sum(ladat[1].data['EINSTEIN_A'][j] *\
                      levpop_rr[up-1])

      sum_dr_out = sum(ladat[1].data['EINSTEIN_A'][j] *\
                      levpop_dr[up-1])
      rate_out = sum(ladat[1].data['EINSTEIN_A'][j])

      print sum_rr_in
      print sum_dr_in
      print sum_rr_out
      print sum_dr_out
      print rate_out















  print "ionbal", ionbal
#  print "levpop = %e"%(lev_pop[up-1])
  print "DI in = %e"%(direction)
  print "Excit-Auto in = %e"%(excitauto)
  print "Rad in = %e"%(lain)
  print "Rad out = %e" %(out)
  print "E_excite in = %e" %(eexcin)
  print "E_excite out = %e"%(eexcout)
  print "E_dexcite in = %e"%(edexin)
  print "E_dexcite out = %e"%(edexout)
  print "P_excite in = %e"%(pexcin)
  print "P_excite out = %e"%(pexcout)
  print "P_dexcite in = %e"%(pdexin)
  print "P_dexcite out = %e"%(pdexout)
  print "RR in = %e"%(sum_rr_out)
  print "DR in = %e"%(sum_dr_out)

  cfactor = 0.8365

  s = "%i %e %e %e %e %e %e %e %e %e %e %e\n"%\
      (linelabel, Te, Ne, eexcin*cfactor, edexout*cfactor, \
       pexcin * cfactor, pdexout*cfactor,\
       lain * cfactor, out, \
       sum_rr_out * cfactor, sum_dr_out*cfactor, \
       (direction+excitauto)*cfactor)

  print s
  return s

def lorentz_levpop(version):
  """
  Calculate the level populating processes for each line in the stronglines
  Files. This will require a significant rerun of APEC. Hmmmmm

  Processes to be tracked: electron excitation, electron de-excitation,
  proton excitation and dexcitation, cascade into the level, radiative out,
  recombination (incl. cascade) in, DR (incl cascade) in, and inner-shell
  ionization in (why only inner shell?)
  """

  # Method:
  # For each line, calculate the population of the level completely (i.e run most
  # of apec) from each of a range of sources (ionization, recombination,
  # excitation). Apply appropriate modifiers to get a plasma with the relevant
  # density, abundance and othe animals

  # this is just a test call for now:


  linelist = {}
  linelist['ID']=[ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17]
  linelist['Z']= [26,26,26,26,26, 8 ,8, 8, 8,26,26,26,26, 8, 8,26,26]
  linelist['z1']=[17,17,17,17,17, 7, 7, 7, 7,25,25,25,25, 8, 8,26,26]
  linelist['up']=[27,23, 5, 3, 2, 7, 6, 5, 2, 7, 6, 5, 2, 3, 4, 3, 4]
  linelist['lo']=[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

  Telist = [1e6, 6e6, 4.642e7]
  Nelist = [1e0, 1e12]
  now = datetime.datetime.now()
  f = open('LevelPop_atomdb_%s.dat'%(version),'w')
  f.write('# Generated %s by %s\n'%(now.strftime("%c"), os.getenv('USER')))
  f.write('# Populating processes for lines listed in table 2 of\n')
  f.write('# the Lorentz standards, at T=1e6, 6e6 and 2.321e7K,\n')
  f.write('# and n_e=1e6 and 1e18m^-3\n')
  f.write('# Rates are given in units of per m^3 s^-1\n')
  f.write('# we include here all direct and excitation-ionization processes\n')
  f.write('# and subsequent cascade under inner shell ionization/excitation\n')
  f.write('# columns are:\n')
  f.write('#   LineID (from Table 2)\n')
  f.write('#   Temperature (K)\n')
  f.write('#   Electron Density (cm^-3)\n')
  f.write('#   Electron excitation into level (cm^-3 s^-1)\n')
  f.write('#   Electron de-excitation out of level (cm^-3 s^-1)\n')
  f.write('#   Proton excitation into level (cm^-3 s^-1)\n')
  f.write('#   Electron excitation out of level (cm^-3 s^-1)\n')
  f.write('#   Excitation-Cascade from higher levels (cm^-3 s^-1)\n')
  f.write('#   Radiative decay out (s^-1)\n')
  f.write('#   Radiative Recombination into level (and cascade from capture to higher levels) (cm^-3 s^-1)\n')
  f.write('#   Dielectronic Recombination into level (and cascade from capture to higher levels) (cm^-3 s^-1)\n')
  f.write('#   Ionization into the level (cm^-3 s^-1)\n')

  f.write('  Num     Te    Ne    EExc    EDeExc    PExc    PDeexc    CascadeTo    RadiativeOut     RRin    DRin    ISIon\n')
  for iNe in range(len(Nelist)):
    for iTe in range(len(Telist)):
      for iline in range(len(linelist['ID'])):

        s=get_lorentz_levpop(linelist['Z'][iline],\
                           linelist['z1'][iline],\
                           linelist['up'][iline],\
                           linelist['lo'][iline],\
                           Telist[iTe],\
                           Nelist[iNe],
                           version, linelist['ID'][iline])
        f.write(s)
  f.close()
