
"""
This module contains methods for creating spectra from the AtomDB files. Some
are more primitive than others...
"""

# Adam Foster
#
# Version 0.1, July 17th 2015
# First Release
# Adam Foster

try:
  import astropy.io.fits as pyfits
except ImportError:
  import pyfits

import numpy, os, hashlib, pickle
# other pyatomdb modules
from . import atomic, util, const, atomdb, apec
from scipy.stats import norm
import time

def make_spectrum(bins, index, linefile="$ATOMDB/apec_line.fits",\
                  cocofile="$ATOMDB/apec_coco.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  elements=False, abund=False, dummyfirst=False,\
                  dolines = True, docont=True, dopseudo=True):

  """
  make_spectrum is the most generic "make me a spectrum" routine.

  It returns the emissivity in counts cm^3 s^-1 bin^-1.

  Parameters
  ----------
  bins : array(float)
       The bin edges for the spectrum to be calculated on, in \
       units of keV or Angstroms. Must be monotonically\
       increasing. Spectrum will return len(bins)-1 values.
  index : int
       The index to plot the spectrum from. note that the AtomDB files\
       the emission starts in hdu number 2. So for the first block, you\
       set index=2
  linefile : str
       The file containing all the line emission. Defaults to \
       "$ATOMDB/apec_line.fits"
  cocofile : str
       The file containing all the continuum emission. Defaults to \
       "$ATOMDB/apec_coco.fits"
  binunits : {'keV','A'}
       The energy units for bins. "keV" or "A". Default keV.
  broadening : float
       Line broadening to be applied
  broadenunits : {'keV','A'}
       Units of line broadening "keV" or "A". Default keV.
  elements : iterable of int
       Elements to include, listed by atomic number. if not set, include all.
  abund : iterable of float, length same as elements.
       If set, and array of length (elements) with the abundances of each\
       element relative to the Andres and Grevesse values. Otherwise, assumed to\
       be 1.0 for all elements
  dummyfirst : bool
       If true, add a "0" to the beginning of the return array so it is of the
       same length as bins (can be useful for plotting results)
  dolines : bool
       Include lines in the spectrum
  docont : bool
       Include the continuum in the spectrum
  dopseudo : bool
       Include the pseudocontinuum in the spectrum.

  Returns
  -------
  array of floats
      Emissivity in counts cm^3 s^-1 bin^-1.

  """
#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015
#
#  Version 0.2
#    Added dummyfirst keyword
#    Adam Foster July 21st 2015
#
#  Version 0.3
#    Fixed bug in angstrom spectrum generation
#    Adam Foster February 28th 2018
#


  # set up the bins
  if (sum((bins[1:]-bins[:-1])<0) > 0):
    print("*** ERROR: bins must be monotonically increasing. Exiting ***")
    return -1

  if binunits.lower()=='kev':
    ebins = bins*1.0
    flipspectrum=False
  elif binunits.lower() in ['a', 'angstrom', 'angstroms']:
    ebins = const.HC_IN_KEV_A/bins[::-1]
    flipspectrum=True
  else:
    print("*** ERROR: unknown binning unit %s, Must be keV or A. Exiting ***"%\
          (binunits))





  if util.keyword_check(linefile):
    # ok, we should do something with this
    # if it is a string, look for the file name
    if isinstance(linefile, str):
      lfile = os.path.expandvars(linefile)
      if not os.path.isfile(lfile):
        print("*** ERROR: no such file %s. Exiting ***" %(lfile))
        return -1
      ldat = pyfits.open(lfile)
    elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      ldat = linefile
    else:
      print("Unknown data type for linefile. Please pass a string or an HDUList")
      return -1

  if util.keyword_check(cocofile):
    if isinstance(cocofile, str):
      cfile = os.path.expandvars(cocofile)
      if not os.path.isfile(cfile):
        print("*** ERROR: no such file %s. Exiting ***" %(cfile))
        return -1
      cdat = pyfits.open(cfile)
    elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      cdat = cocofile
    else:
      print("Unknown data type for cocofile. Please pass a string or an HDUList")
      return







#  lfile = os.path.expandvars(linefile)
#  cfile = os.path.expandvars(cocofile)
#  if not os.path.isfile(lfile):
#    print "*** ERROR: no such file %s. Exiting ***" %(lfile)
#    return -1
#  if not os.path.isfile(cfile):
#    print "*** ERROR: no such file %s. Exiting ***" %(cfile)
#    return -1

  # open the files
#  ldat = pyfits.open(lfile)
#  cdat = pyfits.open(cfile)

  # get the index
  #if ((index < 2) | (index > len(ldat))):
  #  print "*** ERRROR: Index must be in range %i to %i"%(2, len(ldat)-1)
  #  return -1

  lldat = ldat[index].data
  ccdat = cdat[index].data


  if not util.keyword_check(elements):
    Zl = util.unique(lldat['element'])
    Zc = util.unique(ccdat['Z'])
    Zlist = util.unique(numpy.append(Zl,Zc))

  else:
    Zlist = elements

  if not util.keyword_check(abund):
    abund= numpy.ones(len(Zlist))

  lspectrum = numpy.zeros(len(bins)-1, dtype=float)
  cspectrum = numpy.zeros(len(bins)-1, dtype=float)

  if dolines:
    for iZ, Z in enumerate(Zlist):
      # ADD  LINES
      lspectrum += add_lines(Z, abund[iZ], lldat, ebins, broadening=broadening, broadenunits=broadenunits)

  if docont | dopseudo:
    for iZ, Z in enumerate(Zlist):
    # ADD  CONTINUUM
      cspectrum += make_ion_index_continuum(ebins, Z, cocofile=ccdat,\
                                           binunits='keV', no_coco=not(docont),\
                                           no_pseudo=not(dopseudo))*abund[iZ]

  # broaden the continuum if required:
  if broadening:
    cspectrum = broaden_continuum(ebins, cspectrum, binunits = 'keV', \
                      broadening=broadening,\
                      broadenunits=broadenunits)

  # now flip results back around if angstroms
  if flipspectrum:
    cspectrum = cspectrum[::-1]
    lspectrum = lspectrum[::-1]



  if dummyfirst:
    return numpy.append([0],   cspectrum+lspectrum)
  else:
    return cspectrum+lspectrum



def make_ion_spectrum(bins, index, Z,z1, linefile="$ATOMDB/apec_nei_line.fits",\
                  cocofile="$ATOMDB/apec_nei_comp.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  abund=False, dummyfirst=False, nei = True,\
                  dolines = True, docont=True, dopseudo=True):

  """
  make_spectrum is the most generic "make me a spectrum" routine.

  It returns the emissivity in counts cm^3 s^-1 bin^-1.

  Parameters
  ----------
  bins : array(float)
       The bin edges for the spectrum to be calculated on, in \
       units of keV or Angstroms. Must be monotonically\
       increasing. Spectrum will return len(bins)-1 values.
  index : int
       The index to plot the spectrum from. note that the AtomDB files\
       the emission starts in hdu number 2. So for the first block, you\
       set index=2
  Z : int
       Element of spectrum (e.g. 6 for carbon)
  z1 : int
       Ion charge +1 for the spectrum (e.g. 3 for C III)
  linefile : str or HDUList
       The file containing all the line emission. Defaults to \
       "$ATOMDB/apec_line.fits". Can also pass in the opened file, \
       i.e. "linefile = pyatomdb.pyfits.open('apec_nei_line.fits')"
  cocofile : str or HDUList
       The file containing all the continuum emission. Defaults to \
       "$ATOMDB/apec_coco.fits". Can also pass in the opened file, \
       i.e. "cocofile = pyatomdb.pyfits.open('apec_nei_comp.fits')"
  binunits : {'keV','A'}
       The energy units for bins. "keV" or "A". Default keV.
  broadening : float
       Line broadening to be applied
  broadenunits : {'keV','A'}
       Units of line broadening "keV" or "A". Default keV.
  elements : iterable of int
       Elements to include, listed by atomic number. if not set, include all.
  abund : iterable of float, length same as elements.
       If set, and array of length (elements) with the abundances of each\
       element relative to the Andres and Grevesse values. Otherwise, assumed to\
       be 1.0 for all elements
  dummyfirst : bool
       If true, add a "0" to the beginning of the return array so it is of the
       same length as bins (can be useful for plotting results)
  nei : bool
       If set, return the spectrum from the driving ion being Z, rmJ. If not set,
       return the spectrum for the collisional ionization equilibrium *BUT*
       note that the continuum will be wrong, as it is provided for each element
       as a whole.
  dolines : bool
       Include lines in the spectrum
  docont : bool
       Include the continuum in the spectrum
  dopseudo : bool
       Include the pseudocontinuum in the spectrum.

  Returns
  -------
  array of floats
      Emissivity in counts cm^3 s^-1 bin^-1.

  """
#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015
#
#  Version 0.2
#    Added dummyfirst keyword
#    Adam Foster July 21st 2015
#


  # set up the bins
  if (sum((bins[1:]-bins[:-1])<0) > 0):
    print("*** ERROR: bins must be monotonically increasing. Exiting ***")
    return -1

  if binunits.lower()=='kev':
    ebins = bins*1.0
  elif binunits.lower() in ['a', 'angstrom', 'angstroms']:
    ebins = const.HC_IN_KEV_A/bins[::-1]
  else:
    print("*** ERROR: unknown binning unit %s, Must be keV or A. Exiting ***"%\
          (binunits))

  # check the files exist
  # first, check if the line file is set
  if util.keyword_check(linefile):
    # ok, we should do something with this
    # if it is a string, look for the file name
    if isinstance(linefile, str):
      if ((linefile == "$ATOMDB/apec_nei_line.fits") & (nei==False)):
        linefile = "$ATOMDB/apec_line.fits"
      lfile = os.path.expandvars(linefile)
      if not os.path.isfile(lfile):
        print("*** ERROR: no such file %s. Exiting ***" %(lfile))
        return -1
      ldat = pyfits.open(lfile)
    elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      ldat = linefile
    else:
      print("Unknown data type for linefile. Please pass a string or an HDUList")
      return -1

  if util.keyword_check(cocofile):
    if isinstance(cocofile, str):
      if ((cocofile == "$ATOMDB/apec_nei_comp.fits") & (nei==False)):
        cocofile = "$ATOMDB/apec_coco.fits"
      cfile = os.path.expandvars(cocofile)
      if not os.path.isfile(cfile):
        print("*** ERROR: no such file %s. Exiting ***" %(cfile))
        return -1
      cdat = pyfits.open(cfile)
    elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      cdat = cocofile
    else:
      print("Unknown data type for cocofile. Please pass a string or an HDUList")
      return




  # get the index
  #if ((index < 2) | (index > len(ldat))):
  #  print "*** ERRROR: Index must be in range %i to %i"%(2, len(ldat)-1)
  #  return -1

  lldat = ldat[index].data
  ccdat = cdat[index].data

  if not abund:
    abund= 1.0

  lspectrum = numpy.zeros(len(bins)-1, dtype=float)
  cspectrum = numpy.zeros(len(bins)-1, dtype=float)
  pspectrum = numpy.zeros(len(bins)-1, dtype=float)

  if dolines:
    # ADD  LINES
    if nei:
      lspectrum += add_lines(Z, abund, lldat, ebins, broadening=broadening,\
                             broadenunits=broadenunits,z1_drv=z1)
    else:
      lspectrum += add_lines(Z, abund, lldat, ebins,broadening, broadenunits,z1=z1)

  if docont:
    # ADD  LINES
    if nei:
      cspectrum += make_ion_index_continuum(ebins, Z, cocofile=ccdat,\
                                         ion = z1, binunits=binunits,\
                                         no_pseudo=True)*abund
    else:
      cspectrum += make_ion_index_continuum(ebins, Z, cocofile=ccdat,\
                                         ion = 0, binunits=binunits,\
                                         no_pseudo=True)*abund


  if dopseudo:
    # ADD  LINES
    if nei:
      pspectrum += make_ion_index_continuum(ebins, Z, cocofile=ccdat,\
                                         ion = z1, binunits=binunits, \
                                         no_coco=True)*abund
    else:
      pspectrum += make_ion_index_continuum(ebins, Z, cocofile=ccdat,\
                                         ion = 0, binunits=binunits, \
                                         no_coco=True)*abund



  # broaden the continuum if required:
  if broadening:

    cspectrum = broaden_continuum(ebins, cspectrum, binunits = binunits, \
                      broadening=broadening,\
                      broadenunits=broadenunits)
    pspectrum = broaden_continuum(ebins, pspectrum, binunits = binunits, \
                      broadening=broadening,\
                      broadenunits=broadenunits)
  if dummyfirst:
    return numpy.append([0],   cspectrum+lspectrum+pspectrum)
  else:
    return cspectrum+lspectrum+pspectrum


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def add_lines(Z, abund, lldat, ebins, z1=False, z1_drv=False, \
              broadening=False, broadenunits='A'):
  """
  Add lines to spectrum, applying gaussian broadening.

  Add the lines in list lldat, with atomic number Z, to a spectrum delineated
  by ebins (these are the edges, in keV). Apply broadening to the spectrum
  if broadening != False, with units of broadenunits (so can do constant
  wavelength or energy broadening)

  Parameters
  ----------
  Z : int
    Element of interest (e.g. 6 for carbon)
  abund : float
    Abundance of element, relative to AG89 data.
  lldat : dtype linelist
    The linelist to add. Usually the hdu from the apec_line.fits \
    file, often with some filters pre-applied.
  ebins : array of floats
    Energy bins. Will return spectrum with nbins-1 data points.
  z1 : int
    Ion charge +1 of ion to return
  z1_drv : int
    Driving Ion charge +1 of ion to return
  broadening : float
    Apply spectral broadening if > 0. Units of A of keV
  broadenunits : {'A' , 'keV'}
    The units of broadening, Angstroms or keV


  Returns
  -------
  array of float
    broadened emissivity spectrum, in photons cm^3 s^-1 bin^-1. Array has \
    len(ebins)-1 values.
  """
#
#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015
#

  lammax = (const.HC_IN_KEV_A/ebins[0])
  lammin = (const.HC_IN_KEV_A/ebins[-1])
  if broadenunits.lower() in ['a','angstrom','angstroms']:
    bunits = 'a'
  elif broadenunits.lower() =='kev':
    bunits = 'kev'
  else:
    print("Error: unknown broadening unit %s, Must be keV or A. Exiting ***"%\
          (broadenunits))
    return -1
  if broadening:
    if bunits == 'a':
      lammax += broadening
      lammin -= broadening
    else:
      lammax += lammax**2 * broadening/const.HC_IN_KEV_A
      lammin -= lammin**2 * broadening/const.HC_IN_KEV_A

  l = lldat[(lldat['element']==Z) &\
            (lldat['lambda'] <= lammax) &\
            (lldat['lambda'] >= lammin)]

  if z1:
    l = l[l['ion'] ==z1]
  if z1_drv:
    l = l[l['ion_drv'] ==z1_drv]

  spectrum = numpy.zeros(len(ebins)-1, dtype=float)

  if broadening:
    if  bunits == 'a':
      for ll in l:
        spectrum+=atomdb.addline2(ebins, const.HC_IN_KEV_A/ll['lambda'], \
                 ll['epsilon']* abund,\
                 broadening*const.HC_IN_KEV_A/(ll['lambda']**2))
    else:
      for ll in l:
        spectrum+=atomdb.addline2(ebins, const.HC_IN_KEV_A/ll['lambda'], \
                 ll['epsilon']* abund,\
                 broadening)
  else:
    for ll in l:
      spectrum[numpy.argmax(\
               numpy.where(ebins <= const.HC_IN_KEV_A/ll['lambda'])[0])]+=\
               ll['epsilon']* abund

  return spectrum
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def get_index(te, filename='$ATOMDB/apec_line.fits', \
              teunits='keV', logscale=False):
  """
  Finds HDU with kT closest ro desired kT in given line or coco file.

  Opens the line or coco file, and looks for the header unit
  with temperature closest to te. Use result as index input to make_spectrum

  Parameters
  ----------
  te : float
    Temperature in keV or K
  teunits : {'keV' , 'K'}
    Units of te (kev or K, default keV)
  logscale : bool
    Search on a log scale for nearest temperature if set.
  filename : str or hdulist
    line or continuum file, already opened or filename.

  Returns
  -------
  int
    Index in HDU file with nearest temperature to te.
  """

#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015
#
#
#  Version 0.2 - fixed bug so teunits = l works properly
#    Adam Foster Jan 26th 2016

  if teunits.lower() == 'kev':
    teval = te
  elif teunits.lower() == 'k':
    teval = te*const.KBOLTZ
  else:
    print("*** ERROR: unknown temeprature unit %s. Must be keV or K. Exiting ***"%\
          (teunits))
  if type(filename) == pyfits.hdu.hdulist.HDUList:
    a = filename[1].data
  elif type(filename) == pyfits.hdu.table.BinTableHDU:
    a = filename.data
  elif type(filename)==type('somestring'):
    a = pyfits.open(os.path.expandvars(filename))[1].data


  if logscale:
    i = numpy.argmin(numpy.abs(numpy.log(a['kT'])-numpy.log(teval)))
  else:
    i = numpy.argmin(numpy.abs(a['kT']-teval))
  # need to increase the HDU by 2.
  return i+2



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def list_lines(specrange, lldat=False, index=False, linefile=False,\
              units='angstroms', Te=False, teunit='K', minepsilon=1e-20):
  """
  Gets list of the lines in a given spectral range

  Note that the output from this can be passed directly to print_lines


  Parameters
  ----------
  specrange : [float,float]
    spectral range [min,max] to return lines on
  lldat : see notes
    line data
  index : int
    index in lldat, see notes
  linefile : see notes
    line data file, see notes
  units : {'A' , 'keV'}
    units of specrange (default A)
  Te : float
    electron temperature (used if index not set to select appropriate
    data HDU from line file)
  teunit : {'K','keV','eV'}
    units of Te
  minepsilon : float
    minimum epsilon for lines to be returned, in ph cm^3 s^-1

  Notes
  -----
  The actual line list can be defined in one of several ways:

  specrange = [10,100]

  1. lldat as an actual list of lines::

       a = pyfits.open('apec_line.fits')
       llist = a[30].data
       l = list_lines(specrange, lldat=llist)

  2. lldat as a numpy array of lines::

       a = pyfits.open('apec_line.fits')
       llist = numpy.array(a[30].data)
       l = list_lines(specrange, lldat=llist)

  3. lldat is a BinTableHDU from pyfits::

       a = pyfits.open('apec_line.fits')
       llist = numpy.array(a[30])
       l = list_lines(specrange, lldat=llist)

  4. lldat is a HDUList from pyfits. In this case index must also be set::

       a = pyfits.open('apec_line.fits')
       index = 30
       l = list_lines(specrange, lldat=a, index=index)

  5. lldat NOT set, linefile contains apec_line.fits file location, index
     identifies the HDU::

       linefile = 'mydir/apec_v2.0.2_line.fits'
       index = 30
       l = list_lines(specrange, linefile=linefile, index=index)

  6. lldat NOT set & linefile NOT set, linefile is set to
     $ATOMDB/apec_line.fits. index identifies the HDU::

       index = 30
       l = list_lines(specrange, index=index)

  Returns
  -------
  linelist : dtype=([('Lambda', '>f4'), \
           ('Lambda_Err', '>f4'), \
           ('Epsilon', '>f4'), \
           ('Epsilon_Err', '>f4'), \
           ('Element', '>i4'), \
           ('Ion', '>i4'), \
           ('UpperLev', '>i4'), \
           ('LowerLev', '>i4')])
     A line list filtered by the various elements.
  """

#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015
#
#  Version 0.2 - bugfix:
#    corrected lldat test to be if not false, instead of if true.
#    Adam Foster September 16th 2015
#

  # check the units

  if units.lower()=='kev':
    specrange = [const.HC_IN_KEV_A/specrange[1], const.HC_IN_KEV_A/specrange[0]]
  elif units.lower() in ['a', 'angstrom', 'angstroms']:
    specrange = specrange
  else:
    print("*** ERROR: unknown unit %s, Must be keV or A. Exiting ***"%\
          (units))

  # open the line file if only specified by a name
  if lldat==False:
    # use default file if not specified
    if linefile==False:
      linefile=os.path.expandvars("$ATOMDB/apec_line.fits")
    # open file
    lldat = pyfits.open(linefile)

  # convert temperature units

  if Te != False:
    if teunit.lower() == 'kev':
      kT = Te*1.0
    elif teunit.lower() == 'ev':
      kT = Te/1000.0
    elif teunit.lower() == 'k':
      kT = Te*const.KBOLTZ
    else:
      print("*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunit))

  # if the temperature is specified, get the index
  if index != False:

    if Te != False:
      print("Warning: both index and Te specified. Using index")
    else:
      # everything is fine!
      pass

  else:
    if Te != False:
      # open the file and get the index of nearest block
      index =  get_index(kT, filename=lldat, \
              teunits='keV', logscale=True)
    else:
      if not type(lldat) in [pyfits.fitsrec.FITS_rec, numpy.ndarray]:
        print("Error: did not specify index or Te")
        return False
    # index is specified, so we'll use it.
    pass




  if lldat != False:
    #options here:
    # (1) This is a line list, i.e. the ldata[index].data from a file,
    #     either in original pyfits format or a numpy array
    #
    # (2) This is an hdu from a file
    #
    # (3) This is a _line.fits file, and requires an index to make sense of it

    if type(lldat) == pyfits.hdu.hdulist.HDUList:
      if not(index):
        print("*** ERROR. lldat provided as HDUList, but no index specified.")
        print(" Exiting")
      llist = numpy.array(lldat[index].data)
    elif type(lldat) == pyfits.hdu.table.BinTableHDU:
      llist = numpy.array(lldat.data)
    elif type(lldat) in [pyfits.fitsrec.FITS_rec, numpy.ndarray]:
      llist = numpy.array(lldat)
    else:
      print("ERROR: unkonwn llist type!")
  else:
    # there should now be no way to get here. commenting out this section
    print("ERROR: I SHOULD NEVER BE HERE")
    pass
    # no line data supplied.
    #if linefile==False:
      #linefile = os.path.expandvars('$ATOMDB/apec_line.fits')
    #if not os.path.isfile(linefile):
      #print "*** ERROR. No linefile supplied but $ATOMDB/apec_line.fits is"
      #print " not a file. Exiting"
    #else:

      #if not index:
        #print "*** ERROR. No index specified for line list. Exiting"
      #else:
        #lfile = pyfits.open(os.path.expandvars(linefile))
        #llist= numpy.array(lfile[index].data)

  # at this point, we have data. Filter by specrange, minepsilon
  llist = llist[(llist['Lambda']>= specrange[0]) &\
                (llist['Lambda']<= specrange[1]) &\
                (llist['Epsilon']>= minepsilon)]



  return llist


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def list_nei_lines(specrange, Te, tau, Te_init=False,  lldat=False, linefile=False,\
              units='angstroms', teunit='K', minepsilon=1e-20, \
              datacache=False):
  """
  Gets list of the lines in a given spectral range for a given NEI plasma

  For speed purposes, this takes the nearest temperature tabulated in the
  linefile, and applies the exact ionization balance as calculated to this.
  This is not perfect, but should be good enough.

  Note that the output from this can be passed directly to print_lines


  Parameters
  ----------
  specrange : [float,float]
    spectral range [min,max] to return lines on
  Te : float
    electron temperature
  tau : float
    electron density * time (cm^-3 s)
  Te_init : float
    initial ionization balance temperature
  lldat : see notes
    line data
  linefile : see notes
    line data file, see notes
  units : {'A' , 'keV'}
    units of specrange (default A)
  teunit : {'K' , 'keV'}
    units of temperatures (default K)
  minepsilon : float
    minimum emissivity (ph cm^3 s^{-1}) for inclusion in linelist

  Notes
  -----
  The actual line list can be defined in one of several ways:

  specrange = [10,100]

  1. lldat as an actual list of lines::

       a = pyfits.open('apec_nei_line.fits')
       llist = a[30].data
       l = list_nei_lines(specrange, lldat=llist)

  2. lldat as a numpy array of lines::

       a = pyfits.open('apec_nei_line.fits')
       llist = numpy.array(a[30].data)
       l = list_nei_lines(specrange, lldat=llist)

  3. lldat is a BinTableHDU from pyfits::

       a = pyfits.open('apec_nei_line.fits')
       llist = numpy.array(a[30])
       l = list_nei_lines(specrange, lldat=llist)

  4. lldat is a HDUList from pyfits. In this case index must also be set::

       a = pyfits.open('apec_nei_line.fits')
       index = 30
       l = list_nei_lines(specrange, lldat=a, index=index)

  5. lldat NOT set, linefile contains apec_line.fits file location, index
     identifies the HDU::

       linefile = 'mydir/apec_v3.0.2_nei_line.fits'
       index = 30
       l = list_nei_lines(specrange, linefile=linefile, index=index)

  6. lldat NOT set & linefile NOT set, linefile is set to
     $ATOMDB/apec_line.fits. index identifies the HDU::

       index = 30
       l = list_nei_lines(specrange, Te, tau)

  Returns
  -------
  linelist : dtype=([('Lambda', '>f4'), \
           ('Lambda_Err', '>f4'), \
           ('Epsilon', '>f4'), \
           ('Epsilon_Err', '>f4'), \
           ('Element', '>i4'), \
           ('Elem_drv', '>i4'), \
           ('Ion', '>i4'), \
           ('Ion_drv', '>i4'), \
           ('UpperLev', '>i4'), \
           ('LowerLev', '>i4')])
     A line list filtered by the various elements.
  """

#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster November 02nd 2015
#

  # check the units

  if units.lower()=='kev':
    specrange = [const.HC_IN_KEV_A/specrange[1], const.HC_IN_KEV_A/specrange[0]]
  elif units.lower() in ['a', 'angstrom', 'angstroms']:
    specrange = specrange
  else:
    print("*** ERROR: unknown unit %s, Must be keV or A. Exiting ***"%\
          (units))

  # convert Te into keV

  if teunit.lower() == 'kev':
    kT = Te*1.0
  elif teunit.lower() == 'ev':
    kT = Te/1000.0
  elif teunit.lower() == 'k':
    kT = Te*const.KBOLTZ
  else:
    print("*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunit))


  if Te_init != False:
    if teunit.lower() == 'kev':
      kT_init = Te_init*1.0
    elif teunit.lower() == 'ev':
      kT_init = Te_init/1000.0
    elif teunit.lower() == 'k':
      kT_init = Te_init*const.KBOLTZ
    else:
      print("*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunit))
  else:
  # Te_init was not set:
    kT_init = 1e4*const.KBOLTZ



  # sort out the line file...

  if lldat != False:
    #options here:
    # (1) This is a line list, i.e. the ldata[index].data from a file,
    #     either in original pyfits format or a numpy array
    #
    # (2) This is an hdu from a file
    #
    # (3) This is a _line.fits file, and requires an index to make sense of it

    if type(lldat) == pyfits.hdu.hdulist.HDUList:
      # go get the index
      te_index = get_index(kT, filename=lldat, \
              teunits='keV', logscale=True)
      llist = numpy.array(lldat[te_index].data)
    elif type(lldat) == pyfits.hdu.table.BinTableHDU:
      # no need to get index
      llist = numpy.array(lldat.data)
    elif type(lldat) in [pyfits.fitsrec.FITS_rec, numpy.ndarray]:
      llist = numpy.array(lldat)
  else:
    # no line data supplied.
    if linefile==False:
      linefile = os.path.expandvars('$ATOMDB/apec_nei_line.fits')
    if not os.path.isfile(linefile):
      print("*** ERROR. Linefile %s is "%(linefile), end=' ')
      print(" not a file. Exiting")
    else:
      lldat = pyfits.open(os.path.expandvars(linefile))
      te_index = get_index(kT, filename=lldat, \
              teunits='keV', logscale=True)
      llist= numpy.array(lldat[te_index].data)

  # get filtered line list

  llist = llist[(llist['Lambda']>= specrange[0]) &\
                (llist['Lambda']<= specrange[1]) &\
                (llist['Epsilon'] >= minepsilon)]

  # get the index

  # get list of all the elements present
  Zlist = util.unique(llist['Element'])
  # Calculate the ionization balance.
  ionbal ={}
  for Z in Zlist:
    ionbal[Z] = apec.solve_ionbal_eigen(Z, kT, tau=tau, Te_init = kT_init,\
                                    teunit='keV', datacache=datacache)
  # multiply everything by the appropriate ionization fraction
  if 'Elem_drv' in llist.dtype.names:

    for il in llist:
      il['Epsilon'] *= ionbal[il['Elem_drv']][il['Ion_drv']-1]
  else:
    for il in llist:
      il['Epsilon'] *= ionbal[il['Element_drv']][il['Ion_drv']-1]

  # filter again based on new epsilon values
  llist=llist[llist['Epsilon']>minepsilon]
  # at this point, we have data
  return llist
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def print_lines(llist, specunits = 'A', do_cfg=False):
  """
  Prints lines in a linelist to screen

  This routine is very primitive as things stand. Plenty of room for refinement.

  Parameters
  ----------
  llist: dtype(linelist)
    list of lines to print. Typically returned by list_lines.
  specunits: {'A' , 'keV'}
    units to list the line positions by (A or keV, default A)
  do_cfg: bool
    Show full configuration information for each level

  Returns
  -------
  Nothing, though prints data to standard out.
  """

#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015
#

  if specunits.lower()=='kev':
    specunits = 'keV'
    llist['Lambda']=const.HC_IN_KEV_A/llist['Lambda']
  elif specunits.lower() in ['a', 'angstrom', 'angstroms']:
    specunits = 'A'
  else:
    print("*** ERROR: unknown unit %s, Must be keV or A. Exiting ***"%\
          (specunits))


  # now print the header lines

  if specunits == 'keV':
    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Energy','Epsilon','Element','Ion','Ion_drv','UpperLev','LowerLev')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Energy','Epsilon','Element','Ion','UpperLev','LowerLev')
    print(s)

    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('keV','ph cm3 s-1','','','','','')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('keV','ph cm3 s-1','','','','')
    print(s)

  else:
    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Lambda','Epsilon','Element','Ion','Ion_drv','UpperLev','LowerLev')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Lambda','Epsilon','Element','Ion','UpperLev','LowerLev')
    print(s)

    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('A','ph cm3 s-1','','','','','')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('A','ph cm3 s-1','','','','')
    print(s)

  # now print the data
  d={}
  if 'Ion_drv' in llist.dtype.names:

    for il in llist:
      s = "%.4e %.4e %-10i %-10i %-10i %-10i %-10i" %\
         (il['Lambda'],\
          il['Epsilon'],\
          il['Element'],\
          il['Ion'],\
          il['Ion_drv'],\
          il['UpperLev'],\
          il['LowerLev'])

      # get the configuration info
      if do_cfg:
        lvdat = atomdb.get_data(il['Element'],il['Ion'],'LV',\
                                         datacache=d)

        s+= " %40s %3f %2i %2i"%(lvdat[1].data['ELEC_CONFIG'][il['UpperLev']-1],\
                                 lvdat[1].data['S_QUAN'][il['UpperLev']-1],\
                                 lvdat[1].data['L_QUAN'][il['UpperLev']-1],\
                                 lvdat[1].data['LEV_DEG'][il['UpperLev']-1])

        s+= " -> %40s %3f %2i %2i"%(lvdat[1].data['ELEC_CONFIG'][il['LowerLev']-1],\
                                 lvdat[1].data['S_QUAN'][il['LowerLev']-1],\
                                 lvdat[1].data['L_QUAN'][il['LowerLev']-1],\
                                 lvdat[1].data['LEV_DEG'][il['LowerLev']-1])
      else:
        pass
      print(s)
  else:
    for il in llist:
      s = "%.4e %.4e %-10i %-10i %-10i %-10i" %\
         (il['Lambda'],\
          il['Epsilon'],\
          il['Element'],\
          il['Ion'],\
          il['UpperLev'],\
          il['LowerLev'])
      if do_cfg:
        lvdat = pyatomdb.atomdb.get_data(il['Element'],il['Ion'],'LV',\
                                         datacache=d)

        s+= " %40s %3f %2i %2i"%(lvdat[1].data['ELEC_CONFIG'][il['UpperLev']-1],\
                                 lvdat[1].data['S_QUAN'][il['UpperLev']-1],\
                                 lvdat[1].data['L_QUAN'][il['UpperLev']-1],\
                                 lvdat[1].data['LEV_DEG'][il['UpperLev']-1])

        s+= " -> %40s %3f %2i %2i"%(lvdat[1].data['ELEC_CONFIG'][il['LowerLev']-1],\
                                 lvdat[1].data['S_QUAN'][il['LowerLev']-1],\
                                 lvdat[1].data['L_QUAN'][il['LowerLev']-1],\
                                 lvdat[1].data['LEV_DEG'][il['LowerLev']-1])
      else:
        pass
      print(s)






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def make_ion_index_continuum(bins,  element, \
                             index = False,\
                             cocofile='$ATOMDB/apec_coco.fits',\
                             binunits = 'keV', \
                             fluxunits='ph', no_coco=False, no_pseudo=False, \
                             ion=0,\
                             broadening=False,\
                             broadenunits='keV'):
  """
  Creates the continuum for a given ion.

  Parameters
  ----------
  bins : array(float)
    The bin edges for the spectrum to be calculated on, in units of keV \
    or Angstroms. Must be monotonically increasing. Spectrum will return \
    len(bins)-1 values.
  element : int
    Atomic number of element to make spectrum of (e.g. 6 for carbon)
  binunits : {'keV' , 'A'}
    The energy units for bins. "keV" or "A". Default keV.
  fluxunits : {'ph' , 'erg'}
    Whether to return the emissivity in photons ('ph') or ergs ('erg').\
    Defaults to  photons
  no_coco : bool
    If true, do not include the compressed continuum
  no_pseudo : bool
    If true, do not include the pseudo continuum (weak lines)
  ion : int
    Ion to calculate, e.g. 4 for C IV. By default, 0 (whole element).
  index : int
    The index to generate the spectrum from. Note that the AtomDB files
    the emission starts in hdu number 2. So for the first block, you
    set index=2. Only required if cocofile is a filename or an HDULIST
  cocofile : HDUList, HDU or str
    The continuum file, either already open (HDULIST) or filename.
    alternatively, provide the HDU itself, and then do not need
    to define the index
  broadening: float
    Broaden the continuum by gaussians of this width (if False,
    no broadening is applied)
  broadenunits: {'keV' , 'A'}
    Units for broadening (kev or A)

  Returns
  -------
  array(float)
    len(bins)-1 array of continuum emission, in units of \
    photons cm^3 s^-1 bin^-1 (fluxunits = 'ph') or
    ergs cm^3 s^-1 bin^-1 (fluxunits = 'erg')

  """
#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015


  if binunits.lower() in ['kev']:
    angstrom = False
  elif binunits.lower() in ['a','angstrom','angstroms']:
    angstrom = True
  else:
    print("*** ERROR: unknown units %s for continuum spectrum. Exiting" %\
          (binunits))

  if fluxunits.lower() in ['ph', 'photon','photons', 'p']:
    ergs = False
  elif fluxunits.lower() in ['erg','ergs']:
    ergs = true
  else:
    print("*** ERROR: unknown units %s for continuum flux. Exiting" %\
          (fluxunits))

    return -1

  if angstrom:
    # need to work in keV for this bit
    bins = const.HC_IN_KEV_A/bins[::-1]

  # open the data files

  if type(cocofile) == pyfits.hdu.hdulist.HDUList:
    cdat = cocofile[index].data
  elif type(cocofile) == pyfits.hdu.table.BinTableHDU:
    cdat = cocofile.data
  elif type(cocofile) == pyfits.fitsrec.FITS_rec:
    cdat = cocofile
  elif type(cocofile) == type('somestring'):
    cdat = pyfits.open(os.path.expandvars(cocofile))[index].data
  else:
    print("*** ERROR: unable to parse cocofile = %s" %repr(cocofile))
    return -1


  nbins = len(bins)-1
  spectrum = numpy.zeros(nbins, dtype=float)
  Z = element
  rmj = ion

  d = cdat[(cdat['Z']==Z) & (cdat['rmJ'] == rmj)]
  if len(d)==0:
    return spectrum
  else:
    d = d[0]
  if not(no_coco):
    spectrum += expand_E_grid(bins,\
                                d['N_cont'],\
                                d['E_cont'],\
                                d['Continuum'])

  if not(no_pseudo):
    spectrum += expand_E_grid(bins,\
                                d['N_Pseudo'],\
                                d['E_Pseudo'],\
                                d['Pseudo'])


  if (ergs):
    spectrum[iBin] *= const.ERGperKEV*0.5*(bins[:-1]+bins[1:])


  if (angstrom):
    spectrum = spectrum[::-1]  # reverse spectrum to match angstroms

  return spectrum

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def expand_E_grid(eedges, n,Econt_in_full, cont_in_full):

  """
  Code to expand the compressed continuum onto a series of bins.

  Parameters
  ----------
  eedges : float(array)
    The bin edges for the spectrum to be calculated on, in units of keV
  n : int
    The number of good data points in the continuum array
  Econt_in_full: float(array)
    The compressed continuum energies
  cont_in_full: float(array)
    The compressed continuum emissivities

  Returns
  -------
  float(array)
    len(bins)-1 array of continuum emission, in units of \
    photons cm^3 s^-1 bin^-1
  """

#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015

  import scipy.integrate
  cont_in = cont_in_full[:n]
  Econt_in = Econt_in_full[:n]


  # ok. So. Sort this.
  E_all = numpy.append(Econt_in, eedges)

  cont_tmp = numpy.interp(eedges, Econt_in, cont_in)
  C_all = numpy.append(cont_in, cont_tmp)

  iord = numpy.argsort(E_all)

  # order the arrays
  E_all = E_all[iord]
  C_all = C_all[iord]

  ihi = numpy.where(iord>=n)[0]
  cum_cont = scipy.integrate.cumtrapz(C_all, E_all, initial=0)
  C_out = numpy.zeros(len(eedges))
#  for i in range(len(eedges)):
#    arg  = numpy.argwhere(iord==n+i)[0][0]
#    C_out[i] = cum_cont[ihi[i]]

  C_out = cum_cont[ihi]

  cont = C_out[1:]-C_out[:-1]
  return cont



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def broaden_continuum(bins, spectrum, binunits = 'keV', \
                      broadening=False,\
                      broadenunits='keV'):
  """
  Apply a broadening to the continuum

  Parameters
  ----------
  bins : array(float)
    The bin edges for the spectrum to be calculated on, in units of \
    keV or Angstroms. Must be monotonically increasing. Spectrum \
    will return len(bins)-1 values.
  spectrum : array(float)
    The emissivities in each bin in the unbroadened spectrum
  binunits : {'keV' , 'A'}
    The energy units for bins. "keV" or "A". Default keV.
  broadening : float
    Broaden the continuum by gaussians of this width (if False,\
    no broadening is applied)
  broadenunits : {'keV' , 'A'}
    Units for broadening (kev or A)

  Returns
  -------
  array(float)
    spectrum broadened by gaussians of width broadening
  """

#  History
#  -------
#  Version 0.1 - initial release
#    Adam Foster July 17th 2015

  # convert to energy grid
  if binunits.lower() in ['kev']:
    angstrom = False
  elif binunits.lower() in ['a','angstrom','angstroms']:
    angstrom = True
  else:
    print("*** ERROR: unknown units %s for continuum spectrum. Exiting" %\
          (binunits))

  if angstrom:
    bins = const.HC_IN_KEV_A/bins[::-1]

  # broadening
  if broadening:
    if broadenunits.lower() in ['a','angstrom','angstroms']:
      bunits = 'a'
    elif broadenunits.lower() in ['kev']:
      bunits = 'kev'
    else:
      print("*** ERROR: unknown units %s for continuum broadening. Exiting" %\
            (broadenunits))
      return -1
    # do the broadening
    spec = numpy.zeros(len(spectrum))
    emid = (bins[1:]+bins[:-1])/2
    if broadenunits == 'a':
      # convert to keV
      broadenvec = emin**2 *broadening/const.HC_IN_KEV_A
    else:
      broadenvec = numpy.zeros(len(emid))
      broadenvec[:] = broadening
    for i in range(len(spec)):

      spec += atomdb.addline2(bins, emid[i], \
                 spectrum[i],\
                 broadenvec[i])
    spectrum=spec
  if angstrom:
    spectrum=spectrum[::-1]
  return spectrum

def get_response_ebins(rmf):
  """
  Get the energy bins from the rmf file

  Parameters
  ----------
  rmf : string or pyfits.hdu.hdulist.HDUList
    The filename of the rmf or the opened rmf file

  Returns
  -------
  specbins_in : array(float)
    input energy bins used. nbins+1 length, with the last item being the final bin
    This is the array on which the input spectrum should be calculated
  specbins_out : array(float)
    output energy bins used. nbins+1 length, with the last item being the final bin
    This is the array on which the output spectrum will be returned. Often (but not
    always) the same as specbins_in
  """
  #
  # Update 2016-05-25
  #
  # Changed return to be the MATRIX, not EBOUNDS
  # In many instruments these are the same grids, but not all.


  if type(rmf)==str:
    rmfdat = pyfits.open(rmf)
  elif type(rmf) == pyfits.hdu.hdulist.HDUList:
    rmfdat = rmf
  else:
    print("ERROR: unknown rmf type, %s"%(repr(type(rmf))))
    return
  try:
    k=rmfdat.index_of('MATRIX')
    matrixname = 'MATRIX'
  except KeyError:
    try:
      k=rmfdat.index_of('SPECRESP MATRIX')
      matrixname = 'SPECRESP MATRIX'
    except KeyError:
      print("Cannot find index for matrix in this data")
      raise


  specbins_in = rmfdat[matrixname].data['ENERG_LO']
  specbins_in = numpy.append(specbins_in, rmfdat[matrixname].data['ENERG_HI'][-1])


  specbins_out = rmfdat['EBOUNDS'].data['E_MIN']
  specbins_out = numpy.append(specbins_out , rmfdat['EBOUNDS'].data['E_MAX'][-1])


  return specbins_in, specbins_out


def get_effective_area(rmf, arf=False):
  """
  Get the effective area of a response file

  Parameters
  ----------
  rmf : string or pyfits.hdu.hdulist.HDUList
    The filename of the rmf or the opened rmf file
  arf : string or pyfits.hdu.hdulist.HDUList
    The filename of the arf or the opened arf file
  Returns
  -------
  array(float)
    energy grid (keV) for returned response
  array(float)
    effective area for the returned response
  """
#
# Update 2016-05-25
#
# Changed to return the energy grid and the spectrum, as apparently in some
# instruments these are not the same as the input energy grid.

  if arf:
    if type(arf)==str:
      arfdat = pyfits.open(arf)
    elif type(arf) == pyfits.hdu.hdulist.HDUList:
      arfdat = arf
    else:
      print("ERROR: unknown arf type, %s"%(repr(type(arf))))
      return
    arfarea = arfdat['SPECRESP'].data['SPECRESP']
  else:
    arfarea = 1.0


  if type(rmf)==str:
    rmfdat = pyfits.open(rmf)
  elif type(rmf) == pyfits.hdu.hdulist.HDUList:
    rmfdat = rmf
  else:
    print("ERROR: unknown rmf type, %s"%(repr(type(rmf))))
    return

  ebins_in, ebins_out = get_response_ebins(rmf)

  area = numpy.zeros(len(ebins_in)-1, dtype=float)


  try:
    k=rmfdat.index_of('MATRIX')
    matrixname = 'MATRIX'
  except KeyError:
    try:
      k=rmfdat.index_of('SPECRESP MATRIX')
      matrixname = 'SPECRESP MATRIX'
    except KeyError:
      print("Cannot find index for matrix in this data")
      raise

  matname = 'MATRIX'
  if not matname in rmfdat[matrixname].data.names:
    matname ='SPECRESP MATRIX'
    if not matname in rmfdat[matrixname].data.names:
      print("Error: Cannot find Matrix in rmf data")
      return
  for ibin, i in enumerate(rmfdat[matrixname].data):
    area[ibin] = sum(i[matname])

  area *= arfarea
  return ebins_in, area


class CIESession():
  """
  Load and generate a collisional ionization equilibrium spectrum

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : arraylike(int), optional
    The atomic number of elements to include (default all)
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  datacache : dict
    Any Atomdb FITS files which have to be opened are stored here
  spectra : CIESpectra
    Object storing the actual spectral data
  elements : list(int)
    Nuclear charge of elements to include.
  default_abundset : string
    The abundance set used for the original emissivity file calculation
  abundset : string
    The abundance set to be used for the returned spectrum
  abundsetvector : array_like(float)
    The relative abundance between default_abundset and abundset for each element
  response_set : bool
    Have we loaded a response (or set a dummy response)
  dolines : bool
    Calculate line emission
  docont : bool
    Calculate continuum emission
  dopseudo : bool
    Calculate pseudocontinuum emission
  broaden_limit : float
    Apply broadening to lines with epsilon > this value (ph cm3 s-1)
  thermal_broadening : bool
    Apply thermal broadening to lines (default = False)
  velocity_broadening : float
    Apply velocity broadening with this velocity (km/s). If <=0, do not apply.

  Examples
  --------

  Create a session instance:

  >>> s=CIESession()

  Set up the responses, in this case a dummy response from 0.1 to 10 keV

  >>> ebins = numpy.linspace(0.1,10,1000)
  >>> s.set_response(ebins, raw=True)

  Turn on thermal broadening

  >>> s.set_broadening(True)
  Will thermally broaden lines with emissivity > 1.000000e-18 ph cm3 s-1

  Return spectrum at 1.0keV

  >>> spec = s.return_spectrum(1.0)

  spec is in photons cm^3 s^-1 bin^-1; ebins are the bin edges (so spec is
  1 element shorter than ebins)
  """

  def __init__(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits",\
                     elements=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30],\
                     abundset='AG89'):
    """
    Initialization routine. Can set the line and continuum files here

    Input
    -----
    linefile : str or HDUList
      The filename of the line emissivity data, or the opened file.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the opened file.
    elements : array_like(int)
      The atomic numbers of the elements to include. Defaults to all (1-30)
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance
      for list of options.
    """


    self.datacache={}

    # Open up the APEC files
    self.set_apec_files(linefile, cocofile)



    # if elements are specified, use them. Otherwise, use Z=1-30
    if util.keyword_check(elements):
      self.elements = elements
    else:
      self.elements=list(range(1,const.MAXZ_CIE+1))

    # a hold for the spectra
    self.spectra=CIESpectrum(self.linedata, self.cocodata)

    # Set both the current and the default abundances to those that
    # the apec data was calculated on
    self.abundset=self.linedata[0].header['SABUND_SOURCE']
    self.default_abundset=self.linedata[0].header['SABUND_SOURCE']

    self.abundsetvector = numpy.zeros(const.MAXZ_CIE+1)
    for Z in self.elements:
      self.abundsetvector[Z] = 1.0

    #  but if another vector was already specified, use this instead
    if util.keyword_check(abundset):
      self.set_abundset(abundset)

    self.abund = numpy.zeros(const.MAXZ_CIE+1)

    for Z in self.elements:
      self.abund[Z]=1.0

    # Set a range of parameters which can be overwritten later
    self.response_set = False # have we loaded a response file?
    self.dolines=True # Include lines in spectrum
    self.docont=True # Include continuum in spectrum
    self.dopseudo=True # Include pseudo continuum in spectrum

    #self.broaden_limit = 1e-18
    #self.thermal_broadening=True
    #self.velocity_broadening=0.0

    self.set_broadening(False, broaden_limit=1e-18)


  def set_broadening(self, thermal_broadening, broaden_limit=False, \
                           velocity_broadening=0.0, \
                           velocity_broadening_units='km/s'):

    """
    Turn on or off thermal broadening, and the emissivity limit for
    thost lines

    Parameters
    ----------

    thermal_broadening : bool
      If true, turn on broadening. If False, turn it off.
    broaden_limit : float
      The emissivity limit for lines to be broadened. If False, this value
      will not be updated.
    velocity_broadening : float
      velocity broadening to apply. If <=0, not applied
    velocity_broadening_units : string
      Units of velocity_broadening. 'km/s' is default and only value so far.

    Notes
    -----
    Updates attributes thermal_broadening, broaden_limit, velocity_broadening,
    velocity_broadening_units
    """

    self.thermal_broadening = thermal_broadening
    if broaden_limit != False:
      self.broaden_limit = broaden_limit

    if self.thermal_broadening==True:
      print("Will thermally broaden lines with emissivity > %e ph cm3 s-1"%(self.broaden_limit))
    else:
      print("Will not thermally broaden lines")

    self.velocity_broadening=velocity_broadening

    allowed_velocity_broadening_units= ['km/s']

    if not velocity_broadening_units.lower() in ['km/s']:
      print("Error: velocity broadening units of %s is not in allowed set "%\
             (velocity_broadening_units), allowed_velocity_broadening_units)
      return

    self.velocity_broadening_units=velocity_broadening_units

    self.spectra.thermal_broadening = self.thermal_broadening
    self.spectra.broaden_limit = self.broaden_limit

    self.spectra.velocity_broadening=self.velocity_broadening
    self.spectra.velocity_broadening_units=self.velocity_broadening_units




  def set_response(self, rmf, arf=False, raw=False):
    """
    Set the response. rmf, arf can either be the filenames or the
    opened files (latter is faster if called repeatedly)


    Extended Summary
    ----------------
    Amends the following items:

    self.rmffile : string
      The rmf file name
    self.rmf : string
      The response matrix
    self.arffile : string
      The arf file name
    self.arf : string
      The arf data
    self.specbins : array(float)
      The spectral bins on which to calculate the spectrum (keV or A)
    self.specbin_units : string ['A','keV']
      Units of specbins
    self.ebins : array(float)
      The spectral bins on which to calculate the spectrum (keV).
    self.ebins_out : array(float)
      The spectral bins on which to return the spectrum (keV). Can be
      different from specbins depending on the spectrum
    self.response_set : bool
      A resonse has been loaded
    self.specbins_set : bool
      The spectral bins are set
    self.ebins_checksum : string
      The md5checksum of the specbins

    Parameters
    ----------
    rmf: string or HDUlist or array
      The response matrix file or energy bins (see raw k/w)
    arf: string or HDUlist
      The ancillary response file
    raw : bool
      If true, the rmf variable contains the energy bin edges (keV) for all
      the bins, and each bin has a perfect response. This is effectively
      a dummy response

    Returns
    -------
    none

    """

    if raw==True:
      # make a diagonal perfect response
      self.specbins = rmf
      self.ebins_out = rmf
      self.specbin_units='keV'

      self.rmfmatrix = numpy.zeros([len(rmf)-1, len(rmf)-1])
      for i in range(len(rmf)-1):
        self.rmfmatrix[i,i] = 1.0
      self.aeff = self.rmfmatrix.sum(1)

      self.response_set = True
      self.specbins_set = True
      self.arf = False
      self.ebins_checksum =hashlib.md5(self.specbins).hexdigest()
    else:

      if util.keyword_check(arf):
        if type(arf)==str:
          self.arffile = arf
          self.arf = pyfits.open(arf)
        elif type(arf) == pyfits.hdu.hdulist.HDUList:
          self.arf = arf
          self.arffile = arf.filename()
        else:
          print("ERROR: unknown arf type, %s"%(repr(type(arf))))
          return
      else:
        self.arf=False


      if type(rmf)==str:
        self.rmffile = rmf
        self.rmf = pyfits.open(rmf)
      elif type(rmf) == pyfits.hdu.hdulist.HDUList:
        self.rmf = rmf
        self.rmffile = rmf.filename()
      else:
        print("ERROR: unknown rmf type, %s"%(repr(type(rmf))))
        return

      # get the rmf matrix

    # alternate where we do matrix generation?

    # these are the *output* energy bins
      ebins = self.rmf['EBOUNDS'].data['E_MIN']
      if ebins[-1] > ebins[0]:
        ebins = numpy.append(ebins, self.rmf['EBOUNDS'].data['E_MAX'][-1])
      else:

        ebins = numpy.append(self.rmf['EBOUNDS'].data['E_MAX'][0],ebins)
      # find the name of the rmf matrix HDU
      try:
        k=self.rmf.index_of('MATRIX')
        matrixname = 'MATRIX'
      except KeyError:
        try:
          k=self.rmf.index_of('SPECRESP MATRIX')
          matrixname = 'SPECRESP MATRIX'
        except KeyError:
          print("Cannot find index for matrix in this data")
          raise

    # bugfix: not all missions index from 0 (or 1).
    # Use chanoffset to correct for this.
      chanoffset = self.rmf['EBOUNDS'].data['CHANNEL'][0]

      self.rmfmatrix = numpy.zeros([len(self.rmf[matrixname].data),len(self.rmf['EBOUNDS'].data)])
      for ibin, i in enumerate(self.rmf[matrixname].data):

        lobound = 0

        fchan = i['F_CHAN']*1
        nchan = i['N_CHAN']*1

        if numpy.isscalar(fchan):
          fchan = numpy.array([fchan])
        fchan -= chanoffset
        if numpy.isscalar(nchan):
          nchan = numpy.array([nchan])

        for j in range(len(fchan)):
          ilo = fchan[j]
          if ilo < 0: continue

          ihi = fchan[j] + nchan[j]
          self.rmfmatrix[ibin,ilo:ihi]=i['MATRIX'][lobound:lobound+nchan[j]]
          lobound = lobound+nchan[j]

      self.specbins, self.ebins_out = get_response_ebins(self.rmf)
      self.specbin_units='keV'
      self.aeff = self.rmfmatrix.sum(1)
      if self.arf != False:
        self.aeff *=self.arf['SPECRESP'].data['SPECRESP']
      self.response_set = True
      self.specbins_set = True

      self.ebins_checksum =hashlib.md5(self.specbins).hexdigest()


  def return_spectrum(self, te, teunit='keV', nearest=False,\
                      get_nearest_t=False):
    """
    Get the spectrum at an exact temperature.
    Interpolates between 2 neighbouring spectra

    Finds HDU with kT closest to desired kT in given line or coco file.

    Opens the line or coco file, and looks for the header unit
    with temperature closest to te. Use result as index input to make_spectrum

    Parameters
    ----------
    te : float
      Temperature in keV or K
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    raw : bool
      If set, return the spectrum without response applied. Default False.
    nearest : bool
      If set, return the spectrum from the nearest tabulated temperature
      in the file, without interpolation
    get_nearest_t : bool
      If set, and `nearest` set, return the nearest tabulated temperature
      as well as the spectrum.

    Returns
    -------
    spectrum : array(float)
      The spectrum in photons cm^5 s^-1 bin^-1, with the response, or
      photons cm^3 s^-1 bin^-1 if raw is set.
    nearest_T : float, optional
      If `nearest` is set, return the actual temperature this corresponds to.
      Units are same as `teunit`
    """

    # Check that there is a response set
    if not self.response_set:
      raise util.ReadyError("Response not yet set: use set_response to set.")

    # make element and abundance lists
    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]


    self.spectra.ebins = self.specbins
    self.spectra.ebins_checksum=hashlib.md5(self.spectra.ebins).hexdigest()
    s= self.spectra.return_spectrum(te, teunit=teunit, nearest=nearest,elements = el_list, abundances=ab)

    ss = self.apply_response(s)

    return ss


  def apply_response(self, spectrum):

    """
    Apply a response to a spectrum

    Parameters
    ----------
    spectrum : array(float)
      The spectrum, in counts/bin/second, to have the response applied to. Must be
      binned on the same grid as the rmf.
    Returns
    -------
    array(float)
      energy grid (keV) for returned spectrum
    array(float)
      spectrum folded through the response
    """

    arfdat = self.arf

    if arfdat:
      res = spectrum * arfdat['SPECRESP'].data['SPECRESP']
    else:
      res = spectrum*1.0

    try:
      ret = numpy.matmul(res,self.rmfmatrix)
    except ValueError:
      ret = numpy.matmul(res,self.rmfmatrix.transpose())

    return ret







    # if nearest:
      # index = numpy.argmin(numpy.abs(self.linedata[1].data['kT']-teval))+2
      # if not (index in list(self.spectra.keys())):
        # self.spectra[index] = NewSpec(self, index)
        # self.spectra[index].calc_spectrum(teval/const.KBOLTZ)
      # te_nearest = self.linedata[1].data['kT'][index-2]
      # if teunit.lower()=='kev':
        # pass
      # elif teunit.lower() == 'ev':
        # te_nearest /=1000
      # elif teunit.lower() == 'k':
        # te_nearest = te * const.KBOLTZ

      # if not self.response_set:
        # raw=True

# # edit 8/30/19 ARF
# # old code
# #      if raw:
# #        s = self.spectra[index].spectrum
# #      else:
# #        s = self.spectra[index].spectrum_withresp
# #      if util.keyword_check(get_nearest_t):
# #        return s, te_nearest
# #      else:
# #        return s
# # new code
      # s = self.spectra[index].spectrum
# # end new code

    # else:
      # if ((teval > self.linedata[1].data['kT'][-1]) |\
          # (teval < self.linedata[1].data['kT'][0])):
        # print("*** ERROR: temperature %f keV is out of range %f-%f ***" %\
              # (teval, self.linedata[1].data['kT'][0], self.linedata[1].data['kT'][-1]))
        # return
    # # find the 2 nearest temperatures

      # index = numpy.where(self.linedata[1].data['kT'] > teval)[0][0]


      # loind = index+1
      # upind = index+2

      # # get the spectra at these temperatures
      # if not (loind in list(self.spectra.keys())):
        # self.spectra[loind] = NewSpec(self, loind)
        # self.spectra[loind].calc_spectrum(teval,thermal_broadening, self.broaden_limit, velocity_broadening)
      # if not (upind in list(self.spectra.keys())):
        # self.spectra[upind] = NewSpec(self, upind)
        # self.spectra[upind].calc_spectrum(teval,thermal_broadening, self.broaden_limit, velocity_broadening)

      # # now sum the spectra and add as a response

      # t1 = self.linedata[1].data['kT'][loind-2]
      # t2 = self.linedata[1].data['kT'][upind-2]

# # edit 8/30/19 ARF
# # old code
      # # if not self.response_set:
        # # raw=True

      # # if raw:
        # # s1 = self.spectra[loind].spectrum
        # # s2 = self.spectra[upind].spectrum

      # # else:
        # # s1 = self.spectra[loind].spectrum_withresp
        # # s2 = self.spectra[upind].spectrum_withresp
# # new code
      # s1 = self.spectra[loind].spectrum
      # s2 = self.spectra[upind].spectrum
# # end new code

      # # linear interp

      # r1 = 1- (teval-t1)/(t2-t1)
      # r2 = 1- r1

      # s = s1*r1 + s2*r2

# edit 8/30/19 ARF
# new code
    # do the response in here

    if ((raw==True) | (self.response_set==False)):
      return s
    else:
      s = apply_response(s, self.rmf, arf=self.arf)
      return s

  # def set_specbins(self, specbins, specunits='A'):
    # """
    # Set the energy or wavelength bin for the raw spectrum

    # Note that this is overridden if a response is loaded

    # Parameters
    # ----------
    # ebins : array(float)
      # The edges of the spectral bins (for n bins, have n+1 edges)
    # specunits : {'a','kev'}
      # The spectral bin units to use. Default is angstroms

    # Returns
    # -------
    # None

    # Notes
    # -----
    # updates  self.specbins, self.binunits, self.specbins_set
    # """

    # # set the energy bins for this spectrum
    # self.specbins=specbins
    # self.binunits=specunits
    # self.specbins_set=True

    # # reset the spectrum
    # self.spectra={}




  def set_apec_files(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits"):
    """
    Set the apec line and coco files, and load up their data

    Parameters
    ----------
    linefile : str or HDUList
      The filename of the line emissivity data, or the opened file.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the opened file.

    Returns
    -------
    None

    Notes
    -----
    Updates self.linefile, self.linedata, self.cocofile and self.cocodata
    """
    if util.keyword_check(linefile):
      if isinstance(linefile, str):
        lfile = os.path.expandvars(linefile)
        if not os.path.isfile(lfile):
          print("*** ERROR: no such file %s. Exiting ***" %(lfile))
          return -1
        self.linedata = pyfits.open(lfile)
        self.linefile = lfile

      elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.linedata=linefile
        self.linefile=linefile.filename()

      else:
        print("Unknown data type for linefile. Please pass a string or an HDUList")

    if util.keyword_check(cocofile):

      if isinstance(cocofile, str):

        cfile = os.path.expandvars(cocofile)
        if not os.path.isfile(cfile):
          print("*** ERROR: no such file %s. Exiting ***" %(cfile))
          return -1
        self.cocodata=pyfits.open(cfile)
        self.cocofile=cfile

      elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.cocodata=cocofile
        self.cocofile=cocofile.filename()

      else:
        print("Unknown data type for cocofile. Please pass a string or an HDUList")

  def set_abund(self, elements, abund):
    """
    Set the elemental abundance, relative to the abundset. Defaults to
    1.0 for everything

    Parameters
    ----------
    elements : int or array_like(int)
      The elements to change the abundance of
    abund : float or array_like(float)
      The new abundances. If only 1 value, set all `elements` to this abundance
      Otherwise, should be of same length as elements.

    Returns
    -------
    None

    Examples
    --------
    Set the abundance of iron to 0.5

    >>> myspec.set_abund(26, 0.5)

    Set the abundance of iron and nickel to 0.1 and 0.2 respectively

    >>> myspec.set_abund([26, 28], [0.1,0.2])

    Set the abundance of oxygen, neon, magnesium and iron to 0.1

    >>> myspec.set_abund([8,10,12,26],0.1)
    """

    abundvec, aisvec = util.make_vec(abund)
    elementvec, eisvec = util.make_vec(elements)
    if (aisvec):
      if len(abundvec)!= len(elementvec):

        print("abundance vector and element vector must have same number"+\
              " of elements")
      else:

        self.abund[elementvec] = abundvec
    elif (eisvec):
      # set all these eleemtns to the same abundance
      for el in elementvec:
        self.abund[el]=abund

    else:
      self.abund[elements]=abund



  def set_abundset(self, abundstring):
    """
    Set the abundance set.

    Parameters
    ----------
    abundstring : string
      The abundance string (e.g. "AG89", "uniform". Case insensitive.
      See atomdb.get_abundance for list of possible abundances

    Returns
    -------
    none
      updates self.abundset and self.abundsetvector.
    """

    # read in the abundance the raw data was calculated on
    old = atomdb.get_abundance(abundset=self.default_abundset,\
                               datacache=self.datacache)



    # read in the new abundance
    new = atomdb.get_abundance(abundset=abundstring,\
                               datacache=self.datacache)

    # divide the 2, store the replacement ratio to self.abundsetvector

    for Z in range(1,const.MAXZ_CIE+1):
      try:
        self.abundsetvector[Z]=new[Z]/old[Z]
      except ZeroDivisionError:
        self.abundsetvector[Z] = 0.0
      except IndexError:
        print("spectrum IndexError for Z = %i"%(Z))



    # update the current abundance string to represent your input
    self.abundset=abundstring

    #self.recalc()


  def return_line_emissivity(self, Telist, Z, z1, up, lo, \
                             specunit='A', teunit='keV', \
                             apply_aeff=False, apply_abund=True,\
                             log_interp = True):
    """
    Get line emissivity as function of Te.


    Parameters
    ----------
    Telist : float or array(float)
      Temperature(s) in keV or K
    Z : int
      nuclear charge of element
    z1 : int
      ion charge +1 of ion
    up : int
      upper level for transition
    lo : int
      lower level for transition
    specunit : {'Angstrom','keV'}
      Units for wavelength or energy (a returned value)
    teunit : {'keV' , 'K'}
      Units of Telist (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the line emissivity in the
      linelist to modify their intensities.
    apply_abund : bool
      If true, apply the abundance set in the session to the result.
    log_interp : bool
      Interpolate between temperature on a log-log scale (default).
      Otherwise linear

    Returns
    -------
    ret : dict
      Dictionary containing:
      Te, tau, teunit: as input
      wavelength : line wavelength (A)
      energy : line energy (keV)
      epsilon : emissivity in ph cm^3 s-1 (or ph cm^5 s^-1 if apply_aeff=True)
                first index is temperature, second is tau.

    """

    Tevec, Teisvec = util.make_vec(Telist)

    kTlist = convert_temp(Tevec, teunit, 'keV')
    if apply_abund:
      ab = self.abund[Z]*self.abundsetvector[Z]
    else:
      ab = 1.0

    eps = numpy.zeros(len(Tevec))
    ret={}

    for ikT, kT in enumerate(kTlist):
      e, lam = self.spectra.return_line_emissivity(kT, Z, z1, \
                                                   up, lo, \
                                                   specunit='A', \
                                                   teunit='keV', \
                                                   abundance=ab)
      eps[ikT] = e
      if lam != False:
        ret['wavelength'] = lam * 1.0
      else:
        ret['wavelength'] = None

    ret['Te'] = Telist
    ret['teunit'] = teunit
    if ret['wavelength'] != None:
      ret['energy'] = const.HC_IN_KEV_A/ret['wavelength']
    else:
      ret['energy'] = None


    if apply_aeff == True:
      e = ret['energy']
      ibin = numpy.where(self.specbins<e)[0][-1]

      eps = eps*self.aeff[ibin]

    # now correct for vectors

    if not Teisvec:
      eps = eps[0]

    ret['epsilon'] = eps

    return ret


  def return_linelist(self, Te, specrange, specunit='A', \
                               teunit='keV', apply_aeff=False):
    """
    Get the list of line emissivities vs wavelengths


    Parameters
    ----------
    Te : float
      Temperature in keV or K
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Ansgtrom','keV'}
      Units for specrange
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the lines in the linelist to
      modify their intensities.

    Returns
    -------
    linelist : array(dtype)
      The list of lines with lambda (A), energy (keV), epsilon (ph cm3 s-1),\
      epsilon_aeff (ph cm5 s-1) ion (string) and upper & lower levels.

    """
    kT = convert_temp(Te, teunit, 'keV')


    # make element and abundance lists
    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]


    s= self.spectra.return_linelist(kT, specrange=specrange, teunit='keV',\
                                        specunit=specunit, elements=el_list,\
                                        abundances = ab)

    # do the response thing
    #resp  = s.response()

    if apply_aeff == True:
      ibin = numpy.zeros(len(s), dtype=int)
      for i, ss in enumerate(s):
        e = const.HC_IN_KEV_A/ss['Lambda']
        ibin[i] = numpy.where(self.specbins<e)[0][-1]

      s["Epsilon_Err"] = s['Epsilon']*self.aeff[ibin]

    return(s)

def convert_temp(Te, teunit, teunitout):
  """
  Convert temperature (Te) from units teunit to teunitout

  Parameters
  ----------
  Te : float
    The temperature
  teunit : string
    units of Te
  teunitout : string
    output temperature units
  """

  teunit2 = teunit.lower()
  teunitout2 = teunitout.lower()

  if teunitout2==teunit2: return Te

  allowedunits = ['ev','kev','k']

  cfactors = [1000.,1.0,1/const.KBOLTZ]
  try:
    cfac = cfactors[allowedunits.index(teunitout2)]/cfactors[allowedunits.index(teunit2)]
  except ValueError:
    # Assume one of the units is bad
    if not (teunit2 in allowedunits):
      raise util.UnitsError("%s is not a recognized temperature unit %s"%(teunit2, allowedunits))
    elif not (teunitout2 in allowedunits):
      raise util.UnitsError("%s is not a recognized temperature unit %s"%(teunitout2, allowedunits))
    raise
  return cfac*Te


def convert_spec(spec, specunit, specunitout):
  """
  Convert spectral ranges from specunit to specunitout

  Parameters
  ----------
  spec : array
    The units to return
  specunit : string
    The input spectral unit ('keV', 'A')
  specunitout : string
    The output spectral unit ('keV', 'A')

  Returns
  -------
  specout : array
    spec, converted to specunitout
  """

  allowedunits = ['a','ang','angstrom','angstroms','ev','kev']
  cfactors = [1.0, 1.0, 1.0, 1.0, const.HC_IN_KEV_A*1000, const.HC_IN_KEV_A]
  ctype    =  ['w','w','w','w','e','e']

  specunit2 = specunit.lower()
  specunitout2 = specunitout.lower()

  # If units are the same, do nothing
  if specunit2==specunitout2: return spec

  try:
    cfac_in = cfactors[allowedunits.index(specunit2)]
    cfac_out = cfactors[allowedunits.index(specunitout2)]



  except ValueError:
    # Assume one of the units is bad
    if not (specunit2 in allowedunits):
      raise util.UnitsError("%s is not a recognized spectroscopic unit %s"%\
           (specunit2, allowedunits))
    elif not (specunitout2 in allowedunits):
      raise util.UnitsError("%s is not a recognized spectroscopic unit %s"%\
           (specunitout2, allowedunits))
    raise
  ctype_in = ctype[allowedunits.index(specunit2)]
  ctype_out = ctype[allowedunits.index(specunitout2)]


  if ctype_in == 'w':
    try:
      spec_ang = cfac_in * spec
    except TypeError:
      spec = numpy.array(spec)
      spec_ang = cfac_in * spec
  elif ctype_in == 'e':
    try:
      spec_ang = cfac_in / spec
    except TypeError:
      spec = numpy.array(spec)
      spec_ang = cfac_in / spec
    spec_ang = cfac_in / spec


  if ctype_out == 'w':
    spec_out = spec_ang * cfac_out
  elif ctype_out == 'e':
    spec_out = cfac_out/spec_ang

  if ctype_out != ctype_in:
    # invert array if converting between energy and wavelength
    spec_out = spec_out[::-1]

  return spec_out






class CIESpectrum():
  """
  A class holding the emissivity data for CIE emission, and returning
  spectra

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : arraylike(int), optional
    The atomic number of elements to include (default all)
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  session : CIESession
    The parent CIESession
  SessionType : string
    "CIE"
  spectra : dict of ElementSpectra
    a dictionary containing the emissivity data for each HDU,
    subdivided by element (spectra[12][18] is an ElementSpectrum object
    containing the argon data for the 12th HDU)
  kTlist : array
    The temperatures for each emissivity HDU, in keV
  logkTlist : array
    log of kTlist
  """

  def __init__(self, linedata, cocodata):
    """
    Initializes the code. Populates the line and emissivity data in all
    temperature HDUs.

    Parameters
    ----------
    """

    self.SessionType = 'CIE'



    # HERE RELOAD ALL THE SPECTRAL DATA OF FILL IT OUT, YOUR CALL
    picklefname = os.path.expandvars('$ATOMDB/spectra_%s_%s.pkl'%\
                                (linedata[0].header['CHECKSUM'],\
                                 cocodata[0].header['CHECKSUM']))
    if os.path.isfile(picklefname):
      self.spectra = pickle.load(open(picklefname,'rb'))
      self.kTlist = self.spectra['kTlist']
    else:
    # < insert test for existing pkl file for spectra, otherwise populate>
      self.spectra={}
      self.kTlist = numpy.array(linedata[1].data['kT'].data)
      self.spectra['kTlist']=numpy.array(linedata[1].data['kT'].data)

      for ihdu in range(len(self.kTlist)):
        self.spectra[ihdu]={}
        self.spectra[ihdu]['kT'] = self.kTlist[ihdu]
        ldat = numpy.array(linedata[ihdu+2].data.data)
        cdat = numpy.array(cocodata[ihdu+2].data.data)

        Zarr = numpy.zeros([len(ldat), const.MAXZ_CIE+1], dtype=bool)
        Zarr[numpy.arange(len(ldat), dtype=int), ldat['Element']]=True


        for Z in range(1,const.MAXZ_CIE+1):
          ccdat = cdat[(cdat['Z']==Z) & (cdat['rmJ']==0)]

          if len(ccdat)==1:
            c = ccdat[0]
          else:
            c = False

          self.spectra[ihdu][Z]=ElementSpectrum(ldat[Zarr[:,Z]],\
                                              c, \
                                              Z)
      pickle.dump(self.spectra, open(picklefname,'wb'))


    self.logkTlist=numpy.log(self.kTlist)


  def get_nearest_Tindex(self, Te, teunit='keV', nearest=False,\
                         log_interp=True):
    """
    Return the nearest temperature index in the emissivity file, or,
    alternatively, the array of fractions to sum

    Parameters
    ----------
    Te : float
      Temperature (keV by default)
    teunit : string
      Units of kT
    nearest : bool
      If true, return only nearest. Otherwise, return nearest 2 and
      fractions

    Returns
    -------
    ikT : list[int]
      Index of temperature in HDU file (from 0, not 2)
    f : list[float]
      fractional weight to apply to each ikT. Should sum to 1.
    """
    kT = convert_temp(Te, teunit, 'keV')

    # find the nearest temperature
    if kT < self.kTlist[0]:
      print("kT = %f is below minimum range of %f. Returning lowest kT spectrum available"%\
             (kT, self.kTlist[0]))
      ikT = [0]
      f=[1.0]

    elif kT > self.kTlist[-1]:
      print("kT = %f is above maximum range of %f. Returning highest kT spectrum available"%\
             (kT, self.kTlist[-1]))
      ikT = [len(self.kTlist)-1]
      f=[1.0]
    else:

      if log_interp:
        if nearest:
          ikT = [numpy.argmin(numpy.abs(self.logkTlist - numpy.log(kT)))]
          f=[1.0]
        else:
          ikT = numpy.where(self.kTlist < kT)[0][-1]
          ikT = [ikT, ikT+1]

          f = 1- (numpy.log(kT)-self.logkTlist[ikT[0]])/\
                 (self.logkTlist[ikT[1]]-self.logkTlist[ikT[0]])
          f = [f,1-f]

      else:
        if nearest:
          ikT = [numpy.argmin(numpy.abs(self.kTlist - kT))]
          f=[1.0]
        else:
          ikT = numpy.where(self.kTlist < kT)[0][-1]
          ikT = [ikT, ikT+1]

          f = 1- (kT-self.kTlist[ikT[0]])/\
                 (self.kTlist[ikT[1]]-self.kTlist[ikT[0]])
          f = [f,1-f]

    return ikT, f

#-----------------------------------------------------------------------

  def return_spectrum(self, Te, teunit='keV', nearest = False,
                             elements=False, abundances=False, log_interp=True):

    """
    Return the spectrum of the element on the energy bins in
    self.session.specbins

    Parameters
    ----------
    Te : float
      Electron temperature (default, keV)
    teunit : string
      Units of kT (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.

    Returns
    -------
    spec : array(float)
      The element's emissivity spectrum, in photons cm^3 s^-1 bin^-1
    """

    # get kT in keV
    kT = convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_CIE+1)


    if abundances == False:
      abundances = {}
      for Z in elements:
        abundances[Z] = 1.0

    s = 0.0

    for Z in elements:
      abund = abundances[Z]
      if abund > 0:
        epslimit =  self.broaden_limit/abund

           # go caclulate the spectrum, with broadening as assigned.
        sss=0.0
        print("ikT", ikT, 'f', f)
        for i in range(len(ikT)):


          ss = self.spectra[ikT[i]][Z].return_spectrum(self.ebins,\
                                  kT,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening) *\
                                  abund

          if log_interp:
            sss += numpy.log(ss+const.MINEPSOFFSET) *f[i]
          else:
            sss += ss*f[i]
        # get the interp
        if log_interp:
          s += numpy.exp(sss)-const.MINEPSOFFSET*len(ikT)
        else:
          s += sss

    return s



  def return_linelist(self, Te, teunit='keV', nearest = False,\
                      specrange=False, specunit='A', elements=False, abundances=False):

    """
    Return the linelist of the element

    Parameters
    ----------

    Te : float
      Electron temperature (default, keV)
    teunit : string
      Units of kT (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Ansgtrom','keV'}
      Units for specrange (default A)


    """

    # get kT in keV
    kT = convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_CIE+1)


    if abundances == False:
      abundances = {}
      for Z in elements:
        abundances[Z] = 1.0


    linelist = False

    for Z in elements:
      abund = abundances[Z]
      if abund > 0:
        elemlinelist = False
        for i in range(len(ikT)):


          ss = self.spectra[ikT[i]][Z].return_linelist(specrange,\
                                  teunit='keV', specunit=specunit)

          if len(ss) > 0:
            ss['Epsilon']*=abund*f[i]
            if type(elemlinelist)==bool:
              if elemlinelist==False:
                elemlinelist = ss
            else:
              isnew = numpy.zeros(len(ss), dtype=bool)

              for inew, new in enumerate(ss):
                imatch = numpy.where((new['Element']==elemlinelist['Element']) &\
                                     (new['Ion']==elemlinelist['Ion']) &\
                                     (new['UpperLev']==elemlinelist['UpperLev']) &\
                                     (new['LowerLev']==elemlinelist['LowerLev']))[0]
                if len(imatch)==1:
                  elemlinelist[imatch[0]]['Epsilon']+=new['Epsilon']
                else:
                  isnew[inew]=True

              s = sum(isnew)
              if s > 0:
                elemlinelist = numpy.append(elemlinelist, ss[isnew])
        if elemlinelist != False:
          if linelist==False:
            linelist = elemlinelist
          else:
            linelist =  numpy.append(linelist, elemlinelist)

    return linelist


# class NewSpec():
    # """
    # An individual spectrum, from a specifically
    # tabulated temperature in a line/coco file.

    # Attributes
    # ----------
    # temperature : float
      # The temperature of this spectrum, in keV
    # index : int
      # The index in the line file for this spectrum
    # """



    # def __init__(self, session, index):
      # self.temperature = session.linedata[1].data['kT'][index-2]
      # self.index = index
      # self.spectra={}
      # self.session=session
      # self.linelist=False



    # def set_index(T, teunit='K', logscale = False):
      # """
      # Finds HDU with kT closest to desired kT in given line or coco file.

      # Opens the line or coco file, and looks for the header unit
      # with temperature closest to te. Use result as index input to make_spectrum

      # Parameters
      # ----------
      # te : float
        # Temperature in keV or K
      # teunits : {'keV' , 'K', 'eV'}
        # Units of te (kev or K, default keV)
      # logscale : bool
        # Search on a log scale for nearest temperature if set.

      # Returns
      # -------
      # none

      # Notes
      # -----
      # modifies
      # self.index : int
      # Index in HDU file with nearest temperature to te.

      # """

      # if teunit.lower() == 'kev':
        # teval = te
      # elif teunit.lower() == 'ev':
        # teval = te/1000.0
      # elif teunit.lower() == 'k':
        # teval = te*const.KBOLTZ
      # else:
        # print("*** ERROR: unknown temeprature unit %s. Must be eV, keV or K. Exiting ***"%\
              # (teunits))

      # if logscale:
        # i = numpy.argmin(numpy.abs(numpy.log(self.linedata[1].data['kT'])-numpy.log(teval)))
      # else:
        # i = numpy.argmin(numpy.abs(self.linedata[1].data['kT']-teval))
      # # need to increase the HDU by 2.
      # self.index = i+2

    # def calc_spectrum(self, T, themalbroadening=False, broaden_limit=1e-18, velocity_broadening=0.0):

      # """
      # Calculates the spectrum for each element on a single temperature

      # Parameters
      # ----------
      # session : Session
        # The parent Session
      # dolines : bool
        # Include lines in the spectrum
      # docont : bool
        # Include continuum in the spectrum
      # dopseudo : bool
        # Include pseudocontinuum in the spectrum
      # Outputs
      # -------
      # none

      # Notes
      # -----
      # Modifies:\n
      # dict : self.spectrum_by_Z  the spectrum of each element\n

      # Then calls `recalc()` to update the spectra
      # """
      # # now, we shall calculate the spectrum for each individual element

      # # set the linefile


      # self.temperature = self.session.linedata[1].data['kT'][self.index-2]

  # #    if util.keyword_check(elements):
        # #self.set_abund(elements, abund)

      # #if util.keyword_check(abund):
        # #self.set_abund(elements, abund)

      # # to hold the total spectrum
      # spec = numpy.zeros(len(self.session.specbins)-1)
      # ebins_checksum = hashlib.md5(self.session.specbins).hexdigest()
      # t0 = time.time()

      # for Z in self.session.elements:
      # # make the generic spectrum
        # if not Z in self.spectra.keys():
          # self.spectra[Z]=ElementSpectrum(self.session.linedata[self.index].data, \
                                                   # self.session.cocodata[self.index].data,\
                                                   # Z, parent=self)

        # # calculate the element's spectrum
        # self.spectra[Z].calc_spectrum(self.session.specbins, \
                                                       # T*const.KBOLTZ, \
                                                       # ebins_checksum=ebins_checksum,\
                                                       # thermal_broadening=thermal_broadening,\
                                                       # broaden_limit=broaden_limit,\
                                                       # velocity_broadening = velocity_broadening)



      # t1 = time.time()
      # if self.session.response_set:
        # for Z in self.session.elements:
          # self.spectra[Z].apply_response()
      # t2 = time.time()

      # print('Time to calculate spectra: %fs'%(t1-t0))
      # print('Time to apply response: %fs'%(t2-t1))


      # self.recalc()


    # def calc_line_emissivities(self, apply_aeff=False):

      # """
      # Calculates the line emissivities

      # Parameters
      # ----------
      # apply_aeff : bool
        # apply effective area to the results
      # Outputs
      # -------
      # none

      # Notes
      # -----
      # Modifies:\n
      # dict : self.spectrum_by_Z  the spectrum of each element\n
      # Then calls `recalc()` to update the spectra
      # """
      # # now, we shall calculate the spectrum for each individual element

      # # set the linefile


      # self.temperature = self.session.linedata[1].data['kT'][self.index-2]

  # #    if util.keyword_check(elements):
        # #self.set_abund(elements, abund)

      # #if util.keyword_check(abund):
        # #self.set_abund(elements, abund)

      # # to hold the list of lines
      # if self.linelist == False:
        # self.linelist = numpy.zeros(len(self.session.linedata[self.index].data), dtype=\
                                # numpy.dtype({'names':   ['lambda','energy','epsilon','epsilon_aeff','ionsymb','upperlev','lowerlev'],\
                                             # 'formats': [float, float, float, float, '|S10', int, int]}))

        # iline = 0

        # # if required, create the ElementSpectrum objects for each case
        # ebins_checksum = hashlib.md5(self.session.specbins).hexdigest()
        # for Z in self.session.elements:
          # if not Z in self.spectra.keys():
            # self.spectra[Z]=ElementSpectrum(self.session.linedata[self.index].data, \
                                                     # self.session.cocodata[self.index].data,\
                                                     # Z, parent=self)

        # # calculate the element's spectrum
          # nlines = len(self.spectra[Z].lines.lines)
          # self.linelist[iline:iline+nlines]['lambda'] = self.spectra[Z].lines.lines['lambda']
          # self.linelist[iline:iline+nlines]['energy'] = const.HC_IN_KEV_A/self.spectra[Z].lines.lines['lambda']
          # self.linelist[iline:iline+nlines]['epsilon'] = self.spectra[Z].lines.lines['epsilon']
          # self.linelist[iline:iline+nlines]['upperlev'] = self.spectra[Z].lines.lines['upperlev']
          # self.linelist[iline:iline+nlines]['lowerlev'] = self.spectra[Z].lines.lines['lowerlev']
          # for ii in range(iline, iline+nlines):
            # self.linelist[ii]['ionsymb'] = atomic.spectroscopic_name(Z, self.spectra[Z].lines.lines['ion'][ii-iline])
          # iline += nlines
      # nlines = len(self.linelist)

      # # apply the effective area

      # if apply_aeff:

        # if self.session.response_set:
          # self.linelist['epsilon_aeff'] = numpy.zeros(len(self.linelist))

          # igood = (self.linelist['energy'] >= self.session.specbins[0]) &\
                  # (self.linelist['energy'] < self.session.specbins[-1])


          # aeffindex = numpy.digitize(self.linelist[igood]['energy'], self.session.specbins)-1
          # zzz = self.linelist[igood]['epsilon']*self.session.aeff[aeffindex]
          # self.linelist['epsilon_aeff'][igood]= zzz

    # def recalc(self):
      # """
      # Recalculate the spectrum - just for changing abundances etc.
      # Does not recalculate spectrum fully, just changes the multipliers.
      # Does nothing if self.ready is False, should be run after calc_spectrum.

      # Parameters
      # ----------
      # none

      # Returns
      # -------
      # none

      # Notes
      # -----
      # Modifies:\n
      # self.spectrum : array_like (float)\n
      # """
      # print("recalc: ", self.session.ready, self.session.specbins_set)
      # if self.session.ready:
        # if self.session.specbins_set:

          # self.spectrum = numpy.zeros(len(self.session.specbins)-1)
          # for Z in self.session.elements:
            # self.spectrum += self.spectra[Z].spectrum * self.session.abund[Z] * self.session.abundsetvector[Z]




class ElementSpectrum():
  """
  A class holding the emissivity data for an element in one HDU

  Parameters
  ----------
  linedata : array(linedatatype)
    array of line wavelengths and emissivities, from AtomDB files.
    Should already be filtered to only be from one element.
  cocodata : array(cocodatatype)
    array of continuum wavelengths and emissivities, from AtomDB files
    Should already be filtered to only be from one element.
  Z : int
    The atomic number of the element
  z1_drv : int
    The charge + 1 for the ion. 0 = whole element.
  parent : CIESpectrum
    Parent CIESpectrum object

  Attributes
  ----------
  lines : LineData
    A LineData object containing all the line information
  continuum : ContinuumData
    A ContinuumData object containing all the contrinuum information
  parent : CIESpectrum
    Parent CIESpectrum object
  session : CIESession
    Parent Session of parent CIESpectrum object
  """

  def __init__(self, linedata, cocodata, Z, z1_drv=0):
    """
    Initialization

    Parameters
    ----------
    linedata : hdu
      the line data passed in
    cocodata : hdu
      the continuum data passed in
    Z : int
      the atomic number of the element
    z1 : int
      the ion charge of the driving ion (0 for whole element)
    """

    # intialize
    if z1_drv != 0:
      tmp = linedata[(linedata['Element'] == Z) &\
                               (linedata['Ion_drv'] == z1_drv)]
      self.lines = LineData(tmp)

#      tmp = cocodata[(cocodata['Z']==Z) &\
#                     (cocodata['rmJ']==z1_drv)]
#      print(tmp)
      self.continuum = ContinuumData(cocodata)
    else:
#      tmp = linedata[(linedata['Element'] == Z)]
      self.lines = LineData(linedata)

#      tmp = cocodata[(cocodata['Z']==Z) &\
                     #(cocodata['rmJ']==0)]

#      self.continuum = ContinuumData(tmp[0], parentElementSpectrum=self)
      self.continuum = ContinuumData(cocodata)


  def return_spectrum(self, eedges, Te, ebins_checksum=False,\
                    thermal_broadening=False,\
                    broaden_limit=False,\
                    velocity_broadening=0.0,\
                    teunit = 'keV'):
    """
    Calculate the spectrum

    Parameters
    ----------
    eedges : array
      bin edges (keV)
    Te : float
      temperature (keV by defult)
    ebins_checksum : string
      the md5 checksum of eedges
    thermal_broadening : bool
      true to apply thermal broadening
    broaden_limit : float
      only broaden lines stronger than this.
    velocity_broadening : float
      velocity broadening to apply, km/s. Set <=0 for none (default)
    teunit : string
      Temperature unit (K, keV, eV)

    Returns
    -------
    spectrum : array(float)
      The spectrum in ph cm^3 s^-1 bin^-1
    """

    T = convert_temp(Te, teunit, 'K')

    if ebins_checksum == False:
      # check the parent
#      if self.parentElementSpectrum != False:
#        ebins_checksum = self.parentElementSpectrum.ebins_checksum

      # check again, in case there was no parent
#      if ebins_checksum == False:
        # generate the checksum
      ebins_checksum = hashlib.md5(eedges).hexdigest()

    self.ebins_checksum = ebins_checksum
    self.T = T

    spec = self.lines.return_spec(eedges, T, ebins_checksum=ebins_checksum,\
                                  thermal_broadening=thermal_broadening,\
                                  broaden_limit=broaden_limit,\
                                  velocity_broadening=velocity_broadening) +\
           self.continuum.return_spec(eedges, ebins_checksum = ebins_checksum)

    self.spectrum = spec

    return self.spectrum

  def return_linelist(self, specrange, specunit='A',
                      teunit = 'keV'):
    """
    Return the list of lines in specrange

    Parameters
    ----------
    specrange : array
      spectral range to look for lines in
    specunits : string
      units of spectral range (A or keV)
    teunit : string
      units of Te (keV, eV, K)

    Returns
    -------
    linelist : array
      list of lines and epsilons
    """
    wave = convert_spec(specrange, specunit, 'A')

    llist = self.lines.lines[(self.lines.lines['Lambda']>=wave[0]) &\
                       (self.lines.lines['Lambda']<=wave[1])]

    return llist

  def return_linematch(self, Z, z1, up, lo, z1_drv=0):
    """
    Return the line(s) which match the transition

    Parameters
    ----------
    Z : int
      nuclear charge
    z1 : int
      ion charge + 1
    up : int
      upper level of transition
    lo : int
      lower level of transition
    z1_drv : int
      if provided, also filter on driving ion charge (NEI only)
    Returns
    -------
    linelist : array
      list of lines and epsilons
    """

    print(self.lines.lines.dtype.names)
    llist = self.lines.lines[(self.lines.lines['Element']==Z) &\
                             (self.lines.lines['Ion']==z1) &\
                             (self.lines.lines['UpperLev']==up) &\
                             (self.lines.lines['LowerLev']==lo)]
    if z1_drv != 0:
      llist = llist[llist['Ion_drv']==z1_drv]


    return llist

  # def apply_response(self):
    # """
    # Apply a response to a spectrum

    # Parameters
    # ----------
    # spectrum : array(float)
      # The spectrum, in counts/bin/second, to have the response applied to. Must be
      # binned on the same grid as the rmf.
    # rmf : string or pyfits.hdu.hdulist.HDUList
      # The filename of the rmf or the opened rmf file
    # arf : string or pyfits.hdu.hdulist.HDUList
      # The filename of the arf or the opened arf file
    # Returns
    # -------
    # array(float)
      # energy grid (keV) for returned spectrum
    # array(float)
      # spectrum folded through the response
    # """

    # arfdat = self.session.arf

    # if arfdat:
      # #arfdat = arf
      # res = self.spectrum * arfdat['SPECRESP'].data['SPECRESP']
    # else:
      # res = self.spectrum*1.0

    # ret = numpy.matmul(res,self.session.rmfmatrix)

    # self.spectrum_withresp=ret

class LineData():
  """
  A class holding the line data for an element in one HDU

  Parameters
  ----------
  linelist : array(linedatatype)
    array of line wavelengths and emissivities, from AtomDB files.

  Attributes
  ----------
  lines : array(linedatatype)
    List of lines, wavelength and emissivities
  lineenergies : array(float)
    list of line energies
  spectrum_calculated : bool
    True if spectrum has already been calculated, otherwise false
  T: float
    Temperature (K) last spectrum was calculated at (for broadening)
  v : float
    Velocity (km/s) last spectrum was calculated at (for broadening)
  ebins_checksum : string
    md5sum of the ebins the spectrum was last calculated on. Used to
    identify if new calculations are required or can just return the
    previous value.
  """


  def __init__(self, linelist):
    """
    Initialization

    Parameters
    ----------
    linelist : numpy array
      list of lines from AtomDB fits files
    """

    self.lines = linelist

    self.lineenergies = const.HC_IN_KEV_A/self.lines['Lambda']
     # store the parent ElementSpectrum
    self.spectrum_calculated = False
    self.T = 0.0
    self.v = 0.0
    self.ebins_checksum = False


  def return_spec(self, eedges, T, ebins_checksum = False,\
                  thermal_broadening = False, \
                  velocity_broadening = 0.0, \
                  broaden_limit = 1e-18):
    """
    return the line emission spectrum at tempterature T

    Parameters
    ----------

    eedges : array
      energy bin edges, keV
    T : float
      temperature in Kelvin
    ebins_checksum : string
      the md5 checksum of eedges
    thermal_broadening : bool
      true to apply thermal broadening
    velocity_broadening : float
      velocity broadening to apply, km/s. Set <=0 for none (default)
    broaden_limit : float
      only broaden lines stronger than this.

    Returns
    -------
    spectrum : array(float)
      Emissivity on eedges spectral bins of the lines, in ph cm^3 s^-1 bin^-1
    """
    if ebins_checksum == False:
        # generate the checksum
      ebins_checksum = hashlib.md5(eedges).hexdigest()

    if ((thermal_broadening == False) & \
        (velocity_broadening == False)):
      if ((self.ebins_checksum==ebins_checksum) &\
          (self.spectrum_calculated == True)):
        # no need to do anything.
        pass

      else:
        spec,z = numpy.histogram(self.lineenergies, \
                                 bins=eedges, \
                                 weights = self.lines['Epsilon'])
        self.spectrum = spec
        self.spectrum_calculated = True

    else: #if thermal_broadening == True

      if ((thermal_broadening == True) &\
          (velocity_broadening <= 0)):
        # cehck to see if we need ot redo this
        if ((self.ebins_checksum==ebins_checksum) &\
            (self.spectrum_calculated == True) &\
            (self.T == T) &\
            (self.v <=0.0)):
          pass
        else:
          # recalculate!
          recalc = True

      if ((thermal_broadening == False) &\
          (velocity_broadening > 0)):
        # cehck to see if we need ot redo this
        if ((self.ebins_checksum==ebins_checksum) &\
            (self.spectrum_calculated == True) &\
            (self.T == T) &\
            (self.v == velocity_broadening)):
          pass
        else:
          # recalculate!
          recalc = True

      if ((thermal_broadening == True) &\
          (velocity_broadening > 0)):

        if ((self.ebins_checksum==ebins_checksum) &\
            (self.spectrum_calculated == True) &\
            (self.T == T) &\
            (self.v == velocity_broadening)):

          pass
        else:
          # recalculate!
          recalc = True
      if recalc==True:

        # ind = strong line indicies
        # nonind = weak line indicies
        ind = self.lines['Epsilon']>broaden_limit
        nonind = ~ind

        # calculate the widths of the strong lines
        llist = self.lines[ind]

        # get a raw dictionary of masses in amu
        masslist = atomic.Z_to_mass(1,raw=True)

        if thermal_broadening==False:
          T=0.0
          Tb = 0.0
        else:
          Tb = convert_temp(T, 'K','keV')*const.ERG_KEV/(masslist[llist['Element']]*1e3*const.AMUKG)

        if velocity_broadening <0:
          velocitybroadeining = 0.0
          vb=0.0
        else:
          vb = (velocity_broadening * 1e5)**2

        wcoeff = numpy.sqrt(Tb+vb) / (const.LIGHTSPEED*1e2)
        #width = wcoeff*const.HC_IN_KEV_A/llist['Lambda']
        elines = self.lineenergies[ind]
        width = wcoeff*elines


        # Filter out lines more than NSIGMALIMIT sigma outside the range
        NSIGMALIMIT=4
        eplu = elines+NSIGMALIMIT*width
        eneg = elines-NSIGMALIMIT*width
        emax = max(eedges)
        emin = min(eedges)
        # identify all the good lines!
        igood = numpy.where(((elines >= emin) & (eneg < emax))  |\
                  ((elines < emin) & (eplu < emin)))[0]
        spec = numpy.zeros(len(eedges))
        for iline in igood:

          spec += norm.cdf(eedges, loc=const.HC_IN_KEV_A/llist['Lambda'][iline],\
                           scale=width[iline])*llist['Epsilon'][iline]

        spec = spec[1:]-spec[:-1]


        # Then add on the weak lines
        s,z = numpy.histogram(self.lineenergies[nonind], \
                              bins = eedges,\
                              weights = self.lines['Epsilon'][nonind])
        spec+=s
        self.spectrum = spec
        self.spectrum_calculated = True

        if thermal_broadening:
          self.T = T
        else:
          self.T=0.0

        self.v = velocity_broadening
    return self.spectrum


class ContinuumData():
  """
  A class holding the continuum data for an element in one HDU

  Parameters
  ----------
  cocoentry : array(cocodatatype)
    A single row from the continuum data in an AtomDB file.
  parentElementSpectrum : ElementSpectrum
    Parent ElementSpectrum object

  Attributes
  ----------
  ECont : array(float)
    The continuum energies (keV)
  EPseudo : array(float)
    The pseudocontinuum energies (keV)
  Cont : array(float)
    The continuum emissivities (ph cm^3 s^-1 keV^-1)
  Pseudo : array(float)
    The pseudocontinuum energies (ph cm^3 s^-1 keV^-1)
  parentElementSpectrum : ElementSpectrum
    Parent ElementSpectrum object
  spectrum_calculated : bool
    True if spectrum has already been calculated, otherwise false
  ebins_checksum : string
    md5sum of the ebins the spectrum was last calculated on. Used to
    identify if new calculations are required or can just return the
    previous value.
  docont : bool
    Calculate the continuum emission
  dopseudo : bool
    Calculate the pseudocontinuum emission
  """
  def __init__(self, cocoentry, docont=True, dopseudo=True):
    """
    Initialization

    Parameters
    ----------
    cocoentry : numrec
      The data for 1 element/ion
    parentElementSpectrum : ElementSpectrum
      Parent ElementSpectrum object
    """

    if type(cocoentry)==bool:
      if cocoentry == False:
        cocoentry={}
        cocoentry['N_Cont']=2
        cocoentry['N_Pseudo']=2
        cocoentry['E_Cont'] = numpy.array([0.0001,10000])
        cocoentry['Continuum'] = numpy.array([0.0,0.0])
        cocoentry['E_Pseudo'] = numpy.array([0.0001,10000])
        cocoentry['Pseudo'] = numpy.array([0.0,0.0])


    nEC = cocoentry['N_Cont']
    nEP = cocoentry['N_Pseudo']


    self.ECont = cocoentry['E_Cont'][:nEC]
    self.EPseudo = cocoentry['E_Pseudo'][:nEP]

    self.Cont = cocoentry['Continuum'][:nEC]
    self.Pseudo = cocoentry['Pseudo'][:nEP]

    self.spectrum_calculated = False
    self.ebins_checksum = False


    self.docont = docont
    self.dopseudo=dopseudo


  def return_spec(self, eedges, ebins_checksum = False):
    import scipy.integrate


    # get the checksum for the ebins, if not provided
    if ebins_checksum == False:
        # generate the checksum
      ebins_checksum = hashlib.md5(eedges).hexdigest()

    # see if the current checksum matches the stored one: if it does
    # do nothing, return existing data


    if ((self.ebins_checksum==ebins_checksum) &\
        (self.spectrum_calculated == True)):
      pass

    else:
      if self.docont:
        cont = expand_E_grid(eedges, len(self.ECont), self.ECont, self.Cont)
      else:
        cont = 0.0

      if self.dopseudo:
        pseudo = expand_E_grid(eedges, len(self.EPseudo), self.EPseudo, self.Pseudo)
      else:
        pseudo = 0.0

      self.ebins_checksum = ebins_checksum
      self.spectrum = cont+pseudo
      self.spectrum_calculated = True

    return self.spectrum




def apply_response(spectrum, rmf, arf=False):
  """
  Apply a response to a spectrum

  Parameters
  ----------
  spectrum : array(float)
    The spectrum, in counts/bin/second, to have the response applied to. Must be
    binned on the same grid as the rmf.
  rmf : string or pyfits.hdu.hdulist.HDUList
    The filename of the rmf or the opened rmf file
  arf : string or pyfits.hdu.hdulist.HDUList
    The filename of the arf or the opened arf file
  Returns
  -------
  array(float)
    energy grid (keV) for returned spectrum
  array(float)
    spectrum folded through the response
  """
#
# Update 2016-05-25
#
# Changed to return the energy grid and the spectrum, as apparently in some
# instruments these are not the same as the input energy grid.
  if arf:
    if type(arf)==str:
      arfdat = pyfits.open(arf)
    elif type(arf) == pyfits.hdu.hdulist.HDUList:
      arfdat = arf
    else:
      print("ERROR: unknown arf type, %s"%(repr(type(arf))))
      return
    res = spectrum * arfdat['SPECRESP'].data['SPECRESP']
  else:
    res = spectrum*1.0

  sret = numpy.matmul(res,matrix)
  return ebins, ret, sret





class NEISession(CIESession):
  """
  Load and generate a collisional ionization equilibrium spectrum

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : arraylike(int), optional
    The atomic number of elements to include (default all)
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  datacache : dict
    Any Atomdb FITS files which have to be opened are stored here
  spectra : CIESpectra
    Object storing the actual spectral data
  elements : list(int)
    Nuclear charge of elements to include.
  default_abundset : string
    The abundance set used for the original emissivity file calculation
  abundset : string
    The abundance set to be used for the returned spectrum
  abundsetvector : array_like(float)
    The relative abundance between default_abundset and abundset for each element
  response_set : bool
    Have we loaded a response (or set a dummy response)
  dolines : bool
    Calculate line emission
  docont : bool
    Calculate continuum emission
  dopseudo : bool
    Calculate pseudocontinuum emission
  broaden_limit : float
    Apply broadening to lines with epsilon > this value (ph cm3 s-1)
  thermal_broadening : bool
    Apply thermal broadening to lines (default = False)
  velocity_broadening : float
    Apply velocity broadening with this velocity (km/s). If <=0, do not apply.

  Examples
  --------

  Create a session instance:

  >>> s=CIESession()

  Set up the responses, in this case a dummy response from 0.1 to 10 keV

  >>> ebins = numpy.linspace(0.1,10,1000)
  >>> s.set_response(ebins, raw=True)

  Turn on thermal broadening

  >>> s.set_broadening(True)
  Will thermally broaden lines with emissivity > 1.000000e-18 ph cm3 s-1

  Return spectrum at 1.0keV

  >>> spec = s.return_spectrum(1.0)

  spec is in photons cm^3 s^-1 bin^-1; ebins are the bin edges (so spec is
  1 element shorter than ebins)
  """

  def __init__(self, linefile="$ATOMDB/apec_nei_line.fits",\
                     cocofile="$ATOMDB/apec_nei_comp.fits",\
                     elements=False, abundset='AG89'):
    """
    Initialization routine. Can set the line and continuum files here

    Input
    -----
    linefile : str or HDUList
      The filename of the line emissivity data, or the opened file.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the opened file.
    elements : array_like(int)
      The atomic numbers of the elements to include. Defaults to all (1-30)
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance
      for list of options.
    """


    self.datacache={}

    # Open up the APEC files
    self.set_apec_files(linefile, cocofile)


    # if elements are specified, use them. Otherwise, use Z=1-30
    if util.keyword_check(elements):
      self.elements = elements
    else:
      self.elements=list(range(1,const.MAXZ_NEI+1))

    # a hold for the spectra
    self.spectra=NEISpectrum(self.linedata, self.cocodata)


    # Set both the current and the default abundances to those that
    # the apec data was calculated on
    self.abundset=self.linedata[0].header['SABUND_SOURCE']
    self.default_abundset=self.linedata[0].header['SABUND_SOURCE']

    self.abundsetvector = numpy.zeros(const.MAXZ_NEI+1)
    for Z in self.elements:
      self.abundsetvector[Z] = 1.0

    #  but if another vector was already specified, use this instead
    if util.keyword_check(abundset):
      self.set_abundset(abundset)

    self.abund = numpy.zeros(const.MAXZ_NEI+1)

    for Z in self.elements:
      self.abund[Z]=1.0

    # Set a range of parameters which can be overwritten later
    self.response_set = False # have we loaded a response file?
    self.dolines=True # Include lines in spectrum
    self.docont=True # Include continuum in spectrum
    self.dopseudo=True # Include pseudo continuum in spectrum
    self.set_broadening(False, broaden_limit=1e-18)


  def set_apec_files(self, linefile="$ATOMDB/apec_nei_line.fits",\
                     cocofile="$ATOMDB/apec_nei_comp.fits"):
    """
    Set the apec line and coco files, and load up their data

    Parameters
    ----------
    linefile : str or HDUList
      The filename of the line emissivity data, or the opened file.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the opened file.

    Returns
    -------
    None

    Notes
    -----
    Updates self.linefile, self.linedata, self.cocofile and self.cocodata
    """
    if util.keyword_check(linefile):
      if isinstance(linefile, str):
        lfile = os.path.expandvars(linefile)
        if not os.path.isfile(lfile):
          print("*** ERROR: no such file %s. Exiting ***" %(lfile))
          return -1
        self.linedata = pyfits.open(lfile)
        self.linefile = lfile

      elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.linedata=linefile
        self.linefile=linefile.filename()

      else:
        print("Unknown data type for linefile. Please pass a string or an HDUList")

    if util.keyword_check(cocofile):

      if isinstance(cocofile, str):

        cfile = os.path.expandvars(cocofile)
        if not os.path.isfile(cfile):
          print("*** ERROR: no such file %s. Exiting ***" %(cfile))
          return -1
        self.cocodata=pyfits.open(cfile)
        self.cocofile=cfile

      elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.cocodata=cocofile
        self.cocofile=cocofile.filename()

      else:
        print("Unknown data type for cocofile. Please pass a string or an HDUList")


  def return_linelist(self, Te, tau, specrange, specunit='A', \
                               teunit='keV', apply_aeff=False):
    """
    Get the list of line emissivities vs wavelengths


    Parameters
    ----------
    Te : float
      Temperature in keV or K
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the lines in the linelist to
      modify their intensities.

    Returns
    -------
    linelist : array(dtype)
      The list of lines with lambda (A), energy (keV), epsilon (ph cm3 s-1),\
      epsilon_aeff (ph cm5 s-1) ion (string) and upper & lower levels.

    """
    kT = convert_temp(Te, teunit, 'keV')

    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]



    s= self.spectra.return_linelist(kT, tau, specrange=specrange, teunit='keV',\
                                        specunit=specunit, elements=self.elements,\
                                        abundances = ab)

    # do the response thing
    #resp  = s.response()

    if apply_aeff == True:
      ibin = numpy.zeros(len(s), dtype=int)
      for i, ss in enumerate(s):
        e = const.HC_IN_KEV_A/ss['Lambda']
        ibin[i] = numpy.where(self.specbins<e)[0][-1]

      s["Epsilon_Err"] = s['Epsilon']*self.aeff[ibin]

    return(s)


  def return_line_emissivity(self, Telist, taulist, Z, z1, up, lo, \
                             specunit='A', teunit='keV', \
                             apply_aeff=False, apply_abund=True,\
                             log_interp = True):
    """
    Get line emissivity as function of Te, tau. Assumes ionization from neutral.


    Parameters
    ----------
    Telist : float or array(float)
      Temperature(s) in keV or K
    taulist : float
      ionization timescale(s), ne * t (cm^-3 s).
    Z : int
      nuclear charge of element
    z1 : int
      ion charge +1 of ion
    up : int
      upper level for transition
    lo : int
      lower level for transition
    specunit : {'Angstrom','keV'}
      Units for wavelength or energy (a returned value)
    teunit : {'keV' , 'K'}
      Units of Telist (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the line emissivity in the
      linelist to modify their intensities.
    apply_abund : bool
      If true, apply the abundance set in the session to the result.
    log_interp : bool
      Interpolate between temperature on a log-log scale (default).
      Otherwise linear

    Returns
    -------
    ret : dict
      Dictionary containing:
      Te, tau, teunit: as input
      wavelength : line wavelength (A)
      energy : line energy (keV)
      epsilon : emissivity in ph cm^3 s-1 (or ph cm^5 s^-1 if apply_aeff=True)
                first index is temperature, second is tau.

    """

    Tevec, Teisvec = util.make_vec(Telist)
    tauvec, tauisvec = util.make_vec(taulist)


    kTlist = convert_temp(Tevec, teunit, 'keV')
    if apply_abund:
      ab = self.abund[Z]*self.abundsetvector[Z]
    else:
      ab = 1.0

    eps = numpy.zeros([len(Tevec), len(tauvec)])
    ret={}

    for itau, tau in enumerate(tauvec):
      for ikT, kT in enumerate(kTlist):
        e, lam = self.spectra.return_line_emissivity(kT, tau, Z, z1, \
                                                     up, lo, \
                                                     specunit='A', \
                                                     teunit='keV', \
                                                     abundance=ab)
        eps[ikT, itau] = e
        if lam != False:
          ret['wavelength'] = lam * 1.0
        else:
          ret['wavelength'] = None

    ret['Te'] = Telist
    ret['tau'] = taulist
    ret['teunit'] = teunit
    if ret['wavelength'] != None:
      ret['energy'] = const.HC_IN_KEV_A/ret['wavelength']
    else:
      ret['energy'] = None


    if apply_aeff == True:
      e = ret['energy']
      ibin = numpy.where(self.specbins<e)[0][-1]

      eps = eps*self.aeff[ibin]


    # now correct for vectors

    if not tauisvec:
      eps=eps[:,0]
      if not Teisvec:
        eps = eps[0]
    else:
      if not Teisvec:
        eps = eps[0,:]

    ret['epsilon'] = eps

    return ret

  def return_spectrum(self,  Te, tau, init_pop=False, Te_init=False, teunit='keV', nearest=False,\
                      get_nearest_t=False, log_interp=True):
    """
    Get the spectrum at an exact temperature.
    Interpolates between 2 neighbouring spectra

    Finds HDU with kT closest to desired kT in given line or coco file.

    Opens the line or coco file, and looks for the header unit
    with temperature closest to te. Use result as index input to make_spectrum

    Parameters
    ----------
    Te : float
      Temperature in keV or K
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    raw : bool
      If set, return the spectrum without response applied. Default False.
    nearest : bool
      If set, return the spectrum from the nearest tabulated temperature
      in the file, without interpolation
    get_nearest_t : bool
      If set, and `nearest` set, return the nearest tabulated temperature
      as well as the spectrum.

    Returns
    -------
    spectrum : array(float)
      The spectrum in photons cm^5 s^-1 bin^-1, with the response, or
      photons cm^3 s^-1 bin^-1 if raw is set.
    nearest_T : float, optional
      If `nearest` is set, return the actual temperature this corresponds to.
      Units are same as `teunit`
    """

    # Check that there is a response set
    if not self.response_set:
      raise util.ReadyError("Response not yet set: use set_response to set.")

    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]


    self.spectra.ebins = self.specbins
    self.spectra.ebins_checksum=hashlib.md5(self.spectra.ebins).hexdigest()
    s= self.spectra.return_spectrum(Te, tau, init_pop=init_pop, \
                                    Te_init=Te_init, teunit=teunit, \
                                    nearest = nearest,elements = el_list, \
                                    abundances=ab, log_interp=True)

    ss = self.apply_response(s)

    return ss



class NEISpectrum(CIESpectrum):
  """
  A class holding the emissivity data for NEI emission, and returning
  spectra

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : arraylike(int), optional
    The atomic number of elements to include (default all)
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  session : CIESession
    The parent CIESession
  SessionType : string
    "CIE"
  spectra : dict of ElementSpectra
    a dictionary containing the emissivity data for each HDU,
    subdivided by element (spectra[12][18] is an ElementSpectrum object
    containing the argon data for the 12th HDU)
  kTlist : array
    The temperatures for each emissivity HDU, in keV
  logkTlist : array
    log of kTlist
  """

  def __init__(self, linedata, cocodata):
    """
    Initializes the code. Populates the line and emissivity data in all
    temperature HDUs.

    Parameters
    ----------
    linedata :
      The parent CIESession
    """

    self.datacache={}
    self.SessionType = 'NEI'

    picklefname = os.path.expandvars('$ATOMDB/spectra_%s_%s.pkl'%\
                                (linedata[0].header['CHECKSUM'],\
                                 cocodata[0].header['CHECKSUM']))


    if os.path.isfile(picklefname):
      self.spectra = pickle.load(open(picklefname,'rb'))
      self.kTlist = self.spectra['kTlist']
    else:
      self.spectra={}
      self.kTlist = numpy.array(linedata[1].data['kT'].data)
      self.spectra['kTlist'] = numpy.array(linedata[1].data['kT'].data)
      for ihdu in range(len(self.kTlist)):
        self.spectra[ihdu]={}
        self.spectra[ihdu]['kT'] = self.kTlist[ihdu]
        ldat = numpy.array(linedata[ihdu+2].data.data)
        cdat = numpy.array(cocodata[ihdu+2].data.data)


        Zarr = numpy.zeros([len(ldat), const.MAXZ_NEI+1], dtype=bool)
        Zarr[numpy.arange(len(ldat), dtype=int), ldat['Element']]=True


        for Z in range(1,const.MAXZ_NEI+1):

          if not Z in self.spectra[ihdu].keys():
            self.spectra[ihdu][Z] = {}

          for z1 in range(1,Z+2):
            isz1 = (ldat['Ion_drv']==z1)
            isgood = isz1*Zarr[:,Z]
            ccdat = cdat[(cdat['Z']==Z) & (cdat['rmJ']==z1)]

            if len(ccdat)==0:
              ccdat = [False]
            self.spectra[ihdu][Z][z1]=ElementSpectrum(ldat[isgood],\
                                                  ccdat[0], Z, z1_drv=z1)


      pickle.dump(self.spectra, open(picklefname,'wb'))


    self.logkTlist=numpy.log(self.kTlist)




  def return_spectrum(self, Te, tau, init_pop=False, Te_init=False, teunit='keV', nearest = False,
                             elements=False, abundances=False, log_interp=True):

    """
    Return the spectrum of the element on the energy bins in
    self.session.specbins

    Parameters
    ----------
    Te : float
      Electron temperature (default, keV)
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    teunit : string
      Units of kT (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.

    Returns
    -------
    spec : array(float)
      The element's emissivity spectrum, in photons cm^3 s^-1 bin^-1
    """

    # get kT in keV
    kT = convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_NEI+1)


    if abundances == False:
      abundances = {}
      for Z in elements:
        abundances[Z] = 1.0

    s = 0.0

    for Z in elements:

      abund = abundances[Z]
      if abund > 0:
        if Te_init != False:
          kT_init= convert_temp(Te_init, teunit, 'keV')

          # calculate the ion fraction
        ionfrac = apec.solve_ionbal_eigen(Z, kT, init_pop=init_pop, tau=tau, Te_init=Te_init, \
                     teunit='keV', datacache=self.datacache)

        for z1 in range(1, Z+2):
          if ionfrac[z1-1]>1e-10:

              # calculate minimum emissivitiy to broaden, accounting for ion
              # and element abundance.
            epslimit =  self.broaden_limit/(abund*ionfrac[z1-1])



              # return a broadened spectrum
            ss=0.0
            for i in range(len(ikT)):
              ss=0.0
              sss = self.spectra[ikT[i]][Z][z1].return_spectrum(self.ebins,\
                                  kT,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening) *\
                                  abund*ionfrac[z1-1]*f[i]

              # This may slow everythign down and should be moved outside the loop?
              if log_interp:
                ss += numpy.log(sss+const.MINEPSOFFSET)
              else:
                ss += sss
            if log_interp:
              ss = numpy.exp(ss)-const.MINEPSOFFSET*len(ikT)

            s += ss
    return s


  def return_line_emissivity(self, Te, tau, Z, z1, up, lo, specunit='A',
                             teunit='keV', abundance=1.0,
                             log_interp = True):
    """
    Return the emissivity of a line at kT, tau. Assumes ionization from neutral for now


    Parameters
    ----------
    Te : float
      Temperature in keV or K
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    Z : int
      nuclear charge of element
    z1 : int
      ion charge +1 of ion
    up : int
      upper level for transition
    lo : int
      lower level for transition
    specunit : {'Angstrom','keV'}
      Units for wavelength or energy (a returned value)
    teunit : {'keV' , 'K'}
      Units of Telist (kev or K, default keV)
    abundance : float
      Abundance to multiply the emissivity by

    Returns
    -------
    Emissivity : float
      Emissivity in photons cm^3 s^-1
    spec : float
      Wavelength or Energy of line, depending on specunit
    """
    kT = convert_temp(Te, teunit, 'keV')
    print("LOG INTERP", log_interp)
    ikT, f = self.get_nearest_Tindex(kT, \
                                     teunit='keV', \
                                     nearest=False, \
                                     log_interp=log_interp)
    #ikT has the 2 nearest temperature indexes
    # f has the fraction for each

    init_pop = numpy.zeros(Z+1)
    init_pop[0] = 1.0

    ionfrac = apec.solve_ionbal_eigen(Z, \
                                      kT, \
                                      init_pop=init_pop, \
                                      tau=tau, \
                                      teunit='keV', \
                                      datacache=self.datacache)

    eps = 0.0
    lam = 0.0
    print('--------')

      # find lines which match
    for z1_drv in range(1,Z+2):
      # ions which don't exist get skipped
      if ionfrac[z1_drv-1] <= 1e-10: continue
      eps_in = numpy.zeros(len(ikT))

      for i in range(len(ikT)):
        iikT =ikT[i]

        llist = self.spectra[iikT][Z][z1_drv].return_linematch(Z,z1,up,lo)

        for line in llist:
          # add emissivity
          eps_in[i] += line['Epsilon']
          lam = line['Lambda']

      if log_interp:
        eps_out = 0.0
        for i in range(len(ikT)):
          eps_out += f[i]*numpy.log(eps_in[i]+const.MINEPSOFFSET)
        eps += numpy.exp(eps_out-const.MINEPSOFFSET)*abundance * ionfrac[z1_drv-1]
      else:
        eps_out = 0.0
        for i in range(len(ikT)):
          eps_out += f[i]*eps_in[i]
        eps += eps_out*abundance * ionfrac[z1_drv-1]

    print('epsilon',  eps)

    if specunit == 'keV':
      lam = const.HC_IN_KEV_A/lam
    return eps, lam

  def return_linelist(self,  Te, tau, init_pop=False, Te_init=False,
                      teunit='keV', nearest = False, specrange=False,
                      specunit='A', elements=False, abundances=False):

    """
    Return the linelist of the element

    Parameters
    ----------

    Te : float
      Electron temperature (default, keV)
    teunit : string
      Units of kT (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Ansgtrom','keV'}
      Units for specrange (default A)


    """

    # get kT in keV
    kT = convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest)

    if abundances == False:
      abundances = {}
      for Z in elements:
        abundances[Z] = 1.0

    linelist = False

    for Z in elements:
      abund = abundances[Z]
      if abund > 0:
        elemlinelist = False

        print("THIS ISN'T IMPLEMENTED YET ADD IN ION BY ION CALC NOW")

        for i in range(len(ikT)):


          ss = self.spectra[ikT[i]][Z].return_linelist(specrange,\
                                  teunit='keV', specunit=specunit)

          if len(ss) > 0:
            ss['Epsilon']*=abund*f[i]
            if elemlinelist==False:
              elemlinelist = ss
            else:
              isnew = numpy.zeros(len(ss), dtype=bool)

              for inew, new in enumerate(ss):
                imatch = numpy.where((new['Element']==elemlinelist['Element']) &\
                                     (new['Ion']==elemlinelist['Ion']) &\
                                     (new['UpperLev']==elemlinelist['UpperLev']) &\
                                     (new['LowerLev']==elemlinelist['LowerLev']))[0]
                if len(imatch)==1:
                  elemlinelist[imatch[0]]['Epsilon']+=new['Epsilon']
                else:
                  isnew[inew]=True

              s = sum(isnew)
              if s > 0:
                elemlinelist = numpy.append(elemlinelist, ss[isnew])
        if elemlinelist != False:
          if linelist==False:
            linelist = elemlinelist
          else:
            linelist =  numpy.append(linelist, elemlinelist)

    return linelist


class TestSpectrum():
  def __init__(self, linefile='/export1/atomdb_latest/apec_v3.0.9_line.fits', \
                     cocofile='/export1/atomdb_latest/apec_v3.0.9_coco.fits'):
    import pickle
    ldat = pyfits.open(linefile)
    cdat = pyfits.open(cocofile)
    self.kTlist = numpy.array(ldat[1].data.data)['kT']

    self.continuum = {}
    for ikT in range(len(self.kTlist)):
      ihdu=ikT+2
      self.continuum[ikT] = {}
      nd = numpy.array(cdat[ihdu].data.data)
      #print(nd)
      for il in range(len(nd)):
      #  print(l['Z'])
        Z = nd[il]['Z']

        self.continuum[ikT][Z] = TestContinuumData(nd[il])


    pickle.dump(self, open('testpickle1.pkl', 'wb'))

  def makeSpectrum(self, ihdu, ebins):
    i=0
    ikT=20

    s = numpy.zeros(len(ebins)-1)
    for Z in self.continuum[ikT].keys():
      s+=  self.continuum[ikT][Z].return_spec(ebins)
      i+=1


class TestContinuumData():
  """
  A class holding the continuum data for an element in one HDU

  Parameters
  ----------
  cocoentry : array(cocodatatype)
    A single row from the continuum data in an AtomDB file.

  Attributes
  ----------
  ECont : array(float)
    The continuum energies (keV)
  EPseudo : array(float)
    The pseudocontinuum energies (keV)
  Cont : array(float)
    The continuum emissivities (ph cm^3 s^-1 keV^-1)
  Pseudo : array(float)
    The pseudocontinuum energies (ph cm^3 s^-1 keV^-1)
  spectrum_calculated : bool
    True if spectrum has already been calculated, otherwise false
  ebins_checksum : string
    md5sum of the ebins the spectrum was last calculated on. Used to
    identify if new calculations are required or can just return the
    previous value.
  docont : bool
    Calculate the continuum emission
  dopseudo : bool
    Calculate the pseudocontinuum emission
  """
  def __init__(self, cocoentry, docont=True, dopseudo=True):
    """
    Initialization

    Parameters
    ----------
    cocoentry : numrec
      The data for 1 element/ion
    """

    if type(cocoentry)==bool:
      if cocoentry == False:
        cocoentry={}
        cocoentry['N_Cont']=2
        cocoentry['N_Pseudo']=2
        cocoentry['E_Cont'] = numpy.array([0.0001,10000])
        cocoentry['Continuum'] = numpy.array([0.0,0.0])
        cocoentry['E_Pseudo'] = numpy.array([0.0001,10000])
        cocoentry['Pseudo'] = numpy.array([0.0,0.0])


    nEC = cocoentry['N_Cont']
    nEP = cocoentry['N_Pseudo']


    self.ECont = cocoentry['E_Cont'][:nEC]
    self.EPseudo = cocoentry['E_Pseudo'][:nEP]

    self.Cont = cocoentry['Continuum'][:nEC]
    self.Pseudo = cocoentry['Pseudo'][:nEP]

    # store the parent ElementSpectrum

    self.spectrum_calculated = False
    self.ebins_checksum = False


    self.docont = docont
    self.dopseudo=dopseudo


  def return_spec(self, eedges, ebins_checksum = False):
    import scipy.integrate


    # get the checksum for the ebins, if not provided
    if ebins_checksum == False:
      # generate the checksum
      ebins_checksum = hashlib.md5(eedges).hexdigest()

    # see if the current checksum matches the stored one: if it does
    # do nothing, return existing data


    if ((self.ebins_checksum==ebins_checksum) &\
        (self.spectrum_calculated == True)):
      pass

    else:
      if self.docont:
        cont = expand_E_grid(eedges, len(self.ECont), self.ECont, self.Cont)
      else:
        cont = 0.0

      if self.dopseudo:
        pseudo = expand_E_grid(eedges, len(self.EPseudo), self.EPseudo, self.Pseudo)
      else:
        pseudo = 0.0

      self.ebins_checksum = ebins_checksum
      self.spectrum = cont+pseudo
      self.spectrum_calculated = True

    return self.spectrum

