
"""
The modules is separated by the type of spectrum you wish to model. For
now, the CIE, NEI and PShock models are implemented.

Roughly speaking, you load a [Modelname]Session, and you can then
obtain:

   linelist :
     list of lines in a certain wavelength interval)

   line_emissivity :
     the emissivity of a specific line as a function of temperature,
     line intensity and more.

   spectrum :
     the emissivity in photons cm^3 bin^-1 s^-1. This requires a
     response be set first

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
from scipy.sparse import bsr_array
import time
import warnings

def __make_spectrum(bins, index, linefile="$ATOMDB/apec_line.fits",\
                  cocofile="$ATOMDB/apec_coco.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  elements=False, abund=False, dummyfirst=False,\
                  dolines = True, docont=True, dopseudo=True):

  """
  DEPRECATED
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
       Calculate line emission (default True)
  docont : bool
       Calculate Continuum emission (default True)
  dopseudo : bool
       Calculate PseudoContinuum (weak line) emission (default True)

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
  warnings.warn("make_spectrum is a deprecated function and will be removed. Use CIESession.return_spectrum", DeprecationWarning, stacklevel=3)


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



def __make_ion_spectrum(bins, index, Z,z1, linefile="$ATOMDB/apec_nei_line.fits",\
                  cocofile="$ATOMDB/apec_nei_comp.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  abund=False, dummyfirst=False, nei = True,\
                  dolines = True, docont=True, dopseudo=True):

  """
  DEPRECATED
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
    Calculate line emission (default True)
  docont : bool
    Calculate Continuum emission (default True)
  dopseudo : bool
    Calculate PseudoContinuum (weak line) emission (default True)


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
  warnings.warn("make_ion_spectrum is a deprecated function and will be removed. Use NEISession.return_spectrum", DeprecationWarning, stacklevel=3)


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
def __add_lines(Z, abund, lldat, ebins, z1=False, z1_drv=False, \
              broadening=False, broadenunits='A'):
  """
  DEPRECATED

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


  warnings.warn("add_lines is a deprecated function and will be removed", DeprecationWarning, stacklevel=3)
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
        spectrum+=atomdb._addline2(ebins, const.HC_IN_KEV_A/ll['lambda'], \
                 ll['epsilon']* abund,\
                 broadening*const.HC_IN_KEV_A/(ll['lambda']**2))
    else:
      for ll in l:
        spectrum+=atomdb._addline2(ebins, const.HC_IN_KEV_A/ll['lambda'], \
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


def __get_index(te, filename='$ATOMDB/apec_line.fits', \
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
def __list_lines(specrange, lldat=False, index=False, linefile=False,\
              units='angstroms', Te=False, teunit='K', minepsilon=1e-20):
  """
  DEPRECATED
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

  warnings.warn("list_lines is a deprecated function and will be removed. Use CIESession.return_linelist", DeprecationWarning, stacklevel=3)


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
def __list_nei_lines(specrange, Te, tau, Te_init=False,  lldat=False, linefile=False,\
              units='angstroms', teunit='K', minepsilon=1e-20, \
              datacache=False):
  """
  DEPRECATED
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
  warnings.warn("list_nei_lines is a deprecated function and will be removed. Use NEISession.return_linelist", DeprecationWarning, stacklevel=3)

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
    ionbal[Z] = apec.return_ionbal(Z, kT, tau=tau, Te_init = kT_init,\
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def __print_lines(llist, specunits = 'A', do_cfg=False):
  """
  DEPRECATED
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

  warnings.warn("print_lines is a deprecated function and will be removed.", DeprecationWarning, stacklevel=3)

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

def __make_ion_index_continuum(bins,  element, \
                             index = False,\
                             cocofile='$ATOMDB/apec_coco.fits',\
                             binunits = 'keV', \
                             fluxunits='ph', no_coco=False, no_pseudo=False, \
                             ion=0,\
                             broadening=False,\
                             broadenunits='keV'):
  """
  DEPRECATED
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

  warnings.warn("make_ion_index_continuum is a deprecated function and will be removed. Use NEISession.return_spectrum", DeprecationWarning, stacklevel=3)

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
    spectrum += _expand_E_grid(bins,\
                                d['N_cont'],\
                                d['E_cont'],\
                                d['Continuum'])

  if not(no_pseudo):
    spectrum += _expand_E_grid(bins,\
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
def _expand_E_grid(eedges, n,Econt_in_full, cont_in_full):

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
def __broaden_continuum(bins, spectrum, binunits = 'keV', \
                      broadening=False,\
                      broadenunits='keV'):
  """
  DEPRECATED

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
  warnings.warn("broaden_continuum is a deprecated function and will be removed", DeprecationWarning, stacklevel=3)

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

def _get_response_ebins(rmf):
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


def __get_effective_area(rmf, arf=False):
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

class _Gaussian_CDF():
  """
  For fast interpolation, pre-calculate the CDF and interpolate it when
  broadening lines

  Parameters
  ----------
  None

  Examples
  --------

  Create a CDF instance:

  >>> s=_Gaussian_CDF()

  Broaden a line on ebins grid, with centroid and width.

  >>> cdf = a.broaden(centroid, width, ebins)

  Convert to flux in each bin

  >>> flux= cdf[1:]-cdf[:-1]


  """
  def __init__(self):
    from scipy.stats import norm
    self.x = numpy.linspace(-6,6,2400)
    self.cdf = norm.cdf(self.x)
    self.broadentype='Gaussian'

  def broaden(self, centroid, width, ebins):
    """
    Broaden a line, return CDF

    Parameters
    ----------
    centroid : float
      The line energy (keV)
    width : float
      The sigma of the normal distribution, in keV
    ebins : array(float)
      Energy grid to return CDF on

    Returns
    -------
    CDF : array(float)
      cumulative flux distribution of linen at each bin edge.
    """

    # move the energy grid
    etmp = (ebins-centroid)/width

    # interpolate to get the appropriate CDF values
    ret=numpy.interp(etmp, self.x, self.cdf)

    return ret


class _Lorentzian_CDF():
  """
  For fast interpolation, pre-calculate the CDF and interpolate it when
  broadening lines

  Parameters
  ----------
  None

  Examples
  --------

  Create a CDF instance:

  >>> s=_Gaussian_CDF()

  Broaden a line on ebins grid, with centroid and width.

  >>> cdf = a.broaden(centroid, width, ebins)

  Convert to flux in each bin

  >>> flux= cdf[1:]-cdf[:-1]


  """



  def __init__(self):
    from scipy.stats import cauchy
    self.x = numpy.linspace(-12,12,4800)
    self.cdf = cauchy.cdf(self.x)
    self.broadentype='Lorentzian'

  def broaden(self, centroid, width, ebins):
    """
    Broaden a line, return CDF

    Parameters
    ----------
    centroid : float
      The line energy (keV)
    width : float
      The sigma of the normal distribution, in keV
    ebins : array(float)
      Energy grid to return CDF on

    Returns
    -------
    CDF : array(float)
      cumulative flux distribution of linen at each bin edge.
    """

    # move the energy grid
    etmp = (ebins-centroid)/width

    # interpolate to get the appropriate CDF values
    ret=numpy.interp(etmp, self.x, self.cdf)

    return ret



class _Voigt_CDF():
  """
  For fast interpolation, pre-calculate the CDF and interpolate it when
  broadening lines

  Parameters
  ----------
  None

  Examples
  --------

  Create a CDF instance:

  >>> s=_Gaussian_CDF()

  Broaden a line on ebins grid, with centroid and width.

  >>> cdf = a.broaden(centroid, width, ebins)

  Convert to flux in each bin

  >>> flux= cdf[1:]-cdf[:-1]


  """



  def __init__(self, sigma, gamma):
    from scipy.special import voigt_profile
    self.x = numpy.linspace(-12,12,2400)
    self.broadentype='Voigt'



  def __recalc(self, sigma, gamma):

    if ((self.sigma==sigma) & \
        (self.gamma == gamma)): pass



    edges = (self.x[1:]+self.x[:-1])/2
    edges = numpy.append( edges[0]-(edges[1]-edges[0]), \
                          edges, \
                          edges[-1] + (edges[-1]-edges[-2]))


    dx = edges[1:]-edges[:-1]

    pdfmid =  voigt_profile(self.x, sigma, gamma)*dx

    self.cdf = numpy.cumsum(pdfmid)


  def broaden(self, centroid, ebins, sigma, gamma, test_recalc=False):
    """
    Broaden a line, return CDF

    Parameters
    ----------
    centroid : float
      The line energy (keV)
    width : float
      The sigma of the normal distribution, in keV
    ebins : array(float)
      Energy grid to return CDF on

    Returns
    -------
    CDF : array(float)
      cumulative flux distribution of linen at each bin edge.
    """

    # move the energy grid
    etmp = (ebins-centroid)

    # If required, recalculate the parameters
    if test_recalc:
      self.__recalc(sigma, gamma)

    # interpolate to get the appropriate CDF values
    ret=numpy.interp(etmp, self.x, self.cdf)

    return ret


class CIESession():
  """
  Load and generate a collisional ionization equilibrium spectrum

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : iterable of int
       Elements to include, listed by atomic number. if not set, include all.
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  datacache : dict
    Any Atomdb FITS files which have to be opened are stored here
  spectra : CIESpectra
    Object storing the actual spectral data
  elements : iterable of int
    Elements to include, listed by atomic number. if not set, include all.
  default_abundset : string
    The abundance set used for the original emissivity file calculation
  abundset : string
    The abundance set to be used for the returned spectrum
  abundsetvector : array_like(float)
    The relative abundance between default_abundset and abundset for each element
  response_set : bool
    Have we loaded a response (or set a dummy response)
  dolines : bool
    Calculate line emission (default True)
  docont : bool
    Calculate Continuum emission (default True)
  dopseudo : bool
    Calculate PseudoContinuum (weak line) emission (default True)
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

  Set up the responses, in this case a dummy response from 0.1 to 10 keV,
  with area 1cm^2 in each bin

  >>> ebins = numpy.linspace(0.1,10,1000)
  >>> s.set_response(ebins, raw=True)

  (Alternatively, for a real response file, s.set_response(rmffile, arf=arffile)

  Turn on thermal broadening

  >>> s.set_broadening(True)
  Will thermally broaden lines with emissivity > 1.000000e-18 ph cm3 s-1

  Return spectrum at 1.0keV

  >>> spec = s.return_spectrum(1.0)

  spec is in photons cm^5 s^-1 bin^-1; ebins are the bin edges (so spec is
  1 element shorter than ebins)
  """

  def __init__(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits",\
                     elements=False,\
                     abundset='AG89'):
    """
    Initialization routine. Can set the line and continuum files here

    Parameters
    -----
    linefile : str or HDUList
      The filename of the line emissivity data, or the file opened with pyfits.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the file opened with pyfits.
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance
      for list of options.

    Returns
    -------
    None
    """
    self.SessionType='CIE'
    self._session_initialise1(linefile, cocofile, elements, abundset)

    # a hold for the spectra
    self.spectra=_CIESpectrum(self.linedata, self.cocodata)

    self._session_initialise2()


  def _session_initialise1(self, linefile, cocofile, elements, abundset):
    """
    This routine does all the initialization which is the same for all
    the different session classes, separates them out from the
    instance specific ones

    Parameters
    ----------
    linefile : str or HDUList
      The filename of the line emissivity data, or the file opened with pyfits.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the file opened with pyfits.
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance
      for list of options.

    Returns
    -------
    None
    """
    self.datacache={}

    # Open up the APEC files
    self._set_apec_files(linefile, cocofile)

    # if elements are specified, use them. Otherwise, use Z=1-30
    if util.keyword_check(elements):
      self.elements = elements
    else:
      if self.SessionType=='CIE':
        self.elements=list(range(1,const.MAXZ_CIE+1))
      elif self.SessionType in ['NEI','PShock','Kappa']:
        self.elements=list(range(1,const.MAXZ_NEI+1))


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
    self.do_eebrems = True # Do electron-electron bremsstrahlung

    # no response set yet
    self.rmffile=False
    self.arffile=False
    self.raw_response=False

    # verbosity
    self.verbose = False


  def _session_initialise2(self):
    """
    This routine does the remaining initialization which is the same for all
    the different session classes, separates them out from the
    instance specific ones, after the Spectrum class has been generated

    Parameters
    ----------
    None

    Return
    ------
    None
    """

    self.set_broadening(False, broaden_limit=1e-18)
    self.cdf = _Gaussian_CDF()


  def set_broadening(self, thermal_broadening, broaden_limit=False, \
                           velocity_broadening=0.0, \
                           velocity_broadening_units='km/s',\
                           thermal_broaden_temperature=None,\
                           teunit='keV'):

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
    thermal_broaden_temperature : float
      If not None, use this temperature for all line broadening instead of
      plasma temperature (keV)
    Notes
    -----
    Updates attributes thermal_broadening, broaden_limit, velocity_broadening,
    velocity_broadening_units
    """

    self.thermal_broadening = thermal_broadening
    if broaden_limit != False:
      self.broaden_limit = broaden_limit

    if self.verbose:
      if self.thermal_broadening==True:
        print("Will thermally broaden lines with emissivity > %e ph cm3 s-1"%(self.broaden_limit))
      else:
        print("Will not thermally broaden lines")

    self.velocity_broadening=velocity_broadening

    allowed_velocity_broadening_units= ['km/s']

    if not velocity_broadening_units.lower() in ['km/s']:
      raise util.UnitsError("Error: velocity broadening units of %s is not in allowed set %s."%\
             (velocity_broadening_units, repr(allowed_velocity_broadening_units)))
      return

    self.velocity_broadening_units=velocity_broadening_units

    self.spectra.thermal_broadening = self.thermal_broadening
    self.spectra.broaden_limit = self.broaden_limit

    self.spectra.velocity_broadening=self.velocity_broadening
    self.spectra.velocity_broadening_units=self.velocity_broadening_units
    if thermal_broaden_temperature is not None:
      T = util.convert_temp(thermal_broaden_temperature, teunit, 'K')
      self.thermal_broaden_temperature=T
      self.spectra.thermal_broaden_temperature=T
    else:
      self.thermal_broaden_temperature=None
      self.spectra.thermal_broaden_temperature=None




  def set_response(self, rmf, arf=False, raw=False, sparse=False):
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
      A response has been loaded
    self.response_type : string
      'raw', 'standard' or 'sparse' depending on the type implemented
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
    sparse : bool
      If true, the rmf is stored as a sparse matrix and solved using sparse
      matrix algebra. Useful for large RMFs (e.g. XRISM).
      Tests show accuracy to 1 part in 10^{15}.
      Ignored if raw==True
      
    Returns
    -------
    none

    """

    if raw==True:
      if sparse:
        warnings.warn("Sparse matrix requested with raw response. Ignoring.")
      # make a diagonal perfect response
      self.specbins = rmf
      self.ebins_out = rmf
      self.specbin_units='keV'

      self.aeff = numpy.ones(len(rmf)-1)

      self.response_set = True
      self.specbins_set = True
      self.arf = False
      self.ebins_checksum =hashlib.md5(self.specbins).hexdigest()
      self.raw_response=True
      self.response_type = 'raw'
    else:

      if util.keyword_check(arf):
        if type(arf)==str:
          self.arffile = arf
          self.arfdat = pyfits.open(arf)
        elif type(arf) == pyfits.hdu.hdulist.HDUList:
          self.arfdat = arf
          self.arffile = arf.filename()
        else:
          print("ERROR: unknown arf type, %s"%(repr(type(arf))))
          return

        self.arf = numpy.array(self.arfdat['SPECRESP'].data['SPECRESP'])

      else:
        self.arf=1.0


      self.raw_response=False
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
      
      if not sparse:

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

        self.specbins, self.ebins_out = _get_response_ebins(self.rmf)

        if self.ebins_out[-1] < self.ebins_out[0]:
          # need to reverse things
          self.ebins_out=self.ebins_out[::-1]
          self.rmfmatrix = self.rmfmatrix[:,::-1]

        self.specbin_units='keV'
        self.aeff = self.rmfmatrix.sum(1)
        if util.keyword_check(self.arf):
          self.aeff *=self.arf
        self.response_set = True
        self.specbins_set = True
        self.response_type = 'standard'

        self.ebins_checksum =hashlib.md5(self.specbins).hexdigest()

      else:
        # sparse matrix time!
        data=[]
        row=[]
        col=[]
        
        chanoffset = self.rmf['EBOUNDS'].data['CHANNEL'][0]
        
        for ibin, i in enumerate(self.rmf[matrixname].data):
          lobound = 0
          for ngrp in range(i['N_GRP']):
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
              data.extend(i['MATRIX'][lobound:lobound+nchan[j]])
              row.extend(range(ilo,ihi))
              col.extend([ibin]*nchan[j])
              lobound = lobound+nchan[j]
  
        self.specbins, self.ebins_out = _get_response_ebins(self.rmf)
  
        data = numpy.array(data)
        row = numpy.array(row)
        col = numpy.array(col)
        
        if self.ebins_out[-1] < self.ebins_out[0]:
          # need to reverse things
          self.ebins_out=self.ebins_out[::-1]
          col = len(self.ebins_out)-col-1
  
  
        self.rmfmatrix =  bsr_array((data, (row, col)), shape=(len(self.specbins)-1, len(self.ebins_out)-1),\
                          dtype=data.dtype)
  
        
  
        self.specbin_units='keV'
        self.aeff = self.rmfmatrix.sum(1)
        if util.keyword_check(self.arf):
          self.aeff *=self.arf
        self.response_set = True
        self.specbins_set = True
        self.response_type = 'sparse'
  
        self.ebins_checksum =hashlib.md5(self.specbins).hexdigest()
  
  

    # this is now a check for 0 minimums
    if self.specbins_set:
      if self.specbins[0] <=0:
        warnings.warn('Response minimum energy is 0 keV, setting to small finite value (%e keV)'%(self.specbins[1]*1e-6))
        self.specbins[0] = self.specbins[1]*1e-6


  def return_spectrum(self, te, teunit='keV', nearest=False,\
                      get_nearest_t=False, log_interp=True,\
                      dolines=True, docont=True, dopseudo=True):
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
    nearest : bool
      If set, return the spectrum from the nearest tabulated temperature
      in the file, without interpolation
    get_nearest_t : bool
      If set, and `nearest` set, return the nearest tabulated temperature
      as well as the spectrum.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid, instead of linear.
    dolines : bool
      Calculate line emission (default True)
    docont : bool
      Calculate Continuum emission (default True)
    dopseudo : bool
      Calculate PseudoContinuum (weak line) emission (default True)

    Returns
    -------
    spectrum : array(float)
      The spectrum in photons cm^5 s^-1 bin^-1, with the response, or
      photons cm^3 s^-1 bin^-1 if raw is set.
    nearest_T : float, optional
      If `get_nearest_t` is set, return the actual temperature this corresponds to.
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
    s= self.spectra.return_spectrum(te, teunit=teunit, nearest=nearest,\
                                    elements = el_list, abundance=ab, \
                                    broaden_object = self.cdf, \
                                    log_interp=log_interp, dolines=dolines,\
                                    dopseudo=dopseudo, docont=docont, \
                                    do_eebrems = self.do_eebrems)
    ss = self._apply_response(s)

    return ss




  # def _return_test_spectrum(self, spectrum_in):
    # """
    # Get the spectrum at an exact temperature.
    # Interpolates between 2 neighbouring spectra

    # Finds HDU with kT closest to desired kT in given line or coco file.

    # Opens the line or coco file, and looks for the header unit
    # with temperature closest to te. Use result as index input to make_spectrum

    # Parameters
    # ----------
    # te : float
      # Temperature in keV or K
    # teunit : {'keV' , 'K'}
      # Units of te (kev or K, default keV)
    # raw : bool
      # If set, return the spectrum without response applied. Default False.
    # nearest : bool
      # If set, return the spectrum from the nearest tabulated temperature
      # in the file, without interpolation
    # get_nearest_t : bool
      # If set, and `nearest` set, return the nearest tabulated temperature
      # as well as the spectrum.

    # Returns
    # -------
    # spectrum : array(float)
      # The spectrum in photons cm^5 s^-1 bin^-1, with the response, or
      # photons cm^3 s^-1 bin^-1 if raw is set.
    # nearest_T : float, optional
      # If `nearest` is set, return the actual temperature this corresponds to.
      # Units are same as `teunit`
    # """

    # # Check that there is a response set
    # if not self.response_set:
      # raise util.ReadyError("Response not yet set: use set_response to set.")

    # # make element and abundance lists
    # ss = self.apply_response(spectrum_in)

    # return ss



  def _apply_response(self, spectrum):

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
      spectrum folded through the response
    """

    # if the response is raw, no need to matrix multiply (diagonal response, effectively)
    if self.response_type=='raw':
      return spectrum

    elif self.response_type=='standard':
      arfdat = self.arf

      ret = spectrum*self.arf
      
      try:
        ret = numpy.matmul(ret,self.rmfmatrix)
      except ValueError:
        try:
          ret = numpy.matmul(ret,self.rmfmatrix.transpose())
        except ValueError:
          if ret == 0:
            ret = numpy.zeros(len(self.ebins_out)-1)
      return ret
    elif self.response_type=='sparse':
      arfdat = self.arf
      ret = spectrum*self.arf
      ret = (self.rmfmatrix*ret).sum(1)
      pickle
      return(ret)
    
    else:
      raise util.OptionError('Unknown response type %s'%(self.response_type))

  def _set_apec_files(self, linefile, cocofile):
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

  def set_eebrems(self, do_eebrems):
    """
    Set whether to do the electron-electron bremsstrahlung in the spectrum

    Parameters
    ----------

    do_eebrems : bool
      If True, include the electron-electron bremsstrahlung calculation
      If False, don't (this is default, and matches XSPEC behaviour)
    """
    a = bool(do_eebrems)
    self.do_eebrems = a


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
        print('ab', abundvec)
        print('el',elementvec)
      else:

        self.abund[elementvec] = abundvec
    elif (eisvec):
      # set all these eleemtns to the same abundance
      for el in elementvec:
        self.abund[el]=abund

    else:
      self.abund[elements]=abund



  def set_abundset(self, abundstring=None):
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
    if abundstring==None:
      # Get abundance data file
      abunddata = atomdb.get_data(False, False, 'abund',datacache=self.datacache)

      # find the possible strings
      print("Possible abundance sets:")
      for name in abunddata[1].data.field('Source'):
        print("  %s"%(name))
      return

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
        pass


    # update the current abundance string to represent your input
    self.abundset=abundstring

    #self.recalc()


  def return_line_emissivity(self, Te, Z, z1, up, lo, \
                             specunit='A', teunit='keV', \
                             apply_aeff=False, apply_abund=True,\
                             log_interp = True):
    """
    Get line emissivity as function of Te.


    Parameters
    ----------
    Te : float or array(float)
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
      Units of Te (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the line emissivity in the
      linelist to modify their intensities.
    apply_abund : bool
      If true, apply the abundance set in the session to the result.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
      Otherwise linear

    Returns
    -------
    ret : dict
      Dictionary containing:

        Te, tau, teunit : as input

        wavelength : line wavelength (A)

        energy : line energy (keV)

        epsilon : emissivity in ph cm^3 s-1 (or ph cm^5 s^-1 if apply_aeff=True). First index is temperature, second is tau.

    """

    Tevec, Teisvec = util.make_vec(Te)

    kTlist = util.convert_temp(Tevec, teunit, 'keV')
    if apply_abund:
      ab = self.abund[Z]*self.abundsetvector[Z]
    else:
      ab = 1.0

    eps = numpy.zeros(len(Tevec))
    ret={}
    ret['wavelength']=None
    for ikT, kT in enumerate(kTlist):
      e, lam = self.spectra.return_line_emissivity(kT, Z, z1,\
                                                   up, lo,\
                                                   specunit='A',\
                                                   teunit='keV',\
                                                   abundance=ab)

      eps[ikT] = e
      if lam != False:
        ret['wavelength'] = lam * 1.0
      #else:
#        ret['wavelength'] = None

    ret['Te'] = Te
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
                               teunit='keV', apply_aeff=False, nearest=False,\
                               apply_binwidth=False):
    """
    Get the list of line emissivities vs wavelengths


    Parameters
    ----------
    Te : float
      Temperature in keV or K
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the lines in the linelist to
      modify their intensities.
    nearest : bool
      Return spectrum at nearest tabulated temperature, without interpolation
    apply_binwidth : bool
      Divide the line emissivity by the width of the bin they occupy
      to give emissivity per angstrom or per keV.

    Returns
    -------
    linelist : array(dtype)
      The list of lines with lambda (A), energy (keV), epsilon (ph cm3 s-1),\
      epsilon_aeff (ph cm5 s-1) ion (string) and upper & lower levels.

    """
    kT = util.convert_temp(Te, teunit, 'keV')


    # make element and abundance lists
    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]


    s= self.spectra.return_linelist(kT, specrange=specrange, teunit='keV',\
                                        specunit=specunit, elements=el_list,\
                                        abundance = ab, nearest=nearest)

    # do the response thing
    #resp  = s.response()
    if apply_aeff == True:

      epsilon_aeff =  self._apply_linelist_aeff(s, specunit, apply_binwidth)

      s['Epsilon'] = epsilon_aeff
    return(s)

  def _apply_linelist_aeff(self, linelist, specunit, apply_binwidth):
    """
    Apply effective area to the linelist, return in 'Epsilon_Err'.

    Parameters
    ----------
    linelist : array
      List of lines
    specunit : str
      A or keV
    apply_binwidth : bool
      If true, return emissivity per angstrom or per keV. If false, total.

    Returns
    -------
    emiss_aeff : array(float)
      Emissivity * Aeff
    """


    if specunit.lower()=='kev':
      binwidth = self.ebins_out[1:]-self.ebins_out[:-1]
      factor = numpy.zeros(len(s), dtype=float)
      for i, ss in enumerate(linelist):
        e = const.HC_IN_KEV_A/ss['Lambda']
        if e>self.specbins[-1]:
          factor[i] = 0.0
        elif e<self.specbins[0]:
          factor[i] = 0.0
        else:
          ibin = numpy.where(self.specbins<e)[0][-1]
          factor[i]=self.aeff[ibin]
          if apply_binwidth:
            factor[i] /= binwidth[ibin]

      emiss_aeff = linelist['Epsilon']*factor


    elif specunit.lower()=='a':
      wvbins=12.398425/self.ebins_out[::-1]
      binwidth = wvbins[1:]-wvbins[:-1]
      factor = numpy.zeros(len(linelist), dtype=float)
      for i, ss in enumerate(linelist):
        e = ss['Lambda']
        if e>wvbins[-1]:
          factor[i] = 0.0
        elif e<wvbins[0]:
          factor[i] = 0.0
        else:
          ibin = numpy.where(wvbins<e)[0][-1]
          factor[i]=self.aeff[::-1][ibin]
          if apply_binwidth:
            factor[i] /= binwidth[ibin]
      emiss_aeff = linelist['Epsilon']*factor
    return emiss_aeff



  def _adjust_line(self, change, Z=0, z1=0, z1_drv=0, upper=0,lower=0, quantity="Epsilon", method="Replace", trackchanges=False):
    """
    Change the emissivity or wavelength of a line. Integer parameters set to 0 mean "all". Note this all
    happens in memory and does not edit the underlying files.

    Parameters
    ----------
    change : float or str
      If float, set the new value to this.
      If string
    Z : int
      Element
    z1 : int
      Ion
    z1_drv : int
      Driving ion
    upper : int
      Upper level
    lower : int
      Lower level
    quantity : string
      Change "Epsilon" or "Lambda" - emissivity or wavelength - by change
    method : string
      "Replace": replace existing value with change
      "Multiply" : multiply existing value with change
      "Divide" : divide existing value by change
      "Add" : add change to existing value
      "Subtract" : subtract change from existing
    Returns
    -------
    None
    """
    meth = method.lower()

    if Z==0:
      Zlist = self.elements
    else:
      Zlist = [Z]
    for Zt in Zlist:
      if z1==0:
        z1list = range(1,z1+2)
      else:
        z1list=[z1]
      for z1t in z1list:
        if z1_drv==0:
          if self.SessionType=='CIE':
            z1_drvlist=[0]
          else:
            z1_drvlist=range(1,z1+2)
        else:
          if self.SessionType=='CIE':
            z1_drvlist=[0]
          else:
            z1_drvlist=[z1_drv]

        for z1_drvt in z1_drvlist:

          # go through each temperature, see if there is line data for this
          # ion. If so, change it according to uppper, lower

          for ikT in range(len(self.spectra.kTlist)):
             # check if there is any data for this element

            try:
              ldat = self.spectra.spectra[ikT][Zt].lines.lines
            except KeyError:
              # element or temperature data doesn't exist
              continue

            tochange = numpy.ones(len(ldat), dtype=bool)

            if upper!=0:
              tochange[ldat['UpperLev'] != upper] = False
            if lower != 0:
              tochange[ldat['LowerLev'] != lower] = False

            if trackchanges:
              tochange[self.spectra.spectra[ikT][Zt].lines.changed==True] = False

            if sum(tochange) > 0:

              if meth=='replace':
                 ldat[quantity][tochange] = change
              if meth=='add':
                 ldat[quantity][tochange] += change
              if meth=='subtract':
                 ldat[quantity][tochange] -= change
              if meth=='divide':
                 ldat[quantity][tochange] /= change
              if meth=='multiply':
                 ldat[quantity][tochange] *= change

              self.spectra.spectra[ikT][Zt].lines.changed[tochange] = True


    return



  def _adjust_line_lambda(self, change, Z, z1, upper,lower, quantity="Epsilon", method="Replace", trackchanges=False):
    """
    Change the emissivity or wavelength of a line. Integer parameters set to 0 mean "all". Note this all
    happens in memory and does not edit the underlying files.

    Parameters
    ----------
    change : float or str
      If float, set the new value to this.
      If string
    Z : int
      Element
    z1 : int
      Ion
    upper : int
      Upper level
    lower : int
      Lower level
    quantity : string
      Change "Epsilon" or "Lambda" - emissivity or wavelength - by change
    method : string
      "Replace": replace existing value with change
      "Multiply" : multiply existing value with change
      "Divide" : divide existing value by change
      "Add" : add change to existing value
      "Subtract" : subtract change from existing
    Returns
    -------
    None
    """
    meth = method.lower()

    # see if this line is already listed
    try:
      value=self.spectra.fixwavelength[Z][z1][upper][lower]
    except AttributeError:
      self.spectra.fixwavelength={}

    if not Z in self.spectra.fixwavelength.keys():
      self.spectra.fixwavelength[Z]={}

    if not z1 in self.spectra.fixwavelength[Z].keys():
      self.spectra.fixwavelength[Z][z1]={}

    if not upper in self.spectra.fixwavelength[Z][z1].keys():
      self.spectra.fixwavelength[Z][z1][upper]={}

    if not lower in self.spectra.fixwavelength[Z][z1][upper].keys():
      self.spectra.fixwavelength[Z][z1][upper][lower]=quantity





    if Z==0:
      Zlist = self.elements
    else:
      Zlist = [Z]
    for Zt in Zlist:
      if z1==0:
        z1list = range(1,z1+2)
      else:
        z1list=[z1]
      for z1t in z1list:
        if z1_drv==0:
          if self.SessionType=='CIE':
            z1_drvlist=[0]
          else:
            z1_drvlist=range(1,z1+2)
        else:
          if self.SessionType=='CIE':
            z1_drvlist=[0]
          else:
            z1_drvlist=[z1_drv]

        for z1_drvt in z1_drvlist:

          # go through each temperature, see if there is line data for this
          # ion. If so, change it according to uppper, lower

          for ikT in range(len(self.spectra.kTlist)):
             # check if there is any data for this element

            try:
              ldat = self.spectra.spectra[ikT][Zt].lines.lines
            except KeyError:
              # element or temperature data doesn't exist
              continue

            tochange = numpy.ones(len(ldat), dtype=bool)

            if upper!=0:
              tochange[ldat['UpperLev'] != upper] = False
            if lower != 0:
              tochange[ldat['LowerLev'] != lower] = False

            if trackchanges:
              tochange[self.spectra.spectra[ikT][Zt].lines.changed==True] = False

            if sum(tochange) > 0:

              if meth=='replace':
                 ldat[quantity][tochange] = change
              if meth=='add':
                 ldat[quantity][tochange] += change
              if meth=='subtract':
                 ldat[quantity][tochange] -= change
              if meth=='divide':
                 ldat[quantity][tochange] /= change
              if meth=='multiply':
                 ldat[quantity][tochange] *= change

              self.spectra.spectra[ikT][Zt].lines.changed[tochange] = True


    return


class _CIESpectrum():
  """
  A class holding the emissivity data for CIE emission, and returning
  spectra

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)

  Attributes
  ----------
  SessionType : string
    "CIE"
  spectra : dict of _ElementSpectrum
    a dictionary containing the emissivity data for each HDU,
    subdivided by element (spectra[12][18] is an _ElementSpectrum object
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
    linefile : HDUList
      The line emissivity data
    cocofile : HDUList
      The continuum emissivity data
    """

    self.datacache = {}
    self.SessionType = 'CIE'

    picklefname = os.path.expandvars('$ATOMDB/spectra_%s_%s.pkl'%\
                                (linedata[0].header['CHECKSUM'],\
                                 cocodata[0].header['CHECKSUM']))
    havepicklefile = False
    if os.path.isfile(picklefname):
      havepicklefile = True

    if havepicklefile:
      try:
        self.spectra = pickle.load(open(picklefname,'rb'))
        self.kTlist = self.spectra['kTlist']
      except AttributeError:
        havepicklefile=False
        print("pre-stored data in %s is out of date. This can be caused by updates to the data "%(picklefname)+
              "or, more likely, changes to pyatomdb. Regenerating...")

        # delete the old file
        if os.path.isfile(picklefname):
          os.remove(picklefname)

    if not havepicklefile:
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

          self.spectra[ihdu][Z]=_ElementSpectrum(ldat[Zarr[:,Z]],\
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
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
      Otherwise linear

    Returns
    -------
    ikT : list[int]
      Index of temperature in HDU file (from 0, not 2)
    f : list[float]
      fractional weight to apply to each ikT. Should sum to 1.
    """
    kT = util.convert_temp(Te, teunit, 'keV')

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

  def return_spectrum(self, Te, teunit='keV', nearest = False,\
                             elements=False, abundance=False, log_interp=True,\
                             broaden_object=False,\
                             dolines=True, docont=True, dopseudo=True, \
                             do_eebrems = False):

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
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundance : dict(float)
      The abundances of each element, e.g. abund[6]=1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    broaden_object : class
      Object with routine "broaden" which applies line broadening. Usually a Gaussian.
    dolines : bool
      Calculate line emission (default True)
    docont : bool
      Calculate Continuum emission (default True)
    dopseudo : bool
      Calculate PseudoContinuum (weak line) emission (default True)
    do_eebrems : bool
      Calculate electron-electron bremsstrahlung emission (default False)

    Returns
    -------
    spec : array(float)
      The element's emissivity spectrum, in photons cm^3 s^-1 bin^-1
    """

    # get kT in keV
    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_CIE+1)


    if abundance == False:
      abundance = {}
      for Z in elements:
        abundance[Z] = 1.0

    s = 0.0

    # electron-electron bremsstrahlung electron counter
    if do_eebrems:
      nel=0.0
      rawabund = atomdb.get_abundance(datacache=self.datacache)

    # get the line broadening temperature. Typically this is kT unless
    # overridden using XXXSession.set_broadening
    if self.thermal_broaden_temperature is not None:
      Tb= util.convert_temp(self.thermal_broaden_temperature, 'K', 'keV')
    else:
      Tb=kT
    for Z in elements:
      abund = abundance[Z]
      if abund > 0:
        epslimit =  self.broaden_limit/abund

           # go caclulate the spectrum, with broadening as assigned.
        sss=0.0

        if len(ikT) == 1:
          ss = self.spectra[ikT[0]][Z].return_spectrum(self.ebins,\
                                  Tb,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening,\
                                  broaden_object=broaden_object,\
                                  dolines=dolines,\
                                  docont=docont,\
                                  dopseudo=dopseudo) *\
                                  abund

        else:
          ss1 = self.spectra[ikT[0]][Z].return_spectrum(self.ebins,\
                                  Tb,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening,\
                                  broaden_object=broaden_object,\
                                  dolines=dolines,\
                                  docont=docont,\
                                  dopseudo=dopseudo) *\
                                  abund

          ss2 = self.spectra[ikT[1]][Z].return_spectrum(self.ebins,\
                                  Tb,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening,\
                                  broaden_object=broaden_object,\
                                  dolines=dolines,\
                                  docont=docont,\
                                  dopseudo=dopseudo) *\
                                  abund

          ss = self._merge_spectra_temperatures(f,ss1,ss2,log_interp)

        s+=ss
        if do_eebrems:
          ionpop=apec.return_ionbal(Z, kT, datacache=self.datacache, teunit='keV')
          Zabundance = rawabund[Z]*abund

          tmp = sum(ionpop*numpy.arange(Z+1))*Zabundance
          nel +=tmp

    if do_eebrems:
      eespec = calc_ee_brems_spec(self.ebins, kT, nel)
      s+= eespec

    return s



  def return_line_emissivity(self, Te, Z, z1, up, lo, specunit='A',
                             teunit='keV', abundance=1.0,
                             log_interp = True):
    """
    Return the emissivity of a line at temperature Te.


    Parameters
    ----------
    Te : float
      Temperature in keV or K
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
      The abundances of the element, e.g. 1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.

    Returns
    -------
    eps : float
      Emissivity in photons cm^3 s^-1
    lam : float
      Wavelength (A) or Energy (keV) of line, depending on specunit
    """


    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, \
                                     teunit='keV', \
                                     nearest=False, \
                                     log_interp=log_interp)
    #ikT has the 2 nearest temperature indexes
    # f has the fraction for each

      # find lines which match
    eps_in = numpy.zeros(len(ikT))
    eps = 0.0
    lam = 0.0
    for i in range(len(ikT)):
      iikT =ikT[i]

      llist = self.spectra[iikT][Z].return_linematch(Z,z1,up,lo)
      for line in llist:
        # add emissivity
        eps_in[i] += line['Epsilon']
        lam = line['Lambda']


    if log_interp:
      eps_out = 0.0
      for i in range(len(ikT)):
        eps_out += f[i]*numpy.log(eps_in[i]+const.MINEPSOFFSET)
      eps += numpy.exp(eps_out-const.MINEPSOFFSET)*abundance
    else:
      eps_out = 0.0
      for i in range(len(ikT)):
        eps_out += f[i]*eps_in[i]
      eps += eps_out*abundance


    if specunit == 'keV':
      lam = const.HC_IN_KEV_A/lam
    return eps, lam


  def _merge_spectra_temperatures(self, f, spec1, spec2, log_interp):
    """
    Merge spectra 1 and 2, using fractions f, using log or non-log scaling
    as determined by loginterp

    Paramters
    ---------
    f : [float, float]
      The weight to apply to spec1 and spec 2. Should sum to 1.
    spec1 : array(float)
      The 1st spectrum to merge
    spec2 : array(float)
      The 2nd spectrum to merge
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.

    Returns
    -------
    spec : array(float)
      The merged spectrum
    """

    if log_interp:
      out = numpy.log(spec1+const.MINEPSOFFSET) * f[0]+\
            numpy.log(spec2+const.MINEPSOFFSET) * f[1]
      spec = numpy.exp(out)-const.MINEPSOFFSET

    else:
      spec = spec1*f[0] + spec2*f[1]

    spec[spec<0]=0.0

    return spec

  def _merge_linelist_duplicates(self, llist, by_ion_drv=False):
    """
    Look through the linelists and search for duplicates. Duplicated lines
    are summed.

    Parameters
    ----------
    llist : array(linelist)
      The list of lines to check
    by_ion_drv : bool
      If true, keep lines which are the same but have different ion_drv separate.
      Otherwise, merge them.

    Returns
    -------
    llist_out : array(linelist)
      The lines after filtering for duplicates.
    """
    # sort the data
    if by_ion_drv:
      llist =numpy.sort(llist, order=['Ion','Ion_drv','UpperLev','LowerLev'])
    else:
      llist =numpy.sort(llist, order=['Ion','UpperLev','LowerLev'])
    # merge
    keep = numpy.ones(len(llist), dtype=bool)

    # find each level where the next is for the same transition

    if by_ion_drv:
      j = numpy.where((llist[1:]['Ion']==llist[:-1]['Ion']) &\
                      (llist[1:]['Ion_drv']==llist[:-1]['Ion_drv']) &\
                      (llist[1:]['UpperLev']==llist[:-1]['UpperLev']) &\
                      (llist[1:]['LowerLev']==llist[:-1]['LowerLev']))[0]
    else:
      j = numpy.where((llist[1:]['Ion']==llist[:-1]['Ion']) &\
                    (llist[1:]['UpperLev']==llist[:-1]['UpperLev']) &\
                    (llist[1:]['LowerLev']==llist[:-1]['LowerLev']))[0]

    for jj in j:
      # add the emissivity to the second of the 2
      llist['Epsilon'][jj+1] += llist['Epsilon'][jj]
      keep[jj]=False

    # remove all the non-keepers
    llist_out = llist[keep]
    # fix the emissivities if log scaled

    return llist_out

  def _merge_linelists_temperatures(self, f, llist1, llist2, log_interp, by_ion_drv=False):
    """
    Merge linelists 1 and 2, using fractions f, using log or non-log scaling
    as determined by loginterp. Note that only lines which appear in both
    temperatures will be returned.

    Parameters
    ---------
    f : [float, float]
      The weight to apply to spec1 and spec 2. Should sum to 1.
    llist1 : array(float)
      The 1st linelist to merge
    llist2 : array(float)
      The 2nd linelist to merge
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    by_ion_drv : bool
      If true, keep lines which are the same but have different ion_drv separate.
      Otherwise, merge them.

    Returns
    -------
    llist : array(float)
      The merged linelist
    """


    # merge any repeated lines in together.

    llist1 = self._merge_linelist_duplicates(llist1, by_ion_drv=by_ion_drv)
    llist2 = self._merge_linelist_duplicates(llist2, by_ion_drv=by_ion_drv)


    # weight the emissivities and combine into one list

    llist = numpy.append(llist1, llist2)

    if log_interp:
      llist[:len(llist1)]['Epsilon'] = numpy.log(llist[:len(llist1)]['Epsilon'])*f[0]
      llist[len(llist1):]['Epsilon'] = numpy.log(llist[len(llist1):]['Epsilon'])*f[1]
    else:
      llist[:len(llist1)]['Epsilon'] = llist[:len(llist1)]['Epsilon']*f[0]
      llist[len(llist1):]['Epsilon'] = llist[len(llist1):]['Epsilon']*f[1]

    # sort the data

    if by_ion_drv:
      llist =numpy.sort(llist, order=['Ion','Ion_drv', 'UpperLev','LowerLev'])
    else:
      llist =numpy.sort(llist, order=['Ion','UpperLev','LowerLev'])

    # for denoting which levels to keep
    keep = numpy.zeros(len(llist), dtype=bool)

    # find each level where the next is for the same transition


    if by_ion_drv:
      j = numpy.where((llist[1:]['Ion']==llist[:-1]['Ion']) &\
                      (llist[1:]['Ion_drv']==llist[:-1]['Ion_drv']) &\
                      (llist[1:]['UpperLev']==llist[:-1]['UpperLev']) &\
                      (llist[1:]['LowerLev']==llist[:-1]['LowerLev']))[0]
    else:
      j = numpy.where((llist[1:]['Ion']==llist[:-1]['Ion']) &\
                      (llist[1:]['UpperLev']==llist[:-1]['UpperLev']) &\
                      (llist[1:]['LowerLev']==llist[:-1]['LowerLev']))[0]

    for jj in j:
      # add the emissivity to the second of the 2
      llist['Epsilon'][jj+1] += llist['Epsilon'][jj]

      # check if we are at the end of the array. If so, make this a good level
      keep[jj+1]=True

    # remove all the non-keepers
    llist = llist[keep]
    # fix the emissivities if log scaled
    if log_interp:
      llist['Epsilon'] = numpy.exp(llist['Epsilon'])

    return llist

  def return_linelist(self, Te, teunit='keV', nearest = False,\
                      specrange=False, specunit='A', elements=False, abundance=False,\
                      log_interp=True):

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
    specunit : {'Angstrom','keV'}
      Units for specrange (default A)
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundance : dict(float)
      The abundances of each element, e.g. abund[6]=1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.

    Returns
    -------
    linelist : numpy.array(dtype=linelist_cie_spectrum)
      The list of lines at temperature Te.

    """

    # get kT in keV
    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_CIE+1)


    if abundance == False:
      abundance = {}
      for Z in elements:
        abundance[Z] = 1.0

    linelist = numpy.zeros(0, dtype = apec.generate_datatypes('linelist_cie_spectrum'))


    for Z in elements:
      abund = abundance[Z]
      if abund > 0:

        if len(ikT) > 1:


          llist1 = self.spectra[ikT[0]][Z].return_linelist(specrange,\
                                  specunit=specunit)
          llist2 = self.spectra[ikT[1]][Z].return_linelist(specrange,\
                                  specunit=specunit)

          elemlinelist = self._merge_linelists_temperatures(f, llist1, llist2, log_interp)

        else:
          elemlinelist = self.spectra[ikT[0]][Z].return_linelist(specrange,\
                                  specunit=specunit)

        # fix abundance
        elemlinelist['Epsilon'] *= abund

        # append this element's line list onto the total line list

        if len(linelist)==0:
          linelist=elemlinelist
        else:
          linelist =  numpy.append(linelist, elemlinelist)

    return linelist



class _ElementSpectrum():
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

  Attributes
  ----------
  lines : _LineData
    A _LineData object containing all the line information
  continuum : _ContinuumData
    A _ContinuumData object containing all the contrinuum information

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
    z1_drv : int
      the ion charge of the driving ion (0 for whole element)

    Modifies
    --------
    self.lines : _LineData
      The line emission HDU data for this element/ion
    self.continuum : _ContinuumData
      The continuum emission HDU data for this element/ion
    """

    # intialize
    if z1_drv != 0:
      tmp = linedata[(linedata['Element'] == Z) &\
                               (linedata['Ion_drv'] == z1_drv)]
      self.lines = _LineData(tmp)
      self.continuum = _ContinuumData(cocodata)

    else:
      self.lines = _LineData(linedata)
      self.continuum = _ContinuumData(cocodata)


  def return_spectrum(self, eedges, Te, ebins_checksum=False,\
                    thermal_broadening=False,\
                    broaden_limit=False,\
                    velocity_broadening=0.0,\
                    teunit = 'keV',\
                    broaden_object=False,\
                    dolines = True,\
                    docont = True,\
                    dopseudo = True):
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
    broaden_object : class
      Object with routine "broaden" which applies line broadening. Usually a Gaussian.
    dolines : bool
      Calculate line emission (default True)
    docont : bool
      Calculate Continuum emission (default True)
    dopseudo : bool
      Calculate PseudoContinuum (weak line) emission (default True)

    Returns
    -------
    spectrum : array(float)
      The spectrum in ph cm^3 s^-1 bin^-1
    """

    T = util.convert_temp(Te, teunit, 'K')

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
    spec = numpy.zeros(len(eedges)-1)

    if dolines:
      spec+=self.lines.return_spec(eedges, T, ebins_checksum=ebins_checksum,\
                                  thermal_broadening=thermal_broadening,\
                                  broaden_limit=broaden_limit,\
                                  velocity_broadening=velocity_broadening,\
                                  broaden_object=broaden_object)
    if dopseudo+docont > 0:
      spec+=self.continuum.return_spec(eedges, ebins_checksum = ebins_checksum,\
                                       dopseudo=dopseudo, docont=docont)

    self.spectrum = spec

    return self.spectrum

  def return_linelist(self, specrange, specunit='A'):
    """
    Return the list of lines in specrange

    Parameters
    ----------
    specrange : array
      spectral range to look for lines in
    specunit : string
      units of spectral range (A or keV)

    Returns
    -------
    linelist : array
      list of lines and epsilons
    """
    wave = util.convert_spec(specrange, specunit, 'A')

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

class _LineData():
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
    self.spectrum_calculated = False
    self.T = 0.0
    self.v = 0.0
    self.ebins_checksum = False


  def return_spec(self, eedges, T, ebins_checksum = False,\
                  thermal_broadening = False, \
                  velocity_broadening = 0.0, \
                  broaden_limit = 1e-18,\
                  broaden_object=False):
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
    broaden_object : class
      Object with routine "broaden" which applies line broadening. Usually a Gaussian.

    Returns
    -------
    spectrum : array(float)
      Emissivity on eedges spectral bins of the lines, in ph cm^3 s^-1 bin^-1
    """
    if ebins_checksum == False:
        # generate the checksum
      ebins_checksum = hashlib.md5(eedges).hexdigest()
    if velocity_broadening==None:
      velocity_broadening=0.0
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
          Tb = util.convert_temp(T, 'K','keV')*const.ERG_KEV/(masslist[llist['Element']]*1e3*const.AMUKG)

        if velocity_broadening <0:
          velocitybroadeining = 0.0
          vb=0.0
        else:
          vb = (velocity_broadening * 1e5)**2

        wcoeff = numpy.sqrt(Tb+vb) / (const.LIGHTSPEED*1e2)

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
        t0 = time.time()
        for iline in igood:

          spec += broaden_object.broaden(const.HC_IN_KEV_A/llist['Lambda'][iline],\
                           width[iline],eedges)*llist['Epsilon'][iline]

        t1 = time.time()
#        print("Broadeninging %i lines, in %f seconds"%(len(igood), t1-t0))
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


class _ContinuumData():
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
    Calculate Continuum emission (default True)
  dopseudo : bool
    Calculate PseudoContinuum (weak line) emission (default True)
  """
  def __init__(self, cocoentry):
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

    self.spectrum_calculated = False
    self.ebins_checksum = False


  def return_spec(self, eedges, ebins_checksum = False, docont=True, dopseudo=True):
    import scipy.integrate


    # get the checksum for the ebins, if not provided
    if ebins_checksum == False:
        # generate the checksum
      ebins_checksum = hashlib.md5(eedges).hexdigest()

#    else:
    if docont:
      cont = _expand_E_grid(eedges, len(self.ECont), self.ECont, self.Cont)
    else:
      cont = 0.0

    if dopseudo:
      pseudo = _expand_E_grid(eedges, len(self.EPseudo), self.EPseudo, self.Pseudo)
    else:
      pseudo = 0.0

    self.ebins_checksum = ebins_checksum
    self.spectrum = cont+pseudo

    return self.spectrum

# def apply_response(spectrum, rmf, arf=False):
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
# #
# # Update 2016-05-25
# #
# # Changed to return the energy grid and the spectrum, as apparently in some
# # instruments these are not the same as the input energy grid.
  # if arf:
    # if type(arf)==str:
      # arfdat = pyfits.open(arf)
    # elif type(arf) == pyfits.hdu.hdulist.HDUList:
      # arfdat = arf
    # else:
      # print("ERROR: unknown arf type, %s"%(repr(type(arf))))
      # return
    # res = spectrum * arfdat['SPECRESP'].data['SPECRESP']
  # else:
    # res = spectrum*1.0

  # sret = numpy.matmul(res,matrix)
  # return ebins, ret, sret





class NEISession(CIESession):
  """
  Load and generate a Non-equilibrium ionization (NEI) spectrum

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : iterable of int
    Elements to include, listed by atomic number. if not set, include all.
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  datacache : dict
    Any Atomdb FITS files which have to be opened are stored here
  spectra : NEISpectra
    Object storing the actual spectral data
  elements : iterable of int
    Elements to include, listed by atomic number. if not set, include all.
  default_abundset : string
    The abundance set used for the original emissivity file calculation
  abundset : string
    The abundance set to be used for the returned spectrum
  abundsetvector : array_like(float)
    The relative abundance between default_abundset and abundset for each element
  response_set : bool
    Have we loaded a response (or set a dummy response)
  dolines : bool
    Calculate line emission (default True)
  docont : bool
    Calculate Continuum emission (default True)
  dopseudo : bool
    Calculate PseudoContinuum (weak line) emission (default True)
  broaden_limit : float
    Apply broadening to lines with epsilon > this value (ph cm3 s-1)
  thermal_broadening : bool
    Apply thermal broadening to lines (default = False)
  velocity_broadening : float
    Apply velocity broadening with this velocity (km/s). If <=0, do not apply.

  Examples
  --------

  Create a session instance:

  >>> s = NEISession()

  Set up the responses, in this case a dummy response from 0.1 to 10 keV

  >>> ebins = numpy.linspace(0.1,10,1000)
  >>> s.set_response(ebins, raw=True)

  Turn on thermal broadening

  >>> s.set_broadening(True)
  Will thermally broaden lines with emissivity > 1.000000e-18 ph cm3 s-1

  Return spectrum at 1.0keV, Ne *t = 1e11 cm^-3 s

  >>> spec = s.return_spectrum(1.0, 1e11)

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
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance
      for list of options.
    """
    self.SessionType='NEI'

    self._session_initialise1(linefile, cocofile, elements, abundset)

    self.spectra=_NEISpectrum(self.linedata, self.cocodata)

    self._session_initialise2()


  def return_linelist(self,Te, tau, specrange, init_pop='ionizing',specunit='A', \
                               teunit='keV', apply_aeff=False, by_ion_drv = False,\
                               nearest=False, apply_binwidth=False,\
                               log_interp=True, freeze_ion_pop=False):
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
    by_ion_drv : bool
      If true, keep lines which are the same but have different ion_drv separate.
      Otherwise, merge them.
    nearest : bool
      calculate spectrum at nearest temperature in linelist, no
      interpolation. Ionization fraction calculation still exact.
    apply_binwidth :
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.


    Returns
    -------
    linelist : array(linelist)
      The list of lines with lambda (A), energy (keV), epsilon (ph cm3 s-1),\
      epsilon_aeff (ph cm5 s-1) ion (string) and upper & lower levels.

    """

    kT = util.convert_temp(Te, teunit, 'keV')

    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]



    s= self.spectra.return_linelist(kT, tau, init_pop=init_pop, \
                                    specrange=specrange, teunit='keV',\
                                    specunit=specunit, elements=self.elements,\
                                    abundance = ab, by_ion_drv = by_ion_drv,\
                                    log_interp=log_interp, \
                                    freeze_ion_pop = freeze_ion_pop)

    # do the response thing
    #resp  = s.response()

    if apply_aeff == True:

      epsilon_aeff =  self._apply_linelist_aeff(s, specunit, apply_binwidth)

      s['Epsilon'] = epsilon_aeff
    return(s)

  def return_line_emissivity(self, Telist, taulist, Z, z1, up, lo, \
                             specunit='A', teunit='keV', \
                             apply_aeff=False, apply_abund=True,\
                             log_interp = True, init_pop='ionizing',\
                             freeze_ion_pop=False):
    """
    Get line emissivity as function of Te, tau. Assumes ionization from neutral.


    Parameters
    ----------
    Telist : float or array(float)
      Temperature(s) in keV or K
    taulist : float or array(float)
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
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.

    Returns
    -------
    ret : dict
      Dictionary containing:
      Te, tau, teunit: as input
      wavelength : line wavelength (A)
      energy : line energy (keV)
      epsilon : emissivity in ph cm^3 s-1 (or ph cm^5 s^-1 if apply_aeff=True)
                first index is temperature, second is tau. If Te or Tau was
                supplied as a scalar, then that index is removed

    """

    Tevec, Teisvec = util.make_vec(Telist)
    tauvec, tauisvec = util.make_vec(taulist)


    kTlist = util.convert_temp(Tevec, teunit, 'keV')
    if apply_abund:
      ab = self.abund[Z]*self.abundsetvector[Z]
    else:
      ab = 1.0

    eps = numpy.zeros([len(Tevec), len(tauvec)])
    ret={}
    ret['wavelength'] = None
    for itau, tau in enumerate(tauvec):
      for ikT, kT in enumerate(kTlist):
        e, lam = self.spectra.return_line_emissivity(kT, tau, Z, z1, \
                                                     up, lo, \
                                                     specunit='A', \
                                                     teunit='keV', \
                                                     abundance=ab,\
                                                     init_pop=init_pop,\
                                                     freeze_ion_pop = freeze_ion_pop)

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

  def return_spectrum(self,  Te, tau, init_pop='ionizing', teunit='keV', nearest=False,\
                      log_interp=True, freeze_ion_pop=False):
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
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    nearest : bool
      If set, return the spectrum from the nearest tabulated temperature
      in the file, without interpolation
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.

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
                                    teunit=teunit, \
                                    nearest = nearest,elements = el_list, \
                                    abundance=ab, log_interp=True,\
                                    broaden_object=self.cdf, \
                                    freeze_ion_pop = freeze_ion_pop,\
                                    do_eebrems=self.do_eebrems)

    ss = self._apply_response(s)
    self.ionfrac = self.spectra.ionfrac

    return ss



class _NEISpectrum(_CIESpectrum):
  """
  A class holding the emissivity data for NEI emission, and returning
  spectra

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)

  Attributes
  ----------
  SessionType : string
    "NEI"
  spectra : dict of _ElementSpectrum
    a dictionary containing the emissivity data for each HDU,
    subdivided by ion (spectra[12][18][13] is an _ElementSpectrum object
    containing the argon XIII data for the 12th HDU)
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
    linedata : HDUList
      The line emissivity data
    cocodata : HDUList
      The continuum emissivity data
    """

    self.datacache={}
    self.SessionType = 'NEI'

    picklefname = os.path.expandvars('$ATOMDB/spectra_%s_%s.pkl'%\
                                (linedata[0].header['CHECKSUM'],\
                                 cocodata[0].header['CHECKSUM']))

    havepicklefile = False
    if os.path.isfile(picklefname):
      havepicklefile = True

    if havepicklefile:
      try:
        self.spectra = pickle.load(open(picklefname,'rb'))
        self.kTlist = self.spectra['kTlist']
      except AttributeError:
        havepicklefile=False
        print("pre-stored data in %s is out of date. This can be caused by updates to the data "%(picklefname)+
              "or, more likely, changes to pyatomdb. Regenerating...")

        # delete the old file
        if os.path.isfile(picklefname):
          os.remove(picklefname)

    if not havepicklefile:
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
            self.spectra[ihdu][Z][z1]=_ElementSpectrum(ldat[isgood],\
                                                  ccdat[0], Z, z1_drv=z1)


      pickle.dump(self.spectra, open(picklefname,'wb'))


    self.logkTlist=numpy.log(self.kTlist)


  def _calc_ionfrac(self, Te, tau, init_pop='ionizing', teunit='keV', freeze_ion_pop = False,\
                   elements=False):
    """
    Calculate the ion fractions in an NEI case

    Parameters
    ----------
    Te : float
      Electron temperature (default, keV)
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : string
      Units of kT (keV by default, K also allowed)
    freeze_ion_pop : bool
      If True, return the initial ionization fraction as the final
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.


    Returns
    -------
    ionfrac : dict of arrays
      Array of all the ion fractions.
  """
    init_pop_calc={}
    if elements==False:
      elements = self.elements
    # check the format of init_pop
    if isinstance(init_pop, str):
      if init_pop.lower() == 'ionizing':
        for Z in elements:
          init_pop_calc[Z] = numpy.zeros(Z+1)
          init_pop_calc[Z][0] = 1.0
      elif init_pop.lower() == 'recombining':
        for Z in elements:
          init_pop_calc[Z] = numpy.zeros(Z+1)
          init_pop_calc[Z][-1] = 1.0
      else:
        raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
             (init_pop))
    elif isinstance(init_pop, float):
      # this is an initial temperature
      kT_init = util.convert_temp(init_pop, teunit, 'keV')
      for Z in elements:
        init_pop_calc[Z] = apec.return_ionbal(Z, kT_init, \
                                            teunit='keV', \
                                            datacache=self.datacache)
    elif isinstance(init_pop, dict):
      for Z in elements:
        init_pop_calc[Z] = init_pop[Z]
    else:
      raise util.OptionError("Error: invalid type for init_pop")


    ionfrac = {}

    if freeze_ion_pop:
      for Z in elements:
        ionfrac[Z] = init_pop_calc[Z]

    else:
    # no calculate the output
      kT = util.convert_temp(Te, teunit, 'keV')
      for Z in elements:
        ionfrac[Z] = apec.return_ionbal(Z, kT, init_pop=init_pop_calc[Z], \
                                          tau=tau, \
                                          teunit='keV', \
                                          datacache=self.datacache)
    return ionfrac


  def return_spectrum(self, Te, tau, init_pop='ionizing', teunit='keV', nearest = False,
                             elements=False, abundance=False, log_interp=True, broaden_object=False,\
                             freeze_ion_pop = False, do_eebrems = False):

    """
    Return the spectrum of the element on the energy bins in
    self.session.specbins

    Parameters
    ----------
    Te : float
      Electron temperature (default, keV)
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : string
      Units of kT (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundance : dict(float)
      The abundances of each element, e.g. abund[6]=1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    broaden_object : class
      Object with routine "broaden" which applies line broadening. Usually a Gaussian.
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.
    do_eebrems : bool
      Calculate electron-electron bremsstrahlung emission (default False)

    Returns
    -------
    spec : array(float)
      The element's emissivity spectrum, in photons cm^3 s^-1 bin^-1
    """

    # get kT in keV
    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_NEI+1)


    if abundance == False:
      abundance = {}
      for Z in elements:
        abundance[Z] = 1.0

    totspec = 0.0

    # only do the elements where abundance > 0
    el = []
    for Z in elements:
      if abundance[Z]>0.0:
        el.append(Z)

    ionfrac_all = self._calc_ionfrac(kT, tau, init_pop=init_pop, teunit='keV', \
                                    freeze_ion_pop = freeze_ion_pop,\
                                    elements = el)
    self.ionfrac = ionfrac_all
    spec = {}
    if len(ikT)==2:
      spec[0] = 0.0
      spec[1] = 0.0
    else:
      spec[0] = 0.0

    for i in range(len(ikT)):

      ikTspec = 0.0

      for Z in elements:

        abund = abundance[Z]
        if abund > 0:
          ionfrac=ionfrac_all[Z]

          elspec = 0.0

          for z1 in range(1, Z+2):
            ionspec = 0.0

            if ionfrac[z1-1]>1e-10:
              # calculate minimum emissivitiy to broaden, accounting for ion
              # and element abundance.
              epslimit =  self.broaden_limit/(abund*ionfrac[z1-1])

              # return a broadened spectrum

              ionspec = self.spectra[ikT[i]][Z][z1].return_spectrum(self.ebins,\
                                  kT,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening,\
                                  broaden_object=broaden_object) *\
                                  ionfrac[z1-1]



              spec[i] += ionspec*abund

    if len(ikT)==2:
      totspec = self._merge_spectra_temperatures(f, spec[0], spec[1], log_interp)
    else:
      totspec = spec[0]

    if do_eebrems:
      nel = 0.0
      rawabund = atomdb.get_abundance(datacache=self.datacache)
      for Z in ionfrac_all.keys():
        if abundance[Z] > 0.0:
          nel += sum(ionfrac_all[Z]*numpy.arange(Z+1))*rawabund[Z]*abundance[Z]

      eespec = calc_ee_brems_spec(self.ebins, kT, nel)
      totspec += eespec
    return totspec


  def return_line_emissivity(self, Te, tau, Z, z1, up, lo, specunit='A',
                             teunit='keV', abundance=1.0,
                             log_interp = True, init_pop = 'ionizing',\
                             freeze_ion_pop=False):
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
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.

    Returns
    -------
    Emissivity : float
      Emissivity in photons cm^3 s^-1
    spec : float
      Wavelength or Energy of line, depending on specunit
    """

    import collections

    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, \
                                     teunit='keV', \
                                     nearest=False, \
                                     log_interp=log_interp)
    #ikT has the 2 nearest temperature indexes
    # f has the fraction for each

    ionfrac_all = self._calc_ionfrac(kT, tau, init_pop=init_pop, teunit='keV', \
                                    freeze_ion_pop = freeze_ion_pop,\
                                    elements = [Z])


    ionfrac = ionfrac_all[Z]

    eps = 0.0
    lam = 0.0


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


    if specunit == 'keV':
      lam = const.HC_IN_KEV_A/lam
    return eps, lam

  def return_linelist(self,  Te, tau,init_pop='ionizing',
                      teunit='keV', nearest = False, specrange=False,
                      specunit='A', elements=False, abundance=False,\
                      log_interp=True, freeze_ion_pop=False, by_ion_drv = False):

    """
    Return the lines in a spectral range for a given Te, tau

    Parameters
    ----------

    Te : float
      Electron temperature (default, keV)
    tau : float
      ionization timescale, ne * t (cm^-3 s).
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
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
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundance : dict(float)
      The abundances of each element, e.g. abund[6]=1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.
    by_ion_drv : bool
      If true, keep lines which are the same but have different ion_drv separate.
      Otherwise, merge them.

    """
    # get kT in keV
    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_NEI+1)


    if abundance == False:
      abundance = {}
      for Z in elements:
        abundance[Z] = 1.0

    totspec = 0.0

    # only do the elements where abundance > 0
    el = []
    for Z in elements:
      if abundance[Z]>0.0:
        el.append(Z)

    ionfrac_all = self._calc_ionfrac(kT, tau, init_pop=init_pop, teunit='keV', \
                                    freeze_ion_pop = freeze_ion_pop,\
                                    elements = el)

    if by_ion_drv:
      linelist = numpy.zeros(0, dtype=apec.generate_datatypes('linelist_nei_spectrum'))

    else:
      linelist = numpy.zeros(0, dtype=apec.generate_datatypes('linelist_cie_spectrum'))

    for Z in elements:
      abund = abundance[Z]

      #skip if none of element present
      if abund > 0:

        haveelemlist = False
        #elemllist = numpy.zeros(0,dtype=apec.generate_datatypes('linelist_nei_spectrum'))
        #cycle through each driving ion

        ionfrac = ionfrac_all[Z]
        for z1_drv in range(1, Z+2):

          # skip if none of driving ion present
          if ionfrac[z1_drv-1] <1e-10: continue


          # if > 1 temperature, calculate for each
          if len(ikT) > 1:

            # get the linelist
            llist1 = self.spectra[ikT[0]][Z][z1_drv].return_linelist(specrange,\
                                    specunit=specunit)
            llist2 = self.spectra[ikT[1]][Z][z1_drv].return_linelist(specrange,\
                                    specunit=specunit)
            # correct for ion fraction
            llist1['Epsilon'] *= ionfrac[z1_drv-1]
            llist2['Epsilon'] *= ionfrac[z1_drv-1]

            # create linelists for each temperature, or add to them
            if haveelemlist == False:
              elemllist1 = llist1
              elemllist2 = llist2
              haveelemlist = True

            else:
              elemllist1 = numpy.append(elemllist1, llist1)
              elemllist2 = numpy.append(elemllist2, llist2)

          else:
            # only 1 temperature

            # get the linelist
            llist = self.spectra[ikT[0]][Z][z1].return_linelist(specrange,\
                                     specunit=specunit)

            # correct for ion fraction
            llist['Epsilon'] *= ionfrac[z1_drv-1]

            # create linelists for each temperature, or add to them
            if haveelemlist == False:
              elemllist = llist
              haveelemlist = True
            else:
              elemllist = numpy.append(elemllist, llist)

        # now have a complete ionllist for each ion. Merge them and
        # remove duplicates

        if len(ikT) > 1:
          # merge the two temperatures
          elemllist = self._merge_linelists_temperatures(f, elemllist1, elemllist2, \
                                                             log_interp, \
                                                             by_ion_drv=by_ion_drv)
        else:
          # merge duplicate lines
          elemllist = self._merge_linelist_duplicates(elemllist, by_ion_drv=by_ion_drv)


        # fix abundance
        elemllist['Epsilon'] *= abund

        # append this element's line list onto the total line list

        if by_ion_drv:
          # appending all the colums, so easy
          linelist = numpy.append(linelist, elemllist)
        else:
          # Not keep the drv columns, so a little trickier
          linelisttmp = numpy.zeros(len(elemllist), dtype = linelist.dtype)
          for key in linelist.dtype.names:
            linelisttmp[key] = elemllist[key]
          linelist = numpy.append(linelist, linelisttmp)

    return linelist


class PShockSession(NEISession):
  """
  Load and generate a Parallel Shock (PShock) spectrum

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)
  elements : iterable of int
    Elements to include, listed by atomic number. if not set, include all.
  abundset : string
    The abundance set to use. Default AG89.

  Attributes
  ----------
  datacache : dict
    Any Atomdb FITS files which have to be opened are stored here
  spectra : PShockSpectra
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
    Calculate line emission (default True)
  docont : bool
    Calculate Continuum emission (default True)
  dopseudo : bool
    Calculate PseudoContinuum (weak line) emission (default True)
  broaden_limit : float
    Apply broadening to lines with epsilon > this value (ph cm3 s-1)
  thermal_broadening : bool
    Apply thermal broadening to lines (default = False)
  velocity_broadening : float
    Apply velocity broadening with this velocity (km/s). If <=0, do not apply.

  Examples
  --------

  Create a session instance:

  >>> s=PShockSession()

  Set up the responses, in this case a dummy response from 0.1 to 10 keV

  >>> ebins = numpy.linspace(0.1,10,1000)
  >>> s.set_response(ebins, raw=True)

  Turn on thermal broadening

  >>> s.set_broadening(True)
  Will thermally broaden lines with emissivity > 1.000000e-18 ph cm3 s-1

  Return spectrum at 1.0keV, Tau_u = 1e11, Tau_l = 0.0 (default)

  >>> spec = s.return_spectrum(1.0, 1e11)

  Same with Tau_l = 1e9

  >>> spec = s.return_spectrum(1.0, 1e11, tau_l = 1e9)


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
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance
      for list of options.
    """

    self.SessionType='PShock'

    self._session_initialise1(linefile, cocofile, elements, abundset)

    self.spectra=_PShockSpectrum(self.linedata, self.cocodata)

    self._session_initialise2()

    # self.datacache={}

    # # Open up the APEC files
    # self.set_apec_files(linefile, cocofile)


    # # if elements are specified, use them. Otherwise, use Z=1-30
    # if util.keyword_check(elements):
      # self.elements = elements
    # else:
      # self.elements=list(range(1,const.MAXZ_NEI+1))

    # # a hold for the spectra
    # self.spectra=PShockSpectrum(self.linedata, self.cocodata)


    # # Set both the current and the default abundances to those that
    # # the apec data was calculated on
    # self.abundset=self.linedata[0].header['SABUND_SOURCE']
    # self.default_abundset=self.linedata[0].header['SABUND_SOURCE']

    # self.abundsetvector = numpy.zeros(const.MAXZ_NEI+1)
    # for Z in self.elements:
      # self.abundsetvector[Z] = 1.0

    # #  but if another vector was already specified, use this instead
    # if util.keyword_check(abundset):
      # self.set_abundset(abundset)

    # self.abund = numpy.zeros(const.MAXZ_NEI+1)

    # for Z in self.elements:
      # self.abund[Z]=1.0

    # # Set a range of parameters which can be overwritten later
    # self.response_set = False # have we loaded a response file?
    # self.dolines=True # Include lines in spectrum
    # self.docont=True # Include continuum in spectrum
    # self.dopseudo=True # Include pseudo continuum in spectrum
    # self.set_broadening(False, broaden_limit=1e-18)
    # self.cdf = _Gaussian_CDF()



  def return_linelist(self,Te, tau_u, specrange, tau_l = 0.0, init_pop='ionizing',specunit='A', \
                               teunit='keV', apply_aeff=False, by_ion_drv = False,\
                               nearest=False, apply_binwidth=False):
    """
    Get the list of line emissivities vs wavelengths


    Parameters
    ----------
    Te : float
      Temperature in keV or K
    tau_u : float
      Upper limit of ionization timescale, ne * t (cm^-3 s).
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    tau_l :
      lower limit of ionization timescale, ne * t (cm^-3 s) (defaults to 0)
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    specunit : {'Angstrom','keV'}
      Units for specrange
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    apply_aeff : bool
      If true, apply the effective area to the lines in the linelist to
      modify their intensities.
    by_ion_drv : bool
      If true, keep lines which are the same but have different ion_drv separate.
      Otherwise, merge them.
    nearest :
    apply_binwidth :

    Returns
    -------
    linelist : array(dtype)
      The list of lines with lambda (A), energy (keV), epsilon (ph cm3 s-1),\
      epsilon_aeff (ph cm5 s-1) ion (string) and upper & lower levels.

    """

    kT = util.convert_temp(Te, teunit, 'keV')

    el_list = self.elements
    ab = {}
    for Z in el_list:
      ab[Z] = self.abund[Z]*self.abundsetvector[Z]



    s= self.spectra.return_linelist(kT, tau_u, tau_l=tau_l, init_pop=init_pop, \
                                    specrange=specrange, teunit='keV',\
                                    specunit=specunit, elements=self.elements,\
                                    abundance = ab, by_ion_drv = by_ion_drv)

    # do the response thing
    #resp  = s.response()

    if apply_aeff == True:

      epsilon_aeff =  self._apply_linelist_aeff(s, specunit, apply_binwidth)

      s['Epsilon_Err'] = epsilon_aeff
    return(s)


  def return_spectrum(self,  Te, tau_u, tau_l=0.0, init_pop='ionizing',  teunit='keV', nearest=False,\
                      log_interp=True):
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
    tau_u : float
      Upper limit of ionization timescale, ne * t (cm^-3 s).
    tau_l :
      lower limit of ionization timescale, ne * t (cm^-3 s) (defaults to 0)
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : {'keV' , 'K'}
      Units of te (kev or K, default keV)
    nearest : bool
      If set, return the spectrum from the nearest tabulated temperature
      in the file, without interpolation
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.

    Returns
    -------
    spectrum : array(float)
      The spectrum in photons cm^5 s^-1 bin^-1, with the response, or
      photons cm^3 s^-1 bin^-1 if raw is set.
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
    s= self.spectra.return_spectrum(Te, tau_u, tau_l=tau_l, init_pop=init_pop, \
                                    teunit=teunit, \
                                    nearest = nearest,elements = el_list, \
                                    abundance=ab, log_interp=True,\
                                    broaden_object=self.cdf,\
                                    do_eebrems=self.do_eebrems)

    ss = self._apply_response(s)

    return ss



  def return_line_emissivity(self, Telist, tau_ulist, Z, z1, up, lo, \
                             tau_llist = 0.0, specunit='A', teunit='keV', \
                             apply_aeff=False, apply_abund=True,\
                             log_interp = True, init_pop='ionizing'):
    """
    Get line emissivity as function of Te, tau. Assumes ionization from neutral.


    Parameters
    ----------
    Telist : float or array(float)
      Temperature(s) in keV or K
    tau_ulist : float
      ionization timescale(s), ne * t (cm^-3 s).
    Z : int
      nuclear charge of element
    z1 : int
      ion charge +1 of ion
    up : int
      upper level for transition
    lo : int
      lower level for transition
    tau_llist : float
      lower limit of ionization timescale(s), ne * t (cm^-3 s).
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
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.

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
    tau_uvec, tau_uisvec = util.make_vec(tau_ulist)
    tau_lvec, tau_lisvec = util.make_vec(tau_llist)


    kTlist = util.convert_temp(Tevec, teunit, 'keV')
    if apply_abund:
      ab = self.abund[Z]*self.abundsetvector[Z]
    else:
      ab = 1.0

    eps = numpy.zeros([len(Tevec), len(tau_uvec), len(tau_lvec)])
    ret={}
    ret['wavelength'] = None
    for itau_u, tau_u in enumerate(tau_uvec):
      for ikT, kT in enumerate(kTlist):
        for itau_l, tau_l in enumerate(tau_lvec):
          e, lam = self.spectra.return_line_emissivity(kT, tau_u, Z, z1, \
                                                       up, lo, \
                                                       tau_l = tau_l,\
                                                       specunit='A', \
                                                       teunit='keV', \
                                                       abundance=ab,\
                                                       init_pop=init_pop)

          eps[ikT, itau_u, itau_l] = e
          if lam != False:
            ret['wavelength'] = lam * 1.0
#          else:
 #           ret['wavelength'] = None

    ret['Te'] = Telist
    ret['tau_u'] = tau_ulist
    ret['tau_l'] = tau_llist
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

    if not tau_lisvec:
      eps = eps.sum(2)

    if not tau_uisvec:
      eps = eps.sum(1)

    if not Teisvec:
      eps = eps.sum(0)

    ret['epsilon'] = eps

    return ret




class _PShockSpectrum(_NEISpectrum):
  """
  A class holding the emissivity data for NEI emission, and returning
  spectra

  Parameters
  ----------
  linefile : string or HDUList, optional
    The line emissivity data file (either name or already open)
  cocofile : string or HDUList, optional
    The continuum emissivity data file (either name or already open)

  Attributes
  ----------
  session : CIESession
    The parent CIESession
  SessionType : string
    "CIE"
  spectra : dict of _ElementSpectrum
    a dictionary containing the emissivity data for each HDU,
    subdivided by element (spectra[12][18] is an _ElementSpectrum object
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
  linedata : HDUList
    The line emissivity data
  cocodata : HDUList
    The continuum emissivity data
    """

    self.datacache={}
    self.SessionType = 'PShock'

    picklefname = os.path.expandvars('$ATOMDB/spectra_%s_%s.pkl'%\
                                (linedata[0].header['CHECKSUM'],\
                                 cocodata[0].header['CHECKSUM']))


    havepicklefile = False
    if os.path.isfile(picklefname):
      havepicklefile = True

    if havepicklefile:
      try:
        self.spectra = pickle.load(open(picklefname,'rb'))
        self.kTlist = self.spectra['kTlist']
      except AttributeError:
        havepicklefile=False
        print("pre-stored data in %s is out of date. This can be caused by updates to the data "%(picklefname)+
              "or, more likely, changes to pyatomdb. Regenerating...")

        # delete the old file
        if os.path.isfile(picklefname):
          os.remove(picklefname)

    if not havepicklefile:
      self.spectra={}
      self.kTlist = numpy.array(linedata[1].data['kT'].data)
      self.spectra['kTlist'] = numpy.array(linedata[1].data['kT'].data)
      for ihdu in range(len(self.kTlist)):
        self.spectra[ihdu]={}
        self.spectra[ihdu]['kT'] = self.kTlist[ihdu]
        ldat = numpy.array(linedata[ihdu+2].data.data)
        cdat = numpy.array(cocodata[ihdu+2].data.data)


        Zarr = numpy.zeros([len(ldat), const.MAXZ_NEI+1], dtype=bool)
        print(util.unique(ldat['Element']))
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
            self.spectra[ihdu][Z][z1]=_ElementSpectrum(ldat[isgood],\
                                                  ccdat[0], Z, z1_drv=z1)


      pickle.dump(self.spectra, open(picklefname,'wb'))


    self.logkTlist=numpy.log(self.kTlist)



  def _calc_ionfrac(self, Te, tau_u, tau_l=0.0, init_pop='ionizing', teunit='keV', freeze_ion_pop = False,\
                   elements=False):
    """
    Calculate the ion fractions in an NEI case

    Parameters
    ----------
    Te : float
      Electron temperature (default, keV)
    tau_u : float
      Upper limit of ionization timescale, ne * t (cm^-3 s).
    tau_l :
      lower limit of ionization timescale, ne * t (cm^-3 s) (defaults to 0)
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : string
      Units of kT (keV by default, K also allowed)
    freeze_ion_pop : bool
      If True, return the initial ionization fraction as the final
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.

    Returns
    -------
    ionfrac : dict of arrays
      Array of all the ion fractions.
  """
    init_pop_calc={}
    if elements==False:
      elements = self.elements
    # check the format of init_pop
    if isinstance(init_pop, str):
      if init_pop.lower() == 'ionizing':
        for Z in elements:
          init_pop_calc[Z] = numpy.zeros(Z+1)
          init_pop_calc[Z][0] = 1.0
      elif init_pop.lower() == 'recombining':
        for Z in elements:
          init_pop_calc[Z] = numpy.zeros(Z+1)
          init_pop_calc[Z][-1] = 1.0
      else:
        raise util.OptionError("Error: init_pop is set as a string, must be 'ionizing' or 'recombining'. Currently %s."%\
             (init_pop))
    elif isinstance(init_pop, float):
      # this is an initial temperature
      kT_init = util.convert_temp(init_pop, teunit, 'keV')
      for Z in elements:
        init_pop_calc[Z] = apec.return_ionbal(Z, kT_init, \
                                            teunit='keV', \
                                            datacache=self.datacache)
    elif isinstance(init_pop, dict):
      for Z in elements:
        init_pop_calc[Z] = init_pop[Z]
    else:
      raise util.OptionError("Error: invalid type for init_pop: ", init_pop)


    ionfrac = {}
    kT = util.convert_temp(Te, teunit, 'keV')

    if freeze_ion_pop:
      for Z in elements:
        ionfrac[Z] = init_pop_calc[Z]

    else:
      #####
      if tau_l < 1.0:
        nzones = 200

        taulist = numpy.logspace(numpy.log10(1e8), numpy.log10(tau_u),nzones+1)
        taulist = numpy.append(0, taulist[:-1])
        weight = (taulist[1:]-taulist[:-1])/tau_u
        taulist = (taulist[1:]+taulist[:-1])/2

      elif tau_l == tau_u:
        nzones=1
        weight = numpy.array([1.0])
        taulist = numpy.array([tau_u])

      else:
        nzones=40
        taulist = numpy.linspace(tau_l, tau_u, nzones+1)[:-1]
        weight = numpy.zeros(nzones)
        weight[:] = 0.025
        taulist[1:] = 0.5*(taulist[1:]+taulist[:-1])

      for Z in elements:
        # now need to make new ion fractions

        # cycle through the assorted tau values, to get new ionfrac

        # tau_l is zero:

        #now calculate cumulative ionfrac

        ionfractmp = apec.return_ionbal(Z,kT,init_pop = init_pop_calc[Z],\
                                             tau = taulist,teunit='keV', \
                                            datacache=self.datacache)
        for i in range(nzones):
          ionfractmp[i,:]*=weight[i]

        ionfrac[Z] = numpy.sum(ionfractmp,0)
        elspec = 0.0

    return ionfrac


  def return_spectrum(self, Te, tau_u, tau_l=0.0, init_pop='ionizing', teunit='keV', nearest = False,
                             elements=False, abundance=False, log_interp=True, broaden_object=False,\
                             freeze_ion_pop=False, do_eebrems = False):

    """
    Return the spectrum of the element on the energy bins in
    self.session.specbins

    Parameters
    ----------
    Te : float
      Electron temperature (default, keV)
    tau_u : float
      Upper limit of ionization timescale, ne * t (cm^-3 s).
    tau_l :
      lower limit of ionization timescale, ne * t (cm^-3 s) (defaults to 0)
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : string
      Units of kT (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundance : dict(float)
      The abundances of each element, e.g. abund[6]=1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    broaden_object : class
      Object with routine "broaden" which applies line broadening. Usually a Gaussian.
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.
    do_eebrems : bool
      Calculate electron-electron bremsstrahlung emission (default False)

    Returns
    -------
    spec : array(float)
      The element's emissivity spectrum, in photons cm^3 s^-1 bin^-1
    """

    # get kT in keV
    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_NEI+1)


    if abundance == False:
      abundance = {}
      for Z in elements:
        abundance[Z] = 1.0

    totspec = 0.0

    ionfrac_all = self._calc_ionfrac(kT, tau_u, tau_l=tau_l, init_pop=init_pop, teunit=teunit, freeze_ion_pop=freeze_ion_pop, elements=elements)

    for Z in elements:

      abund = abundance[Z]
      if abund > 0:

        # find the initial ionfrac (before anything changes)



        ionfrac = ionfrac_all[Z]
        elspec = 0.0
        for i in range(len(ikT)):

          ikTspec = 0.0

          for z1 in range(1, Z+2):
            ionspec = 0.0

            if ionfrac[z1-1]>1e-10:
              # calculate minimum emissivitiy to broaden, accounting for ion
              # and element abundance.
              epslimit =  self.broaden_limit/(abund*ionfrac[z1-1])

              # return a broadened spectrum

              ionspec = self.spectra[ikT[i]][Z][z1].return_spectrum(self.ebins,\
                                  kT,\
                                  ebins_checksum = self.ebins_checksum,\
                                  thermal_broadening = self.thermal_broadening,\
                                  broaden_limit = epslimit,\
                                  velocity_broadening = self.velocity_broadening,\
                                  broaden_object=broaden_object) *\
                                  ionfrac[z1-1]



              ikTspec += ionspec

          if len(ikTspec) > 1:
            # add appropriately
            if log_interp:
              elspec += numpy.log(ikTspec+const.MINEPSOFFSET)*f[i]
            else:
              elspec += ikTspec*f[i]
          else:
            # if just one then no need to interpolate on log scale
            elspec += ikTspec*f[i]

        # now we have the spectrum of the whole element. Un-log scale if
        # required.

        if log_interp:
          elspec = numpy.exp(elspec)-const.MINEPSOFFSET

        # add the element's spectrum to the whole, including abundance
        totspec += elspec* abund


    if do_eebrems:
      nel = 0.0
      rawabund = atomdb.get_abundance(datacache=self.datacache)
      for Z in ionfrac_all.keys():
        if abundance[Z] > 0.0:
          nel += sum(ionfrac_all[Z]*numpy.arange(Z+1))*rawabund[Z]*abundance[Z]

      eespec = calc_ee_brems_spec(self.ebins, kT, nel)
      totspec += eespec

    return totspec


  def return_line_emissivity(self, Te, tau_u, Z, z1, up, lo, tau_l=0.0, specunit='A',
                             teunit='keV', abundance=1.0,
                             log_interp = True, init_pop = 'ionizing',\
                             freeze_ion_pop=False):
    """
    Return the emissivity of a line at kT, tau. Assumes ionization from neutral for now


    Parameters
    ----------
    Te : float
      Temperature in keV or K
    tau_u : float
      ionization timescale, ne * t (cm^-3 s).
    Z : int
      nuclear charge of element
    z1 : int
      ion charge +1 of ion
    up : int
      upper level for transition
    lo : int
      lower level for transition
    tau_l :
      lower limit of ionization timescale, ne * t (cm^-3 s) (defaults to 0)
    specunit : {'Angstrom','keV'}
      Units for wavelength or energy (a returned value)
    teunit : {'keV' , 'K'}
      Units of Telist (kev or K, default keV)
    abundance : float
      Abundance to multiply the emissivity by
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.

    Returns
    -------
    Emissivity : float
      Emissivity in photons cm^3 s^-1
    spec : float
      Wavelength or Energy of line, depending on specunit
    """

    import collections

    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, \
                                     teunit='keV', \
                                     nearest=False, \
                                     log_interp=log_interp)
    #ikT has the 2 nearest temperature indexes
    # f has the fraction for each

    ionfrac_all = self._calc_ionfrac(kT, tau_u, tau_l=tau_l, init_pop=init_pop, teunit='keV', \
                                    freeze_ion_pop = freeze_ion_pop,\
                                    elements = [Z])


    ionfrac = ionfrac_all[Z]

    eps = 0.0
    lam = 0.0


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


    if specunit == 'keV':
      lam = const.HC_IN_KEV_A/lam
    return eps, lam



  def return_linelist(self,  Te, tau_u, tau_l=0.0, init_pop='ionizing',
                      teunit='keV', nearest = False, specrange=False,
                      specunit='A', elements=False, abundance=False,\
                      log_interp=True, freeze_ion_pop=False, by_ion_drv = False):

    """
    Return the lines in a spectral range for a given Te, tau_u

    Parameters
    ----------

    Te : float
      Electron temperature (default, keV)
    tau_u : float
      Upper limit of ionization timescale, ne * t (cm^-3 s).
    tau_l :
      lower limit of ionization timescale, ne * t (cm^-3 s) (defaults to 0)
    init_pop :    init_pop : string or float
      if 'ionizing' : all ionizing from neutral (so [1,0,0,0...])
      if 'recombining': all recombining from ionized (so[...0,0,1])
      if dict of arrays : the acutal fractional populations (so init_pop[6] is the array for carbon)
      if single float : the temperature (same units as Te)
    teunit : string
      Units of Te (keV by default, K also allowed)
    nearest : bool
      If True, return spectrum for the nearest temperature index.
      If False, use the weighted average of the (log of) the 2 nearest indexes.
      default is False.
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Ansgtrom','keV'}
      Units for specrange (default A)
    elements : iterable of int
      Elements to include, listed by atomic number. if not set, include all.
    abundance : dict(float)
      The abundances of each element, e.g. abund[6]=1.1 means multiply carbon
      abundance by 1.1.
    log_interp : bool
      Perform linear interpolation on a logT/logEpsilon grid (default), or linear.
    freeze_ion_pop : bool
      If True, skip the ion population calculation, use init_pop as the final pop instead.
    by_ion_drv : bool
      If true, keep lines which are the same but have different ion_drv separate.
      Otherwise, merge them.

    """

    # get kT in keV
    kT = util.convert_temp(Te, teunit, 'keV')

    ikT, f = self.get_nearest_Tindex(kT, teunit='keV', nearest=nearest, log_interp=log_interp)

    # check the params:
    if elements==False:
      elements=range(1,const.MAXZ_NEI+1)


    if abundance == False:
      abundance = {}
      for Z in elements:
        abundance[Z] = 1.0

    totspec = 0.0

    # only do the elements where abundance > 0
    el = []
    for Z in elements:
      if abundance[Z]>0.0:
        el.append(Z)
    ionfrac_all = self._calc_ionfrac(kT, tau_u, tau_l=tau_l, init_pop=init_pop, teunit='keV', \
                                    freeze_ion_pop = freeze_ion_pop,\
                                    elements = el)

    if by_ion_drv:
      linelist = numpy.zeros(0, dtype=apec.generate_datatypes('linelist_nei_spectrum'))

    else:
      linelist = numpy.zeros(0, dtype=apec.generate_datatypes('linelist_cie_spectrum'))

    for Z in elements:
      abund = abundance[Z]

      #skip if none of element present
      if abund > 0:

        haveelemlist = False
        #elemllist = numpy.zeros(0,dtype=apec.generate_datatypes('linelist_nei_spectrum'))
        #cycle through each driving ion

        ionfrac = ionfrac_all[Z]
        for z1_drv in range(1, Z+2):

          # skip if none of driving ion present
          if ionfrac[z1_drv-1] <1e-10: continue


          # if > 1 temperature, calculate for each
          if len(ikT) > 1:

            # get the linelist
            llist1 = self.spectra[ikT[0]][Z][z1_drv].return_linelist(specrange,\
                                    specunit=specunit)
            llist2 = self.spectra[ikT[1]][Z][z1_drv].return_linelist(specrange,\
                                    specunit=specunit)
            # correct for ion fraction
            llist1['Epsilon'] *= ionfrac[z1_drv-1]
            llist2['Epsilon'] *= ionfrac[z1_drv-1]

            # create linelists for each temperature, or add to them
            if haveelemlist == False:
              elemllist1 = llist1
              elemllist2 = llist2
              haveelemlist = True

            else:
              elemllist1 = numpy.append(elemllist1, llist1)
              elemllist2 = numpy.append(elemllist2, llist2)

          else:
            # only 1 temperature

            # get the linelist
            llist = self.spectra[ikT[0]][Z][z1].return_linelist(specrange,\
                                     teunit='keV', specunit=specunit)

            # correct for ion fraction
            llist['Epsilon'] *= ionfrac[z1_drv-1]

            # create linelists for each temperature, or add to them
            if haveelemlist == False:
              elemllist = llist
              haveelemlist = True
            else:
              elemllist = numpy.append(elemllist, llist)

        # now have a complete ionllist for each ion. Merge them and
        # remove duplicates

        if len(ikT) > 1:
          # merge the two temperatures
          elemllist = self._merge_linelists_temperatures(f, elemllist1, elemllist2, \
                                                             log_interp, \
                                                             by_ion_drv=by_ion_drv)
        else:
          # merge duplicate lines
          elemllist = self._merge_linelist_duplicates(elemllist, by_ion_drv=by_ion_drv)


        # fix abundance
        elemllist['Epsilon'] *= abund

        # append this element's line list onto the total line list

        if by_ion_drv:
          # appending all the colums, so easy
          linelist = numpy.append(linelist, elemllist)
        else:
          # Not keep the drv columns, so a little trickier
          linelisttmp = numpy.zeros(len(elemllist), dtype = linelist.dtype)
          for key in linelist.dtype.names:
            linelisttmp[key] = elemllist[key]
          linelist = numpy.append(linelist, linelisttmp)

    return linelist

def calc_ee_brems_spec(ebins, Te, dens, teunit='keV'):
  """
  Calculate the electron-electron bremstrahlung per-bin emissivity

  Parameters
  ----------
  ebins : array(float)
    The spectral bin edges in energy order (keV)
  Te : float
    Electron temperature (default, keV)
  dens : float
    The electron density (cm^-3)
  teunit : string
    Units of Te (keV by default, K also allowed)

  Returns
  -------
  eebrems : array(float)
    electron-electron bremsstrahlung in ph cm^3 s^-1 bin^-1

  """

  # convert temperature to keV if required
  kT = util.convert_temp(Te, teunit, 'keV')

  # get the emission at each bin edge
  eespec = apec.calc_ee_brems(ebins, kT, dens)

  # average over bin width & centroid
  #eecent = (eespec[1:]+eespec[:-1])/2
  #binwidth = ebins[1:]-ebins[:-1]
  ee = (ebins[1:]-ebins[:-1]) * (eespec[1:]+eespec[:-1])/2

  return ee


#### LEGACY CODE BEYOND THIS POINT

def __get_nei_line_emissivity(Z, z1, up, lo):
  """
  DEPRECATED
  Return the line emissivity for a single line, separated out by the ion driving it

  PARAMETERS
  ----------
  Z : int
    nuclear charge
  z1 : int
    ion charge +1
  up : int
    the upper level
  lo : int
    the lower level

  RETURNS
  -------
  emiss : dict
    a dictionary containing an array, one for each ion_drv, with the emissivity in it.
    E.g. emiss[6] is an nTe element array, with the emissivity due to z1=6 as a fn of temperature
    Also emiss['Te'] is the tempearture in keV
  """
  warnings.warn("get_nei_line_emissivity is a deprecated function and will be removed. Use NEISession.return_line_emissivity instead", DeprecationWarning, stacklevel=3)

  ldata= pyfits.open(os.path.expandvars("$ATOMDB/apec_nei_line.fits"))

  emiss={}
  datacache = {}
  emiss['Te'] = ldata[1].data['kT']
  for i in range(len(emiss['Te'])):
    j = ldata[i+2].data
    j2 = j[ (j['Element']==Z) &\
            (j['Ion']==z1) &\
            (j['UpperLev']==up) &\
            (j['LowerLev']==lo)]
    if len(j2)>0:
      ionbal = apec.return_ionbal(Z, emiss['Te'][i],  \
                       teunit='keV', \
                       datacache=datacache)


    for jj in j2:
      if not jj['ion_drv'] in emiss.keys():
        emiss[jj['Ion_drv']] = numpy.zeros(len(emiss['Te']))

      emiss[jj['Ion_drv']][i] = jj['Epsilon'] * ionbal[jj['ion_drv']-1]

  return emiss


