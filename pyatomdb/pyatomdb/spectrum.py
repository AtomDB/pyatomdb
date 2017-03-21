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

import numpy, os
# other pystomdb modules
import atomic, util, const, atomdb, apec

def make_spectrum(bins, index, linefile="$ATOMDB/apec_line.fits",\
                  cocofile="$ATOMDB/apec_coco.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  elements=False, abund=False, dummyfirst=False,\
                  dolines = True, docont=True, dopseudo=True):

  r"""
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


  # set up the bins
  if (sum((bins[1:]-bins[:-1])<0) > 0):
    print "*** ERROR: bins must be monotonically increasing. Exiting ***"
    return -1
  
  if binunits.lower()=='kev':
    ebins = bins*1.0
  elif binunits.lower() in ['a', 'angstrom', 'angstroms']:
    ebins = const.HC_IN_KEV_A/bins[::-1]
  else:
    print "*** ERROR: unknown binning unit %s, Must be keV or A. Exiting ***"%\
          (binunits)





  if util.keyword_check(linefile):
    # ok, we should do something with this
    # if it is a string, look for the file name
    if isinstance(linefile, basestring):
      lfile = os.path.expandvars(linefile)
      if not os.path.isfile(lfile):
        print "*** ERROR: no such file %s. Exiting ***" %(lfile)
        return -1
      ldat = pyfits.open(lfile)
    elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      ldat = linefile
    else:
      print "Unknown data type for linefile. Please pass a string or an HDUList"
      return -1
   
  if util.keyword_check(cocofile):
    if isinstance(cocofile, basestring):
      cfile = os.path.expandvars(cocofile)
      if not os.path.isfile(cfile):
        print "*** ERROR: no such file %s. Exiting ***" %(cfile)
        return -1
      cdat = pyfits.open(cfile)
    elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      cdat = cocofile
    else:
      print "Unknown data type for cocofile. Please pass a string or an HDUList"
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
  if ((index < 2) | (index > len(ldat))):
    print "*** ERRROR: Index must be in range %i to %i"%(2, len(ldat)-1)
    return -1
    
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
                                           binunits=binunits, no_coco=not(docont),\
                                           no_pseudo=not(dopseudo))*abund[iZ]
  
  # broaden the continuum if required:
  if broadening:
    cspectrum = broaden_continuum(ebins, cspectrum, binunits = binunits, \
                      broadening=broadening,\
                      broadenunits=broadenunits)
  if dummyfirst:
    return numpy.append([0],   cspectrum+lspectrum)
  else:
    return cspectrum+lspectrum



def make_ion_spectrum(bins, index, Z,z1, linefile="$ATOMDB/apec_nei_line.fits",\
                  cocofile="$ATOMDB/apec_nei_comp.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  abund=False, dummyfirst=False, nei = True,\
                  dolines = True, docont=True, dopseudo=True):

  r"""
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
    print "*** ERROR: bins must be monotonically increasing. Exiting ***"
    return -1
  
  if binunits.lower()=='kev':
    ebins = bins*1.0
  elif binunits.lower() in ['a', 'angstrom', 'angstroms']:
    ebins = const.HC_IN_KEV_A/bins[::-1]
  else:
    print "*** ERROR: unknown binning unit %s, Must be keV or A. Exiting ***"%\
          (binunits)

  # check the files exist
  # first, check if the line file is set
  if util.keyword_check(linefile):
    # ok, we should do something with this
    # if it is a string, look for the file name
    if isinstance(linefile, basestring):
      if ((linefile == "$ATOMDB/apec_nei_line.fits") & (nei==False)):
        linefile = "$ATOMDB/apec_line.fits"
      lfile = os.path.expandvars(linefile)
      if not os.path.isfile(lfile):
        print "*** ERROR: no such file %s. Exiting ***" %(lfile)
        return -1
      ldat = pyfits.open(lfile)
    elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      ldat = linefile
    else:
      print "Unknown data type for linefile. Please pass a string or an HDUList"
      return -1
   
  if util.keyword_check(cocofile):
    if isinstance(cocofile, basestring):
      if ((cocofile == "$ATOMDB/apec_nei_comp.fits") & (nei==False)):
        cocofile = "$ATOMDB/apec_coco.fits"
      cfile = os.path.expandvars(cocofile)
      if not os.path.isfile(cfile):
        print "*** ERROR: no such file %s. Exiting ***" %(cfile)
        return -1
      cdat = pyfits.open(cfile)
    elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
      # no need to do anything, file is already open
      cdat = cocofile
    else:
      print "Unknown data type for cocofile. Please pass a string or an HDUList"
      return
 
 
  
      
  # get the index
  if ((index < 2) | (index > len(ldat))):
    print "*** ERRROR: Index must be in range %i to %i"%(2, len(ldat)-1)
    return -1
    
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
    print "Error: unknown broadening unit %s, Must be keV or A. Exiting ***"%\
          (broadenunits)
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
    print "*** ERROR: unknown temeprature unit %s. Must be keV or K. Exiting ***"%\
          (teunits)
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
    print "*** ERROR: unknown unit %s, Must be keV or A. Exiting ***"%\
          (units)

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
      print "*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunit)

  # if the temperature is specified, get the index
  if index != False:

    if Te != False:
      print "Warning: both index and Te specified. Using index"
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
        print "Error: did not specify index or Te"
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
        print "*** ERROR. lldat provided as HDUList, but no index specified."
        print " Exiting"
      llist = numpy.array(lldat[index].data)
    elif type(lldat) == pyfits.hdu.table.BinTableHDU:
      llist = numpy.array(lldat.data)
    elif type(lldat) in [pyfits.fitsrec.FITS_rec, numpy.ndarray]:
      llist = numpy.array(lldat)
    else:
      print "ERROR: unkonwn llist type!"
  else:
    # there should now be no way to get here. commenting out this section
    print "ERROR: I SHOULD NEVER BE HERE"
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
    print "*** ERROR: unknown unit %s, Must be keV or A. Exiting ***"%\
          (units)

  # convert Te into keV

  if teunit.lower() == 'kev':
    kT = Te*1.0
  elif teunit.lower() == 'ev':
    kT = Te/1000.0
  elif teunit.lower() == 'k':
    kT = Te*const.KBOLTZ
  else:
    print "*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunit)


  if Te_init != False:
    if teunit.lower() == 'kev':
      kT_init = Te_init*1.0
    elif teunit.lower() == 'ev':
      kT_init = Te_init/1000.0
    elif teunit.lower() == 'k':
      kT_init = Te_init*const.KBOLTZ
    else:
      print "*** ERROR: unknown teunit %s, Must be keV or K. Exiting ***"%\
          (teunit)
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
      print "*** ERROR. Linefile %s is "%(linefile),
      print " not a file. Exiting"
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

def print_lines(llist, specunits = 'A'):
  """
  Prints lines in a linelist to screen
  
  This routine is very primitive as things stand. Plenty of room for refinement.
  
  Parameters
  ----------
  llist: dtype(linelist)
    list of lines to print. Typically returned by list_lines.
  specunits: {'A' , 'keV'}
    units to list the line positions by (A or keV, default A)
  
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
    print "*** ERROR: unknown unit %s, Must be keV or A. Exiting ***"%\
          (specunits)


  # now print the header lines
  
  if specunits == 'keV':
    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Energy','Epsilon','Element','Ion','Ion_drv','UpperLev','LowerLev')
    else:  
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Energy','Epsilon','Element','Ion','UpperLev','LowerLev')
    print s

    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('keV','ph cm3 s-1','','','','','')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('keV','ph cm3 s-1','','','','')
    print s
    
  else:
    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Lambda','Epsilon','Element','Ion','Ion_drv','UpperLev','LowerLev')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('Lambda','Epsilon','Element','Ion','UpperLev','LowerLev')
    print s

    if 'Ion_drv' in llist.dtype.names:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('A','ph cm3 s-1','','','','','')
    else:
      s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
         ('A','ph cm3 s-1','','','','')
    print s
    
  # now print the data
  
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
      print s
  else:
    for il in llist:
      s = "%.4e %.4e %-10i %-10i %-10i %-10i" %\
         (il['Lambda'],\
          il['Epsilon'],\
          il['Element'],\
          il['Ion'],\
          il['UpperLev'],\
          il['LowerLev'])
      print s
        
    
           



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
    print "*** ERROR: unknown units %s for continuum spectrum. Exiting" %\
          (binunits)

  if fluxunits.lower() in ['ph', 'photon','photons', 'p']:
    ergs = False
  elif fluxunits.lower() in ['erg','ergs']:
    ergs = true
  else:
    print "*** ERROR: unknown units %s for continuum flux. Exiting" %\
          (fluxunits)
          
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
    print "*** ERROR: unable to parse cocofile = %s" %repr(cocofile)
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
  for i in range(len(eedges)):
#    arg  = numpy.argwhere(iord==n+i)[0][0]
    C_out[i] = cum_cont[ihi[i]]
  
  
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
    print "*** ERROR: unknown units %s for continuum spectrum. Exiting" %\
          (binunits)

  if angstrom:
    bins = const.HC_IN_KEV_A/bins[::-1]
  
  # broadening
  if broadening:
    if broadenunits.lower() in ['a','angstrom','angstroms']:
      bunits = 'a'
    elif broadenunits.lower() in ['kev']:
      bunits = 'kev'
    else:
      print "*** ERROR: unknown units %s for continuum broadening. Exiting" %\
            (broadenunits)
      return -1
    # do the broadening
    spec = numpy.zeros(len(spectrum))
    emid = (bins[1:]+bins[:-1])/2
    if broadenunits == 'a':
      # convert to keV
      broadenvec = const.HC_IN_KEV_A/emid
    else:
      broadenvec = numpy.zeros(len(emid))
      broadenvec[:] = emid
    for i in range(len(spec)):
      
      spec += atomdb.addline2(bins, emid[i], \
                 spectrum[i],\
                 broadenvec[i])
    spectrum=spec
  if angstrom:
    spectrum=spectrum[::-1]
  return spectrum


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
      print "ERROR: unknown arf type, %s"%(repr(type(arf)))
      return
    res = spectrum * arfdat['SPECRESP'].data['SPECRESP']
  else:
    res = spectrum*1.0
  
  
  if type(rmf)==str:
    rmfdat = pyfits.open(rmf)
  elif type(rmf) == pyfits.hdu.hdulist.HDUList:
    rmfdat = rmf
  else:
    print "ERROR: unknown rmf type, %s"%(repr(type(rmf)))
    return
  
  ebins = rmfdat['EBOUNDS'].data['E_MIN']
  ebins = numpy.append(ebins, rmfdat['EBOUNDS'].data['E_MAX'][-1])

  ret = numpy.zeros(len(ebins)-1, dtype=float)

  try:
    k=rmfdat.index_of('MATRIX')
    matrixname = 'MATRIX'
  except KeyError:
    try:
      k=rmfdat.index_of('SPECRESP MATRIX')
      matrixname = 'SPECRESP MATRIX'
    except KeyError:
      print "Cannot find index for matrix in this data"
      raise
    
  for ibin, i in enumerate(rmfdat[matrixname].data):
    if res[ibin]==0.0: continue
    lobound = 0
    
    fchan = i['F_CHAN']
    nchan = i['N_CHAN']
    
    if numpy.isscalar(fchan):
      fchan = numpy.array([fchan])
      
    if numpy.isscalar(nchan):
      nchan = numpy.array([nchan])
    
    for j in range(len(fchan)):
      ilo = fchan[j]
      if ilo < 0: continue
      
      ihi = fchan[j] + nchan[j]
  
      ret[ilo:ihi] += res[ibin]*i['MATRIX'][lobound:lobound+nchan[j]]
      lobound += nchan[j]

#  print spectrum[:100]
#  print ret[:100]
  return ebins, ret


def get_response_ebins(rmf):
  """
  Get the energy bins from the rmf file
  
  Parameters
  ----------
  rmf : string or pyfits.hdu.hdulist.HDUList
    The filename of the rmf or the opened rmf file
  
  Returns
  -------
  array(float)
    input energy bins used. nbins+1 length, with the last item being the final bin
    This is the array on which the input spectrum should be calculated
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
    print "ERROR: unknown rmf type, %s"%(repr(type(rmf)))
    return
#  ret = rmfdat['EBOUNDS'].data['E_MIN']
#  ret = numpy.append(ret, rmfdat['EBOUNDS'].data['E_MAX'][-1])
  try:
    k=rmfdat.index_of('MATRIX')
    matrixname = 'MATRIX'
  except KeyError:
    try:
      k=rmfdat.index_of('SPECRESP MATRIX')
      matrixname = 'SPECRESP MATRIX'
    except KeyError:
      print "Cannot find index for matrix in this data"
      raise


  ret = rmfdat[matrixname].data['ENERG_LO']
  ret = numpy.append(ret, rmfdat[matrixname].data['ENERG_HI'][-1])

  return ret


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
      print "ERROR: unknown arf type, %s"%(repr(type(arf)))
      return
    arfarea = arfdat['SPECRESP'].data['SPECRESP']
  else:
    arfarea = 1.0
  
  
  if type(rmf)==str:
    rmfdat = pyfits.open(rmf)
  elif type(rmf) == pyfits.hdu.hdulist.HDUList:
    rmfdat = rmf
  else:
    print "ERROR: unknown rmf type, %s"%(repr(type(rmf)))
    return
  
  ebins = get_response_ebins(rmf)
  
  area = numpy.zeros(len(ebins)-1, dtype=float)
  
  for ibin, i in enumerate(rmfdat['MATRIX'].data):
    area[ibin] = sum(i['MATRIX'])

  area *= arfarea
  return ebins, area


class Session():
  """
  A session using the same line and coco files, and/or responses
  
  Attributes
  ----------
  linefile : string
    The line emissivity data file
  cocofile : string
    The continuum emissivity data file
  linedata: HDUList  
    The line emissivity data
  cocodata: HDUList  
    The line emissivity data
  elements : array_like, int
    The atomic number of the elements to include. Defaults to all.
  abundset : string
    The elemental abundances to be used. Defaults to Anders and 
    Grevesse 1989.
  ready : bool
    Set when line, continuum and spectral bin data has been
    read in, and a spectrum can be calculated.
  default_abundset : string
    The abundance set used in line and continuum files
  abundset : string
    The abundance set to be used in calculating the spectra.
  response_set : bool
    If a response (rmf & arf) have been loaded, set to true
  spectra : dict of array_like 
    Holds the spectra at each temperature.
  rmffile : string
    Filename of RMF file
  arffile : string
    Filename of ARF file
  rmf : HDUList
    RMF data
  arf : HDUList
    ARF data
    
  """

  def __init__(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits",\
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
    
    
    self.ready=False
    
    self.set_apec_files(linefile, cocofile)
    self.abundset=self.linedata[0].header['SABUND_SOURCE']
    self.default_abundset=self.linedata[0].header['SABUND_SOURCE']
    
    self.specbins_set=False

    # a hold for the spectra
    self.spectra={}

    # if elements are specified, use them. Otherwise, use Z=1-30
    if util.keyword_check(elements):
      self.elements = elements
    else:
      self.elements=range(1,31)
        
    # set the abundances:
    #   (1) the initial vector is whatever set AtomDB was calculated on,
    #       and is therefore 1.0
    self.abundsetvector = {}
    for Z in self.elements:
      self.abundsetvector[Z] = 1.0

    #   (2) but if another vector was already specified, use this instead
    if util.keyword_check(abundset):
      self.set_abundset(abundset)
    
    
    self.abund = {}
    for Z in self.elements:
      self.abund[Z]=1.0
    
    # Open the spectrum files
    #self.linefile = linefile
    #self.linefile = pyfits.open(self.linefile)
    #self.cocofile = cocofile
    #self.cocodata = pyfits.open(self.cocofile)
    
    
    self.response_set = False
    self.ready = True
    #if index >=2:
      #self.index_set = True
      #self.index = index
      #self.ready=True

    self.dolines=True
    self.docont=True
    self.dopseudo=True
    

  def return_spectra(self, te, teunit='keV', raw=False, nearest=False,\
                     get_nearest_t=False):
    """
    Get the spectrum at an exact temperature.
    Interpolates between 2 neighbouring spectra

    Finds HDU with kT closest ro desired kT in given line or coco file.

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

    if teunit.lower() == 'kev':
      teval = te
    elif teunit.lower() == 'ev':
      teval = te/1000.0
    elif teunit.lower() == 'k':
      teval = te*const.KBOLTZ
    else:
      print "*** ERROR: unknown temeprature unit %s. Must be eV, keV or K. Exiting ***"%\
            (teunits)
    
    if ((teval > self.linedata[1].data['kT'][-1]) |\
        (teval < self.linedata[1].data['kT'][0])):
      print "*** ERROR: temperature %f keV is out of range %f-%f ***" %\
            (teval, self.linedata[1].data['kT'][0], self.linedata[1].data['kT'][-1])
      return
    # find the 2 nearest temperatures
    
    if nearest:
      index = numpy.argmin(numpy.abs(self.linedata[1].data['kT']-teval))+2
      if not (index in self.spectra.keys()):
        self.spectra[index] = self.Spec(self, index)
        self.spectra[index].calc_spectrum(self)
      te_nearest = self.linedata[1].data['kT'][index-2]
      if teunit.lower()=='kev':
        pass
      elif teunit.lower() == 'ev':
        te_nearest /=1000
      elif teunit.lower() == 'k':
        te_nearest = te * const.KBOLTZ
        
      if not self.response_set:
        raw=True
    
      if raw:
        s = self.spectra[index].spectrum
      else:  
        s = self.spectra[index].spectrum_withresp
      if util.keyword_check(get_nearest_t):
        return s, te_nearest
      else:
        return s
    else:
      index = numpy.where(self.linedata[1].data['kT'] > teval)[0][0]


      loind = index+1
      upind = index+2
    
      # get the spectra at these temperatures
      if not (loind in self.spectra.keys()):
        self.spectra[loind] = self.Spec(self, loind)
        self.spectra[loind].calc_spectrum(self)
      if not (upind in self.spectra.keys()):
        self.spectra[upind] = self.Spec(self, upind)
        self.spectra[upind].calc_spectrum(self)
    
      # now sum the spectra and add as a response
    
      t1 = self.linedata[1].data['kT'][loind-2]
      t2 = self.linedata[1].data['kT'][upind-2]
    
    
      if not self.response_set:
        raw=True
    
      if raw:
        s1 = self.spectra[loind].spectrum
        s2 = self.spectra[upind].spectrum
    
      else:  
        s1 = self.spectra[loind].spectrum_withresp
        s2 = self.spectra[upind].spectrum_withresp
    
      # linear interp
    
      r1 = 1- (teval-t1)/(t2-t1)
      r2 = 1- r1
    
      s = s1*r1 + s2*r2
    
    return s
    
    
    

  def set_specbins(self, specbins, specunits='A'):
    """
    Set the energy or wavelength bin for the raw spectrum
    
    Note that this is overridden if a response is loaded
    
    Parameters
    ----------
    ebins : array(float)
      The edges of the spectral bins (for n bins, have n+1 edges)
    specunits : {'a','kev'}
      The spectral bin units to use. Default is angstroms
      
    Returns
    -------
    None
    
    Notes
    -----
    updates  self.specbins, self.binunits, self.specbins_set
    """
    
    # set the energy bins for this spectrum
    self.specbins=specbins
    self.binunits=specunits
    self.specbins_set=True
    
    
  def set_response(self, rmf, arf=False):
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
    
    Parameters
    ----------
    rmf: string or HDUlist
      The response matrix file
    arf: string or HDUlist
      The ancillary response file
    
    Returns
    -------
    none

    """
    
    if util.keyword_check(arf):
      if type(arf)==str:
        self.arffile = arf
        self.arf = pyfits.open(arf)
      elif type(arf) == pyfits.hdu.hdulist.HDUList:
        self.arf = arf
        self.arffile = arf.filename()
      else:
        print "ERROR: unknown arf type, %s"%(repr(type(arf)))
        return
#      res   = spectrum * arfdat['SPECRESP'].data['SPECRESP']
    #else:
      #res = spectrum*1.0
  
  
    if type(rmf)==str:
      self.rmffile = rmf
      self.rmf = pyfits.open(rmf)
    elif type(rmf) == pyfits.hdu.hdulist.HDUList:
      self.rmf = rmf
      self.rmffile = rmf.filename()
    else:
      print "ERROR: unknown rmf type, %s"%(repr(type(rmf)))
      return
    
    self.ebins_response = get_response_ebins(self.rmf)
    self.response_set = True

  def set_apec_files(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits"):
    """
    Set the apec line and coco files
    
    Parameters
    ----------
    linefile : str or HDUList
      The filename of the line emissivity data, or the opened file.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the opened file.
    elements : array_like(int)
      The atomic numbers of the elements to include. Defaults to all (1-30)
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance 
    
    Returns
    -------
    None
    
    Notes
    -----
    Updates self.linefile, self.linedata, self.cocofile and self.cocodata
    """
    if util.keyword_check(linefile):
      if isinstance(linefile, basestring):
        lfile = os.path.expandvars(linefile)
        if not os.path.isfile(lfile):
          print "*** ERROR: no such file %s. Exiting ***" %(lfile)
          return -1
        self.linedata = pyfits.open(lfile)
        self.linefile = lfile

      elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.linedata=linefile
        self.linefile=linefile.filename()

      else:
        print "Unknown data type for linefile. Please pass a string or an HDUList"

    if util.keyword_check(cocofile):

      if isinstance(cocofile, basestring):

        cfile = os.path.expandvars(cocofile)
        if not os.path.isfile(cfile):
          print "*** ERROR: no such file %s. Exiting ***" %(cfile)
          return -1
        self.cocodata=pyfits.open(cfile)
        self.cocofile=cfile

      elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.cocodata=cocofile
        self.cocofile=cocofile.filename()

      else:
        print "Unknown data type for cocofile. Please pass a string or an HDUList"
  
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

        print "abundance vector and element vector must have same number"+\
              " of elements"
      else:
        
        self.abund[elementvec] = abundvec
    elif (eisvec):
      # set all these eleemtns to the same abundance
      for el in elementvec:
        self.abund[el]=abund  
      
    else:
      self.abund[elements]=abund  
              
    self.recalc()  
    

  def recalc(self):
    """
    Recalculate the spectrum - just for changing abundances etc. 
    Does not recalculate spectrum fully, just changes the multipliers.
    Does nothing if self.ready is False, should be run after calc_spectrum.
    
    Parameters
    ----------
    none
    
    Returns
    -------
    none
    
    Notes
    --------
    modifies
    self.spectrum
    """
    for index in self.spectra.keys():
      self.spectra[index].recalc(self)
    #if self.ready:
      #self.spectrum = numpy.zeros(len(self.specbins)-1)
      #for Z in self.elements:
        #self.spectrum += self.spectrum_by_Z[Z] * self.abund[Z] * self.abundsetvector[Z]
      #if self.response_set:
        #self.spectrum_withresp = numpy.zeros(len(self.ebins_response)-1)
        #for Z in self.elements:
          #self.spectrum_withresp += self.spectrum_by_Z_withresp[Z] * self.abund[Z] * self.abundsetvector[Z]
      
      
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
    old = atomdb.get_abundance(abundset=self.default_abundset)
    
    # read in the new abundance
    new = atomdb.get_abundance(abundset=abundstring)
  
    # divide the 2, store the replacement ratio to self.abundsetvector
    for Z in self.abundsetvector.keys():
      self.abundsetvector[Z]=new[Z]/old[Z]

    # update the current abundance string to represent your input
    self.abundset=abundstring
    
    self.recalc()
  
  
  class Spec():
      """
      An individual spectrum within a session, from a specifically
      tabulated temperature in a line/coco file.
    
      Attributes
      ----------
      temperature : float
        The temperature of this spectrum, in keV
      index : int
        The index in the line file for this spectrum 
      """
    
    
      def __init__(self, session, index):
        self.temperature = session.linedata[1].data['kT'][index-2]
        self.index = index
      
      
        
      def set_index(T, teunit='K', logscale = False):
        """
        Finds HDU with kT closest to desired kT in given line or coco file.
    
        Opens the line or coco file, and looks for the header unit
        with temperature closest to te. Use result as index input to make_spectrum
    
        Parameters
        ----------
        te : float
          Temperature in keV or K
        teunits : {'keV' , 'K', 'eV'}
          Units of te (kev or K, default keV)
        logscale : bool
          Search on a log scale for nearest temperature if set.
      
        Returns
        -------
        none
      
        Notes
        -----
        modifies
        self.index : int
        Index in HDU file with nearest temperature to te.
            
        """  
    
        if teunit.lower() == 'kev':
          teval = te
        elif teunit.lower() == 'ev':
          teval = te/1000.0
        elif teunit.lower() == 'k':
          teval = te*const.KBOLTZ
        else:
          print "*** ERROR: unknown temeprature unit %s. Must be eV, keV or K. Exiting ***"%\
                (teunits)
      
        if logscale:
          i = numpy.argmin(numpy.abs(numpy.log(self.linedata[1].data['kT'])-numpy.log(teval)))
        else:
          i = numpy.argmin(numpy.abs(self.linedata[1].data['kT']-teval))
        # need to increase the HDU by 2.
        self.index = i+2
      
      def calc_spectrum(self,session,
                        dolines = True, docont=True, dopseudo=True):
      
        """
        Calculates the spectrum for each element on a single temperature
      
        Parameters
        ----------
        session : Session
          The parent Session
        dolines : bool
          Include lines in the spectrum
        docont : bool
          Include continuum in the spectrum
        dopseudo : bool
          Include pseudocontinuum in the spectrum
        Outputs
        -------
        none
      
        Notes
        -----
        Modifies:\n
        dict : self.spectrum_by_Z  the spectrum of each element\n
        dict : self.spectrum_by_Z_withresp  the spectrum of each element, \
                                           folded through response\n
        Then calls `recalc()` to update the spectra
        """
        # now, we shall calculate the spectrum for each individual element
    
        # set the linefile
    
        self.spectrum_by_Z={}
        if session.response_set==True:
          self.spectrum_by_Z_withresp={}
    
        
        
        self.temperature = session.linedata[1].data['kT'][self.index-2]
      
    #    if util.keyword_check(elements):
          #self.set_abund(elements, abund)
      
        #if util.keyword_check(abund):
          #self.set_abund(elements, abund)
      
        for Z in session.elements:
        # make the generic spectrum
          if session.specbins_set:
            self.spectrum_by_Z[Z] = make_spectrum(session.specbins, self.index,\
                                                  session.linedata, session.cocodata,\
                                                  session.binunits,elements=[Z],\
                                                  dolines=dolines,\
                                                  docont = docont,\
                                                  dopseudo = dopseudo)
        # make the spectrum on the response grid
          if session.response_set:
            tmp = make_spectrum(session.ebins_response, self.index,\
                                session.linedata, session.cocodata,\
                                'keV',elements=[Z],\
                                dolines=dolines,\
                                docont = docont,\
                                dopseudo = dopseudo)
                                
            e,self.spectrum_by_Z_withresp[Z] = apply_response(tmp, session.rmf, arf=session.arf)
    
        self.recalc(session)
    
    
      def recalc(self, session):
        """
        Recalculate the spectrum - just for changing abundances etc. 
        Does not recalculate spectrum fully, just changes the multipliers.
        Does nothing if self.ready is False, should be run after calc_spectrum.
        
        Parameters
        ----------
        session : Session
          The parent session
        
        Returns
        -------
        none
        
        Notes
        -----
        Modifies:\n
        self.spectrum : array_like (float)\n
        self.spectrum_withresp : array_like (float)
        """
    
        if session.ready:
          if session.specbins_set:
            
            self.spectrum = numpy.zeros(len(session.specbins)-1)
            for Z in session.elements:
              self.spectrum += self.spectrum_by_Z[Z] * session.abund[Z] * session.abundsetvector[Z]
          if session.response_set:
            self.spectrum_withresp = numpy.zeros(len(session.ebins_response)-1)
            for Z in session.elements:
              self.spectrum_withresp += self.spectrum_by_Z_withresp[Z] * session.abund[Z] * session.abundsetvector[Z]
    
  


class CXSession(Session):
  """
  A Charge Exchange session using the same line and coco files, and/or responses
  
  Attributes
  ----------
  linefile : string
    The line emissivity data file
  cocofile : string
    The continuum emissivity data file
  linedata: HDUList  
    The line emissivity data
  cocodata: HDUList  
    The line emissivity data
  elements : array_like, int
    The atomic number of the elements to include. Defaults to all.
  abundset : string
    The elemental abundances to be used. Defaults to Anders and 
    Grevesse 1989.
  collisionunits : string
    Whether the units are given in energy (kev/amu) or velocity (cm/s)
  ready : bool
    Set when line, continuum and spectral bin data has been
    read in, and a spectrum can be calculated.
  default_abundset : string
    The abundance set used in line and continuum files
  abundset : string
    The abundance set to be used in calculating the spectra.
  response_set : bool
    If a response (rmf & arf) have been loaded, set to true
  spectra : dict of array_like 
    Holds the spectra at each temperature.
  rmffile : string
    Filename of RMF file
  arffile : string
    Filename of ARF file
  rmf : HDUList
    RMF data
  arf : HDUList
    ARF data
  ionbal : dict of array like
    ionization balance of each ion, normalized to 1
    e.g. ionbal[6]=numpy.array([0.5,0.4,0.1,0,0,0,0]) for Carbon
    
    
  """

  def __init__(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits",\
                     elements=False, abundset='AG89',
                     collisionunits = 'kev/amu'):
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
    collisionunits : string
      The units for the particle collision speed. Either 'kev/amu' or 'cm/s'.
    """
    
    
    self.ready=False
    
    self.set_apec_files(linefile, cocofile)
    self.abundset=self.linedata[0].header['SABUND_SOURCE']
    self.default_abundset=self.linedata[0].header['SABUND_SOURCE']
    
    self.specbins_set=False

    # a hold for the spectra
    self.spectra={}

    # if elements are specified, use them. Otherwise, use Z=1-30
    if util.keyword_check(elements):
      self.elements = elements
    else:
      self.elements=range(1,31)
        
    # set the abundances:
    #   (1) the initial vector is whatever set AtomDB was calculated on,
    #       and is therefore 1.0
    self.abundsetvector = {}
    for Z in self.elements:
      self.abundsetvector[Z] = 1.0

    #   (2) but if another vector was already specified, use this instead
    if util.keyword_check(abundset):
      self.set_abundset(abundset)
    
    
    self.abund = {}
    for Z in self.elements:
      self.abund[Z]=1.0
    

    self.collisionunits = collisionunits
    # Open the spectrum files
    #self.linefile = linefile
    #self.linefile = pyfits.open(self.linefile)
    #self.cocofile = cocofile
    #self.cocodata = pyfits.open(self.cocofile)
    
    
    self.response_set = False
    self.ready = True
    #if index >=2:
      #self.index_set = True
      #self.index = index
      #self.ready=True

    self.ionbal_set=False
    
  def set_ionbal_temperature(self, te, teunit='keV'):
    """
    Set the ionization balance to that of a given electron temperature

    Parameters
    ----------
    te : float
      Electron Temperature
    teunit : string
      Units for the temperature. keV or A.

    Return
    ------
    none

    Notes
    -----
    Modifies self.ionbal
    """
    ionbal = {}

    for Z in self.elements:
      ionbal[Z] = apec.solve_ionbal_eigen(Z, te, \
                       teunit=teunit)
    self.ionbal = ionbal
    self.ionbal_set = True
    
    
  def return_spectra(self, collision, raw=False, nearest=False):
    """
    Get the spectrum at an exact temperature.
    Interpolates between 2 neighbouring spectra

    Finds HDU with kT closest ro desired kT in given line or coco file.

    Opens the line or coco file, and looks for the header unit
    with temperature closest to te. Use result as index input to make_spectrum

    Parameters
    ----------
    collision : float
      The energy (kev/amu) or velocity (cm/s) of the collision
    raw : bool
      If set, return the spectrum without response applied. Default False.

    Returns
    -------
    spectrum : array(float)
      The spectrum in photons cm^5 s^-1 bin^-1, with the response, or
      photons cm^3 s^-1 bin^-1 if raw is set.
    """  

    # for each ion, calculate the velocity
    
    
    stot = numpy.zeros(len(self.specbins)-1, dtype=float)
    for Z in self.elements:
      # get convert the velocity, etc, into keV/amu
      if self.collisionunits.lower() == 'cm/s':
        # convert from velocity to energy

        # in joules
        print " requesting velocity = %e cm/s = %e m/s"%(collision, collision/100)
        v_ms = collision/100
        mass_amu=atomic.Z_to_mass(Z)+1
        mass_kg = (atomic.Z_to_mass(Z)+1)* const.AMUKG
        print "mass = %e amu %e kg"%(mass_amu,mass_kg)

        E_J = 0.5*mass_kg * (v_ms)**2
        print "energy = %e J"%(E_J)

        E_keV = E_J/1.602e-16
        print "energy_keV = %e keV"%(E_keV)
        
#e=0.5 * (atomic.Z_to_mass(Z)+1) * const.AMUKG * (collision/100)**2

        e = E_keV /mass_amu
        print "e = %e keV/amu"%(e)

        
      else:
        e = collision*1.0
        
      

      if ((e > self.linedata[1].data['energy'][-1]) |\
          (e < self.linedata[1].data['energy'][0])):
        print "*** WARNING: energy %e keV is out of range %f-%f ***" %\
            (e, self.linedata[1].data['energy'][0], self.linedata[1].data['energy'][-1])
        continue

      # OK, so we have a valid energy
      # Make a spectrum for each ion

      index = numpy.where(self.linedata[1].data['energy'] > e)[0][0]


      loind = index+1
      upind = index+2
      E1 = self.linedata[1].data['Energy'][loind-2]
      E2 = self.linedata[1].data['Energy'][upind-2]
    
      for z1 in range(1,Z+2):
        if self.ionbal[Z][z1-1]>1e-10:
          # we should do things
          for ind in [loind, upind]:
            if not ind in self.spectra:
              self.spectra[ind]={}
            if not Z in self.spectra[ind].keys():
              self.spectra[ind][Z]={}
            if not z1 in self.spectra[ind][Z].keys():
              self.spectra[ind][Z][z1] = self.IonSpec(self, ind, Z, z1)
              self.spectra[ind][Z][z1].calc_spectrum(self)

            
          # OK, so I have made sure that all the spectra are generated
        
          # now evaluate the spectra at those energies.
      
      
          if not self.response_set:
            raw=True
        
          if raw:
            print Z, z1, self.ionbal[Z][z1-1]
            s1 = self.spectra[loind][Z][z1].spectrum * self.ionbal[Z][z1-1]
            s2 = self.spectra[upind][Z][z1].spectrum * self.ionbal[Z][z1-1]
        
          else:  
            s1 = self.spectra[loind][Z][z1].spectrum_withresp * self.ionbal[Z][z1]
            s2 = self.spectra[upind][Z][z1].spectrum_withresp * self.ionbal[Z][z1]
        
          # linear interp
        
          r1 = 1- (e-E1)/(E2-E1)
          r2 = 1- r1
        
          s = s1*r1 + s2*r2
          stot+=s
     
    return stot
  
    
    

  def set_specbins(self, specbins, specunits='A'):
    """
    Set the energy or wavelength bin for the raw spectrum
    
    Note that this is overridden if a response is loaded
    
    Parameters
    ----------
    ebins : array(float)
      The edges of the spectral bins (for n bins, have n+1 edges)
    specunits : {'a','kev'}
      The spectral bin units to use. Default is angstroms
      
    Returns
    -------
    None
    
    Notes
    -----
    updates  self.specbins, self.binunits, self.specbins_set
    """
    
    # set the energy bins for this spectrum
    self.specbins=specbins
    self.binunits=specunits
    self.specbins_set=True
    
    
  def set_response(self, rmf, arf=False):
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
    
    Parameters
    ----------
    rmf: string or HDUlist
      The response matrix file
    arf: string or HDUlist
      The ancillary response file
    
    Returns
    -------
    none

    """
    
    if util.keyword_check(arf):
      if type(arf)==str:
        self.arffile = arf
        self.arf = pyfits.open(arf)
      elif type(arf) == pyfits.hdu.hdulist.HDUList:
        self.arf = arf
        self.arffile = arf.filename()
      else:
        print "ERROR: unknown arf type, %s"%(repr(type(arf)))
        return
#      res   = spectrum * arfdat['SPECRESP'].data['SPECRESP']
    #else:
      #res = spectrum*1.0
  
  
    if type(rmf)==str:
      self.rmffile = rmf
      self.rmf = pyfits.open(rmf)
    elif type(rmf) == pyfits.hdu.hdulist.HDUList:
      self.rmf = rmf
      self.rmffile = rmf.filename()
    else:
      print "ERROR: unknown rmf type, %s"%(repr(type(rmf)))
      return
    
    self.ebins_response = get_response_ebins(self.rmf)
    self.response_set = True

  def set_apec_files(self, linefile="$ATOMDB/apec_line.fits",\
                     cocofile="$ATOMDB/apec_coco.fits"):
    """
    Set the apec line and coco files
    
    Parameters
    ----------
    linefile : str or HDUList
      The filename of the line emissivity data, or the opened file.
    cocofile : str or HDUList
      The filename of the continuum emissivity data, or the opened file.
    elements : array_like(int)
      The atomic numbers of the elements to include. Defaults to all (1-30)
    abundset : string
      The abundance set to use. Defaults to AG89. See atomdb.set_abundance 
    
    Returns
    -------
    None
    
    Notes
    -----
    Updates self.linefile, self.linedata, self.cocofile and self.cocodata
    """
    if util.keyword_check(linefile):
      if isinstance(linefile, basestring):
        lfile = os.path.expandvars(linefile)
        if not os.path.isfile(lfile):
          print "*** ERROR: no such file %s. Exiting ***" %(lfile)
          return -1
        self.linedata = pyfits.open(lfile)
        self.linefile = lfile

      elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.linedata=linefile
        self.linefile=linefile.filename()

      else:
        print "Unknown data type for linefile. Please pass a string or an HDUList"

    if util.keyword_check(cocofile):

      if isinstance(cocofile, basestring):

        cfile = os.path.expandvars(cocofile)
        if not os.path.isfile(cfile):
          print "*** ERROR: no such file %s. Exiting ***" %(cfile)
          return -1
        self.cocodata=pyfits.open(cfile)
        self.cocofile=cfile

      elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
        # no need to do anything, file is already open
        self.cocodata=cocofile
        self.cocofile=cocofile.filename()

      else:
        print "Unknown data type for cocofile. Please pass a string or an HDUList"
  
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

        print "abundance vector and element vector must have same number"+\
              " of elements"
      else:
        
        self.abund[elementvec] = abundvec
    elif (eisvec):
      # set all these eleemtns to the same abundance
      for el in elementvec:
        self.abund[el]=abund  
      
    else:
      self.abund[elements]=abund  
              
    self.recalc()  
    

  def recalc(self):
    """
    Recalculate the spectrum - just for changing abundances etc. 
    Does not recalculate spectrum fully, just changes the multipliers.
    Does nothing if self.ready is False, should be run after calc_spectrum.
    
    Parameters
    ----------
    none
    
    Returns
    -------
    none
    
    Notes
    --------
    modifies
    self.spectrum
    """
    for index in self.spectra.keys():
      self.spectra[index].recalc(self)
    #if self.ready:
      #self.spectrum = numpy.zeros(len(self.specbins)-1)
      #for Z in self.elements:
        #self.spectrum += self.spectrum_by_Z[Z] * self.abund[Z] * self.abundsetvector[Z]
      #if self.response_set:
        #self.spectrum_withresp = numpy.zeros(len(self.ebins_response)-1)
        #for Z in self.elements:
          #self.spectrum_withresp += self.spectrum_by_Z_withresp[Z] * self.abund[Z] * self.abundsetvector[Z]
      
      
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
    old = atomdb.get_abundance(abundset=self.default_abundset)
    
    # read in the new abundance
    new = atomdb.get_abundance(abundset=abundstring)
  
    # divide the 2, store the replacement ratio to self.abundsetvector
    for Z in self.abundsetvector.keys():
      self.abundsetvector[Z]=new[Z]/old[Z]

    # update the current abundance string to represent your input
    self.abundset=abundstring
    
    self.recalc()
  
  
  class IonSpec():
      """
      An individual ion spectrum within a session, from a specifically
      tabulated temperature in a line/coco file.
    
      Attributes
      ----------
      energy : float
        The energy of this spectrum, in keV/amu
      index : int
        The index in the line file for this spectrum
      """
    
    
      def __init__(self, session, index, Z, z1):
        
        self.energy = session.linedata[1].data['Energy'][index-2]
        self.index = index
        self.Z = Z
        self.z1 = z1
      
      
        
      def set_index(E, Eunit='kev/amu', logscale = False):
        """
        Finds HDU with kT closest to desired kT in given line or coco file.
    
        Opens the line or coco file, and looks for the header unit
        with temperature closest to te. Use result as index input to make_spectrum
    
        Parameters
        ----------
        E : float
          Energy in keV/amu
        Eunit : {'keV' , 'K', 'eV'}
          Units of E (kev/amu default)
        logscale : bool
          Search on a log scale for nearest temperature if set.
      
        Returns
        -------
        none
      
        Notes
        -----
        modifies
        self.index : int
        Index in HDU file with nearest temperature to te.
            
        """  
    
        if Eunit.lower() == 'kev/amu':
          Eval = E*1.0
          
#        elif teunit.lower() == 'ev':
#          teval = te/1000.0
#        elif teunit.lower() == 'k':
#          teval = te*const.KBOLTZ
        else:
          print "*** ERROR: unknown energy unit %s. Must be keV/amu. Exiting ***"%\
                (Eunit)
      
        if logscale:
          i = numpy.argmin(numpy.abs(numpy.log(self.linedata[1].data['Energy'])-numpy.log(teval)))
        else:
          i = numpy.argmin(numpy.abs(self.linedata[1].data['Energy']-teval))
        # need to increase the HDU by 2.
        self.index = i+2
      
      def calc_spectrum(self,session,
                        dolines = True, docont=True, dopseudo=True):
      
        """
        Calculates the spectrum for each ion on a single energy
      
        Parameters
        ----------
        session : Session
          The parent Session
        dolines : bool
          Include lines in the spectrum
        docont : bool
          Include continuum in the spectrum
        dopseudo : bool
          Include pseudocontinuum in the spectrum
        Outputs
        -------
        none
      
        Notes
        -----
        Modifies:\n
        dict : self.spectrum  the spectrum of the ion\n
        dict : self.spectrum_withresp the spectrum of the ion, \
                                           folded through response\n
        Then calls `recalc()` to update the spectra
        """
        # now, we shall calculate the spectrum for each individual element
    
        # set the linefile
    
        
        
        
        self.energy = session.linedata[1].data['Energy'][self.index-2]
      
    #    if util.keyword_check(elements):
          #self.set_abund(elements, abund)
      
        #if util.keyword_check(abund):
          #self.set_abund(elements, abund)
      
        # make the generic spectrum
        if session.specbins_set:
          self.spectrum = make_ion_spectrum(session.specbins, self.index,\
                              self.Z, self.z1, \
                              linefile = session.linedata, \
                              cocofile = session.cocodata,\
                              binunits = 'keV',\
                              dolines=dolines,\
                              docont = docont,\
                              dopseudo = dopseudo)
        # make the spectrum on the response grid
        if session.response_set:

          tmp = make_ion_spectrum(session.ebins_response, self.index,\
                              self.Z, self.z1, \
                              linefile = session.linedata, \
                              cocofile = session.cocodata,\
                              binunits = 'keV',\
                              dolines=dolines,\
                              docont = docont,\
                              dopseudo = dopseudo)
                                
          e,self.spectrum = apply_response(tmp, session.rmf, arf=session.arf)
    
        self.recalc(session)
    
    
      def recalc(self, session):
        """
        Recalculate the spectrum - just for changing abundances etc. 
        Does not recalculate spectrum fully, just changes the multipliers.
        Does nothing if self.ready is False, should be run after calc_spectrum.
        
        Parameters
        ----------
        session : Session
          The parent session
        
        Returns
        -------
        none
        
        Notes
        -----
        Modifies:\n
        self.spectrum : array_like (float)\n
        self.spectrum_withresp : array_like (float)
        """
    
        if session.ready:
          if session.specbins_set:
            
            self.spectrum = self.spectrum * session.abund[self.Z] * session.abundsetvector[self.Z]
          if session.response_set:
            self.spectrum_withresp = self.spectrum_withresp * session.abund[self.Z] * session.abundsetvector[self.Z]
    
  
