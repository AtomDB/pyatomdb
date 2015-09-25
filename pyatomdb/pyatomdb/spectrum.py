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
import atomic, util, const, atomdb

def make_spectrum(bins, index, linefile="$ATOMDB/apec_line.fits",\
                  cocofile="$ATOMDB/apec_coco.fits",\
                  binunits='keV', broadening=False, broadenunits='keV', \
                  elements=False, abund=False, dummyfirst=False):

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
  
  Returns
  -------
  array of floats
      eEissivity in counts cm^3 s^-1 bin^-1.
  
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
  lfile = os.path.expandvars(linefile)
  cfile = os.path.expandvars(cocofile)
  if not os.path.isfile(lfile):
    print "*** ERROR: no such file %s. Exiting ***" %(lfile)
    return -1
  if not os.path.isfile(cfile):
    print "*** ERROR: no such file %s. Exiting ***" %(cfile)
    return -1
  
  # open the files
  ldat = pyfits.open(lfile)
  cdat = pyfits.open(cfile)
      
  # get the index
  if ((index < 2) | (index > len(ldat))):
    print "*** ERRROR: Index must be in range %i to %i"%(2, len(ldat)-1)
    return -1
    
  lldat = ldat[index].data
  ccdat = cdat[index].data
  
  if not elements:
    Zl = util.unique(lldat['element'])
    Zc = util.unique(ccdat['Z'])
    Zlist = util.unique(numpy.append(Zl,Zc))
  
  else:
    Zlist = elements
  
  if not abund:
    abund= numpy.ones(len(Zlist))

  lspectrum = numpy.zeros(len(bins)-1, dtype=float)
  cspectrum = numpy.zeros(len(bins)-1, dtype=float)
  
  for iZ, Z in enumerate(Zlist):
    # ADD  LINES
    lspectrum += add_lines(Z, abund[iZ], lldat, ebins, broadening, broadenunits)
    
  for iZ, Z in enumerate(Zlist):
    # ADD  CONTINUUM
    cspectrum += make_ion_index_continuum(ebins, Z, cocofile=ccdat,\
                                         binunits=binunits)*abund[iZ]
  
  # broaden the continuum if required:
  if broadening:
    cspectrum = broaden_continuum(ebins, cspectrum, binunits = binunits, \
                      broadening=broadening,\
                      broadenunits=broadenunits)
  if dummyfirst:
    return numpy.append([0],   cspectrum+lspectrum)
  else:
    return cspectrum+lspectrum
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def add_lines(Z, abund, lldat, ebins, broadening=False, broadenunits='A'):
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
               numpy.where(ebins < const.HC_IN_KEV_A/ll['lambda'])[0])]+=\
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
    i = numpy.argmin(numpy.abs(numpy.log(a['kT'])-numpy.log(te)))
  else:
    i = numpy.argmin(numpy.abs(a['kT']-te))
  # need to increase the HDU by 2.
  return i+2
  
  
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def list_lines(specrange, lldat=False, index=False, linefile=False,\
              units='angstroms'):
  """
  Gets list of the lines in a given spectral range
  
  Note that the output from this can be passed directly to print_lines
  
  
  Parameters
  ----------
  specrange : [float,float]
    spectral range [min,max] to return lines on
  units : {'A' , 'keV'}
    units of specrange (default A)
  lldat : see notes
    line data
  index : int
    index in lldat, see notes
  linefile : see notes
    line data file, see notes

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


  # check that either "lldat" is provided, which is a line list, or
  # that an index in a line file is specified
    
    
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
    # no line data supplied.
    linefile = os.path.expandvars('$ATOMDB/apec_line.fits')
    if not os.path.isfile(linefile):
      print "*** ERROR. No linefile supplied but $ATOMDB/apec_line.fits is"
      print " not a file. Exiting"
    else:

      if not index:
        print "*** ERROR. No index specified for line list. Exiting"
      else:
        lfile = pyfits.open(os.path.expandvars(linefile))
        llist= numpy.array(lfile[index].data)

  # at this point, we have data
  llist = llist[(llist['Lambda']>= specrange[0]) &\
                (llist['Lambda']<= specrange[1])]

               
          
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
    s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
       ('Energy','Epsilon','Element','Ion','UpperLev','LowerLev')
    print s
    s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
       ('keV','ph cm3 s-1','','','','')
    print s
    
  else:
    s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
       ('Lambda','Epsilon','Element','Ion','UpperLev','LowerLev')
    print s
    s= "%-10s %-10s %-10s %-10s %-10s %-10s" %\
       ('A','ph cm3 s-1','','','','')
    print s
    
  # now print the data

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
    cdat = pyfits.open(os.expandvars(cocofile))[index].data
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
    bins = HC_IN_KEV_A/bins[::-1]
  
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
