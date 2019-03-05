import pyatomdb, numpy, os

"""
This code is an example of generating a cooling curve: the total power
radiated in keV cm3 s-1 by each element at each temperature. It will
generate a text file with the emission per element at each temperature
from 1e4 to 1e9K.

This is similar to the atomdb.lorentz_power function, but with a few
normalizations removed to run a little quicker.

Note that the Anders and Grevesse (1989) abundances are built in to
this. These can be looked up using atomdb.get_abundance(abundset='AG89'),
or the 'angr' column of the table at
https://heasarc.nasa.gov/xanadu/xspec/xspec11/manual/node33.html#SECTION00631000000000000000

Adjustable parameters (energy range, element choice) are in the block
marked #### ADJUST THINGS HERE

Usage: python3 calc_power.py
"""



def calc_power_oneelem_oneT(Z, Elo, Ehi, ihdu,linefile="$ATOMDB/apec_line.fits",\
                       cocofile="$ATOMDB/apec_coco.fits"):
  """
  Calculate the radiated power between 2 different energies

  INPUTS
  ------
  Z : int
    The element atomic number
  Elo : float
    The lower energy bound
  Ehi : float
    The upper energy bound
  ihdu : int
    The HDU to use, from 2 to 52. Each is a different temperature from
    10^4 to 10^9K in log space.
  linefile : string/HDUlist
    If a string, filename of the line emission. If HDU list, that file, already open
  cocofile : string/HDUlist
    If a string, filename of the continuum emission. If HDU list, that file, already open

  """

  # make energy bins
  ebins = numpy.linspace(Elo, Ehi, 10000)
  en = (ebins[1:]+ebins[:-1])/2

  #
  spec = pyatomdb.spectrum.make_spectrum(ebins, ihdu, linefile=linefile,\
                       cocofile=cocofile,\
                       elements=[Z])

  # now you have a spectrum in photons. Convert to keV

  E = spec*en # energy in keV cm^3 s^-1

  return sum(E)

def calc_power_oneT(Zlist, Elo, Ehi, ihdu, linefile="$ATOMDB/apec_line.fits",\
                       cocofile="$ATOMDB/apec_coco.fits"):
  """
  Zlist : [int]
    List of element nuclear charges
  Elo : float
    The lower energy bound
  Ehi : float
    The upper energy bound
  ihdu : int
    The HDU to use, from 2 to 52. Each is a different temperature from
    10^4 to 10^9K in log space.
  linefile : string/HDUlist
    If a string, filename of the line emission. If HDU list, that file, already open
  cocofile : string/HDUlist
    If a string, filename of the continuum emission. If HDU list, that file, already open
  """
  E={}
  for Z in Zlist:
    E[Z] = calc_power_oneelem_oneT(Z, Elo, Ehi, ihdu, linefile = linefile, cocofile = cocofile)
  return E

def calc_power(Zlist, Elo, Ehi, hdulist=range(2,53), linefile="$ATOMDB/apec_line.fits",\
                       cocofile="$ATOMDB/apec_coco.fits"):

  """
  Zlist : [int]
    List of element nuclear charges
  Elo : float
    The lower energy bound
  Ehi : float
    The upper energy bound
  hdulist : [int]
    The HDUs to calculate the emission on, from 2 to 52. Each is a different temperature from
    10^4 to 10^9K in log space. If not given, will do for all 51 temperatures.
  linefile : string/HDUlist
    If a string, filename of the line emission. If HDU list, that file, already open
  cocofile : string/HDUlist
    If a string, filename of the continuum emission. If HDU list, that file, already open
  """
  res = {}
  res['power'] = {}
  res['temperature'] = []
  for i, ihdu in enumerate(hdulist):
    # Get temperature (there are more sophisticated ways to do this, this should just work for what you need)
    T = 10**(4+(0.1*(ihdu-2)))


    res['temperature'].append(T)
    res['power'][i] = calc_power_oneT(Zlist, Elo, Ehi, ihdu, linefile = linefile,\
                      cocofile = cocofile)


  return res

if __name__=='__main__':

  #### ADJUST THINGS HERE

  # Elements to include
  #Zlist = range(1,31) <- all the elements
  Zlist = [1,2,6,7,8,10,12,13,14,16,18,20,26,28] #<- just a few

  # specify energy range you want to integrate over (min = 0.001keV, max=100keV)
  Elo = 2 #keV
  Ehi = 10 #

  # specify output file name (default output.txt)
  outfile = 'output.txt'

  #pre-open the emissivity files (not required, but saves a lot of disk access time)
  linedata  = pyatomdb.pyfits.open(os.path.expandvars('$ATOMDB/apec_line.fits'))
  cocodata  = pyatomdb.pyfits.open(os.path.expandvars('$ATOMDB/apec_coco.fits'))

  #### END ADJUST THINGS HERE

  # crunch the numbers
  k = calc_power(Zlist, Elo, Ehi, linefile = linedata, cocofile = cocodata)

  # output generation
  o = open(outfile, 'w')

  # header row
  s = '# Temperature log10(K)'
  for i in range(len(Zlist)):
    s += ' %12i'%(Zlist[i])
  o.write(s+'\n')

  # for each temperature
  for i in range(len(k['temperature'])):
    s = '%22e'%(numpy.log10(k['temperature'][i]))
    for Z in Zlist:
      s+=' %12e'%(k['power'][i][Z])
    o.write(s+'\n')

  # notes
  o.write("# Total Emissivity in keV cm^3 s^-1 for each element with AG89 abundances, between %e and %e keV\n"%(Elo, Ehi))
  o.write("# To get cooling power, multiply by Ne NH")
  o.close


