"""
atomic.py contains routines related to basic atomic data, e.g. converting
integer nuclear charge to element symbols, etc.

Version -.1 - initial release
Adam Foster July 17th 2015
"""

import re, numpy
################################################################################
#
#  Python Module
#
#  Name:        adbatomic.py
#
#  Decription:  Codes for simple atomic data related tasks
#
#  Module contents (and 1 line description: see individual modules for more):
#
#     z0toelsymb
#          Converts z0 to element symbol (eg 6 -> C)
#
#     z0toelname
#          Converts z0 to element name (eg 6 -> Carbon)
#
#     int2roman
#          Converts integer to Roman Numerals (eg 6 -> VI)
#
#     spectroscopic_name
#          Converts z0,ionstage to spectroscopic name (eg 6,3 -> C IV)
#
#  Check individual codes for author and update details
#
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
################################################################################


#*******************************************************************************
#
#  Routine z0toelsymb
#
#  Converts z0 to element symbol
#
#  input: z0 (integer)
#
#  returns: Element symbol (first letter capitalised)
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************

def Ztoelsymb(Z) :
  """
  Returns element symbol of element with nuclear charge Z.
  
  INPUTS
  Z  - nuclear charge of element (e.g. 6 for carbon)
  
  RETURNS
  element symbol (e.g. "C" for carbon)
  
  Version 0.1 28 July 2009
  Adam Foster
  """

  elsymb=('H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
          'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
          'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
          'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
          'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
          'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
          'Pa', 'U')

  if Z < 1 :
      print "Z must be between 1 and 92. You have given Z= " + repr(z0)
      ret=-1
  elif Z > 92 :
      print "Z must be between 1 and 92. You have given Z= " + repr(z0)
      ret=-1
  else :
      ret=elsymb[Z-1]
  return ret

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

def z0toelsymb(z0):
  """
  Returns element symbol of element with nuclear charge z0.
  (wrapper to Ztoelsymb for compatibility purposes)
  
  Parameters
  ----------
  z0 : int
    nuclear charge of element (e.g. 6 for carbon)
  
  Returns
  -------
  str
    element symbol (e.g. "C" for carbon)
  """
#  
#  Version 0.1 28 July 2009
#  Adam Foster
#
  ret = Ztoelsymb(z0)

  return ret

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************


#*******************************************************************************
#
#  Routine z0toelname
#
#  Converts z0 to element name
#
#  input: z0 (integer)
#
#  returns: Element name (first letter capitalised)
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************


def z0toelname(z0):
  """
  Returns element name of element with nuclear charge z0.
  (wrapper to Ztoelname for compatibility purposes)
  
  Parameters
  ----------
  
  z0 :  int
    nuclear charge of element (e.g. 6 for carbon)
  
  Returns
  -------
  str
    element name (e.g. "Carbon" for carbon)
  """
#
#  Version 0.1 28 July 2009
#  Adam Foster
#

  ret = Ztoelname(Z)
  return ret

def Ztoelname(Z):
  """
  Returns element name of element with nuclear charge Z.
  
  Parameters
  ----------
  
  Z :  int
    nuclear charge of element (e.g. 6 for carbon)
  
  Returns
  -------
  str
    element name (e.g. "Carbon" for carbon)
  """
#
#  Version 0.1 28 July 2009
#  Adam Foster
#

  elname=('Hydrogen'    , 'Helium'      , 'Lithium'     , 'Beryllium'   ,
          'Boron'       , 'Carbon'      , 'Nitrogen'    , 'Oxygen'      ,
          'Fluorine'    , 'Neon'        , 'Sodium'      , 'Magnesium'   ,
          'Aluminum'    , 'Silicon'     , 'Phosphorus'  , 'Sulfur'      ,
          'Chlorine'    , 'Argon'       , 'Potassium'   , 'Calcium'     ,
          'Scandium'    , 'Titanium'    , 'Vanadium'    , 'Chromium'    ,
          'Manganese'   , 'Iron'        , 'Cobalt'      , 'Nickel'      ,
          'Copper'      , 'Zinc'        , 'Gallium'     , 'Germanium'   ,
          'Arsenic'     , 'Selenium'    , 'Bromine'     , 'Krypton'     ,
          'Rubidium'    , 'Strontium'   , 'Yttrium'     , 'Zirconium'   ,
          'Niobium'     , 'Molybdenum'  , 'Technetium'  , 'Ruthenium'   ,
          'Rhodium'     , 'Palladium'   , 'Silver'      , 'Cadmium'     ,
          'Indium'      , 'Tin'         , 'Antimony'    , 'Tellurium'   ,
          'Iodine'      , 'Xenon'       , 'Cesium'      , 'Barium'      ,
          'Lanthanum'   , 'Cerium'      , 'Praseodymium', 'Neodymium'   ,
          'Promethium'  , 'Samarium'    , 'Europium'    , 'Gadolinium'  ,
          'Terbium'     , 'Dysprosium'  , 'Holmium'     , 'Erbium'      ,
          'Thulium'     , 'Ytterbium'   , 'Lutetium'    , 'Hafnium'     ,
          'Tantalum'    , 'Tungsten'    , 'Rhenium'     , 'Osmium'      ,
          'Iridium'     , 'Platinum'    , 'Gold'        , 'Mercury'     ,
          'Thallium'    , 'Lead'        , 'Bismuth'     , 'Polonium'    ,
          'Astatine'    , 'Radon'       , 'Francium'    , 'Radium'      ,
          'Actinium'    , 'Thorium'     , 'Protactinium', 'Uranium')

  if Z < 1 :
    print "Z must be between 1 and 92. You have given Z= " + repr(Z)
    ret=-1
  elif Z > 92 :
    print "Z must be between 1 and 92. You have given Z= " + repr(Z)
    ret=-1
  else :
    ret=elname[Z-1]
  return ret

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

#*******************************************************************************
#
#  Routine int2roman
#
#  Converts a number to roman numeral (pilfered off the internet)
#
#  input: number (integer)
#
#  returns: Roman numeral (string)
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************

def int2roman(number):

    numerals = { 1   : "I" , 4   : "IV", 5    : "V" , 9   : "IX", 10  : "X" ,
                 40  : "XL", 50  : "L" , 90   : "XC", 100 : "C" , 400 : "CD",
                 500 : "D" , 900 : "CM", 1000 : "M" }
    result = ""

    for value, numeral in sorted(numerals.items(), reverse=True):
        while number >= value:
            result += numeral
            number -= value
    return result

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#
#  Routine int_to_roman
#
#  Converts a integer to a roman numeral (pilfered off the internet)
#
#  input: Roman numeral (string)
#
#  returns: number (integer)
#
#  First Version:
#       Adam Foster, 03-Nov-2011
#
#*******************************************************************************
def int_to_roman(input):
   """
   Convert an integer to Roman numerals.

   """
   if type(input) != type(1):
      raise TypeError, "expected integer, got %s" % type(input)
   if not 0 < input < 4000:
      raise ValueError, "Argument must be between 1 and 3999"   
   ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   result = ""
   for i in range(len(ints)):
      count = int(input / ints[i])
      result += nums[i] * count
      input -= ints[i] * count
   return result

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#
#  Routine roman2int
#
#  Converts a roman numeral to an integer (pilfered off the internet)
#
#  input: Roman numeral (string)
#
#  returns: number (integer)
#
#  First Version:
#       Adam Foster, 03-Nov-2011
#
#*******************************************************************************
def roman_to_int(input):
   """
   Convert a roman numeral to an integer.
   """
   if type(input) != type(""):
      raise TypeError, "expected string, got %s" % type(input)
   input = input.upper()
   nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
   ints = [1000, 500, 100, 50,  10,  5,   1]
   places = []
   for c in input:
      if not c in nums:
         raise ValueError, "input is not a valid roman numeral: %s" % input
   for i in range(len(input)):
      c = input[i]
      value = ints[nums.index(c)]
      # If the next place holds a larger number, this value is negative.
      try:
         nextvalue = ints[nums.index(input[i +1])]
         if nextvalue > value:
            value *= -1
      except IndexError:
         # there is no next place.
         pass
      places.append(value)
   sum = 0
   for n in places: sum += n
   # Easiest test for validity...
   if int_to_roman(sum) == input:
      return sum
   else:
      raise ValueError, 'input is not a valid roman numeral: %s' % input
#*******************************************************************************
#
#  Routine spectroscopic_name
#
#  Converts z0, ioncharge to element symbol
#
#  input: z0 (integer), ioncharge (integer)
#
#  returns: Spectroscopic name (e.g. C IV or Mg XI)
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************

def spectroscopic_name(Z,z1) :
  """
  Converts Z,z1 to spectroscopic name, e.g. 6,5 to "C V"
  
  Parameters
  ----------
  
  Z : int
    nuclear charge (e.g. 6 for C)
  z1 : int
    ion charge +1 (e.g. 5 for C4+)
   
  Returns
  -------
  str
    spectroscopic symbol for ion (e.g. "C V" for C+4)
  
  """
#
#  Version 0.1 28 July 2009
#  Adam Foster 
#    
    # get element symbol
  elsymb = Ztoelsymb(Z)

    # convert z1 to spectroscopic

  roman = int2roman(z1)
  ret = elsymb + ' ' + roman

  return ret


#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

#*******************************************************************************
#
#  Routine spectroscopictoz0
#
#  Converts z0, ioncharge to element symbol
#
#  input: z0 (integer), ioncharge (integer)
#
#  returns: Spectroscopic name (e.g. C IV or Mg XI)
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************

def spectroscopictoz0(name):
  """
  Converts spectroscopic name to Z, z1, e.g. "C V" to 6,5
  
  Parameters
  ----------
  
  name : str
    Ion name, e.g. "C V" 
  
  Returns
  -------
  int, int
    Z, z1 for the ion. (e.g. 6,5 for C V)
  
  """
#  
#  Version 0.1 28 July 2009
#  Adam Foster 
#
# convert name (e.g. Fe VIII) to z0 & ioncharge (=0 for neutral)

    # get element symbol
  d = name.split()
  
  elsymb = d[0]
  chargesymb = d[1]
  
  z0 = elsymb_to_z0(elsymb)
  
  z1 = roman_to_int(chargesymb)
  z=z1-1
  
  return z0,z



#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#
#  Routine occup_to_cfg
#
#  Converts occupancy vector to configuration string
#                     (e.g. [2,1,0,1] -> 1s2 2s1 3s1)
#
#  input: occupancy (list)
#
#  returns: configuration string
#
#  First Version:
#       Adam Foster, 24-Nov-2009
#
#*******************************************************************************

def occup_to_cfg(occlist) :
  l_list = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r',
            't','u','v','w','x','y','z']

  cfgstr=''
  l=0
  n=0
  for j, i in enumerate(occlist) :
#    print repr(i) +' '+ repr(n) +' '+ repr(l)
    if (l+1 >= n):
      l = 0
      n += 1
    else:
      l += 1
    if (i > 0):
      cfgstr = cfgstr+' '+repr(n)+l_list[l]+repr(i)
    
# return minus leading blank
  
  return cfgstr.strip()
      
    
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

#*******************************************************************************
#
#  Routine elsymb_to_z0
#
#  Converts occupancy element symbol to z0
#                     (e.g. 'He' -> 2) (case insensitive)
#
#  input: occupancy (list)
#
#  returns: configuration string
#
#  First Version:
#       Adam Foster, 24-Nov-2009
#
#*******************************************************************************

def elsymb_to_Z(elsymb) :
  """
  Converts element symbol to nuclear charge, e.g. "C" -> 6
  
  Parameters
  ----------
  
  elsymb : str
    Element symbol, e.g. "C". Case insensitive.

  Returns
  -------
  int
    Z for the ion. (e.g. 6 for C)
  """
#  
#  Version 0.1 28 July 2009
#  Adam Foster 
#
  ellist=('h' , 'he', 'li', 'be', 'b' , 'c' , 'n' , 'o' , 'f' , 'ne',
          'na', 'mg', 'al', 'si', 'p' , 's' , 'cl', 'ar', 'k' , 'ca',
          'sc', 'ti', 'v' , 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn',
          'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y' , 'zr',
          'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn',
          'sb', 'te', 'i ', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd',
          'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb',
          'lu', 'hf', 'ta', 'w' , 're', 'os', 'ir', 'pt', 'au', 'hg',
          'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th',
          'pa', 'u')

  try:
    ind=ellist.index(elsymb.lower().strip())
  except ValueError:
    print "elsymb_to_z0 error: invalid element symbol '"+elsymb+"', returning -1"
    ind=-1
  
  return ind+1
      
def elsymb_to_z0(elsymb) :
  """
  Converts element symbol to nuclear charge, e.g. "C" -> 6
  (wrapper to elsymb_to_Z, retained for consistency)
  
  Parameters
  ----------
  
  elsymb : str
    Element symbol, e.g. "C". Case insensitive.

  Returns
  -------
  int
    Z for the ion. (e.g. 6 for C)
  """
  ret= elsymb_to_Z(elsymb)
  return ret

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

#*******************************************************************************
#
#  Routine z0_to_mass
#
#  Return atomic mass of element with atomic number z0
#
#  input: z0
#
#  returns: atomic mass (float)
#
#  First Version:
#       Adam Foster, 4-Apr-2010
#
#*******************************************************************************
def z0_to_mass(z0):
  """
  Converts element symbol to atomic mass, e.g. "C" -> 12.0107
  
  (wrapper to Z_to_mass, retained for consistency)
  
  Isotope fractions based on those found in earth's crust samples, your
  astrophysical object may vary.
  
  Parameters
  ----------
  z0 : int
    nuclear charge, e.g 6 for C

  Returns
  -------
  float
    mass in a.m.u. for the element. (e.g. 12.0107 for C)

  References
  ----------
  Atomic masses are taken from:
  Pure Appl. Chem. 81 NO 11, 2131-2156 (2009)
  Masses for Technetium, Promethium, Polonium, Astatine, Radon, 
  Francium, Radium & Actinum are estimates. If you need these you
  probably aren't doing astronomy...

  """
#  
#  Version 0.1 28 July 2009
#  Adam Foster 
#
  ret = Z_to_mass(z0)
  return ret

def Z_to_mass(Z):
  """
  Converts element symbol to atomic mass, e.g. "C" -> 12.0107
  
  Isotope fractions based on those found in earth's crust samples, your
  astrophysical object may vary.
  
  Parameters
  ----------
  Z : int
    nuclear charge, e.g 6 for C

  Returns
  -------
  float
    mass in a.m.u. for the element. (e.g. 12.0107 for C)
  
  References
  ----------
  Atomic masses are taken from:
  Pure Appl. Chem. 81 NO 11, 2131-2156 (2009)
  Masses for Technetium, Promethium, Polonium, Astatine, Radon, 
  Francium, Radium & Actinum are estimates. If you need these you
  probably aren't doing astronomy...
  
  """
#  Version 0.1 28 July 2009
#  Adam Foster 
#
  masslist=( 1.00794   ,   4.002602,   6.941    ,   9.012182 ,  10.811    ,
            12.0107    ,  14.0067  ,  15.9994   ,  18.9984032,  20.1797   ,
            22.98976928,  24.3050  ,  26.9815386,  28.0855   ,  30.973762 ,
            32.065     ,  35.453   ,  39.948    ,  39.0983   ,  40.078    ,
            44.955912  ,  47.867   ,  50.9415   ,  51.9961   ,  54.938045 ,
            55.845     ,  58.933195,  58.6934   ,  63.546    ,  65.38     ,
            69.723     ,  72.64    ,  74.92160  ,  78.96     ,  79.904    ,
            83.798     ,  85.4678  ,  87.62     ,  88.90585  ,  91.224    ,
            92.90638   ,  95.96    ,  98.000    , 101.07     , 102.90550  ,
           106.42      , 107.8682  , 112.411    , 114.818    , 118.710    ,
           121.760     , 127.60    , 126.90447  , 131.293    , 132.9054519,
           137.327     , 138.90547 , 140.116    , 140.90765  , 144.242    ,
           145.000     , 150.36    , 151.964    , 157.25     , 158.92535  ,
           162.500     , 164.93032 , 167.259    , 168.93421  , 173.054    ,
           174.9668    , 178.49    , 180.94788  , 183.84     , 186.207    ,
           190.23      , 192.217   , 195.084    , 196.966569 , 200.59     ,
           204.3833    , 207.2     , 208.98040  , 209.000    , 210.000    ,
           222.000     , 223.000   , 226.000    , 227.00     , 232.03806  ,
           231.03588   , 238.02891)

  if Z < 1 :
    print "Z must be between 1 and 92. You have given Z= " + repr(Z)
    ret=-1
  elif Z > 92 :
    print "Z must be between 1 and 92. You have given Z= " + repr(Z)
    ret=-1
  else :
    ret=masslist[Z-1]
  return ret

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def config_to_occup(cfgstr, nel=-1, shlmax=-1, noccup=[-1]):
#  print cfgstr
  if len(cfgstr)==0:
    cfgstr = '1s2'
  cfgsplit = cfgstr.split(' ')
  n = []
  l = []
  o = []
  llist='spdfghiklmnoqrt'
  for cfg in cfgsplit:
#    print cfg
#    print n,l,o
    ntmp = re.search("^[0-9]+",cfg)
    n.append(int(ntmp.group(0)))
    ltmp = re.search("[a-zA-Z]",cfg)
    l.append(llist.index(ltmp.group(0)))
    otmp = re.search("[0-9]+$",cfg)
    o.append(int(otmp.group(0)))
    
  # find the max nl shell
  if shlmax == -1:
    maxshl = -1
#    print 'boo'
    for i in range(len(n)):
      shlind = 0
#      print 'b2'
      shlind=sum(range(1,n[i]+1))
      maxshl = max([maxshl, shlind])
  else:
    maxshl=shlmax
#  print maxshl    
  occup = numpy.zeros(maxshl, dtype=int)
#  print occup
  for i in range(len(n)):
    shlind = 0
    if n[i] > 1:
      for iin in range(1,n[i]):
        shlind = shlind + iin
      shlind = shlind + l[i]
#    print shlind, n[i], l[i], o[i]  
#    print occup, shlind
    occup[shlind] = occup[shlind] + o[i]
       
  inext = 0
  lnext = 0
  nnext = 1
  onext = 2
  if noccup[0]==-1:

    firstoccup = min(numpy.where(occup>0)[0])
    if firstoccup > 0:
#    print "HI2!"
#    print occup
      for i in range(len(occup)):
        if ((occup[i] == 0) &(occup.sum() < nel)):
  #        print "arse",occup.sum()+4*lnext+2, nel
          if (occup.sum()+4*lnext+2 <= nel):
  #          print "PING:"
            occup[i] = occup[i] + 4*lnext+2
  #          print occup
        else:
          break
      
        
        if nnext-lnext == 1:
          nnext += 1
          lnext = 0
        else: 
          lnext += 1

  
    inext = 0
    lnext = 0
    nnext = 1
    onext = 2
  #  print "HI!"
  #  print (sum(occup), nel)
    while (sum(occup) < nel):
  #    print occup
      if occup[inext] == 0:
        if (onext > (nel-sum(occup))):
          occup[inext] += nel-sum(occup)
        else:
          occup[inext] += onext
      
      if nnext-lnext == 1:
        nnext += 1
        lnext = 0
      else: 
        lnext += 1
      onext = 4*lnext+2
      inext += 1
  else:
##    print "WHAT!", occup
    # we have an array defining the number of electrons total in each N shell
    # such as in FAC
    inext = 0
    nnext = 1
    shell_n = numpy.zeros(len(occup), dtype=int)
    shell_l = numpy.zeros(len(occup), dtype=int)
    while inext < len(shell_n):
      shell_n[inext:inext+nnext]=nnext
      shell_l[inext:inext+nnext]=numpy.arange(nnext)
      inext += nnext
      nnext += 1
##    print "shell_n",shell_n
##    print "shell_l",shell_l
    for i_n in range(len(noccup)):
      nnext = i_n+1
      i = numpy.where(shell_n == nnext)[0]
      nel_tot = sum(occup[i])
      nel_targ = noccup[i_n]
##      print "occup: ", occup
##      print "Gives nel_tot=%i for n=%i, while nel_targ=%i" %\
##            (nel_tot, nnext, nel_targ)
      # if the number of electrons match
      if nel_tot == nel_targ: continue
      
      if nel_tot > nel_targ:
        print "ERROR: more electron in n=%i shell than there should be for %s" %\
            (nnext, cfgstr)
        print "   %i vs %i" %(nel_tot, nel_targ)
        
      while nel_tot< nel_targ:
##        print occup
        #find empty l shells
        lposs = []
        for il in i:
##          print "Try n = %i, l=%i, il=%i"%(nnext, shell_l[il], il)
##          print "shell_l[il]*4+2=%i, nel_targ-nel_tot=%i" % \
##                (shell_l[il]*4+2,nel_targ-nel_tot)
##          print "occup[%i] = %i" %(il, occup[il])
          if ((occup[il]==0) &(shell_l[il]*4+2 <= (nel_targ-nel_tot))): 
            lposs.append(shell_l[il])
##            print "Good!"
        # get number of occupancies
        shell_occup = numpy.array(lposs)*4+2
        
        delta_nel = nel_targ - nel_tot
        k = numpy.where(shell_occup == delta_nel)[0]
##        print "k=",k
        if len(k) ==1:
          k = k[0]
          ind = numpy.where((shell_n==nnext) & (shell_l==lposs[k]))[0][0]
##          print "setting occup[%i]=shell_occup[%i]=%i"% \
##                (ind, k, shell_occup[k])
          occup[ind] =shell_occup[k]
##          print occup
        else:
##          print nnext, lposs
          ind = numpy.where((shell_n==nnext) & (shell_l==lposs[0]))[0][0]
          occup[ind] = shell_occup[0]
          
        nel_tot = sum(occup[i])
##        zzz=raw_input()
      
    
  if ((nel > 0) & (sum(occup) != nel)):
#    print occup
    return occup,False
  else:
    return occup, True
    

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def occup_to_config(occup):
  llist='spdfghiklmnopqrstuvwxyz'
  
  s = ''
  lnext = 0
  nnext = 1
  for i,j in enumerate(occup):
    if j > 0:
      s = s+ repr(nnext)+llist[lnext]+repr(j)+' '
    
    if nnext-lnext==1:
      lnext = 0
      nnext += 1
    else:
      lnext=lnext+1
  s = s[:-1]
  return s
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def parse_config(cfgstr):
 # returns n shell, l shell and occupancy for each part of the configuration
 # e.g. [[1,0,2],[2,1,1]] for 1s2 2p1
 
  #split on space
  c = cfgstr.split()
  
  
  llist= 'spdfghiklmnoqrtuvwxyz'
  
  ret=[]
  for ic in c:
    cfg = []
    ntmp = re.search("^[0-9]+",ic)
    cfg.append(int(ntmp.group(0)))
    ltmp = re.search("[a-zA-Z]",ic)
    cfg.append(llist.index(ltmp.group(0)))
    otmp = re.search("[0-9]+$",ic)
    cfg.append(int(otmp.group(0)))
   
    ret.append(cfg)    
  return ret
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def get_parity(cfgstr):
  d = parse_config(cfgstr)
  
  evenparity = True
  
  for i in d:
    if i[1]*i[2] % 2 == 1:
      evenparity = not(evenparity)

  if evenparity:
    return 0
  else:
    return 1      
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
    
def get_maxn(cfgstr):
  d = parse_config(cfgstr)
  
  maxn = max([c[0] for c in d])

  return maxn
