__all__=["atomdb","util","atomic","spectrum","const","apec"]

from .atomdb import *
from . import spectrum
from . import atomic
import ctypes
from . import apec
from . import util
import sys, glob

__version__="0.10.11"

try:
  PATH = os.path.dirname(__file__)
  g = glob.glob("%s/../linear_approx*.dylib"%(PATH))
  if len(g) ==1:
    liblinapprox_file = os.path.join(PATH,g[0])
  else:
    g = glob.glob("%s/../linear_approx*.so"%(PATH))
    liblinapprox_file = os.path.join(PATH,g[0])

  liblinapprox = ctypes.CDLL( liblinapprox_file, ctypes.RTLD_GLOBAL)
except (OSError, IndexError):
  on_rtd=os.environ.get('READTHEDOCS')=='True'
  if on_rtd:
    pass
  else:
    raise

# check initialization
on_rtd=os.environ.get('READTHEDOCS')=='True'
if on_rtd:
  pass
else:
  if os.environ.get('ATOMDB')==None:
    # issues here
    install = ''
    while not install in ['yes','no']:
      install = util.question("ATOMDB environment variable not set. Set up AtomDB files now?", "yes",multichoice=["yes","no","info"])
      if install=='info':
        print("If yes, this will run the pyatomdb.util.initialize() script, which installs the necessary data files to run pyatomdb")
    if install=='yes':
      util.initialize()
          
  elif not os.path.exists(os.path.expandvars("$ATOMDB/userdata")):
    install = ''
    while not install in ['yes','no']:
      install = util.question("ATOMDB Environment Variable Not Set. Set up AtomDB files now?", "yes",multichoice=["yes","no","info"])
      if install=='info':
        print("If yes, this will run the pyatomdb.util.initialize() script, which installs the necessary data files to run pyatomdb")

    if install=='yes':
      util.initialize()
  
    

    
