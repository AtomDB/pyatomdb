__all__=["atomdb","util","atomic","spectrum","const","apec"]

from .atomdb import *
from . import spectrum
from . import atomic
import ctypes
from . import apec
import sys, glob

__version__="0.7.2"

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
