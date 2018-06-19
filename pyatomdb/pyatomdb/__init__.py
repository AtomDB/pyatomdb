__all__=["atomdb","util","atomic","spectrum","const","apec"]

from .atomdb import *
from . import spectrum
from . import atomic
import ctypes
from . import apec
import sys
__version__="0.5.3"

try:
  PATH = os.path.dirname(__file__)
  liblinapprox_file = os.path.join(PATH,"../linear_approx.dylib")
  if not os.path.isfile(liblinapprox_file):
    liblinapprox_file = os.path.join(PATH,"../linear_approx.so")
   
  liblinapprox = ctypes.CDLL( liblinapprox_file, ctypes.RTLD_GLOBAL)
except OSError:
  on_rtd=os.environ.get('READTHEDOCS')=='True'
  if on_rtd:
    pass
  else:
    raise
