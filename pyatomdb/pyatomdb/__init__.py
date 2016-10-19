__all__=["atomdb","util","atomic","spectrum","const","apec"]

from .atomdb import *
import spectrum
import atomic
import ctypes
import apec
__version__="0.0.1.9"
PATH = os.path.dirname(__file__)
liblinapprox_file = os.path.join(PATH,"../linear_approx.dylib")
if not os.path.isfile(liblinapprox_file):
  liblinapprox_file = os.path.join(PATH,"../linear_approx.so")
 
liblinapprox = ctypes.CDLL( liblinapprox_file, ctypes.RTLD_GLOBAL)
