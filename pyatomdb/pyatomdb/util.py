"""
util.py contains a range of miscellaneous helper codes that assist in running
other AtomDB codes but are not in any way part of a physical calculation.

Version -.1 - initial release
Adam Foster July 17th 2015
"""

import numpy, os, errno


################################################################################
#
#  Python Module
#
#  Name:        arflib.py
#
#  Decription:  Python codes that I use.
#
#  Module contents (and 1 line description: see individual modules for more):
#
#     differentiate
#          Differentiate line.
#
#  First Version:
#       Adam Foster, 07-Jun-2010
#
################################################################################


#*******************************************************************************
#
#  Routine fig
#
#  Differentiates vector y wrt x
#
#  input: x,y, lowend=1, highend=1
#
#   x - numpy.array of x values
#   y - numpy.array of y values
#   lowend - 0,1,2: 0 set dy/dx=0 at low end
#                   1 set d2y/dx2=0 at low end
#                   2 set d3y/dx3=0 at low end
#   highend - 0,1,2: 0 set dy/dx=0 at high end
#                   1 set d2y/dx2=0 at high end
#                   2 set d3y/dx3=0 at high end
#
#  returns: numpy array of derivatives
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************

def figcoords(lowxpix, lowypix, highxpix, highypix, 
                  lowxval, lowyval, highxval, highyval, 
                  xpix, ypix, logx=False, logy=False):



  if logx:
    lxval = numpy.log10(lowxval)
    hxval = numpy.log10(highxval)
  else:
    lxval = lowxval*1.0
    hxval = highxval*1.0
  if logy:
    lyval = numpy.log10(lowyval)
    hyval = numpy.log10(highyval)
  else:
    lyval = lowyval*1.0
    hyval = highyval*1.0
  
  dx = (hxval-lxval)/(highxpix-lowxpix)
  dy = (hyval-lyval)/(highypix-lowypix)
  
#  print dx, dy
  xout = lxval + (xpix-lowxpix)*dx
  yout = lyval + (ypix-lowypix)*dy

  if logx:
    xout = 10**xout
  if logy:
    yout = 10**yout
  return xout, yout


def unique(s):
     """Return a list of the elements in s, but without duplicates.

     For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
     unique("abcabc") some permutation of ["a", "b", "c"], and
     unique(([1, 2], [2, 3], [1, 2])) some permutation of
     [[2, 3], [1, 2]].

     For best speed, all sequence elements should be hashable.  Then
     unique() will usually work in linear time.

     If not possible, the sequence elements should enjoy a total
     ordering, and if list(s).sort() doesn't raise TypeError it's
     assumed that they do enjoy a total ordering.  Then unique() will
     usually work in O(N*log2(N)) time.

     If that's not possible either, the sequence elements must support
     equality-testing.  Then unique() will usually work in quadratic
     time.
     """

     n = len(s)
     if n == 0:
         return []

     # Try using a dict first, as that's the fastest and will usually
     # work.  If it doesn't work, it will usually fail quickly, so it
     # usually doesn't cost much to *try* it.  It requires that all the
     # sequence elements be hashable, and support equality comparison.
     u = {}
     try:
         for x in s:
             u[x] = 1
     except TypeError:
         del u  # move on to the next method
     else:
         return u.keys()

     # We can't hash all the elements.  Second fastest is to sort,
     # which brings the equal elements together; then duplicates are
     # easy to weed out in a single pass.
     # NOTE:  Python's list.sort() was designed to be efficient in the
     # presence of many duplicate elements.  This isn't true of all
     # sort functions in all languages or libraries, so this approach
     # is more effective in Python than it may be elsewhere.
     try:
         t = list(s)
         t.sort()
     except TypeError:
         del t  # move on to the next method
     else:
         assert n > 0
         last = t[0]
         lasti = i = 1
         while i < n:
             if t[i] != last:
                 t[lasti] = last = t[i]
                 lasti += 1
             i += 1
         return t[:lasti]

     # Brute force is all that's left.
     u = []
     for x in s:
         if x not in u:
             u.append(x)
     return u
     
     

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise
