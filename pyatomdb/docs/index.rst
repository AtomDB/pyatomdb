.. PyAtomDB documentation master file, created by
   sphinx-quickstart on Fri Sep 18 17:07:36 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyAtomDB
====================================

Contents:

.. toctree::
   :maxdepth: 2
   
   apec
   atomic
   atomdb
   spectrum
   util
   

============
Introduction
============





pyatomdb is a selection of utilities designed to interact with the AtomDB
database. Currently, these utilities are in a far from finished format, however
progress on this is ongoing. These have been taken from a series of codes on my
laptop which were useful. Some produce lots of unhelpful onscreen output.


There are several different modules currently. These are:

-  :doc:`atomdb </atomdb>`  : a series of codes for interacting with the AtomDB atomic database
- :doc:`atomic </atomic>`  : basic atomic data routines - e.g. converting element symbols to atomic number, etc.
- :doc:`const </const>`   : a series of physical constants
- :doc:`spectrum </spectrum>` : routines for generating spectra from the published AtomDB line and continuum emissivity files
- :doc:`util </util>`    : sumple utility codes (sorting etc) that pyatomdb relies on.

Currently, only the spectrum library has been extensively tested. Expect bugs.
Report those bugs!


=======
License
=======
Pyatomdb is released under the Smithsonian License:

Copyright 2015 Smithsonian Institution. Permission is granted to use, copy, 
modify, and distribute this software and its documentation for educational,
research and non-profit purposes, without fee and without a signed
licensing agreement, provided that this notice, including the following
two paragraphs, appear in all copies, modifications and distributions.
For commercial licensing, contact the Office of the Chief Information
Officer, Smithsonian Institution, 380 Herndon Parkway, MRC 1010, Herndon,
VA. 20170, 202-633-5256.

This software and accompanying documentation is supplied "as is" without
warranty of any kind. The copyright holder and the Smithsonian
Institution: (1) expressly disclaim any warranties, express or implied,
including but not limited to any implied warranties of merchantability,
fitness for a particular purpose, title or non-infringement; (2) do not
assume any legal liability or responsibility for the accuracy,
completeness, or usefulness of the software; (3) do not represent that use
of the software would not infringe privately owned rights; (4) do not
warrant that the software is error-free or will be maintained, supported,
updated or enhanced; (5) will not be liable for any indirect, incidental,
consequential special or punitive damages of any kind or nature,
including but not limited to lost profits or loss of data, on any basis
arising from contract, tort or otherwise, even if any of the parties has
been warned of the possibility of such loss or damage.


=====
Usage
=====

--------------------------
Interrogating the database
--------------------------
This makes use of the atomdb module. Before starting, it is a good
idea to download the AtomDB tarballs from www.atomdb.org/Downloads: the latest one is version 3.0.2. This dataset contains the filemap, which identifies the correct files for each different ion and process.

It is also important to set the ATOMDB environment variable to point to where this data was untarred, e.g. the folder with the filemap in it, so $ATOMDB/filemap should exist.
e.g::

 bash: export ATOMDB=/myfolder/atomdb_v3.0.2
 csh: setenv ATOMDB /myfolder/atomdb_v3.0.2

or within python:: 

 import os
 os.environ['ATOMDB']='/myfolder/atomdb_v3.0.2'

Currently, the AtomDB database is more than 10GB of data, so we are avoiding distributing it to all users. You can, however, get the individual data you need using the ``get_data`` routine::

  lvdata = pyatomdb.atomdb.get_data(z0, z1, ftype)

This will try to open the file locally if it exists, and if it does not it will then go to the AtomDB FTP server and download the data for element z0, ion z1, with ftype a 2-character string denoting the type of data to get:

- ``IR``: ionization and recombination
- ``LV``: energy levels
- ``LA``: radiative transition data (lambda and A-values)
- ``EC``: electron collision data
- ``PC``: proton collision data
- ``DR``: dielectronic recombination satellite line data
- ``PI``: XSTAR photoionization data
- ``AI``: autoionization data

So to open the energy levels for oxygen with 2 electrons (O 6+, or O VII)::
  
 lvdata = pyatomdb.atomdb.get_data(8,7,'LV')

Data files are stored in ``$ATOMDB/APED/<elsymb>/<elsymb>_<ionnum>/``

-----------------
Making a Spectrum
-----------------
These functions are in the ``spectrum`` module::
  
  import pyatomdb, numpy, pylab
  
  # set up a grid of energy bins to model the spectrum on:
  ebins=numpy.linspace(0.3,10,1000)
  
  # define a broadening, in keV, for the lines
  de = 0.01
  
  # define the temperature at which to plot (keV)
  te = 3.0
  
  # find the index which is closest to this temperature
  ite = pyatomdb.spectrum.get_index( te, teunits='keV', logscale=False)
  
  # create both a broadened and an unbroadened spectrum
  a = pyatomdb.spectrum.make_spectrum(ebins, ite,dummyfirst=True)
  b = pyatomdb.spectrum.make_spectrum(ebins, ite, broadening=de, \
                                      broadenunits='kev',dummyfirst=True)
  # The dummyfirst argument adds an extra 0 at teh beginning of the
  # returned array so it is the same length as ebins. It allows
  # accurate plotting using the "drawstyle='steps'" flag to plot.
  
  # plot the results
  fig = pylab.figure()
  fig.show()
  ax = fig.add_subplot(111)
  
  ax.loglog(ebins, a, drawstyle='steps', label='Unbroadened')
  ax.loglog(ebins, b, drawstyle='steps', label='sigma = %.2f'%(de))
  ax.set_xlabel('Energy (keV)')
  ax.set_ylabel('Emissivity (ph cm$^{3}$ s$^{-1}$ bin$^{-1}$)')
  ax.legend(loc=0)
  pylab.draw()
  zzz = raw_input("Press enter to continue")
  
  print "Listing lines between 1 and 2 A"
  # now list the lines in a wavelength region
  llist = pyatomdb.spectrum.list_lines([1,2.0], index=ite)
  # print these to screen
  pyatomdb.spectrum.print_lines(llist)
  # print to screen, listing the energy, not the wavelength
  print "Listing lines between 1 and 2 A, using keV."
  
  pyatomdb.spectrum.print_lines(llist, specunits = 'keV')





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

