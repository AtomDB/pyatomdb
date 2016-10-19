.. PyAtomDB documentation master file, created by
   sphinx-quickstart on Fri Sep 18 17:07:36 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyAtomDB
====================================


============
Introduction
============





PyAtomDB is a selection of utilities designed to interact with the `AtomDB
database <https://www.atomdb.org>`_ . Currently, these utilities are in a far from finished format, however progress on this is ongoing. These have been taken from a series of codes on my laptop which were useful. Some produce lots of unhelpful onscreen output.


There are several different modules currently. These are:

- :doc:`atomdb </atomdb>`  : a series of codes for interacting with the AtomDB atomic database
- :doc:`atomic </atomic>`  : basic atomic data routines - e.g. converting element symbols to atomic number, etc.
- :doc:`const </const>`   : a series of physical constants
- :doc:`spectrum </spectrum>` : routines for generating spectra from the published AtomDB line and continuum emissivity files
- :doc:`util </util>`    : sumple utility codes (sorting etc) that pyatomdb relies on.
- :doc:`apec </apec>`  : ultimately, the full apec code. For now, incomplete.

Expect bugs. Report those bugs! Make feature requests! Email the code authors or raise an issue at the `github page <https://github.com/jagophile/atomdb/issues>`_  


Contents
========

.. contents::
   :depth: 2
   :backlinks: entry
   
.. toctree:: 
   :maxdepth: 1
   
   apec
   atomic
   atomdb
   const
   spectrum
   util
   examples
   

=======
License
=======

Pyatomdb is released under the Smithsonian License:

Copyright 2015-16 Smithsonian Institution. Permission is granted to use, copy, 
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

--------
Examples
--------
Note: there are example routines demonstrating use of these features in the examples directory of the package.


------------
Installation
------------
PyAtomDB can be installed from pypi, using the simple ``pip install pyatomdb`` command.

For PyAtomDB to be useful, it requires access to a range of AtomDB database files (these are all `FITS <fits.gsfc.nasa.gov>`_ files). The database has two broad types of files, emissivity files and fundamental atomic data files (APED, the Astrophysical Plasma Emission Database). 

The emissivity files are needed for things such as producing spectra. The APED files are underlying atomic data and are not strictly needed for creating a spectrum, but can be useful for getting later information out.

In order for PyAtomDB to work efficiently, you should choose a location to store all of these files (e.g. /home/username/atomdb). It is strongly recommended that you set the environment variable ATOMDB to point to this, i.e. for bash add the following line to your .bashrc file::

  export ATOMDB=/home/username/atomdb
  
or for csh, add this to your .cshrc or .cshrc.login::
  
  setenv ATOMDB /home/username/atomdb   

If you run the following code, PyAtomDB will download the files you need to get started::

  import pyatomdb
  pyatomdb.util.initialize()

This will prompt you for an install location (defaulting to `$ATOMDB`)and whether to download the emissivity files. It is suggested that you say yes. It will also ask if you mind sharing anonymous download information with us. We would appreciate it if you say yes, but it is not necessary for the functioning of the software.

--------------------------
Example: Making a Spectrum
--------------------------
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






---------------------------------
Interrogating the atomic database
---------------------------------

The atomic database APED contains a range of data for a host of different ions. It contains a host of different files covering a range of different processes. The full database, when uncompressed is more than 10GB of data, so we are avoiding distributing it to all users. You can, however, get the individual data you need using the ``get_data`` routine::

  mydata = pyatomdb.atomdb.get_data(Z, z1, ftype)

This will try to open the file locally if it exists, and if it does not it will then go to the AtomDB FTP server and download the data for element Z, ion z1, with ftype a 2-character string denoting the type of data to get:

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

Downloaded data files are stored in ``$ATOMDB/APED/<elsymb>/<elsymb>_<ionnum>/``. You can delete them if you need to free up space, whenever a code needs the data it will reload them. There are many routines in the atomdb module which relate to extracting the data from the files, i.e. getting collisional excitation rates or line wavelengths. If you have trouble finding a routine to do what you want, please contact us and we'll be happy to write one if we can (this is how this module will grow - through user demand!)



   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

