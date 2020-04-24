.. PyAtomDB documentation master file, created by
   sphinx-quickstart on Fri Sep 18 17:07:36 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyAtomDB & AtomDB
=================
The `AtomDB Project <https://www.atomdb.org>`_ consists of a large atomic database
designed for creating spectra of collisionally excited plasmas in the UV and
X-ray wavebands for use in astronomy and astrophysics research.
It has been successfully used for some time within several
different spectral analysis suites, such as `XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_
and `Sherpa <https://cxc.cfa.harvard.edu/sherpa/>`_.

There are two main parts to AtomDB:

The Astrophysical Plasma Emission Database (APED)
  A series of `FITS <https://heasarc.gsfc.nasa.gov/docs/heasarc/fits.html>`_
  files which store the atomic data necessary for modelling emission from collisional plasmas
  from elements from H through to Zn.

The Astrophysical Plasma Emission Code (APEC)
  Code which takes the APED and uses it to model a collisional
  plasma, creating emissivity files for line and continuum emission. These output files are
  then used in a range of models such as the ``apec``, ``nei`` and ``pshock`` models amongst others.

Starting in 2015, PyAtomDB was developed to achieve a number of goals:

  #. To replace the APEC code (formerly in C) with a more flexible tool for larger data sets
  #. To allow interactive user access to the underlying database, which was always freely
     available but sometimes difficult to use
  #. To enable creating more flexible spectra based on the outputs of AtomDB
  #. To enhance the types of spectra which can be modeled in AtomDB, and astronomer interaction
     with these models
  #. To facililtate inclusion of AtomDB models in other software (often python based).

Contact
=======
This module is still in very active development. If you have feature requests, bug
reports or any other questions, please contact us through the project's
`GitHub <https://github.com/AtomDB/pyatomdb>`_ page.


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
   license




=======
Outline
=======


PyAtomDB is a selection of utilities designed to interact with the `AtomDB
database <https://www.atomdb.org>`_ . These utilities are under constant development. Please get in touch with any issues that arise.


There are several different modules currently. These are:

- :doc:`atomdb </atomdb>`  : a series of codes for interacting with the AtomDB atomic database
- :doc:`atomic </atomic>`  : basic atomic data routines - e.g. converting element symbols to atomic number, etc.
- :doc:`const </const>`   : a series of physical constants
- :doc:`spectrum </spectrum>` : routines for generating spectra from the published AtomDB line and continuum emissivity files
- :doc:`util </util>`    : sumple utility codes (sorting etc) that pyatomdb relies on.
- :doc:`apec </apec>`  : ultimately, the full apec code. For now, incomplete.

Expect bugs. Report those bugs! Make feature requests! Email the code authors or raise an issue at the `github page <https://github.com/jagophile/atomdb/issues>`_



============
Installation
============

.. warning::
  PyAtomDB runs only under Python 3. It will not work on Python 2. If you
  need to install Python 3 in your system, you can use your package manager
  or there are many other sources which
  can help you including `Anaconda <https://www.anaconda.com/>`_.

  Once you have Python 3 installed, you may (depending on your system)
  have to add a ``3`` to many
  of the command line commands, e.g. ``python`` becomes ``python3``, or
  ``pip`` becomes ``pip3``. Commands in this guide omit the ``3``. Note
  that actual Python code is unaffected.



PyAtomDB can be installed in two ways:

  #. From `PyPI <https://pypi.org/>`_ , using the simple ``pip install pyatomdb`` command.
  #. From `GitHub`_ , using the command ``git clone https://github.com/AtomDB/pyatomdb.git``
     to get the source, then ``python setup.py develop`` to install links to the source
     in your Python path.

Note that for both of these options the ``--user`` option can be useful, as it will install
software in your local path if you do not have administrator priviledges on your machine.

You can check that the installation was successful by running:

.. code-block:: python

  >>> import pyatomdb

If it does not immediately throw out an error, it has been successful. It will
then start asking about installing the AtomDB files, see the next section. Note
that there is no longer a need to run the initialize script.

----------------
ATOMDB Directory
----------------
Whenever you import the PyAtomDB module, it performs a check for the $ATOMDB directory.
This directory is where the AtomDB data files will be stored. These are not
distributed with the python package as they are large and most people will only need
a few. PyAtomDB will download the APED data files on demand as you require them, and
they are then stored in this directory until you manually delete them. If you need to
recover disk space, you can delete anything in the $ATOMDB/APED directory without
repercussions - PyAtomDB will re-download the files if it needs them in the future.

You will be asked to select a directory for installation and then whether to download
the emissivity files. It is important that this directory is one where you have write
access, as in the future further files will be added there by the code automatically.

Once installation is complete, ensure that you add the ATOMDB variable to your
shell startup file. Assuming you have installed into /home/username/atomdb,
in bash, add this line to your ~/.bashrc or ~/.bash_profile files
(depends on which one is sourced by your system):

.. code-block:: bash

  export ATOMDB=/home/username/atomdb

or for csh, add this to your ~/.cshrc or ~/.cshrc.login:

.. code-block:: csh

  setenv ATOMDB /home/username/atomdb

Recent version of Mac OS have moved to zsh, in which case modify your ~/.zshrc file as for bash above.

---------------------
Usage Data Collection
---------------------

You will also be asked about anonymous usage data. In order to track roughly how many
people are using PyAtomDB, a randomly generated number is created when you install
PyAtomDB and stored in your ``$ATOMDB/userdata`` file. Whenever PyAtomDB has to fetch
a new file this number, the filename and the current timestamp is stored on our
system so we can estimate how many users there are. We have no way to connect this
to actual individuals, it simply tells us roughly how many unique active users
there are.

If you decline, this number is set to 00000000, and otherwise PyAtomDB functions
as normal.




Note that it requires python 3 to run. On some systems, this will require calling ``pip3`` instead of ``pip``.
In all the examples in this guide, we assume python 3 is running. Again,
depending on your system this can be invoked by either ``python`` or ``python3``.

For PyAtomDB to be useful, it requires access to a range of AtomDB database files. The database has two broad types of files, emissivity files (APEC) and fundamental atomic data files (APED, the Astrophysical Plasma Emission Database).

The emissivity files are needed for things such as producing spectra. The APED files are underlying atomic data and are not strictly needed for creating a spectrum, but can be useful for getting later information out.

In order for PyAtomDB to work efficiently, you should choose a location to store all of these files (e.g. /home/username/atomdb). It is strongly recommended that you set the environment variable ATOMDB to point to this, i.e. for bash add the following line to your .bashrc file::

  export ATOMDB=/home/username/atomdb

or for csh, add this to your .cshrc or .cshrc.login::

  setenv ATOMDB /home/username/atomdb

If you run the following code within a python shell, PyAtomDB will download the files you need to get started::

  import pyatomdb
  pyatomdb.util.initialize()

This will prompt you for an install location (defaulting to `$ATOMDB`) and whether to download the emissivity files. It is suggested that you say yes. It will also ask if you mind sharing anonymous download information with us. We would appreciate it if you say yes, but it is not necessary for the functioning of the software.





--------
Examples
--------
Note: there are example routines demonstrating use of these features in the examples directory of the package.


------------
Installation
------------

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

