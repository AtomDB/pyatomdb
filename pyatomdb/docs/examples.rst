PyAtomDB Example Scripts
====================================

These are examples of using the PyAtomDB module in your projects. They can all be found in the examples subdirectory of your PyAtomDB installation.

The examples here are separated by the module in which they are implemented. Generally speaking, the roles of these modules can be split into the following areas:

spectrum
  Obtaining emissivities from the AtomDB emissivity files (e.g. the ``apec`` model). This module is used to create spectra as required.

apec
  Related to the APEC code - creating the AtomDB emissivity files from the underlying atomic database

atomdb
  Accessing the atomic database - returning rates and coefficients, fetching files

atomic
  Basic atomic functions - getting atomic masses, converting to element symbols etc.

util
  File utilities (not generally interesting to end users)

const
  List of physical constants and code-related constants used throughout PyAtomDB (again, not generally interesting to end users)



.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry


=======
Spectra
=======
The spectrum.py module contains routines for taking the spectrum generated from the
apec model and extracting line emissivities or continuum emission, applying responses,
changing abundances, etc. In these examples, we will use the
`Chandra ACIS-S MEG +1 <https://cxc.cfa.harvard.edu/caldb/prop_plan/grating/index.html>`_
order grating responses as examples where one is required, but you can use any.

Similarly, we use ``pylab`` for plotting. You can of course use whatever system you
like for plotting, this is just what we did. It is included in matplotlib and
therefore hopefully present on most systems.

++++++++++++++++
CIESession Class
++++++++++++++++
The heart of the spectral analysis is the spectrum.py class. This reads in the results of an apec run (by default, $ATOMDB/apec_line.fits and $ATOMDB/apec_coco.fits) and allows the user to obtain spectra at a range of temperatures accounting for instrument responses, thermal and velocity broadening,  abundance changes and other issues.

-----------------
Making a Spectrum
-----------------

.. literalinclude:: ../examples/spectrum_session_examples.py

--------------------
Showing Line Details
--------------------

.. literalinclude:: ../examples/spectrum_session_linelist_examples.py

++++++++++++++++
NEISession Class
++++++++++++++++
Derived from the CIESession class, this handles non-equilibrium spectra. As such,
a few extra parameters should be set. The ionization timescale (tau) and the initial ionization fraction should be specified. This can either be as an initial temperature or an exact specified input distribution of ion populations.

-----------------
Making a spectrum
-----------------

.. literalinclude:: ../examples/spectrum_NEIsession_examples.py

--------------------
Showing Line Details
--------------------

.. literalinclude:: ../examples/spectrum_NEIsession_linelist_examples.py


==============
Make Line List
==============
List the strongest lines in a given temperature and wavelength region:  ``make_line_list.py``

.. literalinclude:: ../examples/make_line_list.py

=====================
Get PI Cross Sections
=====================

Extract the PI cross section data:  ``photoionization_data.py``

.. literalinclude:: ../examples/photoionization_data.py

==================
Make Cooling Curve
==================

Make a cooling curve, total emissivity in keV cm3 s-1, for each element
in a specfied spectral range (e.g. 2 to 10 keV).

.. literalinclude:: ../examples/calc_power.py
