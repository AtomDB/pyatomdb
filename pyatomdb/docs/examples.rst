PyAtomDB Example Scripts
====================================

These are examples of using the pyatomdb module in your projects. They can all be found in the examples subdirectory of the pyatomdb tarball.

.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

====================
Initial installation
====================   
``first_installation.py``

.. literalinclude:: ../examples/first_installation.py

=====================
Examples with Spectra
=====================
The spectrum.py module contains routines for taking the spectrum generated from the
apec model and extracting line emissivities or continuum emission.

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

=====================
Make a Spectrum
=====================

Make a broadened and unbroadened spectrum:  ``make_spectrum.py``

.. literalinclude:: ../examples/make_spectrum.py



===========================
Make a Spectrum version 2.0
===========================

Make a spectrum using the new Session class:  ``new_make_spectrum.py``.
This is significantly faster if you need to make lots of spectra (fitting, 
interpolating between 2 temperatures etc). Note the example requires an
RMF and ARF file - adjust to fit your available response.

.. literalinclude:: ../examples/new_make_spectrum.py

==================
Make Cooling Curve
==================

Make a cooling curve, total emissivity in keV cm3 s-1, for each element
in a specfied spectral range (e.g. 2 to 10 keV).

.. literalinclude:: ../examples/calc_power.py
