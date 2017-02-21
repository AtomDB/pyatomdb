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
