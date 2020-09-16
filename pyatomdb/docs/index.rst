.. PyAtomDB documentation master file, created by
   sphinx-quickstart on Fri Sep 18 17:07:36 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================
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



========
Contents
========

.. contents::
   :depth: 2
   :backlinks: entry

.. toctree::
   :maxdepth: 4

   installation
   examples
   modules
   license
   contact



=======
Outline
=======


PyAtomDB is a selection of utilities designed to interact with the `AtomDB
database <https://www.atomdb.org>`_ . These utilities are under constant development. Please get in touch with any issues that arise.


There are several different modules currently. These are:

  :doc:`atomdb </atomdb>`
     A series of codes for interacting with the AtomDB atomic database
  :doc:`atomic </atomic>`
     Basic atomic data routines - e.g. converting element symbols to atomic number, etc.
  :doc:`const </const>`
     Physical and code related constants
  :doc:`spectrum </spectrum>`
     Routines for generating spectra from the published AtomDB line and continuum emissivity files
  :doc:`util </util>`
     Utility codes (sorting etc) that pyatomdb relies on.
  :doc:`apec </apec>`
     The full APEC code

To report bugs or make feature requests,  email the code authors or raise an issue at the `github page <https://github.com/jagophile/atomdb/issues>`_




==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

