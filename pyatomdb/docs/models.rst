========================
PyAtomDB Spectral Models
========================

There are several different spectral models based on the same atomic data. In this section, we describe the physics and theories behind the models. Examples of usage are given on the :ref:`example-page-label` page.


----------------------------------------
Collisional Ionization Equilibrium (CIE)
----------------------------------------

This is the standard collisional ionization equilibrium model stored on 201
temperatures from 10\ :sup:`4` to 10\ :sup:`9` K. The model assumes an optically thin plasma

STILL UNDER CONSTRUCTION


----------------------------------------
Non-equilibrium Ionization (NEI)
----------------------------------------

STILL UNDER CONSTRUCTION


----------------------------------------
Plane-parallel shock (Pshock)
----------------------------------------

STILL UNDER CONSTRUCTION

.. _models-kappa:

--------------------------------------------
Kappa (non-Maxwellian Electron Distribution)
--------------------------------------------

~~~~~~~~~~~~
Introduction
~~~~~~~~~~~~

The AtomDB Kappa model is designed to create spectra for spectra with non-Maxwellian electron distributions of electron energies. The coefficients of Hahn and Savin [1]_ are used to create similar electron distributions from a sum of Maxwellian electron distributions. 

Given the available atomic data in AtomDB and other similar databases, it is important to ensure that the Maxwellian temperatures required to model the spectrum are present in the underlying atomic spectra database. Important work on this, including exploring the validity of the approximations is in [2]_. Note that there are no warnings included in this package for when you enter into areas where the model is not particularly valid.

~~~~~~~~~~~~
How it Works
~~~~~~~~~~~~

The AtomDB Kappa model uses the non-equilibrium emissivity data from the AtomDB project, combined with their ionization and recombination rates (mostly from [3]_) to determine the spectrum. When a spectrum is required, the following steps are carried out:

1. The Hahn and Savin coefficients are calculated, to seee what Maxwellian temperatures (Tm) are required and with what normalizations to create the kappa distribution

2. The ionization and recombination rates are calculated at each Tm, and then summed to provide kappa ionization and recombination rates. For Tm outside the nominal 10^4 to 10^9 Kelvin range of AtomDB data, the ionization and recombination data can be extrappolated fairly safely.

3. The ionization fraction for each ion is calculated using these kappa ionization and recombination rates

4. For each ion present in the resulting ionization fraction, the spectrum is calculated at each Tm and then summed appropriately to give a final specturm. As spectra cannot easily be extrappolated, and Tm outside 10^4 to 10^9 K are ignored. Again, this can lead to unrealistic results. 


.. [1] Hahn and Savin (2015) `2015ApJ...800...68H <https://ui.adsabs.harvard.edu/abs/2015ApJ...800...68H/abstract>`_
.. [2] Cui, Foster, Yuasa and Smith (2019), `2019ApJ...887..182C <https://ui.adsabs.harvard.edu/abs/2019ApJ...887..182C/abstract>`_
.. [3] Bryans (2009) `2009ApJ...691.1540B <https://ui.adsabs.harvard.edu/abs/2009ApJ...691.1540B/abstract>`_

There is a wrapper for XSPEC, kappa_xspec.py, in the wrappers folder. See this and :ref:`examples-kappa` for a quick examples of syntax and useage.
