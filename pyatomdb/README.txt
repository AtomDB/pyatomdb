========
PYATOMDB
========

This code is designed to interact with the AtomDB database. It is in active
development but is now approaching feature complete for creating spectra
and interacting with fitting codes like XSPEC and Sherpa.

Full documentation can be found at http://atomdb.readthedocs.io

============
Installation
============

Standard python installation:
python setup.py install


===============
Version History
===============

0.0.0.1
July 17th 2015: initial release

0.0.0.2
July 21st 2015: added ``dummyfirst`` keyword to
``pyatomdb.spectrum.make_spectrum``

0.0.0.3
August 12th 2015: fixed some errors in spectrum.py. Updated readme.

0.0.0.4
August 12th 2015: fixed example script in manual.

0.0.0.5
August 17th 2015: added PI cross section graphing data.
Updated get_data to fetch files from remote server

0.0.0.6
August 18th 2015: Fixed bug with documentation

0.0.1.0
September 23rd 2015: Added significant documentation. Converted all 'z0' to 'Z' and 'rmJ' to 'z1' throught the code for consistency.

0.0.1.1
September 25th 2015:   Major documentation revamp.
Separated datacache from settings keyword
Introduced automatic install and update into util package.

0.0.1.2
October 13th 2015:   More documentation improvements
Removed use of "curl" library, as apparently is non-standard. Replaced with 
wget, ftplib and requests.

0.0.1.3
November 5th 2015: Many minor bugfixes. Added "switch_version" to util package
to allow looking at older AtomDB version data. Added use option to get total RRC
emission using integral to infinity. 

0.0.1.5
February 12th 2016: As ever, bugfixes. Also added ability to calculate
ionization balance to apec module.

0.0.1.6
April 22nd 2016: Bugfixes! Also:
Major change to get_maxwell_rate interface
Added ability to apply RMF and ARF to spectra in spectrum module
Fixed hydrogenic RRC calculation

0.0.1.7
May 23rd 2016: Bugfix:update switch_version to get NEI files too.

0.0.1.8
October 07th 2016: Bugfixes! Also:
APEC now included in apec.py. Can create a full apec run from a par file.

0.0.1.9
October 13th 2016: Bugfixes! wget was causing problems during file download due to bugs, Now fixed

0.0.2.0
November 01st 2016: Removed use of sparse algebra. Wasn't accurate enough, causing issues with suprious lines in low ionization states.

Implemented faster ionization balance calculations in "apec.solve_ionbal_eigen"


0.0.2.1
November 16th 2016: Bugfix to the level population calculation for recombination. For ions with level resolved rates and PI cross sections, we were double counting. Fixed now.

0.0.2.2
November 30th 2016: Bugfix to download of new files: urllib.urlcleanup now called more often.

0.0.2.3
December 8th 2016: Updated to handle 3.0.7.

0.0.2.4
January 16th 2017: Bugfix for installation, include joblib as a requirement

0.0.2.5
January 17th 2017: Bugfix for installation, include mock as a requirement

0.0.3.0
February 21st 2017: Introduced new way to calculated spectra, using the Session
and Spec objects in spectrum.py. It is intended to make this the main way to
calculate spectra, as it should lead to significant speed enhancements by
avoiding lots of unnecessary recalculation of the same numbers.

0.0.3.1
February 21st 2017: A series of bug fixes to installation and spectral calculation routines. Closes bugs #1, 2, 5 from github

0.0.3.2
March 21st 2017: Several bug fixes and enhancements to the spectrum module. Updated setup.py to handle import of mock more gracefully

0.0.3.3
May 5th 2017: Several bug fixes and initial inclusion of charge exchange model

0.0.3.4
May 23rd 2017: Bug fix to spectrum.apply_response, now handles response files
indexed from zero (thank you Lia Corrales for finding this issue)

0.5.0
January 9th 2018: Major updates to include CX models. Increased version number to something relevant to 3 decimal style.

0.5.1
February 14th 2018: Bugfix to use 64bit floats for ionization fraction calculation, which was occasionally causing issues.

0.5.2
May 10th 2018: Bugfixes, including one where angstrom spectra where generated backwards.

0.5.3
June 14th 2018: Bugfixes: NEI spectrum session fixed to assign bin units correctly when using a response file

0.5.4
September 03rd 2018: Bugfix: NEI spectrum.Session.set_specbins now correctly resets the spectra when changed. Thanks to Shinya Nakashima for discovering this bug.


0.6.0
November 20th 2018: Switch to python 3 support only. Python 2.7 support available for now through Github

0.7.0
Bug fixes aplenty, mostly due to python 3 conversion.

0.7.1
April 22nd 2019: Fixed bug occuring with spectrum continuum broadening.

0.7.2
May 6th 2019: Add option 'ALL' to get_data to download all the data files for
an ion at once

0.8.0
January 27th 2020: Major overhaul of Spectrum module.
Updated paths for the FTP server due to file move.

0.8.1
February 5th 2020: Corrections to NEISpectrum. Spectrum summation much faster and
bugs in logic fixed.

0.9.0
March 20th 2020: Added PShockSession class to specturm module, mimicking XSPEC pshock (parallel shock) model. Other minor bugfixes.

0.10.0
May 12th 2020: Major update to documentation, spectrum and atomdb modules. Significant rationalization of routines, making many non-public in the documentation. This should help new users understand what is going on. I apologise if it breaks whatever you were doing, please get in touch if the fix is not obvious.

0.10.1
May 12th 2020: Minor packaging fixes

0.10.2
May 27th 2020: Bugfix to read_data

0.10.3
July 14th 2020: Added force keyword to switch_version. Amended make_linelist to handle non-51 temperature files

0.10.4
July 20th 2020: Fixed PShockSession initialization error. Added automatic removal of old pickle files generated from fits files when detected.

0.10.5
September 4th 2020: Fixed issue when handling temperatures above or below the maximum and minimum values in the APEC emissivity files. Added Electron-electron bremsstrahlung emission into APEC model.

0.10.6
September 9th 2020: Fixed missing filename in build - README.txt had been renamed README.rst. Also included pycurl in requirements list

0.10.7
October 12th 2020: Fixed a series of bugs - 28 and 29 on Github, and merged in ionization balance storage information pull request from Jim Slavin.

0.10.8
November 9th 2020: Minor changes for compatibility with changes to kappa module.

0.10.9
June 16th 2022: Several quality of life issues overdue for inclusion

- Spectrum session now has a verbosity switch
- Added more to the wrappers folder
- Added the ability to update wavelengths on the fly in spectrum objects
- Fixed errors when forcing use of older Mazzotta data for ionization

0.10.10
July 22nd 2022: Updated spectrum.py set_response to handle 0keV minimum energy responses

0.10.11
August 2nd 2022: Added "sparse" option to set_response, which uses sparse matrices in the RMF, allowing for large responses (e.g. XRISM) to be handled with limited memory usage.
