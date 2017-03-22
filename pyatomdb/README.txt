========
PYATOMDB
========

This is pre-alpha release code. Many functions have just been ripped from my
hard disk and not yet fully vetted. Possibly the most complete section is the 
"spectrum" module. Sample code is in the docs subfolder.

Full documentation can be found at http://atomdb.readthedocs.org

=======================
Installation
=======================
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
