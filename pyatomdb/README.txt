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
