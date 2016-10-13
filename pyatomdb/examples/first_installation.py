import pyatomdb

"""
This script shows the commands you should run when you first download pyatomdb.

It is recommended that you choose the location you want to install the
AtomDB data files (not the same as the python module) and set your
ATOMDB environment variable to point to it.


Parameters
------
none

Returns
-------
none

"""

# call the setup routine
pyatomdb.util.initialize()

# this routine downloads a bunch of files and sets things up for you. It will 
# take a few minutes, depending on your internet connection.

print "Install complete!"

#and that's it!

# If you want to switch versions of atomdb (in this case to 3.0.2) later, call:
# pyatomdb.util.switch_version('3.0.2')

