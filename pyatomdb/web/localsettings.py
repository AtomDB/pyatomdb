import os

# This should be updated for your particular webserver/site
# DATAFILEPATH is where the pre-calculated spectra are stored 
# DEBUG is whether the server is in debug mode, though typically ignored
# ATOMDB is where the AtomDB data is stored

os.environ['DATAFILEPATH'] = '/var/www/atomdb/apps/newapp/data'
os.environ['DEBUG']='1'
os.environ['ATOMDB'] = '/var/www/atomdb/atomdb_release'
