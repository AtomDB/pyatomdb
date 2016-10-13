"""
util.py contains a range of miscellaneous helper codes that assist in running
other AtomDB codes but are not in any way part of a physical calculation.

Version -.1 - initial release
Adam Foster July 17th 2015
"""

import numpy, os, errno, hashlib
import requests, urllib, time, subprocess, shutil, wget, glob
import datetime
import const, atomic, atomdb
from StringIO import StringIO
import ftplib
try:
  import astropy.io.fits as pyfits
except:
  import pyfits
  
  

################################################################################
#
#  Python Module
#
#  Name:        util.py
#
#  Decription:  Python codes that I use.
#
#  Module contents (and 1 line description: see individual modules for more):
#
#     differentiate
#          Differentiate line.
#
#  First Version:
#       Adam Foster, 07-Jun-2010
#
################################################################################

#-------------------------------------------------------------------------------

def figcoords(lowxpix, lowypix, highxpix, highypix, 
                  lowxval, lowyval, highxval, highyval, 
                  xpix, ypix, logx=False, logy=False):



  if logx:
    lxval = numpy.log10(lowxval)
    hxval = numpy.log10(highxval)
  else:
    lxval = lowxval*1.0
    hxval = highxval*1.0
  if logy:
    lyval = numpy.log10(lowyval)
    hyval = numpy.log10(highyval)
  else:
    lyval = lowyval*1.0
    hyval = highyval*1.0
  
  dx = (hxval-lxval)/(highxpix-lowxpix)
  dy = (hyval-lyval)/(highypix-lowypix)
  
#  print dx, dy
  xout = lxval + (xpix-lowxpix)*dx
  yout = lyval + (ypix-lowypix)*dy

  if logx:
    xout = 10**xout
  if logy:
    yout = 10**yout
  return xout, yout

#-------------------------------------------------------------------------------

def unique(s):
     """Return a list of the elements in s, but without duplicates.

     For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
     unique("abcabc") some permutation of ["a", "b", "c"], and
     unique(([1, 2], [2, 3], [1, 2])) some permutation of
     [[2, 3], [1, 2]].

     For best speed, all sequence elements should be hashable.  Then
     unique() will usually work in linear time.

     If not possible, the sequence elements should enjoy a total
     ordering, and if list(s).sort() doesn't raise TypeError it's
     assumed that they do enjoy a total ordering.  Then unique() will
     usually work in O(N*log2(N)) time.

     If that's not possible either, the sequence elements must support
     equality-testing.  Then unique() will usually work in quadratic
     time.
     
     Parameters
     ----------
     s : list type object
       List to remove the duplicates from
     
     Returns
     -------
     list type object
       ...with all the duplicates removed
     
     References
     ----------
     Taken from Python Cookbook, written by Tim Peters.
     http://code.activestate.com/recipes/52560/
       
     
     """

     n = len(s)
     if n == 0:
         return []

     # Try using a dict first, as that's the fastest and will usually
     # work.  If it doesn't work, it will usually fail quickly, so it
     # usually doesn't cost much to *try* it.  It requires that all the
     # sequence elements be hashable, and support equality comparison.
     u = {}
     try:
         for x in s:
             u[x] = 1
     except TypeError:
         del u  # move on to the next method
     else:
         return u.keys()

     # We can't hash all the elements.  Second fastest is to sort,
     # which brings the equal elements together; then duplicates are
     # easy to weed out in a single pass.
     # NOTE:  Python's list.sort() was designed to be efficient in the
     # presence of many duplicate elements.  This isn't true of all
     # sort functions in all languages or libraries, so this approach
     # is more effective in Python than it may be elsewhere.
     try:
         t = list(s)
         t.sort()
     except TypeError:
         del t  # move on to the next method
     else:
         assert n > 0
         last = t[0]
         lasti = i = 1
         while i < n:
             if t[i] != last:
                 t[lasti] = last = t[i]
                 lasti += 1
             i += 1
         return t[:lasti]

     # Brute force is all that's left.
     u = []
     for x in s:
         if x not in u:
             u.append(x)
     return u
     
#-------------------------------------------------------------------------------

def mkdir_p(path):
  """
  Create a directory. If it already exists, do nothing.
  
  Parameters
  ----------
  path : string
    The directory to make
  
  Returns
  -------
  none
  """
  
  try:
    os.makedirs(path)
  except OSError as exc: # Python >2.5
    if exc.errno == errno.EEXIST:
      pass
    else: raise

#-------------------------------------------------------------------------------

def question(question, default, multichoice  = []):
  """
  Ask question with default answer provided. Return answer
  
  Parameters
  ----------
  question : str
    Question to ask
  default : str
    Default answer to question
  multichoice : str
    if set, answer must be one of these choices
  
  Returns
  -------
  str
    The answer.
  """

#  Version 0.1 - initial release
#  Adam Foster September 23rd 2015

  if len(multichoice) > 0:
    mclower = []
    for i in multichoice:
      mclower.append(i.lower())
    qstring = question+' ('
    for i in multichoice:
      qstring+=i+'/'
    qstring=qstring[:-1]
    qstring+=') [%s]'%(default)
    
    while True:
      ans = raw_input(qstring)
      if ans=='':
        ans=default
      if ans.lower() in mclower:  
        break
  else:      
    ans = raw_input("%s [%s] :" %(question, default))
    if ans =='':
      ans = default
  return ans

#-------------------------------------------------------------------------------

def record_upload(fname):
  """
  Transmits record of a file transfer to AtomDB
  
  This simply transmits the USERID, filename, and time to AtomDB. If USERID=0,
  then the user has chosen not to share this information and this is skipped
  
  Parameters
  ----------
  fname : string
    The file name being downloaded.
  
  Returns
  -------
  None
  
  """
  #
  # Version 0.1 Initial Release
  # Adam Foster 24th September 2015
  #
  try:
    tmp = load_user_prefs()
    userid = tmp['USERID']
  except KeyError:
    print load_user_prefs()
  
  if int(userid) > 0:
    postform = {  'TYPE' : const.FILEDOWNLOAD,\
                  'USERID': int(userid),\
                  'FILE': fname,\
                  'TIME': time.time()}
    r = requests.post('http://www.atomdb.org/util/process_downloads.php',\
                      data = postform)
  return
    
#-------------------------------------------------------------------------------

def md5Checksum(filePath):
  """
  Calculate the md5 checksum of a file
  
  Parameters
  ----------
  filepath : str
    the file to calculate the md5sum of
  
  Returns
  -------
  string
    the hexadecimal string md5 hash of the file
    
  References
  ----------
  Taken from http://joelverhagen.com/blog/2011/02/md5-hash-of-file-in-python/
  
  """
    
  with open(filePath, 'rb') as fh:
    m = hashlib.md5()
    while True:
      data = fh.read(8192)
      if not data:
        break
      m.update(data)
  return m.hexdigest()

#-------------------------------------------------------------------------------

def download_atomdb_emissivity_files(adbroot, userid, version):
  """
  Download the AtomDB equilibrium emissivity files for AtomDB"
  
  This code will go to the AtomDB FTP site and download the necessary files.
  It will then unpack them into a directory adbroot. It will not
  overwrite existing files with the same md5sum (to avoid pointless updates)
  but it will not know this until it has downloaded and unzipped the main 
  file.
  
  Parameters
  ----------
  
  adbroot : string
    The location to install the data. Typically should match $ATOMDB
  userid : string
    An 8 digit ID number. Usually passed as a string, but integer
    is also fine (provided it is all numbers)
  version : string
    The version string for the release, e.g. "3.0.2"
  
  Returns
  -------
  None
  """
  #
  # Version 0.1 Initail Release
  # Adam Foster 24th September 2015
  #
   
  # set up remote file name
  fname = "atomdb_v%s.tar.bz2"%(version)
  
  # set up temporary directory to hold data
  mkdir_p(adbroot)
  if adbroot[-1]=='/':
    tmpdir = adbroot+'installtmp'
  else:
    tmpdir = adbroot+'/installtmp'
  
  mkdir_p(tmpdir)
  
  # get the files
  fnameout = wget.download('ftp://sao-ftp.harvard.edu/AtomDB/releases/%s'%(fname), out="%s/%s"%(tmpdir, fname))

  # collect user statistics if allowed.
  record_upload(fname)
    
  #uncompress
  print "\nUncompressing",
  cwd=os.getcwd()
  os.chdir(tmpdir)
  subprocess.call(["tar", "-xjf", "%s"%(fnameout)])
  print "...done"
  
  print "Moving files to %s" % (adbroot),
  for ifile in glob.glob('./*/*'):
    outfile = adbroot+'/'+ifile.split('/')[-1]
    if os.path.exists(outfile):
      try:
        if md5Checksum(outfile) == md5Checksum(ifile):
          
        # these files are the same, don't bother copying or 
        # asking about copying them.
          continue
      except IOError:
        print "outfile = %s, ifile = %s"%(outfile, ifile)
        raise
      
      overwrite = question("file %s already exists. Overwrite?"%(outfile),"y",["y","n"])
      if overwrite:
        os.remove(outfile)
      else:
        continue
    shutil.move(ifile,outfile)


  os.chdir(cwd)
  shutil.rmtree(tmpdir)
    
  print "...done"
    
#-------------------------------------------------------------------------------
    
def download_atomdb_nei_emissivity_files(adbroot, userid, version):
  """
  Download the AtomDB non-equilibrium emissivity files for AtomDB"
  
  This code will go to the AtomDB FTP site and download the necessary files.
  It will then unpack them into a directory adbroot. It will not
  overwrite existing files with the same md5sum (to avoid pointless updates)
  but it will not know this until it has downloaded and unzipped the main 
  file.
  
  Parameters
  ----------
  
  adbroot : string
    The location to install the data. Typically should match $ATOMDB
  userid : string
    An 8 digit ID number. Usually passed as a string, but integer
    is also fine (provided it is all numbers)
  version : string
    The version string for the release, e.g. "3.0.2"
  
  Returns
  -------
  None
  """
  #
  # Version 0.1 Initail Release
  # Adam Foster 24th September 2015
  #
   
  # set up remote file name
  fname = "atomdb_v%s_nei.tar.bz2"%(version)
  
  # set up temporary directory to hold data
  mkdir_p(adbroot)
  if adbroot[-1]=='/':
    tmpdir = adbroot+'installtmp'
  else:
    tmpdir = adbroot+'/installtmp'
  
  mkdir_p(tmpdir)
  
  # get the files
  fnameout = wget.download('ftp://sao-ftp.harvard.edu/AtomDB/releases/%s'%(fname), out="%s/%s"%(tmpdir, fname))

  # collect user statistics if allowed.
  record_upload(fname)
    
  #uncompress
  print "\nUncompressing",
  cwd=os.getcwd()
  os.chdir(tmpdir)
  subprocess.call(["tar", "-xjf", "%s"%(fnameout)])
  print "... done"
  
  print "Moving files to %s" % (adbroot),
  for ifile in glob.glob('./*/*'):
    outfile = adbroot+'/'+ifile.split('/')[-1]
    if os.path.exists(outfile):
      try:
        if md5Checksum(outfile) == md5Checksum(ifile):
          
        # these files are the same, don't bother copying or 
        # asking about copying them.
          continue
      except IOError:
        print "outfile = %s, ifile = %s"%(outfile, ifile)
        raise
      
      overwrite = question("file %s already exists. Overwrite?"%(outfile),"y",["y","n"])
      if overwrite:
        os.remove(outfile)
      else:
        continue
    shutil.move(ifile,outfile)

  os.chdir(cwd)
  shutil.rmtree(tmpdir)
    
  print "... done"
  
#-------------------------------------------------------------------------------

def load_user_prefs(adbroot="$ATOMDB"):
  """
  Loads user preference data from $ATOMDB/userdata
  
  Parameters
  ----------
  adbroot : string
    The AtomDB root directory. Defaults to environment variable $ATOMDB.
  
  Returns
  -------
  dictionary
    keyword/setting pairs e.g. settings['USERID'] = "12345678"
  
  """
  #
  # Version 0.1 Adam Foster
  # Initial Release 24th Sep 2015
  # 
    
  ret = {}

  fname = os.path.expandvars(adbroot+'/userdata')
  with open(fname,'r') as f:
    for l in f:
      #remove trailing/leading whitespace
      ls = l.strip()
      # remove everything after a "#"
      ls = ls.split('#')[0]
      if len(ls) > 0:
        ls = ls.split('=')
        if len(ls) != 2:
          print "Error: cannot parse line in userdata file: %s\n%s"%(fname, l)
        else:
          ret[ls[0].strip()]=ls[1].strip()
  return ret    

#-------------------------------------------------------------------------------

def write_user_prefs(prefs, adbroot="$ATOMDB"):
    """
    Write user preference data to $ATOMDB/userdata. This will overwrite the 
    entire file. 
    
    Therefore you should use "load_user_prefs", then add in additional 
    keywords, the call write_user_prefs.
    
    Parameters
    ----------
    prefs: dictionary
      keyword/setting pairs e.g. settings['USERID'] = "12345678"
    
    adbroot : string
      The AtomDB root directory. Defaults to environment variable $ATOMDB.
    
    Returns
    -------
    None
    """
    #
    # Version 0.1 Adam Foster
    # Initial Release 24th Sep 2015
    # 
    
    fname = os.path.expandvars(adbroot+'/userdata')
    
    l =  open(fname,'w')
      #remove trailing/leading whitespace
    for i in prefs.keys():
      l.write("%s = %s\n" %(i,prefs[i]))
    l.close()
    return

#-------------------------------------------------------------------------------

def initialize():
  """
  Initialize your AtomDB Setup
  
  This code will let you select where to install AtomDB, get the
  latest version of the filemap, and download the emissivity
  files needed for various functions to work.
  
  Parameters
  ----------
  None. 
  
  Returns
  -------
  None
  
  """
  
  if 'ATOMDB' in os.environ:
    adbroot_init = os.environ['ATOMDB']
  else:
    adbroot_init = os.path.expandvars("$HOME/atomdb")
  adbroot = question("Location to install AtomDB files",\
                          adbroot_init)

  anondat = ''
  while not anondat in ['yes','no']:
    anondat = question("Allow reporting of anonymous usage statistics","yes",multichoice=["yes","no","info"])
    if anondat=='info':
      print "We like to know how many of users use our data. To enable this, we will generate a random number which will be transmitted with your requests to download data files in the future. We will record and transmit no other data beyond this number and what files you are downloading, and we will have no way of connecting this number to you. If you decline, this number will be set to 0"
  if anondat == 'no':
    userid = '00000000'
  else:
    userid = "%08i"%(numpy.random.randint(1e8))
  
  print "You are about to install:"
  print "AtomDB installation location: %s"%(adbroot)
  print "Transmit anonymous user data: %s"%(anondat)
  print "Randomly generated userid: %s"%(userid)
  
  proceed = question("Is this ok?","y",\
                          multichoice=["y","n"])
                          
  if proceed == 'n':
    print "Aborting"
    return
  
  if proceed == 'y':
    
    print "Temporarily setting $ATOMDB environment variable to %s."%(adbroot),
    print " It is *strongly* recommended that you add this to your environment ",
    print " permanently."
    os.environ['ATOMDB'] = adbroot
    
    print "creating directory %s"%(adbroot),
    mkdir_p(adbroot)
    print "...done"
    
    print "creating user data file %s/userdata"%(adbroot),
    userdatafname = "%s/userdata"%(adbroot)
    if os.path.exists(userdatafname):
      print "... user data file already exists. Will update, not overwrite",
      userprefs = load_user_prefs(adbroot=adbroot)
    else:
      userprefs={}
      userprefs['USERID'] = userid
      userprefs['LASTVERSIONCHECK'] = time.time()
    write_user_prefs(userprefs, adbroot=adbroot)

    print "...done"
    
    print "finding current version of AtomDB. ",
#    a=curl.Curl()
#    version=a.get('ftp://sao-ftp.harvard.edu/AtomDB/releases/LATEST')
#    a.close()
    
    ftp = ftplib.FTP('sao-ftp.harvard.edu') 
    x = ftp.login()
    r = StringIO()
    x = ftp.retrbinary('RETR /AtomDB/releases/LATEST', r.write)
    version = r.getvalue()[:-1]
    x = ftp.quit()
    print "Latest version is %s"%(version)
    
    get_new_files=question(\
      "Do you wish to download the emissivity data for these files (recommended)?",\
      "y",multichoice=["y","n"])
    
    if get_new_files=='y':
      download_atomdb_emissivity_files(adbroot, userid, version)

    get_new_nei_files=question(\
      "Do you wish to download the non-equilibrium emissivity data for these files (recommended)?",\
      "y",multichoice=["y","n"])
    
    if get_new_nei_files=='y':
      download_atomdb_nei_emissivity_files(adbroot, userid, version)

#-------------------------------------------------------------------------------

def check_version():
  """
  Checks if there is a more recent version of the database to install.
  
  Parameters
  ----------
  None. 
  
  Returns
  -------
  None
  
  """
  
  try:
    adbroot_init = os.environ['ATOMDB']
  except keyError:
    print "You must set the ATOMDB environment variable for this to work!"
    raise
  
  ftp = ftplib.FTP('sao-ftp.harvard.edu') 
  x = ftp.login()
  r = StringIO()
  x = ftp.retrbinary('RETR /AtomDB/releases/LATEST', r.write)
  newversion = r.getvalue()[:-1]
  x = ftp.quit()

  curversion = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]

  if (curversion != newversion):
    ans = question("New version %s is available. Upgrade?"%(newversion),"y",["y","n"])
  
    if ans=="y":
      get_new_files=question(\
        "Do you wish to download the emissivity data for these files (recommended)?",\
        "y",multichoice=["y","n"])
    
      if get_new_files=='y':
        download_atomdb_emissivity_files(adbroot, userid, version)

      get_new_nei_files=question(\
        "Do you wish to download the non-equilibrium emissivity data for these files (recommended)?",\
        "y",multichoice=["y","n"])
    
      if get_new_nei_files=='y':
        download_atomdb_nei_emissivity_files(adbroot, userid, version)
  else:
    print "Current version %s is up to date" %(curversion)

  # now update the time the last version check happened.
  userprefs = load_user_prefs()
  userprefs['LASTVERSIONCHECK'] = time.time()
  write_user_prefs(userprefs)

#-------------------------------------------------------------------------------

def switch_version(version):
  """
  Changes the AtomDB version. Note this will overwrite several links 
  on your hard disk, and will *NOT* be repaired upon quitting python.

  The files affect are the VERSION file and the soft links
  $ATOMDB/apec_line.fits, $ATOMDB/apec_coco.fits, $ATOMDB/filemap and
  $ATOMDB/apec_linelist.fits
  
  Parameters
  ----------
  version: string
    The version of AtomDB to switch to. Should be of the form "2.0.2"

  Returns
  -------
  None

  """
  #
  # Intial version November 4th 2015
  # Adam Foster
  #
  import re
  try:
    adbroot_init = os.environ['ATOMDB']
  except keyError:
    print "You must set the ATOMDB environment variable for this to work!"
    raise
  
  
  # check the AtomDB version string is a suitable string
  
  if not re.match('^\d\.\d\.\d$',version):
    print "Error: version number must be of format %i.%i.%i, e.g. 3.0.2"
    return

  # check current version
  curversion = open(os.path.expandvars('$ATOMDB/VERSION'),'r').read()[:-1]
  
  if curversion == version:
    print "Already using version %s. Not changing anything!" %(version)
    return
  
  # ok, otherwise we must do things!
  startdir = os.getcwd()

  # check for existing local files. If so we can just change the pointers
  mustdownload = False  
  
  flist =   ['apec_vVERSION_line.fits', 'apec_vVERSION_coco.fits',\
                'filemap_vVERSION', 'apec_vVERSION_linelist.fits']
  if version[0] == '3':
    flist.append('apec_vVERSION_nei_line.fits')
    flist.append('apec_vVERSION_nei_comp.fits')
  
  for f in flist:
    fname_test = re.sub('VERSION',version, os.path.expandvars("$ATOMDB/%s"%(f)))
    
    if os.path.exists(fname_test):
      print "%s already exists"%(fname_test)
      pass
      
    else:
      print "We are missing some files for this version, downloading now"
      mustdownload = True

  if mustdownload:
    # go find the files      
    ftproot = 'sao-ftp.harvard.edu'
    # get the filename
    if version[0] =='2':
      fname = re.sub('VERSION',version,'atomdb_vVERSION_runs.tar.gz')
      dirname = 'AtomDB/releases/old'
    elif version[0] =='3':
      fname = re.sub('VERSION',version,'atomdb_vVERSION.tar.bz2')
      dirname = 'AtomDB/releases'
      
    localfile = os.path.expandvars("$ATOMDB/tmp/%s"%(fname))
      # create temporary folder
    mkdir_p(os.path.expandvars("$ATOMDB/tmp"))
      
     
    print "Attempting to download %s to %s"%('ftp://%s/%s/%s'%(ftproot,dirname,fname), localfile)
    
    if os.path.isfile(localfile):
      os.remove(localfile)
    # get the file
    try:
      wget.download('ftp://%s/%s/%s'%(ftproot,dirname,fname), localfile)
    except IOError:
      print "Cannot find file ftp://%s/%s/%s on server. Please check that version %s is a valid version."%(ftproot,dirname,fname, version)
      return    
    # ok, now open up the relevant file and copy the things we need
    print "\nUncompressing %s..." %(localfile)

    cwd=os.getcwd()
    os.chdir(os.path.expandvars("$ATOMDB/tmp"))
    if fname.split('.')[-1] == 'gz':
      subprocess.call(["tar", "-xvzf", "%s"%(fname)])
    elif fname.split('.')[-1] == 'bz2':
      subprocess.call(["tar", "-xvjf", "%s"%(fname)])

    print "...done"
    urllib.urlcleanup()

    if version[0] =='3':
      fname = re.sub('VERSION',version,'atomdb_vVERSION_nei.tar.bz2')
      dirname = 'AtomDB/releases'
      
      localfile = os.path.expandvars("$ATOMDB/tmp/%s"%(fname))
        # create temporary folder
      mkdir_p(os.path.expandvars("$ATOMDB/tmp"))
      
     
      print "Attempting to download %s to %s"%('ftp://%s/%s/%s'%(ftproot,dirname,fname), localfile)
    # get the file
      if os.path.isfile(localfile):
        os.remove(localfile)

      try:
        wget.download('ftp://%s/%s/%s'%(ftproot,dirname,fname), localfile)
      except IOError:
        print "Cannot find file ftp://%s/%s/%s on server. Please check that version %s is a valid version."%(ftproot,dirname,fname, version)
        return    
    # ok, now open up the relevant file and copy the things we need
      print "\nUncompressing %s..." %(localfile)

      cwd=os.getcwd()
      os.chdir(os.path.expandvars("$ATOMDB/tmp"))
      if fname.split('.')[-1] == 'gz':
        subprocess.call(["tar", "-xvzf", "%s"%(fname)])
      elif fname.split('.')[-1] == 'bz2':
        subprocess.call(["tar", "-xvjf", "%s"%(fname)])

    os.chdir(os.path.expandvars(re.sub('VERSION',version, "$ATOMDB/tmp/atomdb_vVERSION")))

    urllib.urlcleanup()

    for ifile in glob.glob('*%s*'%(version)):
      
      outfile = os.path.expandvars("$ATOMDB/%s"%(ifile))
      if os.path.exists(outfile):
        try:
          if md5Checksum(outfile) == md5Checksum(ifile):
            print "file %s already exists, not overwriting"%(ifile)
        # these files are the same, don't bother copying or 
        # asking about copying them.
            continue
        except IOError:
          print "outfile = %s, ifile = %s"%(outfile, ifile)
          raise
      
        overwrite = question("file %s already exists. Overwrite?"%(outfile),"y",["y","n"])
        if overwrite:
          os.remove(outfile)
          shutil.move(ifile,outfile)
        else:
          continue
      else: 

        shutil.move(ifile,outfile)

    print "...done"
    
    os.chdir(os.path.expandvars("$ATOMDB"))
    # delete temporary directory
    shutil.rmtree('tmp')
    
  os.chdir(os.path.expandvars("$ATOMDB"))
  # OK, download complete. Now to make symlinks
  
  flistlist = ['apec_vVERSION_line.fits', 'apec_vVERSION_coco.fits',\
                'filemap_vVERSION', 'apec_vVERSION_linelist.fits']
  if int(version[0]) >=3:
    flistlist.append('apec_vVERSION_nei_line.fits')
    flistlist.append('apec_vVERSION_nei_comp.fits')
    
  for flist in flistlist:
                  
    # remove existing link if there is one
#    print "CHECKING FOR %s/%s"%(os.getcwd(), re.sub('VERSION',version,flist))
#    if os.path.exists(re.sub('VERSION',version,flist)):
#      os.remove(re.sub('VERSION',version,flist))
#      print "REMOVING LINK1!"
    # add new link
    print flist
    if os.path.islink(re.sub('_vVERSION','',flist)):
      os.remove(re.sub('_vVERSION','',flist))
      
    os.symlink(re.sub('VERSION',version,flist), re.sub('_vVERSION','',flist))
    
    # update version
    a = open('VERSION','w')
    a.write(version+'\n')
    a.close()
  os.chdir(startdir)

#-------------------------------------------------------------------------------

def make_vec(d):
  """
  Create vector version of d, return True or false depending on whether
  input was vector or not

  Parameters
  ----------
  d: any scalar or vector
    The input

  Returns
  -------
  vecd : array of floats
    d as a vector (same as input if already an iterable type)
  isvec : bool
    True if d was a vector, otherwise False.
  

  """
  isvec=True
  try:
    _ = (e for e in d)
  except TypeError:
    isvec=False
    d= numpy.array([d])
    
  return d,isvec

#-------------------------------------------------------------------------------

def write_lv_file(fname, dat, clobber=False):
  """
  Write the data in list dat to fname

  Parameters
  ----------
  fname : string
    The file to write
  dat : list
    The data to write
    Should be a list with the following keywords:
    Z : int
       nuclear charge
    z1 : int
       ion charge + 1
    comments : iterable of strings
       comments to append to the file
    data : numpy.array
      stores all the individual level data, with the following types:
      elec_config : string (40 char max) : Electron configuration strings\n
      energy : float: Level energy (eV)\n
      e_error : float : Energy level error (eV)\n
      n_quan : int : N quantum number\n
      l_quan : int : L quantum number\n
      s_quan : float : S quantum number\n
      lev_deg : int : level degeneracy\n
      phot_type : int : photoionization data type\n
          -1 : none\n
          0  : hydrogenic\n
          1  : Clark\n
          2  : Verner\n
          3  : XSTAR data
      phot_par : float(20) : photoionization paramters (see specific PI type for definition)\n
      Aaut_tot : float (optional) : the total autoionization rate out of the level (s^-1)\n
      Arad_tot : float (optional) : the total radiative rate out of the level (s^-1)\n
      energy_ref : string(20) :  energy reference (usually bibcode)\n
      phot_ref : string(20) : photoionization reference (bibcode)\n
      Aaut_ref : string(20) : total autoionization rate reference (bibcode)\n
      Arad_ref : string(20) : total radiative decay rate reference (bibcode)\n
  clobber : bool
    Overwrite existing file if it exists.
    
  Returns
  -------
  none
  

  """

  # Create primary HDU (it's a dummy one)
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  if 'aaut_tot' in dat['data'].dtype.names:
    fileversion = '1.1.0'
  elif 'AAUT_TOT' in dat['data'].dtype.names:
    fileversion = '1.1.0'
  else:
    fileversion = '1.0.0'

  hdu0.header.update('DATE', now.strftime('%d/%m/%y'))
  hdu0.header.update('E_LEVEL', "Atomic Energy Levels")
  hdu0.header.update('FILENAME', "Python routine")
  hdu0.header.update('ORIGIN', "ATOMDB",comment=os.environ['USER']+", AtomDB project")
  hdu0.header.update('HDUCLASS', "ATOMIC",comment="Atomic Data")
  hdu0.header.update('HDUCLAS1', "E_LEVEL",comment="Atomic Energy levels")
  hdu0.header.update('HDUVERS', fileversion,comment="Version of datafile")


  # do a key case conversion
  datakey=''
  for i in dat.keys():
    if i.lower()=='data':
      datakey=i
  keys={}
  for i in dat[datakey].dtype.names:
    if i.lower()=='elec_config':
      keys['elec_config'] = i
    if i.lower()=='energy':
      keys['energy'] = i
    if i.lower()=='e_error':
      keys['e_error'] = i
    if i.lower()=='n_quan':
      keys['n_quan'] = i
    if i.lower()=='l_quan':
      keys['l_quan'] = i
    if i.lower()=='s_quan':
      keys['s_quan'] = i
    if i.lower()=='lev_deg':
      keys['lev_deg'] = i
    if i.lower()=='phot_type':
      keys['phot_type'] = i
    if i.lower()=='phot_par':
      keys['phot_par'] = i
    if i.lower()=='energy_ref':
      keys['energy_ref'] = i
    if i.lower()=='phot_ref':
      keys['phot_ref'] = i
    if i.lower()=='arad_tot':
      keys['arad_tot'] = i
    if i.lower()=='aaut_tot':
      keys['aaut_tot'] = i
    if i.lower()=='aaut_ref':
      keys['aaut_ref'] = i
    if i.lower()=='arad_ref':
      keys['arad_ref'] = i


  #secondary HDU, hdu1:
  if fileversion=='1.0.0':
    hdu1 = pyfits.new_table(pyfits.ColDefs(
            [pyfits.Column(name='ELEC_CONFIG',
               format='40A',
               array=dat[datakey][keys['elec_config']]),
             pyfits.Column(name='ENERGY',
               format='1E',
               unit='eV',
               array=dat[datakey][keys['energy']]),
             pyfits.Column(name='E_ERROR',
               format='1E',
               unit='eV',
               array=dat[datakey][keys['e_error']]),
             pyfits.Column(name='N_QUAN',
               format='1J',
               array=dat[datakey][keys['n_quan']]),
             pyfits.Column(name='L_QUAN',
               format='1J',
               array=dat[datakey][keys['l_quan']]),
             pyfits.Column(name='S_QUAN',
               format='1E',
               array=dat[datakey][keys['s_quan']]),
             pyfits.Column(name='LEV_DEG',
               format='1J',
               array=dat[datakey][keys['lev_deg']]),
             pyfits.Column(name='PHOT_TYPE',
               format='1J',
               array=dat[datakey][keys['phot_type']]),
             pyfits.Column(name='PHOT_PAR',
               format='20E',
               array=dat[datakey][keys['phot_par']]),
             pyfits.Column(name='ENERGY_REF',
               format='20A',
               array=dat[datakey][keys['energy_ref']]),
             pyfits.Column(name='PHOT_REF',
               format='20A',
               array=dat[datakey][keys['phot_ref']])]
             ))

  elif fileversion=='1.1.0':
    hdu1 = pyfits.new_table(pyfits.ColDefs(
            [pyfits.Column(name='ELEC_CONFIG',
               format='40A',
               array=dat[datakey][keys['elec_config']]),
             pyfits.Column(name='ENERGY',
               format='1E',
               unit='eV',
               array=dat[datakey][keys['energy']]),
             pyfits.Column(name='E_ERROR',
               format='1E',
               unit='eV',
               array=dat[datakey][keys['e_error']]),
             pyfits.Column(name='N_QUAN',
               format='1J',
               array=dat[datakey][keys['n_quan']]),
             pyfits.Column(name='L_QUAN',
               format='1J',
               array=dat[datakey][keys['l_quan']]),
             pyfits.Column(name='S_QUAN',
               format='1E',
               array=dat[datakey][keys['s_quan']]),
             pyfits.Column(name='LEV_DEG',
               format='1J',
               array=dat[datakey][keys['lev_deg']]),
             pyfits.Column(name='PHOT_TYPE',
               format='1J',
               array=dat[datakey][keys['phot_type']]),
             pyfits.Column(name='PHOT_PAR',
               format='20E',
               array=dat[datakey][keys['phot_par']]),
             pyfits.Column(name='AAUT_TOT',
               format='1E',
               unit='s^-1',
               array=dat[datakey][keys['aaut_tot']]),
             pyfits.Column(name='ARAD_TOT',
               format='1E',
               unit='s^-1',
               array=dat[datakey][keys['arad_tot']]),
             pyfits.Column(name='ENERGY_REF',
               format='20A',
               array=dat[datakey][keys['energy_ref']]),
             pyfits.Column(name='PHOT_REF',
               format='20A',
               array=dat[datakey][keys['phot_ref']]),
             pyfits.Column(name='AAUT_REF',
               format='20A',
               array=dat[datakey][keys['aaut_ref']]),
             pyfits.Column(name='ARAD_REF',
               format='20A',
               array=dat[datakey][keys['arad_ref']])] ))

  hdu1.header.update('XTENSION', hdu1.header['XTENSION'],
          comment='Written by '+os.environ['USER']+now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu1.header.update('EXTNAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('HDUCLASS', 'ATOMIC',
          comment='Atomic Data', before="TTYPE1")
  hdu1.header.update('HDUCLAS1', 'E_LEVEL',
          comment='Energy level tables', before="TTYPE1")
  hdu1.header.update('ELEMENT', dat['Z'],
          comment='Numer of protons in element', before="TTYPE1")
  hdu1.header.update('ION_STAT', dat['z1']-1,
          comment='ion state (0 = neutral)', before="TTYPE1")
  hdu1.header.update('ION_NAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('N_LEVELS',len(dat[datakey][keys['elec_config']]) ,
           comment='Number of energy levels', before="TTYPE1")
  hdu1.header.update('HDUVERS1', '1.1.0',
           comment='Version of datafile', before="TTYPE1")

  if 'comments' in dat.keys():
    for icmt in dat['comments']:
      hdu1.header.add_comment(icmt)

  # combine hdus

  hdulist = pyfits.HDUList([hdu0,hdu1])

  # write out file (overwrite any existing file)
  if clobber:
    try:
      os.remove(fname)
    except OSError:
      pass
  
  try:
    hdulist.writeto(fname, checksum=True, clobber=clobber)
  except TypeError:
    hdulist.writeto(fname, clobber=clobber)

  print "file written: "+fname

#-------------------------------------------------------------------------------

def write_la_file(fname, dat, clobber=False):
  """
  Write the data in list dat to fname

  Parameters
  ----------
  fname : string
    The file to write
  dat : list
    The data to write
    Should be a list with the following keywords:
    Z : int
       nuclear charge
    z1 : int
       ion charge + 1
    comments : iterable of strings
       comments to append to the file
    data : numpy.array
      stores all the individual level data, with the following types:
      upper_lev : int
        Upper level of transition
      lower_lev : int
        Lower level of transition
      wavelen : float
        Wavelength of transition (A)
      wave_err : float
        Error in wavelength (A)
      einstein_a : float
        Einstein A coefficient (s-1)
      ein_a_err : float
        Error in A coefficient (s-1)
      wave_ref : string(20)
        wavelength reference (bibcode)
      ein_a_ref : string(20)
        A-value reference (bibcode)
  clobber : bool
    Overwrite existing file if it exists.
    
  Returns
  -------
  none

  """
  # start generation of new HDU:

  #primary HDU, hdu0
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  hdu0.header.update('DATE', now.strftime('%d/%m/%y'))
  hdu0.header.update('EM_LINES', "Emission Line Data")
  hdu0.header.update('FILENAME', "Python routine")
  hdu0.header.update('ORIGIN', "ATOMDB",comment=os.environ['USER']+", AtomDB project")
  hdu0.header.update('HDUCLASS', "ATOMIC",comment="Atomic Data")
  hdu0.header.update('HDUCLAS1', "EM_LINES",comment="Emission Line data")
  hdu0.header.update('HDUVERS', "1.0.0",comment="Version of datafile")

  # do a key case conversion
  datakey=''
  for i in dat.keys():
    if i.lower()=='data':
      datakey=i
  keys={}
  for i in dat[datakey].dtype.names:
    if i.lower()=='upper_lev':
      keys['upper_lev'] = i
    if i.lower()=='lower_lev':
      keys['lower_lev'] = i
    if i.lower()=='upper_lev':
      keys['upper_lev'] = i
    if i.lower()=='wavelen':
      keys['wavelen'] = i
    if i.lower()=='wave_obs':
      keys['wave_obs'] = i
    if i.lower()=='wave_err':
      keys['wave_err'] = i
    if i.lower()=='einstein_a':
      keys['einstein_a'] = i
    if i.lower()=='ein_a_err':
      keys['ein_a_err'] = i
    if i.lower()=='wave_ref':
      keys['wave_ref'] = i
    if i.lower()=='wv_obs_ref':
      keys['wv_obs_ref'] = i
    if i.lower()=='ein_a_ref':
      keys['ein_a_ref'] = i
      
  
  #secondary HDU, hdu1:
  hdu1 = pyfits.new_table(pyfits.ColDefs(
        [pyfits.Column(name='UPPER_LEV',
           format='1J',
           array=dat[datakey][keys['upper_lev']]),
         pyfits.Column(name='LOWER_LEV',
           format='1J',
           array=dat[datakey][keys['lower_lev']]),
         pyfits.Column(name='WAVELEN',
           format='1E',
           unit='Angstrom',
           array=dat[datakey][keys['wavelen']]),
         pyfits.Column(name='WAVE_OBS',
           format='1E',
           unit='Angstrom',
           array=dat[datakey][keys['wave_obs']]),
         pyfits.Column(name='WAVE_ERR',
           format='1E',
           unit='Angstrom',
           array=dat[datakey][keys['wave_err']]),
         pyfits.Column(name='EINSTEIN_A',
           format='1E',
           unit='s**-1',
           array=dat[datakey][keys['einstein_a']]),
         pyfits.Column(name='EIN_A_ERR',
           format='1E',
           unit='s**-1',
           array=dat[datakey][keys['ein_a_err']]),
         pyfits.Column(name='WAVE_REF',
           format='20A',
           array=dat[datakey][keys['wave_ref']]),
         pyfits.Column(name='WV_OBS_REF',
           format='20A',
           array=dat[datakey][keys['wv_obs_ref']]),
         pyfits.Column(name='EIN_A_REF',
           format='20A',
           array=dat[datakey][keys['ein_a_ref']])]
         ))

  hdu1.header.update('XTENSION', hdu1.header['XTENSION'],
          comment='Written by '+os.environ['USER']+now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu1.header.update('EXTNAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('HDUCLASS', 'ATOMIC',
          comment='Atomic Data', before="TTYPE1")
  hdu1.header.update('HDUCLAS1', 'EM_LINES',
          comment='Emission line tables', before="TTYPE1")
  hdu1.header.update('ELEMENT', dat['Z'],
          comment='Numer of protons in element', before="TTYPE1")
  hdu1.header.update('ION_STAT', dat['z1']-1,
          comment='ion state (0 = neutral)', before="TTYPE1")
  hdu1.header.update('ION_NAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('N_LINES',len(dat['data']) ,
           comment='Number of emission lines', before="TTYPE1")
  hdu1.header.update('HDUVERS1', '1.0.0',
           comment='Version of datafile', before="TTYPE1")

  if  'comments' in dat.keys():
    print 'adding comments'
    for icmt in dat['comments']:
      hdu1.header.add_comment(icmt)

  # combine hdus
  print 'combining HDUs'
  hdulist = pyfits.HDUList([hdu0,hdu1])

  # write out file (overwrite any existing file)
  if clobber:
    try:
      os.remove(fname)
    except OSError:
      pass

  print 'writing lafile'
  try:
    hdulist.writeto(fname, checksum=True, clobber=clobber)
  except TypeError:
    hdulist.writeto(fname, clobber=clobber)

  print "file written: "+fname

#-------------------------------------------------------------------------------

def write_ai_file(fname, dat, clobber=False):
  """
  Write the data in list dat to fname

  Parameters
  ----------
  fname : string
    The file to write
  dat : list
    The data to write
    Should be a list with the following keywords:
    Z : int
       nuclear charge
    z1 : int
       ion charge + 1
    comments : iterable of strings
       comments to append to the file
    data : numpy.array
      stores all the individual level data, with the following types:
      ion_init : int
        Inital ion state of transition
      ion_final : int
        Final ion state of transition
      level_init : int
        Initial level of transition
      level_final : int
        Final level of transition
      auto_rate : float
        Autoionization rate (s-1)
      auto_err : float
        Error in autoionization rate (s-1)
      auto_ref : string(20)
        Autoionization rate reference (bibcode)
  clobber : bool
    Overwrite existing file if it exists.
    
  Returns
  -------
  none

  """
  # start generation of new HDU:

  #primary HDU, hdu0
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  hdu0.header.update('DATE', now.strftime('%d/%m/%y'))
  hdu0.header.update('AUTOION', "Autoionization Data")
  hdu0.header.update('FILENAME', "Python routine")
  hdu0.header.update('ORIGIN', "ATOMDB",comment=os.environ['USER']+", AtomDB project")
  hdu0.header.update('HDUCLASS', "ATOMIC",comment="Atomic Data")
  hdu0.header.update('HDUCLAS1', "AUTOION",comment="Autoionization data")
  hdu0.header.update('HDUVERS', "1.0.0",comment="Version of datafile")

  #secondary HDU, hdu1:
  hdu1 = pyfits.new_table(pyfits.ColDefs(
        [pyfits.Column(name='ION_INIT',
           format='1J',
           array=dat['data']['ion_init']),
         pyfits.Column(name='ION_FINAL',
           format='1J',
           array=dat['data']['ion_final']),
         pyfits.Column(name='LEVEL_INIT',
           format='1J',
           array=dat['data']['level_init']),
         pyfits.Column(name='LEVEL_FINAL',
           format='1J',
           array=dat['data']['level_final']),
         pyfits.Column(name='AUTO_RATE',
           format='1E',
           unit='s**-1',
           array=dat['data']['auto_rate']),
         pyfits.Column(name='AUTO_ERR',
           format='1E',
           unit='s**-1',
           array=dat['data']['auto_err']),
         pyfits.Column(name='AUTO_REF',
           format='20A',
           array=dat['data']['auto_ref'])]
         ))

  hdu1.header.update('XTENSION', hdu1.header['XTENSION'],
          comment='Written by '+os.environ['USER']+now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu1.header.update('EXTNAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('HDUCLASS', 'ATOMIC',
          comment='Atomic Data', before="TTYPE1")
  hdu1.header.update('HDUCLAS1', 'AUTOION',
          comment='Autoionization data tables', before="TTYPE1")
  hdu1.header.update('ELEMENT', dat['Z'],
          comment='Numer of protons in element', before="TTYPE1")
  hdu1.header.update('ION_STAT', dat['z1']-1,
          comment='ion state (0 = neutral)', before="TTYPE1")
  hdu1.header.update('ION_NAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('N_LINES',len(dat['data']) ,
           comment='Number of emission lines', before="TTYPE1")
  hdu1.header.update('HDUVERS1', '1.0.0',
           comment='Version of datafile', before="TTYPE1")

  # combine hdus

  
  if 'comments' in dat.keys():
    print 'adding comments'
    for icmt in dat['comments']:
      hdu1.header.add_comment(icmt)

  # combine hdus
  print 'combining HDUs'
  hdulist = pyfits.HDUList([hdu0,hdu1])

  # write out file (overwrite any existing file)
  if clobber:
    try:
      os.remove(fname)
    except OSError:
      pass

  print 'writing lafile'
  try:
    hdulist.writeto(fname, checksum=True, clobber=clobber)
  except TypeError:
    hdulist.writeto(fname, clobber=clobber)

  print "file written: "+fname

#-------------------------------------------------------------------------------

def write_ec_file(fname, dat, clobber=False):
  """
  Write the data in list dat to fname

  Parameters
  ----------
  fname : string
    The file to write
  dat : list
    The data to write
    Should be a list with the following keywords:
    Z : int
       nuclear charge
    z1 : int
       ion charge + 1
    comments : iterable of strings
       comments to append to the file
    data : numpy.array
      stores all the individual level data, with the following types:
      lower_lev : int
        Lower level of transition
      upper_lev : int
        Upper level of transition
      coeff_type : int
        Coefficient type
      min_temp : float
        Minimum temperature in range (K)
      max_temp : float
        Maximum temperature in range (K)
      temperature : float(20)
        List of temperatures (K)
      effcollstrpar : float(20)
        Effective collision strength parameters
      inf_limit  : float (OPTIONAL - if type 1.2.0)
        High temperature limit point, if provided.
      reference : string(20)
        Collisional excitation reference (bibcode)
  clobber : bool
    Overwrite existing file if it exists.
    
  Returns
  -------
  none

  """

  # start generation of new HDU:

  #primary HDU, hdu0
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  # get version type
  if 'inf_limit' in dat['data'].dtype.names:
    version = '1.2.0'
  else:
    version = '1.1.0'
  hdu0.header.update('DATE', now.strftime('%d/%m/%y'))
  hdu0.header.update('COLL_STR', "Collision Strengths")
  hdu0.header.update('FILENAME', "Python routine")
  hdu0.header.update('ORIGIN', "ATOMDB",comment=os.environ['USER']+", AtomDB project")
  hdu0.header.update('HDUCLASS', "ATOMIC",comment="Atomic Data")
  hdu0.header.update('HDUCLAS1', "COLL_STR",comment="e + p Collision Strengths")
  hdu0.header.update('HDUVERS', version,comment="Version of datafile")

  #secondary HDU, hdu1:
  if version == '1.2.0':
    hdu1 = pyfits.new_table(pyfits.ColDefs(
          [pyfits.Column(name='LOWER_LEV',
             format='1J',
             array=dat['data']['lower_lev']),
           pyfits.Column(name='UPPER_LEV',
             format='1J',
             array=dat['data']['upper_lev']),
           pyfits.Column(name='COEFF_TYPE',
             format='1J',
             array=dat['data']['coeff_type']),
           pyfits.Column(name='MIN_TEMP',
             format='1E',
             unit='K',
             array=dat['data']['min_temp']),
           pyfits.Column(name='MAX_TEMP',
             format='1E',
             unit='K',
             array=dat['data']['max_temp']),
           pyfits.Column(name='TEMPERATURE',
             format='20E',
             unit='K',
             array=dat['data']['temperature']),
           pyfits.Column(name='EFFCOLLSTRPAR',
             format='20E',
             array=dat['data']['effcollstrpar']),
           pyfits.Column(name='INF_LIMIT',
             format='E',
             array=dat['data']['inf_limit']),
           pyfits.Column(name='REFERENCE',
             format='20A',
             array=dat['data']['reference'])]
           ))
  else:
    hdu1 = pyfits.new_table(pyfits.ColDefs(
          [pyfits.Column(name='LOWER_LEV',
             format='1J',
             array=dat['data']['lower_lev']),
           pyfits.Column(name='UPPER_LEV',
             format='1J',
             array=dat['data']['upper_lev']),
           pyfits.Column(name='COEFF_TYPE',
             format='1J',
             array=dat['data']['coeff_type']),
           pyfits.Column(name='MIN_TEMP',
             format='1E',
             unit='K',
             array=dat['data']['min_temp']),
           pyfits.Column(name='MAX_TEMP',
             format='1E',
             unit='K',
             array=dat['data']['max_temp']),
           pyfits.Column(name='TEMPERATURE',
             format='20E',
             unit='K',
             array=dat['data']['temperature']),
           pyfits.Column(name='EFFCOLLSTRPAR',
             format='20E',
             array=dat['data']['effcollstrpar']),
           pyfits.Column(name='REFERENCE',
             format='20A',
             array=dat['data']['reference'])]
           ))
    
  hdu1.header.update('XTENSION', hdu1.header['XTENSION'],
          comment='Written by '+os.environ['USER']+now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu1.header.update('EXTNAME', atomic.spectroscopic_name(dat['Z'], dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('HDUCLASS', 'ATOMIC',
          comment='Atomic Data', before="TTYPE1")
  hdu1.header.update('HDUCLAS1', 'COLL_STR',
          comment='Collision Strengths (e,p)', before="TTYPE1")
  hdu1.header.update('ELEMENT', dat['Z'],
          comment='Numer of protons in element', before="TTYPE1")
  hdu1.header.update('ION_STAT', dat['z1']-1,
          comment='ion state (0 = neutral)', before="TTYPE1")
  hdu1.header.update('ION_NAME', atomic.spectroscopic_name(dat['Z'], dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('N_EXCITE',len(dat['data']['upper_lev']) ,
           comment='Number of collisional excitations', before="TTYPE1")
  hdu1.header.update('HDUVERS1', version,
           comment='Version of datafile', before="TTYPE1")

  if 'comments' in dat.keys():
    print 'adding comments'
    for icmt in dat['comments']:
      hdu1.header.add_comment(icmt)


  # combine hdus

  hdulist = pyfits.HDUList([hdu0,hdu1])

  # write out file (overwrite any existing file)
  if clobber:
    try:
      os.remove(fname)
    except OSError:
      pass

  print 'writing lafile'
  try:
    hdulist.writeto(fname, checksum=True, clobber=clobber)
  except TypeError:
    hdulist.writeto(fname, clobber=clobber)
  print "file written: "+fname

#-------------------------------------------------------------------------------

def write_ir_file(fname, dat, clobber=False):
  """
  Write the data in list dat to fname

  Parameters
  ----------
  fname : string
    The file to write
  dat : list
    The data to write
    Should be a list with the following keywords:
    Z : int
       nuclear charge
    z1 : int
       ion charge + 1
    comments : iterable of strings
       comments to append to the file
    ionpot : float
       ionization potential (eV)   
    data : numpy.array
      stores all the individual level data, with the following types:
      element : int
        Nuclear Charge
      ion_init : int
        Initial ion stage
      ion_final : int
        Final ion stage
      level_init : int
        Initial level
      level_final : int
        Final level
      tr_type : string(2)
        Transition type:
        CI = collisional excitaion
        EA = excitation autoionization
        RR = radiative recombination
        DR = dieclectronic recombination
        XI = ionization, excluded from total rate calculation
        XR = recombination, excluded from total rate calculation
        (XR and XI are used to populate level directly)
      tr_index : int
        index within the file
      par_type : int
        parameter type, i.e. how the data is stored
      min_temp : float
        Minimum temperature in range (K)
      max_temp : float
        Maximum temperature in range (K)
      temperature : float(20)
        List of temperatures (K)
      ionrec_par : float(20)
        Ionization and recombination rate parameters
      wavelen : float
        Wavelength of emitted lines (A) [not used]
      wave_obs : float
        Observed wavelength of emitted lines (A) [not used]
      wave_err : float
        Error in these wavelengths (A) [not used]
      br_ratio : float
        Branching ratio of this line [not used]
      br_rat_err : float
        Error in branching ratio [not used]
      label : string(20)
        Label for the transition
      rate_ref : string(20)
        Rate reference (bibcode)
      wave_ref : string(20)
        Wavelength reference (bibcode)
      wv_obs_ref : string(20)
        Observed wavelength reference (bibcode)
      br_rat_ref : string(20)
        Branching ratio reference (bibcode)
  clobber : bool
    Overwrite existing file if it exists.
    
  Returns
  -------
  none

  """
  #primary HDU, hdu0
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  hdu0.header.update('DATE', now.strftime('%d/%m/%y'))
  hdu0.header.update('COLL_STR', "Collision Strengths")
  hdu0.header.update('FILENAME', "Python routine")
  hdu0.header.update('ORIGIN', "ATOMDB",comment=os.environ['USER']+", AtomDB project")
  hdu0.header.update('HDUCLASS', "ATOMIC",comment="Atomic Data")
  hdu0.header.update('HDUCLAS1', "IONREC",comment="Ionization/Recombination rates")
  hdu0.header.update('HDUVERS', "1.0.0",comment="Version of datafile")

  #secondary HDU, hdu1:
  hdu1 = pyfits.new_table(pyfits.ColDefs(
        [pyfits.Column(name='ELEMENT',
           format='1J',
           array=dat['data']['element']),
         pyfits.Column(name='ION_INIT',
           format='1J',
           array=dat['data']['ion_init']),
         pyfits.Column(name='ION_FINAL',
           format='1J',
           array=dat['data']['ion_final']),
         pyfits.Column(name='LEVEL_FINAL',
           format='1J',
           array=dat['data']['level_final']),
         pyfits.Column(name='LEVEL_INIT',
           format='1J',
           array=dat['data']['level_init']),
         pyfits.Column(name='TR_TYPE',
           format='2A',
           array=dat['data']['tr_type']),
         pyfits.Column(name='TR_INDEX',
           format='1J',
           array=numpy.arange(1,len(dat['data']['tr_type'])+1)),
         pyfits.Column(name='PAR_TYPE',
           format='1J',
           array=dat['data']['par_type']),
         pyfits.Column(name='MIN_TEMP',
           format='1E',
           unit='K',
           array=dat['data']['min_temp']),
         pyfits.Column(name='MAX_TEMP',
           format='1E',
           unit='K',
           array=dat['data']['max_temp']),
         pyfits.Column(name='TEMPERATURE',
           format='20E',
           unit='K',
           array=dat['data']['temperature']),
         pyfits.Column(name='IONREC_PAR',
           format='20E',
           array=dat['data']['ionrec_par']),
         pyfits.Column(name='WAVELEN',
           format='1E',
           array=dat['data']['wavelen']),
         pyfits.Column(name='WAVE_OBS',
           format='1E',
           array=dat['data']['wave_obs']),
         pyfits.Column(name='WAVE_ERR',
           format='1E',
           array=dat['data']['wave_err']),
         pyfits.Column(name='BR_RATIO',
           format='1E',
           array=dat['data']['br_ratio']),
         pyfits.Column(name='BR_RAT_ERR',
           format='1E',
           array=dat['data']['br_rat_err']),
         pyfits.Column(name='LABEL',
           format='20A',
           array=dat['data']['label']),
         pyfits.Column(name='RATE_REF',
           format='20A',
           array=dat['data']['rate_ref']),
         pyfits.Column(name='WAVE_REF',
           format='20A',
           array=dat['data']['wave_ref']),
         pyfits.Column(name='WV_OBS_REF',
           format='20A',
           array=dat['data']['wv_obs_ref']),
         pyfits.Column(name='BR_RAT_REF',
           format='20A',
           array=dat['data']['br_rat_ref'])]
         ))

  hdu1.header.update('XTENSION', hdu1.header['XTENSION'],
          comment='Written by '+os.environ['USER']+now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu1.header.update('EXTNAME', 'IONREC',
          comment='Ionization/Recombination rates', before="TTYPE1")
  hdu1.header.update('HDUCLASS', 'ATOMIC',
          comment='Atomic Data', before="TTYPE1")
  hdu1.header.update('HDUCLAS1', 'IONREC',
          comment='Ionization/Recombinatoin rates', before="TTYPE1")
  hdu1.header.update('ELEMENT', dat['Z'],
          comment='Numer of protons in element', before="TTYPE1")
  hdu1.header.update('ION_STAT', dat['z1']-1,
          comment='ion state (0 = neutral)', before="TTYPE1")
  hdu1.header.update('ION_NAME', atomic.spectroscopic_name(dat['Z'],dat['z1']),
          comment='Ion Name', before="TTYPE1")
  hdu1.header.update('N_ION',len(dat['data']['level_init']) ,
           comment='Number of rates', before="TTYPE1")
  hdu1.header.update('HDUVERS1', '1.0.0',
           comment='Version of datafile', before="TTYPE1")
  if 'ionpot' in dat.keys():
    hdu1.header.update('IONPOT', dat['ionpot'],
             comment='Ionization Potential (eV)', before="TTYPE1")
  else:
    print "WARNING: ionpot keyword not found in list"

  if 'comments' in dat.keys():
    print 'adding comments'
    for icmt in dat['comments']:
      hdu1.header.add_comment(icmt)
  # combine hdus

  hdulist = pyfits.HDUList([hdu0,hdu1])

  # write out file (overwrite any existing file)
  if clobber:
    try:
      os.remove(fname)
    except OSError:
      pass

  print 'writing irfile'
  try:
    hdulist.writeto(fname, checksum=True, clobber=clobber)
  except TypeError:
    hdulist.writeto(fname, clobber=clobber)

  print "file written: "+fname

def keyword_check(keyword):
  """
  Returns False is the keyword is in fact false, otherwise returns True
  
  Parameters
  ----------
  keyword: any 
    The keyword value
  
  Returns
  -------
  bool
    True if the keyword is set to not False, otherwise False
  """
  
  # first, check if iterable
  
  isbool=True
  
  try:
    it = iter(keyword)
    return True
  except TypeError: 
    pass


  # so it's not iterable
  
  if keyword==False:
    return False
  else:
    return True  

#-------------------------------------------------------------------------------

def write_develop_data(data, filemapfile, Z, z1, ftype, folder, froot):
  import string
  
  """
  Write the data to the next version of the file. Update the filemap.
  
  """
  
  # ok. This is hard.
  
  # 1 create a new filename
  #
  # filename is of type: "folder/APED/elsymb/elsymb_z1/elsymb_z1_ftype_froot_ITERATION.fits"
  # iteration will be 1 2 3 4 5 etc
  elsymb = atomic.Ztoelsymb(Z)
  fname1 = '%s/APED/%s/%s_%i/%s_%i_%s_%s'%\
           (folder, elsymb.lower(), elsymb.lower(),z1,elsymb.lower(),z1,ftype,froot)

  isunique=False
  i=1
  flocation =string.join(fname1.split('/')[:-1],'/')
  
  if not os.path.isdir(flocation):
    mkdir_p(flocation)
  
  while not isunique:
    fname = "%s_%i.fits"%(fname1, i)
    
    if not os.path.exists(fname):
      isunique=True
    else:
      i+=1
  # ok, we have a unique file
  
  # check we can update the filemap
  ret = atomdb.read_filemap(filemap=filemapfile)
  
  j = numpy.where((ret['Z']==Z) & \
                  (ret['z1']==z1))[0]
  if len(j)==0:
    print "Hmm... this data doesn't exist?"
    
  else:
    ret[ftype.lower()][j[0]] = fname
    atomdb.write_filemap(ret, filemapfile)
    
    data.writeto(fname)

#-------------------------------------------------------------------------------

def generate_xspec_ionbal_files(Z, filesuffix, settings = False):
  """
  Generate the eigen files that XSPEC uses to calculate the ionizatoin
  balances
  
  Parameters
  ----------
  Z : int
    atomic number of element
  filesuffix : string
    the filename will be eigenELSYMB_filesuffix.fits
  settings : dict
    This will let you override some standard inputs for get_data:

    * settings['filemap']: the filemap to use if you do not want to use
      the default $ATOMDB/filemap

    * settings['atomdbroot']: If you have files in non-standard locations
      you can replace $ATOMDB with this value
  
  Returns
  -------
   none
   
  """
  import scipy.linalg
  # get the ion & rec rates

  Telist = numpy.logspace(4,9,1251)
  
  ionlist = numpy.zeros([len(Telist), Z])
  reclist = numpy.zeros([len(Telist), Z])
  
  # outputs:
  feqb = numpy.zeros([len(Telist),Z+1])
  vl_out = numpy.zeros([len(Telist),Z**2])
  vr_out = numpy.zeros([len(Telist),Z**2])
  eig_out = numpy.zeros([len(Telist),Z])
  
  
  for z1 in range(1,Z+1):
    iontmp, rectmp = atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True,
                                            settings=settings)
    ionlist[:,z1-1] = iontmp
    reclist[:,z1-1] = rectmp


  for ite in range(len(Telist)):
    Te = Telist[ite]
    ion = ionlist[ite,:]
    rec = reclist[ite,:]
  

    b = numpy.zeros(Z+1, dtype=numpy.float32)
    a = numpy.zeros((Z+1,Z+1), dtype=numpy.float32)


    for iZ in range(0,Z):
      a[iZ,iZ] -= (ion[iZ])
      a[iZ+1,iZ] += (ion[iZ])

      a[iZ,iZ+1] += (rec[iZ])
      a[iZ+1,iZ+1] -= (rec[iZ])

    # conservation of population
    for iZ in range(0,Z+1):
      a[0,iZ]=1.0
    b[0]=1.0

    c = numpy.linalg.solve(a,b)
    c[0] = 1-sum(c[1:])
    c[c<1e-10]=0.0
    feqb[ite] = c

    ZZ=len(ion)+1
    ndim=ZZ
    AA = numpy.zeros((ndim-1,ndim-1), dtype=numpy.float32)
    # populate with stuff

    for iCol in range(ndim-1):
      for iRow in range(ndim-1):

        if (iRow==0):
          if (iCol==0):
            if (Z>2):
              AA[0,iCol] = -(ion[0] + ion[1] + rec[0])
            else:
              AA[0,iCol] = -(ion[0] + rec[0])

  
          if (iCol==1): AA[0,iCol] = rec[1] - ion[0]
          if (iCol>1):
            AA[0,iCol] = -ion[0]
        else:
          if (iRow==iCol+1):  AA[iRow,iCol]= ion[iRow]
          if (iRow==iCol):
            if (iRow+2<ndim):

              AA[iRow,iCol]=-(ion[iRow+1]+rec[iRow])
            else:
              AA[iRow,iCol]=-rec[iRow]


          if (iRow==iCol-1):
             AA[iRow,iCol]= rec[iRow+1]


    wr_la,wi_la,vl_la,vr_la, info=scipy.linalg.lapack.dgeev(AA)

    leftevec = numpy.zeros(Z**2)
    rightevec = numpy.zeros(Z**2)
  
  # The value VL in which is stored is not actually the left eigenvecotr,
  # but is instead the inverse of vr.
  
    vl = numpy.matrix(vr_la)**-1
  

    for i in range(Z):
      for j in range(Z):
        leftevec[i*Z+j] = vl[i,j]
        rightevec[i*Z+j] = vr_la[j,i]

    vr_out[ite] = rightevec
    vl_out[ite] = leftevec
    eig_out[ite] = wr_la
  

  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  hdu0.header.update('DATE', now.strftime('%d/%m/%y'))
  hdu0.header.update('FILENAME', "Python routine")
  hdu0.header.update('ORIGIN', "ATOMDB",comment=os.environ['USER']+", AtomDB project")
  
  #secondary HDU, hdu1:
  hdu1 = pyfits.new_table(pyfits.ColDefs(
        [pyfits.Column(name='FEQB',
           format='%iD'%(Z+1),
           array=feqb),
         pyfits.Column(name='EIG',
           format='%iD'%(Z),
           array=eig_out),
         pyfits.Column(name='VR',
           format='%iD'%(Z**2),
           array=vr_out),
         pyfits.Column(name='VL',
           format='%iD'%(Z**2),
           array=vl_out)] ))
  
  hdulist = pyfits.HDUList([hdu0,hdu1])
  hdulist[1].header['EXTNAME']='EIGEN'
  
  fname = 'eigen%s_%s.fits'%(atomic.Ztoelsymb(Z).lower(), filesuffix)
  hdulist.writeto(fname, checksum=True, clobber=True)

#-------------------------------------------------------------------------------

def make_release_filetree(filemapfile_in, filemapfile_out, \
                          replace_source, destination, versionname):
  """
  Take an existing filemap, copy the files to the atomdbftp folder as required.
  
  Parameters
  ----------
  filemapfile_in : string
    The existing filemap file for the new release
  filemapfile_out : string
    The filename for the produced filemap
  replace_source : string
    All new files are in this directory.
  destination : string
    The folder to store the files in
  versionname : string
    The version string for the new files (e.g. 3_0_4)
  
  Returns
  -------
  None
  
  Notes
  -----
  This code searches for any files which don't have $ATOMDB in the filename
  and assumes they are new.
  
  It updates the file name to be $ATOMDB/elname/elname_ion/elname_ion_FTYPE_versionname.fits
  
  Versionname will have its last number stsripped and replaced with "a".
  So 3_0_4_2 becomes 3_0_4_a. This reflects that 4-number versions are for revisions of a file
  under development, while 3 number + letter are for released data.
  
  And then copies it to the destination folder, compressing it with gzip.
  """
  import re, gzip
  fmap = atomdb.read_filemap(filemapfile_in, atomdbroot='XXX')
  for i in range(len(fmap['Z'])):
    for key in ['em','ci','pi','la','ai','ir','ec','lv','pc','dr']:
      if fmap[key][i] =='': continue
      if replace_source in fmap[key][i]:
        fin = fmap[key][i]
        fmapf = re.sub(replace_source, 'XXX', fin)
        fmapf = re.sub('%s[^ \t\r\n\v\f]*.fits'%(key.upper()),\
                '%s_v%s_a.fits'%(key.upper(), versionname),fmapf)

        fout = re.sub('XXX',destination, fmapf)
        
        tmp = open(fin, 'rb')
        tmpd = tmp.read()
        print "writing %s ..."%(fout+'.gz'),
        tmpo = gzip.open(fout+'.gz', 'wb')
        tmpo.write(tmpd)
        tmp.close()
        print "done"
        
        fmap[key][i] = fmapf

  atomdb.write_filemap(fmap, filemapfile_out, atomdbroot='XXX')
  print "Wrote filemap %s"%(filemapfile_out)

#-------------------------------------------------------------------------------

def make_linelist(linefile, outfile):
  """
  Create atomdb linelist file from line.fits file
  
  Parameters
  ----------
  
  linefile : string
    The filename of the line file
  outfile : string
    The output filename of the string
  
  Returns
  -------
  none
  
  """
  
  # open the line file
  d = pyfits.open(linefile)

  # get temperature units

  tunit = d[1].header['TUNIT1']

  # get temperatures

  te = d[1].data.field('kT')

  if tunit.lower()=='kev':
    te = te /const.KBOLTZ

  nte = len(te)

  linedattype = numpy.dtype({'names':['Lambda',\
                                      'Lambda_Err',\
                                      'LambdaTh',\
                                      'dLambdaTh',\
                                      'Element',\
                                      'Ion',\
                                      'UpperLev',\
                                      'LowerLev',\
                                      'Arraysize',\
                                      'PeakIndex',\
                                      'Temperature',\
                                      'Density',\
                                      'Emissivity',\
                                      'Emissivity_Err'],\
                             'formats':[numpy.float,\
                                        numpy.float,\
                                        numpy.float,\
                                        numpy.float,\
                                        numpy.int,\
                                        numpy.int,\
                                        numpy.int,\
                                        numpy.int,\
                                        numpy.int,\
                                        numpy.int,\
                                        (numpy.float,51),\
                                        (numpy.float,51),\
                                        (numpy.float,51),\
                                        (numpy.float,51)]})
   
  tmpdattype =  numpy.dtype({'names':['Lambda',\
                                      'Lambda_Err',\
                                      'LambdaTh'],\
                             'formats':[numpy.float,\
                                        numpy.float,\
                                        numpy.float]})
   
  # master data
  ldat_all = numpy.zeros(0, dtype=linedattype)

  for z0 in range(1,31):
    # go by element
    print "Starting Element %s"%(atomic.z0toelsymb(z0))
    ldat_list = {}
    ildat_list = {}
    for z1 in range(1,z0+1):
      ldat_list[z1] = numpy.zeros(1000,dtype=linedattype)
      ildat_list[z1] = 0
    
    tstart = time.time()
    for it in range(0,nte):
      i = it + 2
      print "it = %i, nlines = %i"%(it, len(ldat_list))
      dens = d[i].header['density']
      # filter the data by element
      delem = d[i].data[d[i].data['element']==z0]
      if len(delem)==0: continue
      
      for il in delem:
        z1 = il['Ion']
        ildat = ildat_list[z1]
        imatch = numpy.where((ldat_list[z1]['UpperLev']==il['UpperLev']) &\
                             (ldat_list[z1]['LowerLev']==il['LowerLev']))[0]
        if len(imatch)==0:
          # no match - new line!
#          tmpkT = numpy.zeros(51, dtype=numpy.float)
#          tmpkT[0] = te[it]
#          tmpNe = numpy.zeros(51, dtype=numpy.float)
#          tmpNe[0] = dens
#          tmpEmis = numpy.zeros(51, dtype=numpy.float)
#          tmpEmis[0] = il['Epsilon']
#          tmpEmiserr = numpy.zeros(51, dtype=numpy.float)
#          tmpEmiserr[:] = numpy.nan
#          tmp0 = numpy.float(0)


#          z = numpy.zeros(1,dtype=linedattype)
 #         z['Lambda'] =  il['Lambda']
#          z['Lambda_Err'] = il['Lambda_Err']
#          z['Element'] = il['Element']
#          z['Ion'] = il['Ion']
#          z['UpperLev'] = il['UpperLev']
#          z['LowerLev'] = il['LowerLev']
#          z['Arraysize'] = 1
#          z['Temperature'] = tmpkT
#          z['Density'] = tmpNe
#          z['Emissivity'] = tmpEmis
#          z['Emissivity_Err'] = tmpEmiserr
          
          
          ldat_list[z1]['Lambda'][ildat] = il['Lambda']
          ldat_list[z1]['Lambda_Err'][ildat] = il['Lambda_Err']
          ldat_list[z1]['Element'][ildat] = il['Element']
          ldat_list[z1]['Ion'][ildat] = il['Ion']
          ldat_list[z1]['UpperLev'][ildat] = il['UpperLev']
          ldat_list[z1]['LowerLev'][ildat] = il['LowerLev']
          ldat_list[z1]['Arraysize'][ildat] = 1
          ldat_list[z1]['Temperature'][ildat][0] = d[i].header['temperature']
          ldat_list[z1]['Density'][ildat][0] = dens
          ldat_list[z1]['Emissivity'][ildat][0]  = il['Epsilon']
          ldat_list[z1]['Emissivity_Err'][ildat][:]  = numpy.nan

          ildat_list[z1]+=1
          
          if ildat_list[z1] == len(ldat_list[z1]):
            ldat_list[z1] = numpy.append(ldat_list[z1], numpy.zeros(1000, dtype=linedattype))
          
        else:
          # update existing line
          imatch = imatch[0]
          ind = ldat_list[z1]['Arraysize'][imatch]
          ldat_list[z1]['Temperature'][imatch][ind] = d[i].header['temperature']
          ldat_list[z1]['Density'][imatch][ind] = dens
          ldat_list[z1]['Emissivity'][imatch][ind] = il['Epsilon']
          ldat_list[z1]['Arraysize'][imatch] += 1
    
    for z1 in ldat_list.keys():
      ldat_list[z1] = ldat_list[z1][:ildat_list[z1]]
      ldat_all = numpy.append(ldat_all, ldat_list[z1])
    tend = time.time()
    tdiff = tend-tstart
    print 'Element %s took %d min %d sec'%(atomic.z0toelsymb(z0),\
                                           int(tdiff)/60, int(tdiff)%60)

  # now go through and sort out the peaks
  for line in ldat_all:
    line['PeakIndex'] = numpy.argmax(line['Emissivity'])+1
  
  # get number of points
  maxpoints = numpy.max(ldat_all['Arraysize'])
  hdu1dat = {}
  hdu1dat['lambda']      = ldat_all['Lambda']
  hdu1dat['lambda_err']  = ldat_all['Lambda_Err']
  hdu1dat['lambdath']    = ldat_all['LambdaTh']
  hdu1dat['dlambdath']   = ldat_all['dLambdaTh']
  hdu1dat['element']     = ldat_all['Element']
  hdu1dat['ion']         = ldat_all['Ion']
  hdu1dat['upperlev']    = ldat_all['UpperLev']
  hdu1dat['lowerlev']    = ldat_all['LowerLev']
  hdu1dat['arraysize']   = ldat_all['Arraysize']
  hdu1dat['peakindex']   = ldat_all['PeakIndex']
  hdu1dat['temperature'] = ldat_all['Temperature'][:,:maxpoints]
  hdu1dat['density']     = ldat_all['Density'][:,:maxpoints]
  hdu1dat['emissivity']  = ldat_all['Emissivity'][:,:maxpoints]
  hdu1dat['emissivity_err']= ldat_all['Emissivity_Err'][:,:maxpoints]
  

  hdu2dat = {}
  hdu2dat['lambda']    =  hdu1dat['lambda']
  hdu2dat['lambda_err']=  hdu1dat['lambda_err']
  hdu2dat['lambdath']  =  hdu1dat['lambdath']
  hdu2dat['dlambdath'] =  hdu1dat['dlambdath']
  hdu2dat['element']   =  hdu1dat['element']
  hdu2dat['ion']       =  hdu1dat['ion']
  hdu2dat['upperlev']  =  hdu1dat['upperlev']
  hdu2dat['lowerlev']  =  hdu1dat['lowerlev']
  hdu2dat['peakemissivity'] = numpy.zeros(len(hdu2dat['lambda']), dtype = float)
  hdu2dat['peakemissivityerr'] = numpy.zeros(len(hdu2dat['lambda']), dtype = float)
  hdu2dat['peaktemperature'] = numpy.zeros(len(hdu2dat['lambda']), dtype = float)
  hdu2dat['peakdensity'] = numpy.zeros(len(hdu2dat['lambda']), dtype = float)

  for i in range(len(hdu2dat['lambda'])):
    hdu2dat['peakemissivity'][i] = hdu1dat['emissivity'][i,hdu1dat['peakindex'][i]-1]
    hdu2dat['peakemissivityerr'][i] = hdu1dat['emissivity_err'][i,hdu1dat['peakindex'][i]-1]
    hdu2dat['peaktemperature'][i] = hdu1dat['temperature'][i,hdu1dat['peakindex'][i]-1]
    hdu2dat['peakdensity'][i] = hdu1dat['density'][i,hdu1dat['peakindex'][i]-1]




  # start generation of new HDU:

  #primary HDU, hdu0
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()

  hdu0.header['DATE']= now.strftime('%d/%m/%y')
  hdu0.header['CONTENT']= ("Emissivity", "Line emission output")
  hdu0.header['FILENAME']= (linefile, 'Parent File')
  hdu0.header['ORIGIN']= ("ATOMDB",os.environ['USER']+", AtomDB project")
  hdu0.header['HDUCLASS']= ("EMISSIVITY","Line Emission Output")
  hdu0.header['HDUCLAS1']= ("SHORT_LINE","Line Emission Output")
  hdu0.header['HDUVERS']= ("1.0.0","Version of datafile")


  #secondary HDU, hdu1:
  hdu1 = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(
        [pyfits.Column(name='Lambda',
           format='1E',
           unit='Angstrom',
           array=hdu1dat['lambda']),
         pyfits.Column(name='Lambda_Err',
           format='1E',
           unit='Angstrom',
           array=hdu1dat['lambda_err']),
         pyfits.Column(name='LambdaTh',
           format='1E',
           unit='Angstrom',
           array=hdu1dat['lambdath']),
         pyfits.Column(name='dLambdaTh',
           format='1E',
           unit='Angstrom',
           array=hdu1dat['dlambdath']),
         pyfits.Column(name='Element',
           format='1J',
           array=hdu1dat['element']),
         pyfits.Column(name='Ion',
           format='1J',
           array=hdu1dat['ion']),
         pyfits.Column(name='UpperLev',
           format='1J',
           array=hdu1dat['upperlev']),
         pyfits.Column(name='LowerLev',
           format='1J',
           array=hdu1dat['lowerlev']),
         pyfits.Column(name='ArraySize',
           format='1J',
           array=hdu1dat['arraysize']),
         pyfits.Column(name='PeakIndex',
           format='1J',
           array=hdu1dat['peakindex']),
         pyfits.Column(name='Temperature',
           format=repr(numpy.max(hdu1dat['arraysize']))+'E',
                             unit='K',
           array=hdu1dat['temperature']),
         pyfits.Column(name='Density',
           format=repr(numpy.max(hdu1dat['arraysize']))+'E',
                             unit='cm**-3',
           array=hdu1dat['density']),
         pyfits.Column(name='Emissivity',
           format=repr(numpy.max(hdu1dat['arraysize']))+'E',
                             unit='photons cm**3 s**-1',
           array=hdu1dat['emissivity']),
         pyfits.Column(name='Emissivity_Err',
           format=repr(numpy.max(hdu1dat['arraysize']))+'E',
                             unit='photons cm**3 s**-1',
           array=hdu1dat['emissivity_err'])]
         ))


  hdu1.header['XTENSION']=(hdu1.header['XTENSION'],\
                           'Written by '+os.environ['USER']+\
                           now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu1.header['EXTNAME']='EMISSIVITY'
  hdu1.header['HDUNAME']=('EMISSIVITY','Spectral emission data')
  hdu1.header['HDUCLASS']=('CXC', 'Spectral data')
  hdu1.header['HDUCLAS1']=('PARAMETERS','Line Emissivity per Ion')
  hdu1.header['LOG_TEMP']=(1,'Temperature Space: 1 if log, 0 if linear')
  hdu1.header['NUM_TEMP']=(nte,'Number of temperatures in grid')
  hdu1.header['TEMPSTEP']= (d[0].header['dtemp_step'],\
                            'Temperature step size, in K or dex')
  hdu1.header['TEMPSTRT'] = (d[0].header['dtemp_start'],\
                              'Starting temperature, in K or logK')
  hdu1.header['HDUVERS1'] = ('1.0.0','version of format')







  #tertiary HDU, hdu2:
  hdu2 = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(
        [pyfits.Column(name='Lambda',
           format='1E',
           unit='Angstrom',
           array=hdu2dat['lambda']),
         pyfits.Column(name='Lambda_Err',
           format='1E',
           unit='Angstrom',
           array=hdu2dat['lambda_err']),
         pyfits.Column(name='LambdaTh',
           format='1E',
           unit='Angstrom',
           array=hdu2dat['lambdath']),
         pyfits.Column(name='dLambdaTh',
           format='1E',
           unit='Angstrom',
           array=hdu2dat['dlambdath']),
         pyfits.Column(name='Element',
           format='1J',
           array=hdu2dat['element']),
         pyfits.Column(name='Ion',
           format='1J',
           array=hdu2dat['ion']),
         pyfits.Column(name='UpperLev',
           format='1J',
           array=hdu2dat['upperlev']),
         pyfits.Column(name='LowerLev',
           format='1J',
           array=hdu2dat['lowerlev']),
         pyfits.Column(name='PeakEmissivity',
           format='1E',
                             unit='photons cm**3 s**-1',
           array=hdu2dat['peakemissivity']),
         pyfits.Column(name='PeakEmissivityErr',
           format='1E',
                             unit='photons cm**3 s**-1',
           array=hdu2dat['peakemissivityerr']),
         pyfits.Column(name='PeakTemperature',
           format='1E',
                             unit='K',
           array=hdu2dat['peaktemperature']),
         pyfits.Column(name='PeakDensity',
           format='1E',
                             unit='cm**-3',
           array=hdu2dat['peakdensity'])]
         ))


  hdu2.header['XTENSION']  =(hdu2.header['XTENSION'],\
                             'Written by '+os.environ['USER']+\
                             now.strftime('%a %Y-%m-%d %H:%M:%S')+ 'UTC')
  hdu2.header['EXTNAME']= 'PEAKEMIS'
  hdu2.header['HDUNAME']= ('PEAKEMIS', 'Spectral emission data')
  hdu2.header['HDUCLASS']= ('CXC', 'Spectral data')
  hdu2.header['HDUCLAS1']=('PARAMETERS','Line Emissivity per Ion')
  hdu2.header['HDUVERS1']=( '1.0.0','version of format')



  # combine hdus

  hdulist = pyfits.HDUList([hdu0,hdu1, hdu2])

  # write out file (overwrite any existing file)
  try:
    os.remove(outfile)
  except OSError:
    pass

  hdulist.writeto(outfile, checksum=True)

  print "file written: "+outfile

#-------------------------------------------------------------------------------

def write_ionbal_file(Te, dens, ionpop, filename, Te_linear = False, dens_linear=False):
  """
  Create ionization balance file
  
  Parameters
  ----------
  Te : array(float)
    temperatures (in K)
  dens : array(float)
    electron densities (in cm^-3)
  ionpop : dict of arrays
    one entry for each element:
    ionpop[2] = numpy.array(nion,nte, ndens)
  filename : str
    filename to write to
  Te_linear : bool
    if true, temperature grid is linear
  dens_linear : bool
    if true, density grid is linear
    
  """

  # get list of elements
  Zlist = unique(ionpop.keys())
  Zlist.sort()

  # ok, have a sorted list of elements
  n = len(Te)*len(dens)
  Z_element = numpy.zeros([n,len(Zlist)], dtype=int)
  Z_el1 = numpy.array(Zlist)
  nZ = len(Zlist) + sum(Zlist)
  X_ionpot = numpy.zeros([n, nZ], dtype=float)
  for i in range(n):
    Z_element[i,:] =Z_el1
  
  iZ = 0
  for Z in Zlist:
    X_ionpot[:,iZ:iZ+Z+1] = ionpop[Z]
    iZ+=Z+1
  
  # OK, done!
  Te_vec = numpy.zeros(n)
  Ne_vec = numpy.zeros(n)
  
  for i in range(len(dens)):
    for j in range(len(Te)):
      Te_vec[i*len(Te)+j] = Te[j]
      Ne_vec[i*len(Te)+j] = dens[i]
      
  
  
#  hdu1 = pyfits.TableHDU.from_columns(pyfits.ColDefs(
#        [pyfits.Column(name='TEMPERATURE',
#           format='1E',
#           unit='K',
#           array=Te_vec),
#         pyfits.Column(name='DENSITY',
#           format='1E',
#           unit='cm**-3',
#           array=Ne_vec),
#         pyfits.Column(name='Z_ELEMENT',
#           format='%iI'%(len(Zlist)),
#           array=Z_element),
#         pyfits.Column(name='X_IONPOP',
#           format='%iE'%(nZ),
#           array=X_ionpot)]))

  
  hdu1 = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(
        [pyfits.Column(name='TEMPERATURE',
           format='1E',
           unit='K',
           array=Te_vec),
         pyfits.Column(name='DENSITY',
           format='1E',
           unit='cm**-3',
           array=Ne_vec),
         pyfits.Column(name='Z_ELEMENT',
           format='%iI'%(len(Zlist)),
           array=Z_element),
         pyfits.Column(name='X_IONPOP',
           format='%iE'%(nZ),
           array=X_ionpot)]))
  hdu1.header['EXTNAME'] = ('ION_BAL', 'Ionization Balance table')
  hdu1.header['HDUCLASS'] = ('ION_BAL', 'Ionization Balance table')
  hdu1.header['HDUVERS1'] = ('1.0.0', 'Version of datafile')
  if Te_linear:
    tfirst = Te[0]
    tdelta = (Te[-1]-Te[0])/(len(Te)-1)
  else:
    tfirst = numpy.log10(Te[0])
    tdelta = (numpy.log10(Te[-1])-numpy.log10(Te[0]))/(len(Te)-1)

  if dens_linear:
    nfirst = dens[0]
    if len(dens)==1:
      ndelta = 0.0
    else:
      ndelta = (dens[-1]-dens[0])/(len(dens)-1)
  else:
    nfirst = numpy.log10(dens[0])
    ndelta = (numpy.log10(dens[-1])-numpy.log10(dens[0]))/(len(dens)-1)
    
  
  hdu1.header['T_FIRST'] = (tfirst, '[K] First temperature')
  hdu1.header['T_DELTA'] = (tdelta, '[K] Delta temperature')
  hdu1.header['T_LINEAR'] = (Te_linear, 'Linear (or log) temperature spacing')
  hdu1.header['T_NUMBER'] = (len(Te), 'Number of temperature grid points')
  hdu1.header['N_FIRST'] = (nfirst, '[cm**(-3)] First density')
  hdu1.header['N_DELTA'] = (ndelta, '[cm**(-3)] Delta density')
  hdu1.header['N_LINEAR'] = (dens_linear, 'Linear (or log) density spacing')
  hdu1.header['N_NUMBER'] = (len(dens), 'Number of density grid points')
  
  hdu1.header['N_ELEMEN'] = (len(Zlist), 'Number of elements')
  hdu1.header['N_IONS'] = (nZ, 'Total number of ion stages')
  
  #primary HDU, hdu0
  hdu0 = pyfits.PrimaryHDU()
  now = datetime.datetime.utcnow()
  hdu0.header['DATE']= now.strftime('%d/%m/%y')

  hdulist = pyfits.HDUList([hdu0,hdu1])
  
  hdulist.writeto(filename, checksum=True)


#-------------------------------------------------------------------------------
def make_release_tarballs(ciefileroot, neifileroot, filemap, versionname, \
                          releasenotes, makelinelist=False):
  """
  Create tarball for exmissivity files for a new release.
  
  Parameters
  ----------
  ciefileroot : string
    The path to the CIE line and coco files, with the _line.fits and _coco.fits
    ommitted.
  neifileroot : string
    The path to the NEI line and coco files, with the _line.fits and _coco.fits
    ommitted.
  filemap : string
    The filemap file
  versionname : string
    The version string for the new files (e.g. 3.0.4).
  releasenotes : string
    The file name for the release notes.
  makelinelist : bool
    Remake the line list from the line file. If not specified, assumes linelist
    file already exists.
  Returns
  -------
  None
  """

  import re, gzip, os, shutil, tarfile
  mkdir_p('tmp')
  outdir = 'tmp/atomdb_v%s'%(versionname)
  mkdir_p(outdir)

  shutil.copy2(ciefileroot+'_coco.fits',outdir+'/apec_v%s_coco.fits'%(versionname))
  shutil.copy2(ciefileroot+'_line.fits',outdir+'/apec_v%s_line.fits'%(versionname))
  shutil.copy2(filemap,outdir+'/filemap_v%s'%(versionname))
  shutil.copy2(releasenotes,outdir+'/Release_Notes.txt')
  f=open(outdir+'/VERSION', 'w')
  f.write("%s\n"%(versionname))
  f.close()

  # make the linelist
  if makelinelist:
    make_linelist(outdir+'/apec_v%s_line.fits'%(versionname), outdir+'/apec_v%s_linelist.fits'%(versionname))

  mycwd = os.getcwd()
  os.chdir(outdir)

  # make links
  os.symlink('apec_v%s_coco.fits'%(versionname), 'apec_coco.fits')
  os.symlink('apec_v%s_line.fits'%(versionname), 'apec_line.fits')
  os.symlink('apec_v%s_linelist.fits'%(versionname), 'apec_linelist.fits')
  os.symlink('filemap_v%s'%(versionname), 'filemap')

  # compress
  os.chdir('..')
  tar = tarfile.open(name='%s/atomdb_v%s.tar.bz2'%(mycwd,versionname), mode='w:bz2')
  tar.add('atomdb_v%s'%(versionname))
  tar.close()

  # make a directory for the nei stuff
  os.chdir(mycwd)
  mkdir_p('tmp/nei')
  outdir = 'tmp/nei/atomdb_v%s'%(versionname)
  mkdir_p(outdir)

  shutil.copy2(neifileroot+'_comp.fits',outdir+'/apec_v%s_nei_comp.fits'%(versionname))
  shutil.copy2(neifileroot+'_line.fits',outdir+'/apec_v%s_nei_line.fits'%(versionname))
  shutil.copy2(filemap,outdir+'/filemap_v%s'%(versionname))
  shutil.copy2(releasenotes,outdir+'/Release_Notes.txt')
  f=open(outdir+'/VERSION', 'w')
  f.write("%s\n"%(versionname))
  f.close()

  os.chdir(outdir)

  # make links
  os.symlink('apec_v%s_nei_comp.fits'%(versionname), 'apec_nei_comp.fits')
  os.symlink('apec_v%s_nei_line.fits'%(versionname), 'apec_nei_line.fits')
  #os.symlink('apec_v%s_linelist.fits'%(versionname), 'apec_v%s_linelist.fits')
  os.symlink('filemap_v%s'%(versionname), 'filemap')

  # compress
  os.chdir('..')
  tar = tarfile.open(name='%s/atomdb_v%s_nei.tar.bz2'%(mycwd,versionname), mode='w:bz2')
  tar.add('atomdb_v%s'%(versionname))
  tar.close()

  print "Tarballs written"



def generate_equilibrium_ionbal_files(filename, settings = False):
  """
  Generate the eigen files that XSPEC uses to calculate the ionizatoin
  balances
  
  Parameters
  ----------
  filename : string
    file to write
  settings : dict
    This will let you override some standard inputs for get_data:

    * settings['filemap']: the filemap to use if you do not want to use
      the default $ATOMDB/filemap

    * settings['atomdbroot']: If you have files in non-standard locations
      you can replace $ATOMDB with this value
  
  Returns
  -------
   none
   
  """

  # get the ion & rec rates

  Telist = numpy.logspace(4,9,501)
  Nelist = numpy.array([1.0])
  ionbal = {}
  for Z in range(1,31):
    print Z
    ionlist = numpy.zeros([len(Telist), Z])
    reclist = numpy.zeros([len(Telist), Z])
  
  # outputs:
    feqb = numpy.zeros([len(Telist),Z+1])
  
  
    for z1 in range(1,Z+1):
      iontmp, rectmp = atomdb.get_ionrec_rate(Telist, False, Z=Z, z1=z1, extrap=True,
                                              settings=settings)
      ionlist[:,z1-1] = iontmp
      reclist[:,z1-1] = rectmp


    for ite in range(len(Telist)):
      Te = Telist[ite]
      ion = ionlist[ite,:]
      rec = reclist[ite,:]
  

      b = numpy.zeros(Z+1, dtype=numpy.float32)
      a = numpy.zeros((Z+1,Z+1), dtype=numpy.float32)

 
      for iZ in range(0,Z):
        a[iZ,iZ] -= (ion[iZ])
        a[iZ+1,iZ] += (ion[iZ])

        a[iZ,iZ+1] += (rec[iZ])
        a[iZ+1,iZ+1] -= (rec[iZ])

    # conservation of population
      for iZ in range(0,Z+1):
        a[0,iZ]=1.0
      b[0]=1.0

      c = numpy.linalg.solve(a,b)
      c[0] = 1-sum(c[1:])
      c[c<1e-10]=0.0
      feqb[ite] = c
    ionbal[Z] = feqb
  
  write_ionbal_file(Telist, Nelist, ionbal, filename, Te_linear = False, dens_linear=True)
