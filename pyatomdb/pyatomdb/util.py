"""
util.py contains a range of miscellaneous helper codes that assist in running
other AtomDB codes but are not in any way part of a physical calculation.

Version -.1 - initial release
Adam Foster July 17th 2015
"""

import numpy, os, errno, hashlib
import requests, urllib, time, subprocess, shutil, wget, glob
import const
from StringIO import StringIO
import ftplib

################################################################################
#
#  Python Module
#
#  Name:        arflib.py
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


#*******************************************************************************
#
#  Routine fig
#
#  Differentiates vector y wrt x
#
#  input: x,y, lowend=1, highend=1
#
#   x - numpy.array of x values
#   y - numpy.array of y values
#   lowend - 0,1,2: 0 set dy/dx=0 at low end
#                   1 set d2y/dx2=0 at low end
#                   2 set d3y/dx3=0 at low end
#   highend - 0,1,2: 0 set dy/dx=0 at high end
#                   1 set d2y/dx2=0 at high end
#                   2 set d3y/dx3=0 at high end
#
#  returns: numpy array of derivatives
#
#  First Version:
#       Adam Foster, 28-Jul-2009
#
#*******************************************************************************

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
     
     

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


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

    os.chdir(os.path.expandvars(re.sub('VERSION',version, "$ATOMDB/tmp/atomdb_vVERSION")))


    for ifile in glob.glob('*%s*'%(version)):

      outfile = os.path.expandvars("$ATOMDB/%s"%(ifile))
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
  for flist in ['apec_vVERSION_line.fits', 'apec_vVERSION_coco.fits',\
                'filemap_vVERSION', 'apec_vVERSION_linelist.fits']:
                  
    # remove existing link if there is one
#    print "CHECKING FOR %s/%s"%(os.getcwd(), re.sub('VERSION',version,flist))
#    if os.path.exists(re.sub('VERSION',version,flist)):
#      os.remove(re.sub('VERSION',version,flist))
#      print "REMOVING LINK1!"
    # add new link

    if os.path.islink(re.sub('_vVERSION','',flist)):
      os.remove(re.sub('_vVERSION','',flist))
      
    os.symlink(re.sub('VERSION',version,flist), re.sub('_vVERSION','',flist))
    
    # update version
    a = open('VERSION','w')
    a.write(version+'\n')
    a.close()
  os.chdir(startdir)
