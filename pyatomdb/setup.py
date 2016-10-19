from setuptools import setup, Extension

def get_version(relpath):
    """read version info from file without importing it"""
    from os.path import dirname, join
    for line in open(join(dirname(__file__), relpath)):
        if '__version__' in line:
            if '"' in line:
                # __version__ = "0.9"
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


readme = open('README.txt', 'r')
README_TEXT = readme.read()
readme.close()

linapprox =  Extension('linear_approx',\
                       define_macros = [('MAJOR_VERSION', '1'),\
                                        ('MINOR_VERSION','0')],\
                       sources=['linear_approx.c'])



setup(name='pyatomdb',
      version=get_version('pyatomdb/__init__.py'),
      description='AtomDB python library. This is a very early development version.',
      url='http://www.atomdb.org',
      author='Adam Foster',
      author_email='afoster@cfa.harvard.edu',
      license='Smithsonian',
      packages=['pyatomdb'],
      classifiers=['Development Status :: 2 - Pre-Alpha',\
                   'Environment :: Console',\
                   'Intended Audience :: Developers',\
                   'Intended Audience :: Education',\
                   'Intended Audience :: End Users/Desktop',\
                   'Intended Audience :: Science/Research',\
                   'Topic :: Scientific/Engineering :: Astronomy',\
                   'Topic :: Scientific/Engineering :: Physics',\
                   'Programming Language :: Python',\
                   'Operating System :: POSIX'],
      zip_safe=False,
      long_description = README_TEXT,\
      install_requires=[
      "requests",\
      "wget",\
      "numpy",\
      "scipy"],
      ext_modules = [linapprox])

