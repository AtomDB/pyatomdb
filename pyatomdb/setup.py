from setuptools import setup, Extension
import os

import sys


on_rtd = os.environ.get('READTHEDOCS')=='True'


def get_version(relpath):
    """read version info from file without importing it"""
    from os.path import dirname, join
    for line in open(join(dirname(__file__), relpath)):
        if '__version__' in line:
            if '"' in line:
                # __version__ = "0.10"
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


readme = open('README.rst', 'r')
README_TEXT = readme.read()
readme.close()

linapprox =  Extension("linear_approx",['linear_approx.c'],\
                       define_macros = [('MAJOR_VERSION', '1'),\
                                        ('MINOR_VERSION','0')])





if on_rtd:
  from mock import Mock as MagicMock

  extmos= []

  class Mock(MagicMock):
      @classmethod
      def __getattr__(cls, name):
              return Mock()

  MOCK_MODULES = ['liblinapprox']
  sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

else:
  extmos= [linapprox]

setup(name='pyatomdb',
      version=get_version('pyatomdb/__init__.py'),
      description='AtomDB python library.',
      url='http://www.atomdb.org',
      author='Adam Foster',
      author_email='afoster@cfa.harvard.edu',
      license='Smithsonian',
      packages=['pyatomdb','tests'],
      classifiers=['Development Status :: 4 - Beta',\
                   'Environment :: Console',\
                   'Intended Audience :: Developers',\
                   'Intended Audience :: Education',\
                   'Intended Audience :: End Users/Desktop',\
                   'Intended Audience :: Science/Research',\
                   'Topic :: Scientific/Engineering :: Astronomy',\
                   'Topic :: Scientific/Engineering :: Physics',\
                   'Programming Language :: Python :: 3',\
                   'Operating System :: POSIX'],
      zip_safe=False,
      long_description = README_TEXT,\
      install_requires=[
      "requests",\
      "wget",\
      "numpy>=1.9.0",\
      "scipy>=1.9.0",\
      "joblib",\
      "mock",\
      "astropy",\
      "pycurl"],
      ext_modules = extmos)

