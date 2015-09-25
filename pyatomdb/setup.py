from setuptools import setup

readme = open('README.txt', 'r')
README_TEXT = readme.read()
readme.close()


setup(name='pyatomdb',
      version='0.0.1.1',
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
      long_description = README_TEXT)
