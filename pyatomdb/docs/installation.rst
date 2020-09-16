============
Installation
============

.. warning::
  PyAtomDB runs only under Python 3. It will not work on Python 2. If you
  need to install Python 3 in your system, you can use your package manager
  or there are many other sources which
  can help you including `Anaconda <https://www.anaconda.com/>`_.

  Once you have Python 3 installed, you may (depending on your system)
  have to add a ``3`` to many
  of the command line commands, e.g. ``python`` becomes ``python3``, or
  ``pip`` becomes ``pip3``. Commands in this guide omit the ``3``. Note
  that actual Python code is unaffected.



PyAtomDB can be installed in two ways:

  #. From `PyPI <https://pypi.org/>`_ , using the simple ``pip install pyatomdb`` command.
  #. From `GitHub <https://github.com/AtomDB/pyatomdb>`_ , using the command ``git clone https://github.com/AtomDB/pyatomdb.git``
     to get the source, then ``python setup.py develop`` to install links to the source
     in your Python path.

Note that for both of these options the ``--user`` option can be useful, as it will install
software in your local path if you do not have administrator priviledges on your machine.

You can check that the installation was successful by running:

.. code-block:: python

  >>> import pyatomdb

If it does not immediately throw out an error, it has been successful. It will
then start asking about installing the AtomDB files, see the next section. Note
that there is no longer a need to run the initialize script.


.. warning::
  pycURL issues can arise when installing. If your install above worked without
  errors, then you are fine. If you encountered errors, and they are related to
  pycurl, you will need to consult your system's package managed (or conda if
  that is what you are using) and install pycurl separately - 
  e.g. ``conda install pycurl``.
  
  I'm not sure why this refuses to install directly with pip sometimes, but it seems
  to be a recurring feature.
  
  
----------------
ATOMDB Directory
----------------
Whenever you import the PyAtomDB module, it performs a check for the $ATOMDB directory.
This directory is where the AtomDB data files will be stored. These are not
distributed with the python package as they are large and most people will only need
a few. PyAtomDB will download the APED data files on demand as you require them, and
they are then stored in this directory until you manually delete them. If you need to
recover disk space, you can delete anything in the $ATOMDB/APED directory without
repercussions - PyAtomDB will re-download the files if it needs them in the future.

You will be asked to select a directory for installation and then whether to download
the emissivity files. It is important that this directory is one where you have write
access, as in the future further files will be added there by the code automatically.

Once installation is complete, ensure that you add the ATOMDB variable to your
shell startup file. Assuming you have installed into /home/username/atomdb,
in bash, add this line to your ~/.bashrc or ~/.bash_profile files
(depends on which one is sourced by your system):

.. code-block:: bash

  export ATOMDB=/home/username/atomdb

or for csh, add this to your ~/.cshrc or ~/.cshrc.login:

.. code-block:: csh

  setenv ATOMDB /home/username/atomdb

Recent version of Mac OS have moved to zsh, in which case modify your ~/.zshrc file as for bash above.

---------------------
Usage Data Collection
---------------------

You will also be asked about anonymous usage data. In order to track roughly how many
people are using PyAtomDB, a randomly generated number is created when you install
PyAtomDB and stored in your ``$ATOMDB/userdata`` file. Whenever PyAtomDB has to fetch
a new file this number, the filename and the current timestamp is stored on our
system so we can estimate how many users there are. We have no way to connect this
to actual individuals, it simply tells us roughly how many unique active users
there are.

If you decline, this number is set to 00000000, and otherwise PyAtomDB functions
as normal.
