Installation
============

.. _ref-installation:

Create a conda/anaconda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is not mandatory, but it makes using CRYSTALpytools very smooth. It
is, therefore, very recommended. If you are new to anaconda, please follow
`these steps <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
to install it on your computer.

Create a new conda environment:

.. code-block:: console

   conda create --name crystal python=3.9

In the line above, ``crystal`` is the name of the environment and can be set to
any you like. The ``python=3.9`` ensures that the suitable python distribution
is installed.

Activate the conda environment:

.. code-block:: console

   conda activate crystal

Install CRYSTALpytools
~~~~~~~~~~~~~~~~~~~~~~

By PyPI
-----------------

CRYSTALpytools can be installed from ``pip``. ``pip`` is a package-management
system written in Python and is used to install and manage software packages
(called modules in python).

.. code-block:: console

   pip install --upgrade CRYSTALpytools

Windows users might need to install windows-curses. This can be done by using:

.. code-block:: console

   pip install windows-curses

By Conda Forge
-----------------

CRYSTALpytools is also available on `Conda Forge <https://conda-forge.org/>`_:

.. code-block:: console

   conda install CRYSTALpytools -c conda-forge

Check Installation
------------------

To check that CRYSTALpytools was installed successfully please type:

.. code-block:: console

   conda list | grep -i 'crystalpytools'

For Windows users, please use:

.. code-block:: console

   conda list | findstr /i crystalpytools

This will return to entries containing the text string 'crystalpytools' (case
insensitive). Here there should be 'crystalpytools' with version number (e.g.
'2023.4.4') and source ('pypi'). If this is not the case, something went wrong
during the installation. Please check the location of the environment that is
being displayed. This appears at the beginning of the “conda list” command. The
most common mistake at this stage is that the environment was not activated as
described above.

Please note that both ``pip`` and ``conda`` can only install the functions and
not the example / test cases (provided as Jupyter Notebooks. This decision was
taken in order to reduce the volume of data transferred when installing. If you
are interested in the example notebooks please read the 'Testing' section below.

Set the path to runcry and runprop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you intend to run CRYSTAL on the machine where you are running the
CRYSTALpytools, the path to your local runcry and runprop needs to be specified.
To do so, please run the set_runcry_path and set_runprop_path functions:

.. code-block:: console

   python 3
   >>> from CRYSTALpytools.execute import set_runcry_path, set_runprop_path
   >>>
   >>> set_runcry_path('path_to_your_runcry')
   >>> set_runprop_path('path_to_your_runcry')

Testing
~~~~~~~

To test the CRYSTALpytools please run the test notebook that can be found in
the `examples/ <https://github.com/crystal-code-tools/CRYSTALpytools/tree/main/examples>`_. 
directory. Code examples there also provide a convenient way to learn how the code
works. The tested example code is available in the :ref:`Examples and Test Cases <ref-examples>`
section.
