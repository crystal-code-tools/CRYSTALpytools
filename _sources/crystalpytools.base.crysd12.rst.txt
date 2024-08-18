CRYSTALpytools.base.crysd12 module
==================================

.. _ref-base-crysd12:

The ``Crystal_inputBASE`` object is strictly structured by 'blocks', which, in general, is defined as keywords that are closed by 'END'. It is inherited from the :ref:`BlockBASE <ref-base-inputbase>` object and is inherited by the :ref:`Crystal_input <ref-crystalio>` object. All the blocks are organized in layers and each corresponds to a list of keywords that can be called and set. The current structure of ``Crystal_inputBASE`` is listed below:

Layer 1: ``geom``, ``basisset``, ``scf``  

Layer 2: ``optgeom``, ``freqcalc``, ``hf3c``, ``hfsol3c``, ``dft``, ``dftd3``, ``gcp``, ``geom``, ``base``, ``geba``  

Layer 3: ``preoptgeom``  

For the usages of :ref:`BlockBASE <ref-base-inputbase>` and :ref:`Crystal_input <ref-crystalio>` objects, please refer to the corresponding documentations.


Examples
~~~~~~~~

Note that methods listed below are for :ref:`Crystal_input <ref-crystalio>` objects, which is typically consistent with ``Crystal_inputBASE``, except file read and write functions.

To set force convergence threshold of a optimization run:

.. code-block:: python

    >>> obj = Crystal_input()
    >>> obj.geom.optgeom.toldeg(0.0001)

By calling the 'block-like' attribute, a sub-block object will be automatically generated if no such object is saved in the upper block object. It also be initialized and deleted in a similar way as keyword commands.

.. code-block:: python

    >>> obj = Crystal_input()
    >>> obj.geom.optgeom() # Initialize OPTGEOM
    >>> obj.geom.optgeom(None) # Remove OPTGEOM

To set the values of keywords, their names are called as methods with a really rare exception (1 so far):

.. code-block:: python

    >>> obj.geom.freqcalc.preoptgeom() # Initialize PREOPTGEOM
    >>> obj.geom.freqcalc.preoptgeom.toldeg(0.0003)
    >>> obj.scf.toldee(9) # Set SCF TOLDEE = 9
    >>> obj.scf.toldee(None) # Clean the TOLDEE keyword and value
    >>> obj.scf.ppan() # Print PPAN keyword, without value

The only exception is given below. However, when doing text analysis and printing formatted string to files, the correct keyword is recognized and printed.

.. code-block:: python

    >>> obj.geom.molecule2() # MOLECULE keyword to extract molecule from lattice. Renamed to address the conflict with modelling keyword MOLECULE

Though one can set CRYSTAL input object by manually setting up all the attributes, it is also possible to read a template d12 file and do modifications.

.. code-block:: python

    >>> obj = Crystal_input('opt.d12')
    >>> obj.geom.optgeom(None) # Remove OPTGEOM block
    >>> obj.to_file('scf.d12') # Print it into file

It is also possible to set individual blocks by a string. This is achieved by simply assigning the string variable as input when calling the corresponding method.

.. code-block:: python

    >>> obj.scf.dft('SPIN\nEXCHANGE\nPBE\nCORRELAT\nP86\n')

For basis set, it is not a typical ``BlockBASE`` object (though it inherits ``BlockBASE``). When 'BASISSET' keyword is used, it is called in the same way as other blocks. When explicit definitions of basis set are used, it can be defined via formatted string, file, `Basis Set Exchange (BSE) <https://molssi-bse.github.io/basis_set_exchange/index.html>`_ and :ref:`BasisSetBASE <ref-base-basisset>` object. The ending line '99 0' is required.

.. code-block:: python

    >>> obj.basisset.basisset('def2-SVP')
    >>> obj.basisset.from_file('mybasis.txt')
    >>> obj.basisset.from_bse('6-311G*', [6, 1, 8]) # conventional atomic numbers are supported.

Wrapper methods (``bs_user()`` and ``bs_keyword()``) are added in :ref:`Crystal_input <ref-crystalio>`. For details please checkt the manual :ref:`there <ref-crystalio>`.

.. code-block:: python

    >>> obj.bs_keyword('def2-SVP')
    >>> obj.bs_user('mybasis.txt')
    >>> obj.bs_user('6-311G*', [6, 1, 8]) # bs_user accepts file, string and BSE variables

To examine the data in a block object, including the :ref:`Crystal_input <ref-crystalio>` obj itself, call the ``data`` attribute.

.. code-block:: python

    >>> print(obj.data)


.. automodule:: CRYSTALpytools.base.crysd12
   :members:
   :private-members:
   :undoc-members:
   :show-inheritance:
