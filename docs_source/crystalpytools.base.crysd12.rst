CRYSTALpytools.base.crysd12 module
==================================

.. _ref-base-crysd12:

The ``Crystal_inputBASE`` object is strictly structured by 'blocks', which, in general, is defined as keywords that are closed by 'END'. It is inherited from the :ref:`BlockBASE <ref-base-inputbase>` object and is inherited by the :ref:`Crystal_input <ref-crystalio>` object. All the blocks are organized in layers and each corresponds to a list of keywords that can be called and set. The current structure of ``Crystal_inputBASE`` is listed below:

Layer 1: ``geom``, ``basisset``, ``scf``  

Layer 2: ``optgeom``, ``freqcalc``, ``dft``, ``dftd3``, ``gcp``, ``fixindex``  

Layer 3: ``preoptgeom``, ``geom``, ``base``  

For the usages of :ref:`BlockBASE <ref-base-inputbase>` and :ref:`Crystal_input <ref-crystalio>` objects, please refer to the corresponding documentations.

For example, to set force convergence threshold of a optimization run:

.. code-block::

    >>> obj = Crystal_input()
    >>> obj.geom.optgeom.toldeg(0.0001)

In principle, by calling the 'block-like' attribute, a 'block-like' object will be automatically generated if the attribute is empty. The exception is the 3rd layer attributes, which must be set by ``set_attr()`` method. A warning message is printed to indicate the name of the opened sub-block since it usually does not correspond to CRYSTAL keywords to avoid potential conflicts.

.. code-block::

    >>> obj.geom.freqcalc.set_preoptgeom()
    >>> obj.geom.freqcalc.optgeom.toldeg(0.0003)

Methods and sub-blocks of ``Crystal_inputBASE`` usually have the same name as corresponding keywords. One can setup, change or clean the keyword by calling the corresponding method.

.. code-block::

    >>> obj.scf.toldee(9) # Set SCF TOLDEE = 9
    >>> obj.scf.toldee('') # Clean the TOLDEE keyword and value
    >>> obj.scf.ppan() # Print PPAN keyword, without value

Though one can set CRYSTAL input object by manually setting up all the attributes, it is also possible to read a template d12 file and do modifications.

.. code-block::

    >>> obj.from_file('opt.d12')
    >>> obj.geom.optgeom('') # Remove OPTGEOM block
    >>> obj.to_file('scf.d12') # Print it into file

It is also possible to set individual blocks by a string. The ``set_block`` method should be used. The keyword for the block itself should not be included.

.. code-block::

    >>> obj.scf.set_dft('SPIN\nEXCHANGE\nPBE\nCORRELAT\nP86\n')

For basis set, it is not a typical ``BlockBASE`` object (though it inherits ``BlockBASE``). When 'BASISSET' keyword is used, it is called in the same way as other blocks. When explicit definitions of basis set are used, it can be defined via formatted string, file, `Basis Set Exchange (BSE) <https://molssi-bse.github.io/basis_set_exchange/index.html>`_ and :ref:`BasisSetBASE <ref-base-basisset>` object. The ending line '99 0' is required.

.. code-block::

    >>> obj.basisset.basisset('def2-SVP')
    >>> obj.basisset.from_file('mybasis.txt')
    >>> obj.basisset.from_bse('6-311G*', ['C', 'H', 'O'])

To examine the data in a block object, including Crystal_input obj itself, call the ``data`` attribute.

.. code-block::

    >>> obj.data


.. automodule:: CRYSTALpytools.base.crysd12
   :members:
   :private-members:
   :undoc-members:
   :show-inheritance:
