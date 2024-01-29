CRYSTALpytools.base.basisset module
===================================

.. _ref-base-basisset:

A base module for CRYSTAL basis set, including ``BasisSetBASE``, which is the class for basis sets of the whole system and ``AtomBS``, which is the class of a single atom. Definitions of shells and Gaussian type functions (GTF) are saved as dictionary under ``AtomBS().shells`` attribute. The `Basis Set Exchange (BSE) Python API <https://molssi-bse.github.io/basis_set_exchange/index.html>`_ is used to get and convert basis sets.

To call and modify parameters of a certain GTF, which is usually needed when optimizing basis set:

.. code-block::

    >>> bs = BasisSetBASE.from_bse('6-311G*', [6, 1]) # Download 6-311G* BS from BSE
    >>> print(bs.atoms[0].shells[1]['angular momentum']) # Get the angular momentum
    >>> bs.atoms[0].shells[1]['orbitals'][2][0] = 1.46000 # Change the exponent of the 3rd GTF, 2nd shell, 1st atom
    >>> bs.atoms[0].shells[1]['orbitals'][2][1] = 0.001 # Contraction
    >>> bs.atoms[0].shells[1]['orbitals'][2][2] = 0.815854 # sp coefficient

.. automodule:: CRYSTALpytools.base.basisset
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
