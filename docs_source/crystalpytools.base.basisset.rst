CRYSTALpytools.base.basisset module
===================================

.. _ref-base-basisset:

A base module for CRYSTAL basis set, including objects respectively for the basis set file, atom, shell and Gaussian type functions (GTF). The `Basis Set Exchange (BSE) Python API <https://molssi-bse.github.io/basis_set_exchange/index.html>`_ is used to get and convert basis sets.

To call and modify parameters of a certain GTF, which is usually needed when optimizing basis set:

.. code-block::

    >>> bs = BasisSetBASE.from_bse('6-311G*', ['C', 'H']) # Download 6-311G* BS from BSE
    >>> bs.atom[6].shell[1].gtf[2].exp = 1.46000 # Change the exponent of the 3rd GTF, 2nd shell, C atom (called by conventional atomic number)
    >>> bs.atom[6].shell[1].gtf[2].contr = 0.001 # Contraction
    >>> bs.atom[6].shell[1].gtf[2].pceof = 0.815854 # sp coefficient

.. automodule:: CRYSTALpytools.base.basisset
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
