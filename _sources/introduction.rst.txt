.. image:: _static/CRYSTAL_logo.png
   :width: 200 px
   :alt: CRYSTALpytools
   :align: center

Introduction
============

The CRYSTALpytools python package contains functions that allow the user to
access the `CRYSTAL code <https://www.crystal.unito.it/index.html>`_ input,
output and execution from a python infrastructure, such as Jupyter Notebooks.
Although they can achieve this goal on their own, they achieve their full
potential when used together with `PyMatGen <https://pymatgen.org/index.html>`_.
In this latter scenario, the CRYSTALpytools could be seen as a layer between
CRYSTAL and pymatgen.

In January 2022 the first stable version (v2022.1.10) was released. The latest
release is available on both `PyPI <https://pypi.org/project/CRYSTALpytools/>`_
and `conda forge <https://conda-forge.org/>`_. For details please refer to the
:ref:`Installation <ref-installation>` section.


Structure
=========

The CRYSTALpytools module aims at providing the user a python interface to the
CRYSTAL code. The central data structure, called ``Crystal_object``, is created
by the ``crystal_io`` by parsing CRYSTAL input/output files. The flowchart below
is aimed at showing how different parts of the module interact with the
``Crystal_objects``.

.. figure:: _static/crystal_object.png
   :width: 600 px
   :alt: crystal_object
   :align: center

   The general structure of CRYSTALpytools

