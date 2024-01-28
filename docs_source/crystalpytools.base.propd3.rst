CRYSTALpytools.base.propd3 module
=================================

.. _ref-base-propd3:

The ``Properties_inputBASE`` object is strictly structured by 'blocks', which, in general, is defined as keywords that are closed by 'END'. It is inherited from the :ref:`BlockBASE <ref-base-inputbase>` object and is inherited by the :ref:`Properties_input <ref-crystalio>` object. All the blocks are organized in layers and each corresponds to a list of keywords that can be called and set. The current structure of ``Properties_inputBASE`` is listed below:

Layer 1: Optional, repeated block (same calculation, another time) ``append1`` to ``append5``  

Layer 2: Data grid and DFT/MP2 correlation energy ``ECHG``, ``POTM``, ``CLAS``, ``EDFT/ENECOR``, ``ADFT/ACOR``

For the usages of :ref:`BlockBASE <ref-base-inputbase>` and :ref:`Properties_input <ref-crystalio>` objects, please refer to the corresponding documentations.


Examples
~~~~~~~~

Note that methods listed below are for :ref:`Properties_input <ref-crystalio>` objects, which is typically consistent with ``Properties_inputBASE``, except file read and write functions.

To set a band structure and a projected doss calculation:

.. code-block::

    >>> obj = Properties_input()
    >>> obj.band('Band calc title', 3, 6, 188, 203, 224, 1, 0,
                 [[[0, 0, 0], [3, 0, 0]], # A 3*2*3 list, for 3 line segments, 2 ending points of each segment, xyz for each point
                  [[3, 0, 0], [2, 2, 0]],
                  [[2, 2, 0], [0, 0, 0]]])
    >>> obj.newk(8, 16, 1, 0)
    >>> obj.doss(1, 600, 203, 224, 1, 12, 0, [[-1, 57], [-1, 64]]) # Project to atom 54 and atom 64

By calling the 'block-like' attribute, a sub-block object will be automatically generated if no such object is saved in the upper block object. It also be initialized and deleted in a similar way as keyword commands.

.. code-block::

    >>> obj = Properties_input()
    >>> obj.echg() # Initialize ECHG
    >>> obj.echg(None) # Remove ECHG
    >>> obj.echg(0, 95) # Charge density map, Npoint of MAPNET is 95 (default value 100)

To set the values of keywords, their names are called as methods:

.. code-block::

    >>> obj.echg.coordina([-2.498, 0., 1.696], [-2.498, 0., -1.696], [-1.249, -2.164, -1.696])
    >>> obj.echg.rectangu() # Print RECTANGU keyword, without value
    >>> obj.echg.margins(3, 3, 3, 3)

Though one can set CRYSTAL input object by manually setting up all the attributes, it is also possible to read a template d3 file and do modifications.

.. code-block::

    >>> obj = Properties_input('charge2d.d3')
    >>> obj.echg(None) # Remove ECHG block
    >>> obj.ech3(100) # Define ECH3 block
    >>> obj.ech3.range(-10, 10) # Range of Non-periodic direction
    >>> obj.to_file('charge3d.d3') # Print it into file

It is also possible to set individual blocks by a string. This is achieved by simply assigning the string variable as input when calling the corresponding method.

.. code-block::

    >>> obj.ech3('ECH3\n100\nRANGE\n-10\n10\n')

As stressed in the doc of :ref:`base.inputbase <ref-base-inputbase>`, repeated keywords not protected by sub-blocks are not permitted. That leads to problems when, for example, plotting the charge difference map to analyze bonds, where 'ECHG' is repeated twice in the main block. To address this, the following sub-block is called when the second 'ECHG' is used:

.. code-block::

    >>> obj = Properties_input()
    >>> obj.echg(0) # Set the first ECHG. MAPNET density = 100, default value is used.
    >>> obj.echg.coordina([-2.498, 0., 1.696], [-2.498, 0., -1.696], [-1.249, -2.164, -1.696])
    >>> obj.echg.rectangu()
    >>> obj.echg.margins(3, 3, 3, 3)
    >>> obj.append1() # Initialize the first appended calculation
    >>> obj.append1.pato(1, 0)
    >>> obj.append1.echg(0) # Set the second ECHG. MAPNET density = 100, default value is used.
    >>> obj.append1.echg.coordina([-2.498, 0., 1.696], [-2.498, 0., -1.696], [-1.249, -2.164, -1.696])
    >>> obj.append1.echg.rectangu()
    >>> obj.append1.echg.margins(3, 3, 3, 3)

5 appended calculations (``append1`` to ``append5``) at most can be added and the ``Properties_input`` object is the initial calculation, so the same keyword, which is not protected by subblocks, can repeat 6 times at most. Error is reported if it appears more than 6 times.


To examine the data in a block object, including the :ref:`Properties_input <ref-crystalio>` obj itself, call the ``data`` attribute.

.. code-block::

    >>> print(obj.data)


.. automodule:: CRYSTALpytools.base.propd3
   :members:
   :private-members:
   :undoc-members:
   :show-inheritance:
