CRYSTALpytools.base.inputbase module
====================================

.. _ref-base-inputbase:

Python wrapper for CRYSTAL inputs. Provides basic classes, methods and attributes to read, operate and write CRYSTAL d12/d3 files.

Structure of classes
--------------------

#. CRYSTAL keywords are used as methods. Calling a method modifies the corresponding attributes.

#. Keyword-like attributes (keyword not closed by 'END') are named as ``_attr`` and the corresponding method should be ``attr()``.

#. Block-like attributes (keyword closed by 'END') are named as ``_block_attr``. The corresponding method is ``set_attr()``. Another method with ``@property`` decorator should be set in such way that when called without ``_block_attr`` existing, create a new one, otherwise return to it.

#. ``_block_bg`` attribute has the keyword or undefined 1st line of a block. In class Optgeom, that is 'OPTGEOM\n' and in Geom, that is title line. When it is empty, it should be set as ``None``.

#. Similarly, ``_block_ed`` attribute has the ending line of a block. In class BasisSet, by default that is 'ENDBS\n'. When it is empty, it should be set as ``None``.

#. When ``_block_bg`` and ``_block_ed`` are both ``''`` (note that is not ``None``), the whole block is regarded as empty and its data will not be updated.

#. The ``_block_data`` attribute stored the d12/d3 formatted text of the whole block.

#. The ``_block_dict`` attribute is a dictionary whose keys are CRYSTAL keywords and values are corresponding attributes. The sequence of key-value pairs should follow the requirement of d12/d3 files. For example, in Geom class, the 'OPTGEOM' keyword and corresponding attribute should always be placed after the 'SUPERCEL' keyword.

#. ``_block_key`` and ``_block_value`` are sorted, non-repeated lists of keys and values defined in ``_block_dict``.

Add keyword-related methods
---------------------------

Layer 1: ``Crystal_inputBASE`` or ``Properties_inputBASE`` Objects. 

Layer 2: 3 basic blocks, Geom, BasisSet and SCF. A BlockBASE class should be created.

Layer 3: Sub-blocks closed by END, such as 'OPTGEOM'. A BlockBASE class should be created. 

   Layer 2 and 3 should be created as the value of ``_block_attr`` attribute of the upper layer, with corresponding ``set_attr`` method and ``attr()`` method decorated by ``@property``

Layer 4: Keyword-like inputs.

   * Keywords with 'matrix-like' (ndimen\*ndimen) inputs, such as 'SUPERCEL'. Use ``set_matrix`` + ``assign_keyword``.
   * Keywords with 'list-like' (nline + a list of nline, or 1 line of multiple values) inputs, such as 'ATOMSPIN'. Use ``set_list`` + ``assign_keyword``.  
   * Keywords with 0 or 1 line input, such as 'TOLDEE' and 'CVOLOPT'. Use ``assign_keyword``.  
   * Other irregular keywords, such as 'SHRINK'. Design a special scheme and use ``assign_keyword``.  

Address conflicts
-----------------

To address conflicts between 2 'keyword-like' inputs, one can simply direct 2 different keywords to the same attribute. For example, both 'SUPERCEL' and 'SCELPHONO' keywords directs to the ``_sp_matrix`` attribute of Geom object. By setting either one the other is automatically covered.

To address conflicts involving 'block-like' inputs, the ``clean_conflict`` method should be used. Check the explanations below.

Planned developments
--------------------

#. Add Properties_inputBASE

#. Redo 'FIXINDEX' of SCF block. Make an individual class for GEOM / BASE / GEBA. The current implementation does not support GEBA.

#. Add CPHF / QHA / EOS blocks 

.. automodule:: CRYSTALpytools.base.inputbase
   :members:
   :private-members:
   :undoc-members:
   :show-inheritance:
