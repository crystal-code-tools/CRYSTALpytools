CRYSTALpytools.base.inputbase module
====================================

.. _ref-base-inputbase:

Python wrapper for CRYSTAL inputs. Provides basic classes, methods and attributes to read, operate and write CRYSTAL d12/d3 files.

Structure of classes
--------------------

#. All the information are stored in the dictionary ``obj._block_dict``, with the keywords as keys and lists as values. See below for definitions of list elements. The sequence of key-value pairs should follow the requirement of d12/d3 files. For example, in ``Geom``, the 'OPTGEOM' key-value pair should always be placed after the 'SUPERCEL' key-value pair.

#. Keywords in a block are NON-REPEATABLE. If the same keyword can appear multiple times, please consider to divide the input into multiple subblocks. Both :ref:`Crystal_inputBASE <ref-base-crysd12>` and :ref:`Properties_inputBASE <ref-base-propd3>` use 'dummpy' subblocks. For details please refer to docs there.

#. Keywords are called as methods. For 'keyword-like' keywords (keywords not closed by 'END'), calling the method (``obj.keyword(value)``) modifies the value in ``obj._block_dict``. For 'block-like' keywords (keyword closed by 'END'), calling ``obj.keyword()`` modifies the sub-block attribute ``obj._keyword``, which is a derivative of  ``BlockBASE``. ``obj.keyword`` returns to this attribute.

#. For ``BlockBASE`` objects, ``_block_bg`` attribute has the keyword or title line of a block. For example, ``Optgeom`` has 'OPTGEOM\n' and ``Geom`` has title line.

#. Similarly, ``_block_ed`` attribute has the ending line of a block. In ``BasisSet``, by default that is 'ENDBS\n'. ``''`` is automatically set for ``_block_ed`` if 'BASISSET' keyword is used.

#. The ``_block_data`` attribute is the formatted CRYSTAL input string, which can be called and updated by the 'data' attriute ``obj.data``.

#. ``_block_key`` is the sorted, non-repeated lists of keys defined in ``_block_dict``.

#. ``_block_valid`` is a boolian to indicate whether to print out the data of this block.

Add keyword-related methods
---------------------------

Layer 1: ``Crystal_inputBASE`` or ``Properties_inputBASE`` Objects. 

Layer 2: All of them are inherited from ``BlockBASE``. For d12, 3 basic blocks, ``geom``, ``basisSet`` and ``scf``. For d3, optional, ``append1`` to ``append5``, useful only when the same keyword in main block (not protected by sub-blocks) is set for multiple times, for example, 'ECHG' + 'PATO' + 'ECHG'.

Layer 3: Sub-blocks closed by END, such as 'OPTGEOM'. Inherited from ``BlockBASE``.

Layer 4: Keyword-like inputs or sub-sub-blocks.

   * Keywords with 'matrix-like' (ndimen\*ndimen) inputs, such as 'SUPERCEL'. Use ``set_matrix`` + ``assign_keyword``.
   * Keywords with 'list-like' (nline + a list of nline, or 1 line of multiple values) inputs, such as 'ATOMSPIN'. Use ``set_list`` + ``assign_keyword``.  
   * Keywords with 0 or 1 line input, such as 'TOLDEE' and 'CVOLOPT'. Use ``assign_keyword``.  
   * Other irregular keywords, such as 'SHRINK'. Design a special scheme and use ``assign_keyword``.  
   * Sub-sub-blocks, same as sub-blocks. Change ``_block_bg`` and ``_block_ed`` accordingly.  

All the information, including keywords and subblocks, are stored as a dictionary under ``obj._block_dict`` attribute. The key is CRYSTAL keyword (not all the keywords are supported). The value is a list.

#. The first element
    * For 'keyword-like' keywords ``None``, ``''`` or value(string). ``None`` means skipping the keyword and corresponding values. ``''`` prints out keyword. ``key\nvalue`` is printed out for other formatted string values.
    * For 'block-like' keywords, the first element must not be edited. It is a string that directs to the real attribute of that sub-block, such as ``_optgeom``. Using attribute operation functions and text saved the developer can visit the 'real' attribute without triggering potential commands defined in ``@property`` decorated ``keyword(self)`` functions.

#. The second element: boolian. Whether the keyword requires a block object.

#. The third element: A list of conflicting keywords. See the next section.

#. The fourth element: Valid only if the second element is ``True``. A string that initializes the subblock object. For example, in ``Geom()._block_dict``, the fourth element of the list corresponding to 'OPTGEOM' key is 'crysd12.Optgeom()'. This helps the code to initialize sub-block objects automatically by ``eval()`` command. All the sub-blocks are initialized with ``_block_valid = False`` when the block object is initialized.

Address conflicts
-----------------

To address conflicts of keywords, the ``clean_conflict`` method is automatically called. A list of conflicting keywords should be given. The list can include the input keyword itself. Conflicting keywords, if they are 'keyword-like', the first elements of their 'value list' in ``_block_dict`` are set to ``None``; if they are 'block-like', they are re-initialized with ``_block_valid = False``.

Other comments
--------------------

#. Not all the keywords in CRYSTAL manual have been implemented. Calling the non-existing attributes or methods causes problems. Please raise an issue in GitHub repo to contact the authors if you need any specific keyword immidiately.

#. Please also read documentations and examples for :ref:`Crystal_inputBASE <ref-base-crysd12>`, :ref:`Properties_inputBASE <ref-base-propd3>` and :ref:`crystal_io <ref-crystalio>`.

.. automodule:: CRYSTALpytools.base.inputbase
   :members:
   :private-members:
   :undoc-members:
   :show-inheritance:
