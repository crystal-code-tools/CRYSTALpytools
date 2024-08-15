#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base object of all the input (d12/d3) blocks.
"""
import numpy as np


class BlockBASE():
    """
    The base class of 'block' objects

    Args:
        bg (str): Beginning line. Keyword or title
        ed (str): Ending line.
        dic (dict): Keyword and value (in text) pairs. Format is listed below.

    Returns:
        self (BlockBASE): New attributes listed below
        self._block_bg (str): Keyword of the block
        self._block_ed (str): End of block indicator, 'END'
        self._block_data (str) : Formatted text string for d12.
        self._block_dict (dict): Keyword and value (in text) pairs. Key: CRYSTAL
            keyword in string. Value: A 3\*1 or 4\*1 list, see below.
            * 1st element: For keywords: String, formatted output; ``None``, not including the keyword; ``''``, keyword only
            * 1st element: For subblocks: Name of the 'real' attribute (protected by property decorator)
            * 2nd: bool, whether it is a sub-block or not;
            * 3rd: A list of conflicting keywords.
            * 4th: Subblock only. String that initializes the subblock object.
        self._block_key (list): Allowed keyword list
        self._block_valid (bool): Whether this block is valid and for print.
    """

    def __init__(self, bg, ed, dic):
        from CRYSTALpytools.base import crysd12, propd3

        self._block_bg = bg # Title or keyword of the block
        self._block_ed = ed # End of block indicator
        self._block_data = '' # Formatted text string
        self._block_dict = dic # Data
        self._block_valid = True # Whether this block is valid and for print

        key = list(self._block_dict.keys())
        self._block_key = sorted(set(key), key=key.index)
        for k in key:
            if self._block_dict[k][1] == True: # Initialize all the subblock objs
                obj = eval(self._block_dict[k][3])
                obj._block_valid = False
                setattr(self, self._block_dict[k][0], obj)

    def __call__(self, obj=''):
        if type(obj) == str:
            self.__init__()
            if np.all(obj!=''):
                self.analyze_text(obj)
        elif np.all(obj==None):
            self.__init__()
            self._block_valid = False
        elif type(obj) == type(self):
            self = obj
        else:
            raise ValueError('Unknown data type.')

    @property
    def data(self):
        """
        Settings in all the attributes are summarized here.
        """
        import warnings

        if self._block_valid == False:
            warnings.warn("This block is not visible. Set 'self._block_valid = True' to get data",
                          stacklevel=2)
            return ''

        self.update_block()
        text = ''
        for i in [self._block_bg, self._block_data, self._block_ed]:
            text += i
        return text

    def assign_keyword(self, key, shape, value=''):
        """
        Transform value into string formats.

        Args:
            shape (list[int]): 1D list. Shape of input text. Length: Number of
                lines; Element: Number of values
            value (list | str): List, a 1D list of arguments; ``''`` or a list
                begins with ``''``, return to ``''``; ``None`` or a list
                begins with ``None``, return ``None``.

        Returns:
            text (str): CRYSTAL input
        """
        # Check the validity of key
        if key not in self._block_key:
            raise ValueError("Cannot recognize keyword '{}'.".format(key))

        self.clean_conflict(key)

        if type(value) != list and type(value) != tuple:
            value = [value, ]

        # Keyword only or cleanup: Value is '' or None
        if np.all(value[0]=='') or np.all(value[0]==None):
            self._block_dict[key][0] = value[0]
            return self

        # Wrong input: Number of args defined by shape != input.
        if sum(shape) != len(value):
            raise ValueError(
                "Value does not meet requirement of keyword '{}'.".format(key))

        # Correct number of args
        text = ''
        value_counter = 0
        for nvalue in shape:
            for v in value[value_counter:value_counter + nvalue]:
                text += '{} '.format(v)
            text += '\n'
            value_counter += nvalue

        self._block_dict[key][0] = text
        return self

    @staticmethod
    def set_matrix(mx):
        """
        Set matrix-like data to get assign_keyword inputs. Used for supercell 
        expansion matrix and strain tensor.

        Args:
            mx (list | str): ``ndimen*ndimen`` list, ``None``, or ``''``

        Returns:
            shape (list): ``ndimen*1`` 1D list. All elements are ndimen.
            value (list): ``ndimen*2*1`` 1D list. Flattened matrix.
        """
        if np.all(mx==None):  # Clean data
            return [], None
        elif np.all(mx==''):  # Keyword only
            return [], ''

        matrix = np.array(mx)
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Input matrix is not a square matrix.")

        shape = [matrix.shape[0] for i in range(matrix.shape[0])]
        value = matrix.reshape([1, -1]).tolist()[0]

        return shape, value

    @staticmethod
    def set_list(*args):
        """
        Set list-like data to get assign_keyword inputs. Used for lists with
        known dimensions. Such as atom coordinate list.

        Args:
            \*args : ``''``, Clean data; ``None``, Return keyword only; 
                ``int, list``, int for length of the list, list for list data
        Returns:
            shape (list): 1 + length 1D list or []
            args (list): Flattened list, [] or ''
        """
        if np.all(args[0]==None):  # Clean data
            return [], None
        elif np.all(args[0]==''):  # Keyword only
            return [], ''

        if len(args) != 2 or int(args[0]) != len(args[1]):
            return ValueError('Input format error. Arguments should be int + list')

        shape = [1, ]
        value = [int(args[0]), ]

        if type(args[1][0]) == list or type(args[1][0]) == tuple:  # 2D list (multi-rows)
            for i in args[1]:
                shape += [len(i), ]
                value += i
        else:  # 1D list (single row)
            shape += [len(args[1]), ]
            value += args[1]

        return shape, value

    def clean_conflict(self, key):
        """
        Addressing the conflictions between keywords.

        Args:
            key (str): The keyword specified.
        """
        import warnings
        for cttr in self._block_dict[key][2]:
            cttr = cttr.upper()
            if cttr == key:
                continue

            if self._block_dict[cttr][1] == False: # keyword
                if np.all(self._block_dict[cttr][0]!=None):
                    warnings.warn("'{}' conflicts with the existing '{}'. The old one is deleted.".format(key, cttr),
                                  stacklevel=3)
                    self._block_dict[cttr][0] = None
            else: # subblock
                obj = getattr(self, self._block_dict[cttr][0])
                if obj._block_valid == True:
                    warnings.warn("'{}' conflicts with the existing '{}'. The old one is deleted.".format(key, cttr),
                                  stacklevel=3)
                    obj(None)
                    setattr(self, self._block_dict[cttr][0], obj)

        return self

    def update_block(self):
        """
        Update the ``_block_data`` attribute: Summarizing all the settings to
        ``_block_data`` attribute for inspection and print
        """
        self._block_data = ''
        for key in self._block_key:
            if self._block_dict[key][1] == False: # Keyword-like attributes
                if np.all(self._block_dict[key][0]!=None):
                    self._block_data += key + '\n' + self._block_dict[key][0]
            else: # Block-like attributes, get data from the corresponding attribute
                # It is important to use real attribute here for subblocks
                # To avoid automatically setting objreal._block_valid == True
                objreal = getattr(self, self._block_dict[key][0])
                if objreal._block_valid == True:
                    # It is important to use @property decorated attribute here for sub-sub-blocks
                    # some of them are the same as subblocks but different keywords
                    # The modification method is a decorated attribute in its upper block (self)
                    obj = getattr(self, key.lower())
                    self._block_data += obj.data # update and print subblock
                    setattr(self, self._block_dict[key][0], obj)
        return self

    def analyze_text(self, text):
        """
        Analyze the input text and return to corresponding attributes
        """
        import warnings
        from CRYSTALpytools.base import crysd12, propd3

        if np.all(self._block_ed==None):
            end_block_label = ''
        else:
            end_block_label = 'END'

        textline = text.strip().split('\n')
        attr = ''
        attr_real = ''
        value = ''
        for idx, t in enumerate(textline[::-1]):
            if t.upper() in self._block_key:  # Keyword line: ending point
                t = t.upper()
                if self._block_dict[t][1] == False: # Keyword-like attributes
                    if np.all(self._block_dict[t][0]!=None):
                        warnings.warn("Keyword '{}' exists. The new entry will cover the old one".format(t),
                                      stacklevel=2)
                    self._block_dict[t][0] = value
                else:  # Block-like attributes
                    obj = getattr(self, self._block_dict[t][0])
                    if obj._block_valid == True:
                        warnings.warn("Keyword '{}' exists. The new entry will cover the old one".format(t),
                                      stacklevel=2)
                    # Update obj
                    obj(value)
                    setattr(self, self._block_dict[t][0], obj)
                # Clean values
                value = ''
            elif end_block_label in t.upper(): # End line: starting point
                continue
            else:
                value = t + '\n' + value
        # Last lines if unallocated string exists
        if np.all(value!=''):
            if self._block_bg == value: # for subblocks
                pass
            else: # saved into beginning lines
                textline = value.strip().split('\n')
                for t in textline:
                    if self._block_bg == t + '\n':
                        pass
                    else:
                        self._block_bg = self._block_bg + t + '\n'
        return self
