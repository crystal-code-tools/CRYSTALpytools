#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base object of all the input (d12/d3) blocks.
"""
class BlockBASE():
    """
    The base class of 'block' objects
    """

    def __init__(self):
        self._block_bg = ''
        self._block_ed = ''
        self._block_data = ''
        self._block_dict = {}
        key = list(self._block_dict.keys())
        attr = list(self._block_dict.values())
        self._block_key = sorted(set(key), key=key.index)
        self._block_attr = sorted(set(attr), key=attr.index)

    @property
    def data(self):
        """
        Settings in all the attributes are summarized here.
        """
        self.update_block()
        text = ''
        for i in [self._block_bg, self._block_data, self._block_ed]:
            if i == None:
                continue
            text += i
        return text

    @staticmethod
    def assign_keyword(key, shape, value=None):
        """
        Transform value into string formats.

        Args:
            key (str): CRYSTAL keyword
            shape (list[int]): 1D list. Shape of input text. Length: Number of
                lines; Element: Number of values
            value (list | str): List, a 1D list of arguments; ``None`` or a list
                begins with ``None``, return to keyword only; ``''`` or a list begins
                with ``''``, Clean everything

        Returns:
            text (str): CRYSTAL input
        """
        if type(value) != list and type(value) != tuple:
            value = [value, ]

        # Keyword only : Value is None and and key is not ''
        if value[0] == None and key != '':
            return '{}\n'.format(key)

        # Clean everything : Empty key or value is ''
        if value[0] == '' or key == '':
            return ''

        # Wrong input: Number of args defined by shape != input.
        if sum(shape) != len(value):
            raise ValueError(
                "The number of input parameters '{}' does not meet requirements.".format(value))

        # Correct number of args and valid key. Key = None, no keyword
        if key != None:
            text = '{}\n'.format(key)
        else:
            text = ''
        value_counter = 0
        for nvalue in shape:
            for v in value[value_counter:value_counter + nvalue]:
                text += '{} '.format(v)
            text += '\n'
            value_counter += nvalue

        return text

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
        import numpy as np

        if mx == '':  # Clean data
            return [], ''
        elif mx == None:  # Keyword only
            return [], None

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
        if args[0] == '':  # Clean data
            return [], ''
        elif args[0] == None:  # Keyword only
            return [], None

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

    def clean_conflict(self, newattr, conflict):
        """
        Addressing the conflictions between attributes, usually between blocks
        or block and keywords. For conflictions between keywords, they are set
        to direct the same attribute.

        Args:
            newattr (str): The attribute explicitly specified.
            conflict (list[str]): The list of conflicting attributes including
                the called one. 'Real' attributes (begins with '_') are needed.
        """
        for cttr in conflict:
            if cttr == newattr:
                continue

            if hasattr(self, cttr):
                if '_block' in cttr:  # Block conflicts
                    obj = getattr(self, cttr)
                    obj.clean_block()
                    setattr(self, cttr, obj)
                else:  # Keyword conflict
                    setattr(self, cttr, '')
        return

    def clean_block(self):
        """
        Clean all the keyword-related attributes (accessible attributes).

        .. note::

            This method directly deletes all the attributes. Alternatively, by
            setting an attribute with '', the attribute is kept but its old
            values are erased.

        """
        self._block_bg = ''
        self._block_ed = ''
        self._block_data = ''
        for a in self._block_attr:
            try:
                delattr(self, a)
            except AttributeError:
                continue
        return

    def update_block(self):
        """
        Update the ``_block_data`` attribute: Summarizing all the settings to
        ``_block_data`` attribute for inspection and print
        """
        self._block_data = ''
        if self._block_bg != '' or self._block_ed != '':
            for attr in self._block_attr:
                if attr[0] == '_':  # Keyword-like attributes
                    if hasattr(self, attr):
                        self._block_data += getattr(self, attr)
                else:  # Block-like attributes
                    # If sub-block does not exist, call @property name will create a new one. Call _block_name instead
                    attr_real = '_block_' + attr
                    if hasattr(self, attr_real):
                        obj = getattr(self, attr_real)
                        obj.update_block()
                        for i in [obj._block_bg, obj._block_data, obj._block_ed]:
                            if i == None:
                                continue
                            self._block_data += i
        return

    def analyze_text(self, text):
        """
        Analyze the input text and return to corresponding attributes
        """
        import warnings

        if self._block_ed == None:
            end_block_label = ''
        else:
            end_block_label = self._block_ed

        textline = text.strip().split('\n')
        attr = ''
        attr_real = ''
        value = ''
        for idx, t in enumerate(textline[::-1]):
            if t.upper() in self._block_key:  # Keyword line: ending point
                value += t + '\n'
                attr = self._block_dict[t]
                value = '\n'.join([i for i in value.strip().split('\n')[::-1]]) + '\n' # Reverse the string
                if attr[0] == '_': # Keyword-like attributes
                    attr_real = attr
                    if hasattr(self, attr_real):
                        warnings.warn("Keyword same as or equivalent to '{}' exists. The new entry will cover the old one".format(t),
                                      stacklevel=2)
                    setattr(self, attr, value)
                else:  # Block-like attributes
                    attr_real = '_block_' + attr
                    # If sub-block does not exist, call @property name will create a new one. Call _block_name instead
                    if hasattr(self, attr_real):
                        warnings.warn("Keyword same as or equivalent to '{}' exists. The new entry will cover the old one".format(t),
                                      stacklevel=2)
                    # This step will create an obj if the real attribute does not exist
                    obj = getattr(self, attr)
                    obj.analyze_text(value)
                    # @property does not have setter
                    setattr(self, attr_real, obj)

                # Clean values
                value = ''
            elif t in end_block_label: # End line: starting point
                continue
            else:
                value += t + '\n'
        # Last lines if unallocated string exists
        if value != '':
            value = '\n'.join([i for i in value.strip().split('\n')[::-1]]) + '\n' # Reverse the string
            value = value[1:]
            self._block_data = value + self._block_data

        return
