#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods of keywords used in 'properties' input file (d3).
"""
from CRYSTALpytools.base.inputbase import BlockBASE
import numpy as np


class PInputBlock(BlockBASE):
    """
    The base class of a 'property input block', i.e., a properties calculation,
    with non-repeated keywords.
    """
    def __init__(self):
        # The sequence of keywords should follow rules in the manual
        # Read inputbase.py for the definition of dict values
        dic = {
            # ---- Wavefunction ----
            'RDFMWF'    : [None, False, []],
            # ---- Band ----
            'BAND'      : [None, False, []], # BAND might destory NEWK
            'ANBD'      : [None, False, []],
            'BWIDTH'    : [None, False, []],
            # ---- preliminary ----
            'COMMENS'   : [None, False, []],
            'NOSYMADA'  : [None, False, []],
            'NEWK'      : [None, False, []],
            'PATO'      : [None, False, []],
            'PBAN'      : [None, False, []],
            'PGEOMW'    : [None, False, []],
            'PDIDE'     : [None, False, []],
            # 'PMP2'     : [None, False, []],
            # ---- A posteriori correlation ----
            'ADFT'      : ['_adft', True, [], 'propd3.ADFT()'],
            'ACOR'      : ['_adft', True, [], 'propd3.ADFT()'], # For reading only. The keyword will be changed to 'ADFT'
            'EDFT'      : ['_edft', True, [], 'propd3.EDFT()'],
            'ENECOR'    : ['_edft', True, [], 'propd3.EDFT()'], # For reading only. The keyword will be changed to 'EDFT'
            # ---- Charge and potential ----
            # 'PROPS2COMP'
            'CLAS'      : ['_clas', True, [], 'propd3.CLAS()'],
            'ECHG'      : ['_echg', True, [], 'propd3.ECHG()'],
            'ECH3'      : ['_ech3', True, [], 'propd3.ECH3()'],
            'POTM'      : ['_potm', True, [], 'propd3.POTM()'],
            'POT3'      : ['_pot3', True, [], 'propd3.POT3()'],
            'POTC'      : [None, False, []],
            'HIRSHCHG'  : [None, False, []],
            'PPAN'      : [None, False, []],
            # ---- Electron momentum density ----
            'EMDLDM'    : [None, False, []],
            # 'EMDPDM'   : [None, False, []],
            # 'EMDL'
            # 'EMDP'
            'KINETEMD'  : [None, False, []],
            # ---- Multiple moments ----
            'POLI'      : [None, False, []],
            'POLSPIN'   : [None, False, []],
            # ---- spin ----
            'ANISOTRO'  : [None, False, ['ISOTROPIC', 'ANISOTRO']],
            'ISOTROPIC' : [None, False, ['ISOTROPIC', 'ANISOTRO']],
            # 'SPINCNTM'
            # ---- DOSS-like ----
            'DOSS'      : [None, False, []],
            'COOP'      : [None, False, []],
            'COHP'      : [None, False, []],
            # ---- XRD ----
            # 'XFAC'
            # 'XRDSPEC'
            # ---- Compton Profiles ----
            # 'BIDIERD'  : [None, True, []], # developing
            # ---- Transport properties ----
            # 'BOLTZTRA'
            # ---- spontaneous polarization ----
            # 'SPOLBP'
            # 'SPOLWF'
            # ---- Wannierization ---
            # 'LOCALWF' : [None, True, []], # developing
            # ---- Optical dielectric constant ----
            # 'DIEL'
            # --- Topological analysis
            # 'TOPO'
            # ---- Print ----
            # 'FMWF'
            # 'INFOGUI'
            'PSCF'     : [None, False, []],
        }
        # Initialize the object to empty values
        super().__init__('', '', dic)

    # __call__ method inherited from BlockBASE

    # ---- Wavefunction ----

    def rdfmwf(self, rdfmwf=''):
        super().assign_keyword('RDFMWF', [], rdfmwf); return self

    # ---- Band ----

    def band(self, TITLE='', NLINE='', ISS='', NSUB='', INZB='', IFNB='',
             IPLO='', LPR66='', path=[]):
        """
        ``path`` can either be a NLINE\*2 list of string (labels), or a
        NLINE\*2\*3 list of int (fractional coordinates).
        """
        if path != []:
            shape = [1, 7]
            value = [TITLE, NLINE, ISS, NSUB, INZB, IFNB, IPLO, LPR66]
            if type(path[0][0]) == str: # labels
                for i in range(NLINE):
                    shape.extend([2,])
                    value.extend([path[i][0], path[i][1]])
            else: # values
                for i in range(NLINE):
                    shape.extend([6,])
                    value.extend(['{:<2d}'.format(path[i][0][0]),
                                  '{:<2d}'.format(path[i][0][1]),
                                  '{:<4d}'.format(path[i][0][2]),
                                  '{:<2d}'.format(path[i][1][0]),
                                  '{:<2d}'.format(path[i][1][1]),
                                  '{:<4d}'.format(path[i][1][2])])
            super().assign_keyword('BAND', shape, value); return self
        else:
            super().assign_keyword('BAND', [], TITLE); return self

    def anbd(self, NK='', NB='', TOL='', *args):
        if np.all(NK==None) or np.all(NK==''):
            super().assign_keyword('ANBD', [], NK); return self
        shape = [3,]
        value = [NK, NB, TOL]
        if NK > 0:
            shape.append(NK)
        if NB > 0:
            shape.append(NB)
        for arg in args:
            value.extend(arg)
        super().assign_keyword('ANBD', shape, value); return self

    def bwidth(self, INZB='', IFNB=''):
        super().assign_keyword('BWIDTH', [2,], [INZB, IFNB]); return self

    # ---- preliminary ----

    def commens(self, ICASO=''):
        super().assign_keyword('COMMENS', [1,], ICASO); return self

    def pato(self, IBN='', NPR=''):
        super().assign_keyword('PATO', [2,], [IBN, NPR]); return self

    def nosymada(self, nosymada=''):
        super().assign_keyword('NOSYMADA', [], nosymada); return self

    def newk(self, *args):
        """
        Inputs are separated by comma.
        """
        if len(args) == 0:
            args = ['']
        if np.all(args[0]==None) or np.all(args[0]==''):
            shape = []
            value = args[0]
        else:
            if len(args) == 4:
                shape = [2, 2]
                value = args
            elif len(args) == 7:
                shape = [2, 3, 2]
                value = args
            else:
                raise ValueError('Unknown k mesh definition')
        super().assign_keyword('NEWK', shape, value); return self

    def pban(self, NB='', nbands=[]):
        shape, value = super().set_list(NB, nbands)
        super().assign_keyword('PBAN', shape, value); return self

    def pgeomw(self, pgeomw=''):
        super().assign_keyword('RDFMWF', [], pgeomw); return self

    def pdide(self, EMI='', EMX=''):
        super().assign_keyword('PDIDE', [2,], [EMI, EMX]); return self

    # ---- A posteriori correlation ----

    @property
    def adft(self):
        """
        Subblock object ADFT.
        """
        self._adft._block_valid = True
        return self._adft

    @property
    def acor(self):
        """
        Alias of block ADFT.
        """
        return self.adft

    @property
    def edft(self):
        """
        Subblock object EDFT.
        """
        self._edft._block_valid = True
        return self._edft

    @property
    def enecor(self):
        """
        Alias of block EDFT.
        """
        return self.edft

    # ---- Charge and potential ----

    @property
    def clas(self):
        """
        Subblock object CLAS. Call ``self.clas()`` to set parameters. Call
        ``self.clas.coordina()`` (for example) to set mapnet keywords.

        Args:
            IDER (int): Order of the derivative
            IFOR (int): Number of point multipoles
            charge (list): Formal net charge
            NPY (int): Number of points (MAPNET). Default 100.
        """
        self._clas._block_valid = True
        return self._clas

    @property
    def echg(self):
        """
        Subblock object ECHG. Call ``self.echg()`` to set parameters. Call
        ``self.echg.coordina()`` (for example) to set mapnet keywords.

        Args:
            IDER (int): Order of the derivative
            NPY (int): Number of points (MAPNET). Default 100.
        """
        self._echg._block_valid = True
        return self._echg

    @property
    def ech3(self):
        """
        Subblock object ECH3. Call ``self.ech3()`` to set parameters. Call
        ``self.ech3.range()`` (for example) to set 3D grid keywords.

        Args:
            NP (int): Number of points along the first direction
        """
        self._ech3._block_valid = True
        return self._ech3

    @property
    def potm(self):
        """
        Subblock object POTM. Call ``self.potm()`` to set parameters. Call
        ``self.potm.coordina()`` (for example) to set mapnet keywords.

        Args:
            IDER (int): Order of the derivative
            ITOL (int): Penetration tolerance
            NPY (int): Number of points (MAPNET). Default 100.
        """
        self._potm._block_valid = True
        return self._potm

    @property
    def pot3(self):
        """
        Subblock object POT3. Call ``self.pot3()`` to set parameters. Call
        ``self.pot3.range()`` (for example) to set 3D grid keywords.

        Args:
            NP (int): Number of points along the first direction
            ITOL (int): Penetration tolerance
        """
        self._pot3._block_valid = True
        return self._pot3

    def potc(self, ICA='', NPU='', IPA='', datagrid=[]):
        """
        ``datagrid`` corresponds to X, Y, Z, ZP, ZM and ZD parameters, which
        can be NPU\*3 (``ICA=0``), 2\*1 (``ICA=1``) or 3\*1 (``ICA=2``) lists.
        """
        if np.all(ICA==None) or np.all(ICA==''):
            super().assign_keyword('POTC', [], ICA); return self

        shape = [3,]
        data = [ICA, NPU, IPA]
        datagrid = np.array(datagrid, dtype=float)
        if ICA == 0:
            if NPU > 0:
                for i in range(NPU):
                    shape.extend([3,])
                    data.extend(datagrid[i, :])
            else:
                pass
        elif ICA == 2:
            shape.extend([2,])
            data.extend([datagird[0], datagird[1]])
        elif ICA == 3:
            shape.extend([3,])
            data.extend([datagird[0], datagird[1], datagird[2]])

        super().assign_keyword('POTC', shape, data); return self

    def hirshchg(self, hchg=''):
        super().assign_keyword('HIRSHCHG', [], hchg); return self

    def ppan(self, pchg=''):
        super().assign_keyword('PPAN', [], pchg); return self


    # ---- Electron momentum density ----

    def emdldm(self, N='', PMAX='', STEP='', IPLO='', hkl=[], ICASO='', NSA1='', NSA2=''):
        """
        ``hkl`` is a N\*3 list of int.
        """
        shape = [4,]
        value = [N, PMAX, STEP, IPLO]
        for i in range(N):
            shape.append(3)
            value.extend(hkl[i])
        shape.extend([1, 2])
        value.extend([ICASO, NSA1, NSA2])
        super().assign_keyword('EMDLDM', shape, value); return self

    def kinetemd(self, PMAX='', PINT='', STEP='', STEPDIST='', ICASO=''):
        super().assign_keyword('KINETEMD', [4, 1], [PMAX, PINT, STEP, STEPDIST, ICASO]); return self

    # ---- Multiple moments ----

    def poli(self, IDIPO='', ITENS='', LPR68=''):
        super().assign_keyword('POLI', [1, 2], [IDIPO, ITENS, LPR68]); return self

    def polspin(self, IDIPO='', ITENS='', LPR68=''):
        super().assign_keyword('POLSPIN', [1, 2], [IDIPO, ITENS, LPR68]); return self

    # ---- spin ----

    def anisotro(self, keyword='', N='', IA=''):
        """
        Available keywords: 'ALL', 'UNIQUE', 'SELECT'
        """
        if keyword.upper() == 'SELECT':
            shape = [1, 1, N]
            value = ['SELECT', N] + list(IA[i])
        else:
            shape = [1]
            value = [keyword.upper()]
        super().assign_keyword('ANISOTRO', shape, value); return self

    def isotropic(self, keyword='', N='', IA=''):
        """
        Available keywords: 'ALL', 'UNIQUE', 'SELECT'
        """
        if keyword.upper() == 'SELECT':
            shape = [1, 1, N]
            value = ['SELECT', N] + list(IA[i])
        else:
            shape = [1]
            value = [keyword.upper()]
        super().assign_keyword('ISOTROPIC', shape, value); return self

    # ---- DOSS like ----

    def doss(self, NPRO='', NPT='', INZB='', IFNB='', IPLO='', NPOL='', NPR='', *args):
        """
        Energy range (BMI, BMA, 2\*1 list) and projections (N, NDM,
        Nproj\*(NDM+1) list) are defined by extendable list variables. If both
        are defined, projections should be the last entry.
        """
        shape, value = self._doss_like(NPRO, NPT, INZB, IFNB, IPLO, NPOL, NPR, args)
        super().assign_keyword('DOSS', shape, value); return self

    def coop(self, NINT='', NPT='', INZB='', IFNB='', IPLO='', NPOL='', NPR='', *args):
        """
        Energy range (BMI, BMA, 2\*1 list) and projections (N, NDM,
        Nproj\*(NDM+1) list) are defined by extendable list variables. If both
        are defined, projections should be the last entry.
        """
        shape, value = self._doss_like(NINT, NPT, INZB, IFNB, IPLO, NPOL, NPR, args)
        super().assign_keyword('COOP', shape, value); return self

    def cohp(self, NINT='', NPT='', INZB='', IFNB='', IPLO='', NPOL='', NPR='', *args):
        """
        Energy range (BMI, BMA, 2\*1 list) and projections (N, NDM,
        Nproj\*(NDM+1) list) are defined by extendable list variables. If both
        are defined, projections should be the last entry.
        """
        shape, value = self._doss_like(NINT, NPT, INZB, IFNB, IPLO, NPOL, NPR, args)
        super().assign_keyword('COHP', shape, value); return self

    # ---- Print ----

    def pscf(self, pscf=''):
        super().assign_keyword('PSCF', [], pscf); return self

    # ---- methods developed for properties ----
    @staticmethod
    def _doss_like(NPRO, NPT, INZB, IFNB, IPLO, NPOL, NPR, arg):
        """
        Get formatted input text for 'DOSS'-like keywords, including 'DOSS',
        'COHP' and 'COOP'.
        """
        if np.all(NPRO==None) or np.all(NPRO==''):
            return [], NPRO

        shape = [7,]
        value = [NPRO, NPT, INZB, IFNB, IPLO, NPOL, NPR]
        if INZB < 0 and IFNB < 0:
            shape.append(2)
            value.extend(arg[0])
        if NPRO > 0:
            if INZB < 0 and IFNB < 0:
                proj = arg[1]
            else:
                proj = arg[0]

            for line in proj:
                if abs(line[0]) != len(line) - 1:
                    raise ValueError(
                        'Inconsistent size of data. N AO/atom = {}, but {} entries are given. N AO/atom is changed.'.format(line[0], len(line)-1))
                    line[0] = len(line) - 1

                shape.append(len(line))
                value.extend(line)

        return shape, value


class Properties_inputBASE(PInputBlock):
    """
    The base class of Properties_input class. At maximum 5 same 'properties'
    calculations can be appended
    """
    def __init__(self):
        super().__init__()
        self._block_ed = 'END\n'
        # These are dummy keys, not printed out
        # Add 10 blocks
        for i in range(1, 6):
            key = 'APPEND{}'.format(i)
            value='_append{}'.format(i)
            self._block_dict.update({
                # Do not use the last element - it might cause error.
                # Call PInputBlock directly in this script
                key : [value, True, [], 'PInputBlock()'],
            })
            obj = PInputBlock()
            obj._block_valid = False
            setattr(self, value, obj)

        key = list(self._block_dict.keys())
        self._block_key = sorted(set(key), key=key.index)
        # Partition lines / keywords. Important for charge difference map
        self._block_label = ['PATO', 'PBAN', 'PDIDE']

    @property
    def append1(self):
        self._append1._block_valid = True
        return self._append1

    @property
    def append2(self):
        self._append2._block_valid = True
        return self._append2

    @property
    def append3(self):
        self._append3._block_valid = True
        return self._append3

    @property
    def append4(self):
        self._append4._block_valid = True
        return self._append4

    @property
    def append5(self):
        self._append5._block_valid = True
        return self._append5

    def analyze_text(self, data):
        """
        Analyze formatted d3 file.

        Args:
            data (str): Formatted string.
        """
        import copy

        nblock = 1
        inpkeys = [[]] # No repeated keyword in the same subblock
        repkeyline = [0,]
        divkeyline = [] # Partition lines / keywords

        datalines = data.strip().split('\n')
        for nline, line in enumerate(datalines):
            if line.upper() in self._block_key:
                key = line.upper()
                if key in self._block_label:
                    divkeyline.append(nline)
                if key in inpkeys[nblock-1]: # A new subblock is needed
                    nblock += 1
                    inpkeys.append([key,])
                    repkeyline.append(nline)
                else:
                    inpkeys[nblock-1].append(key)
            else:
                continue

        if nblock == 1: # No repeat
            super().analyze_text(data)
        elif nblock > 6: # Too many repeats
            raise Exception('''The same keyword has been repeated for over 5 times.
Please consider to split the calculation into multiples to reduce repeats.''')
        else:
            # No partition line
            if len(divkeyline) == 0:
                blockline = copy.deepcopy(repkeyline)
            # Insert partition line into repeated keywords.
            else:
                blockline = []
                for div in divkeyline:
                    for rep in repkeyline:
                        if div <= rep: # echg // PATO ECHG
                            blockline.append(div)
                        else: # echg // ECHG // PATO echg
                            blockline.append(rep)
            blockline.append(len(datalines))
            for i in range(nblock):
                text = ''.join('{}\n'.format(j) for j in datalines[blockline[i]:blockline[i+1]])
                if i == 0:
                    super().analyze_text(text)
                else:
                    obj = PInputBlock()
                    obj.analyze_text(text)
                    setattr(self, '_append{}'.format(i), obj)

        return self


class EDFT(BlockBASE):
    """
    EDFT sub-block
    """
    def __init__(self):
        bg = 'EDFT\n'
        ed = 'END\n'
        # The sequence of keywords should follow rules in the manual
        # Read inputbase.py for the definition of dict values
        dic = {
            'BECKE'    : [None, False, []],
            'SAVIN'    : [None, False, []],
            'RADIAL'   : [None, False, []],
            'ANGULAR'  : [None, False, []],
            'PRINT'    : [None, False, []],
            'PRINTOUT' : [None, False, []],
            'TOLLDENS' : [None, False, []],
            'TOLLGRID' : [None, False, []],
        }
        super().__init__(bg, ed, dic)

    # __call__ method inherited from BlockBASE

    def becke(self, becke=''):
        super().assign_keyword('BECKE', [], becke); return self

    def savin(self, savin=''):
        super().assign_keyword('SAVIN', [], savin); return self

    def radial(self, NR='', RL=[], IL=[]):
        if np.all(NR!=None) and np.all(NR!=''):
            RL = list(RL)
            IL = list(IL)
            if NR != len(RL) and NR != len(IL):
                raise ValueError('Inconsistent definition. NR is not equal to lengths of RL or IL')
        super().assign_keyword('RADIAL', [1, NR, NR], [NR,] + RL + IL); return self

    def angular(self, NI='', AL=[], IA=[]):
        if np.all(NI!=None) and np.all(NI!=''):
            AL = list(AL)
            IA = list(IA)
            if NI != len(AL) and NI != len(IA):
                raise ValueError('Inconsistent definition. NI is not equal to lengths of AL or IA')
        super().assign_keyword('ANGULAR', [1, NI, NI], [NI,] + AL + IA); return self

    def print(self, prt=''):
        super().assign_keyword('PRINT', [], prt); return self

    def printout(self, prt=''):
        super().assign_keyword('PRINTOUT', [], prt); return self

    def tolldens(self, ID=''):
        super().assign_keyword('TOLLDENS', [1,], ID); return self

    def tollgrid(self, IG=''):
        super().assign_keyword('TOLLGRID', [1,], IG); return self


class ADFT(EDFT):
    """
    ADFT sub-block
    """
    def __init__(self):
        bg = 'ADFT\n'
        ed = 'END\n'
        # The sequence of keywords should follow rules in the manual
        # Read inputbase.py for the definition of dict values
        dic = {
            'BECKE'    : [None, False, []],
            'SAVIN'    : [None, False, []],
            'RADIAL'   : [None, False, []],
            'ANGULAR'  : [None, False, []],
            'PRINT'    : [None, False, []],
            'PRINTOUT' : [None, False, []],
            'TOLLDENS' : [None, False, []],
            'TOLLGRID' : [None, False, []],
            'NEWBASIS' : ['_newbasis', True, [], 'crysd12.BasisSet()'],
        }
        self = BlockBASE(bg, ed, dic)

    # __call__ method inherited from BlockBASE

    @property
    def newbasis(self):
        """
        Subblock object NEWBASIS
        """
        self._newbasis._block_valid = True
        self._newbasis._block_bg = 'NEWBASIS\n'
        self._newbasis._block_ed = 'END\n'
        return self._newbasis


class Grid2D(BlockBASE):
    """
    A base class for 2D data grid (ECHG, POTM, CLAS).

    Args:
        header (str): Formatted headers
    """
    def __init__(self, header=''):
        bg = header
        ed = 'END\n'
        # The sequence of keywords should follow rules in the manual
        # Read inputbase.py for the definition of dict values
        dic = {
            'COORDINA' : [None, False, ['COORDINA', 'ATOMS']],
            'ATOMS'    : [None, False, ['COORDINA', 'ATOMS']],
            'RECTANGU' : [None, False, []],
            'MARGINS'  : [None, False, []],
            'PRINT'    : [None, False, []],
            'ANGSTROM' : [None, False, []],
            'BOHR'     : [None, False, []],
            'FRACTION' : [None, False, []],
        }
        super().__init__(bg, ed, dic)

    # __call__ method inherited from BlockBASE

    def coordina(self, crda='', crdb='', crdc=''):
        if np.all(crda==None) or np.all(crda==''):
            shape = []
            value = crda
        else:
            shape = [1, 1, 1]
            value = [''.join('{:< 10.6f}'.format(i) for i in crda),
                     ''.join('{:< 10.6f}'.format(i) for i in crdb),
                     ''.join('{:< 10.6f}'.format(i) for i in crdc)]
        super().assign_keyword('COORDINA', shape, value); return self

    def atoms(self, IA='', IB='', IC=''):
        """
        IA IB IC are 4\*1 lists
        """
        super().assign_keyword('ATOMS', [4, 4, 4], [IA, IB, IC]); return self

    def rectangu(self, rec=''):
        super().assign_keyword('RECTANGU', [], rec); return self

    def margins(self, ABM='', CDM='', ADM='', BCM=''):
        super().assign_keyword('MARGINS', [4,], [ABM, CDM, ADM, BCM]); return self

    def print(self, prt=''):
        super().assign_keyword('PRINT', [], prt); return self

    def angstrom(self, ang=''):
        super().assign_keyword('ANGSTROM', [], ang); return self

    def bohr(self, br=''):
        super().assign_keyword('BOHR', [], br); return self

    def fraction(self, frac=''):
        super().assign_keyword('FRACTION', [], frac); return self


class ECHG(Grid2D):
    """
    Class for ECHG. Initialization can be done by formatted string only.
    Otherwise call this object from upper block.

    Args:
        IDER (int): Order of the derivative
        NPY (int): Number of points (MAPNET)
    """
    def __call__(self, IDER='', NPY=100):
        if type(IDER) == str:
            super().__init__('ECHG\n')
            if np.all(IDER!=''):
                self.analyze_text(IDER)
        elif np.all(IDER==None):
            super().__init__('ECHG\n')
            self._block_valid = False
        elif type(IDER) == type(self):
            self = IDER
        elif type(IDER) == int: # valid data
            self.__init__('ECHG\n{}\n{}\n'.format(IDER, NPY))
        else:
            raise ValueError('Unknown input type')


class POTM(Grid2D):
    """
    Class for POTM. Initialization can be done by formatted string only.
    Otherwise call this object from upper block.

    Args:
        IDER (int): Order of the derivative
        ITOL (int): Penetration tolerance
        NPY (int): Number of points (MAPNET)
    """
    def __call__(self, IDER='', ITOL=5, NPY=100):
        if type(IDER) == str:
            super().__init__('POTM\n')
            if np.all(IDER!=''):
                self.analyze_text(IDER)
        elif np.all(IDER==None):
            super().__init__('POTM\n')
            self._block_valid = False
        elif type(IDER) == type(self):
            self = IDER
        elif type(IDER) == int: # valid data
            super().__init__('POTM\n{}{:4d}\n{}\n'.format(IDER, ITOL, NPY))
        else:
            raise ValueError('Unknown input type')


class CLAS(Grid2D):
    """
    Class for CLAS. Initialization can be done by formatted string only.
    Otherwise call this object from upper block.

    Args:
        IDER (int): Order of the derivative
        IFOR (int): Number of point multipoles
        charge (list): Formal net charge
        NPY (int): Number of points (MAPNET)
    """
    def __call__(self, IDER='', IFOR='', charge=[], NPY=100):
        if type(IDER) == str:
            super().__init__('CLAS\n')
            if np.all(IDER!=''):
                self.analyze_text(IDER)
        elif np.all(IDER==None):
            super().__init__('CLAS\n')
            self._block_valid = False
        elif type(IDER) == type(self):
            self = IDER
        elif type(IDER) == int: # valid data
            if IFOR == 0:
                super().__init__('CLAS\n{:<2d}{:2d}\n{}\n'.format(IDER, IFOR, NPY))
            elif IFOR == 1:
                string = ''.join('{:<10.6f} '.format(i) for i in charge)
                super().__init__('CLAS\n{:<2d}{:2d}\n{}\n{}\n'.format(IDER, IFOR, string, NPY))
            else:
                raise ValueError('Unknon IFOR value. Check your input.')
        else:
            raise ValueError('Unknown input type')


class Grid3DBASE(BlockBASE):
    """
    A base class for 3D data grid (ECH3, POT3). Not for 'GRID3D' keyword

    Args:
        header (str): Formatted headers
    """
    def __init__(self, header=''):
        bg = header
        ed = ''
        # The sequence of keywords should follow rules in the manual
        # Read inputbase.py for the definition of dict values
        dic = {
            'RANGE' : [None, False, ['RANGE', 'SCALE']],
            'SCALE' : [None, False, ['RANGE', 'SCALE']],
        }
        super().__init__(bg, ed, dic)

    def scale(self, scale1='', scale2='', scale3=''):
        if np.all(scale1=='') or np.all(scale1==None): # others
            super().assign_keyword('SCALE', [], scale1)
        elif np.all(scale2==''): # 2D
            super().assign_keyword('SCALE', [1,], scale1)
        elif np.all(scale3==''): # 1D
            super().assign_keyword('SCALE', [2,], [scale1, scale2])
        else: # 0D
            super().assign_keyword('SCALE', [3,], [scale1, scale2, scale3])
        return self

    def range(self, rg1='', rg2='', rg3=''):
        """
        Inputs are 2\*1 lists. ``[minvalue, maxvalue]``. The sequence of rg123
        is consistent with CRYSTAL (2D: z; 1D: y, z; 0D: x, y, z)
        """
        if np.all(rg1=='') or np.all(rg1==None): # others
            super().assign_keyword('RANGE', [], rg1)
        elif np.all(rg2==''): # 2D
            super().assign_keyword('RANGE', [1, 1], [rg1[0], rg1[1]])
        elif np.all(rg3==''): # 1D
            super().assign_keyword('RANGE', [2, 2], [rg1[0], rg2[0], rg1[1], rg2[1]])
        else: # 0D
            super().assign_keyword('RANGE', [3, 3], [rg1[0], rg2[0], rg3[0], rg1[1], rg2[1], rg3[1]])
        return self


class ECH3(Grid3DBASE):
    """
    Class for ECH3. Initialization can be done by formatted string only.
    Otherwise call this object from upper block.

    Args:
        NP (int): Number of points along the first direction
    """
    def __call__(self, NP=''):
        if type(NP) == str:
            super().__init__('ECH3\n')
            if np.all(NP!=''):
                self.analyze_text(NP)
        elif np.all(NP==None):
            super().__init__('ECH3\n')
            self._block_valid = False
        elif type(NP) == type(self):
            self = NP
        elif type(NP) == int: # valid data
            super().__init__('ECH3\n{}\n'.format(NP))
        else:
            raise ValueError('Unknown input type')


class POT3(Grid3DBASE):
    """
    Class for POT3. Initialization can be done by formatted string only.
    Otherwise call this object from upper block.

    Args:
        NP (int): Number of points along the first direction
        ITOL (int): Penetration tolerance
    """
    def __call__(self, NP='', ITOL=5):
        if type(NP) == str:
            super().__init__('POT3\n')
            if np.all(NP!=''):
                self.analyze_text(NP)
        elif np.all(NP==None):
            super().__init__('POT3\n')
            self._block_valid = False
        elif type(NP) == type(self):
            self = NP
        elif type(NP) == int: # valid data
            super().__init__('POT3\n{}{:4d}\n'.format(NP, ITOL))
        else:
            raise ValueError('Unknown input type')

