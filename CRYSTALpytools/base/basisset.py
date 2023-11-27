class BS_GTF():
    """
    A Gaussian type function object. Random Text

    Args:
        info (str): GTF radial coefficients in CRYSTAL format
    """

    def __init__(self, info):
        info = info.replace('D', 'E').split()
        if len(info) == 2:
            self.exp = eval(info[0])
            self.contr = eval(info[1])
            self.pcoef = None
        elif len(info) == 3:
            self.exp = eval(info[0])
            self.contr = eval(info[1])
            self.pcoef = eval(info[2])
        else:
            raise ValueError('Input format error.')

    @property
    def data(self):
        """
        Print basis set information in CRYSTAL format
        """
        return self._print_data()

    def _print_data(self):
        """
        Return to a string in CRYSTAL format
        """
        if self.pcoef == None:
            bs_str = '  {: 16.10E}    {: 16.10E}\n'.format(
                self.exp, self.contr)
        else:
            bs_str = '  {: 16.10E}    {: 16.10E}    {: 16.10E}\n'.format(
                self.exp, self.contr, self.pcoef)

        return bs_str


class BS_Shell():
    """
    The Shell object to set a general basis set.

    Args:
        info (list[str]): String in Crystal format

    Other arguments are consistent with CRYSTAL manual.
    """

    def __init__(self, ITYB, LAT, NG, CHE, SCAL, info):
        from CRYSTALpytools.base.basisset import BS_GTF

        if ITYB != 0:
            raise ValueError(
                'BS_Shell class is limited to general definitions of basis set.')

        self.type = ITYB
        self.angular = LAT
        self.ngtf = NG
        self.chg = CHE
        self.scale = SCAL
        self.gtf = []
        for line in info:
            self.gtf.append(BS_GTF(line))

    @property
    def data(self):
        """
        Print basis set information in CRYSTAL format
        """
        return self._print_data()

    def _print_data(self):
        """
        Return to a string in CRYSTAL format
        """
        bs_str = ''
        bs_str += format('%-2i%-2i%-3i%-6.2f%-6.2f\n' %
                         (self.type, self.angular, self.ngtf, self.chg, self.scale))
        for g in self.gtf:
            bs_str += g.data

        return bs_str


class BS_Atom():
    """
    The atom object to set a basis set.

    .. note::

        The type of basis set should be consistent for all the shells.
        Otherwise unexpected errors might happen.

        Free effective core pseudopotential (ECP) definition not supported.

    Args:
        z (int): Convential atomic number Z.
        info (list(str)): String in Crystal BS format
    """

    def __init__(self, z, info):
        import re
        from CRYSTALpytools.base.basisset import BS_Shell

        if z > 99:
            self.z = z % 100
        else:
            self.z = z
        self.conventional_atomic_number(z)
        self.nshell = 0
        self.shell = []
        self.ecp = None
        set_stdbs = False
        info.append('')  # Avoid index out of range error
        for idx, line in enumerate(info):
            if re.match(r'^[0-2]\s+[0-5]\s+[0-9]+\s+[0-9,\.,\-,\+]+\s+[0-9,\.,\-,\+]+$', line):
                self.nshell += 1
                data = line.split()
                ITYB = int(data[0])
                LAT = int(data[1])
                NG = int(data[2])
                CHE = float(data[3])
                SCAL = float(data[4])
                if ITYB == 0:  # General definition of basis sets
                    self.shell.append(
                        BS_Shell(ITYB, LAT, NG, CHE, SCAL,
                                 info[idx + 1: idx + NG + 1])
                    )
                else:  # Standard STO-3G / 3(6)-21G* basis sets
                    set_stdbs = True
                    break
            elif re.match(r'^[A-Z]+', line): # ECP line
                if re.match(r'^INPUT', line):
                    raise ValueError('Free-format ECP not supported')
                self.ecp = line
            else:
                continue

        if set_stdbs == True:
            self._set_stdbs(info)

    def conventional_atomic_number(self, zconv):
        """
        Set convential atomic number.

        Args:
            zconv (int): Z(real atomic number) + n\*100. N is an integer.
        """
        if zconv % 100 != self.z:
            raise ValueError(
                'Convential atomic number should be Z + n*100. N is an integer.')

        self.zconv = zconv
        return

    def _set_stdbs(self, z, info):
        """
        Set the standard STO-3G / 3(6)-21G* basis sets. All the settings will
        be converted to general BS definitions.

        Args:
            info (list(str)): String in Crystal BS format
        """
        import warnings
        import copy
        import basis_set_exchange as bse
        from CRYSTALpytools.base.basisset import Basisset

        warnings.warn('''The built-in standard STO-3G / 3(6)-21G basis set is used.
CRYSTALpytools will substract basis set definition from Basis Set Exchange.
That might lead to discripencies in basis set definitions.''')

        # set a general 1-atom BS template
        ITYB_init = info[1].split()[0]
        NG_init = info[1].split()[2]
        if ITYB_init == 1:  # STO-3G
            bs = Basisset('STO-3G', [self.z, ])
        elif ITYB_init == 2 and NG_init == 3:  # 3-21G
            bs = Basisset('3-21G', [self.z, ])
        elif ITYB_init == 2 and NG_init == 6:  # 6-21G
            bs = Basisset('6-21G', [self.z, ])
        else:
            raise ValueError('Basis set input format error.')

        atom = bs.atom[self.z]
        self.nshell = copy.deepcopy(atom.nshell)
        self.shell = copy.deepcopy(atom.shell)
        del atom, bs
        # Update object with user defined value. Every line defines a shell
        for idx, line in enumerate(info):
            data = line.split()
            if len(data) == 2:
                continue
            if self.shell[idx].angular != int(data[1]) or self.shell[idx].ngtf != int(data[2]):
                warnings.warn(
                    'The sequence of input BS is not consistent with BSE data. Trying to find the matched shell...')
                found_flag = False
                for j in range(len(info)):
                    data2 = info[j].split()
                    if self.shell[j].angular == int(data2[1]) or self.shell[j].ngtf != int(data2[2]):
                        self.shell[j].chg = float(data2[3])
                        self.shell[j].scale = float(data2[4])
                        found_flag = True
                        break
                if found_flag == False:
                    raise Exception(
                        'The input BS is not a correct STO-3G / 3(6)-21G basis set.')
            else:
                self.shell[idx].chg = float(data[3])
                self.shell[idx].scale = float(data[4])
        return

    @property
    def data(self):
        """
        Print basis set information in CRYSTAL format
        """
        return self._print_data()

    def _print_data(self):
        """
        Return to a string in CRYSTAL format

        Returns:
            bs_str (str): Formatted basis set string of an atom.
        """
        bs_str = ''
        bs_str += format('%-5i%-5i\n' % (self.zconv, self.nshell))
        if self.ecp != None:
            bs_str += format('%s\n' % self.ecp)

        for s in self.shell:
            bs_str += s.data

        return bs_str


class BasisSetBASE():
    """
    The basisset object in CRYSTAL format. When called, one can pass the
    following arguments:

    Args:
        name (str | None): When specified, download the corresponding basis set
            from `Basis Set Exchange(BSE) <https://www.basissetexchange.org/>`_
        element (list[str] | list[int]): Elements to download. Either as a list
            of atomic numbers or labels.

    .. note::

        For basis sets from BSE, reference and roles are omitted.

    """

    def __init__(self):
        pass

    @classmethod
    def from_bse(cls, name, element, zconv=None):
        """
        Download basis set definitions from BSE.

        Args:
            name (str): Basis set's name.
            element (list[str] | list[int]): List of elements.
        """
        import basis_set_exchange as bse

        bs_str = bse.get_basis(name, elements=element, fmt='crystal')

        return cls()._set_atom(bs_str, zconv)

    @classmethod
    def from_string(cls, bs_str, fmt='crystal'):
        """
        Define basis set from a string.

        Args:
            bs_str (str)
            fmt (str): Format string. Consistent with BSE python API.
        """
        import re
        import basis_set_exchange as bse

        if not re.match(r'^crystal$', fmt, re.IGNORECASE):
            bs_str = bse.convert_formatted_basis_str(bs_str, fmt, 'crystal')

        return cls()._set_atom(bs_str)

    @classmethod
    def from_file(cls, file, fmt='crystal'):
        """
        Define a basis set from a file.
        """
        import re
        import basis_set_exchange as bse

        bs_file = open(file, 'r')
        bs_str = bs_file.read()
        bs_file.close()
        if not re.match(r'^crystal$', fmt, re.IGNORECASE):
            bs_str = bse.convert_formatted_basis_str(bs_str, fmt, 'crystal')

        return cls()._set_atom(bs_str)

    def _set_atom(self, info, zconv=None):
        """
        Assign basis set parameters at atom level.

        Args:
            info (str): String in Crystal format
            zconv (list[int]): If not none, use the conventional atomic number.
                Its length must be nelement. Its sequence must be consistent
                with basis set's
        """
        import re
        from CRYSTALpytools.base.basisset import BS_Atom

        info = info.strip().split('\n')
        info = [i.strip() for i in info]

        atom_list = []
        atom_range = []
        for idx, line in enumerate(info):
            if re.match(r'^[0-9]+\s+[0-9]+$', line):
                if re.match(r'^99\s+0$', line):
                    atom_range.append(idx)
                    break
                if zconv == None:
                    z = int(line.split()[0])  # Conventional atomic number
                else:
                    z = zconv[len(atom_range) - 1]
                atom_list.append(z)
                atom_range.append(idx)

        self.atom = {}
        for idx, z in enumerate(atom_list):
            atom_obj = BS_Atom(z, info[atom_range[idx]: atom_range[idx + 1]])
            self.atom[z] = atom_obj

        return self

    @property
    def data(self):
        """
        Print basis set information in CRYSTAL format
        """
        bs_str = ''
        for a in list(self.atom.values()):
            bs_str += a.data
        bs_str += '99   0\n'

        return bs_str

    def to_file(self, file='BASISSET.DAT', fmt='crystal'):
        """
        Print formatted data into a text file.

        Args:
            file (str)
            fmt (str): Output format
        """
        import os
        import re
        import warnings
        import basis_set_exchange as bse

        f = open(file, 'w+')
        if os.path.isfile(file):
            warnings.warn('File exists. New entry will be attached to the end.')

        bs_str = self.data
        if not re.match(r'^crystal$', fmt, re.IGNORECASE):
            bs_str = bse.convert_formatted_basis_str(bs_str, 'crystal', fmt)

        f.write('%s' % bs_str)
        f.close()
