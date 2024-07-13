"""
Basis set definitions.
"""

class AtomBS():
    """
    Basis set class for CRYSTALpytools

    Args:
        z (int): Conventional atomic number
        nshell (int): Number of shells
        ECP (str): Effective core potential, printed before basis sets
    """
    def __init__(self, z, nshell, ECP=None):
        from mendeleev import element

        self.z = int(z)
        self.nshell = int(nshell)
        self.ECP = ECP
        self.element = element(z % 100)
        self.shells = []
        # self.charge = self.element.atomic_number # For charge neutrality check - all-electron BSs only

    def define_a_shell(self, ITYB, LAT, NG, CHE, SCAL, orbitals=[]):
        """
        Define shells of atomic basis sets. Arguments are consistent with
        CRYSTAL manual.

        Args:
            ITYB (int): Type of Gaussain orbitals
            LAT (int): Angular momentum
            NG (int): Number of Gaussian orbitals
            CHE (float): Charge
            SCAL (float): Scaling factor
            orbitals (list): nGTO\*2 or nGTO\*3 (sp type) list of exponent
                and contraction factors
        Returns:
            obj (AtomBS)
        """
        shell = {}
        if int(ITYB) not in [0, 1, 2]:
            raise ValueError('Unknown basis set type.')
        shell['type'] = int(ITYB)

        if int(LAT) not in range(6): # Maximum g shell
            raise ValueError('Unknown angular momentum')
        shell['angular momentum'] = int(LAT)

        if int(NG) != len(orbitals) and ITYB == 0:
            raise ValueError('Inconsistent definition of number of Gaussian orbitals.')
        shell['n orbital'] = int(NG)
        shell['charge'] = float(CHE)
        shell['scale factor'] = float(SCAL)
        shell['orbitals'] = orbitals

        self.shells.append(shell)
        # self.charge -= CHE
        return self

    @classmethod
    def read_bse(cls, bs, z):
        """
        Read basis sets saved in `Basis Set Exchange(BSE) <https://www.basissetexchange.org/>`_ 
        dictionary format and assign charge.

        .. note::
            All electron basis sets only. Atoms are charge neutral.

        Args:
            bs (str): Name of basis set
            z (int): Conventional atomic number
        Returns:
            obj (AtomBS)
        """
        from mendeleev import element
        import numpy as np
        import warnings, copy
        try:
            import basis_set_exchange as bse
        except ImportError:
            raise ImportError('basis_set_exchange package is not available or incompatible with your hardware (a known issue for arm64 chips).')

        zreal = int(z)%100
        dic = bse.get_basis(bs, elements=[element(zreal).symbol])

        try:
            nshell = len(dic['elements'][str(zreal)]['electron_shells'])
        except KeyError:
            raise ValueError("Specified basis set '{}' does not have GTO definitions. Is it a pesudopotential name?".format(bs))
        # some basis set definitions do not use 'electron_shells' entry but repeats coefficients
        nshell = 0
        dicold = copy.deepcopy(dic)
        dic['elements'][str(zreal)]['electron_shells'] = []
        for ishell, shell in enumerate(dicold['elements'][str(zreal)]['electron_shells']):
            # Convert them into Non-repeated coefficients, 1 entry for 1 shell
            nshell_ang = len(shell['coefficients'])
            ## Do not covert sp orbitals or non-repeated definitions
            if len(shell['angular_momentum']) == nshell_ang: 
                nshell += 1
                dic['elements'][str(zreal)]['electron_shells'].append(
                        {
                            'function_type'    : shell['function_type'],
                            'region'           : shell['region'],
                            'angular_momentum' : shell['angular_momentum'],
                            'exponents'        : shell['exponents'],
                            'coefficients'     : shell['coefficients']
                        }
                    )
            ## Convert
            else:
                nshell += nshell_ang
                for jshell in range(nshell_ang):
                    coeff = np.array(shell['coefficients'][jshell], dtype=float)
                    idx_real = np.where(coeff!=0.)[0]
                    dic['elements'][str(zreal)]['electron_shells'].append(
                        {
                            'function_type'    : shell['function_type'],
                            'region'           : shell['region'],
                            'angular_momentum' : shell['angular_momentum'],
                            'exponents'        : [shell['exponents'][int(k)] for k in idx_real],
                            'coefficients'     : [[str(k) for k in coeff[idx_real]]]
                        }
                    )
        obj = cls(z, nshell)
        chg = obj._assign_charge(obj.element, dic) # automatically assign charge

        # Add GTO shell
        for ishell in range(nshell):
            angshell = dic['elements'][str(zreal)]['electron_shells'][ishell]['angular_momentum']
            if len(angshell) == 1: # s,p,d,f,g angular momentum in CRYSTAL format
                angshell = angshell[0]
                if angshell == 0:
                    angcrys = angshell
                else:
                    angcrys = angshell + 1
            else: # sp angular momentum in CRYSTAL format
                angcrys = 1

            if angcrys > 5:
                raise ValueError('H orbital and beyond is not available.')

            ngto = len(dic['elements'][str(zreal)]['electron_shells'][ishell]['exponents'])
            orbs = np.vstack([
                np.array(dic['elements'][str(zreal)]['electron_shells'][ishell]['exponents'], dtype=float),
                np.array(dic['elements'][str(zreal)]['electron_shells'][ishell]['coefficients'], dtype=float)
            ]).transpose()
            obj.define_a_shell(0, angcrys, ngto, chg[ishell], 1.0, orbitals=orbs)

        return obj

    @staticmethod
    def _assign_charge(element, bs):
        """
        Assign charge to basis sets from BSE, where the charge info is missing.

        BSE sorts shells by angular momentum acending order, i.e., 1s, 2s, 3s, 2p, 3p, 3d...

        .. note::
            All electron basis sets only. Atoms are charge neutral.

        Args:
            element (element): `Mendeleev <https://mendeleev.readthedocs.io/en/stable/>`_
                Element object.
            bs (dict): BSE basis set dictionary
        Returns:
            chg (list): nshell\*1 list of charge assigned to every shell.
        """
        import numpy as np

        z = str(element.atomic_number)

        # Get orbitals defined in the BS
        econf = [i for i in element.ec.conf]
        # ECP basis sets
        if np.all(bs['elements'][z].get('ecp_potentials')!=None):
            orbitals = []
            corechg = bs['elements'][z]['ecp_electrons']
            totchg = element.atomic_number
            for i in econf[::-1]:
                if corechg == totchg: # all valence electron got
                    break
                elif corechg + element.ec.conf[i] > totchg: # filled d/f/g orbitals
                    continue
                corechg = corechg + element.ec.conf[i]
                orbitals.append(i)

            orbitals = [i for i in orbitals[::-1]]
        else: # all electron BS
            orbitals = econf

        # reorder orbitals by s, p, d...
        ang_symbol = {0 : 's', 1 : 'p', 2 : 'd', 3 : 'f', 4 : 'g'}
        re_orbitals = []
        for i in ['s', 'p', 'd', 'f', 'g']: # s to g orbitals
            for j in range(1, 10): # arbitrary
                if (j, i) in orbitals:
                    re_orbitals.append((j, i))

        chg = []
        nshell = len(bs['elements'][z]['electron_shells'])
        iorbital = 0 # counter for orbitals
        for ishell in range(nshell):
            # Assign charges of physical shells to GTO shells
            angshell = bs['elements'][z]['electron_shells'][ishell]['angular_momentum']
            if len(angshell) == 1: # s, p, d, f, g orbitals
                angshell = angshell[0]
                try:
                    if ang_symbol[angshell] == re_orbitals[iorbital][1]: # same symbol, fill electrons in
                        chg.append(float(element.ec.conf[re_orbitals[iorbital]]))
                        iorbital += 1
                    else: # split-valance orbitals
                        chg.append(0.)
                except IndexError: # polarization orbitals
                    chg.append(0.)
            else: # sp orbitals
                try:
                    n = re_orbitals[iorbital][0]
                    sporbs = []
                    for i in re_orbitals:
                        if i[0] != n:
                            continue
                        if i[1] == 'd': # remove d, f.. orbitals
                            break
                        sporbs.append(i)
                    chgsp = 0.
                    for i in sporbs:
                        chgsp += float(element.ec.conf[i])
                        re_orbitals.remove(i) # current index must corresponds to 's'. Not +1
                    chg.append(chgsp)
                except IndexError: # split-valance or polarization orbitals
                    chg.append(0.)
        return chg

    @classmethod
    def read_crystal(cls, text):
        """
        Analyze text in CRYSTAL format.
        """
        import re

        lines = text.strip().split('\n')
        if re.match(r'^[0-9]+\s+[0-9]+$', lines[0].strip()):
            z = int(lines[0].strip().split()[0])
            nshell = int(lines[0].strip().split()[1])
        else:
            raise ValueError('Input string is not a CRYSTAL formatted basis set. Atom and n shells are unknown.')

        # ECP
        if not re.match(r'^[0-9]+', lines[1]):
            if 'INPUT' not in lines[1]:
                ecp = lines[1].strip()
                nl = 2
            else: # manual definition of ECP
                ecp = lines[1].strip()
                for nl, l in enumerate(lines[2:]):
                    if len(l.strip().split()) != 5:
                        ecp = ecp + '\n' + l
                    else:
                        break
                nl += 1
        else:
            ecp = None
            nl = 1

        obj = cls(z, nshell, ecp)
        title = []
        gto = []
        while nl < len(lines):
            data = lines[nl].strip().split()
            if len(data) == 5:
                title.append([int(data[0]), int(data[1]), int(data[2]), float(data[3]), float(data[4])])
                gto_shell = []
                ngto = int(data[2])
                if int(data[0]) == 0:
                    for i in range(ngto):
                        data2 = lines[nl + i + 1].strip().split()
                        gto_shell.append([float(i.replace('D', 'E')) for i in data2])
                    nl += ngto + 1
                    gto.append(gto_shell)
                else:
                    nl += 1
            else:
                nl += 1

        for nshell, shell in enumerate(title):
            obj.define_a_shell(shell[0], shell[1], shell[2], shell[3], shell[4], gto[nshell])

        return obj

    def print_crystal(self):
        """
        Print basis set into CRYSTAL format
        """
        import numpy as np

        if self.nshell != len(self.shells):
            raise ValueError('Number of shells is not consistent with shells defined.')

        txt_out = '{:d} {:d}\n'.format(self.z, self.nshell)
        if np.all(self.ECP!=None):
            txt_out = txt_out + '{}\n'.format(self.ECP)
        for s in self.shells:
            txt_out = txt_out + '{:d} {:d} {:d} {:.2f} {:.2f}\n'.format(
                s['type'], s['angular momentum'], s['n orbital'], s['charge'], s['scale factor']
            )
            if s['type'] == 0:
                for o in s['orbitals']:
                    if len(o) == 3: # sp orbitals
                        txt_out = txt_out + '{: 18.10f}{: 18.10f}{: 18.10f}\n'.format(o[0], o[1], o[2])
                    elif len(o) == 2: # others
                        txt_out = txt_out + '{: 18.10f}{: 18.10f}\n'.format(o[0], o[1])
                    else:
                        raise ValueError('Input error: Unknown GTO definition.')

        return txt_out


class BasisSetBASE():
    """
    The basisset object base object, as the class of basis sets defined for the
    whole system. Basis sets from string, file or `Basis Set Exchange(BSE) <https://www.basissetexchange.org/>`_,
    can be read and saved as a list of ``AtomBS`` objects. 
    """

    def __init__(self, atoms=[]):
        self.atoms = atoms

    @classmethod
    def from_bse(cls, bs, z):
        """
        Get or append basis sets from `Basis Set Exchange(BSE) <https://www.basissetexchange.org/>`_.

        Args:
            bs (str): Name of basis set. Only one name is accepted.
            z (int | list[int]): Conventional atomic number.
        Returns:
            cls
        """
        # from CRYSTALpytools.base.basisset import AtomBS

        if type(z) == int:
            z = [z]
        elif type(z) != list and type(z) != tuple:
            raise ValueError('Conventional atomic number must be either list or int.')

        return cls([AtomBS.read_bse(bs, onez) for onez in z])

    @classmethod
    def from_string(cls, bs, fmt='crystal'):
        """
        Parse basis set strings.

        Args:
            bs (str): Basis set string.
            fmt (str): Format string. Consistent with BSE python API. For non-
                CRYSTAL formats, only all-electron basis sets are supported.
                Charge of each shell will be automatically assigned to get
                charge neutral atoms if ``fmt`` is not 'crystal'.
        Returns:
            cls
        """
        from mendeleev import element
        import re
        try:
            import basis_set_exchange as bse
        except ImportError:
            if fmt.lower() != 'crystal':
                raise ValueError("""'basis_set_exchange' python module is probably not installed or incompatible with your hardware (arm64 chips).
In this case, only 'fmt=crystal' is accepted.""")
        # from CRYSTALpytools.base.basisset import AtomBS

        atoms = []
        if fmt.lower() != 'crystal':
            bs = bse.read_formatted_basis_str(bs, fmt)
            # get definitions of inidival atoms
            atomlist = []
            elementlist = []
            for i in bs['elements'].keys():
                atom = {
                    'molssi_bse_schema' : bs['molssi_bse_schema'],
                    'elements'          : {i : bs['elements'][i]},
                    'function_types'    : bs['function_types'],
                    'name'              : bs['name'],
                    'description'       : bs['description'],
                }
                atomlist.append(atom)
                elementlist.append(element(int(i)))

            # Get basis sets
            for i in range(len(atomlist)):
                atom = AtomBS.read_crystal(bse.write_formatted_basis_str(atomlist[i], 'crystal'))
                # charge info might be missing
                chg = atom._assign_charge(elementlist[i], atomlist[i])
                for sh in range(len(atom.shells)):
                    atom.shells[sh]['charge'] = chg[sh]
                atoms.append(atom)

        else:
            bs = bs.strip().split('\n')
            block = ''
            endflag = False
            for line in bs:
                if re.match(r'^\s*[0-9]+\s+[0-9]+$', line):
                    if block != '':
                        atoms.append(AtomBS.read_crystal(block))
                    data = line.strip().split()
                    if data[0] == '99' and data[1] == '0':
                        block = ''
                        break
                    else:
                        block = line + '\n'

                elif re.match(r'^\s*$', line):
                    continue
                else:
                    block += line + '\n'
            # When '99 0' is not added
            if block != '':
                atoms.append(AtomBS.read_crystal(block))

        return cls(atoms)

    @classmethod
    def from_file(cls, bs, fmt='crystal'):
        """
        Parse basis set files.

        Args:
            bs (str): Basis set file.
            fmt (str): Format string. Consistent with BSE python API. For non-
                CRYSTAL formats, only all-electron basis sets are supported.
                Charge of each shell will be automatically assigned to get
                charge neutral atoms.
        Returns:
            cls
        """
        file = open(bs, 'r')
        bstr = file.read()
        file.close()
        return cls.from_string(bs=bstr, fmt=fmt)

    def print_crystal(self):
        """
        Print the information into CRYSTAL basis set format
        """
        bstr = ''
        for atom in self.atoms:
            bstr += atom.print_crystal()
        bstr += '99 0\n'

        return bstr
