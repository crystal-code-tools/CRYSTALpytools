"""
Class and methods to deal with geometries, including Pymatgen Structrue and
Molecule geometries and CRYSTAL gui geometries
"""
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Molecule

import numpy as np

# very slow to call mendeleev instances. Define a periodic table here for some simple cases
ptable = {
    0 : 'X',
    1 : 'H', 2 : 'He',
    3 : 'Li', 4 : 'Be', 5 : 'B', 6 : 'C', 7 : 'N', 8 : 'O', 9 : 'F', 10 : 'Ne',
    11 : 'Na', 12 : 'Mg', 13 : 'Al', 14 : 'Si', 15 : 'P', 16 : 'S', 17 : 'Cl', 18 : 'Ar',
    19 : 'K', 20 : 'Ca',
    21 : 'Sc', 22 : 'Ti', 23 : 'V', 24 : 'Cr', 25 : 'Mn', 26 : 'Fe', 27 : 'Co', 28 : 'Ni', 29 : 'Cu', 30 : 'Zn',
    31 : 'Ga', 32 : 'Ge', 33 : 'As', 34 : 'Se', 35 : 'Br', 36 : 'Kr',
    37 : 'Rb', 38 : 'Sr',
    39 : 'Y', 40 : 'Zr', 41 : 'Nb', 42 : 'Mo', 43 : 'Tc', 44 : 'Ru', 45 : 'Rh', 46 : 'Pd', 47 : 'Ag', 48 : 'Cd',
    49 : 'In', 50 : 'Sn', 51 : 'Sb', 52 : 'Te', 53 : 'I', 54 : 'Xe',
    55 : 'Cs', 56 : 'Ba',
    57 : 'La', 58 : 'Ce', 59 : 'Pr', 60 : 'Nd', 61 : 'Pm', 62 : 'Sm', 63 : 'Eu', 64 : 'Gd', 65 : 'Tb', 66 : 'Dy', 67 : 'Ho', 68 : 'Er', 69 : 'Tm', 70 : 'Yb', 71 : 'Lu',
    72 : 'Hf', 73 : 'Ta', 74 : 'W', 75 : 'Re', 76 : 'Os', 77 : 'Ir', 78 : 'Pt', 79 : 'Au', 80 : 'Hg',
    81 : 'Tl', 82 : 'Pb', 83 : 'Bi', 84 : 'Po', 85 : 'At', 86 : 'Rn',
    87 : 'Fr', 88 : 'Ra',
    89 : 'Ac', 90 : 'Th', 91 : 'Pa', 92 : 'U', 93 : 'Np', 94 : 'Pu', 95 : 'Am', 96 : 'Cm', 97 : 'Bk', 98 : 'Cf', 99 : 'Es'
}

ptable_inv = dict(zip(ptable.values(), ptable.keys()))


class Crystal_gui():
    """
    This class can read a CRYSTAL gui file into an object or substrate
    information of the object to generate a gui file.

    Args:
        dimensionality (int): Number of dimensions
        lattice (array): 3\*3 lattice matrix in Angstrom
        symmops (array): n_symmops\*4\*3 matrices of symmetry operators
        atom_number (array): natom\*1 int array of atomic numbers
        atom_positions (array): natom\*3 array of Cartesian coordinates
        space_group (int): CRYSTAL space group number
    """
    def __init__(self, dimensionality=None, lattice=None, symmops=None,
                 atom_number=None, atom_positions=None, space_group=None):
        self.dimensionality = dimensionality
        self.lattice = lattice
        self.symmops = symmops
        if np.all(symmops!=None):
            self.n_symmops = len(symmops)
        else:
            self.n_symmops = 0
        self.atom_number = atom_number
        self.atom_positions = atom_positions
        if np.all(atom_number!=None):
            self.n_atoms = len(atom_number)
        else:
            self.n_atoms = 0
        self.space_group = space_group

    def read_pmg(self, struc, pbc=None, vacuum=500., symmetry=True, zconv=None, **kwargs):
        """
        Read a pymatgen Structure object into a ``CRYSTAL_gui`` object. Vacuum
        layer is set to 500 Angstrom as the default of CRYSTAL for low symmstry
        systems.

        Args:
            struc (Structure|Molecule): Pymatgen Structure / Molecule object.
            pbc (list): 1\*3 boolian list. Implements periodicity along x, y and z
                directions. If none, the code will read it from input structure.
            vacuum (float): Vacuum distance. Unit: Angstrom. Low dimensional
                systems only.
            symmetry (bool): Do symmetry analysis.
            zconv (list[list[int, int]]): 1st element: The **index** of atom;
                    2nd element: The new conventional atomic number.
            **kwargs: Passed to Pymatgen SpacegroupAnalyzer object. Valid only
                if ``symmetry=True``.
        """
        from CRYSTALpytools.geometry import CStructure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.core.surface import center_slab
        from pymatgen.core.lattice import Lattice
        import numpy as np
        import warnings
        import copy

        if type(struc) == Structure:
            struc = CStructure.from_pmg(struc)
        elif type(struc) == Molecule or type(struc) == CMolecule:
            latt = Lattice(np.eye(3)*500, pbc=[False, False, False])
            struc = CStructure(lattice=latt, species=struc.species,
                               coords=struc.cart_coords, coords_are_cartesian=True)
        # dimensionality
        if np.all(pbc==None):
            pbc = struc.pbc

        self.dimensionality = pbc.count(True)

        # Vacuum distance
        latt_mx = np.eye(3)
        thickness_x = np.amax(struc.cart_coords[:, 0]) - np.amin(struc.cart_coords[:, 0])
        thickness_y = np.amax(struc.cart_coords[:, 1]) - np.amin(struc.cart_coords[:, 1])
        thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])
        latt_mx[0, 0] = thickness_x + vacuum
        latt_mx[1, 1] = thickness_y + vacuum
        latt_mx[2, 2] = thickness_z + vacuum

        if self.dimensionality == 0: # 0D
            self.lattice = latt_mx
            self.n_atoms = struc.num_sites
            self.space_group = 1
            self.symmops = []
            self.n_symmops = 1
            self.symmops = np.vstack([np.eye(3), [0.0,0.0,0.0]])
            self.atom_number = list(struc.atomic_numbers)
            self.atom_positions = struc.cart_coords
        else: # 1-3D, rotation and add vacuum layer
            if self.dimensionality == 2:
                if pbc[0] == False: # A no periodicity
                    struc = struc.rot_cel(struc.lattice.matrix[0, :], [0, 0, 1])

                elif pbc[1] == False: # B no periodicity
                    struc = struc.rot_cel(struc.lattice.matrix[1, :], [0, 0, 1])

                else: # C no periodicity
                    struc = struc.rot_cel(struc.lattice.matrix[2, :], [0, 0, 1])

            elif gui.dimensionality == 1:
                if pbc[0] == True: # A periodic
                    struc = struc.rot_cel(struc.lattice.matrix[0, :], [1, 0, 0])

                elif pbc[1] == True: # B periodic
                    struc = struc.rot_cel(struc.lattice.matrix[1, :], [1, 0, 0])

                else: # C periodic
                    struc = struc.rot_cel(struc.lattice.matrix[2, :], [1, 0, 0])

            struc = Structure(lattice=struc.lattice.matrix, species=struc.species,
                              coords=struc.cart_coords, coords_are_cartesian=True)

            self.lattice = struc.lattice.matrix
            self.n_atoms = struc.num_sites

            if symmetry == True:
                if self.dimensionality == 3:
                    self.space_group, self.n_symmops, self.symmops = struc.get_sg_symmops(**kwargs)
                elif self.dimensionality == 2:
                    # Get group number before editing- inheriated from previous version
                    self.space_group = SpacegroupAnalyzer(struc, **kwargs).get_space_group_number()
                    #center the slab first
                    struc = center_slab(struc)
                    # Then center at z=0.0
                    translation = np.array([0.0, 0.0, -0.5])
                    struc.translate_sites(list(range(structure.num_sites)),
                                          translation, to_unit_cell=False)
                    _, self.n_symmops, self.symmops = struc.get_sg_symmops(**kwargs)
                else:
                    warnings.warn('Polymer symmetry currently not examined.')
            else:
                self.space_group = 1
                self.n_symmops = 1
                self.symmops = np.vstack([np.eye(3), [0.0,0.0,0.0]])
                self.symmops = np.reshape(np.array(gui.symmops, dtype=float),
                                         [gui.n_symmops, 4, 3])

            self.atom_number = list(struc.atomic_numbers)
            self.atom_positions = struc.cart_coords

        if np.all(zconv!=None):
            for atom in zconv:
                self.atom_number[atom[0]] = atom[1]

        return self

    def read_gui(self, gui_file):
        """
        Read CRYSTAL gui file and genreate a ``Crystal_gui`` object.

        Args:
            gui_file (str): The CRYSTAL structure (gui) file
        """
        import numpy as np

        file = open(gui_file, 'r')
        data = file.readlines()
        file.close()

        self.dimensionality = int(data[0].split()[0])
        self.lattice = []
        self.symmops = []
        for i in range(1, 4):
            self.lattice.append([float(x) for x in data[i].split()])
        self.n_symmops = int(data[4].split()[0])
        for i in range(5, 5+self.n_symmops*4):
            self.symmops.append(data[i].split())
        self.symmops = np.reshape(np.array(self.symmops, dtype=float),
                                  [self.n_symmops, 4, 3])
        self.n_atoms = int(data[5+self.n_symmops*4].split()[0])
        self.atom_number = []
        self.atom_positions = []
        for i in range(6+self.n_symmops*4, 6+self.n_symmops*4+self.n_atoms):
            atom_line = data[i].split()
            self.atom_number.append(int(atom_line[0]))
            self.atom_positions.append([float(x) for x in atom_line[1:]])
        self.space_group = int(data[-1].split()[0])

        return self

    def write_gui(self, gui_file, symm=True, pseudo_atoms=[]):
        """
        Write a CRYSTAL gui file (to file)

        Args:
            gui_file (str): The name of the gui that is going to be written (
                including .gui).
            symm (bool): Whether to include symmetry operations.
            pseudo_atoms (list[int]): A list of atoms whose core is described
                by a pseudopotential (conventional atomic number = atomic
                number + 200)
        """
        import numpy as np

        if symm == False:
            self.n_symmops = 1
            self.symmops = np.vstack([np.eye(3), [0.0, 0.0, 0.0]])

        file = open(gui_file, 'w')
        # First line
        file.writelines('%4s   1   1\n' % self.dimensionality)
        # Cell vectors
        for vector in self.lattice:
            file.writelines('{}\n'.format(
                ''.join(['{0: 20.12E}'.format(np.round(n, 12)) for n in vector])
            ))
        # N symm ops
        file.writelines('{:5d}\n'.format(self.n_symmops))

        # symm ops
        sym_list = np.reshape(self.symmops, [self.n_symmops*4, 3])
        for symmops in sym_list:
            file.writelines('{}\n'.format(
                ''.join(['{0: 20.12f}'.format(np.round(n, 12)) for n in symmops])
            ))
        # N atoms
        file.writelines('{:5d}\n'.format(self.n_atoms))

        # atom number (including pseudopotentials) + coordinates cart
        for i in range(self.n_atoms):
            if self.atom_number[i] in pseudo_atoms:
                file.writelines('{:5d}{}\n'.format(
                    int(self.atom_number[i])+200,
                    ''.join(['{0: 20.12E}'.format(np.round(x, 12)) for x in self.atom_positions[i]])
                ))
            else:
                file.writelines('{:5d}{}\n'.format(
                    int(self.atom_number[i]),
                    ''.join(['{0: 20.12E}'.format(np.round(x, 12)) for x in self.atom_positions[i]])
                ))

        # space group + n symm ops
        if symm == True:
            file.writelines('{:5d}{:5d}\n'.format(
                self.space_group, self.n_symmops
            ))
        else:
            file.writelines('{:5d}{:5d}\n'.format(1, 1))

        file.close()


class CStructure(Structure):
    """
    Inherited from `Pymatgen Structure <https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.Structure>`_
    object with added methods.

    Arguments not listed are the same as `Pymatgen Structure <https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.Structure>`_

    Args:
        species (list[int]): Same as pymatgen or a 1\*nAtom list of
            **conventional** atomic numbers.
        symmetry_group (int): Symmetry group number or symbol in CRYSTAL
            convention.
        pbc (list|tuple): Periodicity.
        standarize (bool): Whether to use the CRYSTAL convention of periodic
            boundaries for low dimensional materials. It calls the
            ``standarize_pbc()`` method.
        \*\*kwargs: Other arguments passed to pymatgen Structure.
    Returns:
        self (CStructure): ``CStructure`` object.
    """
    def __init__(self, lattice, species, coords, symmetry_group=1, pbc=None,
                 standarize=False, **kwargs):
        import numpy as np
        from pymatgen.core.lattice import Lattice

        # conventional atomic number
        if isinstance(species[0], int) or isinstance(species[0], float) \
        or isinstance(species[0], np.int64) or isinstance(species[0], np.float64):
            zconv = [int(i) for i in species]
            species = [ptable[int(i % 100)] for i in species]
        else:
            zconv = []
        # PBC
        if isinstance(lattice, Lattice):
            if np.all(pbc==None): pbc = lattice.pbc
            latt = Lattice(lattice.matrix, pbc=pbc)
        else:
            if np.all(pbc==None): pbc = (True, True, True)
            latt = Lattice(lattice, pbc=pbc)
        # structure
        kwargs['lattice'] = latt
        kwargs['species'] = species
        kwargs['coords'] = coords
        super().__init__(**kwargs)
        # standarization
        ## Orientation
        if standarize == True: self.standarize_pbc()
        # Extras
        self._symmetry_group = symmetry_group
        if zconv == []:
            self._species_Z = [i.Z for i in self.species]
        else:
            self._species_Z = zconv

    @classmethod
    def from_pmg(cls, struc):
        """
        Get a ``CStructure`` object from Pymatgen structure. ``symmetry_group``
        currently not available.
        """
        if not isinstance(struc, Structure):
            raise ValueError('Not a Pymatgen Structure object')
        obj = cls(lattice=struc.lattice, species=struc.species,
                  coords=struc.cart_coords, coords_are_cartesian=True)
        return obj.standarize_pbc()

    @property
    def ndimen(self):
        """
        Dimensionality (1~3) of the system.
        """
        return self.pbc.count(True)

    @property
    def species_symbol(self):
        """
        Atom symbols.
        """
        return [ptable[int(i % 100)] for i in self.species_Z]

    @property
    def species_Z(self):
        """
        Conventional atomic numbers
        """
        return self._species_Z

    @property
    def crys_coords(self):
        """
        Composite fractional / Cartesian atomic coordinates. Consistent
        with CRYSTAL conventions.

        * 3D: Fractional
        * 2D: Frac, Frac, Cart
        * 1D: Frac, Cart, Cart
        * 0D: Cart, Cart, Cart
        """
        import numpy as np

        if self.ndimen == 0:
            crys_coords = self.cart_coords
        elif self.ndimen == 1:
            crys_coords = np.zeros([natoms, 3], dtype=float)
            crys_coords[:, 0] = self.frac_coords[:, 0]
            crys_coords[:, 1:] = self.cart_coords[:, 1:]
        elif self.ndimen == 2:
            crys_coords = np.zeros([natoms, 3], dtype=float)
            crys_coords[:, :2] = self.frac_coords[:, :2]
            crys_coords[:, 2] = self.cart_coords[:, 2]
        else:
            crys_coords = self.frac_coords
        return crys_coords

    def standarize_pbc(self):
        """
        Use the CRYSTAL standard periodic boundary for low dimensional materials.
        """
        import numpy as np
        from pymatgen.core.lattice import Lattice

        latt = self.lattice.matrix
        if self.ndimen == 3:
            return self
        elif self.ndimen == 0:
            cart_coords = self.cart_coords
            self = CStructure(Lattice(latt, pbc=(False, False, False)),
                              species=self.species, coords=cart_coords,
                              coords_are_cartesian=True)
            return self

        def rotate(v1, v2): # A quick rotation function. v1: old, v2: new
            from scipy.spatial.transform import Rotation
            import numpy as np

            v1 = v1 / np.linalg.norm(v1); v2 = v2 / np.linalg.norm(v2)
            vn = np.cross(v1,v2)
            if np.linalg.norm(vn) < 1e-4:
                return Rotation.from_rotvec([0., 0., 0.])
            else:
                vn = vn / np.linalg.norm(vn)
                ang = np.arccos(np.dot(v1,v2) / np.linalg.norm(v1) / np.linalg.norm(v2))
                return Rotation.from_rotvec(vn*ang)

        if self.ndimen == 2:
            if self.pbc[0] == False:
                oldv = latt[0, :]
            elif self.pbc[1] == False:
                oldv = latt[1, :]
            else: # in case of tilted z
                oldv = latt[2, :]

            cart_coords = self.cart_coords
            rot = rotate(oldv, [0, 0, 1])
            latt = rot.apply(latt)
            latt[0, 2] = 0
            latt[1, 2] = 0
            latt[2] = [0, 0, 500]
            cart_coords = rot.apply(cart_coords)
            self = CStructure(Lattice(latt, pbc=(True, True, False)),
                              species=self.species, coords=cart_coords,
                              standarize=False, coords_are_cartesian=True)
        elif self.ndimen == 1:
            if self.pbc[0] == True: # in case of tilted x
                oldv = latt[0, :]
            elif self.pbc[1] == True:
                oldv = latt[1, :]
            else:
                oldv = latt[2, :]

            cart_coords = self.cart_coords
            latt = self.lattice.matrix
            rot = rotate(oldv, [1, 0, 0])
            latt = rot.apply(latt)
            latt[0, 1:] = 0
            latt[1] = [0, 500, 0]
            latt[2] = [0, 0, 500]
            cart_coords = rot.apply(cart_coords)
            self = CStructure(Lattice(latt, pbc=(True, False, False)),
                              species=self.species, coords=cart_coords,
                              standarize=False, coords_are_cartesian=True)
        return self

    def refine_geometry(self, **kwargs):
        """
        Get refined geometry. Useful when reducing the cell to the irrducible
        one. 3D only.

        Args:
            **kwargs: Passed to Pymatgen `SpacegroupAnalyzer <https://pymatgen.org/pymatgen.symmetry.html#pymatgen.symmetry.analyzer.SpacegroupAnalyzer>`_ object.
        Returns:
            self (CStructure): New attributes listed below
            sg (int): Space group number
            pstruc (Structure): Irrducible structure that is consistent with
                International Crystallographic Table
            platt (list): minimal set of crystallographic cell parameters
            natom_irr (int): number of irrducible atoms
            atom (list): natom\*4 array. 1st element: atomic number; 2-4:
                fractional coordinates
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        import numpy as np

        ndimen = self.pbc.count(True)
        if ndimen < 3:
            raise Exception('This method is for 3D systems only.')

        analyzer = SpacegroupAnalyzer(self, **kwargs)
        # Analyze the refined geometry
        struc1 = analyzer.get_refined_structure()
        analyzer2 = SpacegroupAnalyzer(struc1, **kwargs)
        struc2 = analyzer2.get_primitive_standard_structure()
        analyzer3 = SpacegroupAnalyzer(struc2, **kwargs)

        struc3 = analyzer3.get_symmetrized_structure()
        struc4 = analyzer3.get_refined_structure()
        sg = analyzer3.get_space_group_number()

        latt = []
        if sg >= 1 and sg < 3:  # trilinic
            for i in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
                latt.append(getattr(struc4.lattice, i))
        elif sg >= 3 and sg < 16:  # monoclinic
            for i in ['a', 'b', 'c', 'beta']:
                latt.append(getattr(struc4.lattice, i))
        elif sg >= 16 and sg < 75:  # orthorhombic
            for i in ['a', 'b', 'c']:
                latt.append(getattr(struc4.lattice, i))
        elif sg >= 75 and sg < 143:  # tetragonal
            for i in ['a', 'c']:
                latt.append(getattr(struc4.lattice, i))
        elif sg >= 143 and sg < 168:  # trigonal, converted to hexagonal
            struc5 = analyzer3.get_conventional_standard_structure()
            analyzer4 = SpacegroupAnalyzer(struc5, **kwargs)
            struc3 = analyzer4.get_symmetrized_structure()
            struc4 = analyzer4.get_refined_structure()
            for i in ['a', 'c']:
                latt.append(getattr(struc4.lattice, i))
        elif sg >= 168 and sg < 195:  # hexagonal
            for i in ['a', 'c']:
                latt.append(getattr(struc4.lattice, i))
        else:  # cubic
            latt.append(struc4.lattice.a)

        self.sg = sg
        self.pstruc = struc4
        self.platt = latt
        self.natom_irr = len(struc3.equivalent_sites)
        self.atom = []
        invlatt = np.linalg.inv(struc4.lattice.matrix)
        for i in np.unique(np.array(struc3.site_labels, dtype=int)):
            frac_in_struc4 = struc3.cart_coords[i] @ invlatt
            self.atom.append([struc3.species[i].Z,
                              frac_in_struc4[0],
                              frac_in_struc4[1],
                              frac_in_struc4[2]])
        return self

    def get_sg_symmops(self, **kwargs):
        """
        Get space group number and corresponding symmetry operations. To keep
        consistency with International Crystallographic Table, refined geometry
        is suggested.

        Args:
            **kwargs: Passed to Pymatgen SpacegroupAnalyzer object.
        Returns:
            self (CStructure): New attributes are listed below
            sg (int): Space group number
            n_symmops (int): number of symmetry operations
            symmops (array): n_symmops\*4\*3 array of symmetry operations
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        import numpy as np

        struc = SpacegroupAnalyzer(self, **kwargs).get_refined_structure()
        self.sg = SpacegroupAnalyzer(struc, **kwargs).get_space_group_number()
        all_symmops = SpacegroupAnalyzer(struc, **kwargs).get_symmetry_operations(cartesian=True)
        self.symmops = []
        ops_tmp = []
        self.n_symmops = 0
        # For symmetry operations with same rotation matrix, save the one with 0
        # tranlation vector.
        for symmop in all_symmops:
            if self.n_symmops == 0:
                self.n_symmops += 1
                self.symmops.append(np.vstack([symmop.rotation_matrix, symmop.translation_vector]))
                ops_tmp = [symmop]
            else:
                save = None
                for nop, op in enumerate(ops_tmp):
                    if np.array_equal(op.rotation_matrix, symmop.rotation_matrix):
                        if np.all(op.translation_vector == 0.):
                            save = False
                            break
                        else:
                            save = True
                            save_id = nop
                            break
                    else:
                        continue

                if save == True: # Same rotation, choose the one with no translation
                    self.symmops[save_id] = np.vstack([symmop.rotation_matrix, symmop.translation_vector])
                    ops_tmp[save_id] = symmop
                elif np.all(save==None): # New rotation
                    self.symmops.append(np.vstack([symmop.rotation_matrix, symmop.translation_vector]))
                    ops_tmp.append(symmop)
                    self.n_symmops += 1
                else:
                    continue

        self.symmops = np.reshape(np.array(self.symmops, dtype=float), [self.n_symmops, 4, 3])

        return self

    def get_pcel(self, smx):
        """
        Restore the supercell to primitive cell, with the origin shifted to the
        middle of lattice to utilize symmetry (as the default of CRYSTAL).

        Args:
            smx (array): 3\*3 array of *supercell expansion matrix*. Inverse
                will be taken automatically.
        Returns:
            pcel (CStructure): Pymatgen structure of primitive cell with
                CRYSTALpytools methods.
        """
        from pymatgen.core.lattice import Lattice
        import numpy as np

        ndimen = self.pbc.count(True)
        natom = self.num_sites
        pbc = self.pbc

        # That forces origin back to (0.5,0.5,0.5), but makes pbc to be 3D
        super().make_supercell([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        shrink_mx = np.linalg.inv(smx)
        scel_mx = self.lattice.matrix
        # Shift origin to (0,0,0), consistent with CRYSTAL
        all_coords = self.cart_coords
        for i in range(natom):
            for j in range(ndimen):
                all_coords[i, 0:ndimen] -= 0.5 * scel_mx[j, 0:ndimen]

        pcel_mx = shrink_mx @ scel_mx
        pcel_latt = Lattice(pcel_mx, pbc=pbc)
        # Fractional coords of pcel: Both periodic and no periodic sites
        all_coords = all_coords @ np.linalg.inv(pcel_mx)
        pcel_coords = []
        pcel_species = []
        for i, coord in enumerate(all_coords.round(12)): # Slightly reduce the accuracy
            if np.any(coord[0:ndimen] >= 0.5) or np.any(coord[0:ndimen] < -0.5):
                continue
            else:
                pcel_coords.append(coord)
                pcel_species.append(self.species[i])

        # For low dimen systems, this restores the non-periodic vecter length
        pcel = CStructure(lattice=pcel_latt, species=pcel_species,
                          coords=pcel_coords, coords_are_cartesian=False)

        return pcel

    def get_scel(self, smx):
        """
        Get the supercell from primitive cell, with the origin shifted to the
        middle of lattice to utilize symmetry (as the default of CRYSTAL).

        Args:
            smx (array): 3\*3 array of supercell expansion matrix
        Returns:
            scel (CStructure): Pymatgen structure of supercell
        """
        from pymatgen.core.lattice import Lattice

        ndimen = self.pbc.count(True)
        pbc = self.pbc
        natom = self.num_sites

        super().make_supercell(smx)
        scel_mx = self.lattice.matrix
        # Shift origin to (0,0,0), consistent with CRYSTAL
        all_coords = self.cart_coords
        for i in range(natom):
            for j in range(ndimen):
                all_coords[i, 0:ndimen] -= 0.5 * scel_mx[j, 0:ndimen]

        scel_latt = Lattice(struc.lattice.matrix, pbc=pbc)

        scel = CStructure(lattice=scel_latt, species=self.species,
                          coords=all_coords, coords_are_cartesian=True)
        return scel

    def rot_cel(self, vec1, vec2):
        """
        Rotate the geometry according to 2 vectors. A `rotation vector <https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector>`_
        is defined.

        Args:
            vec1 (array): A Cartesian vector before rotation
            vec2 (array): A Cartesian vector after rotation

        Returns:
            rcel (CStructure): Pymatgen structure of rotated cell
        """
        import numpy as np
        from scipy.spatial.transform import Rotation as Rot
        from pymatgen.core.lattice import Lattice

        vec1 = vec1 / np.linalg.norm(vec1)
        vec2 = vec2 / np.linalg.norm(vec2)
        # define a rotation
        rotvec = np.cross(vec1, vec2)
        if np.all(np.abs(rotvec) < 1e-4): # vec1 and 2 along the same direction
            rotvec = np.zeros([3,])
        else:
            rotvec = rotvec / np.linalg.norm(rotvec) * np.arccos(np.dot(vec1, vec2))
        rot = Rot.from_rotvec(rotvec)

        # lattice
        latt_mx = rot.apply(self.lattice.matrix)
        latt = Lattice(latt_mx, pbc=self.pbc)

        # coordinates
        coords = rot.apply(self.cart_coords)

        rcel = CStructure(lattice=latt, species=self.species,
                          coords=coords, coords_are_cartesian=True)

        return rcel

    def miller_norm(self, miller):
        """
        Find the norm vector of a specified Miller plane

        Args:
            miller (array | list): Miller indices. 3\*1 1D array or nmiller\*3 3D array

        Returns:
            vec (array): Norm vector, normalized to 1. 3\*1 1D array or nmiller\*3 3D array
        """
        import numpy as np

        if len(np.shape(miller)) == 1:
            dimen = 1
            miller = np.array([miller], dtype=int)
        elif len(np.shape(miller)) > 1:
            dimen = 2
            miller = np.array(miller, dtype=int)
        else:
            raise ValueError('Unknown Miller index format.')

        if np.shape(miller)[1] != 3:
            raise ValueError('Unknown Miller index format.')
        nmiller = np.shape(miller)[0]

        vec = np.zeros([nmiller, 3])
        for im, m in enumerate(miller):
            cart_vec = self.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(m)
            vec[im, :] = cart_vec / np.linalg.norm(cart_vec)

        if dimen == 1:
            vec = vec[0, :]

        return vec

    def write_gui(self, gui_file=None, pbc=None, vacuum=500., symmetry=True,
                  zconv=None, **kwargs):
        """
        Read a pymatgen Structure object into a ``CRYSTAL_gui`` object. Vacuum
        layer is set to 500 Angstrom as the default of CRYSTAL for low symmstry
        systems.

        *Developing*

        Args:
            struc (Structure|Molecule): Pymatgen Structure / Molecule object.
            pbc (list): 1\*3 boolian list. Implements periodicity along x, y and z
                directions. If none, the code will read it from input structure.
            vacuum (float): Vacuum distance. Unit: Angstrom. Low dimensional
                systems only.
            symmetry (bool): Do symmetry analysis.
            zconv (list[list[int, int]]): 1st element: The **index** of atom;
                    2nd element: The new conventional atomic number.
            **kwargs: Passed to Pymatgen SpacegroupAnalyzer object. Valid only
                if ``symmetry=True``.
        """

    # ---------------------------- Connectivity ------------------------------#

    def get_bonds(self, scale=1.2, special_bonds={}):
        """
        .. _ref-CStrucGetBonds:

        Get bond connectivity based on distances between atoms.

        Args:
            scale (float): Scale the sum of atomic radius A and radius B.
            special_bonds(dict): Dictionary of bond lengths that are not
                compliant to ``scale``. Should be defined as ``{'A:B' : len}``,
                with ``A`` ``B`` being element symbols and ``len`` being bond
                length in :math:`\\AA`. The sequence of ``A`` ``B`` is arbitrary.
        Returns:
            self (CStructure): New attributes listed below.
            self.bonds (list): nBond\*3 list. The first and second elements are
                atom indices (starting from 1) of the bond. The third element
                is a 1\*3 integer array lattice vector of the second element.
                The first atom always has the lattice vector of \[0,0,0\].
            self.bond_matrix (bsr_array): nAtom\*nAtom sparse matrix of bond
                connectivity. In `scipy.sparse.bsr_array <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.bsr_array.html>`_
                format.
        """
        import numpy as np
        import warnings
        from scipy.sparse import bsr_array

        # get MAX bond lengths
        uspec = np.unique(self.species)
        bondMX = {}
        for i in range(len(uspec)):
            elei = uspec[i]
            # between same atoms
            bondMX['{}:{}'.format(elei.symbol, elei.symbol)] = elei.data['Atomic radius'] * 2 * scale
            for j in range(i+1, len(uspec)):
                elej = uspec[j]
                bondMX['{}:{}'.format(elei.symbol, elej.symbol)] = (elei.data['Atomic radius'] + elej.data['Atomic radius']) * scale
                bondMX['{}:{}'.format(elej.symbol, elei.symbol)] = bondMX['{}:{}'.format(elei.symbol, elej.symbol)]

        if len(special_bonds) != 0:
            for k in special_bonds.keys():
                ka = k.split(':')[0].capitalize()
                kb = k.split(':')[1].capitalize()
                key1 = '{}:{}'.format(ka, kb)
                if key not in bondMX.keys():
                    warnings.warn("The specified pair '{}' does not exist in the structure. The definition is ignored.".format(k),
                                  stacklevel=2)
                    continue

                bondMX[key1] = special_bonds[k]
                if ka != kb:
                    key2 = '{}:{}'.format(kb, ka)
                    bondMX[key2] = special_bonds[k]


        # neighboring cells
        self.standarize_pbc()
        if self.ndimen == 3:
            nbr_cell = np.array([
                [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1],
                [0, 1, 1], [1, 1, 1], [-1,0, 0], [0,-1, 0], [0, 0,-1], [1,-1, 0],
                [-1,1, 0], [1, 0,-1], [-1,0, 1], [0, 1,-1], [0,-1, 1], [-1,1, 1],
                [1,-1, 1], [1, 1,-1], [-1,-1,0], [-1,0,-1], [0,-1,-1], [1,-1,-1],
                [-1,1,-1], [-1,-1,1], [-1,-1,-1]], dtype=int)
        elif self.ndimen == 2:
            nbr_cell = np.array([
                [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [-1,0, 0], [0,-1, 0],
                [1,-1, 0], [-1,1, 0], [-1, -1, 0]], dtype=int)
        elif self.ndimen == 1:
            nbr_cell = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], dtype=int)
        else:
            nbr_cell = np.array([[0, 0, 0]], dtype=int)
        nbr_cell_c = nbr_cell @ self.lattice.matrix

        self.bonds = []
        self.bond_matrix = np.zeros([self.num_sites, self.num_sites])
        # make the code faster
        sbol = [i.symbol for i in self.species]
        crds = self.cart_coords
        nrepeat = nbr_cell.shape[0]
        nsite = self.num_sites
        for i in range(nsite):
            v1 = crds[i]
            s1 = sbol[i]
            for j in range(i+1,nsite):
                v2 = np.repeat([crds[j]], nrepeat, axis=0) + nbr_cell_c
                s2 = sbol[j]
                dists = np.linalg.norm(v1 - v2, axis=1)
                bondTOL = bondMX['{}:{}'.format(s1, s2)]
                for k in np.where(dists<bondTOL)[0]: # only loop 0 or 1 times
                    self.bonds.append([i+1, j+1, nbr_cell[k]])
                    self.bond_matrix[i, j] = 1
                    self.bond_matrix[j, i] = 1
                    break

        self.bond_matrix = bsr_array(self.bond_matrix)
        return self

    def get_molecules(self, **kwargs):
        """
        Substract molecule from molecular crystals or clusters based on bond
        connectivity. For ionic or covalent structures where bonds form an
        integrated network, 1 molecule of all atoms are returned.

        Args:
            \*\*kwargs: See the ``get_bonds()`` method.
        Returns:
            self (CStructure): New attributes listed below.
            self.molecules (list): nMolecule\*nAtom_per_molecule int list of
                atom indices (from 1).
        """
        from scipy.sparse.csgraph import connected_components

        if not hasattr(self, 'bond_matrix'): self.get_bonds(**kwargs)

        nmolecule, label = connected_components(self.bond_matrix)
        self.molecules = [[] for i in range(nmolecule)]
        for iat, imole in enumerate(label):
            self.molecules[imole].append(iat+1)
        return self

    # ---------------------------- Visualization -----------------------------#
    def visualize(self,
                  atom_color='cpk',
                  atom_data=[],
                  atom_null=(0., 0., 0.),
                  atom_cbar_label='Atom Data',
                  bond_color=(0.5, 0.5, 0.5),
                  bond_data=[],
                  bond_null=(0., 0., 0.),
                  bond_cbar_label='Bond Data',
                  atom_bond_ratio='medium',
                  cell_display=True,
                  cell_color=(0., 0., 0.),
                  cell_linewidth=0.1,
                  display_matrix=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                  display_origin=[0., 0., 0.],
                  display_range=[[0., 1.], [0., 1.], [0., 1.]],
                  show_the_scene=True,
                  **kwargs):
        """
        Visualize 3D atomic structure with `MayaVi <https://docs.enthought.com/mayavi/mayavi/>`_
        (*not installed by default*).

        2 display modes are available:

        #. Normal display. Distinguishing atoms by elements.  
        #. Data display. Coloring atoms/bonds/both by the data assigned.

        The displayed stucture is manipulated by the following equation:

        .. math::

            \\mathbf{L} = \\mathbf{M}\\mathbf{L_{0}} + \\mathbf{p}\\mathbf{M}\\mathbf{L_{0}}

        where :math:`\\mathbf{M}` is defined by ``display_matrix``, :math:`\\mathbf{p}`
        is defined by ``display_origin``. :math:`\\mathbf{L_{0}}` is the input
        structure. The object itself will not be changed.

        **A known issue**

        Due to settings of MayaVi, bonds plotted in bond data display mode are
        plotted as lines rather than tubes. That might make bonds look very
        thick in the default rendering window, which can be addressed by
        zooming in or expand the rendering window.

        Args:
            atom_color (str): Color map of atoms. 'jmol' or 'cpk' for normal
                display (by elements). Or MayaVi colormaps for data display.
            atom_data (array): nAtom\*2 array. The first column is atom indices
                (from 1) and the second is data. Atoms without assigned data
                are plotted in the color specified by ``atom_null``.
            atom_null (turple): *Useful only for data-display mode*. Color of
                atoms without data assigned.
            atom_cbar_label (str): Title of atom data colorbar.
            bond_color (turple|str): Color of bonds for normal display or
                MayaVi colormaps for data display.
            bond_data (array): nBond\*3 array. The first 2 columns are atom
                indices (from 1) connecting atom A and B. The 3rd column is
                data.
            bond_null (turple): *Useful only for data-display mode*. Color of
                bonds without data assigned.
            bond_cbar_label (str): Title of bond data colorbar.
            atom_bond_ratio (str): 'balls', 'large', 'medium', 'small' or
                'sticks'. The relative sizes of balls and sticks.
            cell_display (bool): Display lattice boundaries (at \[0., 0., 0.\] only).
            cell_color (turple): Color of lattice boundaries.
            cell_linewidth (float): Linewidth of plotted lattice boundaries.
            display_matrix (array): nDimen\*nDimen **integer** matrix to
                manipulate the structure. Matrices with determinant :math:`\\geq`
                1 are meaningful, otherwise error will be displayed.
            display_origin (array): 1\*3 array of structure origin in fractional
                coordinates of the **expanded supercell**.
            display_range (array): 3\*2 array defining the displayed region.
                Fractional coordinates a, b, c are used but only the periodic
                directions are applied. Defined according to the **expanded
                and translated supercell**.
            show_the_scene (bool): Display the scene by ``mlab.show()`` and
                return None. Otherwise return the scene object.
            \*\*kwargs: Other keywords passed to MayaVi or ``self.get_bonds()``.
                See below.
            scale (float): Scale the sum of atomic radius A and radius B.
            special_bonds(dict): Dictionary of bond lengths that are not
                compliant to ``scale``. Should be defined as ``{'A:B' : len}``,
                with ``A`` ``B`` being element symbols and ``len`` being bond
                length in :math:`\\AA`. The sequence of ``A`` ``B`` is arbitrary.
            azimuth: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            elevation: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            distance: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
                By default set to 'auto'.
            focalpoint: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            roll: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
        Returns:
            None
        """
        try:
            from mayavi import mlab
        except ModuleNotFoundError:
            raise ModuleNotFoundError("MayaVi is required for this method.")
        from ase.data import colors
        import numpy as np
        import copy, warnings

        # Colormap by original indices
        if atom_color.lower() in ['jmol', 'cpk']:
            at_data = False
            if len(atom_data)>0:
                warnings.warn("'atom_data' available but colormap not specified. Using the default display mode.",
                              stacklevel=2)
        else:
            at_data = True
            if len(atom_data)==0:
                raise Exception("'atom_data' not specified or incorrect 'atom_color' input.")

        if not isinstance(bond_color, str):
            bd_data = False
            if len(bond_data)>0:
                warnings.warn("'bond_data' available but colormap not specified. Using the default display mode.",
                              stacklevel=2)
        else:
            bd_data = True
            if len(bond_data)==0:
                raise Exception("'bond_data' not specified or incorrect 'bond_color' input.")

        # Ball and Stick ratio
        if atom_bond_ratio.lower() not in ['balls', 'large', 'medium', 'small', 'sticks']:
            raise ValueError("Invalid atom_bond_ratio: {}".format(atom_bond_ratio))
        scale_factors = {'balls' : 2.0,
                         'large' : 1.0,
                         'medium': 0.8,
                         'small' : 0.6,
                         'sticks': 0.0}
        atscale = scale_factors[atom_bond_ratio.lower()]

        # Bond connectivity
        bondkwd = ['scale', 'special_bonds']; newkwd = dict()
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in bondkwd: newkwd[k] = v

        self.get_bonds(**newkwd)
        # Supercell
        ndimen = self.pbc.count(True)
        smx = np.eye(3, dtype=int)
        smx[0:ndimen, 0:ndimen] = np.array(display_matrix, dtype=int, ndmin=2)[0:ndimen, 0:ndimen]
        if np.round(abs(np.linalg.det(smx)), 12) < 1:
            raise ValueError('Determinant of input display matrix must >= 1.')
        struc = self.make_supercell(smx, in_place=False) # this will get a 3D structure

        # atoms to plot, map supercell indices back, when data display mode is on
        if bd_data == True or at_data == True:
            atplt = np.zeros([struc.num_sites, 4])
            pfrac = np.round(struc.cart_coords @ np.linalg.inv(self.lattice.matrix) % 1, 6)
            reffrac = np.round(self.frac_coords % 1, 6)
            for iat, pf in enumerate(pfrac):
                iref = np.where((pf==reffrac).all(axis=1))[0] # must have 1 entry
                atplt[iat, 0] = iref
                atplt[iat, 1:] = struc.frac_coords[iat]
        else:
            atplt = np.round(np.hstack([
                [[i] for i in range(struc.num_sites)],
                struc.frac_coords
            ]), 12)

        # Display range and origin
        origin = np.array(display_origin, ndmin=1, dtype=float)
        dispbg = np.round([i[0] for i in display_range], 12) + origin
        disped = np.round([i[1] for i in display_range], 12) + origin

        idx = np.where(disped-dispbg<1e-12)[0]
        if len(idx) > 0:
            direct = ['a', 'b', 'c'][idx[0]]
            raise Exception("Display range error along {} axis!\n{} min = {:.2f}, {} max = {:.2f}. No structure is displayed.".format(
                direct, direct, dispbg[idx[0]], direct, disped[idx[0]]))

        ## Atoms in display range
        for i, nbg, ned in zip([1,2,3], np.floor(dispbg), np.ceil(disped)):
            tmp = copy.deepcopy(atplt)
            for s in np.arange(nbg, ned, 1):
                if s == 0: continue
                ttmp = copy.deepcopy(tmp)
                ttmp[:,i] = ttmp[:,i] + s
                atplt = np.vstack([atplt, ttmp])
                del ttmp
            del tmp
        # reduce to fractional boundary
        atplt = atplt[np.where((atplt[:,1]>=dispbg[0])&(atplt[:,1]<disped[0]))]
        atplt = atplt[np.where((atplt[:,2]>=dispbg[1])&(atplt[:,2]<disped[1]))]
        atplt = atplt[np.where((atplt[:,3]>=dispbg[2])&(atplt[:,3]<disped[2]))]
        atplt[:, 1:] = np.round(atplt[:, 1:] @ struc.lattice.matrix, 12)
        if len(atplt) == 0:
            raise Exception("No atom exist in the defined boundary. Check your input.")
        ##  Data associated to atoms
        if at_data == True:
            saved_atoms = []
            atom_data = np.array(atom_data, dtype=float, ndmin=2)
            atom_data[:,0] -= 1
            dat = np.empty([atplt.shape[0],]); dat[:] = np.nan
            for i in range(atplt.shape[0]):
                idx = np.where(atom_data[:,0]==atplt[i, 0])[0]
                if idx.shape[0] == 0: continue
                elif idx.shape[0] == 1: dat[i] = atom_data[idx, 1]
                else: raise Exception("Multiple values defined for atom {:.0f}.".format(
                    atom_data[idx[0], 0]+1))
                saved_atoms.append(idx[0])
            # not assigned data
            saved_atoms = np.unique(saved_atoms)
            if len(saved_atoms) < atom_data.shape[0]:
                warnings.warn("Not all atoms defined by 'atom_data' are available in the visualized structure.",
                              stacklevel=2)
            del saved_atoms, idx, atom_data

        # Bonds to plot, new connectivity based on atplt is needed.
        latt = np.multiply(struc.lattice.matrix, (disped-dispbg).reshape([3,1]))
        species = [struc.species[int(i[0])].symbol for i in atplt]
        obj = CStructure(latt, species, atplt[:,1:4], pbc=self.pbc, coords_are_cartesian=True)
        obj.get_bonds(**newkwd)
        del species, latt

        nbond = len(obj.bonds)
        idx1 = np.array([i[0] for i in obj.bonds], dtype=int) - 1
        idx2 = np.array([i[1] for i in obj.bonds], dtype=int) - 1
        lattpt = np.array([i[2] for i in obj.bonds], dtype=float)
        # align bonds to + directions
        idx_neg = np.unique(np.where(lattpt<0)[0])
        lattpt[idx_neg] = np.abs(lattpt[idx_neg])
        tmp = copy.deepcopy(idx1[idx_neg])
        idx1[idx_neg] = copy.deepcopy(idx2[idx_neg])
        idx2[idx_neg] = copy.deepcopy(tmp)
        del idx_neg, tmp

        if bd_data == False:
            bdplt = np.hstack([ # bond idx, cart coord A, cart coord B
                np.linspace(0, nbond-1, nbond).reshape([-1, 1]),
                atplt[idx1, 1:],
                atplt[idx2, 1:] + lattpt@obj.lattice.matrix,
            ])
            del idx1, idx2, lattpt, obj
        else:
            # Explicit matching of bonds
            bdplt = []
            bdold = np.array([[i[0], i[1]] for i in self.bonds], dtype=int) - 1

            bond_data = np.array(bond_data, dtype=float, ndmin=2)
            bond_data[:,0] -= 1; bond_data[:,1] -= 1
            dbd = np.empty([nbond,]); dbd[:] = np.nan # data to plot

            saved_bonds = []
            for i in range(nbond):
                at1 = atplt[idx1[i]]; at2 = atplt[idx2[i]]
                oldidx = np.sort([at1[0], at2[0]])
                ibd = np.where((bdold==oldidx).all(axis=1))[0]
                if ibd.shape[0] == 0: continue
                bdplt.append(
                    np.hstack([ibd[0], at1[1:], at2[1:]])
                )

                idx = np.where((bond_data[:,0:2]==oldidx).all(axis=1))[0]
                if idx.shape[0] == 0: continue
                elif idx.shape[0] == 1: dbd[i] = bond_data[idx[0], 2]
                else: raise Exception(
                    "Multiple values defined for the bond between atom {:.0f} and {:.0f}.".format(
                    self.bonds[ibd[0]][0], self.bonds[ibd[0]][1])
                )
                saved_bonds.append(idx[0])

            bdplt = np.array(bdplt)
            bdplt[:, 4:7] += lattpt @ obj.lattice.matrix
            # not assigned data
            saved_bonds = np.unique(saved_bonds)
            if len(saved_bonds) < bond_data.shape[0]:
                warnings.warn("Not all bonds defined by 'bond_data' are available in the visualized structure.",
                              stacklevel=2)
            del idx1, idx2, lattpt, obj, bdold, bond_data, oldidx, saved_bonds

        # plot
        ## fig passed, developer only
        if 'fig' in kwargs.keys():
            del kwargs[k]
        else:
            fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
        ## Atoms
        radii = np.array([struc.species[int(i)].data['Atomic radius'] for i in atplt[:,0]])
        radii[np.where(radii<0.4)[0]] = 0.4 # H too small

        if at_data == False:
            cat = np.zeros(atplt.shape, dtype=np.uint8) + 255
            if atom_color.lower() == 'jmol':
                cat[:, 0:3] = np.array([
                    colors.jmol_colors[struc.species[int(i)].Z]*255 for i in atplt[:,0]
                ], dtype=np.uint8)
            else:
                cat[:, 0:3] = np.array([
                    colors.cpk_colors[struc.species[int(i)].Z]*255 for i in atplt[:,0]
                ], dtype=np.uint8)
            pts = mlab.pipeline.scalar_scatter(atplt[:,1], atplt[:,2], atplt[:,3],
                                               figure=fig)
            pts.add_attribute(cat, 'colors')
            pts.data.point_data.set_active_scalars('colors')
            pts.data.point_data.vectors = np.tile(radii, (3,1)).T
            g = mlab.pipeline.glyph(pts)
            g.glyph.glyph.scale_factor = atscale
            g.glyph.scale_mode = 'scale_by_vector'
        else:
            pts = mlab.points3d(atplt[:,1], atplt[:,2], atplt[:,3], dat,
                                figure=fig,
                                colormap=atom_color,
                                scale_factor=atscale)
            pts.glyph.scale_mode = 'scale_by_vector'
            pts.mlab_source.dataset.point_data.vectors = np.tile(radii, (3,1)).T
            cat = np.zeros([4,], dtype=np.uint8) + 255
            cat[0:3] = np.array([i*255 for i in atom_null], dtype=np.uint8)
            pts.module_manager.scalar_lut_manager.lut.nan_color = tuple(cat)
        # Bonds
        if atom_bond_ratio.lower() == 'sticks':
            bond_radius = 0.2
        else:
            bond_radius = 0.15

        if bd_data == False:
            for ib, b in enumerate(bdplt):
                bds = mlab.plot3d([b[1], b[4]], [b[2], b[5]], [b[3], b[6]],
                                  figure=fig,
                                  color=bond_color,
                                  tube_radius=bond_radius)
        else: # use lines rather than tubes
            bond_radius = 10
            src = mlab.pipeline.scalar_scatter(bdplt[:, [1,4]].flatten(),
                                               bdplt[:, [2,5]].flatten(),
                                               bdplt[:, [3,6]].flatten(),
                                               np.vstack([dbd, dbd]).T.flatten())
            src.mlab_source.dataset.lines = [[i, i+1] for i in range(0, 2*nbond, 2)]
            src.update()
            lines = mlab.pipeline.stripper(src); del src
            bds = mlab.pipeline.surface(lines,
                                        figure=fig,
                                        colormap=bond_color,
                                        line_width=bond_radius)
            cbd = np.zeros([4,], dtype=np.uint8) + 255
            cbd[0:3] = np.array([i*255 for i in bond_null], dtype=np.uint8)
            bds.module_manager.scalar_lut_manager.lut.nan_color = tuple(cbd)

        ## Lattice
        if cell_display == True:
            if ndimen == 3:
                lattpts = np.array([[0, 0, 0], [1, 0, 0],
                                    [1, 0, 0], [1, 1, 0],
                                    [1, 1, 0], [0, 1, 0],
                                    [0, 1, 0], [0, 0, 0],
                                    [0, 0, 1], [1, 0, 1],
                                    [1, 0, 1], [1, 1, 1],
                                    [1, 1, 1], [0, 1, 1],
                                    [0, 1, 1], [0, 0, 1],
                                    [0, 0, 0], [0, 0, 1],
                                    [1, 0, 0], [1, 0, 1],
                                    [1, 1, 0], [1, 1, 1],
                                    [0, 1, 0], [0, 1, 1]]) @ struc.lattice.matrix
                lattpts = np.add(lattpts, origin@struc.lattice.matrix)
                for i in range(0, 24, 2):
                    mlab.plot3d([lattpts[i,0], lattpts[i+1,0]],
                                [lattpts[i,1], lattpts[i+1,1]],
                                [lattpts[i,2], lattpts[i+1,2]],
                                figure=fig,
                                color=cell_color,
                                line_width=cell_linewidth)
            elif ndimen == 2:
                lattpts = np.zeros([8, 3], dtype=float)
                lattpts[:, 0:2] = np.array([[0, 0], [1, 0],
                                            [1, 0], [1, 1],
                                            [1, 1], [0, 1],
                                            [0, 1], [0, 0]]) @ struc.lattice.matrix[0:2, 0:2]
                lattpts = np.add(lattpts, origin@struc.lattice.matrix)
                for i in range(0, 8, 2):
                    mlab.plot3d([lattpts[i,0], lattpts[i+1,0]],
                                [lattpts[i,1], lattpts[i+1,1]],
                                [lattpts[i,2], lattpts[i+1,2]],
                                figure=fig,
                                color=cell_color,
                                line_width=cell_linewidth)
            elif ndimen == 1:
                lattpts = np.zeros([2, 3], dtype=float)
                lattpts[:, 0:1] = np.array([[0], [1]]) @ struc.lattice.matrix[0:1, 0:1]
                lattpts = np.add(lattpts, origin@struc.lattice.matrix)
                mlab.plot3d([lattpts[0,0], lattpts[1,0]],
                            [lattpts[0,1], lattpts[1,1]],
                            [lattpts[0,2], lattpts[1,2]],
                            figure=fig,
                            color=cell_color,
                            line_width=cell_linewidth)

        ## Add color bar
        if at_data == True:
            mlab.colorbar(object=pts, orientation='horizontal', title=atom_cbar_label, label_fmt='%.2f')
        if bd_data == True:
            mlab.colorbar(object=bds, orientation='vertical', title=bond_cbar_label, label_fmt='%.2f')

        ## Final setups
        keys = ['azimuth', 'elevation', 'distance', 'focalpoint', 'roll']
        keywords = dict(figure=fig, distance='auto')
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v
        mlab.view(**keywords)

        if show_the_scene == False:
            return fig
        else:
            mlab.gcf().scene.parallel_projection = True
            mlab.show()
            return


class CMolecule(Molecule):
    """Deprecated, use CStructure"""
    def __init__(species, coords, symmetry_group=1, **kwargs):
        import warnings
        import numpy as np

        if isinstance(species[0], int) or isinstance(species[0], float) \
        or isinstance(species[0], np.int64) or isinstance(species[0], np.float64):
            zconv = [int(i) for i in species]
            species = [ptable[int(i % 100)] for i in zconv]

        self = Molecule(species=species, coords=coords, **kwargs)
        self._symmetry_group = symmetry_group
        self._species_Z = zconv

        warnings.warn("This is a deprecated method and returns to a pymatgen Molecule class only. Use CStructure with pbc=(False, False, False) instead.",
                      stacklevel=2)










