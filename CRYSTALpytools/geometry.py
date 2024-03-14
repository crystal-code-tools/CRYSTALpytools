"""
Class and methods to deal with geometries, including Pymatgen Structrue and
Molecule geometries and CRYSTAL gui geometries
"""
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Molecule


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
        if symmops != None:
            self.n_symmops = len(symmops)
        else:
            self.n_symmops = 0
        self.atom_number = atom_number
        self.atom_positions = atom_positions
        if atom_number != None:
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
        if pbc == None:
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

        if zconv != None:
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
    object with added methods. Instantiation method is the same as Pymatgen Structure.
    """
    @classmethod
    def from_pmg(cls, struc):
        """
        Get a ``CStructure`` object from Pymatgen structure.
        """
        if type(struc) != Structure:
            raise ValueError('Not a Pymatgen Structure object')

        return cls(lattice=struc.lattice, species=struc.species,
                   coords=struc.cart_coords, coords_are_cartesian=True)

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
        natom_eq = int(struc3.num_sites / self.natom_irr)
        for i in range(self.natom_irr):
            idx_eq = int(i * natom_eq)
            z = struc4.species[idx_eq].Z
            self.atom.append([z, struc4.sites[idx_eq].a, struc4.sites[idx_eq].b,
                             struc4.sites[idx_eq].c])

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
                elif save == None: # New rotation
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

    def Miller_norm(self, miller, d=1.0):
        """
        Find the norm vector of a specified Miller plane

        Args:
            miller (array | list): 3\*1 list of Miller index
            d (fload): Length of norm vector

        Returns:
            vec (array): 3\*1 norm vector, normalized to 1.
        """
        import numpy as np

        zeros = np.argwhere(abs(miller) < 1)
        if len(zeros) == 0:
            vec1 = np.array([1/miller[0], 0, 0]) - np.array([0, 0, 1/miller[2]])
            vec2 = np.array([0, 1/miller[1], 0]) - np.array([0, 0, 1/miller[2]])
            vec1 = np.dot(vec1, self.lattice.matrix)
            vec2 = np.dot(vec2, self.lattice.matrix)
            vec = np.cross(vec1, vec2)
        elif len(zeros) == 1:
            if zeros[0][0] == 0:
                vec1 = [1, 0, 0]
                vec2 = np.array([0, 1/miller[1], 0]) - np.array([0, 0, 1/miller[2]])
            elif zeros[0][0] == 1:
                vec1 = np.array([1/miller[0], 0, 0]) - np.array([0, 0, 1/miller[2]])
                vec2 = [0, 1, 0]
            else:
                vec1 = [0, 0, 1]
                vec2 = np.array([1/miller[0], 0, 0]) - np.array([0, 1/miller[1], 0])
            vec1 = np.dot(vec1, self.lattice.matrix)
            vec2 = np.dot(vec2, self.lattice.matrix)
            vec = np.cross(vec1, vec2)
        else:
            if zeros[0][0] == 0 and zeros[1][0] == 1:
                vec = np.array([0., 0., 1.])
            elif zeros[0][0] == 0 and zeros[1][0] == 2:
                vec = np.array([0., 1., 0.])
            else:
                vec = np.array([1., 0., 0.])
            vec = np.dot(vec, self.lattice.matrix)

        vec = vec / np.linalg.norm(vec) * d
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


class CMolecule(Molecule):
    """
    Inherited from `Pymatgen Molecule <https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.Molecule>`_
    object with added methods. Instantiation method is the same as Pymatgen Molecule.

    Developing
    """










