#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Objects of input / output files of CRYSTAL. Methods to edit or substract data
from corresponding files are provided.
"""
from CRYSTALpytools import units
from CRYSTALpytools.base.crysd12 import Crystal_inputBASE
from CRYSTALpytools.base.propd3 import Properties_inputBASE
from CRYSTALpytools.base.output import POutBASE
from CRYSTALpytools.geometry import Crystal_gui

import numpy as np


class Crystal_input(Crystal_inputBASE):
    """
    Crystal input object inherited from the :ref:`Crystal_inputBASE <ref-base-crysd12>`
    object. For the basic set-ups of keywords, please refer to manuals :ref:`there <ref-base-crysd12>`.

    Args:
        source (str): Template file or string. An empty object is generated for
            ``''``.
    """
    def __init__(self, source=''):
        import os

        super().__init__()
        if type(source) != str:
            raise ValueError('Unknown input value. Only string is allowed')

        if source != '':
            if os.path.exists(source):
                inp = open(source, 'r')
                source = inp.read()
                inp.close()
            super().analyze_text(source)

    def read_file(self, source):
        """
        Initialize the object if not yet done.

        Args:
            source (str): The name of the input file.
        Returns:
            self (Crystal_input)
        """
        self.__init__(source)

    def write_file(self, file):
        """
        Write data to a file
        """
        out = open(file, 'w')
        out.write('%s' % self.data)
        out.close()
        return self

    def geom_from_cif(self, file, zconv=None, keyword='EXTERNAL',
                      pbc=[True, True, True], gui_name='fort.34', **kwargs):
        """
        Read geometry from cif file and put infomation to geom block, either as
        'EXTERNAL' or 'CRYSTAL'. CIF files with a single geometry only.

        CIF with Symmetry is required.

        .. note::

            CIF files should have the '.cif' extension.

        Args:
            file (str): CIF file.
            keyword (str): 'EXTERNAL' or 'CRYSTAL'.
            zconv (list[list[int, int]]): 1st element: The **index** of atom;
                2nd element: The new conventional atomic number. Atoms of the
                irreducible unit is required.
            pbc (list[bool]): Valid only if ``keyword = EXTERNAL``. Force to remove
                periodic boundary conditions along x, y, z axis.
            gui_name (str): Valid only if ``keyword = EXTERNAL``. Gui file's name.
            **kwargs: Passed to Pymatgen `SpacegroupAnalyzer <https://pymatgen.org/pymatgen.symmetry.html#pymatgen.symmetry.analyzer.SpacegroupAnalyzer>`_ object.
        """
        import re
        from CRYSTALpytools.geometry import CStructure
        from pymatgen.core.lattice import Lattice

        struc = CStructure.from_file(file)
        if keyword != 'CRYSTAL': # PBC might vary
            latt_mx = struc.lattice.matrix
            latt = Lattice(latt_mx, pbc=pbc)
            struc = CStructure(lattice=latt, species=struc.species,
                               coords=struc.cart_coords, coords_are_cartesian=True)

        self.geom_from_pmg(struc, zconv, keyword, gui_name, **kwargs)

        return self

    def geom_from_pmg(self, struc, zconv=None, keyword='EXTERNAL',
                      gui_name='fort.34', **kwargs):
        """
        Read geometry defined by PyMatGen structure object and put infomation
        into geom block, either as 'EXTERNAL' or 'CRYSTAL'.

        See ``geom_from_cif`` for definition of arguments.

        .. note::

            Coordinates of corresponding atoms may not consistent with the
            original CIF file if 'CRYSTAL' is used, in which case
            coordinates of another symmetry equivalent atom is used.
        """
        from CRYSTALpytools.convert import cry_pmg2gui
        from CRYSTALpytools.geometry import CStructure

        if type(struc) != CStructure:
            struc = CStructure.from_pmg(struc)

        if keyword.upper() == 'EXTERNAL':
            super(Crystal_input, self).geom.external()
            gui = cry_pmg2gui(struc, gui_file=gui_name, symmetry=True,
                              zconv=zconv, **kwargs)
        elif keyword.upper() == 'CRYSTAL':
            struc.refine_geometry(**kwargs)
            if np.all(zconv!=None):
                z_atom_index = [i[0] for i in zconv]
                for i in range(struc.natom_irr):
                    try:
                        atom_to_sub = z_atom_index.index(i)
                        z_input = zconv[atom_to_sub][1]
                    except ValueError:
                        z_input = struc.atom[i][0]
                    struc.atom[i][0] = z_input
            super().geom.crystal(IGR=struc.sg, latt=struc.platt, atom=struc.atom)
        else:
            raise ValueError("Input keyword format error: {}".format(keyword))

        return self

    def bs_user(self, bs_str, z=[], fmt='crystal', append=False, charge=None,
                ECP=None, BSfile=None, title='MYBASIS'):
        """
        A shortcut to get user defined (free format) basis sets from `Basis Set
        Exchange (BSE) <https://www.basissetexchange.org/>`_, formatted string
        or file.

        In some rare cases the user might want to manually change the charge
        assigned to the atoms to get ionic initial guess.

        Args:
            bs_str (str): Name of basis set(BSE), formatted string or file.
            z (list[int]): List of elements, specified by conventional atomic
                numbers. BSE only.
            fmt (str): Format string. Consistent with `BSE python API
                       <https://molssi-bse.github.io/basis_set_exchange/usage.html#versioning>`_.
                       For string and files.
            append (bool): Whether to cover old entries. If the old entry
                contains 'BASISSET', it will be removed anyway. Useful when
                different series of basis sets are defined
            charge (dict): To define charge. Keys: z, conventional atomic
                number; Values: charge of shells, whose sequence must be
                consistent with BS definition.
            ECP (dict): Definition of effective core pseudopotentials. Keys:
                z, conventional atomic number; Values: ECP string.
            BSfile (str): If not None, print basis set definitions into
                ``BSfile`` and use ``title`` as the entry to 'BASISSET' keyword
            title (str): Useful when ``BSfile`` is not None.
        """
        if z == []:
            try:
                file = open(bs_str, 'r') # file
                self.basisset.from_file(bs_str, fmt, append)
            except FileNotFoundError: # string
                self.basisset.from_string(bs_str, fmt, append)
        else: # BSE
            self.basisset.from_bse(bs_str, z, append)

        # set charge
        if np.all(charge!=None):
            for i in list(charge.keys()):
                found_atom = False
                for j, at in enumerate(self.basisset._bs_obj.atoms):
                    if at.z != i:
                        continue
                    else:
                        found_atom = True
                        break
                if found_atom == False:
                    raise Exception('Unknown z value.')
                if len(at.shells) != len(charge[i]):
                    raise Exception('Charge definition of element {} is inconsistent with its basis set.'.format(i))
                for ishell, chg in enumerate(charge[i]):
                    self.basisset._bs_obj.atoms[j].shells[ishell]['charge'] = chg

        # set ECP
        if np.all(ECP!=None):
            for i in list(ECP.keys()):
                found_atom = False
                for j, at in enumerate(self.basisset._bs_obj.atoms):
                    if at.z != i:
                        continue
                    else:
                        found_atom = True
                        break
                if found_atom == False:
                    raise Exception('Unknown z value.')
                self.basisset._bs_obj.atoms[j].ECP = ECP[i]

        if np.all(BSfile!=None):
            if append == False:
                file = open(BSfile, 'w')
            else:
                file = open(BSfile, 'a')
            file.write('%s\n' % title)
            file.write('%s' % self.basisset._bs_obj.print_crystal())
            file.close()
            self.basisset._bs_obj = None
            self.basisset._block_dict['BASISSET'][0] = None
            self.basisset.basisset(title)

        return self

    def bs_keyword(self, keyword):
        """
        A shortcut to define basis sets by keywords. Equivalent to ``self.basisset.basisset``

        Args:
            keyword (str): Name of basis set
        """
        self.basisset.basisset(keyword)
        return self


class Crystal_output:
    """
    This class reads a CRYSTAL output and generates an object.

    Args:
        output (str): Filename
    """

    def __init__(self, output=None):
        if np.all(output!=None):
            self.read_file(output)

    def read_file(self, output_name):
        """
        Reads a CRYSTAL output file.

        Args:
            output_name (str): Name of the output file.
        Returns:
            self (Crystal_output)
        """
        import re

        self.name = output_name

        # Check if the file exists
        try:
            if output_name[-3:] != 'out' and output_name[-4:] != 'outp':
                output_name = output_name+'.out'
            file = open(output_name, 'r', errors='ignore')
            self.data = file.readlines()
            file.close()
        except:
            raise FileNotFoundError('EXITING: a .out file needs to be specified')

        # Check the calculation terminated correctly
        self.terminated = False

        for i, line in enumerate(self.data[::-1]):
            if re.match(r'^ EEEEEEEEEE TERMINATION', line):
                self.terminated = True
                # This is the end of output
                self.eoo = len(self.data)-1-i
                break

        if self.terminated == False:
            self.eoo = len(self.data)

        return self

    def get_dielectric_tensor(self):
        """Extracts the dielectric tensor from the output.

        Returns:
            list: Dielectric tensor values.
        """
        import re

        for i, line in enumerate(self.data):
            if re.match(r'^ TENSOR IN PRINCIPAL AXES SYSTEM', line):
                # This is the end of output
                self.dielectric_tensor = [
                    float(x) for x in self.data[i+1].split()[1::2]]
                return self.dielectric_tensor
        return None

    def get_eigenvectors(self):
        """Extracts eigenvectors from the output."""
        import re

        for i, line in enumerate(self.data):
            if re.match(r'\s NUMBER OF AO', line) != None:
                self.num_ao = int(line.split()[3])

            if re.match(r'\s SHRINK. FACT.(MONKH.)', line) != None:
                self.num_k = int(line.split()[13])

            if re.match(r'\s SHRINK. FACT.(MONKH.)', line) != None:
                self.num_k = int(line.split()[13])

    #### Geometry ####

    def get_dimensionality(self):
        """
        Gets the dimensionality of the system.

        Returns:
            self.dimensionality (int): Dimensionality of the system.
        """
        import re

        self.dimensionality = None
        for line in self.data:
            if re.match(r'\sGEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM', line):
                self.dimensionality = int(line.split()[9])
                break

        if np.all(self.dimensionality==None):
            raise Exception('Invalid file. Dimension information not found.')

        return self.dimensionality

    def get_symmops(self):
        """
        Return the symmetry operators

        Returns:
            self.symmops (numpy.ndarray): Symmetry operators
        """
        import re
        import numpy as np

        self.n_symmops = 0
        self.symmops = np.array([])

        symmops = []
        for i, line in enumerate(self.data):
            if re.match(r'^ \*\*\*\*   \d+ SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS', line):
                self.n_symmops = int(line.split()[1])
                for j in range(0, self.n_symmops):
                    symmops.append(self.data[i+3+j].split()[2:])
                break

        if self.n_symmops > 0:
            self.symmops = np.reshape(np.array(symmops, dtype=float), [self.n_symmops, 4, 3])

        return self.symmops

    def get_geometry(self, initial=True, write_gui=False,
                     gui_name=None, symmetry='pymatgen', **kwargs):
        """
        Get the geometry.

        Args:
            initial (bool): Read the initial or last gemetry. Useful in
                case of geometry optimization.
            write_gui (bool): If True, write .gui file
            gui_name (str): Valid only if ``write_gui = True``. Gui file
                is named as 'gui_name'. If None, use 'basename.gui'. The
                basename is the same as output file.
            symmetry (str): Valid only if ``write_gui = True``. 'pymatgen'
                to use symmetry info from a pymatgen SpacegroupAnalyzer;
                'initial' to use symmstry information on output file. If
                None, no symmstry. Otherwise it is taken from the existing
                gui file.
            **kwargs: Valid only if ``write_gui = True`` and
                ``symmetry = 'pymatgen'``.  Passed to Pymatgen
                SpacegroupAnalyzer object.

        Returns:
            self.geometry (CStructure | CMolecule): A modified pymatgen Structure
                or molecule object.
        """
        import os
        import re
        import warnings
        import numpy as np
        from CRYSTALpytools.base.output import GeomBASE
        from CRYSTALpytools.convert import cry_pmg2gui
        from CRYSTALpytools.crystal_io import Crystal_gui

        # Get geometry
        bg_line = -1
        if initial == True:
            for nline, line in enumerate(self.data[:self.eoo]):
                # Use atom coords to read molecule geometries. Go 4 lines up for periodic systems
                if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                    bg_line = nline - 4
                    break
        else:
            for nline, line in enumerate(self.data[self.eoo::-1]):
                # Use atom coords to read molecule geometries. Go 4 lines up for periodic systems
                if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                    bg_line = len(self.data[:self.eoo]) - nline - 4
                    break

        if bg_line < 0:
            raise Exception('Geometry information not found.')

        output = GeomBASE.read_geom(self.data, bg_line)
        struc = output[1]

        # Get the last lattice matrix: structure obtained by GeomBASE might be rotated.
        ndimen = self.get_dimensionality()
        lattice_line = -1
        if ndimen != 0:
            if initial == True:
                for nline, line in enumerate(self.data[:self.eoo]):
                    if re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN', line):
                        lattice_line = nline
                        break
            else:
                for nline, line in enumerate(self.data[self.eoo::-1]):
                    if re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN', line):
                        lattice_line = len(self.data[:self.eoo]) - nline
                        break

        if lattice_line != -1:
            a_crys = np.array(self.data[lattice_line + 2].strip().split(), dtype=float)
            a_pmg = struc.lattice.matrix[0, :]
            # Rotate the geometry back
            struc = struc.rot_cel(a_pmg, a_crys)

        self.geometry = struc

        # Write gui files
        if write_gui == True:
            # Conventional atomic numbers
            zconv = [[i, self.atom_numbers[i]] for i in range(self.n_atoms)]
            if np.all(gui_name==None):
                gui_name = os.path.splitext(self.name)[0]
                gui_name = '{}.gui'.format(gui_name)

            if symmetry == 'pymatgen':
                gui = cry_pmg2gui(struc, gui_file=gui_name, symmetry=True, zconv=zconv, **kwargs)
            elif np.all(symmetry==None):
                gui = cry_pmg2gui(struc, gui_file=gui_name,  symmetry=False, zconv=zconv)
            elif symmetry == 'initial':
                self.get_symmops()
                gui = cry_pmg2gui(struc, gui_file=None, symmetry=False, zconv=zconv)
                gui.symmops = self.symmops
                gui.n_symmops = self.n_symmops
                gui.space_group = self.sg_number
                gui.write_gui(gui_name, symm=True)
            else:
                warnings.warn('Symmetry adapted from reference geometry. Make sure that is desired.',
                              stacklevel=2)
                gui_ref = Crystal_gui().read_gui(symmetry)
                gui = cry_pmg2gui(struc, gui_file=None, symmetry=False, zconv=zconv)
                # Replace the symmops with the reference file
                gui.symmops = gui_ref.symmops
                gui.n_symmops = gui_ref.n_symmops
                gui.space_group = gui_ref.space_group
                gui.write_gui(gui_list[idx_s], symm=True)

        return self.geometry

    def get_last_geom(self, write_gui_file=True, symm_info='pymatgen'):
        """
        Return the last optimised geometry.
        """
        struc = self.get_geometry(initial=False, write_gui=write_gui_file, symm_info=symm_info)
        if 'Molecule' in str(type(struc)):
            self.last_geom = [[[500., 0., 0.], [0., 500., 0.], [0., 0., 500.]],
                              self.atom_numbers,
                              self.atom_positions_cart.tolist()]
        else:
            self.last_geom = [struc.lattice.matrix.tolist(),
                              self.atom_numbers,
                              self.atom_positions_cart.tolist()]
        return self.last_geom

    def get_lattice(self, initial=True):
        """
        Returns the lattice of the system. Unit: Angstrom.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.lattice (np.ndarray): Lattice of the system.
        """
        import re
        import warnings
        import numpy as np

        ndimen = self.get_dimensionality()
        self.lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.lattice

        self.get_geometry(initial=initial, write_gui=False)
        self.lattice = self.geometry.lattice.matrix

        return self.lattice

    def get_reciprocal_lattice(self, initial=True):
        """
        Returns the reciprocal lattice of the system. Unit: Angstrom^-1.
        Reciprocal lattice is consistent with CRYSTAL: 2pi is added.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.reciprocal_lattice (np.ndarray): Lattice of the system.
        """
        import warnings

        ndimen = self.get_dimensionality()
        self.reciprocal_lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.reciprocal_lattice

        self.get_lattice(initial=initial)
        self.reciprocal_lattice = self.geometry.lattice.reciprocal_lattice.matrix

        return self.reciprocal_lattice

    @property
    def sg_number(self):
        """
        CRYSTAL 0~3D space group number. Before geometry editing.
        """
        import re
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        ndimen = self.get_dimensionality()
        sg = 'unknown'
        if ndimen == 0:
            rexp = r'^\s*POINT  GROUP  N\.'
        elif ndimen == 1:
            rexp = r'^POLYMER GROUP N\.'
        elif ndimen == 2:
            rexp = r'^\s*TWO\-SIDED PLANE GROUP N\.'

        if ndimen < 3:
            for nline, line in enumerate(self.data):
                if re.match(rexp, line):
                    if ndimen == 0:
                        sg = int(line.strip().split()[3])
                    elif ndimen == 1:
                        sg = int(line.strip()[16:19].strip())
                    elif ndimen == 2:
                        sg = int(line.strip().split()[3])
                    break
        else:
            struc = self.get_geometry(initial=True, write_gui=False)
            sg = SpacegroupAnalyzer(struc).get_space_group_number()
        return sg

    @property
    def sg_symbol(self):
        """
        CRYSTAL 0~3D 0~3D space group symbol. Before geometry editing.
        """
        import re
        import warnings
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        ndimen = self.get_dimensionality()
        sg = 'unknown'
        if ndimen == 0:
            rexp = r'^\s*POINT  GROUP  N\.'
        elif ndimen == 1:
            rexp = r'^POLYMER GROUP N\.'
        elif ndimen == 2:
            rexp = r'^\s*TWO\-SIDED PLANE GROUP N\.'
        elif ndimen == 3:
            rexp = r'^\s*SPACE GROUP \(CENTROSYMMETRIC\)'

        for nline, line in enumerate(self.data):
            if re.match(rexp, line):
                if ndimen == 0:
                    sg = line.strip().split(':')[1].split('OR')[1].strip()
                elif ndimen == 1:
                    sg = line.strip().split(':')[1].strip()
                elif ndimen == 2:
                    sg = line.strip().split(':')[1].strip()
                else:
                    sg = line.strip().split(':')[1].strip()
                break

        if sg == 'unknown':
            warnings.warn('Symmstry information lost. Trying to get from pymatgen...',
                          stacklevel=2)
            struc = self.get_geometry(initial=True, write_gui=False)
            sg = SpacegroupAnalyzer(struc).get_space_group_symbol()
        return sg

    @property
    def n_atoms(self):
        """
        Number of atoms. After geometry editing.
        """
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.num_sites

    @property
    def atom_symbols(self):
        """
        Atom symbols. After geometry editing.
        """
        from mendeleev import element

        symbol = []
        for i in self.atom_numbers:
            symbol.append(element(int(i % 100)).symbol)

        return symbol

    @property
    def atom_numbers(self):
        """
        Conventional atom numbers. After geometry editing.
        """
        from CRYSTALpytools.base.output import GeomBASE

        _, conv_z = GeomBASE.read_conv_z(self.data, 0)
        return conv_z

    @property
    def atom_positions(self):
        """
        Composite fractional / Cartesian atomic coordinates. Consistent
        with CRYSTAL definitions. 3D: Fractional; 2D: Frac, Frac, Cart; 1D
        Frac, Cart, Cart; 0D: Cart, Cart, Cart. After geometry editing
        (before optimization).
        """
        import numpy as np

        ndimen = self.get_dimensionality()
        struc = self.get_geometry(initial=True, write_gui=False)

        if ndimen == 0:
            composite_coord = struc.cart_coords.tolist()
        elif ndimen == 3:
            composite_coord = struc.frac_coords.tolist()
        else:
            cart_coord = struc.cart_coords.tolist()
            frac_coord = struc.frac_coords.tolist()
            composite_coord = []
            for i in range(struc.num_sites):
                composite_coord.append(frac_coord[i][:ndimen] + cart_coord[i][ndimen:])

        return np.array(composite_coord, dtype=float)

    @property
    def atom_positions_cart(self):
        """
        Cartesian atomic coordinates. After geometry editing (before optimization).
        """
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.cart_coords

    @property
    def atom_positions_frac(self):
        """
        Fractional atomic coordinates. After geometry editing (before optimization).
        """
        import warnings

        ndimen = self.get_dimensionality()
        if ndimen != 3:
            warnings.warn("Low dimension systems. Property 'atom_positions' is called instead.",
                          stacklevel=2)
            return self.atom_positions
        else:
            struc = self.get_geometry(initial=True, write_gui=False)
            return struc.frac_coords

    def get_trans_matrix(self):
        """
        Get cell transformation matrix

        Returns:
            self.trans_matrix (np.ndarray): 3\*3 array of supercell
                expansion matrix
        """
        import re
        import numpy as np

        ndimen = self.get_dimensionality()
        self.trans_matrix = np.eye(3, dtype=float)
        for i, line in enumerate(self.data[:self.eoo]):
            if re.match(r'^\s+EXPANSION MATRIX OF PRIMITIVE CELL', line):
                mx = np.array([
                    self.data[i + 1].strip().split()[1:],
                    self.data[i + 2].strip().split()[1:],
                    self.data[i + 3].strip().split()[1:],
                ], dtype=float)
                self.trans_matrix[:ndimen, :ndimen] = mx[:ndimen, :ndimen]
                break

        return self.trans_matrix

    def get_primitive_geometry(self, initial=True, write_gui=False, gui_name=None,
                               symmetry='pymatgen', **kwargs):
        """
        Get the primitive geometry, reduced by cell transformation matrix
        inversed.

        .. note::
            This is not the standard 'primitive cell'. This method returns
            geometry before supercell expansion keywords such as 'SUPERCEL'
            or 'SCELPHONO'.

            Conventional atomic numbers are not available.

        Args:
            initial (bool): Read the initial or last geometry. Useful in
                case of geometry optimization.
            write_gui (bool): If True, write .gui file
            gui_name (str): Valid only if ``write_gui = True``. Gui file
                is named as 'gui_name'. If None, use 'basename.gui'. The
                basename is the same as output file.
            symmetry (str): Valid only if ``write_gui = True``. 'pymatgen'
                to use symmetry info from a pymatgen SpacegroupAnalyzer;
                'initial' to use symmstry information on output file. If
                None, no symmstry. Otherwise it is taken from the existing
                gui file.
            **kwargs: Valid only if ``write_gui = True`` and
                ``symmetry = 'pymatgen'``.  Passed to Pymatgen
                SpacegroupAnalyzer object.

        Returns:
            self.primitive_geometry (CStructure | CMolecule): A modified
                pymatgen Structure or molecule object.
        """
        import os
        import re
        import warnings
        import numpy as np
        from CRYSTALpytools.convert import cry_pmg2gui
        from CRYSTALpytools.crystal_io import Crystal_gui

        ndimen = self.get_dimensionality()
        self.get_geometry(initial=initial, write_gui=False)
        self.get_trans_matrix()

        if ndimen == 0:
            warnings.warn('0D system. Nothing to reduce.', stacklevel=2)
            self.primitive_geometry = self.geometry
            return self.primitive_geometry

        shrink_mx = np.linalg.inv(self.trans_matrix)
        pstruc = self.geometry.get_pcel(self.trans_matrix)

        self.primitive_geometry = pstruc
        # Write gui files
        if write_gui == True:
            if np.all(gui_name==None):
                gui_name = os.path.splitext(self.name)[0]
                gui_name = '{}.gui'.format(gui_name)

            if symmetry == 'pymatgen':
                gui = cry_pmg2gui(pstruc, gui_file=gui_name, symmetry=True, **kwargs)
            elif np.all(symmetry==None):
                gui = cry_pmg2gui(pstruc, gui_file=gui_name, symmetry=False)
            elif symmetry == 'initial':
                self.get_symmops()
                gui = cry_pmg2gui(pstruc, gui_file=None, symmetry=False)
                gui.symmops = self.symmops
                gui.n_symmops = self.n_symmops
                gui.space_group = self.sg_number
                gui.write_gui(gui_name, symm=True)
            else:
                warnings.warn('Symmetry adapted from reference geometry. Make sure that is desired.',
                              stacklevel=2)
                gui_ref = Crystal_gui().read_gui(symmetry)
                gui = cry_pmg2gui(pstruc, gui_file=None, symmetry=False)
                # Replace the symmops with the reference file
                gui.symmops = gui_ref.symmops
                gui.n_symmops = gui_ref.n_symmops
                gui.space_group = gui_ref.space_group
                gui.write_gui(gui_name, symm=True)

        return self.primitive_geometry

    def get_primitive_lattice(self, initial=True):
        """
        Returns the primitive lattice of the system reduced by cell
        transformation matrix inverse. Unit: Angstrom.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.primitive_lattice (np.ndarray): Lattice of the system.
        """
        import warnings

        ndimen = self.get_dimensionality()
        self.primitive_lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.primitive_lattice

        self.get_primitive_geometry(initial=initial, write_gui=False)
        self.primitive_lattice = self.primitive_geometry.lattice.matrix

        return self.primitive_lattice

    def get_primitive_reciprocal_lattice(self, initial=False):
        """
        Returns the primitive reciprocal lattice of the system before
        expansion by cell transformation matrix inverse. Unit: Angstrom^-1.

        Args:
            initial (bool): Read the initial or last lattice. Useful in
                case of geometry optimization.
        Returns:
            self.primitive_reciprocal_lattice (np.ndarray): Lattice of the system.
        """
        import warnings

        ndimen = self.get_dimensionality()
        self.primitive_reciprocal_lattice = None
        if ndimen == 0:
            warnings.warn('0D system. No lattice.')
            return self.primitive_reciprocal_lattice

        self.get_primitive_lattice(initial=initial)
        self.primitive_reciprocal_lattice = self.primitive_geometry.lattice.reciprocal_lattice.matrix

        return self.primitive_reciprocal_lattice

    def get_config_analysis(self, return_multiplicity=False):
        """
        Return the configuration analysis for solid solutions (CONFCON keyword in input)

        Args:
            return_multiplicity (bool, optional): Return multiplicity information. Defaults to False.
        Returns:
            list or str: Configuration analysis if available, or a warning message
        """
        import re

        import numpy as np

        # Check this is a configuration analysis calculation
        try:
            begin = self.data.index(
                '                             CONFIGURATION ANALYSIS\n')
        except:
            return "WARNING: this is not a CONFCNT analysis."

        for line in self.data[::-1]:
            if '----' in line:
                dash_line = line.rstrip().lstrip()
                break

        for i, line in enumerate(self.data[begin:]):
            if re.match(r'^ COMPOSITION', line):
                self.n_classes = line.split()[9]
                original_atom = str(line.split()[2])
                begin = begin+i
        config_list = []

        # Read all the configurations
        for line in self.data[begin:]:
            if not re.match(r'^   WARNING', line):
                config_list.extend(line.split())

        '''multiplicity = []

        for i in range(len(config_list)):
            if config_list[i] == 'MULTIPLICITY' and i < len(config_list) - 1:
                number = re.findall(r'\d+', config_list[i+1])
                if number:
                    config_list.append(int(number[0]))'''

        config_list = np.array(config_list)
        warning = np.where(config_list == 'WARNING')
        config_list = np.delete(config_list, warning)
        atom1_begin = np.where(config_list == original_atom)[0]
        atom1_end = np.where(
            config_list == dash_line)[0]
        atom2_begin = np.where(config_list == 'XX')[0]
        atom2_end = np.where(
            config_list == '<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')[0]
        end = np.where(
            config_list == '===============================================================================')[0][-1]
        atom2_end = np.append(atom2_end, end)
        atom_type1 = []
        atom_type2 = []
        if return_multiplicity == True:
            input_string = ' '.join(config_list.tolist())
            matches = re.findall(r'MULTIPLICITY\s*:\s*(\d+)', input_string)

            # multiplicity_tmp = config_list[np.where(config_list == 'MULTIPLICITY')[0]+1]
            multiplicity = [int(x) for x in matches]
            self.multiplicity = multiplicity
        config_list = config_list.tolist()
        for i in range(len(atom1_end)):
            atom_type1.append(
                [int(x) for x in config_list[atom1_begin[i+1]+1:atom1_end[i]]])
            atom_type2.append(
                [int(x) for x in config_list[atom2_begin[i]+1:atom2_end[i]]])

        self.atom_type1 = atom_type1
        self.atom_type2 = atom_type2

        if return_multiplicity == True:
            return [self.atom_type1, self.atom_type2, multiplicity]
        else:
            return [self.atom_type1, self.atom_type2]

    #### SCF / OPT Convergence ####

    def get_convergence(self, history=False):
        """
        The upper level of ``get_scf_convergence`` and ``get_opt_convergence``.
        For analysing the geometry and energy convergence.

        .. note::

            It might not work well with SCF / OPT cycles of multiple systems
            such as PREOPTGEOM + EOS calculations

        Args:
            history (bool): If true, the convergence history of optimisation
                (energy,gradient, displacement) / SCF is returned.

        Returns:
            self (Crystal_output): New attributes listed below
            self.final_energy (float) The converged energy of SCF / Opt. Unit: eV

        For other attributes, see ``get_scf_convergence`` and
        ``get_opt_convergence`` methods on the same page.
        """
        import re
        import warnings

        self.final_energy = None
        for line in self.data[self.eoo::-1]:
            if re.match(r'\s\W OPT END - CONVERGED', line) != None:
                data = line.strip().split()
                self.final_energy = units.H_to_eV(float(data[7]))
                self.opt_cycles = int(data[-2])
                if re.search('CONVERGED', line):
                    self.opt_status = 'converged'
                elif re.search('FAILED', line):
                    self.opt_status = 'failed'
                else:
                    warnings.warn('Unknown termination.', stacklevel=2)
                    self.opt_status = 'unknown'
                is_scf = False
                break

            elif re.match(r'^ == SCF ENDED', line) != None:
                data = line.strip().split()
                self.final_energy = units.H_to_eV(float(data[8]))
                self.scf_cycles = int(data[-1])
                if re.search('CONVERGENCE', line):
                    self.scf_status = 'converged'
                elif re.search('TOO MANY CYCLES', line):
                    self.scf_status = 'too many cycles'
                else:
                    warnings.warn('Unknown termination.', stacklevel=2)
                    self.scf_status = 'unknown'
                is_scf = True
                break

        if np.all(self.final_energy==None):
            warnings.warn('No final energy found in the output file. self.final_energy = None',
                          stacklevel=2)

        if history == True:
            if is_scf == True:
                self.get_scf_convergence(all_cycles=False)
            else:
                self.get_opt_convergence()

        return self

    def get_final_energy(self):
        """
        Get the final energy of the system. A wrapper of ``self.get_convergence``.

        Returns:
            self.final_energy (float): The final energy of the system.
        """
        self.get_convergence(history=False)

        return self.final_energy

    #### SCF ####

    def get_scf_convergence(self, all_cycles=False):
        """
        Returns the scf convergence energy and energy difference. A wrapper of
        ``CRYSTALpytools.base.SCFBASE.read_convergence``.

        Args:
            all_cycles (bool, optional): Return all SCF steps for a geometry
                opt. The 'ONELOG' CRYSTAL keyword is needed.

        Returns:
            self (Crystal_output): New attributes listed below
            self.scf_cycles (int | array): Number of cycles. Array if
                ``all_cycles=True``.
            self.scf_status (str | list): 'terminated', 'converged',
                'too many cycles' and 'unknown'. List if ``all_cycles=True``.
            self.scf_energy (array): SCF energy convergence. Unit: eV
            self.scf_deltae (array): Energy difference. Unit: eV
        """
        import numpy as np
        from CRYSTALpytools.base.output import SCFBASE

        if all_cycles == True:
            self.scf_cycles = []
            self.scf_status = []
            self.scf_energy = []
            self.scf_deltae = []

        countline = 0
        while countline < self.eoo:
            output = SCFBASE.read_convergence(self.data[:self.eoo], countline)
            if all_cycles == False:
                self.scf_cycles = output[1]
                self.scf_status = output[2]
                self.scf_energy = output[3]
                self.scf_deltae = output[4]
                break
            else:
                countline = output[0]
                self.scf_cycles.append(output[1])
                self.scf_status.append(output[2])
                self.scf_energy.append(output[3])
                self.scf_deltae.append(output[4])
                countline += 1

        if all_cycles == True:
            self.scf_cycles = np.array(self.scf_cycles, dtype=int)
            self.scf_energy = np.array(self.scf_energy)
            self.scf_deltae = np.array(self.scf_deltae)

        return self

    def get_fermi_energy(self, history=False):
        """
        Returns the system Fermi energy.

        Args:
            history (bool): Whether to read the convergence history of Fermi energy.

        Returns:
            self.fermi_energy (float | array): Fermi energy of the system. For
                spin-polarized insulating systems, ``self.fermi_energy`` would
                be either a 2\*1 array (``history=False``) or a nCYC\*2 array
                (``history=True``).
        """
        from CRYSTALpytools.base.output import SCFBASE

        output = SCFBASE.read_fermi_energy(self.data[:self.eoo], self.eoo - 1, history=history)
        self.spin_pol = output[1]
        self.fermi_energy = output[2]

        return self.fermi_energy

    def get_band_gap(self, history=False):
        """
        Returns the system band gap.

        Args:
            history (bool): Whether to read the convergence history of band gap.

        Returns:
            self.band_gap (float | array): Band gap of the system. For
                spin-polarized systems, ``self.band_gap`` would be either a
                2\*1 array (``history=False``) or a nCYC\*2 array
                (``history=True``).
        """
        from CRYSTALpytools.base.output import SCFBASE

        output = SCFBASE.read_band_gap(self.data[:self.eoo], self.eoo - 1, history=history)
        self.spin_pol = output[1]
        self.band_gap = output[2]

        return self.band_gap

    def get_mulliken_charges(self):
        """
        Return the atomic Mulliken charges (PPAN keyword in input).

        Returns:
            self.mulliken_charges (array): natom\*1 for non spin-polarised systems.
                natom\*3 for spin-polarised systems. [total, :math:`\\alpha`, :math:`\\beta`].
        """
        import re
        import warnings
        import numpy as np

        mulliken = []  # empty, 1*1 or 2*1 list
        countline = 0
        countm = 0
        while countline < self.eoo:
            line = self.data[countline]
            if re.match(r'\s*MULLIKEN POPULATION ANALYSIS', line):
                mulliken_charge = []  # natom*1
                countline += 4
                line2 = self.data[countline]
                while len(line2.strip()) != 0:
                    if re.match(r'^\s+[0-9]+\s+[A-Z, a-z]+\s+[0-9+]', line2):
                        data = line2.strip().split()
                        mulliken_charge.append(data[3])
                    countline += 1
                    line2 = self.data[countline]

                mulliken.append(mulliken_charge)
                continue
            else:
                countline += 1

        if len(mulliken) == 0:
            warnings.warn('Mulliken analysis not found.', stacklevel=2)
            self.mulliken_charges = np.array([], dtype=float)
        elif len(mulliken) == 1:
            self.mulliken_charges = np.array(mulliken[0], dtype=float)
            self.spin_pol = False
        else:
            apb = np.array(mulliken[0], dtype=float)
            amb = np.array(mulliken[1], dtype=float)
            self.mulliken_charges = np.array([apb, (apb + amb) / 2, (apb - amb) / 2])
            self.spin_pol = True

        return self.mulliken_charges

    #### Optimization ####

    def get_opt_convergence(self, primitive=False, scf_history=False,
                            write_gui=False, gui_name=None,
                            symmetry='pymatgen', **kwargs):
        """
        Returns optimisation convergence. A wrapper of
        ``CRYSTALpytools.base.OptBASE.read_convergence``.

        Args:
            primitive (bool): Restore the primitive cell (multiply by the
                cell transform matrix inversed)
            scf_history (bool): Read SCF history of each optimisation step.
                Keyword 'ONELOG' is needed. Please refer to
                ``self.get_scf_convergence(all_cycles=True)`` method.
            write_gui (bool): If True, write .gui file of each step
            gui_name (str): Valid only if ``write_gui = True``. Gui file
                is named as 'gui_name-optxxx.gui'. If None, use
                'basename-optxxx.gui'. The basename is the same as output.
            symmetry (str): Valid only if ``write_gui = True``. 'pymatgen'
                to use symmetry info from a pymatgen SpacegroupAnalyzer;
                'initial' to use symmstry information on output file. If
                None, no symmstry. Otherwise it is taken from the existing
                gui file.
            **kwargs: Valid only if ``write_gui = True`` and
                ``symmetry = 'pymatgen'``.  Passed to Pymatgen
                SpacegroupAnalyzer object.

        Returns:
            self (Crystal_output): New attributes listed below
            self.opt_cycles (int): Number of cycles.
            self.opt_status (str): 'terminated', 'converged', 'failed' and 'unknown'
            self.opt_energy (array): Total energy convergence. Unit: eV
            self.opt_deltae (array): Total energy difference. Unit: eV
            self.opt_geometry (list): Modified CStructure/CMolecule at each step.
            self.opt_maxgrad (array): Maximum gradient convergence. Unit: Hartree/Bohr
            self.opt_rmsgrad (array): RMS gradient convergence. Unit: Hartree/Bohr
            self.opt_maxdisp (array): Maximum displacement convergence. Unit: Bohr
            self.opt_rmsdisp (array): RMS displacement convergence. Unit: Bohr
        """
        from CRYSTALpytools.crystal_io import Crystal_gui
        from CRYSTALpytools.convert import cry_pmg2gui
        from CRYSTALpytools.base.output import OptBASE
        from pymatgen.core.lattice import Lattice
        import numpy as np
        import re
        import os

        ndimen = self.get_dimensionality()
        countline = 0

        # Initial geometry
        struc0 = self.get_geometry(initial=True, write_gui=False)

        # Initial SCF energy. No SCF when optimisation is restarted
        e0 = None
        self.opt_cycles = 0
        lattice_line = -1
        while countline < self.eoo:
            line = self.data[countline]
            # Initial step SCF
            if re.match(r'^\s*== SCF ENDED', line):
                line_data = line.strip().split()
                e0 = float(line_data[8])
                countline += 1
            # Initial step SCF empirical corrections
            elif re.match(r'^\s*TOTAL ENERGY \+', line):
                line_data = line.strip().split()
                e0 = float(line_data[-1])
                countline += 1
            # Initial Gradient
            elif re.match(r'^\s+MAX GRADIENT', line):
                line_data = line.strip().split()
                maxg0 = float(line_data[2])
                countline += 1
            elif re.match(r'^\s+RMS GRADIENT', line):
                line_data = line.strip().split()
                rmsg0 = float(line_data[2])
                countline += 1
            # Enter the Opt block
            elif re.match(r'^\s*OPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPT', line):
                output = OptBASE.read_optblock(self.data[:self.eoo], countline)
                self.opt_cycles = output[1]
                self.opt_status = output[2]
                self.opt_energy = output[3]
                self.opt_deltae = output[4]
                self.opt_geometry = output[5]
                self.opt_maxgrad = output[6]
                self.opt_rmsgrad = output[7]
                self.opt_maxdisp = output[8]
                self.opt_rmsdisp = output[9]
                break
            # Lattice matrix rotation issue
            elif re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN', line):
                lattice_line = countline
                countline += 1
            else:
                countline += 1

        # Converged at initial step, or wrong file
        if self.opt_cycles == 0:
            if np.all(e0==None):
                raise Exception('Valid data not found.')
            else:
                self.opt_cycles += 1
                self.opt_status = 'converged at initial step'
                self.opt_energy = np.array([e0,])
                self.opt_deltae = np.array([e0,])
                self.opt_geometry = [struc0,]
                self.opt_maxgrad = np.array([maxg0,])
                self.opt_rmsgrad = np.array([rmsg0,])
                self.opt_maxdisp = np.array([])
                self.opt_rmsdisp = np.array([])
        else:
            # restarted from previous opt
            if np.all(e0==None):
                pass
            else:
                self.opt_cycles += 1
                self.opt_energy = np.concatenate([[e0,], self.opt_energy])
                self.opt_deltae = np.concatenate([[e0,], self.opt_deltae])
                self.opt_geometry = [struc0,] + self.opt_geometry
                self.opt_maxgrad = np.concatenate([[maxg0,], self.opt_maxgrad])
                self.opt_rmsgrad = np.concatenate([[rmsg0,], self.opt_rmsgrad])
                self.opt_maxdisp = np.concatenate([[0.,], self.opt_maxdisp])
                self.opt_rmsdisp = np.concatenate([[0.,], self.opt_rmsdisp])

        # Lattice matrix rotation issue
        if ndimen != 0 and lattice_line != -1:
            # Reference lattice vector a
            a_crys = np.array(self.data[lattice_line + 2].strip().split(), dtype=float)
            # Pmg lattice vector a
            a_pmg = struc0.lattice.matrix[0, :]
            # Rotate lattice: CStructure.rot_cel
            for idx_s, s in enumerate(self.opt_geometry):
                self.opt_geometry[idx_s] = s.rot_cel(a_crys, a_pmg)

        # Get primitive cell
        if primitive == True:
            ndimen = self.get_dimensionality()
            if ndimen == 0:
                warnings.warn('0D system. Nothing to reduce.', stacklevel=2)
            else:
                self.get_trans_matrix()
                shrink_mx = np.linalg.inv(self.trans_matrix)
                for idx_s, s in enumerate(self.opt_geometry):
                    self.opt_geometry[idx_s] = s.get_pcel(self.trans_matrix)

        # SCF history
        if scf_history == True:
            self.get_scf_convergence(all_cycles=True)

        # Write gui files
        if write_gui == True:
            if np.all(gui_name==None):
                gui_name = os.path.splitext(self.name)[0]
            gui_list = ['{}-opt{:0=3d}.gui'.format(gui_name, i+1) for i in range(self.opt_cycles)]

            if symmetry == 'pymatgen':
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=gui_list[idx_s],
                                      symmetry=True, **kwargs)
            elif np.all(symmetry==None):
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=gui_list[idx_s], symmetry=False)
            elif symmetry == 'initial':
                self.get_symmops()
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=None, symmetry=False)
                    gui.symmops = self.symmops
                    gui.n_symmops = self.n_symmops
                    gui.space_group = self.sg_number
                    gui.write_gui(gui_list[idx_s], symm=True)
            else:
                warnings.warn('Symmetry adapted from reference geometry. Make sure that is desired.',
                              stacklevel=2)
                gui_ref = Crystal_gui().read_gui(symmetry)
                for idx_s, s in enumerate(self.opt_geometry):
                    gui = cry_pmg2gui(s, gui_file=None, symmetry=False)
                    # Replace the symmops with the reference file
                    gui.symmops = gui_ref.symmops
                    gui.n_symmops = gui_ref.n_symmops
                    gui.space_group = gui_ref.space_group
                    gui.write_gui(gui_list[idx_s], symm=True)

        return self

    def get_forces(self, initial=True, grad=False):
        """
        Read forces.

        Args:
            initial (bool): Return forces from the initial calculation. If
                ``initial=False``, return to the last forces, which is valid
                only for geometry optimizations and the keyword 'ONELOG' is
                needed in d12 file.
            grad (bool): Return gradient convergence history. For optimizations
                only.

        Returns:
            self (Crystal_output): New attributes listed below
            self.forces_atoms (array): natom\*3 array. Atomic forces. Unit: Hartree/Bohr
            self.forces_cell (array): 3\*3 array. Cell forces, 3D only. Unit: Hartree/Bohr
            self.opt_maxgrad (array): Maximum gradient convergence. Unit: Hartree/Bohr
            self.opt_rmsgrad (array): RMS gradient convergence. Unit: Hartree/Bohr
        """
        import warnings, re
        import numpy as np

        if initial == False or grad == True:
            if ' OPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPT\n' not in self.data:
                warnings.warn('Not a geometry optimisation: Set initial = True and grad = False', stacklevel=2)
                initial = True
                grad = False

        self.forces_atoms = []
        self.forces_cell = []

        if grad == True:
            self.get_opt_convergence()

        if initial == True:
            for i, line in enumerate(self.data):
                if re.match(r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)', line):
                    for j in range(i+2, i+2+self.n_atoms):
                        self.forces_atoms.append(self.data[j].strip().split()[2:])
                    self.forces_atoms = np.array(self.forces_atoms, dtype=float)

                if re.match(r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR', line):
                    for j in range(i+4, i+7):
                        self.forces_cell.append(self.data[j].strip().split())
                    self.forces_cell = np.array(self.forces_cell, dtype=float)

        else:
            for i, line in enumerate(self.data[::-1]):
                if re.match(r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR', line):
                    for j in range(len(self.data)-i+3, len(self.data)-i+6):
                        self.forces_cell.append(self.data[j].strip().split())
                    self.forces_cell = np.array(self.forces_cell, dtype=float)

                if re.match(r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)', line):
                    for j in range(len(self.data)-i+1, len(self.data)-i+1+self.n_atoms):
                        self.forces_atoms.append(self.data[j].strip().split()[2:])
                    self.forces_atoms = np.array(self.forces_atoms, dtype=float)

        return self

    #### Lattice Dynamics ####

    def get_phonon(self, read_eigvt=False, rm_imaginary=True, rm_overlap=True,
                   imaginary_tol=-1e-4, q_overlap_tol=1e-4, eigvt_amplitude=1.):
        """
        Read phonon-related properties from output file.

        Args:
            read_eigvt (bool): Whether to read phonon eigenvectors and
                normalize it to 1.
            rm_imaginary (bool): Remove the modes with negative frequencies and
                set all the related properties to NaN.
            rm_overlap (bool): *For dispersion calculations* Remove repeated q
                points and recalculate their weights.
            imaginary_tol (float): *``rm_imaginary`` = True only* The threshold
                of negative frequencies.
            q_overlap_tol (float): *``rm_overlap`` = True only* The threshold of
                overlapping points, defined as the 2nd norm of the difference
                of fractional q vectors
            eigvt_amplitude (float|str): Normalize the eigenvector to a certian
                amplitude. Either a number or 'classical' (classical amplitude).

        Returns:
            self (Crystal_output): New attributes listed below
            self.edft (array[float]): :math:`E_{0}` Energy with empirical
                correction. Unit: kJ/mol.
            self.nqpoint (int): Number of q points
            self.qpoint (list[list[array[float], float]]): A nqpoint\*1 list of
                2\*1 list whose first element is a 3\*1 array of fractional
                coordinates and the second is its weight.
            self.nmode (array[int]): Number of modes at q point. nqpoint\*1
                array.
            self.frequency (array[float]): nqpoint\*nmode array ofvibrational
                frequency. Unit: THz
            self.intens (array[float]): nqpoint\*nmode array of harmonic
                intensiy. Unit: km/mol
            self.IR (array[bool]): nqpoint\*nmode array of boolean values
                specifying whether the mode is IR active
            self.Raman (array[bool]): nqpoint\*nmode array of boolean values
                specifying whether the mode is Raman active
            self.eigenvector (array[complex]): *``read_eigvt = True only``*
                nqpoint\*nmode\*natom\*3 array of eigenvectors. Normalized to 1.

        .. note::

            In QHA calculations, ``self.nqpoint`` refer to harmonic phonons
            computed. In other cases it refers to actual q points in reciprocal
            space.
        """
        import re
        import numpy as np
        from CRYSTALpytools.base.output import PhononBASE
        from CRYSTALpytools.units import H_to_kjmol

        is_freq = False
        found_anti = True
        self.edft = []
        self.nqpoint = 0
        self.qpoint = []
        self.nmode = []
        self.frequency = []
        self.intens = []
        self.IR = []
        self.Raman = []
        self.eigenvector = []

        countline = 0
        while countline < self.eoo:
            line = self.data[countline]
            # Whether is a frequency file
            if re.match(r'^\s*\+\+\+\sSYMMETRY\sADAPTION\sOF\sVIBRATIONAL\sMODES\s\+\+\+', line):
                is_freq = True
                countline += 1
                continue
            # E_0 with empirical corrections
            elif re.match(r'^\s+CENTRAL POINT', line):
                self.edft.append(float(line.strip().split()[2]))
                countline += 1
                continue
            # Q point info + frequency
            ## Dispersion
            elif re.match(r'^.+EXPRESSED IN UNITS\s+OF DENOMINATOR', line):
                shrink = int(line.strip().split()[-1])
                countline += 1
                continue
            elif re.match(r'\s+DISPERSION K POINT NUMBER', line):
                coord = np.array(line.strip().split()[7:10], dtype=float)
                weight = float(line.strip().split()[-1])
                self.qpoint.append([coord / shrink, weight])
                self.nqpoint += 1
                countline += 2
                # Read phonons
                phonon = PhononBASE.readmode_basic(self.data[:self.eoo], countline)
                countline = phonon[0]
                self.frequency.append(phonon[1])
                self.intens.append(phonon[2])
                self.IR.append(phonon[3])
                self.Raman.append(phonon[4])
            ## Gamma point
            elif re.match(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+IRREP', line) and self.nqpoint == 0:
                countline += 2
                # Read phonons
                phonon = PhononBASE.readmode_basic(self.data[:self.eoo], countline)
                countline = phonon[0]
                self.frequency.append(phonon[1])
                self.intens.append(phonon[2])
                self.IR.append(phonon[3])
                self.Raman.append(phonon[4])
            ## Phonon eigenvector
            ### Gamma point: real numbers. Imaginary = 0
            elif re.match(r'^\s+NORMAL MODES NORMALIZED', line):
                if read_eigvt == False:
                    countline += 1
                    continue
                countline += 2
                eigvt = PhononBASE.readmode_eigenvector(self.data[:self.eoo], countline)
                countline = eigvt[0]
                self.eigenvector.append(eigvt[1] + 0.j)
            ### Dispersion: complex numbers
            elif re.match(r'^\s+MODES IN PHASE', line):
                if read_eigvt == False:
                    countline += 1
                    continue
                if found_anti == False:  # Real k point
                    self.eigenvector.append(tmp_eigvt)
                countline += 2
                found_anti = False
                eigvt = PhononBASE.readmode_eigenvector(self.data[:self.eoo], countline)
                countline = eigvt[0]
                tmp_eigvt = eigvt[1] + 0.j
            elif re.match(r'^\s+MODES IN ANTI\-PHASE', line):
                if read_eigvt == False:
                    countline += 1
                    continue
                countline += 2
                found_anti = True
                eigvt_anti = PhononBASE.readmode_eigenvector(self.data[:self.eoo], countline)
                countline = eigvt_anti[0]
                self.eigenvector.append(tmp_eigvt + eigvt_anti[1] * 1.j)
            # Other data
            else:
                countline += 1
                continue

        if is_freq == False:
            raise Exception('Not a frequency calculation.')
        if found_anti == False and read_eigvt == True:  # The last real k point
            self.eigenvector.append(tmp_eigvt)

        # Format data
        # HA/QHA Gamma point calculation
        if self.nqpoint == 0:
            self.nqpoint = len(self.edft)
            self.qpoint = [[np.zeros([3,]), 1.] for i in range(self.nqpoint)]
            if len(self.edft) == 1:
                self.edft = [self.edft[0] for i in range(self.nqpoint)]
        # Dispersion
        else:
            self.qpoint = [[i[0], i[1] / self.nqpoint] for i in self.qpoint]
            self.edft = [self.edft[0] for i in range(self.nqpoint)]

        self.edft = H_to_kjmol(np.array(self.edft))

        self.frequency = np.array(self.frequency)
        self.nmode = np.array([len(i) for i in self.frequency], dtype=int)
        if self.intens[0] == []:
            self.intens = []
            self.IR = []
            self.Raman = []
        else:
            self.intens = np.array(self.intens)

        if self.eigenvector != []:
            self.eigenvector = np.array(self.eigenvector)
            # already normalised to classical amplitude
            if str(eigvt_amplitude).lower() == 'classical':
                pass
            # remove classical amplitude
            # elif str(eigvt_amplitude).lower() == 'classical-rev':
            #     struc = self.get_geometry(initial=False, write_gui=False)
            # To a specific value
            else:
                for idx_q in range(self.nqpoint):
                    self.eigenvector[idx_q] = PhononBASE.normalize_eigenvector(
                        self.eigenvector[idx_q],
                        amplitude=float(eigvt_amplitude),
                    )

        if rm_imaginary == True:
            self = PhononBASE.clean_imaginary(self, threshold=imaginary_tol)

        if rm_overlap == True and self.nqpoint > 1:
            self = PhononBASE.clean_q_overlap(self, threshold=q_overlap_tol)

        return self

    def get_anh_spectra(self):
        """
        This method reads data from a CRYSTAL output file and processes it to 
        extract anharmonic (VSCF and VCI) IR and Raman spectra.

        Returns:
            self.IR_HO_0K (array[float]): 2D numpy array containing harmonic IR
                frequency and intensities computed at 0 K. 
            self.IR_HO_T (array[float]): 2D numpy array containing harmonic IR
                frequency and intensities computed at temperature T. 
            self.IR_VSCF_0K (array[float]): 2D numpy array containing VSCF IR
                frequency and intensities computed at 0 K. 
            self.IR_VSCF_T (array[float]): 2D numpy array containing VSCF IR 
                frequency and intensities computed at temperature T.
            self.IR_VCI_0K (array[float]): 2D numpy array containing VCI IR
                frequency and intensities computed at 0 K. 
            self.IR_VCI_T (array[float]): 2D numpy array containing VCI IR 
                frequency and intensities computed at temperature T.
            self.Ram_HO_0K_tot (array[float]): 2D numpy array containing harmonic
                Raman frequency and intensities (total) computed at 0 K. 
            self.Ram_HO_0K_per (array[float]): 2D numpy array containing harmonic
                Raman frequency and intensities (perpendicular component ) computed at 
                temperature 0 K. 
            self.Ram_HO_0K_par (array[float]): 2D numpy array containing harmonic
                Raman frequency and intensities (parallel component ) computed at 
                temperature 0 K. 
            self.Ram_HO_T_tot (array[float]): 2D numpy array containing harmonic
                Raman frequency and intensities (total) computed at temperature T. 
            self.Ram_HO_T_per (array[float]): 2D numpy array containing harmonic
                Raman frequency and intensities (perpendicular component ) computed at 
                temperature T. 
            self.Ram_HO_T_par (array[float]): 2D numpy array containing harmonic
                Raman frequency and intensities (parallel component ) computed at 
                temperature T. 
            self.Ram_VSCF_0K_tot (array[float]): 2D numpy array containing VSCF
                Raman frequency and intensities (total) computed at 0 K. 
            self.Ram_VSCF_0K_per (array[float]): 2D numpy array containing VSCF
                Raman frequency and intensities (perpendicular component) computed at 
                0 K. 
            self.Ram_VSCF_0K_par (array[float]): 2D numpy array containing VSCF
                Raman frequency and intensities (parallel component) computed at 
                0 K. 
            self.Ram_VSCF_T_tot (array[float]): 2D numpy array containing VSCF
                Raman frequency and intensities (total) computed at temperature T. 
            self.Ram_VSCF_T_per (array[float]): 2D numpy array containing VSCF
                Raman frequency and intensities (perpendicular component) computed at 
                temperature T. 
            self.Ram_VSCF_T_par (array[float]): 2D numpy array containing VSCF
                Raman frequency and intensities (parallel component) computed at 
                temperature T. 
            self.Ram_VCI_0K_tot (array[float]): 2D numpy array containing VCI
                Raman frequency and intensities (total) computed at 0 K. 
            self.Ram_VCI_0K_per (array[float]): 2D numpy array containing VCI
                Raman frequency and intensities (perpendicular component) computed at 
                0 K. 
            self.Ram_VCI_0K_par (array[float]): 2D numpy array containing VCI
                Raman frequency and intensities (parallel component) computed at 
                0 K. 
            self.Ram_VCI_T_tot (array[float]): 2D numpy array containing VCI
                Raman frequency and intensities (total) computed at temperature T. 
            self.Ram_VCI_T_per (array[float]): 2D numpy array containing VCI
                Raman frequency and intensities (perpendicular component) computed at 
                temperature T. 
            self.Ram_VCI_T_par (array[float]): 2D numpy array containing VCI
                Raman frequency and intensities (parallel component) computed at 
                temperature T.
            self.Ram_HO_0K_comp_xx (array[float]): 2D numpy array containing
                harmonic Raman frequency and intensities (xx component) computed at 0
                K.
            self.Ram_HO_T_comp_xx (array[float]): 2D numpy array containing
                harmonic Raman frequency and intensities (xx component) computed at
                temperature T.
            self.Ram_VSCF_0K_comp_xx (array[float]): 2D numpy array containing
                VSCF Raman frequency and intensities (xx component) computed at 0 K.
            self.Ram_VSCF_T_comp_xx (array[float]): 2D numpy array containing
                VSCF Raman frequency and intensities (xx component) computed at
                temperature T.
            self.Ram_VCI_0K_comp_xx (array[float]): 2D numpy array containing
                VCI Raman frequency and intensities (xx component) computed at 0 K.
            self.Ram_VCI_T_comp_xx (array[float]): 2D numpy array containing
                VCI Raman frequency and intensities (xx component) computed at
                temperature T.

        .. note::

            Please, note that for the sake of brevity, only the xx Raman component
            attributes have been listed here, but the yy, zz, xy, xz, yz components
            are available as well.  
        """
        import re
        import numpy as np

        # Initialize some logical variables
        save = False
        HO_IR = False
        HO_Ram_0K_tot = False
        HO_Ram_T_tot = False
        HO_Ram_0K_comp = False
        HO_Ram_T_comp = False
        VSCF_IR = False
        VSCF_Ram_0K_tot = False
        VSCF_Ram_0K_comp = False
        VSCF_Ram_T_tot = False
        VSCF_Ram_T_comp = False
        VCI_IR = False
        VCI_Ram_0K_tot = False
        VCI_Ram_0K_comp = False
        VCI_Ram_T_tot = False
        VCI_Ram_T_comp = False

        # Initialize some member variables
        IR_HO = []
        IR_VSCF = []
        IR_VCI = []
        Ram_HO_0K_tot = []
        Ram_HO_T_tot = []
        Ram_HO_0K_comp = []
        Ram_HO_T_comp = []
        Ram_VSCF_0K_tot = []
        Ram_VSCF_T_tot = []
        Ram_VSCF_0K_comp = []
        Ram_VSCF_T_comp = []
        Ram_VCI_0K_tot = []
        Ram_VCI_T_tot = []
        Ram_VCI_0K_comp = []
        Ram_VCI_T_comp = []

        # Initialize some buffers
        bufferHO_IR = []
        bufferHO_Ram_0K_tot = []
        bufferHO_Ram_T_tot = []
        bufferHO_Ram_0K_comp = []
        bufferHO_Ram_T_comp = []
        bufferVSCF_IR = []
        bufferVSCF_Ram_0K_tot = []
        bufferVSCF_Ram_T_tot = []
        bufferVSCF_Ram_0K_comp = []
        bufferVSCF_Ram_T_comp = []
        bufferVCI_IR = []
        bufferVCI_Ram_0K_tot = []
        bufferVCI_Ram_T_tot = []
        bufferVCI_Ram_0K_comp = []
        bufferVCI_Ram_T_comp = []

        # Big loop over lines of CRYSTAL output file -->
        for i, line in enumerate(self.data):

            if re.match(r'\s*HARMONIC IR SPECTRUM', line):
                if re.match(r'\s*HO', self.data[i-1]):
                    HO_IR = True
                else:
                    print('Something went wrong with your CRYSTAL output file.')
                save = True

            if re.match(r'\s*ANHARMONIC IR SPECTRUM', line):
                if re.match(r'\s*VSCF', self.data[i-1]):
                    VSCF_IR = True
                elif re.match(r'\s*VCI*', self.data[i-1]):
                    VCI_IR = True
                else:
                    print('Something went wrong with your CRYSTAL output file.')
                save = True

            if re.match(r'\s*HARMONIC RAMAN SPECTRUM', line):
                if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                    if re.match(r'\s*I_TOT', self.data[i+4]):
                        HO_Ram_0K_tot = True
                    else:
                        HO_Ram_0K_comp = True
                else:
                    if re.match(r'\s*I_TOT', self.data[i+4]):
                        HO_Ram_T_tot = True
                    else:
                        HO_Ram_T_comp = True
                save = True

            if re.match(r'\s*ANHARMONIC RAMAN SPECTRUM', line):
                if re.match(r'\s*VSCF', self.data[i-1]):
                    if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VSCF_Ram_0K_tot = True
                        else:
                            VSCF_Ram_0K_comp = True
                    else:
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VSCF_Ram_T_tot = True
                        else:
                            VSCF_Ram_T_comp = True
                save = True
                if re.match(r'\s*VCI*', self.data[i-1]):
                    if re.match(r'\s*\[ 0 K \]', self.data[i+1]):
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VCI_Ram_0K_tot = True
                        else:
                            VCI_Ram_0K_comp = True
                    else:
                        if re.match(r'\s*I_TOT', self.data[i+4]):
                            VCI_Ram_T_tot = True
                        else:
                            VCI_Ram_T_comp = True
                save = True

            if re.match(r'\s*HHHHHHHHHHHHH', line):
                save = False
                HO_IR = False
                HO_Ram_0K_tot = False
                HO_Ram_T_tot = False
                HO_Ram_0K_comp = False
                HO_Ram_T_comp = False
                VSCF_IR = False
                VSCF_Ram_0K_tot = False
                VSCF_Ram_0K_comp = False
                VSCF_Ram_T_tot = False
                VSCF_Ram_T_comp = False
                VCI_IR = False
                VCI_Ram_0K_tot = False
                VCI_Ram_0K_comp = False
                VCI_Ram_T_tot = False
                VCI_Ram_T_comp = False

            if save:
                if HO_IR:
                    bufferHO_IR.append(line)
                if HO_Ram_0K_tot:
                    bufferHO_Ram_0K_tot.append(line)
                if HO_Ram_T_tot:
                    bufferHO_Ram_T_tot.append(line)
                if HO_Ram_0K_comp:
                    bufferHO_Ram_0K_comp.append(line)
                if HO_Ram_T_comp:
                    bufferHO_Ram_T_comp.append(line)
                if VSCF_IR:
                    bufferVSCF_IR.append(line)
                if VSCF_Ram_0K_tot:
                    bufferVSCF_Ram_0K_tot.append(line)
                if VSCF_Ram_T_tot:
                    bufferVSCF_Ram_T_tot.append(line)
                if VSCF_Ram_0K_comp:
                    bufferVSCF_Ram_0K_comp.append(line)
                if VSCF_Ram_T_comp:
                    bufferVSCF_Ram_T_comp.append(line)
                if VCI_IR:
                    bufferVCI_IR.append(line)
                if VCI_Ram_0K_tot:
                    bufferVCI_Ram_0K_tot.append(line)
                if VCI_Ram_T_tot:
                    bufferVCI_Ram_T_tot.append(line)
                if VCI_Ram_0K_comp:
                    bufferVCI_Ram_0K_comp.append(line)
                if VCI_Ram_T_comp:
                    bufferVCI_Ram_T_comp.append(line)
        # <--

        # Save and parse VSCF data for IR spectrum
        n_VSCF_ir = len(bufferVSCF_IR)
        if n_VSCF_ir > 0:
            for i, line in enumerate(bufferVSCF_IR[5:n_VSCF_ir-1]):
                IR_VSCF.append(line.split()[3:6])
                for j in range(3):
                    IR_VSCF[i][j] = float(IR_VSCF[i][j])
            IR_VSCF = np.array(IR_VSCF)
            self.IR_VSCF_T = IR_VSCF[:, 0:3:2]
            self.IR_VSCF_0K = IR_VSCF[:, 0:2]

        # Save and parse VCI data for IR spectrum
        n_VCI_ir = len(bufferVCI_IR)
        if n_VCI_ir > 0:
            for i, line in enumerate(bufferVCI_IR[5:n_VCI_ir-1]):
                IR_VCI.append(line.split()[3:6])
                for j in range(3):
                    IR_VCI[i][j] = float(IR_VCI[i][j])
            IR_VCI = np.array(IR_VCI)
            self.IR_VCI_T = IR_VCI[:, 0:3:2]
            self.IR_VCI_0K = IR_VCI[:, 0:2]

        # Save and parse HO data for IR spectrum
        n_HO_ir = len(bufferHO_IR)
        if n_HO_ir > 0:
            for i, line in enumerate(bufferHO_IR[5:n_HO_ir-1]):
                IR_HO.append(line.split()[3:6])
                for j in range(3):
                    IR_HO[i][j] = float(IR_HO[i][j])
            IR_HO = np.array(IR_HO)
            self.IR_HO_T = IR_HO[:, 0:3:2]
            self.IR_HO_0K = IR_HO[:, 0:2]

        # Save and parse HO data for Raman spectrum (0K, tot)
        n_HO_Ram_0K_tot = len(bufferHO_Ram_0K_tot)
        if n_HO_Ram_0K_tot > 0:
            for i, line in enumerate(bufferHO_Ram_0K_tot[6:n_HO_Ram_0K_tot-1]):
                Ram_HO_0K_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_HO_0K_tot[i][j] = float(Ram_HO_0K_tot[i][j])
            Ram_HO_0K_tot = np.array(Ram_HO_0K_tot)
            self.Ram_HO_0K_tot = Ram_HO_0K_tot[:, 0:2]
            self.Ram_HO_0K_per = Ram_HO_0K_tot[:, 0:3:2]
            self.Ram_HO_0K_par = Ram_HO_0K_tot[:, 0:4:3]

        # Save and parse HO data for Raman spectrum (T, tot)
        n_HO_Ram_T_tot = len(bufferHO_Ram_T_tot)
        if n_HO_Ram_T_tot > 0:
            for i, line in enumerate(bufferHO_Ram_T_tot[6:n_HO_Ram_T_tot-1]):
                Ram_HO_T_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_HO_T_tot[i][j] = float(Ram_HO_T_tot[i][j])
            Ram_HO_T_tot = np.array(Ram_HO_T_tot)
            self.Ram_HO_T_tot = Ram_HO_T_tot[:, 0:2]
            self.Ram_HO_T_per = Ram_HO_T_tot[:, 0:3:2]
            self.Ram_HO_T_par = Ram_HO_T_tot[:, 0:4:3]

        # Save and parse HO data for Raman spectrum (0K, comp)
        n_HO_Ram_0K_comp = len(bufferHO_Ram_0K_comp)
        if n_HO_Ram_0K_comp > 0:
            for i, line in enumerate(bufferHO_Ram_0K_comp[6:n_HO_Ram_0K_comp-1]):
                Ram_HO_0K_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_HO_0K_comp[i][j] = float(Ram_HO_0K_comp[i][j])
            Ram_HO_0K_comp = np.array(Ram_HO_0K_comp)
            self.Ram_HO_0K_comp_xx = Ram_HO_0K_comp[:, 0:2]
            self.Ram_HO_0K_comp_xy = Ram_HO_0K_comp[:, 0:3:2]
            self.Ram_HO_0K_comp_xz = Ram_HO_0K_comp[:, 0:4:3]
            self.Ram_HO_0K_comp_yy = Ram_HO_0K_comp[:, 0:5:4]
            self.Ram_HO_0K_comp_yz = Ram_HO_0K_comp[:, 0:6:5]
            self.Ram_HO_0K_comp_zz = Ram_HO_0K_comp[:, 0:7:6]

        # Save and parse HO data for Raman spectrum (T, comp)
        n_HO_Ram_T_comp = len(bufferHO_Ram_T_comp)
        if n_HO_Ram_T_comp > 0:
            for i, line in enumerate(bufferHO_Ram_T_comp[6:n_HO_Ram_T_comp-1]):
                Ram_HO_T_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_HO_T_comp[i][j] = float(Ram_HO_T_comp[i][j])
            Ram_HO_T_comp = np.array(Ram_HO_T_comp)
            self.Ram_HO_T_comp_xx = Ram_HO_T_comp[:, 0:2]
            self.Ram_HO_T_comp_xy = Ram_HO_T_comp[:, 0:3:2]
            self.Ram_HO_T_comp_xz = Ram_HO_T_comp[:, 0:4:3]
            self.Ram_HO_T_comp_yy = Ram_HO_T_comp[:, 0:5:4]
            self.Ram_HO_T_comp_yz = Ram_HO_T_comp[:, 0:6:5]
            self.Ram_HO_T_comp_zz = Ram_HO_T_comp[:, 0:7:6]

        # Save and parse VSCF data for Raman spectrum (0K, tot)
        n_VSCF_Ram_0K_tot = len(bufferVSCF_Ram_0K_tot)
        if n_VSCF_Ram_0K_tot > 0:
            for i, line in enumerate(bufferVSCF_Ram_0K_tot[6:n_VSCF_Ram_0K_tot-1]):
                Ram_VSCF_0K_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VSCF_0K_tot[i][j] = float(Ram_VSCF_0K_tot[i][j])
            Ram_VSCF_0K_tot = np.array(Ram_VSCF_0K_tot)
            self.Ram_VSCF_0K_tot = Ram_VSCF_0K_tot[:, 0:2]
            self.Ram_VSCF_0K_per = Ram_VSCF_0K_tot[:, 0:3:2]
            self.Ram_VSCF_0K_par = Ram_VSCF_0K_tot[:, 0:4:3]

        # Save and parse VSCF data for Raman spectrum (T, tot)
        n_VSCF_Ram_T_tot = len(bufferVSCF_Ram_T_tot)
        if n_VSCF_Ram_T_tot > 0:
            for i, line in enumerate(bufferVSCF_Ram_T_tot[6:n_VSCF_Ram_T_tot-1]):
                Ram_VSCF_T_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VSCF_T_tot[i][j] = float(Ram_VSCF_T_tot[i][j])
            Ram_VSCF_T_tot = np.array(Ram_VSCF_T_tot)
            self.Ram_VSCF_T_tot = Ram_VSCF_T_tot[:, 0:2]
            self.Ram_VSCF_T_per = Ram_VSCF_T_tot[:, 0:3:2]
            self.Ram_VSCF_T_par = Ram_VSCF_T_tot[:, 0:4:3]

        # Save and parse VSCF data for Raman spectrum (0K, comp)
        n_VSCF_Ram_0K_comp = len(bufferVSCF_Ram_0K_comp)
        if n_VSCF_Ram_0K_comp > 0:
            for i, line in enumerate(bufferVSCF_Ram_0K_comp[6:n_VSCF_Ram_0K_comp-1]):
                Ram_VSCF_0K_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VSCF_0K_comp[i][j] = float(Ram_VSCF_0K_comp[i][j])
            Ram_VSCF_0K_comp = np.array(Ram_VSCF_0K_comp)
            self.Ram_VSCF_0K_comp_xx = Ram_VSCF_0K_comp[:, 0:2]
            self.Ram_VSCF_0K_comp_xy = Ram_VSCF_0K_comp[:, 0:3:2]
            self.Ram_VSCF_0K_comp_xz = Ram_VSCF_0K_comp[:, 0:4:3]
            self.Ram_VSCF_0K_comp_yy = Ram_VSCF_0K_comp[:, 0:5:4]
            self.Ram_VSCF_0K_comp_yz = Ram_VSCF_0K_comp[:, 0:6:5]
            self.Ram_VSCF_0K_comp_zz = Ram_VSCF_0K_comp[:, 0:7:6]

        # Save and parse VSCF data for Raman spectrum (T, comp)
        n_VSCF_Ram_T_comp = len(bufferVSCF_Ram_T_comp)
        if n_VSCF_Ram_T_comp > 0:
            for i, line in enumerate(bufferVSCF_Ram_T_comp[6:n_VSCF_Ram_T_comp-1]):
                Ram_VSCF_T_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VSCF_T_comp[i][j] = float(Ram_VSCF_T_comp[i][j])
            Ram_VSCF_T_comp = np.array(Ram_VSCF_T_comp)
            self.Ram_VSCF_T_comp_xx = Ram_VSCF_T_comp[:, 0:2]
            self.Ram_VSCF_T_comp_xy = Ram_VSCF_T_comp[:, 0:3:2]
            self.Ram_VSCF_T_comp_xz = Ram_VSCF_T_comp[:, 0:4:3]
            self.Ram_VSCF_T_comp_yy = Ram_VSCF_T_comp[:, 0:5:4]
            self.Ram_VSCF_T_comp_yz = Ram_VSCF_T_comp[:, 0:6:5]
            self.Ram_VSCF_T_comp_zz = Ram_VSCF_T_comp[:, 0:7:6]

        # Save and parse VCI data for Raman spectrum (0K, tot)
        n_VCI_Ram_0K_tot = len(bufferVCI_Ram_0K_tot)
        if n_VCI_Ram_0K_tot > 0:
            for i, line in enumerate(bufferVCI_Ram_0K_tot[6:n_VCI_Ram_0K_tot-1]):
                Ram_VCI_0K_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VCI_0K_tot[i][j] = float(Ram_VCI_0K_tot[i][j])
            Ram_VCI_0K_tot = np.array(Ram_VCI_0K_tot)
            self.Ram_VCI_0K_tot = Ram_VCI_0K_tot[:, 0:2]
            self.Ram_VCI_0K_per = Ram_VCI_0K_tot[:, 0:3:2]
            self.Ram_VCI_0K_par = Ram_VCI_0K_tot[:, 0:4:3]

        # Save and parse VCI data for Raman spectrum (T, tot)
        n_VCI_Ram_T_tot = len(bufferVCI_Ram_T_tot)
        if n_VCI_Ram_T_tot > 0:
            for i, line in enumerate(bufferVCI_Ram_T_tot[6:n_VCI_Ram_T_tot-1]):
                Ram_VCI_T_tot.append(line.split()[3:7])
                for j in range(4):
                    Ram_VCI_T_tot[i][j] = float(Ram_VCI_T_tot[i][j])
            Ram_VCI_T_tot = np.array(Ram_VCI_T_tot)
            self.Ram_VCI_T_tot = Ram_VCI_T_tot[:, 0:2]
            self.Ram_VCI_T_per = Ram_VCI_T_tot[:, 0:3:2]
            self.Ram_VCI_T_par = Ram_VCI_T_tot[:, 0:4:3]

        # Save and parse VCI data for Raman spectrum (0K, comp)
        n_VCI_Ram_0K_comp = len(bufferVCI_Ram_0K_comp)
        if n_VCI_Ram_0K_comp > 0:
            for i, line in enumerate(bufferVCI_Ram_0K_comp[6:n_VCI_Ram_0K_comp-1]):
                Ram_VCI_0K_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VCI_0K_comp[i][j] = float(Ram_VCI_0K_comp[i][j])
            Ram_VCI_0K_comp = np.array(Ram_VCI_0K_comp)
            self.Ram_VCI_0K_comp_xx = Ram_VCI_0K_comp[:, 0:2]
            self.Ram_VCI_0K_comp_xy = Ram_VCI_0K_comp[:, 0:3:2]
            self.Ram_VCI_0K_comp_xz = Ram_VCI_0K_comp[:, 0:4:3]
            self.Ram_VCI_0K_comp_yy = Ram_VCI_0K_comp[:, 0:5:4]
            self.Ram_VCI_0K_comp_yz = Ram_VCI_0K_comp[:, 0:6:5]
            self.Ram_VCI_0K_comp_zz = Ram_VCI_0K_comp[:, 0:7:6]

        # Save and parse VCI data for Raman spectrum (T, comp)
        n_VCI_Ram_T_comp = len(bufferVCI_Ram_T_comp)
        if n_VCI_Ram_T_comp > 0:
            for i, line in enumerate(bufferVCI_Ram_T_comp[6:n_VCI_Ram_T_comp-1]):
                Ram_VCI_T_comp.append(line.split()[3:10])
                for j in range(7):
                    Ram_VCI_T_comp[i][j] = float(Ram_VCI_T_comp[i][j])
            Ram_VCI_T_comp = np.array(Ram_VCI_T_comp)
            self.Ram_VCI_T_comp_xx = Ram_VCI_T_comp[:, 0:2]
            self.Ram_VCI_T_comp_xy = Ram_VCI_T_comp[:, 0:3:2]
            self.Ram_VCI_T_comp_xz = Ram_VCI_T_comp[:, 0:4:3]
            self.Ram_VCI_T_comp_yy = Ram_VCI_T_comp[:, 0:5:4]
            self.Ram_VCI_T_comp_yz = Ram_VCI_T_comp[:, 0:6:5]
            self.Ram_VCI_T_comp_zz = Ram_VCI_T_comp[:, 0:7:6]

        return self

    def get_elatensor(self, *thickness):
        """
        Extracts the elastic tensor from the data.

        Args:
            \*thickness (float): Effective thickness of low dimensional
                 materials, in :math:`\\AA`.
        Returns:
            self.tensor (numpy.ndarray): Symmetrized elastic tensor in Voigt
                notation. For 3D systems, 6\*6; for 2D, 3\*3; for 1D, 1\*1.
                The matrix always has 2 dimensions. Unit: GPa for 3D and 1D, 2D
                if effective thickness of materials are specified. Otherwise
                GPa.m for 2D and GPa.m:math:`^{2}` for 1D (might lead to very
                small numbers).
        """
        import re, warnings

        title = ''
        buffer = []
        countline = 0
        while countline < self.eoo:
            line = self.data[countline]
            if re.match(r'^\s*SYMMETRIZED ELASTIC CONSTANTS', line):
                title = line
                countline += 2
                while self.data[countline].strip() != '':
                    buffer.append(self.data[countline])
                    countline += 1
                break
            else:
                countline += 1
                continue

        # Build tensor
        if len(buffer) == 0:
            raise Exception('Elastic tensor not found. Check your output file.')
        else:
            dimen = len(buffer)
            self.tensor = np.zeros([dimen,dimen], dtype=float)

        for i in range(dimen):
            # Clean buffer and copy it in strtensor
            self.tensor[i, i:] = np.array(
                buffer[i].strip().split()[1:-1], dtype=float
            )
        buffer.clear()

        # Symmetrize tensor
        for i in range(dimen):
            for j in range(i,dimen):
                self.tensor[j][i] = self.tensor[i][j]

        # Unit conversion
        if re.search('gpa', title.lower()):
            pass
        elif re.search('hartree', title.lower()):
            if dimen == 1:
                length = units.angstrom_to_au(self.get_lattice(initial=False)[0, 0]) * 1e-10 # AA --> m
                self.tensor =  units.au_to_GPa(self.tensor) * (units.au_to_angstrom(1.)*1e-10)**3 / length # Eh --> GJ/m
            elif dimen == 3:
                area = np.linalg.norm(np.cross(
                    self.get_lattice(initial=False)[0, :2],
                    self.get_lattice(initial=False)[1, :2]
                )) * 1e-20 # AA^2 --> m^2
                self.tensor = units.au_to_GPa(self.tensor) * (units.au_to_angstrom(1.)*1e-10)**3 / area # Eh --> GJ/m
        else:
            warnings.warn('Unknown unit identified, return to CRYSTAL default unit.',
                          stacklevel=2)

        # Effective thickness
        if len(thickness) > 0:
            if dimen == 1:
                if len(thickness) == 1:
                    self.tensor = self.tensor / (thickness[0]*1e-10)**2
                else:
                    self.tensor = self.tensor / (thickness[0]*1e-10) / (thickness[1]*1e-10)
            elif dimen == 3:
                self.tensor = self.tensor / (thickness[0]*1e-10)
        else:
            if dimen != 6:
                warngings.warn('Low dimensional materials without effective thickness! Output units have extra dimensions of length.',
                               stacklevel=2)
        return self.tensor


class Properties_input(Properties_inputBASE):
    """
    Properties input object inherited from the :ref:`Properties_inputBASE <ref-base-propd3>`
    object. For the basic set-ups of keywords, please refer to manuals :ref:`there <ref-base-propd3>`.

    Args:
        source (str): Template file or string. An empty object is generated for
            ``''``.
    """
    def __init__(self, source=''):
        import os

        super().__init__()
        if type(source) != str:
            raise ValueError('Unknown input value. Only string is allowed')

        if source != '':
            if os.path.exists(source):
                inp = open(source, 'r')
                source = inp.read()
                inp.close()
            super().analyze_text(source)

    def read_file(self, source):
        """
        Initialize the object if not yet done.

        Args:
            source (str): The name of the input file.
        Returns:
            self (Properties_input)
        """
        self.__init__(source)

    def write_file(self, file):
        """
        Writes the properties input to a file.

        Args:
            file (str): The name of the d3 file.
        """
        out = open(file, 'w')
        out.write('%s' % self.data)
        out.close()
        return self

    def make_bands_block(self, k_path, n_kpoints, first_band, last_band,
                         print_eig=0, print_option=1, precision=5,
                         title='BAND STRUCTURE CALCULATION'):
        """
        A shortcut for ``self.band()``. Also accepts pymatgen `HighSymmKpath <https://pymatgen.org/pymatgen.symmetry.html#pymatgen.symmetry.bandstructure.HighSymmKpath>`_
        objects for k path.

        Args:
            k_path (list | HighSymmKpath): The k-path for the bands calculation.
            n_kpoints (int): The number of k-points along the path.
            first_band (int): The index of the first band.
            last_band (int): The index of the last band.
            print_eig (int): Printing options for eigenvalues (default is 0).
            print_option (int): Properties printing options (default is 1).
            precision (int): Precision in the calculation of ``numpy.gcd``
            title (str): The title of the calculation.
        """
        import numpy as np

        bands_block = []

        # path from a pymatgen k_path object
        if 'HighSymmKpath' in str(type(k_path)):
            k_path_flat = [item for sublist in k_path.kpath['path']
                           for item in sublist]
            k_path_pmg = []
            for i in k_path_flat:
                # This is a pmg HighSymmKpath object
                k_path_pmg.append(k_path.kpath['kpoints'][i].tolist())
            k_path = np.array(k_path_pmg)
        elif type(k_path[0]) == list:
            # This is a list of lists
            k_path = np.array(k_path)
        else:
            raise ValueError('k_path must be a list of list (k coordinates) or a pymatgen HighSymmKpath object. {} used'.format(type(k_path)))

        k_unique = np.unique(k_path)
        # Find the shrinking factor
        k_unique = np.array(np.around(k_unique, precision)*10**precision, dtype=int)
        if len(k_unique) > 2:
            gcd = np.gcd.reduce(k_unique)
        else:
            gcd = np.gcd(k_unique[0], k_unique[1])

        k_path = (k_path/gcd)*10**precision
        k_path = np.array(k_path, dtype=int)
        shrink = int(10**precision/gcd)

        # Add the symmetry lines
        k_path_prt = []
        for i in range(len(k_path)-1):
            # The same K points repeat multiple times for HighSymmKpath
            if np.linalg.norm(k_path[i] - k_path[i+1]) < 1e-4:
                continue
            k_path_prt.append([k_path[i], k_path[i+1]])

        return self.band(title, len(k_path_prt)-1, shrink, n_kpoints, first_band,
                         last_band, print_option, print_eig, k_path_prt)

    def make_doss_block(self, n_points=200, band_range=None, e_range=None,
                        plotting_option=2, poly=12, print_option=1,
                        projections=[], proj_type='atom', output_file=None):
        """
        A shortcut for ``self.doss()``.

        Args:
            n_points (int): The number of points in the DOS plot.
            band_range (list or tuple): The range of bands to include in the DOS calculation.
            e_range (list or tuple): The energy range for the DOS calculation. Unit: eV
            plotting_option (int): DOS plotting options.
            poly (int): Degree of the polynomial for smearing.
            print_option (int): Properties printing options.
            projections (list): nProj\*nElement list of ints specifying the
                projections for the pdoss calculation, or nProj\*1 list of
                strings specifying element names to project on all atoms of the
                same element (``output_file`` is needed).
            proj_type (str): Type of projection ('atom' or 'site').
            output_file (str): Output file of 'crystal' calculation.
        """
        # either band_range or e_range needs to be specified
        if np.all(band_range==None) and np.all(e_range==None):
            raise ValueError('Either band_range or e_range should be specified. None specified.')
        elif np.all(band_range!=None) and np.all(e_range!=None):
            raise ValueError('Either band_range or e_range should be specified. 2 specified')
        elif type(band_range) == list and len(band_range) == 2:
            pass
        elif type(e_range) == list and len(e_range) == 2:
            e_range = [units.eV_to_H(i) for i in e_range]
        else:
            raise ValueError('Format error for band_range or the e_range (2 item list)')

        prj = []
        if projections != []:
            if proj_type != 'atom' and proj_type != 'ao':
                raise ValueError("proj_type should be either be 'ao' or 'atom'")

            if type(projections[0]) == str:
                if np.all(output_file==None):
                    raise ValueError('Outut file is needed.')
                else:
                    if proj_type == 'ao':
                        warnings.warn("When element name is specified, proj_type must be 'atom'.",
                                      stacklevel=2)

                    output = Crystal_output(output_file)
                    for prje in projections:
                        index = [i+1 for i, ele in enumerate(output.atom_symbols) if prje.upper() == ele.upper()]
                        prj.append([-len(index),] + index)
            else:
                if proj_type == 'atom':
                    for p in projections:
                        prj.append([-len(p),] + list(p))
                else:
                    for p in projections:
                        prj.append([len(p),] + list(p))

        if np.all(e_range==None):
            return self.doss(len(prj), n_points, band_range[0], band_range[1], plotting_option, poly, print_option, prj)
        else:
            return self.doss(len(prj), n_points, -1, -1, plotting_option, poly, print_option, e_range, prj)


    def make_pdoss_block(self, projections, proj_type='atom', output_file=None,
                         n_points=200, band_range=None, e_range=None,
                         plotting_option=2, poly=12, print_option=1):
        """
        Deprecated. Use either ``self.doss()`` or ``self.make_doss_block()``
        """
        import warnings

        warnings.warn('Deprecated. Use make_doss_block() method instead.', stacklevel=2)
        return self.make_doss_block(n_points, band_range, e_range, plotting_option,
                                    poly, print_option, projections, proj_type, output_file)

    def make_newk_block(self, shrink1, shrink2, Fermi=1, print_option=0):
        """
        Deprecated. Use ``self.newk()``.
        """
        import warnings

        warnings.warn('Deprecated. Use newk() method instead.', stacklevel=2)
        return self.newk(shrink1, shrink2, Fermi, print_option)

    def from_file(self, input_name):
        """
        Deprecated. Use ``read_file()``
        """
        import warnings

        warnings.warn('Deprecated. Use read_file() method instead.', stacklevel=2)
        return self.__init__(input_name)

    def write_properties_input(self, input_name):
        """
        Deprecated. Use ``write_file()``
        """
        import warnings

        warnings.warn('Deprecated. Use write_file() method instead.', stacklevel=2)
        return self.write_file(input_name)


class Properties_output(POutBASE):
    """
    Creates a Properties_output object. Since not all results from 'properties'
    executable creates are summarized in output, methods of this class usually
    requires multiple files.

    .. note::

        Most of methods directly dealing with output file itself are saved in
        ``base.output.POutBASE`` object. Users typically do not need to call
        methods there since they are called internally.

    Args:
        properties_output (str): Properties output file
    """

    def __init__(self, properties_output=None):
        if np.all(properties_output!=None):
            self.read_file(properties_output)
        else:
            pass

    def read_file(self, properties_output):
        """
        Parse the properties output file.

        Args:
            properties_output (str): The properties output file.
        Returns:
            Properties_output: The updated Properties_output object.
        """
        import os

        self.file_name = properties_output

        try:
            file = open(self.file_name, 'r')
            self.data = file.readlines()
            file.close()

            # directory
            dir_name = os.path.split(properties_output)[0]
            self.abspath = os.path.join(dir_name)

            # title (named "title" only to distinguish from "file_name" which means another thing)
            self.title = os.path.split(properties_output)[1]

        except:
            raise FileNotFoundError('EXITING: a CRYSTAL properties file needs to be specified')

    def read_vecfield(self, properties_output, which_prop):
        """Reads the fort.25 file to return data arrays containing one or more vectiorial density properties.

        Args:
            properties_output (str): The properties output file.
            which_prop (str): The density property selected by the user.
            'm' (magnetization), 'j' (spin current), 'J' (spin current density)
        Returns:
            Properties_output (str): The fort.25 output file.
        """

        import numpy as np

        self.read_file(properties_output)

        data = self.data

        # Reads the header information
        nrow = int(data[0].split()[1])
        ncol = int(data[0].split()[2])
        stepx = float(data[0].split()[3])
        stepy = float(data[0].split()[4])
        cosxy = float(data[0].split()[5])

        A = np.array([float(data[1].split()[0]), float(
            data[1].split()[1]), float(data[1].split()[2])])
        B = np.array([float(data[1].split()[3]), float(
            data[1].split()[4]), float(data[1].split()[5])])

        C = np.array([float(data[2].split()[0]), float(
            data[2].split()[1]), float(data[2].split()[2])])
        naf = int(data[2].split()[3])
        ldim = int(data[2].split()[4])

        self.header = (nrow, ncol, stepx, stepy, cosxy, A, B, C, naf, ldim)

        # Elaborates the header data
        skip = 6 + naf

        for i in range(2, 20):
            if (nrow % i) == 0:
                nrow_split = int(nrow/i)

        for i in range(2, 20):
            if (ncol % i) == 0:
                ncol_split = int(ncol/i)

        bline = (nrow*ncol)/6
        if (bline % 6) == 0:
            bline = int(bline)
        else:
            bline = int(bline) + 1

        # Reads the types of density properties requested by the user and initializes the data arrays
        check = np.zeros(3, dtype=int)
        if 'm' in which_prop:
            check[0] = 1
            self.dens_m = np.zeros((nrow, ncol, 3), dtype=float)
        if 'j' in which_prop:
            check[1] = 1
            self.dens_j = np.zeros((nrow, ncol, 3), dtype=float)
        if 'J' in which_prop:
            check[2] = 1
            self.dens_JX = np.zeros((nrow, ncol, 3), dtype=float)
            self.dens_JY = np.zeros((nrow, ncol, 3), dtype=float)
            self.dens_JZ = np.zeros((nrow, ncol, 3), dtype=float)
        if (not check[0]) and (not check[1]) and (not check[2]):
            print('Error: Invalid Entry. Only the m, j, and J characters are supported')
            sys.exit(1)

        # Gathers the data
        iamhere = 0

        if check[0]:
            iamhere = 3
            r = 0
            s = 0
            for i in range(0, bline):
                for j in range(0, len(data[i+iamhere].split())):
                    self.dens_m[r, s, 0] = data[i+iamhere].split()[j]
                    self.dens_m[r, s, 1] = data[i +
                                                iamhere+bline+skip].split()[j]
                    self.dens_m[r, s, 2] = data[i+iamhere +
                                                (2*bline)+(2*skip)].split()[j]
                    if s == (ncol - 1):
                        r += 1
                        s = 0
                    else:
                        s += 1
            iamhere = iamhere + 3*bline + 2*skip
        if check[1]:
            if iamhere == 0:
                iamhere = 3
            else:
                iamhere = iamhere + skip
            r = 0
            s = 0
            for i in range(0, bline):
                for j in range(0, len(data[i+iamhere].split())):
                    self.dens_j[r, s, 0] = data[i+iamhere].split()[j]
                    self.dens_j[r, s, 1] = data[i +
                                                iamhere+bline+skip].split()[j]
                    self.dens_j[r, s, 2] = data[i +
                                                iamhere+2*bline+2*skip].split()[j]
                    if s == (ncol - 1):
                        r += 1
                        s = 0
                    else:
                        s += 1
            iamhere = iamhere + 3*bline + 2*skip
        if check[2]:
            if iamhere == 0:
                iamhere = 3
            else:
                iamhere = iamhere + skip
            r = 0
            s = 0
            for i in range(0, bline):
                for j in range(0, len(data[i+iamhere].split())):
                    self.dens_JX[r, s, 0] = data[i+iamhere].split()[j]
                    self.dens_JX[r, s, 1] = data[i +
                                                 iamhere+bline+skip].split()[j]
                    self.dens_JX[r, s, 2] = data[i+iamhere +
                                                 (2*bline)+(2*skip)].split()[j]
                    self.dens_JY[r, s, 0] = data[i+iamhere +
                                                 (3*bline)+(3*skip)].split()[j]
                    self.dens_JY[r, s, 1] = data[i+iamhere +
                                                 (4*bline)+(4*skip)].split()[j]
                    self.dens_JY[r, s, 2] = data[i+iamhere +
                                                 (5*bline)+(5*skip)].split()[j]
                    self.dens_JZ[r, s, 0] = data[i+iamhere +
                                                 (6*bline)+(6*skip)].split()[j]
                    self.dens_JZ[r, s, 1] = data[i+iamhere +
                                                 (7*bline)+(7*skip)].split()[j]
                    self.dens_JZ[r, s, 2] = data[i+iamhere +
                                                 (8*bline)+(8*skip)].split()[j]
                    if s == (ncol - 1):
                        r += 1
                        s = 0
                    else:
                        s += 1
        return self

    def read_electron_band(self, band_file):
        """
        Generate bands object from CRYSTAL BAND.DAT or fort.25 file. Energy
        unit: eV. E Fermi is aligned to 0.

        Args:
            band_file (str): Name of BAND.DAT or fort.25 file

        Returns:
            self.bands (ElectronBand): An ``CRYSTALpytools.electronics.ElectronBand`` object
        """
        from CRYSTALpytools.electronics import ElectronBand
        import warnings

        if not hasattr(self, 'file_name'):
            warnings.warn('Properties output file not found: 3D k path not available',
                          stacklevel=2)
            self.bands = ElectronBand.from_file(band_file)
        else:
            self.bands = ElectronBand.from_file(band_file, self.file_name)

        return self.bands

    def read_electron_dos(self, dos_file):
        """
        Generate doss object from CRYSTAL DOSS.DAT or fort.25 file. Energy
        unit: eV. E Fermi is aligned to 0.
        Args:
            dos_file (str): Name of DOSS.DAT or fort.25 file

        Returns:
            self.doss (ElectronDOS): An ``CRYSTALpytools.electronics.ElectronDOS`` object
        """
        from CRYSTALpytools.electronics import ElectronDOS

        self.doss = ElectronDOS.from_file(dos_file)

        return self.doss

    def read_cry_contour(self, properties_output):
        """Read the CRYSTAL contour files to create the contour objects.

        Args:
            properties_output (str): The properties output file.
        Returns:
            Properties_output: The updated Properties_output object.
        """
        import re
        import sys

        import numpy as np
        import pandas as pd

        self.read_file(properties_output)

        filename = str(properties_output)

        tipo = ''

        if (filename.endswith('.SURFRHOO')):
            self.tipo = 'SURFRHOO'
            self.path = filename
        elif (filename.endswith('.SURFLAPP')):
            self.tipo = 'SURFLAPP'
            self.path = filename
        elif (filename.endswith('.SURFLAPM')):
            self.tipo = 'SURFLAPM'
            self.path = filename
        elif (filename.endswith('.SURFGRHO')):
            self.tipo = 'SURFGRHO'
            self.path = filename
        elif (filename.endswith('.SURFELFB')):
            self.tipo = 'SURFELFB'
            self.path = filename
        elif (filename.endswith('.SURFVIRI')):
            self.tipo = 'SURFVIRI'
            self.path = filename
        elif (filename.endswith('.SURFGKIN')):
            self.tipo = 'SURFGKIN'
            self.path = filename
        elif (filename.endswith('.SURFKKIN')):
            self.tipo = 'SURFKKIN'
            self.path = filename
        else:
            sys.exit('Please choose a valid file')

        l_dens = self.data

        n_punti_x = int(l_dens[1].strip().split()[0])
        n_punti_y = int(l_dens[1].strip().split()[1])

        self.npx = n_punti_x

        x_min = units.au_to_angstrom(float(l_dens[2].strip().split()[0]))
        x_max = units.au_to_angstrom(float(l_dens[2].strip().split()[1]))
        x_step = units.au_to_angstrom(float(l_dens[2].strip().split()[2]))

        y_min = units.au_to_angstrom(float(l_dens[3].strip().split()[0]))
        y_max = units.au_to_angstrom(float(l_dens[3].strip().split()[1]))
        y_step = units.au_to_angstrom(float(l_dens[3].strip().split()[2]))

        l_dens = l_dens[5:]

        m_dens = []
        for i in l_dens:
            m_dens.append(re.sub("\s\s+", " ", i))

        n_dens = []
        for i in m_dens:
            n_dens.append(i.replace('\n', '').split())

        self.df = pd.DataFrame(n_dens)

        self.x_points = np.linspace(x_min, x_max, n_punti_x)
        self.y_points = np.linspace(y_min, y_max, n_punti_y)

        a = x_max - x_min
        b = y_max - y_min
        r = a/b

        self.x_graph_param = 10
        self.y_graph_param = 10 / r

        ctr1 = np.array([0.002, 0.004, 0.008, 0.02, 0.04,
                         0.08, 0.2, 0.4, 0.8, 2, 4, 8, 20])
        colors1 = ['r', 'r', 'r', 'r', 'r', 'r',
                   'r', 'r', 'r', 'r', 'r', 'r', 'r']
        ls1 = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        ctr2 = np.array([-8, -4, -2, -0.8, -0.4, -0.2, -0.08, -0.04, -0.02, -0.008, -0.004, -0.002, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                         0.2, 0.4, 0.8, 2, 4, 8])
        colors2 = ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b',
                   'b', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
        ls2 = ['--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--',
               '--', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        ctr3 = np.array([0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
                         0.95, 1])
        colors3 = ['k', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b',
                   'b', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
        ls3 = ['dotted', '--', '--', '--', '--', '--', '--', '--', '--',
               '--', '--', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        if (self.tipo == 'SURFRHOO') or (self.tipo == 'SURFGRHO') or (self.tipo == 'SURFGKIN'):
            self.levels = ctr1
            self.colors = colors1
            self.linestyles = ls1
            self.fmt = '%1.3f'
        elif (self.tipo == 'SURFLAPP') or (self.tipo == 'SURFLAPM') or (self.tipo == 'SURFVIRI') or (self.tipo == 'SURFKKIN'):
            self.levels = ctr2
            self.colors = colors2
            self.linestyles = ls2
            self.fmt = '%1.3f'
        elif (self.tipo == 'SURFELFB'):
            self.levels = ctr3
            self.colors = colors3
            self.linestyles = ls3
            self.fmt = '%1.2f'

        return self

    def read_cry_xrd_spec(self, properties_output):
        """
        Read XRD spectrum data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted XRD spectrum data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        if filename.endswith('.outp'):
            pass
        else:
            sys.exit('please, choose a valid file or rename it properly')

        spectrum = re.compile(
            '2THETA    INTENS  INTENS-LP INTENS-LP-DW', re.DOTALL)

        match = []

        a = 0

        for line in data:
            if spectrum.search(line):
                match.append('WRITE LINE:' + line)
                a = 1
            else:
                match.append('WRONG LINE:' + line)

        if (a == 0):
            sys.exit('please, choose a valid file or rename it properly')

        df = pd.DataFrame(match)

        num_riga = (df[df[0].str.contains(u'WRITE')].index.values)

        num_riga = num_riga[0]

        match = match[num_riga:]

        pattern = re.compile(
            '\s+ \d+\.\d+ \s+ \d+\.\d+ \s+ \d+\.\d+ \s+ \d+\.\d+\n', re.DOTALL)

        match_2 = []

        for line in match:
            if pattern.search(line):
                # pulisco dalle scritte di prima
                line = line.replace('WRONG LINE:', '')
                match_2.append(line)

        df = pd.DataFrame([i.strip().split() for i in match_2])

        for i in range(0, 4):
            df[i] = df[i].astype(float)

        df = df.rename(columns={0: '2THETA', 1: 'INTENS',
                                2: 'INTENS-LP', 3: 'INTENS-LP-DW'})

        self.x = df['2THETA']
        self.y = df['INTENS-LP']

        self.title = title[:-1]

        return self

    def read_cry_rholine(self, properties_output):
        """
        Read density line data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted density line data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        l_dens = self.data
        filename = self.abspath
        title = self.title

        if filename.endswith('.RHOLINE'):
            pass
        else:
            sys.exit('please, choose a valid file or rename it properly')

        m_dens = []
        for i in l_dens:
            m_dens.append(re.sub("\s\s+", " ", i))

        n_dens = []
        for i in m_dens:
            n_dens.append(i.replace('\n', '').split())

        df_dens = pd.DataFrame(n_dens)
        df_dens = df_dens.dropna()

        for i in range(0, len(df_dens.columns)):
            df_dens[i] = pd.to_numeric(df_dens[i])

        self.x = units.au_to_angstrom((df_dens[0]-5.55))
        self.y = df_dens[1]/0.148184743

        self.title = title[:-4]

        return self

    def read_cry_seebeck(self, properties_output):
        """
        Read Seebeck coefficient data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted Seebeck coefficient data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        spectrum = re.compile('Npoints', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE:' + line)
            else:
                match.append('WRONG LINE:' + line)

        df = pd.DataFrame(match)
        indx = list(df[df[0].str.contains("RIGHT")].index)

        lin = []
        for i in indx:
            lin.append(i+1)

        diffs = [abs(x - y) for x, y in zip(lin, lin[1:])]

        length = diffs[0] - 1

        lif = []
        for i in lin:
            lif.append(i+length)

        c = []
        for i in range(len(lin)):
            c.append(lin[i])
            c.append(lif[i])

        d = [c[i:i + 2] for i in range(0, len(c), 2)]

        l = []
        for i in range(0, len(d)):
            pd.DataFrame(l.append(df[d[i][0]:d[i][1]]))

        right = df[df[0].str.contains("RIGHT")]
        right = right.reset_index().drop('index', axis=1)

        self.temp = []

        for i in range(0, len(right)):
            self.temp.append(float(str(right[0][i])[20:24]))

        ll = []
        for k in range(0, len(l)):
            ll.append(l[k].reset_index().drop('index', axis=1))

        self.all_data = []
        for k in range(0, len(ll)):
            for i in ll[k]:
                self.all_data.append(ll[k][i].apply(
                    lambda x: x.replace('WRONG LINE:', '')))

        self.volume = (float(str(match[2:3])[-13:-4]))

        self.title = title

        return self

    def read_cry_sigma(self, properties_output):
        """
        Read electrical conductivity data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted electrical conductivity data.
        """
        import re
        import sys

        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        spectrum = re.compile('Npoints', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE:' + line)
            else:
                match.append('WRONG LINE:' + line)

        df = pd.DataFrame(match)
        indx = list(df[df[0].str.contains("RIGHT")].index)

        lin = []
        for i in indx:
            lin.append(i+1)

        diffs = [abs(x - y) for x, y in zip(lin, lin[1:])]

        length = diffs[0] - 1

        lif = []
        for i in lin:
            lif.append(i+length)

        c = []
        for i in range(len(lin)):
            c.append(lin[i])
            c.append(lif[i])

        d = [c[i:i + 2] for i in range(0, len(c), 2)]

        l = []
        for i in range(0, len(d)):
            pd.DataFrame(l.append(df[d[i][0]:d[i][1]]))

        right = df[df[0].str.contains("RIGHT")]
        right = right.reset_index().drop('index', axis=1)

        self.temp = []

        for i in range(0, len(right)):
            self.temp.append(float(str(right[0][i])[20:24]))

        ll = []
        for k in range(0, len(l)):
            ll.append(l[k].reset_index().drop('index', axis=1))

        self.all_data = []
        for k in range(0, len(ll)):
            for i in ll[k]:
                self.all_data.append(ll[k][i].apply(
                    lambda x: x.replace('WRONG LINE:', '')))

        self.volume = (float(str(match[2:3])[-13:-4]))

        self.title = title

        return self

    def read_cry_lapl_profile(self, properties_output):
        """
        Read Laplacian profile data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted Laplacian profile data.
        """
        import re

        import numpy as np
        import pandas as pd

        data = self.data
        filename = self.abspath
        title = self.title

        self.read_file(properties_output)

        spectrum = re.compile('PROFILE ALONG THE POINTS', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE: ' + line)
            else:
                match.append('WRONG LINE: ' + line)

        df = pd.DataFrame(match)
        num_riga = (df[df[0].str.contains(u'RIGHT')].index.values)
        num_in = num_riga + 8

        spectrum_fin = re.compile('EEEEEEEEEE TERMINATION  DATE', re.DOTALL)

        match_fin = []

        for line in data:
            if spectrum_fin.search(line):
                match_fin.append('RIGHT LINE: ' + line)
            else:
                match_fin.append('WRONG LINE: ' + line)

        df_fin = pd.DataFrame(match_fin)

        num_fin = (df_fin[df_fin[0].str.contains(u'RIGHT')].index.values)
        match = match[num_in[0]:num_fin[0]-2]

        match_2 = []
        for line in match:
            line = line.replace('WRONG LINE:', '')
            match_2.append(line)

        df = pd.DataFrame([i.strip().split() for i in match_2])
        df = df.rename({0: 'x', 1: 'y', 2: 'z', 3: 'dist',
                        4: 'rho', 5: 'lapl', 6: 'L3', 7: 'ellip'}, axis=1)
        df = df.drop(df[df.lapl == '********'].index)
        df = df.reset_index().drop('index', axis=1)
        df['lapl'] = df['lapl'].apply(pd.to_numeric)
        df['dist'] = df['dist'].apply(pd.to_numeric)
        self.datax = df.dist
        self.datay = df.lapl

        return self

    def read_cry_density_profile(self, properties_output):
        """
        Read density profile data from a file.

        Args:
            properties_output (str): Path to the properties output file.
        Returns:
            self: The modified object with extracted density profile data.
        """
        import re

        import numpy as np
        import pandas as pd

        self.read_file(properties_output)

        data = self.data
        filename = self.abspath
        title = self.title

        spectrum = re.compile('PROFILE ALONG THE POINTS', re.DOTALL)

        match = []

        for line in data:
            if spectrum.search(line):
                match.append('RIGHT LINE: ' + line)
            else:
                match.append('WRONG LINE: ' + line)

        df = pd.DataFrame(match)
        num_riga = (df[df[0].str.contains(u'RIGHT')].index.values)
        num_in = num_riga + 8

        spectrum_fin = re.compile('EEEEEEEEEE TERMINATION  DATE', re.DOTALL)

        match_fin = []

        for line in data:
            if spectrum_fin.search(line):
                match_fin.append('RIGHT LINE: ' + line)
            else:
                match_fin.append('WRONG LINE: ' + line)

        df_fin = pd.DataFrame(match_fin)

        num_fin = (df_fin[df_fin[0].str.contains(u'RIGHT')].index.values)
        match = match[num_in[0]:num_fin[0]-2]

        match_2 = []
        for line in match:
            line = line.replace('WRONG LINE:', '')
            match_2.append(line)

        df = pd.DataFrame([i.strip().split() for i in match_2])
        df = df.rename({0: 'x', 1: 'y', 2: 'z', 3: 'dist',
                        4: 'rho', 5: 'lapl', 6: 'L3', 7: 'ellip'}, axis=1)
        df = df.drop(df[df.lapl == '********'].index)
        df = df.reset_index().drop('index', axis=1)
        df['rho'] = df['rho'].apply(pd.to_numeric)
        df['dist'] = df['dist'].apply(pd.to_numeric)
        self.datax = df.dist
        self.datay = df.rho

        return self

    def read_cry_ECHG(self, f25_file):
        """
        Read density profile data from a file.

        Args:
            f25_file (str): Path to the fort.25 file.
        Returns:
            self.echg (ChargeDensity): ``electronics.ChargeDensity`` object.
        """
        from CRYSTALpytools.electronics import ChargeDensity

        self.echg = ChargeDensity.read_ECHG(f25_file)
        return self.echg

    def read_cry_ECHG_delta(self, f25_file1, f25_file2):
        """
        Read density profile data from two files and plots the difference. It
        is important to have consistent grid definitions, otherwise error is
        raised.

        Args:
            f25_file1 (str): Path to first fort.25 file.
            f25_file2 (str): Path to second fort.25 file.
        Returns:
            self.echg (ChargeDensity): ``electronics.ChargeDensity`` object.
        """
        from CRYSTALpytools.electronics import ChargeDensity

        self.echg = ChargeDensity.read_ECHG(f25_file1, f25_file2)
        return self.echg

#------------------------------------------------------------------------------#
#--------------------------------obsolete methods------------------------------#
#------------------------------------------------------------------------------#
    def read_cry_band(self, band_file):
        """
        Obsolete. Use ``read_electron_band``.
        """
        import warnings

        warnings.warn("You are calling an obsolete function. Use 'read_electron_band' instead.",
                      stacklevel=2)
        return self.read_electron_band(band_file)

    def read_cry_doss(self, dos_file):
        """
        Obsolete. Use ``read_electron_dos``.
        """
        import warnings

        warnings.warn("You are calling an obsolete function. Use 'read_electron_dos' instead.",
                      stacklevel=2)
        return self.read_electron_dos(dos_file)



class Crystal_gui(Crystal_gui):
    """
    Inherited from ``geometry.Crystal_gui``. See documentations :ref:`there <ref-geometry>`.
    """


class Crystal_density():
    """
    Read density data from a .f98 file.
    """

    def __init__(self):

        pass

    def read_cry_irr_density(self, fort98_unit):
        """
        Read density profile data from a CRYSTAL .f98 file.
        Args:
            fort98_unit (str): The file containing the formatted density matrix.
        Returns:
            None
        Note:
        This is a work in progress. If you are interested in this functionality,
        please open an Issue on GitHub.
        """
        self.file_name = fort98_unit
        self.is_irr = True

        import re
        import sys

        import numpy as np

        try:
            file = open(self.file_name, 'r')
            data = file.readlines()
            file.close()
        except:
            print('EXITING: a CRYSTAL .f98 file needs to be specified')
            sys.exit(1)

        self.all_file = data

        # the keyword BASATO appears twice, this is a check to see which basato
        # is being read
        basato1 = False

        for i, line in enumerate(data):

            if re.match(r'^LIMINF LIMTOL LIMPAR', line):
                inf_vec_len, tol_vec_len, par_vec_len = [
                    int(x) for x in data[i+1].split()]

            elif re.match(r'^INF', line):
                self.inf_vec = []
                inf_n_line = int(np.ceil(inf_vec_len/8))
                for j in range(inf_n_line):
                    self.inf_vec.extend([int(x) for x in data[i+1+j].split()])
                n_symmops = self.inf_vec[0]
                n_atoms = self.inf_vec[23]
                n_shells = self.inf_vec[19]
                '''if self.inf_vec[26] == 0:
                    self.spin_pol == False
                elif self.inf_vec[26] == 1:
                    self.spin_pol == True'''
                n_prim_gto = self.inf_vec[74]
                f_irr_len = (self.inf_vec[63]+1)*self.inf_vec[227]
                p_irr_len = (self.inf_vec[63]+1)*self.inf_vec[18]
                nnnc_len = self.inf_vec[190]
                la3_len = self.inf_vec[55]
                # n_symmops_noinv = inf_vec[1]

            elif re.match(r'^TOL', line):
                self.tol_vec = []
                tol_n_line = int(np.ceil(tol_vec_len/8))
                for j in range(tol_n_line):
                    self.tol_vec.extend([int(x) for x in data[i+1+j].split()])

            elif re.match(r'^PAR', line):
                self.par_vec = []
                par_n_line = int(np.ceil(par_vec_len/4))
                for j in range(par_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+1+j])/20)):
                        self.par_vec.append(
                            float(data[i+1+j][(item)*20:(item+1)*20]))

            elif re.match(r'^XYVGVE', line):
                # This vector contains the rotations, translation,
                # lattice vectors and transformation matrix from primitive to
                # crystallographic cell
                # Read all of it first and separate later
                xyvgve_n_line = int(np.ceil((n_symmops*12+18)/4))
                xyvgve_vec = []
                for j in range(xyvgve_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+1+j])/20)):
                        xyvgve_vec.append(
                            float(data[i+1+j][(item)*20:(item+1)*20]))
                # Now let's split the xyvgve_vec
                self.rotations_vec = xyvgve_vec[0:n_symmops*9]
                self.translations_vec = xyvgve_vec[n_symmops *
                                                   9:n_symmops*9+n_symmops*3]
                self.direct_lattice_vec = xyvgve_vec[n_symmops *
                                                     12:n_symmops*12+9]
                self.transf_matrix = xyvgve_vec[-9:]

            elif re.match(r'^BASATO', line):
                if basato1 == False:
                    basato_n_line = int(
                        np.ceil((n_atoms*4+n_shells*5+n_prim_gto*7)/4))
                    basato_vec = []
                    for j in range(basato_n_line):
                        # The negative elements appear connected to the previous one
                        # eg:  0.0000000000000E+00-1.0000000000000E+00
                        # The line below fixes that issue
                        for item in range(0, int(len(data[i+1+j])/20)):
                            basato_vec.append(
                                float(data[i+1+j][(item)*20:(item+1)*20]))
                    # Extract the iformation we need from basato

                    # Atom coordinates
                    self.atom_coord = []
                    for j in range(0, 3*n_atoms, 3):
                        self.atom_coord.append(
                            basato_vec[(n_atoms+j):(n_atoms+j+3)])
                    # self.atom_coord = np.array(self.atom_coord)

                    # Assign the shell to the atom
                    self.shell_coord = []
                    # The beginning of the part of BASATO I need here
                    init = 4*n_atoms + 2*n_shells
                    for j in range(0, 3*n_shells, 3):
                        self.shell_coord.append(
                            basato_vec[(init+j):(init+j+3)])
                    # self.shell_coord = np.array(self.shell_coord)

                    # Array that defines which atom a shell belongs to
                    self.shell_to_atom = []
                    for coord in self.shell_coord:
                        self.shell_to_atom.append(self.atom_coord.index(coord))
                    basato1 = True

                elif basato1 == True:
                    self.basato2 = []
                    j = i + 1
                    while 'SPINOR' not in data[j].split()[0]:
                        # The negative elements appear connected to the previous one
                        # eg:  0.0000000000000E+00-1.0000000000000E+00
                        # As opposite to the loops above where the float read was 20
                        # characters long, this ones are 21
                        # The line below fixes that issue
                        self.basato2.extend([int(x) for x in data[j].split()])
                        j += 1
                    self.atom_shell = self.basato2[-n_shells:]

            elif re.match(r'^SPINOR', line):
                self.f_irr = []
                self.charges = []
                # self.spin = [0]*(2*n_atoms) #tmp
                self.spin = []
                self.ghost = []
                n_ghost = 0
                n_spin_line = int(np.ceil((n_atoms*2)/8))
                n_basold = 0
                if 'BASOLD' in data[i+n_spin_line+1]:
                    n_basold = 9 + 3 * \
                        self.inf_vec[1] + n_shells + 3*n_atoms+3 * \
                        n_shells+self.inf_vec[4]+1+3*self.inf_vec[78]
                n_basold_line = int(np.ceil((n_basold)/4))
                n_charge_line = int(np.ceil(n_atoms/4))
                skip = n_spin_line + n_charge_line + 1 + n_basold_line
                for j in range(n_spin_line):
                    self.spin.extend([int(x) for x in data[i + j + 1].split()])
                if 'IGHOST' in data[i+n_spin_line+1]:
                    n_ghost = int(np.ceil((n_atoms)/8)) + 1
                    skip = skip + n_ghost
                    for j in range(n_ghost-1):
                        self.ghost.extend(
                            [float(x) for x in data[i + j + n_spin_line + 2].split()])
                f_irr_n_line = int(np.ceil(f_irr_len/4))
                for j in range(n_charge_line):
                    self.charges.extend(
                        [float(x) for x in data[i+j+n_spin_line+n_basold+n_ghost+1].split()])
                for j in range(f_irr_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+skip+j])/21)):
                        self.f_irr.append(
                            float(data[i+skip+j][(item)*21:(item+1)*21]))
                self.p_irr = []
                p_irr_n_line = int(np.ceil(p_irr_len/4))
                skip += 1
                for k in range(i+skip+j, i+skip+j+p_irr_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    for item in range(0, int(len(data[k])/21)):
                        self.p_irr.append(
                            float(data[k][(item)*21:(item+1)*21]))

            elif re.match(r'^   NCF', line):
                # The ncf vector contains the pointers to the symmetry irerducible
                # shell couples la3, la4
                self.ncf = []
                j = i+1
                while 'NSTATG' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.ncf.extend([int(x) for x in data[j].split()])
                    j += 1

            elif re.match(r'^NSTATG', line):
                # The nstatg vector contains the pointers to the starting point
                # of each couple set in P_irr and F_irr
                self.nstatg = []
                j = i+1
                while 'NSTAFG' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.nstatg.extend([int(x) for x in data[j].split()])
                    j += 1

            elif re.match(r'^  NNNC', line):
                # The nnnc points the starting position in the P matrix for
                # each couple, and its size corresponds to the total number
                # of shell couple in the shell couple sets
                self.nnnc = []
                nnnc_n_line = int(np.ceil(nnnc_len/8))
                for j in range(nnnc_n_line):
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    for item in range(0, int(len(data[i+1+j])/10)):
                        self.nnnc.append(
                            int(data[i+1+j][(item)*10:(item+1)*10]))
            elif re.match(r'^   LA3', line):
                # The nnnc points the starting position in the P matrix for
                # each couple, and its size corresponds to the total number
                # of shell couple in the shell couple sets
                self.la3 = []
                # nnnc_n_line = int(np.ceil(nnnc_len/8))
                j = i+1
                while 'LA4' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.la3.extend([int(x) for x in data[j].split()])
                    j += 1
            elif re.match(r'^   LA4', line):
                # The nnnc points the starting position in the P matrix for
                # each couple, and its size corresponds to the total number
                # of shell couple in the shell couple sets
                self.la4 = []
                # nnnc_n_line = int(np.ceil(nnnc_len/8))
                j = i+1
                while 'IROF' not in data[j].split()[0]:
                    # The negative elements appear connected to the previous one
                    # eg:  0.0000000000000E+00-1.0000000000000E+00
                    # As opposite to the loops above where the float read was 20
                    # characters long, this ones are 21
                    # The line below fixes that issue
                    self.la4.extend([int(x) for x in data[j].split()])
                    j += 1


def cry_combine_density(density1, density2, density3, new_density='new_density.f98', spin_pol=False):
    """
    Combine density matrix files.

    Args:
        density1 (str): The first density matrix file.
        density2 (str): The second density matrix file.
        density3 (str): The density matrix file for the whole system.
        new_density (str, optional): The name of the new density matrix. Defaults to 'new_density.f98'.
        spin_pol (bool, optional): Specifies if the system is spin polarized. Defaults to False.
    Returns:
        None
    Note:
        This is a work in progress. If you are interested in this functionality,
        please open an Issue on GitHub.
    """

    import sys

    import numpy as np

    try:
        density1_data = Crystal_density(density1)  # substrate
        density2_data = Crystal_density(density2)  # adsorbate

        density3_data_obj = Crystal_density(density3)
        density3_data = density3_data_obj.all_file
    except:
        print('EXITING: a CRYSTAL .f98 file needs to be specified')
        sys.exit(1)

    # Find P_irr <-> atom correspondence
    fragment_1 = []
    fragment_2 = []

    for i, j in enumerate(density1_data.ncf):
        # density1_data.la3[j] is the shell number
        # density1_data.atom_shell[density1_data.la3[j]] is the atom position number (1-6)
        # density1_data.ghost[density1_data.atom_shell[density1_data.la3[j]]] is either 0 or atomic number depending on ghost or not
        # This tells me if the shell belongs to this fragment
        n_elements = density1_data.nstatg[i] - density1_data.nstatg[i - 1]

        if density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1] == 0 and \
                density1_data.ghost[density1_data.atom_shell[density1_data.la4[j-1]-1]-1] == 0:
            fragment_1.extend([True]*n_elements)
        else:
            fragment_1.extend([False]*n_elements)

        if density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1] != 0 and \
                density1_data.ghost[density1_data.atom_shell[density1_data.la4[j-1]-1]-1] != 0:
            fragment_2.extend([True]*n_elements)
        else:
            fragment_2.extend([False]*n_elements)

    if spin_pol == True:
        spin_p1 = fragment_1.copy()
        spin_p2 = fragment_2.copy()
        fragment_1.extend(spin_p1)
        fragment_2.extend(spin_p2)

    beginning = density3_data.index('SPINOR\n')
    end = density3_data.index('   NCF\n')
    sum_density = (np.array(density1_data.p_irr) +
                   np.array(density2_data.p_irr))/2
    sum_fock = np.array(density1_data.f_irr)+np.array(density2_data.f_irr)
    sum_charges = np.array(density1_data.charges) + \
        np.array(density2_data.charges)
    spinor = ['SPINOR\n']
    charges = []
    fock = []
    density = []
    new_fock = sum_fock  # TMP
    new_fock = [0] * len(density3_data_obj.f_irr)
    new_p = []

    for i in range(len(fragment_1)):
        if fragment_1[i] == True and fragment_2[i] == False:
            # new_fock.append(density1_data.f_irr[i])
            new_p.append(density1_data.p_irr[i])
        elif fragment_1[i] == False and fragment_2[i] == True:
            # new_fock.append(density2_data.f_irr[i])
            new_p.append(density2_data.p_irr[i])
        elif fragment_1[i] == False and fragment_2[i] == False:
            # new_fock.append(0.)
            new_p.append(sum_density[i])
            # new_p.append(0.)
            # new_p.append(density3_data_obj.p_irr[i])

    for i in range(0, len(density3_data_obj.spin), 8):
        spinor.append(' '.join([str(x)
                                for x in density3_data_obj.spin[i:i+8]])+'\n')
    for i in range(0, len(sum_charges), 4):
        charges.append(' '.join(["{:.13e}".format(x)
                                 for x in sum_charges[i:i+4]])+'\n')
    for i in range(0, len(new_fock), 4):
        fock.append(' '.join(["{:.13e}".format(x)
                              for x in new_fock[i:i+4]])+'\n')
    for i in range(0, len(new_p), 4):
        density.append(' '.join(["{:.13e}".format(x)
                                 for x in new_p[i:i+4]])+'\n')

    final_fort98 = density3_data[0:beginning] + \
        spinor+charges+fock+density+density3_data[end:]
    with open(new_density, 'w') as file:
        for line in final_fort98:
            file.writelines(line)


def write_cry_density(fort98_name, new_p, new_fort98):
    """
    Write the formatted density matrix.

    Args:
        fort98_name (str): The name of the previous density matrix file.
        new_p (list): The new density matrix.
        new_fort98 (str): The name of the new density matrix file.

        Returns:
        None
    Note:
        This is a work in progress. If you are interested in this functionality,
        please open an Issue on GitHub.
    """

    import numpy as np

    file = open(fort98_name, 'r')
    data = file.readlines()
    file.close()

    density = Crystal_density(fort98_name)

    n_spin_line = int(np.ceil((density.inf_vec[23] * 2) / 8))
    n_charges_line = int(np.ceil((density.inf_vec[23]) / 4))
    beginning = data.index('SPINOR\n') + n_spin_line + n_charges_line + 1
    end = data.index('   NCF\n')

    new_fock_vect = [0] * len(density.f_irr)

    new_fock = []
    for i in range(0, len(new_fock_vect), 4):
        new_fock.append(' '.join(["{:.13e}".format(x)
                                  for x in new_fock_vect[i:i + 4]]) + '\n')

    new_density = []
    for i in range(0, len(new_p), 4):
        new_density.append(' '.join(["{:.13e}".format(x)
                                     for x in new_p[i:i+4]])+'\n')

    final_fort98 = data[0:beginning]+new_fock+new_density+data[end:]
    with open(new_fort98, 'w') as file:
        for line in final_fort98:
            file.writelines(line)


class External_unit:
    # WORK IN PROGRESS
    # This class will generate an object from the CRYSTAL external units like: IRSPEC.DAT
    # RAMSPEC.DAT tha can be plotted through the corresponding function in plot.py

    def __init__(self):

        pass

    def read_external_unit(self, external_unit):
        import os

        self.file_name = external_unit
        try:
            file = open(external_unit, 'r')
            self.data = file.readlines()
            file.close()

            # directory
            dir_name = os.path.split(external_unit)[0]
            self.abspath = os.path.join(dir_name)

            # title (named "title" only to distinguish from "file_name" which means another thing)
            self.title = os.path.split(external_unit)[1]

        except:
            raise FileNotFoundError(
                'EXITING: a CRYSTAL generated .DAT unit file needs to be specified')

    def read_cry_irspec(self, external_unit):
        import numpy as np

        self.read_external_unit(external_unit)

        data = self.data

        for index, line in enumerate(data):
            data[index] = line.split()

        columns = len(data[0])
        no_points = len(data)

        if columns > 3:
            self.calculation = 'solid'
        else:
            self.calculation = 'molecule'

        irspec = np.zeros((no_points, columns))

        for i, line in enumerate(data):
            for j, element in enumerate(line):
                irspec[i, j] = float(element)

        self.irspec = irspec

        return self

    def read_cry_ramspec(self, external_unit):
        import numpy as np

        self.read_external_unit(external_unit)

        data = self.data

        for index, line in enumerate(data):
            data[index] = line.split()

        columns = len(data[0])
        no_points = len(data)

        ramspec = np.zeros((no_points, columns))

        for i, line in enumerate(data):
            for j, element in enumerate(line):
                ramspec[i, j] = float(element)

        self.ramspec = ramspec

        return self
