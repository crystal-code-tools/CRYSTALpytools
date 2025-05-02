#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Objects of input / output files of CRYSTAL. Methods to edit or subtract data
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

    @classmethod
    def read_file(cls, source):
        """
        Instantiate class object from file.

        Args:
            source (str): The name of the input file.
        Returns:
            cls (Crystal_input)
        """
        return cls(source)

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
                       <https://molssi-bse.github.io/basis_set_exchange/bse_cli.html#list-formats>`_.
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
                    raise Exception('Charge definition of element {} is inconsistent with its basis set. {} entries required.'.format(i, len(at.shells)))
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
        output_name (str): Filename
    """

    def __init__(self, output_name=None):
        import re
        import pandas as pd

        if np.all(output_name!=None):
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
            self.df = pd.DataFrame(self.data)
            idx_end = self.df[self.df[0].str.contains(r'^ EEEEEEEEEE TERMINATION')].index

            if len(idx_end) == 0:
                self.terminated = False
                self.eoo = len(self.data)
            else:
                self.terminated = True
                self.eoo = idx_end[0]

    @classmethod
    def read_file(cls, output_name):
        """
        Reads a CRYSTAL output file.

        Args:
            output_name (str): Name of the output file.
        Returns:
            cls (Crystal_output)
        """
        return cls(output_name=output_name)

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
        import pandas as pd

        dimen_line = self.df[self.df[0].str.contains(
            r'\sGEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM'
        )].index

        if len(dimen_line) == 0:
            raise Exception('Invalid file. Dimension information not found.')

        ndimen = self.df[0][dimen_line].map(lambda x: x.strip().split()).tolist()
        self.dimensionality = int(ndimen[0][9])
        return self.dimensionality

    def get_symmops(self, initial=True):
        """
        Return the symmetry operators

        Args:
            initial (bool): *Optimization Only* Read symmetry operators of the
                initial or the final geometry.

        Returns:
            self.symmops (numpy.ndarray): Symmetry operators
        """
        import numpy as np

        self.n_symmops = 0
        self.symmops = np.array([])

        symmops = []
        symmtitle = self.df[self.df[0].str.contains(
            r'^ \*\*\*\*\s+[0-9]+\s*SYMMOPS \- TRANSLATORS IN FRACTIONAL UNITS')].index
        if initial == True: titleidx = symmtitle[0]
        else: titleidx = symmtitle[-1]

        self.n_symmops = int(self.df[0][titleidx].strip().split()[1])
        if self.n_symmops > 0:
            block = self.df[0].loc[titleidx+3:titleidx+3+self.n_symmops-1].map(
                lambda x: x.strip().split()[2:]
            ).tolist()
            self.symmops = np.array(block, dtype=float).reshape([self.n_symmops, 4, 3])
        else:
            self.symmops = []
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
        import os, warnings
        import pandas as pd
        import numpy as np
        from CRYSTALpytools.base.output import GeomBASE
        from CRYSTALpytools.convert import cry_pmg2gui
        from CRYSTALpytools.crystal_io import Crystal_gui

        # Get geometry
        # Use atom coords to read molecule geometries. Go 4 lines up for periodic systems
        coord_lines = self.df[self.df[0].str.contains(r'^\s*ATOMS IN THE ASYMMETRIC UNIT')].index
        if len(coord_lines) == 0:
            raise Exception('Geometry information not found.')

        ## only for developers. If initial is an integer, read the i th geometry from output
        if isinstance(initial, bool):
            if initial == True: bg_line = coord_lines[0] - 4
            else: bg_line = coord_lines[-1] - 4
        else:
            bg_line = coord_lines[int(initial)] - 4

        struc = GeomBASE.read_geom(self.df[0][bg_line:])

        # Get the last lattice matrix: structure obtained by GeomBASE might be rotated.
        ndimen = self.get_dimensionality()

        if ndimen != 0:
            lattice_line = self.df[self.df[0].str.contains(r'^ DIRECT LATTICE VECTORS CARTESIAN')].index
            # always use the last one. Lattice vector is only used to indicate the direction.
            lattice_line = lattice_line[-1]
            a_crys = np.array(self.df[0][lattice_line+2].strip().split(), dtype=float)
            a_pmg = struc.lattice.matrix[0, :]
            # Rotate the geometry back
            struc = struc.rot_cel(a_pmg, a_crys)

        self.geometry = struc

        # Write gui files
        if write_gui == True:
            # Conventional atomic numbers
            zconv = [[i, self.geometry.species_Z[i]] for i in range(self.geometry.num_sites)]
            if np.all(gui_name==None):
                gui_name = os.path.splitext(self.name)[0]
                gui_name = '{}.gui'.format(gui_name)

            if symmetry == 'pymatgen':
                gui = cry_pmg2gui(struc, gui_file=gui_name, symmetry=True, zconv=zconv, **kwargs)
            elif np.all(symmetry==None):
                gui = cry_pmg2gui(struc, gui_file=gui_name,  symmetry=False, zconv=zconv)
            elif symmetry == 'initial':
                self.get_symmops(initial=True)
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

    def get_trans_matrix(self):
        """
        Get cell transformation matrix

        Returns:
            self.trans_matrix (np.ndarray): 3\*3 array of supercell
                expansion matrix
        """
        import pandas as pd
        import numpy as np

        ndimen = self.get_dimensionality()
        self.trans_matrix = np.eye(3, dtype=float)

        mx_line = self.df[self.df[0].str.contains(r'^\s+EXPANSION MATRIX OF PRIMITIVE CELL')].index
        mx = self.df[0][mx_line[0]+1:mx_line[0]+4].map(lambda x: x.strip().split()).tolist()
        mx = np.array([i[1:] for i in mx], dtype=float)

        self.trans_matrix[:ndimen, :ndimen] = mx[:ndimen, :ndimen]
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
        The wrapper of ``get_scf_convergence`` and ``get_opt_convergence``.
        For analysing the geometry and energy convergence.

        .. note::

            It might not work well with SCF / OPT cycles of multiple systems
            such as PREOPTGEOM + EOS calculations

        .. note::

            For optimizations, it only returns to optimization convergence. SCF
            convergence is not available. Call ``get_opt_convergence()`` if needed.

        Args:
            history (bool): Deprecated. Always return to history.

        Returns:
            self (Crystal_output): New attributes listed below
            self.final_energy (float) The converged SCF / Opt energy. Unit: eV

        For other attributes, see ``get_scf_convergence()`` and
        ``get_opt_convergence()`` methods on the same page.
        """
        import re
        import warnings

        opttitle = self.df[self.df[0].str.contains(r'^\s*[A-Z]+ OPTIMIZATION - POINT')].index
        if len(opttitle) == 0: # SCF only
            self.get_scf_convergence(scflog=None)
        else:
            self.get_opt_convergence(scf_history=False, scflog=None)
        return self

    def get_final_energy(self):
        """
        Get the final energy of the system. A wrapper of ``self.get_convergence``.

        Returns:
            self.final_energy (float): The final energy of the system in eV.
        """
        self.get_convergence()
        return self.final_energy

    #### SCF ####

    def get_scf_convergence(self, all_cycles=False, scflog=None):
        """
        Returns the scf convergence energy and energy difference.

        .. note::

            With empirical corrections (DFT-D or gCP), the ``scf_energy`` refers
            to DFT energy only. The ``final_energy`` attribute includes
            corrections at the final step. ``final_energy`` is valid only if
            ``scf_status='converged'``.

        .. note::

            For optimizations, using ``get_opt_convergence()`` with
            ``scf_history=True`` is suggested.

        Args:
            all_cycles (bool): Deprecated.
            scflog (str): Path to 'SCFOUT.LOG'. For optimizations, ``scflog=None``
                returns to either the initial SCF convergence or every SCF step
                if 'ONELOG' keyword is enabled. This option helps to get SCF history
                from 'SCFOUT.LOG' file.

        Returns:
            self (Crystal_output): New attributes listed below
            self.scf_cycles (int|list): Number of cycles. Array if
                multiple SCF cycles are read.
            self.scf_status (str|list): 'terminated', 'converged',
                'too many cycles' and 'unknown'. List if multiple SCF cycles are read.
            self.scf_energy (array|list): SCF energy convergence. Unit: eV
            self.scf_deltae (array|list): Energy difference. Unit: eV
            self.final_energy (float): Last step energy with corrections. Unit: eV
        """
        import numpy as np
        import pandas as pd
        import copy
        from CRYSTALpytools.base.output import SCFBASE
        from CRYSTALpytools.units import H_to_eV

        self.scf_cycles = []
        self.scf_status = []
        self.scf_energy = []
        self.scf_deltae = []

        nSCF, SCFrange = SCFBASE.get_SCF_blocks(self.df[0])
        for i in SCFrange:
            output = SCFBASE.read_convergence(self.df[0].loc[i[0]:i[1]])
            self.scf_cycles.append(output[0])
            self.scf_status.append(output[1])
            self.scf_energy.append(output[2])
            self.scf_deltae.append(output[3])

        if np.all(scflog!=None):
            df = pd.DataFrame(open(scflog, 'r'))
            nSCFlog, SCFrangelog = SCFBASE.get_SCF_blocks(df[0])
            nSCF += nSCFlog
            for i in SCFrangelog:
                output = SCFBASE.read_convergence(df[0].loc[i[0]:i[1]])
                self.scf_cycles.append(output[0])
                self.scf_status.append(output[1])
                self.scf_energy.append(output[2])
                self.scf_deltae.append(output[3])

        if nSCF == 1:
            self.scf_cycles = int(self.scf_cycles[-1])
            self.scf_status = str(self.scf_status[-1])
            self.scf_energy = self.scf_energy[-1]
            self.scf_deltae = self.scf_deltae[-1]
        # final energy, with corrections
        status = copy.deepcopy(self.scf_status)
        status = np.array(status, ndmin=1)
        if status[-1] == 'converged':
            finalline = self.df[self.df[0].str.contains(r'^\s*TOTAL ENERGY \+ ')].index
            optline = self.df[self.df[0].str.contains(r'^\s*\* OPT END')].index
            if len(optline) == 0 and len(finalline) == 0:
                self.final_energy = self.scf_energy[-1]
            elif len(optline) == 0:
                self.final_energy = H_to_eV(float(self.df[0][finalline[-1]].strip().split()[-1]))
            else:
                self.final_energy = H_to_eV(float(self.df[0][optline[-1]].strip().split()[-4]))
        else:
            self.final_energy = 0.
        return self

    def get_fermi_energy(self, history=False, scflog=None):
        """
        Returns the system Fermi energy / valance band maximum (insulators).
        For ``history=True``, Fermi energy at step 0 is always 0 as no Fermi
        energy is solved. Similary, steps with SPINLOCK or level shifting also
        return to 0 Fermi energy.

        Args:
            history (bool): Whether to read the convergence history of Fermi energy.
            scflog (str): Path to 'SCFOUT.LOG'. For optimizations, ``scflog=None``
                returns to either the initial Fermi energy or energies of every SCF
                if 'ONELOG' keyword is enabled. This option helps to get SCF history
                from 'SCFOUT.LOG' file.

        Returns:
            self.fermi_energy (float|array): Fermi energy of the system. For
                spin-polarized insulators, the highest VBM is returned.
            self.spin_pol (bool): *Not returned but defined* Whether the
                calculation is spin-polarized.
        """
        from CRYSTALpytools.base.output import SCFBASE
        import numpy as np

        nSCF, SCFrange = SCFBASE.get_SCF_blocks(self.df[0])
        if history == True:
            self.fermi_energy = []
            for i in SCFrange:
                output = SCFBASE.read_convergence(self.df[0].loc[i[0]:i[1]])
                self.fermi_energy.append(output[5])
            self.spin_pol = output[4]
        else:
            output = SCFBASE.read_convergence(self.df[0].loc[SCFrange[-1,0] : SCFrange[-1,1]])
            self.spin_pol = output[4]
            self.fermi_energy = output[5][-1]

        if np.all(scflog!=None):
            df = pd.DataFrame(open(scflog, 'r'))
            nSCFlog, SCFrangelog = SCFBASE.get_SCF_blocks(df[0])
            nSCF += nSCFlog
            if history == True:
                for i in SCFrangelog:
                    output = SCFBASE.read_convergence(df.loc[i[0]:i[1]])
                    self.fermi_energy.append(output[5])
                self.spin_pol = output[4]
            else:
                output = SCFBASE.read_convergence(df.loc[SCFrangelog[-1,0] : SCFrangelog[-1,1]])
                self.spin_pol = output[4]
                self.fermi_energy = output[5][-1]

        if history == True:
            self.fermi_energy = np.array(self.fermi_energy)
            if nSCF == 1: self.fermi_energy = self.fermi_energy[0]
        return self.fermi_energy

    def get_band_gap(self, history=False, scflog=None):
        """
        Returns the system band gap. For ``history=True``, gap at step 0 is
        always 0 as no energy gap is solved. Similary, steps with SPINLOCK or
        level shifting also return to 0 gap.

        Args:
            history (bool): Whether to read the convergence history of band gap.
            scflog (str): Path to 'SCFOUT.LOG'. For optimizations, ``scflog=None``
                returns to either the initial band gap or gaps of every SCF if
                'ONELOG' keyword is enabled. This option helps to get SCF history
                from 'SCFOUT.LOG' file.

        Returns:
            self.band_gap (float|array): Band gap of the system. For
                spin-polarized systems, ``self.band_gap`` would be either a 2\*1
                array (``history=False``) or a nCYC\*2 array (``history=True``).
            self.spin_pol (bool): *Not returned but defined* Whether the
                calculation is spin-polarized.
        """
        from CRYSTALpytools.base.output import SCFBASE
        import numpy as np

        nSCF, SCFrange = SCFBASE.get_SCF_blocks(self.df[0])
        if history == True:
            self.band_gap = []
            for i in SCFrange:
                output = SCFBASE.read_convergence(self.df[0].loc[i[0]:i[1]])
                self.band_gap.append(output[6])
            self.spin_pol = output[4]
        else:
            output = SCFBASE.read_convergence(self.df[0].loc[SCFrange[-1,0] : SCFrange[-1,1]])
            self.spin_pol = output[4]
            self.band_gap = output[6][-1]

        if np.all(scflog!=None):
            df = pd.DataFrame(open(scflog, 'r'))
            nSCFlog, SCFrangelog = SCFBASE.get_SCF_blocks(df[0])
            nSCF += nSCFlog
            if history == True:
                for i in SCFrangelog:
                    output = SCFBASE.read_convergence(df.loc[i[0]:i[1]])
                    self.band_gap.append(output[6])
                self.spin_pol = output[4]
            else:
                output = SCFBASE.read_convergence(df.loc[SCFrangelog[-1,0] : SCFrangelog[-1,1]])
                self.spin_pol = output[4]
                self.band_gap = output[6][-1]

        if history == True:
            self.band_gap = np.array(self.band_gap)
            if nSCF == 1: self.band_gap = self.band_gap[0]
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
            self.mulliken_charges = np.vstack([apb, (apb+amb)/2, (apb-amb)/2]).transpose()
            self.spin_pol = True
        return self.mulliken_charges

    #### Optimization ####

    def get_opt_convergence(self, primitive=False, scf_history=False, scflog=None,
                            write_gui=False, gui_name=None, symmetry='pymatgen',
                            **kwargs):
        """
        Returns optimisation convergence.

        Args:
            primitive (bool): Deprecated. Use ``get_primitive_geometry()`` for
                initial and last geometries. The method here always returns to
                the actual structure optimized.
            scf_history (bool): Read SCF history of each optimisation step. Also
                refer to the ``get_scf_convergence()`` method.
            scflog (str): Path to 'SCFOUT.LOG'. For optimizations, ``scflog=None``
                returns to either the initial SCF convergence or every SCF step
                if 'ONELOG' keyword is enabled. This option helps to get SCF history
                from 'SCFOUT.LOG' file.
            write_gui (bool): If True, write .gui file of each step
            gui_name (str): Valid only if ``write_gui = True``. Gui file is
                named as 'gui_name-optxxx.gui'. If None, use 'basename-optxxx.gui'.
                The basename is the same as output.
            symmetry (str): Valid only if ``write_gui = True``. 'pymatgen' to
                use symmetry info from a pymatgen SpacegroupAnalyzer; 'initial'
                to use symmstry information on output file. If None, no symmstry.
                Otherwise it is taken from the existing gui file.
            **kwargs: Valid only if ``write_gui = True`` and ``symmetry = 'pymatgen'``.
                Passed to Pymatgen SpacegroupAnalyzer object.

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
            self.final_energy (float): Last step energy with corrections. Unit: eV
        """
        from CRYSTALpytools.crystal_io import Crystal_gui
        from CRYSTALpytools.convert import cry_pmg2gui
        from CRYSTALpytools.base.output import OptBASE, GeomBASE
        import numpy as np
        import warnings

        # Initial geometry
        ndimen = self.get_dimensionality()
        struc0 = self.get_geometry(initial=True, write_gui=False)

        # steps
        self.opt_cycles, OPTrange, self.opt_status = OptBASE.get_opt_block(self.df[0])
        self.opt_energy = np.zeros([self.opt_cycles,], dtype=float)
        self.opt_deltae = np.zeros([self.opt_cycles,], dtype=float)
        self.opt_geometry = [struc0,]
        self.opt_maxgrad = np.zeros([self.opt_cycles,], dtype=float)
        self.opt_rmsgrad = np.zeros([self.opt_cycles,], dtype=float)
        self.opt_maxdisp = np.zeros([self.opt_cycles,], dtype=float)
        self.opt_rmsdisp = np.zeros([self.opt_cycles,], dtype=float)

        for i in range(self.opt_cycles):
            bg = OPTrange[i, 0]; ed = OPTrange[i, 1]
            try:
                output = OptBASE.read_opt_block(self.df[0].loc[bg:ed])
            except IndexError: # Unfinished Opt step
                warnings.warn('Hit the unfinished step. Properties except geometry are set to 0.',
                              stacklevel=2)
                self.opt_geometry.append(GeomBASE.read_geom(self.df[0].loc[bg:ed]))
                continue

            self.opt_energy[i] = output[0]
            self.opt_deltae[i] = output[1]
            if i > 0: self.opt_geometry.append(output[2])
            self.opt_maxgrad[i] = output[3]
            self.opt_rmsgrad[i] = output[4]
            self.opt_maxdisp[i] = output[5]
            self.opt_rmsdisp[i] = output[6]

        # SCF history
        if np.all(scflog!=None): scf_history = True
        if scf_history == True: self.get_scf_convergence(scflog=scflog)

        # final energy
        self.final_energy = self.opt_energy[-1]

        # Geometry
        ## Lattice matrix rotation issue
        if ndimen != 0:
            lattice_line = self.df[self.df[0].str.contains(r'^ DIRECT LATTICE VECTORS CARTESIAN')].index
            # always use the last one. Lattice vector is only used to indicate the direction.
            lattice_line = lattice_line[-1]
            a_crys = np.array(self.df[0][lattice_line+2].strip().split(), dtype=float)
            a_pmg = struc0.lattice.matrix[0, :]
            # Rotate the geometry back
            for idx_s, s in enumerate(self.opt_geometry):
                self.opt_geometry[idx_s] = s.rot_cel(a_crys, a_pmg)
        ## Write gui files
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

    def get_forces(self, initial=True, grad=False, scflog=None):
        """
        Read forces.

        .. note::

            To get forces of the last optimization step, 'SCFOUT.LOG' file or
            the 'ONELOG' keyword are needed.

        Args:
            initial (bool): Deprecated.
            grad (bool): Return gradient convergence history. For optimizations
                only.
            scflog (str): Path to 'SCFOUT.LOG'. For optimizations, ``scflog=None``
                returns forces from the initial SCF convergence or forces of
                every SCF step if 'ONELOG' keyword is enabled. This option helps
                to get forces of every step from 'SCFOUT.LOG' file.

        Returns:
            self (Crystal_output): New attributes listed below
            self.forces_atoms (array): natom\*3 array. Atomic forces. Unit: Hartree/Bohr
            self.forces_cell (array): 3\*3 array. Cell forces, 3D only. Unit: Hartree/Bohr
            self.opt_maxgrad (array): Maximum gradient convergence. Unit: Hartree/Bohr
            self.opt_rmsgrad (array): RMS gradient convergence. Unit: Hartree/Bohr
        """
        import warnings, re
        import numpy as np
        import pandas as pd

        title = self.df[self.df[0].str.contains(r'^\s*OPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPTOPT')].index
        if initial == False or grad == True:
            if len(title) == 0:
                warnings.warn('Not a geometry optimisation: Set initial = True and grad = False', stacklevel=2)
                initial = True; grad = False

        if grad == True:
            self.get_opt_convergence()

        forcetitle = self.df[self.df[0].str.contains(r'^\s*\*\s+FORCE CALCULATION')].index
        if len(forcetitle) == 0: raise Exception('Force calculation not found.')

        # dimensionalities
        atomtitle = self.df[self.df[0].str.contains(
            r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)')].index.to_numpy(dtype=int)
        celltitle = self.df[self.df[0].str.contains(
            r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR')].index.to_numpy(dtype=int)
        natomf = len(atomtitle)
        ncellf = len(celltitle)
        struc = self.get_geometry(initial=True, write_gui=False)
        natom = struc.num_sites
        if np.all(scflog!=None):
            df = pd.DataFrame(open(scflog))
            atomtitleLOG = df[df[0].str.contains(
                r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)')].index.to_numpy(dtype=int)
            celltitleLOG = df[df[0].str.contains(
                r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR')].index.to_numpy(dtype=int)
            natomf += len(atomtitleLOG)
            ncellf += len(celltitleLOG)

        if natomf != 0:
            self.forces_atoms = np.zeros([natomf, natom, 3])
            for i in range(len(atomtitle)):
                atomf = self.df[0].loc[atomtitle[i]+2 : atomtitle[i]+1+natom].map(lambda x: x.strip().split()).tolist()
                self.forces_atoms[i, :, :] = np.array([j[2:] for j in atomf], dtype=float)
            if np.all(scflog!=None):
                for i in range(natomf-len(atomtitle)):
                    atomf = df[0].loc[atomtitleLOG[i]+2 : atomtitleLOG[i]+1+natom].map(lambda x: x.strip().split()).tolist()
                    self.forces_atoms[i+len(atomtitle), :, :] = np.array([j[2:] for j in atomf], dtype=float)
        else:
            self.forces_atoms = []

        if ncellf != 0:
            self.forces_cell = np.zeros([ncellf, 3, 3])
            for i in range(len(celltitle)):
                cellf = self.df[0].loc[celltitle[i]+4 : celltitle[i]+6].map(lambda x: x.strip().split()).tolist()
                self.forces_cell[i, :, :] = np.array(cellf, dtype=float)
            if np.all(scflog!=None):
                for i in range(ncellf-len(celltitle)):
                    cellf = df[0].loc[celltitleLOG[i]+4 : celltitleLOG[i]+6].map(lambda x: x.strip().split()).tolist()
                    self.forces_cell[i+len(celltitle), :, :] = np.array(cellf, dtype=float)
        else:
            self.forces_cell = []

        if natomf == 1: self.forces_atoms = self.forces_atoms[0]
        if ncellf == 1: self.forces_cell = self.forces_cell[0]
        return self

    #### Lattice Dynamics ####

    def get_phonon(self, imaginary_tol=-1e-4, q_overlap_tol=1e-4, eigvec='unsaved',
                   **kwargs):
        """
        Read phonon-related properties from output file. This method was
        developed as a basic I/O function for phonon vibration and thermodynamic
        analysis. For phonon band / dos or IR / Raman intensities or spectra,
        please refer to the ``get_phonon_band``, ``get_phonon_dos`` and
        ``get_spectra`` methods.

        .. note::

            In QHA calculations, ``self.nqpoint`` refer to harmonic phonons
            computed at :math:`\\Gamma`, as the current version (CRYSTAL23)
            supports direct-space approach only.

        Available eigenvector formats for ``eigvec``:

        * 'unsaved', Do not save eigenvector.  
        * 'original', Mass-weighted phased eigenvectors normalized to 1, , i.e.,
           eigenvectors from diagonalizing dynamic matrices.  
        * 'classical'，Mass-weighted phased eigenvectors normalized to classical
           amplitude in :math:`\\AA`.  
        * 'mass unweight', Mass-unweighted phased eigenvectors normalized to 1.

        Args:
            imaginary_tol (float): The threshold of negative frequencies to
                remove. 'None' for keeping the original data.
            q_overlap_tol (float): The threshold of overlapped q points to remove,
                defined as the 2nd norm of the difference of fractional q vectors.
                'None' for keeping the original data.
            eigvec (str): Eigenvectors saved. 'unsaved', 'original', 'classical',
                'mass unweight'. See above.
            \*\*kwargs: Developers only. See below.
            read_spec (bool): Read spectra information at Gamma. Extra 1\*nmode
                arrays of 'IR', 'intens' and 'Raman' are returned.
            read_eigvt : Deprecated. Use ``eigvec`` instead.
            eigvt_amplitude: Deprecated. Use ``eigvec`` instead.
            rm_imaginary: Deprecated. Use ``imaginary_tol`` instead.
            rm_overlap: Deprecated. Use ``q_overlap_tol=None`` instead.

        Returns:
            self (Crystal_output): New attributes listed below
            self.phonon (Phonon): ``phonons.Phonon`` object. In practice calling
                this attribute is suggested. The followings are for compatbility.
            self.edft (float|arra): :math:`E_{0}` Energy with empirical
                corrections. An 1D array for QHA output. Unit: kJ/mol.
            self.nqpoint (int): Number of q points. For QHA, number of calculations.
            self.qpoint (array[float): nQpoint\*4 array. The first three
                elements are fractional coordinates and the last is weight, i.e.,
                number of equivalent q points. For QHA, \[0, 0, 0, 1\] repeated
                for nqpoint times.
            self.nmode (int): Number of modes at q point.
            self.frequency (array[float]): nQpoint\*nMode array ofvibrational
                frequency. Unit: THz
            self.mode_symm (array[str]): nQpoint\*nMode array of the
                irreducible representations. In Mulliken symbols.
            self.eigenvector (array[complex]): *``read_eigvt = True only``*
                nqpoint\*nmode\*natom\*3 array of eigenvectors.
        """
        import re
        import pandas as pd
        from CRYSTALpytools.phonons import Phonon

        is_freq = False
        found_anti = True
        read_spec = False
        scale = 1. # scale is for compatibility only
        edft = []
        nqpoint = 0
        qpoint = []
        nmode = []
        frequency = []
        mode_symm = []
        eigenvector = []

        # Developer and deprecated arguments:
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k == 'rm_imaginary' and v == False: imaginary_tol = None; continue
            if k == 'rm_overlap' and v == False: q_overlap_tol = None; continue
            if k == 'read_eigvt':
                if v == False: eigvec = 'unsaved'; continue
                elif v == 'classical': eigvec = 'classical'; continue
                else: eigvec = 'original'; scale=v; continue
            if 'read_spec' in kwargs:
                read_spec = True
                save_spec = False # flag for gamma point
                intens = []
                IR = []
                Raman = []

        # Whether is a frequency file
        title = self.df[self.df[0].str.contains(
            r'^\s*\+\+\+\sSYMMETRY\sADAPTION\sOF\sVIBRATIONAL\sMODES\s\+\+\+')].index
        if len(title) == 0: raise Exception('Not a frequency calculation.')

        # E_0 with empirical corrections
        edft_idx = self.df[self.df[0].str.contains(r'^\s+CENTRAL POINT')].index
        edft = [i[2] for i in self.df[0][edft_idx].map(lambda x: x.strip().split()).tolist()]
        edft = units.H_to_eV(np.array(edft, dtype=float))

        # Geometry: Reduce SCELPHONO for dispersions. Not for Gamma point.
        band_title = self.df[self.df[0].str.contains(r'^\s*\*\s+LIST OF THE K POINTS USED FO[F,R] PHONON DISPERSION\.\s+\*\s*$')].index
        scelphono = self.df[self.df[0].str.contains(r'\s*\*+\s+ATOMS IN THE SUPERCELL REORDERED FOR PHONON CALCULATION\s*$')].index
        if len(scelphono) > 0 and len(band_title) > 0:
            struc = self.get_primitive_geometry(initial=False)
            sstruc = self.get_geometry(initial=False)
            edft /= sstruc.num_sites / struc.num_sites; del sstruc
        else:
            struc = self.get_geometry(initial=False)

        # Q point info + frequency
        empty_line = self.df[self.df[0].map(lambda x: x.strip() == '')].index.to_numpy(dtype=int)
        if len(band_title) > 0:
            ## Dispersion q point info
            kheader = self.df[self.df[0].str.contains(r'^\s*\*\s+K\s+WEIGHT\s+COORD\s+\*\s*$')].index[0]
            kend = self.df[self.df[0].str.contains(r'^\s*\*\s+WITH SHRINKING FACTORS: IS1 =\s+[0-9]+')].index[0]
            nqpoint = int(kend - kheader - 1)
            qpoint = np.array(self.df[0][kheader+1:kend].map(lambda x: x.strip().split()[1:-1]).tolist(),
                              dtype=float)
            is1 = float(self.df[0][kend].strip().split()[6])
            is2 = float(self.df[0][kend].strip().split()[9])
            is3 = float(self.df[0][kend].strip().split()[12])
            qpoint = [[i[2]/is1, i[3]/is2, i[4]/is3, i[1]] for i in qpoint]
            freq_header = self.df[self.df[0].str.contains(r'\s+DISPERSION K POINT NUMBER')].index
            # Generate list for IRREP symbols
            IRREP = []
            irreptitle = self.df[self.df[0].str.contains(r'^\s+IRREP/CLA')].index
            bg = irreptitle[0]
            ed = empty_line[np.where(empty_line>bg)[0][0]]
            dfmini = self.df[0][bg+4:ed]
            IRREP = dfmini.map(lambda x: x[5:9].strip()).tolist()
        else:
            ## Gamma point / Gamma point QHA.
            nqpoint = len(edft)
            qpoint = [[0., 0., 0., 1.0] for i in range(nqpoint)]
            ## frequency
            freq_header = self.df[self.df[0].str.contains(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+')].index
            IRREP = []
        qpoint = np.array(qpoint)

        # Frequency
        for bg in freq_header:
            ed = empty_line[np.where(empty_line>=bg+2)[0][0]]
            dfmini = self.df[0][bg:ed]
            block = dfmini[dfmini.str.contains(r'^\s*[0-9]+\-\s+')]
            idx = 0; empty = empty_line[np.where(empty_line>=ed)[0]]
            while len(block) == 0:
                ed = empty[idx]
                dfmini = self.df[0][bg:ed]
                block = dfmini[dfmini.str.contains(r'^\s*[0-9]+\-\s+')]
                idx += 1

            mode1 = block.map(lambda x: x[0:5]).to_numpy(dtype=int)
            mode2 = block.map(lambda x: x[6:10]).to_numpy(dtype=int)
            # in cm-1 for more decimal places
            freq_tmp = block.map(lambda x: x[24:36]).to_numpy(dtype=float)
            symm_tmp = block.map(lambda x: x[49:52].strip()).tolist()
            # for phonon dispersions, convert indices into symbols
            if len(re.findall(r'[A-Z,a-z]', symm_tmp[0])) == 0:
                symm_tmp = [IRREP[int(i)-1] for i in symm_tmp]
            if read_spec == True:
                if 'I' in data[data.index[0]] or 'A' in data[data.index[0]]:
                    save_spec = True
                    intens_tmp = data.map(lambda x: x[59:68]).to_numpy(dtype=float)
                    IR_tmp = data.map(lambda x: x[56]=='A').to_numpy(dtype=bool)
                    Raman_tmp = data.map(lambda x: x[73]=='A').to_numpy(dtype=bool)
            # repeat
            freq = []; symm=[]
            count = 0
            for m1, m2 in zip(mode1, mode2):
                freq.append(freq_tmp[count])
                symm.append(symm_tmp[count])
                for i in range(m1, m2):
                    freq.append(freq_tmp[count])
                    symm.append(symm_tmp[count])
                ## Raman and IR
                if read_spec == True and save_spec == True:
                    intens.append(intens_tmp[count])
                    IR.append(IR_tmp[count])
                    Raman.append(Raman_tmp[count])
                    for i in range(m1, m2):
                        intens.append(intens_tmp[count])
                        IR.append(IR_tmp[count])
                        Raman.append(Raman_tmp[count])
                count += 1
            # cm to thz
            freq = units.cm_to_thz(np.array(freq, dtype=float))
            symm = np.array([i for i in symm])
            if read_spec == True and save_spec == True:
                intens = np.array(intens, dtype=float)
                IR = np.array(IR, dtype=bool)
                Raman = np.array(Raman, dtype=bool)
                save_spec == False
            frequency.append(freq)
            mode_symm.append(symm)

        frequency = np.array(frequency, dtype=float)
        mode_symm = np.array(mode_symm, dtype=str)
        nmode = frequency.shape[1]

        ## eigenvector
        if eigvec != 'unsaved':
            eigvt_header = self.df[self.df[0].str.contains(r'^\s*NORMAL MODES NORMALIZED TO')].index.tolist()
            eigvt_header += self.df[self.df[0].str.contains(r'^\s*MODES IN PHASE')].index.tolist()
            eigvt_header = np.sort(np.array(eigvt_header, dtype=int))

            if len(band_title) > 0: # DISP
                eigvt_end = freq_header.tolist()[1:]
                eigvt_end.append(self.eoo)
            else: # Gamma frequency
                if nqpoint == 1:
                    eigvt_end = [self.eoo]
                else: # QHA
                    dftitle = self.df[self.df[0].str.contains(r'\s*\*\s+CELL DEFORMATION\s*$')].index
                    eigvt_end = dftitle.tolist()[1:] + [self.eoo]

            countblock = 0
            for bg, ed in zip(eigvt_header, eigvt_end):
                dfmini = self.df[0][bg:ed]
                img = dfmini[dfmini.str.contains(r'MODES IN ANTI\-PHASE')].index
                dfmini = dfmini[dfmini.str.contains(r'[X,Y,Z]\s+\-*[0-9]\.')]
                nlines = len(dfmini)
                if len(img) == 1: # complex
                    rvt = dfmini[0:int(nlines/2)].map(lambda x: x[13:].strip().split()).tolist()
                    ivt = dfmini[int(nlines/2):].map(lambda x: x[13:].strip().split()).tolist()
                else: # real
                    rvt = dfmini.map(lambda x: x[13:].strip().split()).tolist()
                    ivt = []

                nblock = int(np.ceil(nmode / 6))
                nline = int(len(rvt) / nblock) # nline per block of 6 modes

                rvttmp = np.array(rvt[0:nline], dtype=float)
                for i in range(1, nblock):
                    rvttmp = np.hstack(
                        [rvttmp, np.array(rvt[int(nline*i):int(nline*(i+1))], dtype=float)]
                    )
                rvttmp = rvttmp.T.reshape([nmode, int(nline/3), 3], order='C')

                if len(ivt) != 0:
                    ivttmp = np.array(ivt[0:nline], dtype=float)
                    for i in range(1, nblock):
                        ivttmp = np.hstack(
                            [ivttmp, np.array(ivt[int(nline*i):int(nline*(i+1))], dtype=float)]
                        )
                    ivttmp = ivttmp.T.reshape([nmode, int(nline/3), 3], order='C')
                    eigvt = rvttmp + ivttmp * 1j
                else:
                    eigvt = rvttmp + 0j
                countblock += 1
                eigenvector.append(eigvt)
            ## Remove classical amplitude
            eigenvector = units.au_to_angstrom(np.array(eigenvector))
            eigenvector = Phonon.get_eigvec(struc, frequency, eigenvector, 'remove classical')
        else:
            eigenvector = np.array([])

        if len(edft) == 1: edft = edft[0]
        self.phonon = Phonon(struc, edft, qpoint, frequency, mode_symm, eigenvector)
        if imaginary_tol != None:
            self.phonon.clean_imaginary(threshold=imaginary_tol)
        if q_overlap_tol != None:
            self.phonon.clean_q_overlap(threshold=q_overlap_tol)

        self.edft = units.H_to_kjmol(units.eV_to_H(self.phonon.u_0))
        self.nqpoint = self.phonon.nqpoint
        self.qpoint = self.phonon.qpoint
        self.nmode = self.phonon.nmode
        self.frequency = self.phonon.frequency
        self.mode_symm = self.phonon.mode_symm
        if eigvec == 'original':
            self.eigenvector = self.phonon.eigenvector * scale
        elif eigvec == 'classical':
            self.eigenvector = self.phonon.classical_eigvec()
        elif eigvec == 'mass unweight':
            self.eigenvector = self.phonon.unweight_eigvec() * scale
        elif eigvec == 'unsaved':
            self.eigenvector = self.phonon.eigenvector
        if read_spec == True:
            return self, intens, IR, Raman
        else:
            return self

    def get_phonon_band(self, q_overlap_tol=1e-4):
        """
        Read phonon bands from CRYSTAL output file.

        Args:
            q_overlap_tol (float): The threshold for overlapped k points. Only
                used for getting tick positions.

        Returns:
            self.pband (PhononBand): ``phonons.PhononBand`` object.
        """
        from CRYSTALpytools.phonons import PhononBand
        import pandas as pd
        import numpy as np

        band_title = self.df[self.df[0].str.contains(r'^\s*\*\s+PHONON BANDS\s+\*\s*$')].index
        scelphono = self.df[self.df[0].str.contains(r'\s*\*+\s+ATOMS IN THE SUPERCELL REORDERED FOR PHONON CALCULATION\s*$')].index
        if len(band_title) == 0:
            raise Exception('Not a phonon band file')

        self.get_phonon(read_eigvt=False, rm_imaginary=False, rm_overlap=False)
        if len(scelphono) > 0:
            struc = self.get_primitive_geometry(initial=False)
        else:
            struc = self.get_geometry(initial=False)

        k_path3d = np.vstack([i[0] for i in self.qpoint])
        bands = np.reshape(self.frequency.transpose(), [self.nmode, self.nqpoint, 1])
        recp_latt = struc.lattice.reciprocal_lattice.matrix

        # get labels and 1D k path
        headline = self.df[self.df[0].str.contains(r'\s*FROM K\s+\(.*WITH DENOMINATOR')].index
        headline = self.df[0][headline].map(lambda x: x.strip().split()).tolist()

        steps = self.df[self.df[0].str.contains(r'\s*PHONONS ALONG PATH\:\s+[0-9]+\s+NUMBER OF K POINTS\:\s+[0-9]+')].index
        steps = self.df[0][steps].map(lambda x: x.strip().split()[-1]).to_numpy(dtype=int)

        k_path = np.array([], dtype=float); tick_pos = []; tick_label = []
        tick_pos3d = []; totlen = 0.
        for iline, line in enumerate(headline):
            pt1 = np.array(line[3:6], dtype=int)
            pt2 = np.array(line[10:13], dtype=int)
            denominator = int(line[-1])
            seglen = np.linalg.norm((pt2-pt1) @ recp_latt)
            tick_pos3d.append(pt1 / denominator)
            tick_pos3d.append(pt2 / denominator)
            tick_label.append('({:d} {:d} {:d})/{:d}'.format(pt1[0], pt1[1], pt1[2], denominator))
            tick_label.append('({:d} {:d} {:d})/{:d}'.format(pt2[0], pt2[1], pt2[2], denominator))
            k_path = np.hstack([k_path, np.linspace(totlen, totlen+seglen, steps[iline])])
            tick_pos.append(totlen)
            tick_pos.append(totlen+seglen)
            totlen += seglen

        # overlapped points
        for i in range(len(tick_label)-1, 1, -1):
            dist = tick_pos[i]-tick_pos[i-1]
            if dist <= q_overlap_tol:
                del tick_pos[i], tick_label[i], tick_pos3d[i]

        # instantiation
        self.pband = PhononBand(tick_pos, tick_label, bands, k_path, struc,
                                tick_pos3d=tick_pos3d, k_path3d=k_path3d)
        return self.pband

    def get_phonon_dos(self, read_INS=False, atom_prj=[], element_prj=[]):
        """
        Read phonon density of states from CRYSTAL output file.

        Args:
            read_INS (bool): Read the inelastic neutron scattering spectra,
                instead of the PDOS.
            atom_prj (list): Read the projections of atoms with specified labels.
            element_prj (list): Read projections of elements with specified
                conventional atomic numbers.

        Returns:
            self.pdos (PhononDOS): ``phonons.PhononDOS`` object.
        """
        import pandas as pd
        import numpy as np
        from CRYSTALpytools.units import cm_to_thz, thz_to_cm
        from CRYSTALpytools.phonons import PhononDOS

        # read from .out file
        if read_INS == False:
            header = self.df[self.df[0].str.contains(r'^\s*FREQUENCY \(CM\*\*\-1\)    TOTAL PDOS')].index
        else:
            header = self.df[self.df[0].str.contains(r'^\s*FREQUENCY \(CM\*\*\-1\)    TOTAL NW\-PDOS')].index
            if len(header) == 0:
                raise Exception('File does not have INS phonon DOS.')

        empty_line = self.df[self.df[0].map(lambda x: x.strip() == '')].index.tolist()
        empty_line = np.array(empty_line, dtype=int)
        # total prj
        totprj = []
        for h in header:
            end = empty_line[np.where(empty_line>h+2)[0][0]]
            totprj.append(np.array(
                self.df[0][h+2:end].map(lambda x: x.strip().split()).tolist(),
                dtype=float
            ))
        totprj = np.array(totprj, dtype=float)

        ## projections
        atom_prj = np.array(atom_prj, ndmin=1)
        element_prj = np.array(element_prj, ndmin=1)
        nprj = len(atom_prj) + len(element_prj) + 1
        doss = np.zeros([nprj, totprj.shape[1], 1], dtype=float)
        freq = cm_to_thz(totprj[0, :, 0])
        doss[0, :, 0] = np.sum(totprj[:, :, 1], axis=0)

        ## atomic
        if len(atom_prj) > 0:
            for ia, a in enumerate(atom_prj):
                header = self.df[self.df[0].str.contains(r'^\s\s\sATOM\s+{:d}$'.format(a))].index
                if len(header) == 0:
                    raise ValueError("The specified atom label: '{:d}' is not found.".format(a))
                for h in header:
                    end = empty_line[np.where(empty_line>h+4)[0][0]]
                    doss[ia+1, :, 0] += np.array(
                        self.df[0][h+4:end].map(lambda x: x.strip().split()).tolist(),
                        dtype=float
                    )[:, 1]
        ## species
        if len(element_prj) > 0:
            for ie, e in enumerate(element_prj):
                header = self.df[self.df[0].str.contains(r'^\s\s\sATOMIC NUMBER\s+{:d}$'.format(e))].index
                if len(header) == 0:
                    raise ValueError("The specified atomic number: '{:d}' is not found.".format(e))
                for h in header:
                    end = empty_line[np.where(empty_line>h+4)[0][0]]
                    doss[ie+len(atom_prj)+1, :, 0] += np.array(
                        self.df[0][h+4:end].map(lambda x: x.strip().split()).tolist(),
                        dtype=float
                    )[:, 1]

        doss = thz_to_cm(doss)
        self.pdos = PhononDOS(doss, freq, unit='THz')
        return self.pdos

    def get_spectra(self, specfile, type='infer'):
        """
        Read spectra from IRSPEC.DAT / RAMSPEC.DAT files.

        .. note::

            In principle, the code cannot tell the type of the files. It is
            important to assign the appropriate type to the ``type`` keyword.
            Currently the available options are 'IRSPEC' and 'RAMSPEC'.

        Args:
            specfile (str): File name of spectra data.
            type (str): Type of the file. See above. If 'infer', type is
                inferred from input file name.

        Returns:
            self.\* (IR|Raman): Dependending on the ``type`` keyword, return to
                ``spectra.IR`` or ``spectra.Raman`` objects. Attribute names
                same as ``type`` are set.
        """
        import numpy as np
        from CRYSTALpytools.spectra import IR, Raman
        import warnings

        # sanity check
        accepted_types = ['IRSPEC', 'RAMSPEC']
        if type.lower() == 'infer':
            name = specfile.upper()
            for a in accepted_types:
                if a in name: type = a; break
            if type.lower() == 'infer':
                raise Exception("Type of file cannot be inferred from input file name: '{}'.".format(specfile))

        if type.upper() not in accepted_types:
            raise ValueError("The specified type: '{}' is not vaild.".format(type))
        type = type.upper()

        if not hasattr(self, 'df'):
            warnings.warn('Output file not available. Geometry information missing.',
                          stacklevel=2)
        else:
            if type == 'IRSPEC':
                title = self.df[self.df[0].str.contains(
                    r'^\s*\*\s+CALCULATION OF INFRARED ABSORBANCE \/ REFLECTANCE SPECTRA'
                )].index
                if len(title) == 0:
                    warnings.warn(
                        'IR spectra block is not found in the screen output. Are files from the same calculation?',
                         stacklevel=2
                    )
            elif type == 'RAMSPEC':
                # mute pandas warning
                warnings.filterwarnings("ignore", 'This pattern is interpreted as a regular expression, and has match groups.')
                title = self.df[self.df[0].str.contains(
                    r'^\s*(\<RAMAN\>){11}'
                )].index
                if len(title) == 0:
                    warnings.warn(
                        'Raman spectra block is not found in the screen output. Are files from the same calculation?',
                        stacklevel=2
                    )
        # read file and instantiation
        data = np.loadtxt(specfile)
        if type == 'IRSPEC':
            if data.shape[1] > 3: # crystal IR
                obj = IR(freq=data[:, 0].T, absorb=data[:, 2:6].T, reflec=data[:, 6:].T, type='crystal')
            else: # molecule IR
                obj = IR(freq=data[:, 0].T, absorb=data[2].T, reflec=[], type='molecule')
        elif type == 'RAMSPEC':
            obj = Raman(freq=data[:, 0].T, poly=data[:, 1:4].T, single=data[:, 4:].T)

        setattr(self, type, obj)
        return getattr(self, type)

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
                GPa.m for 2D and GPa.m :math:`^{2}` for 1D (might lead to very
                small numbers).
        """
        import re, warnings

        etitle = self.df[
            self.df[0].str.contains(r'^\s*SYMMETRIZED ELASTIC CONSTANTS')
        ].index.to_numpy(dtype=int)

        empty_line = self.df[
            self.df[0].map(lambda x: x.strip() == '')
        ].index.to_numpy(dtype=int)

        if len(etitle) == 0:
            raise Exception('Elastic tensor not found. Check your output file.')

        bg = etitle[0] + 2
        ed = empty_line[np.where(empty_line>bg)[0][0]]

        tens = self.df[0].loc[bg:ed-1].map(lambda x: x.strip().split()).tolist()
        # print(self.df[0].loc[bg:ed-1])
        ndimen = len(tens)
        self.tensor = np.zeros([ndimen,ndimen], dtype=float)
        for i in range(ndimen):
            self.tensor[i, i:] = np.array(tens[i][1:-1], dtype=float)
        # Symmetrize tensor
        for i in range(ndimen):
            for j in range(i,ndimen):
                self.tensor[j][i] = self.tensor[i][j]

        # Unit conversion
        title = self.df[0][etitle[0]].lower()
        if re.search('gpa', title):
            pass
        elif re.search('hartree', title):
            if ndimen == 1:
                length = units.angstrom_to_au(self.get_lattice(initial=False)[0, 0]) * 1e-10 # AA --> m
                self.tensor = units.au_to_GPa(self.tensor) * (units.au_to_angstrom(1.)*1e-10)**3 / length # Eh --> GJ/m
            elif ndimen == 3:
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
            if ndimen == 1:
                if len(thickness) == 1:
                    self.tensor = self.tensor / (thickness[0]*1e-10)**2
                else:
                    self.tensor = self.tensor / (thickness[0]*1e-10) / (thickness[1]*1e-10)
            elif ndimen == 3:
                self.tensor = self.tensor / (thickness[0]*1e-10)
        else:
            if ndimen != 6:
                warngings.warn('Low dimensional materials without effective thickness! Output units have extra dimensions of length.',
                               stacklevel=2)
        return self.tensor


    def get_anh_spectra(self):
        """
        .. note::

            **This is not for the released feature of CRYSTAL23 v1.0.1**

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
        import re, warnings
        import numpy as np

        warnings.warn('This is not a released feature of CRYSTAL23 v1.0.1, make sure that you know what you are doing.',
                      stacklevel=2)

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

#------------------------------Deprecated-------------------------------------#
    @property
    def n_atoms(self):
        """
        Deprecated. Get structure object by 'get_geometry' and call the
        'num_sites' attribute
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry' and call the 'natoms' attribute of output.",
                      stacklevel=2)
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.num_sites

    @property
    def atom_symbols(self):
        """
        Deprecated. Get structure object by 'get_geometry' and call the
        'species_symbol' attribute
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry' and call the 'species_symbol' attribute of output.",
                      stacklevel=2)
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.species_symbol

    @property
    def atom_numbers(self):
        """
        Deprecated. Get structure object by 'get_geometry' and call the
        'species_Z' attribute
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry' and call the 'species_Z' attribute of output.",
                      stacklevel=2)
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.species_Z

    @property
    def atom_positions(self):
        """
        Deprecated. Get structure object by 'get_geometry' and call the
        'crys_coords' attribute
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry' and call the 'crys_coords' attribute of output.",
                      stacklevel=2)
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.crys_coords

    @property
    def atom_positions_cart(self):
        """
        Deprecated. Get structure object by 'get_geometry' and call the
        'cart_coords' attribute
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry' and call the 'cart_coords' attribute of output.",
                      stacklevel=2)
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.cart_coords

    @property
    def atom_positions_frac(self):
        """
        Deprecated. Get structure object by 'get_geometry' and call the
        'crys_coords' attribute
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry' and call the 'crys_coords' attribute of output.",
                      stacklevel=2)
        struc = self.get_geometry(initial=True, write_gui=False)
        return struc.crys_coords

    def get_last_geom(self, write_gui_file=True, symm_info='pymatgen'):
        """
        Deprecated. Use ``get_geometry(initial=False)``. 
        """
        import warnings

        warnings.warn("You are calling a deprecated property. Use 'get_geometry(initial=False)'.",
                      stacklevel=2)

        struc = self.get_geometry(initial=False, write_gui=write_gui_file, symm_info=symm_info)
        self.last_geom = [struc.lattice.matrix.tolist(), struc.species_Z, struc.cart_coords.tolist()]
        return self.last_geom


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

    @classmethod
    def read_file(cls, source):
        """
        Instantiate the object from a file.

        Args:
            source (str): The name of the input file.
        Returns:
            cls (Properties_input)
        """
        return cls(source)

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

    def make_band_block(self, k_path, n_kpoints, first_band, last_band,
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
                    struc = output.get_geometry(initial=False, write_gui=False)
                    for prje in projections:
                        index = [i+1 for i, ele in enumerate(struc.species_symbol) if prje.upper() == ele.upper()]
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


#-------------------------------obsolete methods-------------------------------#
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

    def make_bands_block(self, k_path, n_kpoints, first_band, last_band,
                         print_eig=0, print_option=1, precision=5,
                         title='BAND STRUCTURE CALCULATION'):
        """
        Deprecated. Use ``self.make_band_block()``.
        """
        import warnings

        warnings.warn('Deprecated. Use make_band_block() method instead.', stacklevel=2)
        return self.make_band_block(k_path, n_kpoints, first_band, last_band, print_eig,
                                    print_option, precision, title)

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
        import os

        if np.all(properties_output!=None):
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

    @classmethod
    def read_file(cls, properties_output):
        """
        Parse the properties output file.

        Args:
            properties_output (str): The properties output file.
        Returns:
            cls (Properties_output)
        """
        return cls(properties_output=properties_output)

#----------------------------electronic structure------------------------------#

    def read_electron_band(self, band_file):
        """
        Generate bands object from CRYSTAL BAND.DAT / fort.25 file. Energy
        unit: eV. E Fermi is aligned to 0.

        Args:
            band_file (str): Name of BAND.DAT / fort.25 file.

        Returns:
            self.bands (ElectronBand): The :ref:`electronics.ElectronBand <ref-ElectronBand>` object.
        """
        from CRYSTALpytools.base.extfmt import CrgraParser, DLVParser
        from CRYSTALpytools.units import H_to_eV, angstrom_to_au
        from CRYSTALpytools.electronics import ElectronBand
        import warnings

        file = open(band_file)
        flag = file.readline()
        file.close()
        if '-%-' in flag: # fort.25
            bandout = CrgraParser.band(band_file)
        elif '#' in flag: # BAND.DAT
            bandout = DLVParser.band(band_file)
        else:
            raise Exception("Pattern not found in '{}'. Is it a band structure file?".format(band_file))

        if not hasattr(self, 'file_name'):
            warnings.warn('Properties output file not found: 3D k path not available.',
                          stacklevel=2)
            struc = None; t3d = None; k3d = None
        else:
            struc = super().get_geometry()
            t3d, k3d = super().get_3dkcoord()
        self.bands = ElectronBand(
            bandout[0], bandout[1], bandout[2], bandout[3], bandout[4],
            bandout[5], struc, None, t3d, k3d, bandout[6])
        return self.bands


    def read_Fermi_surface(self, f35_file):
        """
        Generate the Fermi surface, i.e., band energy across the first brillouin
        zone :math:`E(k)` from CRYSTAL fort.35 file (with 'DLV_BAND' keyword).
        Energy unit: eV.

        .. note::

            When output is available, energies are aligned to :math:`E_{F}=0`.
            Otherwise values reported in fort.35 file are used.

        Args:
            f35_file (str): Name of fort.35 file.

        Returns:
            self.FermiSurf (FermiSurface): The :ref:`electronics.FermiSurface <ref-FermiSurface>` object
        """
        from CRYSTALpytools.base.extfmt import DLVParser
        from CRYSTALpytools.units import H_to_eV, angstrom_to_au
        from CRYSTALpytools.electronics import FermiSurface
        import warnings

        rlatt, band, _ = DLVParser.fort35(f35_file)
        rlatt = angstrom_to_au(rlatt) # A^-1 to Bohr^-1
        band = H_to_eV(band)
        if not hasattr(self, 'file_name'):
            warnings.warn('Properties output file not found: Fermi energy not available.',
                          stacklevel=2)
            efermi = 0.
        else:
            efermi = super().get_Fermi()
        self.FermiSurf = FermiSurface(rlatt, band-efermi, efermi=efermi, unit='eV')
        return self.FermiSurf


    def read_electron_dos(self, dos_file):
        """
        Get density of states from CRYSTAL DOSS.DAT or fort.25 file. Energy
        unit: eV. E Fermi is aligned to 0.

        Args:
            dos_file (str): Name of DOSS.DAT or fort.25 file
        Returns:
            self.doss (ElectronDOS): The :ref:`electronics.ElectronDOS <ref-ElectronDOS>` object.
        """
        from CRYSTALpytools.base.extfmt import CrgraParser, DLVParser
        from CRYSTALpytools.electronics import ElectronDOS

        file = open(dos_file)
        flag = file.readline()
        file.close()
        if '-%-' in flag:  # fort.25 file format
            dosout = CrgraParser.dos(dos_file)
        elif '#' in flag: # DOSS.DAT
            dosout = DLVParser.dos(dos_file)
        else:
            raise Exception("Pattern not found in '{}'. Is it a DOSS file?".format(dos_file))

        self.doss = ElectronDOS(spin=dosout[0], efermi=dosout[1], doss=dosout[2],
                   energy=dosout[3], unit=dosout[4])
        return self.doss

#-----------------------------2D scalar field----------------------------------#

    def read_topond(self, topondfile, type='infer'):
        """
        Read the 2D scalar plot files ('SURF\*.DAT') or trajectory files
        (TRAJ\*.DAT) written by `TOPOND <https://www.crystal.unito.it/topond.html>`_.

        Geometry information is printed in the standard ouput, which is not
        mandatory for 'SURF\*.DAT' but is mandatory for 'TRAJ\*.DAT'

        .. note::

            For the convenience of analysis and plotting, it is important to select
            the correct type for your input file. By default `type='infer'` will
            search for (case insensitive) the following strings:

            'SURFRHOO', 'SURFSPDE', 'SURFLAPP', 'SURFLAPM', 'SURFGRHO',
            'SURFKKIN', 'SURFGKIN', 'SURFVIRI', 'SURFELFB', 'TRAJGRAD',
            'TRAJMOLG'.

            For their meanings, please refer `TOPOND manual <https://www.crystal.unito.it/include/manuals/topond.pdf>`_.

        Args:
            topondfile (str): TOPOND formatted 2D plot file
            type (str): 'infer' or specified. Otherwise warning will be given.

        Returns:
            self.\* (ChargeDensity|SpinDensity|Gradient|Laplacian|HamiltonianKE|LagrangianKE|VirialField|ELF|GradientTraj|ChemicalGraph):
                Return to ``topond`` property classes, depending on input file types. The
                attribute name is upper case type names. If unknown, return to
                ``self.TOPOND``. A ``topond.ChargeDensity`` or
                ``topond.GradientTraj`` class is generated.
        """
        import warnings
        from CRYSTALpytools.base.extfmt import TOPONDParser
        from CRYSTALpytools.topond import \
            ChargeDensity, SpinDensity, Gradient, Laplacian, HamiltonianKE, \
            LagrangianKE, VirialField, ELF, GradientTraj, ChemicalGraph

        surflist = ['SURFRHOO', 'SURFSPDE', 'SURFLAPP', 'SURFLAPM', 'SURFGRHO',
                    'SURFKKIN', 'SURFGKIN', 'SURFVIRI', 'SURFELFB']
        trajlist = ['TRAJGRAD', 'TRAJMOLG']

        if type.lower() == 'infer':
            type = topondfile.upper()
        else:
            type = type.upper()

        issurf = False; istraj = False
        for t in surflist:
            if t in type: issurf = True; type = t; break
        for t in trajlist:
            if t in type: istraj = True; type = t; break

        # still need to distinguish surf and traj
        if issurf==False and istraj==False:
            warnings.warn("Unknown type string / filename does not contian type string.",
                          stacklevel=2)
            type = 'unknown'
            file = open(topondfile, 'r')
            header = file.readline()
            file.close()
            if 'DSAA' in header: issurf = True
            else: istraj = True

        if issurf == True:
            _, a, b, c, _, _, map, unit = TOPONDParser.contour2D(topondfile)
            if not hasattr(self, 'file_name'):
                warnings.warn('Properties output file not found: Geometry not available',
                              stacklevel=2)
                struc = None
                # The a, b, c by are dummy base vectors. Info in 3D space lost.
                base = np.vstack([a, b, c])
            else:
                struc = super().get_geometry()
                # no atom plot currently, though read method is given
                _, base = super().get_topond_geometry()
            # class instantiation
            if type == 'SURFRHOO' or type == 'unknown':
            # map from base method has spin dimension
                obj = ChargeDensity(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFSPDE':
                obj = SpinDensity(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFLAPP':
                obj = Laplacian(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFLAPM':
                obj = Laplacian(-map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFGRHO':
                obj = Gradient(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFKKIN':
                obj = HamiltonianKE(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFGKIN':
                obj = LagrangianKE(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFVIRI':
                obj = VirialField(map[:,:,0], base, 2, struc, unit)
            elif type == 'SURFELFB':
                obj = ELF(map[:,:,0], base, 2, struc, unit)
            obj._set_unit('Angstrom')

        elif istraj == True:
            if not hasattr(self, 'file_name'):
                raise Exception("Properties output file is mandatory for 'TRAJ' files.")
            wtraj, traj, unit = TOPONDParser.traj(topondfile)
            struc = super().get_geometry()
            # no atom plot currently, though read method is given
            _, base = super().get_topond_geometry()
            # class instantiation
            if type == 'TRAJGRAD' or type == 'unknown':
                obj = GradientTraj(wtraj, traj, base, struc, unit)
            elif type == 'TRAJMOLG':
                obj = ChemicalGraph(wtraj, traj, base, struc, unit)
            obj._set_unit('Angstrom')

        if type == 'unknown': type = 'TOPOND'

        setattr(self, type, obj)
        return getattr(self, type)

    def read_ECHG(self, *f25_files, method='normal', index=None):
        """
        Read charge / spin density data from a file. Unit: :math:`e.\\AA^{-3}`.

        Available methods are:

        * 'normal': Normal 1-file reading.  
        * 'subtract': Subtracting data from the first entry based on following
            entries. Multiple entries or 1 entry with 'PATO' keyword enabled.
            For multiple entries, make sure the charge map is in the first (and
            ideally the only) 'MAPN' data block, otherwise the code won't get
            the correct data. For 1 entry with 'PATO', data will be subtracted
            from the 'normal' system.  
        * 'alpha_beta': Save spin-polarized data in :math:`\\alpha` /
            :math:`\\beta` states, rather than charge(:math:`\\alpha+\\beta`)
            / spin(:math:`\\alpha-\\beta`). Single entry only.

        .. note::

            The standard screen output is highly recommended to add, which
            identifies the indices of corresponding 2D data maps. Otherwise the
            ``index`` input can be specified. Otherwise, the code only reads
            the first 2D data map for spin =1 and first 2 maps for spin=2.

        Args:
            \*f25_files (str): Path to the fort.25 file(s).
            method (str): Data processing method. See above.
            index (int|list): Sequence number of headers with the '-%-MAPN'
                pattern. Useful only if the standard screen output is not
                available. Starting from 0.
        Returns:
            self.echg (ChargeDensity): ``electronics.ChargeDensity`` object.
        """
        from CRYSTALpytools.base.extfmt import CrgraParser
        from CRYSTALpytools.electronics import ChargeDensity
        import numpy as np
        import pandas as pd
        import warnings

        method = method.lower()
        if method == 'substract': method = 'subtract'# an old typo
        if method != 'subtract' and method != 'alpha_beta' and method != 'normal':
            raise ValueError("Unknown method: '{}'.".format(method))

        pato = [] # used for method check
        if not hasattr(self, 'file_name'):
            if np.all(index==None):
                warnings.warn('Properties output file not found: Only the first 1 (2) density map(s) will be read for spin=1(2).',
                              stacklevel=2)
                index = None
            else:
                index = np.array(index, dtype=int, ndmin=1)
        else:
            df = pd.DataFrame(self.data)
            headers = df[df[0].str.contains(r'^\s*-%-[0-4]MAPN')].index.to_numpy(dtype=int)
            chg = df[df[0].str.contains(r'^\s*ALPHA\+BETA ELECTRONS DENSITY')].index.to_numpy(dtype=int)
            spin = df[df[0].str.contains(r'^\s*ALPHA\-BETA ELECTRONS DENSITY')].index.to_numpy(dtype=int)
            pato = df[df[0].str.contains(r'^\s*DENSITY MATRIX AS SUPERPOSITION')].index.to_numpy(dtype=int)

            if len(chg) == 0:
                raise Exception('Charge density calculation not found in the output file.')
            elif len(chg) == 1:
                if len(spin) == 1: # find the header closest to keywords
                    index = np.array([np.where(headers>chg[0])[0][0],
                                      np.where(headers>spin[0])[0][0]], dtype=int)
                    index = np.sort(index)
                else:
                    index = np.where(headers>chg[0])[0][0]
            else:
                if len(pato) == 1 and len(chg) == 2:
                    if chg[0] < pato[0] and chg[1] > pato[0]: use_idx = 0
                    else: use_idx = 1

                    if method != 'subtract': # Normal read of charge densities (No PATO)
                        if len(spin) != 0:
                            index = np.array([np.where(headers>chg[use_idx])[0][0],
                                              np.where(headers>spin[use_idx])[0][0]], dtype=int)
                            index = np.sort(index)
                        else:
                            index = np.where(headers>chg[0])[0][0]
                    else: # spin dimension is not kept
                        index = np.array([np.where(headers>chg[use_idx])[0][0],
                                          np.where(headers>chg[1-use_idx])[0][0]], dtype=int)
                        index = np.sort(index)
                else:
                    warnings.warn('Multiple charge densities exist in the calculation. Only the first density map will be read.',
                                  stacklevel=2)
                    index = 0

        # read file 0
        spin, a, b, c, cosxy, struc, map, unit = CrgraParser.mapn(f25_files[0], index)
        # methods
        if len(f25_files) == 1 and len(pato) == 1 and method == 'subtract': # PATO in the same file
            self.echg = ChargeDensity(map[use_idx], np.vstack([a[0],b[0],c[0]]), spin, 2, struc[0], unit)
            obj = ChargeDensity(map[1-use_idx], np.vstack([a[1],b[1],c[1]]), spin, 2, struc[1], unit)
            self.echg = self.echg.subtract(obj)
            self.echg._set_unit('Angstrom')
            self.echg.data = self.echg.data[::-1] # base vector use BA rather than AB
        else: # others
            if spin == 1:
                self.echg = ChargeDensity(map, np.vstack([a,b,c]), 1, 2, struc, unit)
            else:
                if isinstance(cosxy, float):
                    raise Exception('Broken file: charge / spin density missing for spin polarized systems.')
                self.echg = ChargeDensity(np.dstack([map[0], map[1]]), np.vstack([a[0],b[0],c[0]]), spin, 2, struc[0], unit)
            self.echg._set_unit('Angstrom')
            self.echg.data = self.echg.data[::-1] # base vector use BA rather than AB
            if method == 'alpha_beta':
                if len(f25_files) > 1:
                    warnings.warn("The 'alpha_beta' method is used only for a single entry. Nothing is done to other entries.",
                                  stacklevel=2)
                elif spin != 2:
                    warnings.warn("Not a spin-polarized system, do nothing", stacklevel=2)
                else:
                    self.echg.alpha_beta()
            elif method == 'subtract':
                if len(f25_files) > 1:
                    self.echg = self.echg.subtract(*[f for f in f25_files[1:]])
                else:
                    warnings.warn("Nothing to subtract.", stacklevel=2)

        return self.echg

    def read_ECH3(self, *cubefiles, method='normal'):
        """
        Read 3D charge / spin density data from CUBE files. Unit: :math:`e.\\AA^{-3}`.

        .. note::

            Only compatible with CRYSTAL cube outputs. Lattice constants are
            annotated in the comment line.

        Available methods are:

        * 'normal': Normal reading, 1 or 2 entries for charge and spin
            densities.  
        * 'subtract': Subtracting data from the first entry based on following
            entries.  
        * 'alpha_beta': Save spin-polarized data in :math:`\\alpha` /
            :math:`\\beta` states, rather than charge(:math:`\\alpha+\\beta`)
            / spin(:math:`\\alpha-\\beta`). 1 or 2 entries for charge and spin
            densities.

        Args:
            \*cubefiles (str): Path to the CUBE file(s).
            method (str): Data processing method. See above.

        Returns:
            self.ech3 (ChargeDensity): ``electronics.ChargeDensity`` object.
        """
        from CRYSTALpytools.base.extfmt import CUBEParser
        from CRYSTALpytools.electronics import ChargeDensity
        import numpy as np
        import pandas as pd
        import warnings

        method = method.lower()
        if method == 'substract': method = 'subtract'# an old typo
        if method != 'subtract' and method != 'alpha_beta' and method != 'normal':
            raise ValueError("Unknown method: '{}'.".format(method))
        if len(cubefiles) > 2 and method != 'subtract':
            raise ValueError("Only 1 or 2 entries are permitted for method: '{}'.".format(method))
        if (method=='subtract' or method=='alpha_beta') and len(cubefiles) < 2:
            warings.warn("At least 2 files are needed for the specified method. Using 'normal' now.",
                         stacklevel=2)
        # The first entry
        origin, a, b, c, struc, data, _ = CUBEParser.read_cube(cubefiles[0])
        # Structure from output if provided, to keep periodicity settings.
        ## compare structures of CUBE and output
        def compare_struc(struc0, struc1):
            if np.linalg.norm(struc0.lattice.matrix - struc1.lattice.matrix) > 1e-4:
                return False
            if struc0.num_sites != struc1.num_sites:
                return False
            if np.linalg.norm(struc0.frac_coords%1-struc1.frac_coords%1)>1e-2:
                return False
            return True
        if hasattr(self, 'file_name'):
            struc1 = super().get_geometry()
            if compare_struc(struc, struc1) == False:
                warnings.warn('Inconsistent geometries are given in output and CUBE files, using the one from output.',
                              stacklevel=2)
            struc = struc1

        # Other entries
        if len(cubefiles) > 1:
            for f in cubefiles[1:]:
                o1, a1, b1, c1, struc1, data1, _ = CUBEParser.read_cube(f)
                if np.linalg.norm(origin-o1)>1e-4 or np.linalg.norm(a-a1)>1e-4 \
                or np.linalg.norm(b-b1)>1e-4 or np.linalg.norm(c-c1)>1e-4 \
                or np.linalg.norm(np.array(data1.shape)-np.array(data.shape))>1e-4:
                    raise Exception("Inconsistent data grid between the initial and the file: '{}'.".format(f))
                if method == 'subtract':
                    data -= data1
                else:
                    if compare_struc(struc, struc1) == False and not hasattr(self, 'file_name'):
                        # only raise when no reliable structure from output is available
                        raise Exception("Inconsistent structure between the initial and the file: '{}'.".format(f))

        if method == 'subtract':
            self.ech3 = ChargeDensity(np.expand_dims(data, axis=3),
                                      [origin, a, b, c], 1, 3, struc=struc, unit='a.u.')
            del data, data1
        else:
            spin = len(cubefiles)
            datanew = np.zeros([data.shape[0], data.shape[1], data.shape[2], spin])
            datanew[:, :, :, 0] = data; del data
            if spin > 1:
                datanew[:, :, :, 1] = data1; del data1
            self.ech3 = ChargeDensity(datanew, [origin, a, b, c], spin, 3,
                                      struc=struc, unit='a.u.')
            del datanew

        if method == 'alpha_beta':
            self.ech3.alpha_beta()
        self.ech3._set_unit('Angstrom')
        return self.ech3

#-----------------------------2D vector field----------------------------------#

    def read_relativistics(self, f25_file, type, index=None):
        """
        Read 2D scalar / vector fields from 2c-SCF calculations, generated by
        the 'PROPS2COMP' keyword.

        .. note ::

            The standard screen output is highly recommended to add, which
            identifies the indices of corresponding 2D data maps. Otherwise the
            ``index`` must be specified by integer or list of integers.

        Args:
            f25_file (str): File name of the fort.25 file.
            type (str|int|list): Property to calculate, 'DENSITY', 'MAGNETIZ',
                'ORBCURDENS', or 'SPICURDENS'.
            index (int|list): Sequence number of headers with the '-%-MAPN'
                pattern. Useful only if the standard screen output is not
                available. Starting from 0.

        Returns:
            self.\* (ChargeDensity|Magnetization|OrbitalCurrentDensity|SpinCurrentDensity):
                Return to classes defined in the ``relativisitcs`` module,
                depending on ``type``. The attribute name is upper case type
                names. Unit: charge densities, :math:`\\AA^{-3}`; magnetization,
                A/m; Orbital/spin densities, A/m :math:`^{2}`.
        """
        import numpy as np
        import pandas as pd
        from CRYSTALpytools.base.extfmt import CrgraParser
        from CRYSTALpytools.relativistics import (ChargeDensity, Magnetization,
            OrbitalCurrentDensity, SpinCurrentDensity)

        type_avail = ['DENSITY', 'MAGNETIZ', 'ORBCURDENS', 'SPICURDENS']
        type = type.upper()
        if type not in type_avail:
            raise ValueError("Unknown type: '{}'.".format(type))

        if not hasattr(self, 'file_name'):
            if np.all(index==None):
                raise Exception("Properties output file is mandatory here, otherwise specify the index.")
            index = np.array(index, dtype=int, ndmin=1)
        else:
            struc = super().get_geometry()
            # get corresponding MPNET entries from output file
            df = pd.DataFrame(self.data)
            headers = df[df[0].str.contains(r'^\s*-%-[0-4]MAPN')].index.to_numpy(dtype=int)
            dens = df[df[0].str.contains(r'^\s+PARTICLE\-NUMBER DENSITY')].index.to_numpy(dtype=int)
            mag = df[df[0].str.contains(r'^\s+[X,Y,Z]\-COMP MAGNETIZATION')].index.to_numpy(dtype=int)
            orbt = df[df[0].str.contains(r'^\s+[X,Y,Z]\-COMP ORB\-CURR DENS')].index.to_numpy(dtype=int)
            spinc = df[df[0].str.contains(r'^\s+[X,Y,Z]\-COMP [X,Y,Z] S\-CURR DENS')].index.to_numpy(dtype=int)

            indices = np.concatenate([dens, mag, orbt, spinc])
            indices = np.sort(indices)
            if len(indices) == 0:
                raise Exception('2c-SCF calculation not found in the output file.')
            headers_2c = headers[np.where(headers>indices[-1])] # headers for 2c scf block
            if type == 'DENSITY':
                if len(dens) != 1:
                    raise Exception('Charge density not found, or found more than once, in the calculation.')
                # index in headers_2c
                idx = np.where(indices==dens[0])[0]
            elif type == 'MAGNETIZ':
                if len(mag) != 3:
                    raise Exception('Magnetization not found, or found more than once, in the calculation.')
                 # index in headers_2c
                idx = np.array([np.where(indices==mag[0])[0],
                                np.where(indices==mag[1])[0],
                                np.where(indices==mag[2])[0]], dtype=int)
            elif type == 'ORBCURDENS':
                if len(orbt) != 3:
                    raise Exception('Orbital-current density not found, or found more than once, in the calculation.')
                # index in headers_2c
                idx = np.array([np.where(indices==orbt[0])[0],
                                np.where(indices==orbt[1])[0],
                                np.where(indices==orbt[2])[0]], dtype=int)
            elif type == 'SPICURDENS':
                if len(spinc) != 9:
                    raise Exception('Spin-current density not found, or found more than once, in the calculation.')
                # index in headers_2c
                idx = np.concatenate([
                    np.where(indices==spinc[0])[0], np.where(indices==spinc[1])[0],
                    np.where(indices==spinc[2])[0], np.where(indices==spinc[3])[0],
                    np.where(indices==spinc[4])[0], np.where(indices==spinc[5])[0],
                    np.where(indices==spinc[6])[0], np.where(indices==spinc[7])[0],
                    np.where(indices==spinc[8])[0]
                ], dtype=int)
            # index in the full fort.25 file.
            index = [np.where(headers==headers_2c[i])[0][0] for i in idx]

        # read file
        spin, a, b, c, cosxy, struc, map, unit = CrgraParser.mapn(f25_file, index)

        if type == 'DENSITY':
            obj = ChargeDensity(map, np.vstack([a,b,c]), 1, 2, struc, unit)
            obj.data = obj.data[::-1] # base vector use BA rather than AB
            obj._set_unit('Angstrom')
        elif type == 'MAGNETIZ':
            obj = Magnetization(np.dstack([map[0], map[1], map[2]]),
                                np.vstack([a[0],b[0],c[0]]), 2, struc, unit)
            obj.data = obj.data[::-1] # base vector use BA rather than AB
            obj._set_unit('SI')
        elif type == 'ORBCURDENS':
            obj = OrbitalCurrentDensity(np.dstack([map[0], map[1], map[2]]),
                                        np.vstack([a[0],b[0],c[0]]), 2, struc, unit)
            obj.data = obj.data[::-1] # base vector use BA rather than AB
            obj._set_unit('SI')
        elif type == 'SPICURDENS':
            obj = SpinCurrentDensity(np.dstack([map[0], map[1], map[2]]),
                                     np.dstack([map[3], map[4], map[5]]),
                                     np.dstack([map[6], map[7], map[8]]),
                                     np.vstack([a[0],b[0],c[0]]), 2, struc, unit)
            obj.data_x = obj.data_x[::-1] # base vector use BA rather than AB
            obj.data_y = obj.data_y[::-1] # base vector use BA rather than AB
            obj.data_z = obj.data_z[::-1] # base vector use BA rather than AB
            obj._set_unit('SI')

        setattr(self, type, obj)
        return getattr(self, type)

#-------------------------------------XRD--------------------------------------#

    def read_XRDspec(self, option='LP'):
        """
        Read XRD spectra from standard screen output of properties calculation.
        It is envoked by keyword 'XRDSPEC'.

        Args:
            option (str): 'NC' for no correction (The 'INTENS' col); 'LP' for
                Lorentz and polarization effects ('INTENS-LP') and 'DW' for LP
                with Debye-Waller thermal factors ('INTENS-LP-DW').
        Returns:
            self.XRDspec (XRD): The ``spectra.XRD`` object with spectra
                information.
        """
        import numpy as np
        from CRYSTALpytools.spectra import XRD

        if not hasattr(self, 'file_name'):
            raise Exception("The screen output file is not defined. Use the 'from_file()' method to add it.")

        spec = super().get_XRDSPEC()

        option = option.upper()
        if option == 'NC':
            self.XRDspec = XRD(theta=spec[:, 0], spectra=spec[:, 1])
        elif option == 'LP':
            self.XRDspec = XRD(theta=spec[:, 0], spectra=spec[:, 2])
        elif option == 'DW':
            self.XRDspec = XRD(theta=spec[:, 0], spectra=spec[:, 3])
        else:
            raise ValueError("Unknown XRD spectra read option: '{}'.".format(option))
        return self.XRDspec

#--------------------------------1D line profile-------------------------------#

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

#-----------------------------transport properties-----------------------------#

    def read_transport(self, boltztra_out):
        """
        Read electron transport properties by the BOLTZTRA keyword, including
        'KAPPA', 'SIGMA', 'SIGMAS', 'SEEBECK' and 'TDF'. Though currently the
        geometry information is not required, it is saved if the standard
        output file is given.

        .. note::

            For 'SEEBECK', all the 9 elements of the tensor was printed. As far
            as the developers have been aware of, it is symmetrized. Therefore
            the redundant 'yx', 'zx' and 'zy' dimensions are removed to keep
            consistent with other outputs.

        Args:
            boltztra_out (str): 'DAT' files by CRYSTAL BOLTZTRA keyword.

        Returns:
            self.\* (Tensor|Distribution): ``transport.Tensor`` ('KAPPA',
                'SIGMA', 'SIGMAS', 'SEEBECK') or ``transport.Distribution``
                (TDF) classes, depending on the input file. The attribute name
                is upper case types.
        """
        from CRYSTALpytools.base.extfmt import BOLTZTRAParaser
        from CRYSTALpytools.transport import Kappa, Sigma, Seebeck, SigmaS, TDF

        if hasattr(self, 'file_name'):
            struc = super().get_geometry()
        else:
            struc = None

        file = open(boltztra_out)
        header = file.readline()
        file.close()
        if 'Transport distribution function' in header:
            out = BOLTZTRAParaser.distribution(boltztra_out)
            if out[1] == 'TDF':
                obj = TDF(out[2], out[3], struc)
            else:
                raise TypeError("Unknown distribution function. Please contact the developers for updates.")
        else:
            out = BOLTZTRAParaser.tensor(boltztra_out)
            if out[1] == 'KAPPA':
                obj = Kappa(out[3], out[4], out[5], out[6], struc)
            elif out[1] == 'SIGMA':
                obj = Sigma(out[3], out[4], out[5], out[6], struc)
            elif out[1] == 'SIGMAS':
                obj = SigmaS(out[3], out[4], out[5], out[6], struc)
            elif out[1] == 'SEEBECK':
                obj = Seebeck(out[3], out[4], out[5], out[6], struc)

        setattr(self, out[1], obj)
        return getattr(self, out[1])

#------------------------------------------------------------------------------#
#--------------------------------obsolete methods------------------------------#
#------------------------------------------------------------------------------#
    def read_cry_band(self, band_file):
        """
        Deprecated. Use ``read_electron_band``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_electron_band' instead.",
                      stacklevel=2)
        return self.read_electron_band(band_file)

    def read_cry_doss(self, dos_file):
        """
        Deprecated. Use ``read_electron_dos``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_electron_dos' instead.",
                      stacklevel=2)
        return self.read_electron_dos(dos_file)

    def read_cry_ECHG(self, f25_file):
        """
        Deprecated. Use ``read_ECHG``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_ECHG' instead.",
                      stacklevel=2)
        return self.read_ECHG(f25_file, method='normal')

    def read_cry_ECHG_delta(self, f25_file1, f25_file2):
        """
        Deprecated. Use ``read_ECHG``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_ECHG' instead.",
                      stacklevel=2)
        return self.read_ECHG(f25_file1, f25_file2, method='subtract')

    def read_cry_contour(self, properties_output):
        """
        Deprecated. Use ``read_topond``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_topond' instead.",
                      stacklevel=2)
        return self.read_topond(properties_output)

    def read_cry_seebeck(self, properties_output):
        """
        Deprecated. Use ``read_transport``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_transport' instead.",
                      stacklevel=2)
        obj = self.read_transport(properties_output)
        if obj.type != 'SEEBECK':
            raise Exception('Input is not a SEBECK coefficient file.')
        return obj

    def read_cry_sigma(self, properties_output):
        """
        Deprecated. Use ``read_transport``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_transport' instead.",
                      stacklevel=2)
        obj = self.read_transport(properties_output)
        if obj.type != 'SIGMA':
            raise Exception('Input is not a conductivity file.')
        return obj

    def read_vecfield(self, properties_output, which_prop):
        """
        Deprecated. Use ``read_relativistics``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_relativistics' instead.",
                      stacklevel=2)
        index = [0]
        for i in which_prop:
            if 'm' in which_prop:
                index.extend([index[-1]+1, index[-1]+2, index[-1]+3])
                type = 'MAGNETIZ'
            if 'j' in which_prop:
                index.extend([index[-1]+1, index[-1]+2, index[-1]+3])
                type = 'ORBCURDENS'
            if 'J' in which_prop:
                index.extend([index[-1]+1, index[-1]+2, index[-1]+3,
                              index[-1]+4, index[-1]+5, index[-1]+6,
                              index[-1]+7, index[-1]+8, index[-1]+9])
                type = 'SPICURDENS'
        return self.read_relativistics(properties_output, type, index[1:])

    def read_cry_xrd_spec(self, properties_output):
        """
        Deprecated. Use ``read_XRDspec``.
        """
        import warnings

        warnings.warn("You are calling a deprecated function. Use 'read_XRDspec' instead.",
                      stacklevel=2)
        if not hasattr(self, 'file_name'): self = Properties_output(properties_output)
        self.read_XRDspec()
        return self.XRDspec



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

