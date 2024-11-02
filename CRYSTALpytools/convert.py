#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions that do conversion between data / file formats
"""

def cry_ase2gui(structure, gui_file=None, vacuum=None, symmetry=True):
    """
    Transform an ASE Structure object into a Pymatgen structure object and then
    a CRYSTAL structure (gui) object. Vacuum layer is set to 500 Angstrom
    as the default of CRYSTAL for low symmstry systems

    Args:
        structure (ASE Structure): ASE Structure object.
        gui_file (str): CRYSTAL gui / fort.34 file.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            length of non-periodic direction to 500 Angstrom. Low dimensional systems only.
        symmetry (bool): Perform symmetry analysis.
    """
    # First transform into pmg and then write the gui

    from pymatgen.io.ase import AseAtomsAdaptor
    from CRYSTALpytools.convert import cry_pmg2gui

    pmg_structure = AseAtomsAdaptor().get_structure(structure)
    gui = cry_pmg2gui(pmg_structure, gui_file, vacuum=vacuum, symmetry=symmetry)

    return gui


def cry_bands2pmg(band, output, labels=None):
    """
    Transform a CRYSTAL bands object/file into a Pymatgen ``BandStructureSymmLine``
    object (inherited from ``BandStructure``). No projection is available for now.
    Output is mandatory here.

    Args:
        band (str|ElectronBand): CRYSTAL fort.25 or BAND.DAT file or
            CRYSTALpytools ElectronBand object
        output (str): CRYSTAL properties output file
        labels (list[str]): K point labels to display in the band structure.

    Returns:
        BandStructureSymmLine: Pymatgen band structure object.
    """
    from CRYSTALpytools.electronics import ElectronBand
    from CRYSTALpytools.crystal_io import Properties_output

    # Generate the Band object
    if type(band) == str:
        band = ElectronBand.from_file(band, output)
    else:
        pout = Properties_output(output)
        band._set_unit('eV')
        band.geometry = pout.get_geometry()
        band.reciprocal_latt = pout.get_reciprocal_lattice()
        band.tick_pos3d, band.k_path3d = pout.get_3dkcoord()
    return band.to_pmg(labels)


def cry_gui2ase(gui, vacuum=None, **kwargs):
    """
    Transform a CRYSTAL structure (gui) file into an ASE atoms object.

    Args:
        gui (str|geometry.Crystal_gui): CRYSTAL gui / fort.34 file or CRYSTALpytools gui object
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of ASE atoms object. Low dimensional systems only.
        **kwargs: Passed to ASE Atoms constructor
    Returns:
        Atoms: ASE atoms object.
    """
    from CRYSTALpytools.convert import cry_gui2pmg
    from pymatgen.io.ase import AseAtomsAdaptor

    struc = cry_gui2pmg(gui, vacuum=vacuum)

    return AseAtomsAdaptor().get_atoms(struc, **kwargs)


def cry_gui2cif(gui, cif_file_name, vacuum=None, **kwargs):
    """
    Read a CRYSTAL structure (gui) file and save a cif file. The `CifWriter <https://pymatgen.org/pymatgen.io.html#pymatgen.io.cif.CifWriter>`_
    object of Pymatgen is called. By default, ``symprec = 0.01`` is used.

    Args:
        gui (str|geometry.Crystal_gui): CRYSTAL gui / fort.34 file or CRYSTALpytools gui object
        cif_file_name (str): Name (including path) of the cif file to be saved
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of Pymatgen atoms object. Low dimensional systems only.
        **kwargs: Passed to Pymatgen CifWriter.
    """
    from CRYSTALpytools.convert import cry_gui2pmg
    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(gui, vacuum=vacuum, molecule=False)
    if len(kwargs) == 0:
        CifWriter(structure, symprec=0.01, **kwargs).write_file(cif_file_name)
    else:
        CifWriter(structure, **kwargs).write_file(cif_file_name)

    return


def cry_gui2pmg(gui, vacuum=None, molecule=True):
    """
    Transform a CRYSTAL structure (gui) object into a Pymatgen Structure object.

    Args:
        gui (str|geometry.Crystal_gui): CRYSTAL gui / fort.34 file or CRYSTALpytools gui object
        vacuum (float): Vacuum distance. Unit: :math:`\\AA`. If none, set the
            ``pbc`` attribute of Pymatgen object. Low dimensional systems only.
        molecule (bool): Generate a Molecule Pymatgen object for 0D structures.

    Returns:
        Structure or Molecule: Pymatgen Structure or Molecule object.
    """
    from CRYSTALpytools.geometry import Crystal_gui
    from pymatgen.core.structure import Structure, Molecule
    from pymatgen.core.lattice import Lattice
    import numpy as np

    # Generate the gui object
    if type(gui) == str:
        gui = Crystal_gui().read_gui(gui)

    if gui.dimensionality == 0:
        if molecule == True:
            return Molecule(gui.atom_number, gui.atom_positions)

        elif molecule == False:
            if np.all(vacuum!=None):
                pbc = (True, True, True)
#                gui.lattice.setflags(write=1)
                thickness_x = np.amax(np.array(gui.atom_positions)[:, 0]) - \
                        np.amin(np.array(gui.atom_positions)[:, 0])
                thickness_y = np.amax(np.array(gui.atom_positions)[:, 1]) - \
                        np.amin(np.array(gui.atom_positions)[:, 1])
                thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                        np.amin(np.array(gui.atom_positions)[:, 2])

                gui.lattice[0][0] = thickness_x + vacuum
                gui.lattice[1][1] = thickness_y + vacuum
                gui.lattice[2][2] = thickness_z + vacuum
            else:
                pbc = (False, False, False)

    if gui.dimensionality == 1:
        if np.all(vacuum!=None):
            pbc = (True, True, True)
#            gui.lattice.setflags(write=1)
            thickness_y = np.amax(np.array(gui.atom_positions)[:, 1]) - \
                        np.amin(np.array(gui.atom_positions)[:, 1])
            thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                        np.amin(np.array(gui.atom_positions)[:, 2])

            gui.lattice[1][1] = thickness_y + vacuum
            gui.lattice[2][2] = thickness_z + vacuum
        else:
            pbc = (True, False, False)

    if gui.dimensionality == 2:
        if np.all(vacuum!=None):
            pbc = (True, True, True)
#            gui.lattice.setflags(write=1)
            thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                        np.amin(np.array(gui.atom_positions)[:, 2])

            gui.lattice[2][2] = thickness_z + vacuum
        else:
            pbc = (True, True, False)

    #Convert the pseudopotential / different basis set atoms to their original atomic number
    atomic_numbers = np.array(gui.atom_number)%100
    if gui.dimensionality == 3:
        pbc = (True, True, True)

    latt = Lattice(gui.lattice, pbc=pbc)
    return Structure(latt, gui.atom_number, gui.atom_positions, coords_are_cartesian=True)


def cry_gui2xyz(gui, xyz_file_name, **kwargs):
    """
    Transform a CRYSTAL structure (gui) file into an XYZ file.

    Args:
        xyz_file_name (str): Name of the XYZ file to be saved.
        gui (str): CRYSTAL gui / fort.34 file.
        **kwargs: Passed to Pymatgen XYZ object.
    """

    from pymatgen.io.xyz import XYZ
    from CRYSTALpytools.convert import cry_gui2pmg

    structure = cry_gui2pmg(gui, molecule=True) #this returns a pmg Molecule object
    XYZ(structure, **kwargs).write_file(xyz_file_name)


def cry_out2ase(output, vacuum=None, initial=False, **kwargs):
    """
    Transform a CRYSTAL output object into an ASE atoms object.

    Args:
        output (str|Crystal_output): Crystal output file or Crystal_output object
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of ASE atoms object. Low dimensional systems only.
        initial (bool): Read the last geometry of the output file.
        **kwargs: Passed to ASE Atoms constructor

    Returns:
        Atoms: ASE atoms object.
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from CRYSTALpytools.convert import cry_out2pmg

    struc = cry_out2pmg(output, vacuum=vacuum, initial=initial)

    return AseAtomsAdaptor().get_atoms(struc, **kwargs)


def cry_out2cif(output, cif_file_name, vacuum=None, initial=False, **kwargs):
    """
    Read geometry from a CRYSTAL output file and save it as a CIF file. The
    `CifWriter <https://pymatgen.org/pymatgen.io.html#pymatgen.io.cif.CifWriter>`_
    object of Pymatgen is called. By default, ``symprec = 0.01`` is used.

    Args:
        output (str|Crystal_output): Crystal output file or Crystal_output object
        cif_file_name (str): Name (including path) of the CIF file to be saved.
        vacuum (float): Vacuum distance. Unit: :math:`\\AA`. If none, set the
            ``pbc`` attribute of Pymatgen atoms object. Low dimensional systems only.
        initial (bool): Read the last geometry of the output file.
        **kwargs: Passed to Pymatgen CifWriter.
    """
    from CRYSTALpytools.convert import cry_out2pmg
    from pymatgen.io.cif import CifWriter

    structure = cry_out2pmg(output, vacuum=vacuum, initial=initial, molecule=False)
    if len(kwargs) == 0:
        CifWriter(structure, symprec=0.01, **kwargs).write_file(cif_file_name)
    else:
        CifWriter(structure, **kwargs).write_file(cif_file_name)


def cry_out2pmg(output, vacuum=None, initial=False, molecule=True):
    """
    Read geometry from a CRYSTAL output file and save it as a pymatgen structure
    object.

    Args:
        output (str|Crystal_output): Crystal output file or Crystal_output object
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of Pymatgen object. Low dimensional systems only.
        initial (bool): Read the last geometry of the output file.
        molecule (bool): Generate a Molecule Pymatgen object for 0D structures.

    Returns:
        Structure: Pymatgen Structure object.
    """
    from CRYSTALpytools.crystal_io import Crystal_output
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure, Molecule
    import numpy as np
    import copy

    #Extract information from the output file
    if type(output) == str:
        out = Crystal_output(output)
    else:
        out = output

    ndimen = out.get_dimensionality()
    struc = out.get_geometry(initial=initial, write_gui=False)
    if ndimen != 0:
        latt_mx = copy.deepcopy(struc.lattice.matrix)
    else:
        latt_mx = np.eye(3) * 500.

    if ndimen == 0:
        if molecule == True:
            struc = Molecule(species=struc.species, coords=struc.cart_coords)
            return struc
        else:
            if np.all(vacuum!=None):
                pbc = (True, True, True)
                thickness_x = np.amax(struc.cart_coords[:, 0]) - np.amin(struc.cart_coords[:, 0])
                thickness_y = np.amax(struc.cart_coords[:, 1]) - np.amin(struc.cart_coords[:, 1])
                thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])

                latt_mx[0, 0] = thickness_x + vacuum
                latt_mx[1, 1] = thickness_y + vacuum
                latt_mx[2, 2] = thickness_z + vacuum
            else:
                pbc = (False, False, False)

    if ndimen == 1:
        if np.all(vacuum!=None):
            pbc = (True, True, True)
            thickness_y = np.amax(struc.cart_coords[:, 1]) - np.amin(struc.cart_coords[:, 1])
            thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])

            latt_mx[1, 1] = thickness_y + vacuum
            latt_mx[2, 2] = thickness_z + vacuum
        else:
            pbc = (True, False, False)

    if ndimen == 2:
        if np.all(vacuum!=None):
            pbc = (True, True, True)
            thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])

            latt_mx[2, 2] = thickness_z + vacuum
        else:
            pbc = (True, True, False)

    if ndimen == 3:
        pbc = (True, True, True)

    latt = Lattice(latt_mx, pbc=pbc)
    return Structure(latt, struc.atomic_numbers, struc.cart_coords, coords_are_cartesian=True)


def cry_out2xyz(output, xyz_file_name, initial=False, **kwargs):
    """
    Read geometry from a CRYSTAL output file and save it as an XYZ file.

    Args:
        output (str|Crystal_output): Crystal output file or Crystal_output object
        xyz_file_name (str): Name (including path) of the XYZ file to be saved.
        initial (bool): Read the last geometry of the output file.
        **kwargs: Passed to Pymatgen XYZ object.
    """
    from pymatgen.io.xyz import XYZ
    from CRYSTALpytools.convert import cry_out2pmg

    structure = cry_out2pmg(output, initial=initial, molecule=True) #this returns a pmg Molecule object
    XYZ(structure, **kwargs).write_file(xyz_file_name)


def cry_pmg2gui(structure, gui_file=None, pbc=None, vacuum=None, symmetry=True,
                zconv=None, **kwargs):
    """
    Save a pymatgen Structure object into a CRYSTAL gui file. Vacuum layer is
    set to 500 Angstrom as the default of CRYSTAL for low symmstry systems.

    Args:
        structure (Structure | Molecule): Pymatgen Structure / Molecule object.
        gui_file (str): CRYSTAL gui / fort.34 file.
        pbc (list): 1\*3 boolian list. Implements periodicity along x, y and z
            directions. If none, the code will read it from input structure.
        vacuum (float): Vacuum distance. Unit: :math:`\\AA`. If none, set the
            length of non-periodic direction to 500 Angstrom. Low dimensional
            systems only.
        symmetry (bool): Do symmetry analysis.
        zconv (list[list[int, int]]): 1st element: The **index** of atom;
                2nd element: The new conventional atomic number.
        **kwargs: Passed to Pymatgen SpacegroupAnalyzer object. Valid only
            if ``symmetry=True``.
    """
    from CRYSTALpytools.geometry import Crystal_gui
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer,PointGroupAnalyzer
    from CRYSTALpytools.geometry import CStructure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import center_slab
    from pymatgen.core.structure import Structure, Molecule

    import numpy as np
    import warnings
    import copy

    # dimensionality
    if np.all(pbc==None):
        if isinstance(structure, Molecule):
            pbc = (False, False, False)
            structure = Molecule(structure.species,
                                 [i.coords for i in structure.sites])
        else:
            pbc = structure.pbc

    gui = Crystal_gui()
    dimensionality = pbc.count(True)

    if dimensionality == 0 and isinstance(structure, Molecule):
        molecule = structure
        is_molecule = True # 0D object called as molecule
    elif dimensionality > 0 or isinstance(structure, CStructure):
        is_molecule = False # >0D object called as structure
    elif dimensionality == 0 and not isinstance(structure, Molecule) and not isinstance(structure, CStructure):
        warnings.warn('Dimensionality is set to 0, but the structure is not a molecule. Periodicity will be removed.')
        molecule = Molecule(structure.species,
                            [i.coords for i in structure.sites])
        is_molecule = True # 0D object called as molecule
    elif dimensionality > 0 and isinstance(structure, Molecule):
        warnings.warn('Dimensionality is set to 1-3, but the structure is a molecule. Periodicity will be added.')
        structure = structure.get_boxed_structure(500., 500., 500.)
        is_molecule = False # >0D object called as structure

    gui.dimensionality = dimensionality

    if is_molecule == True: # 0D
        lattice_vectors = np.identity(3)*500.
        gui.lattice = lattice_vectors
        gui.n_atoms = molecule.num_sites
        if symmetry == True:
            symmops = PointGroupAnalyzer(structure).get_symmetry_operations()
            n_symmops = 0
            gui.symmops = []
            for symmop in symmops:

                if np.all(symmop.translation_vector == 0.):

                    n_symmops += 1
                    gui.symmops.extend(symmop.rotation_matrix.tolist())
                    gui.symmops.append(symmop.translation_vector.tolist())

            gui.n_symmops = n_symmops
            gui.space_group = 1
        else:
            gui.space_group = 1
            gui.symmops = []
            gui.n_symmops = 1
            gui.symmops.extend(np.identity(3).tolist())
            gui.symmops.append([0.0,0.0,0.0])
        gui.atom_number = list(molecule.atomic_numbers)
        gui.atom_positions = molecule.cart_coords.tolist()
    else: # 1-3D
        latt_mx = structure.lattice.matrix
        if dimensionality == 2:
            if pbc[0] == False: # X no periodicity
                warnings.warn('The non-periodic direction will be rotated to z axis.')
                rot = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
                lattice_vectors = copy.deepcopy(structure.lattice.matrix) @ rot
                lattice_vectors[2, 1] = latt_mx[0, 0] # Perpendicular elements are not on diagonal after rotation
                atom_coords = copy.deepcopy(structure.cart_coords) @ rot

            elif pbc[1] == False: # Y no periodicity
                warnings.warn('The non-periodic direction will be rotated to z axis.')
                rot = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
                lattice_vectors = copy.deepcopy(structure.lattice.matrix) @ rot
                lattice_vectors[2, 0] = latt_mx[1, 1] # Perpendicular elements are not on diagonal after rotation
                atom_coords = copy.deepcopy(structure.cart_coords) @ rot

            else: # Z no periodicity
                lattice_vectors = copy.deepcopy(structure.lattice.matrix)
                lattice_vectors[2, 2] = latt_mx[2, 2]
                atom_coords = copy.deepcopy(structure.cart_coords)

        elif gui.dimensionality == 1:
            if pbc[0] == True: # X periodic
                lattice_vectors = copy.deepcopy(structure.lattice.matrix)
                lattice_vectors[1, 1] = latt_mx[1, 1]
                lattice_vectors[2, 2] = latt_mx[2, 2]
                atom_coords = copy.deepcopy(structure.cart_coords)

            elif pbc[1] == True: # Y periodic
                warnings.warn('The periodic direction will be rotated to x axis.')
                rot = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
                lattice_vectors = copy.deepcopy(structure.lattice.matrix) @ rot
                lattice_vectors[1, 0] = latt_mx[2, 2] # Perpendicular elements are not on diagonal after rotation
                lattice_vectors[2, 1] = latt_mx[0, 0]
                atom_coords = copy.deepcopy(structure.cart_coords) @ rot

            else: # Z periodic
                warnings.warn('The periodic direction will be rotated to x axis.')
                rot = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
                lattice_vectors = copy.deepcopy(structure.lattice.matrix) @ rot
                lattice_vectors[1, 2] = latt_mx[0, 0] # Perpendicular elements are not on diagonal after rotation
                lattice_vectors[2, 0] = latt_mx[1, 1]
                atom_coords = copy.deepcopy(structure.cart_coords) @ rot

        else:
            lattice_vectors = copy.deepcopy(structure.lattice.matrix)
            atom_coords = copy.deepcopy(structure.cart_coords)

        structure = CStructure(lattice_vectors, structure.atomic_numbers,
                               atom_coords, coords_are_cartesian=True)

        gui.lattice = structure.lattice.matrix
        gui.n_atoms = structure.num_sites

        if symmetry == True:
            if gui.dimensionality == 3:
                structure.get_sg_symmops(**kwargs)
                gui.space_group = structure.sg
                gui.n_symmops = structure.n_symmops
                gui.symmops = structure.symmops
            elif gui.dimensionality == 2:
                # Get group number before editing- inheriated from previous version
                gui.space_group = SpacegroupAnalyzer(structure, **kwargs).get_space_group_number()
                #center the slab first
                structure = center_slab(structure)
                # Then center at z=0.0
                translation = np.array([0.0, 0.0, -0.5])
                structure.translate_sites(list(range(structure.num_sites)),
                                          translation, to_unit_cell=False)

                sg = SpacegroupAnalyzer(structure)
                ops = sg.get_symmetry_operations(cartesian=True)
                n_symmops = 0; gui.symmops = []
                for op in ops:
                    if np.all(op.translation_vector == 0.):
                        n_symmops += 1
                        gui.symmops.extend(op.rotation_matrix.tolist())
                        gui.symmops.append(op.translation_vector.tolist())

                gui.n_symmops = n_symmops


            else:
                warnings.warn('Check the polymer is correctly centered in the cell and that the correct symmops are used.')
        else:
            gui.space_group = 1
            gui.n_symmops = 1
            gui.symmops = np.vstack([np.eye(3), [0.0,0.0,0.0]])
            gui.symmops = np.reshape(np.array(gui.symmops, dtype=float),
                                     [gui.n_symmops, 4, 3])

        gui.atom_number = list(structure.atomic_numbers)
        gui.atom_positions = structure.cart_coords

    if np.all(zconv!=None):
        for atom in zconv:
            gui.atom_number[atom[0]] = atom[1]

    if np.all(gui_file!=None):
        gui.write_gui(gui_file, symm=symmetry)

    return gui

# Need to check unit and periodicity.
# def cry_cube2xsf(cube, xsf, periodicity='CRYSTAL'):
#     """
#     The shortcut function to convert CRYSTAL CUBE file into XCrySDen XSF format.
#     Useful for 3D scalar field data.

#     .. note::

#         #. CUBE file must follow CRYSTAL convention. The lattice constants are
#         annotated in the second line of the file. Non-periodic directions must
#         along the Z (YZ) axis (axes) for 2D (1D) systems. Data grid is aligned
#         with lattice.

#         #. It's important to specify periodicity of the system with CRYSTAL
#         keywords ('CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE'). In some rare cases
#         the default 3D perodicity might cause file loading error for low
#         dimensional systems. See :ref:`electronics.ChargeDensity.to_xsf() <ref-ChgDensToXSF>`
#         for more information.

#     Args:
#         cube (ChargeDensity | str): `electronics.ChargeDensity` class or name
#             of the CUBE file.
#         xsf (str): Name of the output XSF file.
#         periodicity (str): 'CRYSTAL', 'SLAB', 'POLYMER' or 'MOLECULE'

#     Returns:
#         None
#     """
#     from CRYSTALpytools.electronics import ChargeDensity
#     from CRYSTALpytools.geometry import CStructure
#     from pymatgen.core.lattice import Lattice

#     if not isinstance(cube, ChargeDensity):
#         cube = ChargeDensity.from_file(cube)

#     periodicity = periodicity.upper()
#     if periodicity not in ['CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE']:
#         raise ValueError("Unknown periodicity specification: '{}'.".format(periodicity))
#     pbc = {'CRYSTAL' : (True, True, True),
#            'SLAB'    : (True, True, False),
#            'POLYMER' : (True, False, False),
#            'MOLECULE': (False, False, False)}
#     # Update structure with periodicity
#     cube.structure = CStructure(cube.structure.lattice, cube.structure.species,
#                                 cube.structure.frac_coords, pbc=pbc[periodicity])
#     cube.to_xsf(xsf)
#     return


