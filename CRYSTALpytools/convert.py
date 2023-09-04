#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions that do conversion between data / file formats
"""

def cry_ase2gui(structure, vacuum=None, symmetry=True):
    """
    Transform an ASE Structure object into a Pymatgen structure object and then
    a CRYSTAL structure (gui) object. Vacuum layer is set to 500 Angstrom
    as the default of CRYSTAL for low symmstry systems

    Args:
        structure (ASE Structure): ASE Structure object.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            length of non-periodic direction to 500 Angstrom.
        symmetry (bool): Perform symmetry analysis.

    Returns:
        Crystal_gui: CRYSTAL structure (gui) object.
    """
    # First transform into pmg and then write the gui

    from pymatgen.io.ase import AseAtomsAdaptor
    from CRYSTALpytools.convert import cry_pmg2gui

    pmg_structure = AseAtomsAdaptor().get_structure(structure)

    return cry_pmg2gui(pmg_structure, vacuum=vacuum, symmetry=symmetry)


def cry_bands2pmg(output, bands, labels=None):
    """
    Transform a CRYSTAL bands object into a Pymatgen bands object.
    
    Args:
        output: Crystal output object.
        bands: Crystal bands object.
        labels (list): K point labels to display in the band structure.
        
    Returns:
        BandStructureSymmLine: Pymatgen band structure object.
    """
    import numpy as np
    from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
    
    from pymatgen.core.lattice import Lattice
    from pymatgen.electronic_structure.core import Spin
    
    # Read the reciprocal lattice from the output file
    output.get_reciprocal_lattice()
    labels_dict = {}
    
    if labels is not None:
        for i, j in enumerate(bands.n_points):
            labels_dict[labels[i]] = bands.k_point_coordinates[j-1]
    
    
    # List of k points coordinates as symmetry lines
    band_energy = bands.bands 
    
    # pymatgen will plot the bands wrt to the Fermi Energy
    band_energy[1:, :, :] = band_energy[1:, :, :] + output.get_fermi_energy()
    k_points_coordinates = []
    
    for i, coord in enumerate(bands.k_point_coordinates):
        k_points_coordinates.append(bands.k_point_coordinates)
        if len(bands.n_points) > 1:
            if i+1 in bands.n_points[1:-1]:
                k_points_coordinates.append(bands.k_point_coordinates)
                
    k_points_coordinates = bands.k_point_coordinates
    if len(bands.n_points) > 1:
        for i, point in enumerate(bands.n_points[1:-1]):
            k_points_coordinates.insert(point-1+i, k_points_coordinates[point-1+i])
            band_energy = np.insert(band_energy, point-1+i, band_energy[:, point-1+i, :], axis=1)
    eigenvals = {Spin.up: band_energy[:, :, 0]}

    if len(bands.bands[0, 0, :]) > 1:
        eigenvals[Spin.down] = band_energy[:, :, 1]

    return BandStructureSymmLine(k_points_coordinates, eigenvals, 
                                 Lattice(output.reciprocal_lattice), 
                                 bands.efermi, labels_dict,
                                 coords_are_cartesian=False)
    

def cry_gui2ase(gui_file, vacuum=None, **kwargs):
    """
    Transform a CRYSTAL structure (gui) file into an ASE atoms object.

    Args:
        gui_file (str): Path to the CRYSTAL structure (gui) file.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of ASE atoms object.
        **kwargs: Passed to ASE Atoms constructor
    Returns:
        Atoms: ASE atoms object.
    """
    from CRYSTALpytools.convert import cry_gui2pmg
    from pymatgen.io.ase import AseAtomsAdaptor

    struc = cry_gui2pmg(gui_file, vacuum=vacuum)

    return AseAtomsAdaptor().get_atoms(struc, **kwargs)


def cry_gui2cif(cif_file_name, gui, vacuum=None, **kwargs):
    """
    Read a CRYSTAL structure (gui) file and save a cif file. The `CifWriter <https://pymatgen.org/pymatgen.io.html#pymatgen.io.cif.CifWriter>`_
    object of Pymatgen is called. By default, ``symprec = 0.01`` is used.

    Args:
        cif_file_name (str): Name (including path) of the cif file to be saved
        gui (Crystal_gui): CRYSTALpytools gui object
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of Pymatgen atoms object.
        **kwargs: Passed to Pymatgen CifWriter.
    """
    from CRYSTALpytools.convert import cry_gui2pmg
    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(gui, vacuum=vacuum, molecule=False)
    if len(kwargs) == 0:
        CifWriter(structure, symprec=0.01, **kwargs).write_file(cif_file_name)
    else:
        CifWriter(structure, **kwargs).write_file(cif_file_name)


def cry_gui2pmg(gui, vacuum=None, molecule=True):
    """
    Transform a CRYSTAL structure (gui) object into a Pymatgen Structure object.

    Args:
        gui: CRYSTAL structure (gui) object.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of Pymatgen object.
        molecule (bool): Generate a Molecule Pymatgen object for 0D structures.

    Returns:
        Structure or Molecule: Pymatgen Structure or Molecule object.
    """

    from pymatgen.core.structure import Structure, Molecule
    from pymatgen.core.lattice import Lattice
    import numpy as np

    if gui.dimensionality == 0:
        if molecule == True:
            return Molecule(gui.atom_number, gui.atom_positions)

        elif molecule == False:
            if vacuum != None:
                pbc = (True, True, True)
                gui.lattice.setflags(write=1)
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
        if vacuum != None:
            pbc = (True, True, True)
            gui.lattice.setflags(write=1)
            thickness_y = np.amax(np.array(gui.atom_positions)[:, 1]) - \
                        np.amin(np.array(gui.atom_positions)[:, 1])
            thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                        np.amin(np.array(gui.atom_positions)[:, 2])

            gui.lattice[1][1] = thickness_y + vacuum
            gui.lattice[2][2] = thickness_z + vacuum
        else:
            pbc = (True, False, False)

    if gui.dimensionality == 2:
        if vacuum != None:
            pbc = (True, True, True)
            gui.lattice.setflags(write=1)
            thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                        np.amin(np.array(gui.atom_positions)[:, 2])

            gui.lattice[2][2] = thickness_z + vacuum
        else:
            pbc = (True, True, False)

    latt = Lattice(gui.lattice, pbc=pbc)

    return Structure(latt, gui.atom_number, gui.atom_positions, coords_are_cartesian=True)


def cry_gui2xyz(xyz_file_name, gui, **kwargs):
    """
    Transform a CRYSTAL structure (gui) file into an XYZ file.

    Args:
        xyz_file_name (str): Name of the XYZ file to be saved.
        gui (Crystal_gui): CRYSTAL structure (gui) object.
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
        output: Crystal output object.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of ASE atoms object.
        initial (bool): Read the last geometry of the output file.
        **kwargs: Passed to ASE Atoms constructor

    Returns:
        Atoms: ASE atoms object.
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from CRYSTALpytools.convert import cry_out2pmg

    struc = cry_out2pmg(output, vacuum=vacuum, initial=initial)

    return AseAtomsAdaptor().get_atoms(struc, **kwargs)


def cry_out2cif(cif_file_name, output, vacuum=None, initial=False, **kwargs):
    """
    Save a CRYSTAL output object as a CIF file. The `CifWriter <https://pymatgen.org/pymatgen.io.html#pymatgen.io.cif.CifWriter>`_
    object of Pymatgen is called. By default, ``symprec = 0.01`` is used.

    Args:
        cif_file_name (str): Name (including path) of the CIF file to be saved.
        output: Crystal output object.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of Pymatgen atoms object.
        initial (bool): Read the last geometry of the output file.
        **kwargs: Passed to Pymatgen CifWriter.
    """
    from CRYSTALpytools.convert import cry_out2pmg
    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(output, vacuum=vacuum, initial=initial, molecule=False)
    if len(kwargs) == 0:
        CifWriter(structure, symprec=0.01, **kwargs).write_file(cif_file_name)
    else:
        CifWriter(structure, **kwargs).write_file(cif_file_name)


def cry_out2pmg(output, vacuum=None, initial=False, molecule=True):
    """
    Transform a CRYSTAL output object into a pymatgen structure object.

    Args:
        output (CRYSTAL output object): CRYSTAL output object.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            ``pbc`` attribute of Pymatgen object.
        initial (bool): Read the last geometry of the output file.
        molecule (bool): Generate a Molecule Pymatgen object for 0D structures.

    Returns:
        Structure: Pymatgen Structure object.
    """
    from CRYSTALpytools.crystal_io import Crystal_output
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure
    import numpy as np

    #Extract information from the output file
    out = Crystal_output().read_cry_output(output)
    ndimen = out.get_dimensionality()
    struc = out.get_geometry(initial=initial, write_gui=False)
    latt_mx = struc.lattice.matrix

    if ndimen == 0:
        if molecule == True:
            return struc

        elif molecule == False:
            if vacuum != None:
                pbc = (True, True, True)
                thickness_x = np.amax(struc.cart_coords[:, 0]) - np.amin(struc.cart_coords[:, 0])
                thickness_y = np.amax(struc.cart_coords[:, 1]) - np.amin(struc.cart_coords[:, 1])
                thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])

                latt_mx[0, 0] = thickness_x + vacuum
                latt_mx[1, 1] = thickness_y + vacuum
                latt_mx[2, 2] = thickness_z + vacuum
            else:
                pbc = (False, False, False)

    if gui.dimensionality == 1:
        if vacuum != None:
            pbc = (True, True, True)
            thickness_y = np.amax(struc.cart_coords[:, 1]) - np.amin(struc.cart_coords[:, 1])
            thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])

            latt_mx[1, 1] = thickness_y + vacuum
            latt_mx[2, 2] = thickness_z + vacuum
        else:
            pbc = (True, False, False)

    if gui.dimensionality == 2:
        if vacuum != None:
            pbc = (True, True, True)
            thickness_z = np.amax(struc.cart_coords[:, 2]) - np.amin(struc.cart_coords[:, 2])

            latt_mx[2, 2] = thickness_z + vacuum
        else:
            pbc = (True, True, False)

    latt = Lattice(latt_mx, pbc=pbc)

    return Structure(latt, struc.atomic_numbers, struc.cart_coords, coords_are_cartesian=True)


def cry_out2xyz(xyz_file_name, output, initial=False, **kwargs):
    """
    Transform a CRYSTAL output object into an XYZ file.

    Args:
        xyz_file_name (str): Name (including path) of the XYZ file to be saved.
        output: CRYSTAL output object.
        initial (bool): Read the last geometry of the output file.
        **kwargs: Passed to Pymatgen XYZ object.
    """
    from pymatgen.io.xyz import XYZ
    from CRYSTALpytools.convert import cry_out2pmg

    structure = cry_out2pmg(output, initial=initial, molecule=True) #this returns a pmg Molecule object
    XYZ(structure, **kwargs).write_file(cif_file_name)


def cry_pmg2gui(structure, vacuum=None, symmetry=True, zconv=None, **kwargs):
    """
    Transform a pymatgen Structure object into a CRYSTAL structure (gui) object.
    Vacuum layer is set to 500 Angstrom as the default of CRYSTAL for low
    symmstry systems

    Args:
        structure (Structure | Molecule): Pymatgen Structure / Molecule object.
        symmetry (bool): Do symmetry analysis.
        vacuum (float): Vacuum distance. Unit: Angstrom. If none, set the
            length of non-periodic direction to 500 Angstrom.
        zconv (list[list[int, int]]): 1st element: The **index** of atom;
                2nd element: The new conventional atomic number.
        **kwargs: Passed to Pymatgen SpacegroupAnalyzer object. Valid only
            if ``symmetry=True``.
    """
    from CRYSTALpytools.crystal_io import Crystal_gui
    from CRYSTALpytools.geometry import get_sg_symmops
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import center_slab
    from pymatgen.core.structure import Structure, Molecule

    import numpy as np
    import warnings

    # dimensionality
    if 'Molecule' in str(type(structure)):
        pbc = (False, False, False)
        structure = Structure(lattice=np.eye(3)*500,
                              species=list(structure.atomic_numbers),
                              coords=structure.cart_coords.tolist(),
                              coords_are_cartesian=True)
    else:
        pbc = structure.pbc

    gui = Crystal_gui()
    gui.dimensionality = pbc.count(True)

    # Vacuum distance
    latt_mx = np.eye(3)*500
    if vacuum != None:
        thickness_x = np.amax(structure.cart_coords[:, 0]) - np.amin(structure.cart_coords[:, 0])
        thickness_y = np.amax(structure.cart_coords[:, 1]) - np.amin(structure.cart_coords[:, 1])
        thickness_z = np.amax(structure.cart_coords[:, 2]) - np.amin(structure.cart_coords[:, 2])
        latt_mx[0, 0] = thickness_x + vacuum
        latt_mx[1, 1] = thickness_y + vacuum
        latt_mx[2, 2] = thickness_z + vacuum

    if gui.dimensionality == 0: # 0D
        gui.lattice = latt_mx
        gui.n_atoms = structure.num_sites
        gui.space_group = 1
        gui.symmops = []
        gui.n_symmops = 1
        gui.symmops = np.vstack([np.eye(3), [0.0,0.0,0.0]])
        gui.atom_number = list(structure.atomic_numbers)
        gui.atom_positions = structure.cart_coords.tolist()
    else: # 1-3D, rotation and add vacuum layer
        if gui.dimensionality == 2:
            if pbc[0] == False: # X no periodicity
                warnings.warn('The non-periodic direction will be rotated to z axis.')
                rot = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
                lattice_vectors = structure.lattice.matrix @ rot
                lattice_vectors[2, 1] = latt_mx[0, 0] # Perpendicular elements are not on diagonal after rotation
                atom_coords = structure.cart_coords @ rot

            elif pbc[1] == False: # Y no periodicity
                warnings.warn('The non-periodic direction will be rotated to z axis.')
                rot = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
                lattice_vectors = structure.lattice.matrix @ rot
                lattice_vectors[2, 0] = latt_mx[1, 1] # Perpendicular elements are not on diagonal after rotation
                atom_coords = structure.cart_coords @ rot

            else: # Z no periodicity
                lattice_vectors = structure.lattice.matrix
                lattice_vectors[2, 2] = latt_mx[2, 2]
                atom_coords = structure.cart_coords

        elif gui.dimensionality == 1:
            if pbc[0] == True: # X periodic
                lattice_vectors = structure.lattice.matrix
                lattice_vectors[1, 1] = latt_mx[1, 1]
                lattice_vectors[2, 2] = latt_mx[2, 2]
                atom_coords = structure.cart_coords

            elif pbc[1] == True: # Y periodic
                warnings.warn('The periodic direction will be rotated to x axis.')
                rot = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
                lattice_vectors = structure.lattice.matrix @ rot
                lattice_vectors[1, 0] = latt_mx[2, 2] # Perpendicular elements are not on diagonal after rotation
                lattice_vectors[2, 1] = latt_mx[0, 0]
                atom_coords = structure.cart_coords @ rot

            else: # Z periodic
                warnings.warn('The periodic direction will be rotated to x axis.')
                rot = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)
                lattice_vectors = structure.lattice.matrix @ rot
                lattice_vectors[1, 2] = latt_mx[0, 0] # Perpendicular elements are not on diagonal after rotation
                lattice_vectors[2, 0] = latt_mx[1, 1]
                atom_coords = structure.cart_coords @ rot

        else:
            lattice_vectors = structure.lattice.matrix
            atom_coords = structure.cart_coords

        structure = Structure(lattice_vectors, structure.atomic_numbers,
                              atom_coords,  coords_are_cartesian=True)

        gui.lattice = structure.lattice.matrix
        gui.n_atoms = structure.num_sites

        if symmetry == True:
            if gui.dimensionality == 3:
                gui.space_group, gui.n_symmops, gui.symmops = get_sg_symmops(structure, **kwargs)
            elif gui.dimensionality == 2:
                # Get group number before editing- inheriated from previous version
                gui.space_group = SpacegroupAnalyzer(structure, **kwargs).get_space_group_number()
                #center the slab first
                structure = center_slab(structure)
                # Then center at z=0.0
                translation = np.array([0.0, 0.0, -0.5])
                structure.translate_sites(list(range(structure.num_sites)),
                                          translation, to_unit_cell=False)
                # Analyze the refined geometry
                _, gui.n_symmops, gui.symmops = get_sg_symmops(structure, **kwargs)
            else:
                warnings.warn('Check the polymer is correctly centered in the cell and that the correct symmops are used.')
        else:
            gui.n_symmops = 1
            gui.symmops = np.vstack([np.eye(3), [0.0,0.0,0.0]])
            gui.symmops = np.reshape(np.array(gui.symmops, dtype=float),
                                     [gui.n_symmops, 4, 3])

        gui.atom_number = list(structure.atomic_numbers)
        gui.atom_positions = structure.cart_coords

    if zconv != None:
        for atom in zconv:
            gui.atom_number[atom[0]] = atom[1]

    return gui

