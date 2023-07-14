#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions that do conversion between data / file formats
"""

def cry_ase2gui(structure, pbc=[True, True, True], symmetry=True):
    """
    Transform an ASE Structure object into a Pymatgen structure object and then
    a CRYSTAL structure (gui) object.

    Args:
        structure (ASE Structure): ASE Structure object.
        symmetry (bool): Perform symmetry analysis.
        
    Returns:
        Crystal_gui: CRYSTAL structure (gui) object.
    """
    # First transform into pmg and then write the gui

    from pymatgen.io.ase import AseAtomsAdaptor

    pmg_structure = AseAtomsAdaptor().get_structure(structure)

    return cry_pmg2gui(pmg_structure, pbc=pbc, symmetry=symmetry)


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
    

def cry_gui2ase(gui_file):
    """
    Transform a CRYSTAL structure (gui) file into an ASE atoms object.
    
    Args:
        gui_file (str): Path to the CRYSTAL structure (gui) file.
        
    Returns:
        Atoms: ASE atoms object.
    """

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_gui2pmg(gui_file))


def cry_gui2cif(cif_file_name, gui, symprec=0.01, angle_tolerance=5.0):
    """
    Read a CRYSTAL structure (gui) file and save a cif file. The CifWriter
    object of PyMatGen is called.

    Args:
        cif_file_name (str): Name (including path) of the cif file to be saved
        gui (Crystal_gui): CRYSTALpytools gui object
        symprec (float): Refer to `CifWriter's manual <https://pymatgen.org/pymatgen.io.cif.html#pymatgen.io.cif.CifWriter>`_.
            If none, CIF output without symmetry.
        angle_tolerance (float): Refer to `CifWriter's manual <https://pymatgen.org/pymatgen.io.cif.html#pymatgen.io.cif.CifWriter>`_
    """
    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(gui)

    CifWriter(structure,
              symprec=symprec,
              angle_tolerance=angle_tolerance).write_file(cif_file_name)


def cry_gui2pmg(gui, vacuum=10, molecule = True):
    """
    Transform a CRYSTAL structure (gui) object into a Pymatgen Structure object.
    
    Args:
        gui: CRYSTAL structure (gui) object.
        vacuum (float): Vacuum distance.
        molecule (bool): Generate a Molecule Pymatgen object for 0D structures.
        
    Returns:
        Structure or Molecule: Pymatgen Structure or Molecule object.
    """
    
    from pymatgen.core.structure import Structure, Molecule
    import numpy as np

    if gui.dimensionality == '0':
        if molecule == True:
            return Molecule(gui.atom_number, gui.atom_positions)
        
        elif molecule == False:
            thickness_x = np.amax(np.array(gui.atom_positions)[:, 0]) - \
                    np.amin(np.array(gui.atom_positions)[:, 0])
            thickness_y = np.amax(np.array(gui.atom_positions)[:, 1]) - \
                    np.amin(np.array(gui.atom_positions)[:, 1])
            thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                    np.amin(np.array(gui.atom_positions)[:, 2])
            
            gui.lattice[0][0] = thickness_x + vacuum
            gui.lattice[1][1] = thickness_y + vacuum
            gui.lattice[2][2] = thickness_z + vacuum

    if gui.dimensionality == '1':
        thickness_y = np.amax(np.array(gui.atom_positions)[:, 1]) - \
                    np.amin(np.array(gui.atom_positions)[:, 1])
        thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                    np.amin(np.array(gui.atom_positions)[:, 2])
        
        gui.lattice[1][1] = thickness_y + vacuum
        gui.lattice[2][2] = thickness_z + vacuum

    if gui.dimensionality == '2':       
        thickness_z = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                    np.amin(np.array(gui.atom_positions)[:, 2])
                
        gui.lattice[2][2] = thickness_z + vacuum

    return Structure(gui.lattice, gui.atom_number, gui.atom_positions, coords_are_cartesian=True)


def cry_gui2xyz(xyz_file_name, gui):
    """
    Transform a CRYSTAL structure (gui) file into an XYZ file.
    
    Args:
        xyz_file_name (str): Name (including path) of the XYZ file to be saved.
        gui: CRYSTAL structure (gui) object.
    """

    from pymatgen.io.xyz import XYZ
    import sys

    if gui.dimensionality != 0:
        print('WARNING: the structure is periodic, please use cry_gui2cif()')
        sys.exit(1)
    else:
        structure = cry_gui2pmg(gui, molecule=True) #this returns a pmg Molecule object
    
    XYZ(structure).write_file(xyz_file_name)


def cry_out2ase(output, initial=False, vacuum=10):
    """
    Transform a CRYSTAL output object into an ASE atoms object.

    Args:
        output: Crystal output object.
        initial (bool): Read the last geometry of the output file.
        vacuum (float): Vacuum distance.
        
    Returns:
        Atoms: ASE atoms object.
    """

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_out2pmg(output,initial=initial,vacuum=vacuum))


def cry_out2cif(cif_file_name, output):
    """
    Save a CRYSTAL output object as a CIF file.
    from pymatgen.io.cif import CifWriter

    Args:
        cif_file_name (str): Name (including path) of the CIF file to be saved.
        output: Crystal output object.
    """


    structure = cry_gui2pmg(output)
    
    CifWriter(structure).write_file(cif_file_name)


def cry_out2pmg(output, vacuum=10, initial = False, molecule = True):
    """
    Transform a CRYSTAL output object into a pymatgen structure object.

    Args:
        output (CRYSTAL output object): CRYSTAL output object.
        vacuum (float): Vacuum distance.
        initial (bool): Read the last geometry of the output file.
        molecule (bool): Generate a Molecule Pymatgen object for 0D structures.
        
    Returns:
        Structure: Pymatgen Structure object.
    """
    
    from pymatgen.core.structure import Structure, Molecule
    import numpy as np

    #Extract information from the output file
    dimensionality = output.get_dimensionality()
    output.get_last_geom(write_gui_file=False)
    atom_positions = output.atom_positions_cart
    vectors = output.get_primitive_lattice(initial=initial)
    

    # Add vacuum for lower dimensionality structures    
    if dimensionality == 0:
        if molecule == True:
            if molecule == True:
                return Molecule(output.atom_numbers, atom_positions)
        elif molecule == False:
            thickness_x = np.amax(np.array(atom_positions)[:, 0]) - \
                        np.amin(np.array(atom_positions)[:, 0])
            thickness_y = np.amax(np.array(atom_positions)[:, 1]) - \
                    np.amin(np.array(atom_positions)[:, 1])
            thickness_z = np.amax(np.array(atom_positions)[:, 2]) - \
                    np.amin(np.array(atom_positions)[:, 2])
        
        vectors[0, 0] = thickness_x + vacuum
        vectors[1, 1] = thickness_y + vacuum
        vectors[2, 2] = thickness_z + vacuum
        
        
    elif dimensionality == 1:
        thickness_y = np.amax(np.array(atom_positions)[:, 1]) - \
                    np.amin(np.array(atom_positions)[:, 1])
        thickness_z = np.amax(np.array(atom_positions)[:, 2]) - \
                    np.amin(np.array(atom_positions)[:, 2])
        
        vectors[1][1] = thickness_y + vacuum
        vectors[2][2] = thickness_z + vacuum

    elif dimensionality == 2:       
        thickness_z = np.amax(np.array(atom_positions)[:, 2]) - \
                    np.amin(np.array(atom_positions)[:, 2])
                
        vectors[2][2] = thickness_z + vacuum

    structure = Structure(vectors, output.atom_numbers, 
                              atom_positions, coords_are_cartesian=True)

    return structure


def cry_pmg2gui(structure, pbc=[True, True, True], symmetry=True, zconv=None):
    """
    Transform a pymatgen Structure object into a CRYSTAL structure (gui) object.
    Vacuum needs to be included because pymatgen only includes 3D symmetry.

    Args:
        structure (Structure | Molecule): Pymatgen Structure / Molecule object.
        pbc (list[bool]): Periodic boundary conditions along the x, y, z directions.
        symmetry (bool): Do symmetry analysis.
        zconv (list[list[int, int]]): 1st element: The **index** of atom;
                2nd element: The new conventional atomic number.
    """
    from CRYSTALpytools.crystal_io import Crystal_gui
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import center_slab
    from pymatgen.core.structure import Molecule

    import numpy as np
    import warnings
    import copy

    gui = Crystal_gui()
    dimensionality = pbc.count(True)

    if dimensionality == 0 and 'Molecule' in str(type(structure)):
        molecule = structure
        is_molecule = True # 0D object called as molecule
    elif dimensionality > 0:
        is_molecule = False # >0D object called as structure
    elif dimensionality == 0 and 'Molecule' not in str(type(structure)):
        warnings.warn('Dimensionality is set to 0, but the structure is not a molecule. Periodicity will be removed.')
        molecule = Molecule(structure.species,
                            [i.coords for i in structure.sites],
                            structure.charge,
                            [i.properties for i in structure.sites]) 
        is_molecule = True # 0D object called as molecule
    elif dimensionality > 0 and 'Molecule' in str(type(structure)):
        warnings.warn('Dimensionality is set to 1-3, but the structure is a molecule. Periodicity will be added.')
        structure = structure.get_boxed_structure(500., 500., 500.)
        is_molecule = False # >0D object called as structure

    gui.dimensionality = dimensionality

    if is_molecule == True: # 0D
        lattice_vectors = np.identity(3)*500.
        gui.lattice = lattice_vectors
        gui.n_atoms = molecule.num_sites
        gui.space_group = 1
        gui.symmops = []
        gui.n_symmops = 1
        gui.symmops.extend(np.identity(3).tolist())
        gui.symmops.append([0.0,0.0,0.0])
        gui.atom_number = list(molecule.atomic_numbers)
        gui.atom_positions = molecule.cart_coords.tolist()
    else: # 1-3D
        if dimensionality == 2:
            if pbc[0] == False: # X no periodicity
                warnings.warn('The non-periodic direction will be rotated to z axis.')
                mx = structure.lattice.matrix
                lattice_vectors = np.array([[mx[1, 1], mx[1, 2], 0.],
                                            [mx[2, 1], mx[2, 2], 0.],
                                            [0., 0., 500.]])
            elif pbc[1] == False: # Y no periodicity
                warnings.warn('The non-periodic direction will be rotated to z axis.')
                mx = structure.lattice.matrix
                lattice_vectors = np.array([[mx[0, 2], mx[0, 0], 0.],
                                            [mx[2, 2], mx[2, 0], 0.],
                                            [0., 0., 500.]])
            else: # Z no periodicity
                lattice_vectors = copy.deepcopy(structure.lattice.matrix)

        elif dimensionality == 1:
            if pbc[0] == True: # X periodic
                lattice_vectors = copy.deepcopy(structure.lattice.matrix)
            elif pbc[1] == True: # Y periodic
                warnings.warn('The periodic direction will be rotated to x axis.')
                mx = structure.lattice.matrix
                lattice_vectors = np.array([[mx[1, 1], 0., 0.],
                                            [0., 500., 0.],
                                            [0., 0., 500.]])
            else: # Z periodic
                warnings.warn('The periodic direction will be rotated to x axis.')
                mx = structure.lattice.matrix
                lattice_vectors = np.array([[mx[2, 2], 0., 0.],
                                            [0., 500., 0.],
                                            [0., 0., 500.]])
        else:
            lattice_vectors = copy.deepcopy(structure.lattice.matrix)

        gui.lattice = lattice_vectors
        gui.n_atoms = structure.num_sites
        gui.space_group = SpacegroupAnalyzer(structure).get_space_group_number()
        gui.symmops = []

        if symmetry == True:
            n_symmops = 0
            if dimensionality == 3:
                symmops = SpacegroupAnalyzer(structure).get_symmetry_operations(cartesian=True)

                for symmop in symmops:
                    if np.all(symmop.translation_vector == 0.):
                        n_symmops += 1
                        gui.symmops.extend(symmop.rotation_matrix.tolist())
                        gui.symmops.append(symmop.translation_vector.tolist())

                gui.n_symmops = n_symmops

            elif dimensionality == 2:

                #center the slab first
                structure = center_slab(structure)

                # Then center at z=0.0
                translation = np.array([0.0, 0.0, -0.5])
                structure.translate_sites(list(range(structure.num_sites)),
                                          translation, to_unit_cell=False)

                sg = SpacegroupAnalyzer(structure)
                ops = sg.get_symmetry_operations(cartesian=True)
                for op in ops:
                    if np.all(op.translation_vector == 0.):
                        n_symmops += 1
                        gui.symmops.extend(op.rotation_matrix.tolist())
                        gui.symmops.append(op.translation_vector.tolist())

                gui.n_symmops = n_symmops

            else:
                warnings.warn('Check the polymer is correctly centered in the cell and that the correct symmops are used.')
        else:
            gui.n_symmops = 1
            gui.symmops.extend(np.identity(3).tolist())
            gui.symmops.append([0.0,0.0,0.0])

        gui.atom_number = list(structure.atomic_numbers)
        gui.atom_positions = structure.cart_coords.tolist()

    if zconv != None:
        for atom in zconv:
            gui.atom_number[atom[0]] = atom[1]

    return gui


def cry_out2xyz(xyz_file_name, output):
    """
    Transform a CRYSTAL output object into an XYZ file.
    
    Args:
        xyz_file_name (str): Name (including path) of the XYZ file to be saved.
        output: CRYSTAL output object.
    """
    from pymatgen.io.xyz import XYZ
    import sys

    if output.dimensionality != 0:
        print('WARNING: the structure is periodic, please use cry_out2cif()')
        sys.exit(1)
    else:
        structure = cry_gui2pmg(output, molecule=True) #this returns a pmg Molecule object
    
    XYZ(structure).write_file(cif_file_name)


