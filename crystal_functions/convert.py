#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:16 2021

"""

from curses.ascii import CR


def cry_out2pmg(output, vacuum=10, initial = False, molecule = True):
    #Transform a CRYSTAL output object into a pymatgen structure object

    # output_file is a crystal output object
    # initial == False reads the last geometry of the output file
    # vacuum needs to be specified because pymatgen does not have 2D symmetry tools
    
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
                return Molecule(output.atom_number, atom_positions)
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


def cry_bands2pmg(output, bands, labels=None):
    # WORK IN TRANSFORMATION
    #Transform a CRYSTAL bands object into a pymatgen bands object

    # output_file is a crystal output object
    # bands is a crystal bands object
        # classBandStructure(kpoints, eigenvals, lattice, efermi, labels_dict=None, coords_are_cartesian=False,
        # structure=None, projections=None)
    # labels are the k point labels to display in the band structure
    
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
    

def cry_gui2pmg(gui, vacuum=10, molecule = True):
    #Transform a CRYSTAL structure (gui) object into a pymatgen Structure object
    #Vacuum needs to be included because pymatgen only includes 3D symmetry
    # molecule = True generates a Molecule pymatgen object for 0D structures
    # molecule = False generates a Molecule pymatgen with vacuum object for 0D structures
    
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
        

def cry_pmg2gui(structure, dimensionality = 3, symmetry = True):
    #Transform a CRYSTAL structure (gui) object into a pymatgen Structure object
    #Vacuum needs to be included because pymatgen only includes 3D symmetry
    # molecule = True generates a Molecule pymatgen object for 0D structures
    # molecule = False generates a Molecule pymatgen with vacuum object for 0D structures

    from file_readwrite import Crystal_gui
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import center_slab
    
    import numpy as np
    import sys

    gui = Crystal_gui()

    if dimensionality == 0 and 'Molecule' not in str(type(structure)):
        print('WARNING: dimensionality is set to 0, but the structure is not a molecule')
        sys.exit(1)
    elif dimensionality == 0 and 'Molecule' in str(type(structure)):
        lattice_vectors = np.identity(3)*500.
    elif dimensionality > 0:
        lattice_vectors = structure.lattice.matrix
    
    gui.dimensionality = dimensionality

    if dimensionality == 2:
        lattice_vectors[2][2] = 500.
    
    if dimensionality == 1:
        lattice_vectors[1][1] = 500.
        lattice_vectors[2][2] = 500.

    gui.lattice = lattice_vectors
    gui.n_atoms = structure.num_sites
    gui.space_group = SpacegroupAnalyzer(structure).get_space_group_number()
    gui.symmops = []
    
    if symmetry == True:
        if dimensionality == 3:
            symmops = SpacegroupAnalyzer(structure).get_symmetry_operations(cartesian=True)
        
        elif dimensionality == 2:

            #center the slab first
            structure = center_slab(structure)

            # Then center at z=0.0
            translation = np.array([0.0, 0.0, -0.5])
            structure.translate_sites(list(range(structure.num_sites)),
                                          translation, to_unit_cell=False)
            sg = SpacegroupAnalyzer(structure)
            ops = sg.get_symmetry_operations(cartesian=True)
            symmops = []
            for op in ops:
                if op.translation_vector[2] == 0.:
                    symmops.extend([op.rotation_matrix.tolist()])
                    symmops.append([op.translation_vector.tolist()])

        elif dimensionality == 1:
            print('WARNING: check the polymer is correctly centered in the cell and that the correct symmops are used.')      

        elif dimensionality == 0:
            print('WARNING: 0D in development')
            sys.exit(1)

        gui.n_symmops = len(symmops)

        for symmop in symmops:
            gui.symmops.extend(symmop.rotation_matrix.tolist())
            gui.symmops.append(symmop.translation_vector.tolist())
    else:
        gui.n_symmops = 1
        gui.symmops.extend(np.identity(3).tolist())
        gui.symmops.append([0.0,0.0,0.0])
    
    gui.atom_number = list(structure.atomic_numbers)
    gui.atom_positions = structure.cart_coords.tolist()
    
    return gui

def cry_ase2gui(structure, dimensionality = 3, symmetry = True):
    # First transform into pmg and then write the gui

    from pymatgen.io.ase import AseAtomsAdaptor

    pmg_structure = AseAtomsAdaptor().get_structure(structure)

    return cry_pmg2gui(pmg_structure, dimensionality=dimensionality, symmetry=symmetry)
    
    
def cry_out2ase(output, initial=False, vacuum=10):
    #Transform a CRYSTAL output object into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    # output_file is a crystal output object
    # initial == False reads the last geometry of the output file
    # dimensionality is the dimensionality of the system  
    # vacuum needs to be specified because pymatgen does not have 2D symmetry tools

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_out2pmg(output,initial=initial,vacuum=vacuum))

def cry_gui2ase(gui_file):
    #Transform a CRYSTAL structure (gui) file into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_gui2pmg(gui_file))


def cry_gui2cif(cif_file_name, gui):
    #Read a CRYSTAL structure (gui) file and save a cif file
    # cif_file_name: name (including path) of the cif file to be saved
    # gui_file: CRYSTAL gui object

    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(gui)
    
    CifWriter(structure).write_file(cif_file_name)

def cry_out2cif(cif_file_name, output):
    #Save a CRYSTAL structure (gui) object as a cif file
    #The gui file is first transfomed into a pymatgen object

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(output)
    
    CifWriter(structure).write_file(cif_file_name)

def cry_gui2xyz(xyz_file_name, gui):
    #Transform a CRYSTAL structure (gui) file into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.cif import XYZ
    import sys

    if gui.dimensionality != 0:
        print('WARNING: the structure is periodic, please use cry_gui2cif()')
        sys.exit(1)
    else:
        structure = cry_gui2pmg(gui, molecule=True) #this returns a pmg Molecule object
    
    XYZ(structure).write_file(xyz_file_name)

def cry_out2xyz(xyz_file_name, output):
    #Transform a CRYSTAL structure (gui) file into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.cif import XYZ
    import sys

    if output.dimensionality != 0:
        print('WARNING: the structure is periodic, please use cry_out2cif()')
        sys.exit(1)
    else:
        structure = cry_gui2pmg(output, molecule=True) #this returns a pmg Molecule object
    
    XYZ(structure).write_file(cif_file_name)


