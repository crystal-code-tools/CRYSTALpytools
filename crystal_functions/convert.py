#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:16 2021

"""

from ctypes import Structure


def cry_out2pmg(output, vacuum=10, initial = False):
    #Transform a CRYSTAL output object into a pymatgen structure object

    # output_file is a crystal output object
    # initial == False reads the last geometry of the output file
    # vacuum needs to be specified because pymatgen does not have 2D symmetry tools
    
    from pymatgen.core import Structure  
    import numpy as np

    #Extract information from the output file
    dimensionality = output.get_dimensionality()
    output.get_last_geom(write_gui_file=False)
    output.get_primitive_lattice(initial=initial)
    
    if dimensionality == 3:
        structure = Structure(output.primitive_lattice, output.atom_numbers, 
                              output.atom_positions_cart, coords_are_cartesian=True)
    # Add vacuum for 2D structures    
    elif dimensionality == 2:
        thickness = np.amax(np.array(output.atom_positions_cart)[:, 2]) - \
                    np.amin(np.array(output.atom_positions_cart)[:, 2])
        
        vectors = output.primitive_vectors
        vectors[2, 2] = thickness + vacuum
        
        structure = Structure(vectors, output.atom_numbers, 
                              output.atom_positions_cart, coords_are_cartesian=True)
    
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
    

def cry_gui2pmg(gui, vacuum=10):
    #Transform a CRYSTAL structure (gui) object into a pymatgen Structure object
    
    from pymatgen.core.structure import Structure
    import numpy as np

    if gui.dimensionality == '2'
        thickness = np.amax(np.array(gui.atom_positions)[:, 2]) - \
                    np.amin(np.array(gui.atom_positions)[:, 2])
            
        gui.lattice[2][2] = thickness + vacuum
        
    return Structure(gui.lattice, gui.atom_number, gui.atom_positions, coords_are_cartesian=True)
        

    
def cry_out2ase(output, initial=False, dimensionality=3, vacuum=10):
    #Transform a CRYSTAL output object into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    # output_file is a crystal output object
    # initial == False reads the last geometry of the output file
    # dimensionality is the dimensionality of the system  
    # vacuum needs to be specified because pymatgen does not have 2D symmetry tools

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_out2pmg(output,initial=initial,dimensionality=dimensionality,vacuum=vacuum))

def cry_gui2ase(gui_file):
    #Transform a CRYSTAL structure (gui) file into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_gui2pmg(gui_file))


def cry_gui2cif(cif_file_name, gui_file):
    #Read a CRYSTAL structure (gui) file and save a cif file
    # cif_file_name: name (including path) of the cif file to be saved
    # gui_file: name (including path) of the gui file to be read

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(gui_file)
    
    CifWriter().write_file(cif_file_name)

def cry_out2cif(cif_file_name, ):
    #Transform a CRYSTAL structure (gui) file into an ASE bands object
    #The gui file is firt transfomed into a pymatgen object

    #gui_file is the CRYSTAL structure (gui) file

    from pymatgen.io.cif import CifWriter

    structure = cry_gui2pmg(gui_file)
    
    CifWriter().write_file(cif_file_name)

