#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:16 2021

"""

def cry_out2pmg(output, initial=False, vacuum=10):
    # output_file is a crystal output object
    
    from pymatgen.core import Structure  
    import numpy as np

    dimensionality = output.get_dimensionality()
    output.get_last_geom(write_gui_file=False)
    output.get_primitive_lattice(initial=initial)
    
    if dimensionality == 3:
        structure = Structure(output.primitive_lattice, output.atom_numbers, 
                              output.atom_positions_cart, coords_are_cartesian=True)
        
    elif dimensionality == 2:
        thickness = np.amax(np.array(output.atom_positions_cart)[:, 2]) - \
                    np.amin(np.array(output.atom_positions_cart)[:, 2])
        
        vectors = output.primitive_vectors
        vectors[2, 2] = thickness + vacuum
        
        structure = Structure(vectors, output.atom_numbers, 
                              output.atom_positions_cart, coords_are_cartesian=True)
    
    return structure

    
###TESTING
'''from crystal_functions.file_readwrite import Crystal_output
cry_output = Crystal_output('../examples/data/mgo.out')
print(cry_output)
cry_out2pmg(cry_output)'''


def cry_bands2pmg(output, bands, labels=None):
    # WORK IN PROGRESS Function to transform a crystal band file into a pmg band object Format of the pmg object:
    # classBandStructure(kpoints, eigenvals, lattice, efermi, labels_dict=None, coords_are_cartesian=False,
    # structure=None, projections=None)
    
    import numpy as np
    from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
    
    from pymatgen.core.lattice import Lattice
    from pymatgen.electronic_structure.core import Spin

    output.get_reciprocal_lattice()
    labels_dict = {}
    
    if labels is not None:
        for i, j in enumerate(bands.n_points):
            labels_dict[labels[i]] = bands.k_point_coordinates[j-1]
    
    '''#This defines the Kpoint objects. Not needed at the moment, but might be useful in the future
    k_points = []
    for i, coord in enumerate(bands.k_point_coordinates):
        k_points.append(Kpoint(np.array(coord), 
                               Lattice(output.reciprocal_vectors)))
        if len(bands.n_points) > 1:
            if i+1 in bands.n_points[1:-1]:
                print(i+1)
                k_points.append(Kpoint(np.array(coord), 
                               Lattice(output.reciprocal_vectors)))'''
    
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
    


'''cry_output = Crystal_output('../examples/data/mgo_optgeom.out')
cry_bands = Crystal_bands('../examples/data/mgo_BAND_dat.BAND')
bs = cry_bands2pmg(cry_output,cry_bands,labels=['\\Gamma','B','C','\\Gamma','E'])'''

'''###TESTING
from cry_file_readwrite import Crystal_bands, Crystal_output
from pymatgen.core.lattice import Lattice

mgo_bands = Crystal_bands('examples/data/mgo_SPIN_BAND_dat.BAND')
mgo_out = Crystal_output('examples/data/mgo_SPIN.out')
#print(mgo_bands.n_points,labels=['A','B','C','D','E'])
#mgo_out.reciprocal_lattice()


#cry_bands2pmg(bands,output_name)
#mgo_output_name = 'data/mgo_BAND.outp'
cry_bands2pmg(mgo_out,mgo_bands,labels=['A','B','C','D','E'])'''


def cry_gui2pmg(gui_file):

    from pymatgen.core.structure import Structure
    import sys

    try:
        if gui_file[-3:] != 'gui' and gui_file[-3:] != 'f34':
            gui_file = gui_file + '.gui'
        file = open(gui_file, 'r')
        data = file.readlines()
        file.close()
    except:
        print('EXITING: a .gui file needs to be specified')
        sys.exit(1)
    if data[0].split()[0] == '3':
        lattice = []
        for i in range(1, 4):
            lattice.append([float(x) for x in data[i].split()])
        n_symmops = int(data[4].split()[0])
        n_atoms = int(data[5+n_symmops*4].split()[0])
        atom_number = []
        atom_positions = []
        for i in range(6+n_symmops*4, 6+n_symmops*4+n_atoms):
            atom_line = data[i].split()
            atom_number.append(str(atom_line[0]))
            atom_positions.append([float(x) for x in atom_line[1:]])

        return Structure(lattice, atom_number, atom_positions, coords_are_cartesian=True)

    else:
        return 'Lower dimensionality not yet implemented'

    
def cry_out2ase(output, initial=False, dimensionality=3, vacuum=10):

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_out2pmg(output,initial=initial,dimensionality=dimensionality,vacuum=vacuum))

def cry_gui2ase(gui_file):

    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor().get_atoms(cry_gui2pmg(gui_file))
