#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:16 2021

@author: brunocamino
"""

def cry_out2pmg(output,initial=False):
    #output_file is a crystal output object
    
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core import Structure
   
    output.extract_last_geom(write_gui_file=False,print_cart=False)
    output.primitive_lattice(initial=initial)
    
    structure = Structure(output.primitive_vectors, output.atom_numbers, output.atom_positions)
    structure_conv = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    
    return structure_conv
    
    
'''###TESTING
cry_out2pmg('examples/data/mgo.out')'''

def cry_bands2pmg(output,bands):
    #WORK IN PROGRESS
    #Function to transform a crystal band file into a pmg band object
    #Format of the pmg object:
    #classBandStructure(kpoints, eigenvals, lattice, efermi, labels_dict=None, coords_are_cartesian=False, structure=None, projections=None)
    import sys
    import re
    import numpy as np
    from pymatgen.electronic_structure.bandstructure import Kpoint, BandStructureSymmLine
    
    from pymatgen.core.lattice import Lattice
    from pymatgen.electronic_structure.core import Spin
    

    output.reciprocal_lattice()
    
    #This defines the Kpoint objects. Not needed at the moment, but might be useful in the future
    k_points = []
    for coord in bands.k_point_coordinates:
        k_points.append(Kpoint(np.array(coord), Lattice(output.reciprocal_vectors)))
    
    eigenvals = {}
    eigenvals[Spin.up] = bands.bands[:,:,0]
    
    if len(bands.bands[0,0,:]) > 1:
        eigenvals[Spin.down] = bands.bands[:,:,1]
    
    return BandStructureSymmLine(bands.k_point_coordinates, eigenvals, Lattice(output.reciprocal_vectors), bands.efermi,labels_dict={})
    
'''###TESTING
from cry_file_readwrite import Crystal_bands, Crystal_output
from pymatgen.core.lattice import Lattice

mgo_bands = Crystal_bands('examples/data/mgo_SPIN_BAND_dat.BAND')
mgo_out = Crystal_output('examples/data/mgo_SPIN.out')
mgo_out.reciprocal_lattice()


#cry_bands2pmg(bands,output_name)
#mgo_output_name = 'data/mgo_BAND.outp'
cry_bands2pmg(mgo_out,mgo_out)'''