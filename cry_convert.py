#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:16 2021

@author: brunocamino
"""

def cry_out2pmg(output_file,initial=False):
    #output_file is a crystal output file
    
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from cry_file_readwrite import Crystal_output
    from pymatgen.core import Structure
    
    output = Crystal_output(output_file)
    
    output.extract_last_geom(write_gui_file=False,print_cart=False)
    output.primitive_lattice(initial=initial)
    
    structure = Structure(output.primitive_vectors, output.atom_numbers, output.atom_positions)
    structure_conv = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    
    return structure_conv
    
    
'''###TESTING
cry_out2pmg('examples/data/mgo.out')'''