#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:54 2021

@author: brunocamino
"""
def sub_ads_indices(structure):
    
    substrate_atom_index = []
    adsorbate_atom_index = []
    for i, site in enumerate(structure.sites):
        if "adsorbate" in site.surface_properties:
            print('Adsorbate',i,structure.atomic_numbers[i])
            adsorbate_atom_index.append(i+1)
        elif "surface" in site.surface_properties:
            print('Substrate',i,structure.atomic_numbers[i])
            substrate_atom_index.append(i+1)
    
    indices = {'adsorbate':adsorbate_atom_index,
               'substrate':substrate_atom_index}
    
    return indices
            

###TESTING
from pymatgen.core import Structure, Lattice, Molecule             
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer,PointGroupAnalyzer
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

o = Molecule('O',[[0.0, 0.0, 0.0]])

substrate = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.61491), ["Cu"], [[0, 0, 0]])
substrate = SpacegroupAnalyzer(substrate).get_conventional_standard_structure()

slab = SlabGenerator(substrate, (1,0,0), 5., 10., center_slab=True,max_normal_search=1).get_slab()

system = AdsorbateSiteFinder(slab).adsorb_both_surfaces(o,repeat=[1,1,1])[0]


print(sub_ads_indices(system))
