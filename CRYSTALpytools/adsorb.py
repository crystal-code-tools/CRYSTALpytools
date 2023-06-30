#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:54 2021

"""
def sub_ads_indices(structure):
    """
    Identify the indices of substrate and adsorbate atoms in the given structure.

    Args:
        structure (pymatgen.core.structure.Structure): The input structure.

    Returns:
        dict: A dictionary containing the indices of adsorbate and substrate atoms.
              The dictionary has two keys:
              - 'adsorbate': A list of indices corresponding to the adsorbate atoms.
              - 'substrate': A list of indices corresponding to the substrate atoms.
    """
    
    substrate_atom_index = []
    adsorbate_atom_index = []
    for i, site in enumerate(structure.sites):
        if "adsorbate" in site.surface_properties:
            adsorbate_atom_index.append(i+1)
        elif "surface" in site.surface_properties:
            substrate_atom_index.append(i+1)
    
    indices = {'adsorbate':adsorbate_atom_index,
               'substrate':substrate_atom_index}
    
    return indices

