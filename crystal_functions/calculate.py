#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:35 2021

@author: brunocamino
"""


def cry_ads_energy(e_full_system, e_substrate, e_adsorbate):
    
    return e_full_system-(e_substrate+e_adsorbate)


def cry_shrink(structure, spacing=0.2):
    # structure is a pymatgen Structure object
    
    import numpy as np

    #short_vector = np.min(np.array(structure.lattice.reciprocal_lattice.lengths))
    vector = np.average(structure.lattice.reciprocal_lattice.lengths)

    return int(np.ceil(vector/spacing))


'''###TESTING
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
bulk = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.597), ["Cu"], [[0, 0, 0]])
cry_shrink(bulk)'''