#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:35 2021

@author: brunocamino
"""


def cry_ads_energy(e_full_system, e_substrate, e_adsorbate):
    """
    Calculate the adsorption energy of a system.

    Args:
        e_full_system (float): Total energy of the full system.
        e_substrate (float): Energy of the substrate.
        e_adsorbate (float): Energy of the adsorbate.

    Returns:
        float: Adsorption energy calculated as the difference between the
               total energy of the full system and the sum of the energies
               of the substrate and adsorbate.
    """
    return e_full_system-(e_substrate+e_adsorbate)


def cry_shrink(structure, spacing=0.2):
    """
    Determine the number of unit cells needed to achieve a desired spacing.

    Args:
        structure (pymatgen.core.structure.Structure): The input structure.
        spacing (float): The desired spacing between unit cells. Default is 0.2.

    Returns:
        int: The number of unit cells (rounded up) needed to achieve the desired spacing.
    """
    
    import numpy as np

    vector = np.average(structure.lattice.reciprocal_lattice.lengths)

    return int(np.ceil(vector/spacing))


