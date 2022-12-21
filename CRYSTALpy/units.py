#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21/12/2022
"""

from scipy import constants
import numpy as np

def H_to_eV(energy):
    # Conversion from Hartree to eV
    return energy/constants.physical_constants['electron volt-hartree relationship'][0]

def eV_to_H(energy):
    # Conversion from eV to Hartree
    return energy*constants.physical_constants['electron volt-hartree relationship'][0]

def au_to_angstrom(length):
    return length*constants.physical_constants['atomic unit of length'][0]/10**-10

def angstrom_to_au(length):
    return length/constants.physical_constants['atomic unit of length'][0]/10**-10
