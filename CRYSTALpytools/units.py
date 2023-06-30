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

def H_to_kjmol(energy):
    # Conversion from Hartree to kJ / mol
    return energy*(constants.physical_constants['Hartree energy'][0] * 1e-3 * constants.Avogadro)

def kjmol_to_H(energy):
    # Conversion from kJ / mol to Hartree
    return energy/(constants.physical_constants['Hartree energy'][0] * 1e-3 * constants.Avogadro)

def au_to_angstrom(length):
    return length*(constants.physical_constants['atomic unit of length'][0] * 1e10)

def angstrom_to_au(length):
    return length/(constants.physical_constants['atomic unit of length'][0] * 1e10)

def cm_to_thz(freq):
    # Conversion from cm^-1 to THz
    return freq*(constants.physical_constants['speed of light in vacuum'][0] * 1e-10)

def thz_to_cm(freq):
    # Conversion from cm^-1 to THz
    return freq/(constants.physical_constants['speed of light in vacuum'][0] * 1e-10)
