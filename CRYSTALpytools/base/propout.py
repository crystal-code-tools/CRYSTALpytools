#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to phrase output files by 'properties' calculations.
"""

class DOSBASE():
    """
    Base object for doensity of states.

    Args:
        n_energy (int): Number of energy points
        n_proj (int): Number of projections
        spin (int): 1 or 2, restricted or open shell calculations
        efermi (float): Fermi energy. Unit: eV
        doss (array[float]): n_energy\*(n_proj+1)\*spin array of data. Unit: eV / cm^-1
    """
    def __init__(self, n_energy, n_proj, spin, efermi, doss):
        self.n_energy = n_energy
        self.n_proj = n_proj
        self.spin = spin
        self.efermi = efermi
        self.doss = doss

    @classmethod
    def DOSS_parser(cls, data):
        """
        Parse DOSS.DAT file for electron DOS.

        Args:
            data (list[str]): A list of string. Divided by lines.
        """
        import re
        import numpy as np
        from CRYSTALpytools import units

        # Read the information about the file
        n_energy = int(data[0].split()[2])
        n_proj = int(data[0].split()[4])
        spin = int(data[0].split()[6])
        efermi = units.H_to_eV(float(data[-1].split()[3]))

        first_energy = 4

        # Allocate the doss as np arrays
        doss = np.zeros((n_energy, n_proj+1, spin), dtype=float)
        # Read the doss and store them into a numpy array
        for i, line in enumerate(data[first_energy:first_energy+n_energy]):
            doss[i, :n_proj+1, 0] = np.array([float(n) for n in line.split()])

        if spin == 2:
            # line where the first beta energy is. Written this way to help identify
            first_energy_beta = first_energy + n_energy + 3
            for i, line in enumerate(data[first_energy_beta:-1]):
                doss[i, :n_proj+1, 1] = np.array([float(n) for n in line.split()])

        # Convert all the energy to eV
        doss[:, 0, :] = units.H_to_eV(doss[:, 0, :])

        return cls(n_energy, n_proj, spin, efermi, doss)

    @classmethod
    def f25_parser(cls, data):
        """
        Parse fort.25 file for electron / phonon DOS

        Args:
            data (list[str]): A list of string. Divided by lines.
        """
        import numpy as np
        from CRYSTALpytools.units import H_to_eV

        data_in_block = []
        energy_in_block = []
        countline = 0
        while countline < len(data):
            line = data[countline]
            if '-%-' in line:
                line = line.strip().split()
                ihferm = int(line[0][3])
                type = line[0][4:]
                npt = int(line[2])
                # Format issue: there might be no space between dx, dy and fermi
                dy = float(data[countline][30:42])
                efermi = float(data[countline][42:54])

                miny = float(data[countline + 1][12:24])

                countline += 3
                countpt = 0
                data_per_block = []
                while countpt < npt:
                    line = data[countline]
                    value = line.strip().split()
                    data_per_block += value
                    countline += 1
                    countpt += len(value)
                data_per_block = np.array(data_per_block, dtype=float)
                energy_per_block = np.linspace(miny, miny + dy * (npt - 1), npt)
                # Align Fermi energy to 0, consistent with DOSS file
                data_per_block = data_per_block - efermi
                energy_per_block = energy_per_block - efermi
                data_in_block.append(data_per_block)
                energy_in_block.append(energy_per_block)
                countline += 1
            else:
                countline += 1

        nblock = len(data_in_block)
        n_energy = npt
        if ihferm % 2 == 0:
            spin = 2
            n_proj = int(nblock / 2)
        else:
            spin = 1
            n_proj = nblock
        efermi = H_to_eV(efermi)

        doss = np.zeros([n_energy, n_proj + 1, spin], dtype=float)
        for idx_block, block in data_in_block:
            if idx_block < n_proj: # alpha state
                idx_spin = 0
            else:
                idx_spin = 1
            doss[:, 0, idx_spin] = energy_in_block[idx_block]
            doss[:, idx_block + 1, idx_spin] = block

        # Convert all the energy to eV
        doss[:, 0, :] = H_to_eV(doss[:, 0, :])

        return cls(n_energy, n_proj, spin, efermi, doss)
