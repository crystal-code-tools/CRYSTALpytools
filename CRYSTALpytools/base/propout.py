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
        efermi (float): Fermi energy. Unit: eV / THz
        doss (array[float]): n_energy\*(n_proj+1)\*spin array of data.
            The first entry of 2nd dimension is energy / frequency axis.
            Unit: eV / THz. Other entries are data in states / eV or states / THz.
    """
    def __init__(self, n_energy, n_proj, spin, efermi, doss, unit):
        self.n_energy = n_energy
        self.n_proj = n_proj
        self.spin = spin
        self.efermi = efermi
        self.doss = doss
        self.unit = unit

    @classmethod
    def DOSS_parser(cls, data):
        """
        Parse DOSS.DAT / PHONDOS file for electron / phonon DOS. Unit: eV / THz.

        Args:
            data (list[str]): A list of string. Divided by lines.
        """
        import re
        import numpy as np
        from CRYSTALpytools.units import H_to_eV, eV_to_H, cm_to_thz, thz_to_cm

        # Read the information about the file
        n_energy = int(data[0].split()[2])
        n_proj = int(data[0].split()[4])
        spin = int(data[0].split()[6])
        if 'EFERMI' in data[-1]:
            efermi = H_to_eV(float(data[-1].split()[3]))
            is_electron = True
            unit = 'eV'
        else:
            efermi = 0.
            is_electron = False
            unit = 'THz'

        if n_proj > 16: # 16 entries per line at most. A problem for PHONDOS
            raise Exception('Too many projects. Use fort.25 or output file.')

        first_energy = 4
        # Allocate the doss as np arrays
        doss = np.zeros((n_energy, n_proj+1, spin), dtype=float)
        # Read the doss and store them into a numpy array
        for i, line in enumerate(data[first_energy:first_energy + n_energy]):
            doss[i, :n_proj+1, 0] = np.array([float(n) for n in line.split()])

        if spin == 2:
            # line where the first beta energy is. Written this way to help identify
            first_energy_beta = first_energy + n_energy + 3
            for i, line in enumerate(data[first_energy_beta:-1]):
                doss[i, :n_proj+1, 1] = np.array([float(n) for n in line.split()])

        # Convert all the energy to eV / THz
        if is_electron == True:
            doss[:, 0, :] = H_to_eV(doss[:, 0, :])
            doss[:, 1:, :] = eV_to_H(doss[:, 1:, :]) # states/Hartree to states/eV
        else:
            doss[:, 0, :] = cm_to_thz(doss[:, 0, :])
            doss[:, 1:, :] = thz_to_cm(doss[:, 1:, :])

        return cls(n_energy, n_proj, spin, efermi, doss, unit)

    @classmethod
    def f25_parser(cls, data):
        """
        Parse fort.25 file for electron / phonon DOS. Unit: eV / THz.

        Args:
            data (list[str]): A list of string. Divided by lines.
        """
        import numpy as np
        from CRYSTALpytools.units import H_to_eV, eV_to_H, cm_to_thz, thz_to_cm
        import re

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
                    value = re.findall(r'.{12}', line)
                    data_per_block += value
                    countline += 1
                    countpt += len(value)
                data_per_block = np.array(data_per_block, dtype=float)
                # Align Fermi energy to 0, consistent with DOSS file
                energy_per_block = np.linspace(miny, miny + dy * (npt - 1), npt)
                energy_per_block = energy_per_block - efermi
                data_in_block.append(data_per_block)
                energy_in_block.append(energy_per_block)
            else:
                countline += 1

        nblock = len(data_in_block)
        n_energy = npt
        if ihferm % 2 == 0:
            spin = 1
            n_proj = nblock
        else:
            spin = 2
            n_proj = int(nblock / 2)
        efermi = H_to_eV(efermi)

        doss = np.zeros([n_energy, n_proj + 1, spin], dtype=float)
        for idx_block, block in enumerate(data_in_block):
            if idx_block < n_proj: # alpha state
                idx_proj = idx_block + 1
                idx_spin = 0
            else:
                idx_proj = idx_block - n_proj + 1
                idx_spin = 1
            doss[:, 0, idx_spin] = energy_in_block[idx_block]
            doss[:, idx_proj, idx_spin] = block

        # Convert all the energy to eV
        if type == 'DOSS':
            doss[:, 0, :] = H_to_eV(doss[:, 0, :])
            doss[:, 1:, :] = eV_to_H(doss[:, 1:, :]) # states/Hartree to states/eV
            unit = 'eV'
        elif type == 'PDOS':
            doss[:, 0, :] = cm_to_thz(doss[:, 0, :])
            doss[:, 1:, :] = thz_to_cm(doss[:, 1:, :])
            unit = 'THz'

        return cls(n_energy, n_proj, spin, efermi, doss, unit)


class BandsBASE():
    """
    Base object for electron / phonon band structure.

    Args:
        spin (int): 1 or 2. Closed or open shell calculation
        n_tick (int): Number of high symmetric k points
        tick_position (array): n_tick\*1, 1D coordinates of k points along the path
        tick_label (list[str]): Label of k points
        efermi (str): Fermi energy
        n_bands (int): Number of bands
        bands (array): n_bands\*n_kpoints\*spin array of band data. Unit: eV / THz
        n_kpoints (int): Number of k points along the path
        k_point_plot (float): 1D k coordinates along the path
    """
    def __init__(self, spin, n_tick, tick_position, tick_label, efermi,
                 n_bands, bands, n_kpoints, k_point_plot, unit):
        self.spin = spin
        self.n_tick = n_tick
        self.tick_position = tick_position
        self.tick_label = tick_label
        self.efermi = efermi
        self.n_bands = n_bands
        self.bands = bands
        self.n_kpoints = n_kpoints
        self.k_point_plot = k_point_plot
        self.unit = unit
        # Empty and never called by the old Properties_output.read_cry_bands method.
        # But kept anyway
        self.k_point_inp_coordinates = []
        self.n_points = []

    @classmethod
    def BAND_parser(cls, data):
        """
        Parse BAND.DAT / PHONBANDS.DAT file for electron / phonon band structure.
        Unit: eV / THz.

        Args:
            data (list[str]): A list of string. Divided by lines.
        """
        import numpy as np
        from CRYSTALpytools.units import H_to_eV, cm_to_thz
        import re
        # Read the information about the file
        # number of k points in the calculation
        n_kpoints = int(data[0].split()[2])
        # number of bands in the calculation
        n_bands = int(data[0].split()[4])
        spin = int(data[0].split()[6])  # number of spin
        # number of tick in the band plot
        n_tick = int(data[1].split()[2])+1
        """
        # finds all the coordinates of the ticks and the k points
        k_point_inp_coordinates = []
        n_points = []
        for i in range(n_tick):
            n_points.append(int(data[2+i].split()[1]))
            coord = []
            for j in range(3):
                l = re.findall('\d+', data[2+i].split()[2])
                coord.append(float(l[j])/float(l[3]))
            k_point_inp_coordinates.append(coord)
        k_point_inp_coordinates = np.array(k_point_inp_coordinates)
        k_point_coordinates = [k_point_inp_coordinates[0]]
        for i in range(1, n_tick):
            step = (k_point_inp_coordinates[i]- k_point_inp_coordinates[i-1])/float(
                n_points[i]-n_points[i-1])
            for j in range(n_points[i]-n_points[i-1]):
                # coordinates of the k_points in the calculation
                k_point_coordinates.append(
                    (k_point_inp_coordinates[i-1]+step*float(j+1)).tolist())
        """
        tick_position = []  # positions of the ticks
        tick_label = []  # tick labels
        for i in range(n_tick):
            tick_position.append(float(data[16+n_tick+i*2].split()[4]))
            tick_label.append(str(data[17+n_tick+i*2].split()[3][2:]))

        if 'EFERMI' in data[-1]:
            efermi = H_to_eV(float(data[-1].split()[3]))
            is_electron = True
            unit = 'eV'
        else:
            efermi = 0.
            is_electron = False
            unit = 'THz'

        # Allocate the bands as np arrays
        bands = np.zeros([n_bands, n_kpoints, spin], dtype=float)

        # Allocate the k_points a one dimensional array
        k_point_plot = np.zeros([n_kpoints,])

        # line where the first band is. Written this way to help identify
        # where the error might be if there are different file lenghts
        first_k = 2 + n_tick + 14 + 2*n_tick + 2

        # Read the bands and store them into a numpy array
        for i, line in enumerate(data[first_k:first_k+n_kpoints]):
            bands[:n_bands+1, i, 0] = np.array([float(n) for n in line.split()[1:]])
            k_point_plot[i] = float(line.split()[0])

        if spin == 2:
            # line where the first beta band is. Written this way to help identify
            first_k_beta = first_k + n_kpoints + 15 + 2*n_tick + 2
            for i, line in enumerate(data[first_k_beta:-1]):
                bands[:n_bands+1, i, 1] = np.array([float(n) for n in line.split()[1:]])

        if is_electron == True: # Convert all the energy to eV
            bands[:, :, :] = H_to_eV(bands[:, :, :])
        else: # Convert all the frequency to THz
            bands[:, :, :] = cm_to_thz(bands[:, :, :])

        return cls(spin, n_tick, tick_position, tick_label, efermi,
                   n_bands, bands, n_kpoints, k_point_plot, unit)

    @classmethod
    def f25_parser(cls, data):
        """
        Parse fort.25 file for electron / phonon band structure. Unit: eV / THz.

        .. note::

            If Fermi energy is 0, the file is read as phonon band file.

        Args:
            data (list[str]): A list of string. Divided by lines.
        """
        import numpy as np
        from CRYSTALpytools.units import H_to_eV, cm_to_thz
        import re

        data_in_block = []
        k_in_block = []
        n_kpoints = 0
        tick_position = [0.,]
        tick_label = []
        countline = 0
        while countline < len(data):
            line = data[countline]
            if '-%-' in line:
                line = line.strip().split()
                ihferm = int(line[0][3])
                n_bands = int(line[1])
                npt = int(line[2])
                n_kpoints += npt
                # Format issue: there might be no space between dx, dy and fermi
                dk = float(data[countline][30:42])
                efermi = float(data[countline][42:54])

                tick_line = data[countline + 2].strip().split()
                tick_bg = '({:1d},{:1d},{:1d})'.format(int(tick_line[0]), int(tick_line[1]), int(tick_line[2]))
                tick_ed = '({:1d},{:1d},{:1d})'.format(int(tick_line[3]), int(tick_line[4]), int(tick_line[5]))
                if tick_label == []:
                    tick_label = [tick_bg, tick_ed]
                else:
                    tick_label.append(tick_ed)

                countline += 3
                countpt = 0
                data_per_block = []
                while countpt < int(npt * n_bands):
                    line = data[countline]
                    value = re.findall(r'.{12}', line)
                    data_per_block += value
                    countline += 1
                    countpt += len(value)
                # Align Fermi energy to 0, consistent with BAND.DAT file
                data_per_block = np.array(data_per_block, dtype=float) - efermi
                if k_in_block == []: # Initial k path
                    k_per_block = np.linspace(0, dk * (npt - 1), npt)
                else:
                    bg = k_in_block[-1][-1] + dk
                    k_per_block = np.linspace(bg, bg + dk * (npt - 1), npt)
                data_in_block.append(data_per_block)
                k_in_block.append(k_per_block)
                tick_position.append(k_per_block[-1])
            else:
                countline += 1

        if abs(efermi) < 1e-5:
            is_electron = False
            unit = 'THz'
        else:
            is_electron = True
            efermi = H_to_eV(efermi)
            unit = 'eV'

        if ihferm % 2 == 0 or is_electron == False:
            spin = 1
            nblock = len(data_in_block)
        else:
            spin = 2
            nblock = int(len(data_in_block) / 2)
            n_kpoints = int(n_kpoints / 2)

        k_in_block = k_in_block[:nblock]
        n_tick = nblock + 1
        tick_position = tick_position[:n_tick]
        tick_label = tick_label[:n_tick]

        k_point_plot = np.array([], dtype=float)
        bands = np.array([], dtype=float)
        for idx_block, block in enumerate(data_in_block):
            if idx_block < nblock: # alpha state
                k_point_plot = np.concatenate([k_point_plot, k_in_block[idx_block]])
                bands = np.concatenate([bands, block])
            else:
                bands = np.concatenate([bands, block])

        bands = np.reshape(bands, [n_bands, n_kpoints, spin], order='F')

        if is_electron == True:
            bands[:, :, :] = H_to_eV(bands[:, :, :])
        else:
            bands[:, :, :] = cm_to_thz(bands[:, :, :])

        return cls(spin, n_tick, tick_position, tick_label, efermi,
                   n_bands, bands, n_kpoints, k_point_plot, unit)

