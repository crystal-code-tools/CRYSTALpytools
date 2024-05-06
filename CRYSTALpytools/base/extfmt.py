#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to parse multiple external output formats by 'crystal' and
'properties' executables, such as 'BAND.DAT' and 'fort.25' formats.
"""


class CrgraParser():
    """
    A collection of functions to parse Crgra fort.25 files. Instantiation of
    this object is not recommaneded.
    """
    @classmethod
    def band(cls, filename):
        """
        Parse fort.25 file for electron / phonon band structure('-%-BAND').
        Unit: eV / THz. E Fermi is aligned to 0.

        Args:
            filename (str): File name

        Returns:
            spin (int): 1, closed shell; 2, open shell
            tick_pos (array): n_tick\*1 array of 1D tick coordinates. Unit: Angstrom
            tick_label (list): n_tick\*1 of default tick labels
            efermi (float): Fermi energy. Unit: eV. 0 for phonon bands.
            bands (array): n_bands\*n_kpoints\*spin array of energy / frequency.
                Unit: eV / THz
            k_path (array): 1D coordinates of k points. Unit: Angstrom
            unit (str): 'eV' or 'THz'
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, au_to_angstrom, cm_to_thz

        file = open(filename, 'r')
        data = file.readlines()
        file.close

        if '-%-' not in data[0] or 'BAND' not in data[0]:
            raise Exception(
                "File '{}' is not a Crgra fort.25 BAND format file.".format(filename))

        data_in_block = []
        k_in_block = []
        n_kpoints = 0
        tick_pos = [0.,]
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
                tick_bg = '({:1d},{:1d},{:1d})'.format(
                    int(tick_line[0]), int(tick_line[1]), int(tick_line[2]))
                tick_ed = '({:1d},{:1d},{:1d})'.format(
                    int(tick_line[3]), int(tick_line[4]), int(tick_line[5]))
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
                if k_in_block == []:  # Initial k path
                    k_per_block = np.linspace(0, dk * (npt - 1), npt)
                else:
                    bg = k_in_block[-1][-1] + dk
                    k_per_block = np.linspace(bg, bg + dk * (npt - 1), npt)
                data_in_block.append(data_per_block)
                k_in_block.append(k_per_block)
                tick_pos.append(k_per_block[-1])
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
        tick_pos = tick_pos[:n_tick]
        tick_label = tick_label[:n_tick]

        k_path = np.array([], dtype=float)
        bands = np.array([], dtype=float)
        for idx_block, block in enumerate(data_in_block):
            if idx_block < nblock:  # alpha state
                k_path = np.concatenate([k_path, k_in_block[idx_block]])
                bands = np.concatenate([bands, block])
            else:
                bands = np.concatenate([bands, block])

        bands = np.reshape(bands, [n_bands, n_kpoints, spin], order='F')

        if is_electron == True:
            bands[:, :, :] = H_to_eV(bands[:, :, :])
        else:
            bands[:, :, :] = cm_to_thz(bands[:, :, :])

        # k coordinates unit. Typically that does not matter
        tick_pos = au_to_angstrom(np.array(tick_pos, dtype=float))
        k_path = au_to_angstrom(k_path)

        return spin, tick_pos, tick_label, efermi, bands, k_path, unit

    @classmethod
    def dos(cls, filename):
        """
        Parse fort.25 file for electron / phonon density of states ('-%-DOSS'
        and '-%-PDOS'). Unit: eV^-1 / THz^-1. E Fermi is aligned to 0. All
        projections must have the same energy / frequency range

        Args:
            filename (str): File name

        Returns:
            spin (array): 1, closed shell; 2, open shell
            efermi (float): Fermi energy. Unit: eV. 0 for phonon bands.
            doss (array): n_proj\*n_energy\*spin array of DOS. Positive values
                for both spin up and spin down states
            energy (int): Number of sampling points (energy or frequency).
            unit (str): 'eV' or 'THz'
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, cm_to_thz, eV_to_H, thz_to_cm

        file = open(filename, 'r')
        data = file.readlines()
        file.close

        if '-%-' not in data[0] or 'DOS' not in data[0]:
            raise Exception(
                "File '{}' is not a Crgra fort.25 PDOS/DOSS format file.".format(filename))
        else:
            # Assuming all projections have the same energy/frequency range
            line = data[0].strip().split()
            npt = int(line[2])
            # Format issue: there might be no space between dx, dy and fermi
            dy = float(data[0][30:42])
            efermi = float(data[0][42:54])
            miny = float(data[1][12:24])
            # Align Fermi energy to 0, consistent with DOSS.DAT file
            energy = np.linspace(miny, miny + dy * (npt - 1), npt) - efermi

        data_in_block = []
        countline = 0
        while countline < len(data):
            line = data[countline]
            if '-%-' in line:
                line = line.strip().split()
                ihferm = int(line[0][3])
                ftype = line[0][4:]
                npt = int(line[2])

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
                data_in_block.append(data_per_block)
            else:
                countline += 1

        nblock = len(data_in_block)
        n_energy = len(energy)
        if ihferm % 2 == 0:
            spin = 1
            n_proj = nblock
        else:
            spin = 2
            n_proj = int(nblock / 2)
        efermi = H_to_eV(efermi)

        doss = np.zeros([n_proj, n_energy, spin], dtype=float)
        for idx_block, block in enumerate(data_in_block):
            if idx_block < n_proj:  # alpha state
                idx_proj = idx_block
                idx_spin = 0
                doss[idx_proj, :, idx_spin] = block
            else:
                idx_proj = idx_block - n_proj
                idx_spin = 1
                doss[idx_proj, :, idx_spin] = -block

        # Convert all the energy to eV
        if ftype == 'DOSS':
            energy = H_to_eV(energy)
            doss = eV_to_H(doss)  # states/Hartree to states/eV
            unit = 'eV'
        elif ftype == 'PDOS':
            energy = cm_to_thz(energy)
            doss = thz_to_cm(doss)
            unit = 'THz'

        return spin, efermi, doss, energy, unit

    @classmethod
    def mapn(cls, filename):
        """
        Parse fort.25 file for 2D isovalue maps of classic electrostatic
        potential (CLAS), charge/spin density (ECHG), electrostatic potential
        (POTM), electron momentum density (EMDPDM) ('-%-MAPN'). Unit: a.u.

        Args:
            filename (str): File name

        Returns:
            spin (array): 1, closed shell; 2, open shell
            a (array): 3D Cartesian coordinates of MAPNET point A
            b (array): 3D Cartesian coordinates of MAPNET point B
            c (array): 3D Cartesian coordinates of MAPNET point C
            cosxy (float): Cosine of vector AB and BC
            struc (CStructure): Extended Pymatgen Structure object.
            map1 (array): 2D scalar field map commensurate with MAPNET defined above.
            map2 (array): *Valid for spin* 2D scalar field map commensurate with MAPNET defined above.
            unit (str): 'a.u.'
        """
        import re

        import numpy as np

        from CRYSTALpytools.geometry import CStructure

        file = open(filename, 'r')
        data = file.readlines()
        file.close

        if '-%-' not in data[0] or 'MAPN' not in data[0]:
            raise Exception(
                "File '{}' is not a Crgra fort.25 2D isovalue map file.".format(filename))

        # Spin
        ihferm = int(data[0][3])
        spin = ihferm % 2 + 1

        points_ab = int(data[0][8:13])
        points_bc = int(data[0][13:18])
        cosxy = float(data[0][42:54])

        a = np.array([data[1][0:12], data[1][12:24],
                     data[1][24:36]], dtype=float)
        b = np.array([data[1][36:48], data[1][48:60],
                     data[1][60:72]], dtype=float)
        c = np.array([data[2][0:12], data[2][12:24],
                     data[2][24:36]], dtype=float)
        pbc_dict = {3: [True, True, True],
                    2: [True, True, False],
                    1: [True, False, False],
                    0: [False, False, False]}
        pbc = pbc_dict[int(data[2][45:])]
        npt = points_ab * points_bc

        countspin = 0
        countpt = 0
        countline = 3
        atom_spec = []
        atom_coord = []
        latt = []
        density_maps = [[], []]
        while countline < len(data):
            line = data[countline]
            if countpt < npt:
                value = re.findall(r'.{12}', line)
                density_maps[countspin].extend(value)
                countpt += len(value)
            else:
                if '-%-' in line:  # Spin block
                    countspin += 1
                    countpt = 0
                    countline += 3
                    continue
                elif re.match(r'^\s+[0-9]+\s+[A-Z]+', line):
                    atom_spec.append(int(line[0:4]))
                    atom_coord.append(
                        [float(line[-60:-40]), float(line[-40:-20]), float(line[-20:])])
                else:
                    latt.append([float(line[0:20]), float(
                        line[20:40]), float(line[40:])])

            countline += 1

        struc = CStructure(latt, atom_spec, atom_coord,
                           coords_are_cartesian=True)
        struc.lattice._pbc = pbc
        map1 = np.array(density_maps[0], dtype=float)
        map1 = np.reshape(map1, [points_ab, points_bc], order='F')
        if countspin != 0:
            map2 = np.array(density_maps[1], dtype=float)
            map2 = np.reshape(map2, [points_ab, points_bc], order='F')
        else:
            map2 = None

        return spin, a, b, c, cosxy, struc, map1, map2, 'a.u.'


class XmgraceParser():
    """
    A collection of functions to parse Xmgrace files (also used for DLV).
    Instantiation of this object is not recommaneded.
    """
    @classmethod
    def band(cls, filename):
        """
        Parse BAND.DAT / PHONBANDS.DAT file for electron / phonon band structure.
        Unit: eV / THz. E Fermi is aligned to 0.

        Args:
            filename (str): BAND.DAT or PHONBANDS.DAT.
        Returns:
            spin (int): 1, closed shell; 2, open shell
            tick_pos (array): n_tick\*1 array of 1D tick coordinates. Unit: Angstrom
            tick_label (list): n_tick\*1 of default tick labels
            efermi (float): Fermi energy. Unit: eV. 0 for phonon bands.
            bands (array): n_bands\*n_kpoints\*spin array of energy / frequency.
                Unit: eV / THz
            k_path (array): 1D coordinates of k points. Unit: Angstrom
            unit (str): 'eV' or 'THz'
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, au_to_angstrom, cm_to_thz

        file = open(filename, 'r')
        data = file.readlines()
        file.close

        if '#' not in data[0] or 'NBND' not in data[0]:
            raise Exception(
                "File '{}' is not a BAND.DAT / PHONBANDS.DAT file.".format(filename))

        # Read the information about the file
        # number of k points in the calculation
        n_kpoints = int(data[0].split()[2])
        # number of bands in the calculation
        n_bands = int(data[0].split()[4])
        spin = int(data[0].split()[6])  # number of spin
        # number of tick in the band plot
        n_tick = int(data[1].split()[2])+1
        tick_pos = []  # positions of the ticks
        tick_label = []  # tick labels
        for i in range(n_tick):
            tick_pos.append(data[16+n_tick+i*2].split()[4])
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
        k_path = np.zeros([n_kpoints,])

        # line where the first band is. Written this way to help identify
        # where the error might be if there are different file lenghts
        first_k = 2 + n_tick + 14 + 2*n_tick + 2

        # Read the bands and store them into a numpy array
        for i, line in enumerate(data[first_k:first_k+n_kpoints]):
            bands[:n_bands+1, i,
                  0] = np.array([float(n) for n in line.split()[1:]])
            k_path[i] = float(line.split()[0])

        if spin == 2:
            # line where the first beta band is. Written this way to help identify
            first_k_beta = first_k + n_kpoints + 15 + 2*n_tick + 2
            for i, line in enumerate(data[first_k_beta:-1]):
                bands[:n_bands+1, i,
                      1] = np.array([float(n) for n in line.split()[1:]])

        if is_electron == True:  # Convert all the energy to eV
            bands[:, :, :] = H_to_eV(bands[:, :, :])
        else:  # Convert all the frequency to THz
            bands[:, :, :] = cm_to_thz(bands[:, :, :])

        # k coordinates unit. Typically that does not matter
        tick_pos = au_to_angstrom(np.array(tick_pos, dtype=float))
        k_path = au_to_angstrom(k_path)

        return spin, tick_pos, tick_label, efermi, bands, k_path, unit

    @classmethod
    def dos(cls, filename):
        """
        Parse DOSS.DAT / PHONDOS.DAT file for electron / phonon density of
        states. Unit: eV^-1 / THz^-1. E Fermi is aligned to 0. All projections
        must have the same energy / frequency range

        Args:
            filename (str): File name

        Returns:
            spin (array): 1, closed shell; 2, open shell
            efermi (float): Fermi energy. Unit: eV. 0 for phonon bands.
            doss (array): n_proj\*n_energy\*spin array of DOS. Positive values
                for both spin up and spin down states
            energy (int): Number of sampling points (energy or frequency).
            unit (str): 'eV' or 'THz'
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, cm_to_thz, eV_to_H, thz_to_cm

        file = open(filename, 'r')
        data = file.readlines()
        file.close
        if '#' not in data[0] or 'NPROJ' not in data[0]:
            raise Exception(
                "File '{}' is not a DOSS.DAT / PHONDOS.DAT file.".format(filename))

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

        if n_proj > 16:  # 16 entries per line at most. A problem for PHONDOS
            raise Exception('Too many projects. Use fort.25 or output file.')

        first_energy = 4
        # Allocate the doss as np arrays
        energy = np.zeros([n_energy,], dtype=float)
        doss = np.zeros([n_energy, n_proj, spin], dtype=float)
        # Read the doss and store them into a numpy array
        for i, line in enumerate(data[first_energy:first_energy + n_energy]):
            line_data = np.array(line.strip().split(), dtype=float)
            energy[i] = line_data[0]
            doss[i, :, 0] = line_data[1:]

        if spin == 2:
            # line where the first beta energy is. Written this way to help identify
            first_energy_beta = first_energy + n_energy + 3
            for i, line in enumerate(data[first_energy_beta:-1]):
                line_data = np.array(line.strip().split(), dtype=float)
                doss[i, :, 1] = -line_data[1:]

        # Convert all the energy to eV / THz
        if is_electron == True:
            energy = H_to_eV(energy)
            doss = eV_to_H(doss)  # states/Hartree to states/eV
            # print(np.shape(doss))
        else:
            energy = cm_to_thz(energy)
            doss = thz_to_cm(doss)

        return spin, efermi, doss, energy, unit
