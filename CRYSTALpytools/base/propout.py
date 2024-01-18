#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to phrase output files by 'properties' calculations.
"""


class OutBASE():
    """
    Base object for Properties output file. Auxiliary information is
    substracted. Other data is read from formatted files respectively.

    Args:
        filename (str): Properties output file name.
    """

    def __init__(self, filename):
        try:
            file = open(filename, 'r', errors='ignore')
            self.data = file.readlines()
            file.close()
        except:
            raise FileNotFoundError(
                'EXITING: an output file needs to be specified')

    @classmethod
    def get_geometry(cls, filename):
        """
        Get geometry from properties output calculation. A 3D geometry is
        generated since no dimensionality information is provided.

        Args:
            filename (str): Properties output file name.
        Returns:
            struc (Structure): Pymatgen structure
        """
        import re

        import numpy as np
        from pymatgen.core.structure import Structure

        data = cls(filename=filename).data
        countline = 0
        lattice = []
        cart_coord = []
        species = []
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+DIRECT LATTICE VECTOR COMPONENTS', line):
                lattice = [data[countline+1].strip().split(),
                           data[countline+2].strip().split(),
                           data[countline+3].strip().split()]
                countline += 4
                continue
            elif re.match(r'^\s+ATOM N\.AT\.\s+SHELL\s+X\(A\)', line):
                countline += 2
                line = data[countline]
                while not re.match(r'^\s*\*+\s*$', line):
                    line_data = line.strip().split()
                    species.append(line_data[2].capitalize())
                    cart_coord.append(line_data[4:7])
                    countline += 1
                    line = data[countline]

                break
            else:
                countline += 1
                continue

        lattice = np.array(lattice, dtype=float)
        cart_coord = np.array(cart_coord, dtype=float)
        if len(lattice) == 0 or len(cart_coord) == 0:
            raise Exception('Valid geometry not found.')

        return Structure(lattice=lattice, species=species, coords=cart_coord,
                         coords_are_cartesian=True)

    @classmethod
    def get_lattice(cls, filename):
        """
        Get lattice matrix from properties output calculation. A 3D lattice is
        generated since no dimensionality information is provided.

        Args:
            filename (str): Properties output file name.
        Returns:
            matrix (array): 3\*3 lattice matrix
        """
        from CRYSTALpytools.base.propout import OutBASE

        struc = OutBASE.get_geometry(filename)
        return struc.lattice.matrix

    @classmethod
    def get_reciprocal_lattice(cls, filename):
        """
        Get reciprocal lattice matrix from properties output calculation. A 3D
        lattice is generated since no dimensionality information is provided.

        Args:
            filename (str): Properties output file name.
        Returns:
            matrix (array): 3\*3 reciprocal lattice matrix
        """
        from CRYSTALpytools.base.propout import OutBASE

        struc = OutBASE.get_geometry(filename)
        return struc.lattice.reciprocal_lattice.matrix

    @classmethod
    def get_3dkcoord(cls, filename):
        """
        BANDS calculation only. Get 3D coordinates of k points and shrinking
        factors from output file.

        Args:
            filename (str): Properties output file name.
        Returns:
            tick_pos3d (array): ntick\*3 array of fractional coordinates of
                high symmetry k points
            k_pos3d(array): nkpoint\*3 fractional coordinates of k points
        """
        import re

        import numpy as np

        data = cls(filename=filename).data
        is_band = False
        tick_pos3d = []
        k_pos3d = np.array([np.nan, np.nan, np.nan], dtype=float)
        for nline, line in enumerate(data):
            if re.match(r'^\s*\*\s+BAND STRUCTURE\s+\*$', line):
                is_band = True
            elif re.match(r'^\s*LINE\s+[0-9]+\s+\(', line):
                bg = np.array(line[10:25].strip().split(), dtype=float)
                ed = np.array(line[26:41].strip().split(), dtype=float)
                if len(tick_pos3d) > 0:
                    # do not repeat the same point in the middle
                    if np.array_equal(tick_pos3d[-1], bg):
                        tick_pos3d.append(ed)
                    else:
                        tick_pos3d.append(bg)
                        tick_pos3d.append(ed)
                else:
                    tick_pos3d.append(bg)
                    tick_pos3d.append(ed)
            elif re.match(r'^\s*[0-9]+ POINTS \- SHRINKING', line):
                nkp = int(line.strip().split()[0])
                kpos = np.concatenate([np.linspace(bg[0], ed[0], nkp),
                                       np.linspace(bg[1], ed[1], nkp),
                                       np.linspace(bg[2], ed[2], nkp)])
                kpos = np.reshape(kpos, [3, nkp], order='C')
                k_pos3d = np.vstack([k_pos3d, kpos.transpose()])
            elif re.match(r'^\s*[0-9]+ DATA WRITTEN ON UNIT 25', line):
                break

        if is_band == False:
            raise Exception('Not a valid band calculation.')

        tick_pos3d = np.array(tick_pos3d)
        k_pos3d = k_pos3d[1:, :]

        return tick_pos3d, k_pos3d


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
        E Fermi is aligned to 0.

        Args:
            data (list[str]): A list of string (DOSS.DAT).
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, cm_to_thz, eV_to_H, thz_to_cm

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
        doss = np.zeros((n_energy, n_proj+1, spin), dtype=float)
        # Read the doss and store them into a numpy array
        for i, line in enumerate(data[first_energy:first_energy + n_energy]):
            doss[i, :n_proj+1, 0] = np.array([float(n) for n in line.split()])

        if spin == 2:
            # line where the first beta energy is. Written this way to help identify
            first_energy_beta = first_energy + n_energy + 3
            for i, line in enumerate(data[first_energy_beta:-1]):
                doss[i, :n_proj+1, 1] = np.array([float(n)
                                                 for n in line.split()])

        # Convert all the energy to eV / THz
        if is_electron == True:
            doss[:, 0, :] = H_to_eV(doss[:, 0, :])
            # states/Hartree to states/eV
            doss[:, 1:, :] = eV_to_H(doss[:, 1:, :])
        else:
            doss[:, 0, :] = cm_to_thz(doss[:, 0, :])
            doss[:, 1:, :] = thz_to_cm(doss[:, 1:, :])

        return cls(n_energy, n_proj, spin, efermi, doss, unit)

    @classmethod
    def f25_parser(cls, data):
        """
        Parse fort.25 file for electron / phonon DOS. Unit: eV / THz. E Fermi
        is aligned to 0.

        Args:
            data (list[str]): A list of string (fort.25).
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, cm_to_thz, eV_to_H, thz_to_cm

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
                energy_per_block = np.linspace(
                    miny, miny + dy * (npt - 1), npt)
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
            if idx_block < n_proj:  # alpha state
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
            # states/Hartree to states/eV
            doss[:, 1:, :] = eV_to_H(doss[:, 1:, :])
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

    **Uninitialized attributes**
    * ``self.geometry``: Pymatgen structure  
    * ``self.reciprocal_latt``: array, matrix of reciprocal lattice  
    * ``self.tick_pos3d``: array, 3D fractional coordinates of tick labels in reciprocal space  
    * ``self.k_point_pos3d``: array, 3D fractional coordinates of k points in reciprocal space  
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
        # To set the following attributes use 'CRYSTALpytools.base.propout.OutBASE.get_geometry'
        self.geometry = None
        self.reciprocal_latt = None
        # To set the following attributes use 'CRYSTALpytools.base.propout.OutBASE.get_3dkcoord'
        self.k_point_pos3d = None
        self.tick_pos3d = None

    @classmethod
    def BAND_parser(cls, data):
        """
        Parse BAND.DAT / PHONBANDS.DAT file for electron / phonon band structure.
        Unit: eV / THz. E Fermi is aligned to 0.

        Args:
            data (list[str]): A list of string (BAND.DAT).
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, cm_to_thz

        # Read the information about the file
        # number of k points in the calculation
        n_kpoints = int(data[0].split()[2])
        # number of bands in the calculation
        n_bands = int(data[0].split()[4])
        spin = int(data[0].split()[6])  # number of spin
        # number of tick in the band plot
        n_tick = int(data[1].split()[2])+1
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
            bands[:n_bands+1, i,
                  0] = np.array([float(n) for n in line.split()[1:]])
            k_point_plot[i] = float(line.split()[0])

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

        return cls(spin, n_tick, tick_position, tick_label, efermi,
                   n_bands, bands, n_kpoints, k_point_plot, unit)

    @classmethod
    def f25_parser(cls, data):
        """
        Parse fort.25 file for electron / phonon band structure. Unit: eV / THz.
        E Fermi is aligned to 0.

        .. note::

            If Fermi energy is 0, the file is read as phonon band file.

        Args:
            data (list[str]): A list of string (fort.25).
        """
        import re

        import numpy as np

        from CRYSTALpytools.units import H_to_eV, cm_to_thz

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
            if idx_block < nblock:  # alpha state
                k_point_plot = np.concatenate(
                    [k_point_plot, k_in_block[idx_block]])
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
