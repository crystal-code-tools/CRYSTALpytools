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
        file.close()

        if '-%-' not in data[0]:
            raise Exception("File '{}' is not in Crgra fort.25 format.".format(filename))
        # Band and DOS data might be written into the same f25 file.
        isband = False
        bgline = None
        for nline, line in enumerate(data):
            if 'BAND' in line:
                isband = True
                bgline = nline
                break
            else:
                continue
        if isband != True:
            raise Exception("'BAND' keyword is not found in file '{}'.".format(filename))

        data_in_block = []
        k_in_block = []
        n_kpoints = 0
        tick_pos = [0.,]
        tick_label = []
        countline = bgline
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\-\%\-', line):
                if not re.match(r'^\-\%\-.*BAND', line): # Other data
                    break
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
        file.close()

        if '-%-' not in data[0]:
            raise Exception("File '{}' is not in Crgra fort.25 format.".format(filename))
        # Band and DOS data might be written into the same f25 file.
        isdos = False
        bgline = None
        for nline, line in enumerate(data):
            if 'DOS' in line:
                isdos = True
                bgline = nline
                break
            else:
                continue
        if isdos != True:
            raise Exception("'*DOS*' keyword is not found in file '{}'.".format(filename))

        # Assuming all projections have the same energy/frequency range
        line = data[bgline].strip().split()
        npt = int(line[2])
        # Format issue: there might be no space between dx, dy and fermi
        dy = float(data[bgline][30:42])
        efermi = float(data[bgline][42:54])
        miny = float(data[bgline+1][12:24])
        # Align Fermi energy to 0, consistent with DOSS.DAT file
        energy = np.linspace(miny, miny + dy * (npt - 1), npt) - efermi

        data_in_block = []
        countline = bgline
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\-\%\-', line):
                if not re.match(r'^\-\%\-.*DOS', line): # Other data
                    break
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
        file.close()

        if '-%-' not in data[0]:
            raise Exception("File '{}' is not in Crgra fort.25 format.".format(filename))
        # Multiple data might be written into the same f25 file.
        is2d = False
        bgline = None
        for nline, line in enumerate(data):
            if 'MAPN' in line:
                is2d = True
                bgline = nline
                break
            else:
                continue
        if is2d != True:
            raise Exception("'*MAPN*' keyword is not found in file '{}'.".format(filename))

        # Spin
        ihferm = int(data[bgline][3])
        spin = ihferm % 2 + 1

        points_ab = int(data[bgline][8:13]) # nrow
        points_bc = int(data[bgline][13:18]) # ncol
        cosxy = float(data[bgline][42:54])

        a = np.array([data[bgline+1][0:12], data[bgline+1][12:24],
                     data[bgline+1][24:36]], dtype=float)
        b = np.array([data[bgline+1][36:48], data[bgline+1][48:60],
                     data[bgline+1][60:72]], dtype=float)
        c = np.array([data[bgline+2][0:12], data[bgline+2][12:24],
                     data[bgline+2][24:36]], dtype=float)
        pbc_dict = {3: [True, True, True],
                    2: [True, True, False],
                    1: [True, False, False],
                    0: [False, False, False]}
        pbc = pbc_dict[int(data[bgline+2][45:])]
        npt = points_ab * points_bc

        countspin = 0
        countpt = 0
        countline = bgline+3
        atom_spec = []
        atom_coord = []
        latt = []
        density_maps = [[], []]
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\-\%\-', line) and 'MAPN' not in line: # Other data
                break
            if countpt < npt:
                value = re.findall(r'.{12}', line)
                density_maps[countspin].extend(value)
                countpt += len(value)
            else:
                if re.match(r'^\-\%\-.*MAPN', line):  # Spin block
                    if spin == 1 or countspin == 1: # at most 2 mapnets
                        break
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

        if spin == 2:
            latt = latt[0:3]
            atom_spec = atom_spec[0:int(len(atom_spec)/2)]
            atom_coord = atom_spec[0:int(len(atom_coord)/2)]

        struc = CStructure(latt, atom_spec, atom_coord, coords_are_cartesian=True)
        struc.lattice._pbc = pbc
        map1 = np.array(density_maps[0], dtype=float)
        map1 = np.reshape(map1, [points_ab, points_bc], order='F') # nrow*ncol
        map1 = map1[::-1] # Use BC, BA base vectors, rather than BC, AB.
        if countspin != 0:
            map2 = np.array(density_maps[1], dtype=float)
            map2 = np.reshape(map2, [points_ab, points_bc], order='F')
            map2 = map2[::-1] # Use BC, BA base vectors, rather than BC, AB.
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
        file.close()

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
            energy (int): Number of sampling points (energy or frequency).
            unit (str): 'eV' or 'THz'
        """
        import re
        import numpy as np
        from CRYSTALpytools.units import H_to_eV, cm_to_thz, eV_to_H, thz_to_cm

        file = open(filename, 'r')
        data = file.readlines()
        file.close()
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
        doss = np.zeros([n_proj, n_energy, spin], dtype=float)
        # Read the doss and store them into a numpy array
        for i, line in enumerate(data[first_energy:first_energy + n_energy]):
            line_data = np.array(line.strip().split(), dtype=float)
            energy[i] = line_data[0]
            doss[:, i, 0] = line_data[1:]

        if spin == 2:
            # line where the first beta energy is. Written this way to help identify
            first_energy_beta = first_energy + n_energy + 3
            for i, line in enumerate(data[first_energy_beta:-1]):
                line_data = np.array(line.strip().split(), dtype=float)
                doss[:, i, 1] = -line_data[1:]

        # Convert all the energy to eV / THz
        if is_electron == True:
            energy = H_to_eV(energy)
            doss = eV_to_H(doss)  # states/Hartree to states/eV
        else:
            energy = cm_to_thz(energy)
            doss = thz_to_cm(doss)

        return spin, efermi, doss, energy, unit



class TOPONDParser():
    """
    A collection of functions to parse TOPOND output files. Instantiation of
    this object is not recommaneded.
    """
    @classmethod
    def contour2D(cls, filename):
        """
        Parse TOPOND 2D scalar contour plot files (SURF*.DAT). Unit: a.u.

        Args:
            filename (str)
        Returns:
            spin (array): Always 1
            a (array): 3D Cartesian coordinates of MAPNET point A (xmin, ymax)
            b (array): 3D Cartesian coordinates of MAPNET point B (xmin, ymin)
            c (array): 3D Cartesian coordinates of MAPNET point C (xmax, ymin)
            cosxy (float): Always 0
            struc (None): Always None
            map1 (array): 2D scalar field map commensurate with MAPNET defined above.
            map2 (None): Always None
            unit (str): 'a.u.'
        """
        import numpy as np
        import pandas as pd

        file = open(filename, 'r')
        tmp = file.readline()
        tmp = file.readline()
        npt_x, npt_y = tmp.strip().split()
        tmp = file.readline()
        x_min, x_max, _ = tmp.strip().split()
        tmp = file.readline()
        y_min, y_max, _ = tmp.strip().split()
        file.close()
        npt_x = int(npt_x); npt_y = int(npt_y)

        # To be commensurate with CrgraParser.mapn
        spin = 1
        # Use BC, BA base vectors
        a = np.array([x_min, y_max, 0.], dtype=float)
        b = np.array([x_min, y_min, 0.], dtype=float)
        c = np.array([x_max, y_min, 0.], dtype=float)
        cosxy = 0.
        struc = None
        map1 = np.zeros([npt_y, npt_x], dtype=float)
        map2 = None

        tabtmp = pd.read_table(filename, sep='\s+', skiprows=5, header=None)
        tabtmp = tabtmp.to_numpy(dtype=float)
        nline_per_y = np.ceil(npt_x/np.shape(tabtmp)[1])
        last_line_entry = npt_x % np.shape(tabtmp)[1]
        if last_line_entry == 0:
            last_line_entry = np.shape(tabtmp)[1]

        regular_entries = npt_x-last_line_entry
        for i in range(npt_y):
            tabbg = int(i * nline_per_y)
            tabed = int((i + 1) *nline_per_y)
            map1[i, :regular_entries] = tabtmp[tabbg:tabed-1, :].flatten()
            map1[i, regular_entries:] = tabtmp[tabed-1, 0:last_line_entry]

        return spin, a, b, c, cosxy, struc, map1, map2, 'a.u.'

    @classmethod
    def traj(cls, filename):
        """
        Parse TOPOND trajectory plot files (TRAJ*.DAT). Unit: a.u.

        Args:
            filename (str)
        Returns:
            wtraj (list[int]): 1\*nPath, weight of the path
            traj (list[array]): 1\*nPath, list of critical paths. Every array
                is the nPoint\*3 3D ref framework coordinates of points on the
                path.
            unit (str): 'a.u.'
        """
        import numpy as np
        import re
        import pandas as pd

        wtraj = []; traj = []
        tab = pd.read_fwf(filename, header=None)
        tab = tab.to_numpy(dtype=float)

        countline = 0
        while countline < len(tab):
            # header lines
            line = tab[countline]
            wtraj.append(line[1])
            npt_line = int(line[0])
            traj.append(tab[countline+1:countline+npt_line+1, 1:])
            countline += npt_line+1

        return wtraj, traj, 'a.u.'


