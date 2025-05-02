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
    def mapn(cls, filename, index):
        """
        Parse fort.25 file for 2D isovalue maps generated by the 'MAPNAT'
        formatted keyword block ('-%-MAPN'). Unit: a.u.

        Args:
            filename (str): File name
            index (int|list|None): The sequence of '-%-MAPN' headers starting
                from 0. List inputs are accepted. 'None' to read the first 1(2)
                MAPN entries for spin=1(2).

        Returns:
            spin (array): 1, closed shell; 2, open shell
            all_a (array|list): 3D Cartesian coordinates of MAPNET point A.
                1\*3 array of ``index`` is an integer or 1\*1 list. 1\*nIndex
                list otherwise. Same below.
            all_b (array|list): 3D Cartesian coordinates of MAPNET point B
            all_c (array|list): 3D Cartesian coordinates of MAPNET point C
            all_cosxy (float|list): Cosine of vector AB and BC
            all_struc (CStructure|list): Extended Pymatgen Structure object.
            all_map (array|list): 2D scalar field map commensurate with MAPNET
                defined above, nY\*nX\*nSpin.
            unit (str): 'a.u.'
        """
        import pandas as pd
        import numpy as np
        import re, warnings
        from CRYSTALpytools.geometry import CStructure

        df = pd.DataFrame(open(filename))
        # get the head of data
        head_lines = df[df[0].str.contains(r'^-%-[0-4]MAPN')].index.to_numpy(dtype=int)
        if len(head_lines) == 0:
            raise Exception("The 2D grid data header '-%-MAPN' is not found.")
        # get the end of data: element + coordinate lines
        tail_lines = df[df[0].str.contains(r'^\s*[0-9]+\s*[A-Z]+')].index.to_numpy(dtype=int)
        # pair them
        mapn_lines = []
        for head in head_lines:
            tail = tail_lines[np.where(tail_lines>head)[0][0]]
            mapn_lines.append([head, tail])

        # append a very last line for the last entry
        end_lines = df[df[0].str.contains(r'^-%-')].index.to_numpy(dtype=int)
        if end_lines[-1] == head_lines[-1]: # MAPN is the last entry
            mapn_lines.append([df.shape[0], df.shape[0]])
        else: # not the last entry
            endl = df.shape[0]
            for l in end_lines[::-1]:
                if l == head_lines[-1]:
                    break
                endl = l
            mapn_lines.append([endl, endl])

        # spin
        spin = int(df[0][head_lines[0]][3]) % 2 + 1

        # index
        if np.all(index==None):
            if spin == 1: index = 0
            else: index = [0, 1]

        index = np.array(index, ndmin=1, dtype=int)
        if len(index) > len(mapn_lines)-1:
            warnings.warn(
                "{:d} 2D map entries found, while {:d} 2D data grids are required. Only the first {:d} indices are read".format(len(mapn_lines)-1, len(index), len(mapn_lines)-1),
                stacklevel=2
            )
            mapn_lines = mapn_lines[:len(mapn_lines)-1]

        # dimension
        dimen = {0 : (False, False, False),
                 1 : (True, False, False),
                 2 : (True, True, False),
                 3 : (True, True, True)}

        all_map = []; all_a = []; all_b = []; all_c = []; all_cosxy = []; all_struc = []
        for io, o in enumerate(index):
            # base
            points_ab = int(df[0][mapn_lines[o][0]][8:13])
            points_bc = int(df[0][mapn_lines[o][0]][13:18])
            all_cosxy.append(float(df[0][mapn_lines[o][0]][42:54]))
            ab = re.findall(r'\-?[0-9]\.[0-9]{5}E[+,-][0-9][0-9]',
                            df[0][mapn_lines[o][0]+1])
            c = re.findall(r'\-?[0-9]\.[0-9]{5}E[+,-][0-9][0-9]',
                           df[0][mapn_lines[o][0]+2])
            all_a.append(np.array(ab[0:3], dtype=float))
            all_b.append(np.array(ab[3:], dtype=float))
            all_c.append(np.array(c, dtype=float))
            ndimen = int(df[0][mapn_lines[o][0]+2].strip().split()[-1])
            # map
            map = df[0][mapn_lines[o][0]+3:mapn_lines[o][1]].str.findall(
                r'\-?[0-9]\.[0-9]{5}E[+,-][0-9][0-9]').tolist()
            map = np.array([c for row in map for c in row], dtype=float)
            all_map.append(map.reshape([points_ab, points_bc, 1], order='F'))
            # struc
            next_o = o + 1
            # if io + 1 >= len(index):
            #     next_o = index[io] + 1
            # else:
            #     next_o = index[io+1]
            atoms = df[0][mapn_lines[o][1]:mapn_lines[next_o][0]-3].map(
                lambda x: x.strip().split()).tolist()
            species = np.array([i[0] for i in atoms], dtype=int)
            coords = np.array([i[2:] for i in atoms], dtype=float)
            latt = df[0][mapn_lines[next_o][0]-3:mapn_lines[next_o][0]].map(
                    lambda x: x.strip().split()).tolist()
            latt = np.array(latt, dtype=float)
            all_struc.append(CStructure(latt, species, coords, pbc=dimen[ndimen],
                                        coords_are_cartesian=True))

        if len(all_map) == 1:
            all_map = all_map[0]; all_a = all_a[0]; all_b = all_b[0];
            all_c = all_c[0]; all_cosxy = all_cosxy[0]; all_struc = all_struc[0]

        return spin, all_a, all_b, all_c, all_cosxy, all_struc, all_map, 'a.u.'


class DLVParser():
    """
    A collection of functions to parse DLV / Xmgrace files. Instantiation of
    this object is not recommaneded.
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

    @classmethod
    def fort35(cls, filename):
        """
        Keyword 'DLV_BAND'. For 3D band structures in reciprocal space. 3D and
        2D systems only.

        Returns:
            rlatt (array): 3\*3 reciprocal lattice matrix in Bohr.
            bandnew (array): nBand\*nZ\*nY\*nX\*nSpin array of 3D band in Hartree.
            unit (str): 'a.u.'
        """
        import numpy as np
        import pandas as pd

        df = pd.DataFrame(open(filename))
        # remove all the empty lines
        df = df.replace(r'^\s*#*\s*$', np.nan, regex=True).dropna()
        df.reset_index(drop=True, inplace=True)

        iband = df[df[0].str.contains(r'^\s*Band\s+[0-9]+\s*$')].index.tolist()
        nband = 0
        if len(iband) > 0: # no-spin
            iband.append(len(df[0]))
            iband = np.array(iband, ndmin=2, dtype=int)
            nband = iband.shape[1] - 1
            nspin = 1
        else:
            iaband = df[df[0].str.contains(r'^\s*Alpha band\s+[0-9]+\s*$')].index.tolist()
            ibband = df[df[0].str.contains(r'^\s*Beta band\s+[0-9]+\s*$')].index.tolist()
            iaband.append(ibband[0])
            ibband.append(len(df[0]))
            iband = np.vstack([iaband, ibband], dtype=int)
            nband = iband.shape[1] - 1
            nspin = 2
        if nband == 0:
            raise Exception("Band info not found. Check your file name: '{}'.".format(filename))

        npt = np.array(df[0].loc[0].strip().split(), dtype=int)
        rlatt = np.array(
            df[0].loc[2:4].map(lambda x: x.strip().split()).tolist(),
            dtype=float
        )
        rlatt[0] = rlatt[0] * ((npt[0]-1)*0.5+1)
        rlatt[1] = rlatt[1] * ((npt[1]-1)*0.5+1)
        rlatt[2] = rlatt[2] * ((npt[2]-1)*0.5+1)
        band = np.zeros([nband, npt[2], npt[1], npt[0], nspin], dtype=float)
        for i in range(nband):
            for j in range(nspin):
                bg = iband[j, i] + 1
                ed = iband[j, i+1] - 1
                aband = df[0].loc[bg:ed].map(lambda x: x.strip().split()).tolist()
                band[i, :, :, :, j] = np.array(
                    [c for row in aband for c in row], dtype=float
                ).reshape(npt[::-1], order='C')
        del aband
        # Redefine band in recepical cell rather than 1BZ
        # na: points along re-latt vector
        na = int((band.shape[3]-1)/2+1)
        nb = int((band.shape[2]-1)/2+1)
        nc = int((band.shape[1]-1)/2+1)
        if na > 1: newna = na+1
        else: newna = na
        if nb > 1: newnb = nb+1
        else: newnb = nb
        if nc > 1: newnc = nc+1
        else: newnc = nc
        bandnew = np.zeros([nband, newnc, newnb, newna, nspin]) # general grid
        idx = np.where(band<9999)
        bandnew[:, (idx[1]+1)%nc, (idx[2]+1)%nb, (idx[3]+1)%na, :] = band[:, idx[1], idx[2], idx[3], :]
        del band, idx
        if nc > 1: bandnew[:, -1, :, :, :] = bandnew[:, 0, :, :, :]
        if nb > 1: bandnew[:, :, -1, :, :] = bandnew[:, :, 0, :, :]
        if na > 1: bandnew[:, :, :, -1, :] = bandnew[:, :, :, 0, :]
        return rlatt, bandnew, 'a.u.'


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
            map (array): 2D scalar field map commensurate with MAPNET defined
                above. nY\*nX\*1 (spin dimension kept but no use).
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
        map = np.zeros([npt_y, npt_x, 1], dtype=float)

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
            map[i, :regular_entries, :] = tabtmp[tabbg:tabed-1, :].flatten().reshape([-1,1])
            map[i, regular_entries:, :] = tabtmp[tabed-1, 0:last_line_entry].flatten().reshape([-1,1])

        return spin, a, b, c, cosxy, struc, map, 'a.u.'

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


class BOLTZTRAParaser():
    """
    A collection of functions to parse BOLTZTRA output files. Instantiation of
    this object is not recommaneded.
    """
    @classmethod
    def tensor(cls, filename):
        """
        Read properties in tensor forms, include KAPPA, SEEBECK, SIGMA, SIGMAS.

        Returns:
            spin (int): 1 for closed shell and 2 for open shell.
            type (str): Type of tensor output.
            v (float): Volume in cm:math:`^{3}`.
            t (array): Temperature in K.
            mu (array): Chemical potential in eV.
            dc (array): nT\*nMu\*nspin array of carrier density. Unit:
                cm:math:`^{-3}`.
            tensor (array): nT\*nMu\*ndimen\*nspin array of tensor elements.
            unit (str): Unit of tensor.
        """
        import pandas as pd
        import numpy as np

        df = pd.DataFrame(open(filename))
        # remove all the empty lines
        df = df.replace(r'^\s*#*\s*$', np.nan, regex=True).dropna()
        df.reset_index(drop=True, inplace=True)
        spin = int(df[0][0].strip().split()[1]) + 1
        v = float(df[0][2].strip().split()[-1])
        if 'in W/m/K' in df[0][0]:
            type = 'KAPPA'; unit = 'W/m/K'
        elif 'in V/K' in df[0][0]:
            type = 'SEEBECK'; unit = 'V/K'
        elif 'in 1/Ohm/m' in df[0][0]:
            type = 'SIGMA'; unit = '1/Ohm/m'
        elif 'in A/m/K' in df[0][0]:
            type = 'SIGMAS'; unit = 'A/m/K'
        else:
            raise Exception('Unknown property. Check your input BOLTZTRA file.')

        titles = df[df[0].str.contains('# ')].index.tolist()
        mu = df[0][titles[2]+1:titles[3]].map(lambda x: x.strip().split()[0])
        mu = mu.to_numpy(dtype=float)
        ncol = len(df[0][titles[2]+1].strip().split())

        if spin == 1 or type == 'SEEBECK': # no beta state
            t = df[0][titles[2:]].map(lambda x: x.strip().split()[3])
            t = t.to_numpy(dtype=float)
            newtitle = [titles + [len(df[0])]]
        else: # open shell
            betablock = df[df[0].str.contains('Beta')].index.tolist()[0]
            t = df[0][titles[2:int(len(titles)/2)]].map(lambda x: x.strip().split()[3])
            t = t.to_numpy(dtype=float)

            tmp = []; newtitle = []
            for i in titles:
                if i < betablock or i > betablock:
                    tmp.append(i)
                else:
                    newtitle.append(tmp + [betablock])
                    tmp = [betablock]
            newtitle.append(tmp + [len(df[0])])

        dc = np.zeros([spin, len(t), len(mu)], dtype=float)
        tensor = np.zeros([spin, len(t), len(mu), ncol-3], dtype=float)

        for ispin, titles_block in enumerate(newtitle):
            for nt in range(2, len(titles_block)-1):
                block = df[0].loc[titles_block[nt]+1:titles_block[nt+1]-1].map(lambda x: x.strip().split())
                block = np.array(block.tolist(), dtype=float)
                dc[ispin, nt-2, :] = block[:, 2] / v # N carrier to density
                tensor[ispin, nt-2, :, :] = block[:, 3:]

        dc = np.transpose(dc, axes=[1,2,0])
        tensor = np.transpose(tensor, axes=[1,2,3,0])
        if type == 'SEEBECK': # SEEBECK output is not symmetrized
            if tensor.shape[2] == 9: # 3D
                tensor = tensor[:, :, [0,1,2,4,5,8], :]
            elif tensor.shape[2] == 4: # 2D
                tensor = tensor[:, :, [0,1,3], :]
        return spin, type, v, t, mu, dc, tensor, unit

    @classmethod
    def distribution(cls, filename):
        """
        Read transport distribution function.

        Returns:
            spin (int): 1 for closed shell and 2 for open shell.
            type (str): Type. Currently only 'TDF'
            energy (array): nEnergy\*nspin array, Energy in eV.
            distr (array): nEnergy\*nDimen\*nspin array of distribution. Unit:
                :math:`\\hbar^{-2}.eV.fs.\\AA^{-2}`.
        """
        import pandas as pd
        import numpy as np

        df = pd.DataFrame(open(filename))
        # remove all the empty lines
        df = df.replace(r'^\s*#*\s*$', np.nan, regex=True).dropna()
        df.reset_index(drop=True, inplace=True)
        spin = int(df[0][0].strip().split()[1]) + 1
        if 'Transport distribution function' not in df[0][0]:
            raise ValueError('Not a TDF file.')

        # open shell
        data = []
        titles = df[df[0].str.contains('# ')].index
        if spin != 1:
            data.append(df[0][titles[1]+1:titles[2]].map(lambda x: x.strip().split()).tolist())
            data.append(df[0][titles[3]+1:len(df[0])].map(lambda x: x.strip().split()).tolist())
        else:
            data.append(df[0][titles[1]+1:len(df[0])].map(lambda x: x.strip().split()).tolist())

        data = np.array(data, dtype=float)
        data = np.transpose(data, axes=[1,2,0])
        return spin, 'TDF', data[:, 0, :], data[:, 1:, :], '1/hbar^2*eV*fs/angstrom'


class CUBEParser():
    """
    A collection of functions to parse **CRYSTAL** CUBE files. Instantiation of
    this object is not recommaneded.

    .. note::

        Developed especially for CUBE formatted output of CRYSTAL. Lattice
        parameters must be included in the comment line in Bohr and degree.

    """
    @classmethod
    def read_cube(cls, filename):
        """
        Read CUBE formatted files. For base vectors, they are defined by 4 points,
        `origin`, `a`, `b` and `c`.

        Args:
            filename (str):

        Returns:
            origin (array): 1\*3 Cartesian coordinates of origin. Unit: :math:`\\AA`.
            a (array): 1\*3 Cartesian coordinates of grid x base vector end. Unit: :math:`\\AA`.
            b (array): 1\*3 Cartesian coordinates of grid y base vector end. Unit: :math:`\\AA`.
            c (array): 1\*3 Cartesian coordinates of grid z base vector end. Unit: :math:`\\AA`.
            struc (CStructure): Extended Pymatgen Structure object.
            grid (array): nZ\*nY\*nX array of 3D data grid.
            unit (str): Data grid unit, 'a.u.'
        """
        import pandas as pd
        import numpy as np
        from CRYSTALpytools.geometry import CStructure
        from CRYSTALpytools.units import au_to_angstrom
        from pymatgen.core.lattice import Lattice
        from scipy.spatial.transform import Rotation as Rot

        df = pd.DataFrame(open(filename))
        latt = np.array(df[0][1].strip().split(), dtype=float)
        latt[0:3] = au_to_angstrom(latt[0:3])

        # get geom
        ## lattice
        gridv = np.array(df[0].loc[2:5].map(lambda x: x.strip().split()).tolist(),
                         dtype=float)
        natom = int(gridv[0, 0])
        origin = au_to_angstrom(gridv[0, 1:])
        a = au_to_angstrom(gridv[1, 1:] * (gridv[1, 0]-1)) + origin
        b = au_to_angstrom(gridv[2, 1:] * (gridv[2, 0]-1)) + origin
        c = au_to_angstrom(gridv[3, 1:] * (gridv[3, 0]-1)) + origin
        atoms = df[0].loc[6:5+natom].map(lambda x: x.strip().split()).tolist()
        species = np.array([i[0] for i in atoms], dtype=int)
        site_properties = {'charges' : species - np.array([i[1] for i in atoms], dtype=float)}
        coords = np.array([i[2:] for i in atoms], dtype=float)

        ## align lattice vector a with grid base vector a
        lattice = Lattice.from_parameters(a=latt[0], b=latt[1], c=latt[2],
                                          alpha=latt[3], beta=latt[4], gamma=latt[5])
        lattmx = lattice.matrix
        vec1 = lattmx[0] / np.linalg.norm(lattmx[0])
        vec2 = (a-origin) / np.linalg.norm(a-origin)
        rotvec = np.cross(vec1, vec2)
        if np.all(np.abs(rotvec) < 1e-4): rotvec = np.zeros([3,])
        else: rotvec = rotvec / np.linalg.norm(rotvec) * np.arccos(np.dot(vec1, vec2))
        rot = Rot.from_rotvec(rotvec)
        lattice = Lattice(rot.apply(lattmx))
        struc = CStructure(lattice, species, au_to_angstrom(coords), pbc=(True, True, True),
                           site_properties=site_properties, coords_are_cartesian=True)
        # read data
        na = int(gridv[1, 0]); nb = int(gridv[2, 0]); nc = int(gridv[3, 0])
        grid = []
        df[0].loc[6+natom:].map(lambda x: grid.extend(x.strip().split()))
        grid = np.array(grid, dtype=float).reshape([nc, nb, na], order='F')
        return origin, a, b, c, struc, grid, 'a.u.'


