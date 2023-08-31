#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to phrase 'crystal' output file.
"""
class GeomBASE():
    """
    A container of basic methods for SCF geometry.
    """
    @classmethod
    def read_geom(cls, data, countline):
        """
        Read geometry from 'ATOMS IN THE ASYMMETRIC UNIT' block. It
        terminates at the first empty line.

        .. note::
            ``data`` must not start with 'ATOMS IN THE ASYMMETRIC UNIT' -
            The code refers to 2 lines before it.

        Args:
            data (list[str]): output file read by readlines()
            countline (int): The starting line number

        Returns:
            countline (int): Line number of output file.
            struc (Pymatgen Structure)
        """
        import re
        import numpy as np
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.structure import Structure, Molecule

        pbc = {0 : (False, False, False),
               1 : (True, False, False),
               2 : (True, True, False),
               3 : (True, True, True)}

        species = []
        coords = []
        latt_mx = []
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s*ATOMS IN THE ASYMMETRIC UNIT', line):
                countline += 1
                ndimen = 3 - len(re.findall('\(ANGSTROM\)', data[countline]))

                if ndimen != 0:
                    countline -= 3
                    line_data = data[countline].strip().split()
                    latt = Lattice.from_parameters(
                        a=float(line_data[0]), b=float(line_data[1]),
                        c=float(line_data[2]), alpha=float(line_data[3]),
                        beta=float(line_data[4]), gamma=float(line_data[5]),
                        pbc=pbc[ndimen]
                    )
                    latt_mx = latt.matrix
                else:
                    latt_mx = np.eye(3)

                countline += 5
                line = data[countline].strip()
                coords = []
                species = []
                while line != '':
                    line_data = line.split()
                    coords.append(line_data[-3:])
                    species.append(line_data[3].capitalize())
                    countline += 1
                    line = data[countline].strip()
                countline += 1
                coords = np.array(coords, dtype=float)
                coords[:, 0:ndimen] = coords[:, 0:ndimen] @ latt_mx[0:ndimen, 0:ndimen] # to cartesian coords
                if ndimen != 0:
                    struc = Structure(lattice=latt, species=species, coords=coords,
                                      coords_are_cartesian=True, )
                else:
                    struc = Molecule(species=species, coords=coords)
                break
            else:
                countline += 1

        if len(latt_mx) == 0 or len(species) == 0 or len(coords) == 0:
            raise Exception('Geometry information not found.')

        return countline, struc


class SCFBASE():
    """
    A container of basic methods for SCF loop.
    """
    @classmethod
    def read_convergence(cls, data, countline):
        """
        Read SCF convergence.

        Returns:
            countline (int)
            ncyc (int): Number of cycles
            endflag (str): 'terminated', 'converged', 'too many cycles' and 'unknown'
            e (array): nCYC\*1 array of SCF energy. Unit: eV
            de (array): nCYC\*1 array of SCF energy difference. Unit: eV
        """
        import re
        import warnings
        import numpy as np
        from CRYSTALpytools.units import H_to_eV

        e = []
        de = []
        endflag = 'terminated'
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s*CYC', line):
                line_data = line.strip().split()
                e.append(line_data[3])
                de.append(line_data[5])
                countline += 1
            elif re.match(r'^\s*== SCF ENDED', line):
                if re.search('CONVERGENCE', line):
                    endflag = 'converged'
                elif re.search('TOO MANY CYCLES', line):
                    endflag = 'too many cycles'
                else:
                    warnings.warn('Unknown termination. Check line {}.'.format(countline + 1),
                                  stacklevel=3)
                    endflag = 'unknown'
                break
            else:
                countline += 1

        if endflag != 'converged':
            warnings.warn('SCF convergence not achieved.', stacklevel=3)

        ncyc = len(e)
        e = np.array(e, dtype=float)
        de = np.array(de, dtype=float)

        return countline, ncyc, endflag, H_to_eV(e), H_to_eV(de)

    @classmethod
    def read_fermi_energy(cls, data, countline, history=False):
        """
        Read Fermi energy.

        Args:
            history (bool): Whether to read e fermi of all steps

        Returns:
            countline (int)
            spin (bool): Whether the system is spin-polarised.
            efermi (float | array): Fermi energy. Unit: eV
        """
        import re
        import warnings
        import numpy as np
        from CRYSTALpytools.units import H_to_eV

        if history == True:
            efermi = []
        else:
            efermi = None

        spin = False
        tmp_efermi = []
        while countline >= 0:
            line = data[countline]
            # Metal spin or no spin
            if re.match(r'^ POSSIBLY CONDUCTING STATE - EFERMI', line):
                line_data = line.strip().split()
                if history == False:
                    efermi = float(line_data[5])
                    break
                else:
                    efermi.append(line_data[5])
                    countline -= 1 # Note the reversed sequence here
            # spin flag
            elif re.match(r'^\s*SUMMED SPIN DENSITY', line):
                spin = True
                tmp_efermi = []
                countline -= 1
            # insulator, use top of valence bands
            elif re.match(r'^\s*TOP OF VALENCE BANDS', line):
                line_data = line.strip().split()
                tmp_efermi.append(float(line_data[10]))
                if spin == True and len(tmp_efermi) == 2:
                    if history == False:
                        # Note the reversed sequence here
                        efermi = [tmp_efermi[1], tmp_efermi[0]]
                        break
                    else:
                        efermi.append([tmp_efermi[1], tmp_efermi[0]])
                        tmp_efermi = []
                elif spin == False:
                    if history == False:
                        efermi = tmp_efermi[0]
                        break
                    else:
                        efermi.append(tmp_efermi[0])
                        tmp_efermi = []

                countline -= 1
            # Before the SCF block
            elif re.match(r'^\s*A+$', line):
                break
            else:
                countline -= 1

        if spin == True or history == True:
            efermi = np.array(efermi, dtype=float)
        if history == True:
            # Note the reversed sequence here
            efermi = efermi[::-1]

        return countline, spin, H_to_eV(efermi)

    @classmethod
    def read_band_gap(cls, data, countline, history=False):
        """
        Read band gap.

        Args:
            history (bool): Whether to read band gap of all steps

        Returns:
            countline (int)
            spin (bool): Whether the system is spin-polarised.
            gap (float | array): Band gap. Unit: eV
        """
        import re
        import warnings
        import numpy as np

        if history == True:
            gap = []
        else:
            gap = None

        spin = False
        tmp_gap = []
        while countline >= 0:
            line = data[countline]
            # Metal spin or no spin
            if re.match(r'^ POSSIBLY CONDUCTING STATE', line):
                warnings.warn('Conducting state is identified.', stacklevel=3)
                if history == False:
                    gap = 0.
                    break
                else:
                    gap.append(0.)
                    countline -= 1
            # spin flag
            elif re.match(r'^\s*SUMMED SPIN DENSITY', line):
                spin = True
                tmp_gap = []
                countline -= 1
            # insulator
            elif re.match(r'^\s*.*DIRECT ENERGY BAND GAP', line):
                line_data = line.strip().split()
                tmp_gap.append(float(line_data[4]))
                if spin == True and len(tmp_gap) == 2:
                    if history == False:
                        # Note the reversed sequence here
                        gap = [tmp_gap[1], tmp_gap[0]]
                        break
                    else:
                        gap.append([tmp_gap[1], tmp_gap[0]])
                        tmp_gap = []
                elif spin == False:
                    if history == False:
                        gap = tmp_gap[0]
                        break
                    else:
                        gap.append(tmp_gap[0])
                        tmp_gap = []

                countline -= 1
            # Before the SCF block
            elif re.match(r'^\s*A+$', line):
                break
            else:
                countline -= 1

        if spin == True or history == True:
            gap = np.array(gap, dtype=float)
        if history == True:
            # Note the reversed sequence here
            gap = gap[::-1]

        return countline, spin, gap


class OptBASE():
    """
    A container of basic methods for Opt loop.
    """
    @classmethod
    def read_optblock(cls, data, countline):
        """
        Read optimisation blocks.

        Returns:
            countline (int): Line number of output file.
            ncyc (int): Number of cycles
            endflag (str): 'terminated', 'converged', 'failed' and 'unknown'
            e (array): nCYC\*1 array of total energy. Unit: eV
            de (array): nCYC\*1 array of total energy difference. Unit: eV
            struc (list[Structure]): nCYC\*1 list of pymatgen structures
            maxg (array): nCYC\*1 array of max energy gradient convergence.
                Unit: Hartree / Bohr
            rmsg (array): nCYC\*1 array of RMS energy gradient convergence.
                Unit: Hartree / Bohr
            maxd (array): nCYC\*1 array of max displacement convergence.
                Unit: Bohr
            rmsd (array): nCYC\*1 array of RMS displacement convergence.
                Unit: Bohr
        """
        import re
        import warnings
        import numpy as np
        from CRYSTALpytools.units import H_to_eV
        from CRYSTALpytools.base.crysout import GeomBASE

        e = []
        de = []
        struc = []
        maxg = []
        rmsg = []
        # Displacement: the initial step is 0
        maxd = []
        rmsd = []

        endflag = 'terminated'
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+TOTAL ENERGY\(DFT\)\(AU\)\(.+DE \(AU\)', line):
                line_data = line.strip().split()
                e.append(line_data[3])
                de.append(line_data[6])
                countline += 1
            elif re.match(r'^\s*PRIMITIVE CELL \- CENTRING CODE', line):
                struc.append(GeomBASE.read_geom(data, countline))
            elif re.match(r'^\s+MAX GRADIENT', line):
                line_data = line.strip().split()
                maxg.append(line_data[2])
                countline += 1
            elif re.match(r'^\s+RMS GRADIENT', line):
                line_data = line.strip().split()
                rmsg.append(line_data[2])
                countline += 1
            elif re.match(r'^\s+MAX DISPLAC\.', line):
                line_data = line.strip().split()
                maxd.append(line_data[2])
                countline += 1
            elif re.match(r'^\s+RMS DISPLAC\.', line):
                line_data = line.strip().split()
                rmsd.append(line_data[2])
                countline += 1
            elif re.match(r'\s*\*\s*OPT END', line):
                if re.search('CONVERGED', line):
                    endflag = 'converged'
                elif re.search('FAILED', line):
                    endflag = 'failed'
                else:
                    warnings.warn('Unknown termination. Check line {}.'.format(countline + 1),
                                  stacklevel=3)
                    endflag = 'unknown'
                break
            else:
                countline += 1

        if endflag != 'converged':
            warnings.warn('Optimisation convergence not achieved.', stacklevel=3)

        e = np.array(e, dtype=float)
        de = np.array(de, dtype=float)
        maxg = np.array(maxg, dtype=float)
        rmsg = np.array(rmsg, dtype=float)
        maxd = np.array(maxd, dtype=float)
        rmsd = np.array(rmsd, dtype=float)

        return countline, ncyc, endflag, H_to_eV(e), H_to_eV(de), struc, maxg, rmsg, maxd, rmsd


class PhononBASE():
    """
    A container of basic methods for phonon information.
    """
    @classmethod
    def readmode_basic(cls, data, countline):
        """
        Read basic frequency information.

        Returns:
            countline (int): Line number of output file.
            frequency (array[float]): nmode \* 1
            intens (array[float]): nmode \* 1
            IR (array[bool]): nmode \* 1
            Raman (array[bool]): nmode \* 1
        """
        import re
        import numpy as np

        frequency = []
        intens = []
        IR = []
        Raman = []

        has_spec = False
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+IRREP', line):
                countline += 2
                continue
            elif re.match(r'^\s+\(HARTREE\*\*2\)\s+\(CM\*\*\-1\)', line):
                countline += 1
                continue
            elif line.strip() == '':
                countline += 1
                break
            line_data = line.split()
            nm_a = int(line_data[0].strip('-'))
            nm_b = int(line_data[1])
            freq = float(line_data[4])
            # IR/Raman analysis, closed by default in dispersion calcs
            if 'A' in line or 'I' in line:
                has_spec = True
                intensity = float(line_data[-2].strip(')')) / (nm_b - nm_a + 1)
                ir = line_data[-4] == 'A'
                ram = line_data[-1] == 'A'

            for mode in range(nm_a, nm_b + 1):
                frequency.append(freq)
                if has_spec:
                    intens.append(intensity)
                    IR.append(ir)
                    Raman.append(ram)
            countline += 1

        return countline, frequency, intens, IR, Raman

    @classmethod
    def readmode_eigenvector(cls, data, countline):
        """
        Get mode eigenvectors.

        Returns:
            countline (int): Line number of output file.
            eigvt (array[float]): nmode\*natom\*3 array.
        """
        import numpy as np
        import re

        # Read the eigenvector region as its original shape
        total_data = []
        nmode = 0
        while countline < len(data):
            line = data[countline]
            # Found a block
            if re.match(r'^\s+FREQ\(CM\*\*\-1\)', line):
                countline += 2
                block_data = []
                while data[countline].strip() != '':
                    # Trim annotation part (12 characters)
                    line_data = re.findall(r'\-*[\d\.]+[E\d\-\+]*',
                                           data[countline][13:])
                    if line_data:
                        block_data.append(line_data)
                    countline += 1

                nmode += len(line_data)
                total_data.append(block_data)
                countline += 1
                continue
            # Empty lines
            elif line.strip() == '':
                countline += 1
                continue
            # Other lines
            else:
                break

        eigvt = np.hstack([i for i in total_data]) # (3*natom) * nmode
        eigvt = np.array(eigvt, dtype=float)
        # 1st dimension, nmode
        eigvt = np.transpose(eigvt) # nmode * (3*natom)
        # 2nd dimension, natom
        natom = int(nmode / 3)
        eigvt = np.reshape(eigvt, [nmode, natom, 3], order='C')

        return countline, eigvt

    @classmethod
    def normalize_eigenvector(cls, eigvt, amplitude=1., freq=None, struc=None):
        """
        Normalize the mode of eigenvectors.

        Args:
            eigvt (array[complex]): nmode\*natom\*3 array.
            amplitude (float | str): Amplitude of normalization, Or
                'classical', 'classical-rev', classical amplitude and remove
                classical amplitude.
            freq (array[float]): nmode\*1 array of frequency. 'classical' and 'classical-rev' only.
            struc (Pymatgen Structure): Geometry. 'classical' and 'classical-rev' only.
        Returns:
            eigvt (array[complex]): Normalized eigenvector.
        """
        import numpy as np
        import re
        from scipy import constants
        from CRYSTALpytools.units import amu_to_me, thz_to_hartree

        if re.match('^classical$', amplitude, re.IGNORECASE) or re.match('^classical\-rev$', amplitude, re.IGNORECASE):
            if freq == None or struc == None:
                raise ValueError('Frequency / mass unknown. Cannot generate classical amplitude.')

            nmode = len(freq)
            mass_rev = np.zeros([nmode,]) # In AU
            for i in range(0, nmode, 3):
                atid = i // 3
                atmass = amu_to_me(float(struc.species[atomid].atomic_mass))
                mass_rev[i] = 1 / np.sqrt(atmass)
                mass_rev[i+1] = 1 / np.sqrt(atmass)
                mass_rev[i+2] = 1 / np.sqrt(atmass)

            freq_rev = 1 / np.sqrt(thz_to_hartree(freq))

            if re.match('^classical$', amplitude, re.IGNORECASE):
                for idx_m, m in enumerate(eigvt):
                    eigvt[idx_m] = eigvt[idx_m] * freq_rev[idx_m] * mass_rev[m]
            else:
                for idx_m, m in enumerate(eigvt):
                    eigvt[idx_m] = eigvt[idx_m] / freq_rev[idx_m] / mass_rev[m]

        elif type(amplitude) == float:
            for idx_m, m in enumerate(eigvt):
                eigvt[idx_m] = eigvt[idx_m] / np.linalg.norm(m) * amplitude
        else:
            raise ValueError('Unknown amplitude value: {}'.format(amplitude))

        return eigvt

    @classmethod
    def clean_q_overlap(cls, crysout, threshold):
        """
        Remove the repeated q points at both ends of line segment when
        dispersion is read. The weight of q points will be updated here.

        Args:
            crysout (Crystal_output): :code:`CRYSTALpytools.crystal_io.Crystal_output` object
            threshold (float): The q point overlap threshold.
        """
        import numpy as np
        import warnings

        if crysout.nqpoint <= 1:
            return crysout

        overlap = []
        for idx_q1 in range(crysout.nqpoint - 1):
            qvec1 = crysout.qpoint[idx_q1][0]
            for idx_q2 in range(idx_q1 + 1, crysout.nqpoint): # q1 < q2
                qvec2 = crysout.qpoint[idx_q2][0]
                if np.linalg.norm(qvec1 - qvec2) < threshold:
                    warnings.warn('Overlap of q points is detected between q points {:3d} and {:3d}'.format(idx_q1, idx_q2),
                                  stacklevel=2)
                    overlap.append([idx_q1, idx_q2])

        overlap = np.array(overlap, dtype=int)
        weight = np.array([i[1] for i in crysout.qpoint])
        for idx_q in range(crysout.nqpoint, 1, -1):
            if len(overlap) == 0: # No overlap
                break
            for idx_o in np.where(overlap[:, 1] == idx_q)[0]:
                crysout.nqpoint -= 1
                del crysout.qpoint[overlap[idx_o, 1]]
                weight = np.delete(weight, overlap[idx_o, 1])
                crysout.nmode = np.delete(crysout.nmode, overlap[idx_o, 1])
                crysout.frequency  = np.delete(crysout.frequency, overlap[idx_o, 1], axis=0)
                if len(crysout.intens) != 0:
                    crysout.intens = np.delete(crysout.intens, overlap[idx_o, 1], axis=0)
                    crysout.IR = np.delete(crysout.IR, overlap[idx_o, 1], axis=0)
                    crysout.Raman = np.delete(crysout.Raman, overlap[idx_o, 1], axis=0)
                if len(crysout.eigenvector) != 0:
                    crysout.eigenvector = np.delete(crysout.eigenvector, overlap[idx_o, 1], axis=0)
                break # Avoid repeatly deleting

        # Update weight
        weight = weight * 1 / np.sum(weight)
        for i in range(crysout.nqpoint):
            crysout.qpoint[i][1] = weight[i]

        return crysout

    @classmethod
    def clean_imaginary(cls, crysout, threshold):
        """
        Substitute imaginary modes and corresponding eigenvectors with numpy
        NaN format and print warning message.

        Args:
            crysout (Crystal_output): :code:`CRYSTALpytools.crystal_io.Crystal_output` object
            threshold (float): The threshold to identify a phonon mode as negative.
        """
        import numpy as np
        import warnings

        for q, freq in enumerate(crysout.frequency):
            if np.isnan(freq[0]) or freq[0] > threshold:
                continue

            warnings.warn('Negative frequencies detected.\n  Calculated thermodynamics might be inaccurate. Negative frequencies will be substituted by NaN.',
                          stacklevel=2)

            neg_rank = np.where(freq <= threshold)[0]
            crysout.frequency[q, neg_rank] = np.nan

            if len(crysout.eigenvector) != 0:
                natom = int(crysout.nmode[q] / 3)
                nan_eigvt = np.full([natom, 3], np.nan)
                crysout.eigenvector[q, neg_rank] = nan_eigvt

            if len(crysout.intens) != 0:
                crysout.intens[q, neg_rank] = np.nan
                for n in neg_rank:
                    crysout.IR[q][n] = False
                    crysout.Raman[q][n] = False

        return crysout
