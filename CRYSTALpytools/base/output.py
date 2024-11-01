#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods to parse the output files (screen output, or .out and
.outp files) of 'crystal' and 'properties' executables.
"""
class GeomBASE():
    """
    A container of basic methods for SCF geometry.
    """
    @classmethod
    def read_geom(cls, data):
        """
        Read lattice from 'A B C ALPHA BETA GAMMA' block, periodic boundary
        condition and atom positions from 'ATOMS IN THE ASYMMETRIC UNIT'
        block. It terminates at the first empty line after that block.

        Args:
            data (DataFrame): Pandas DataFrame of the output

        Returns:
            struc (CStructure): Extended Pymatgen Structure
        """
        import pandas as pd
        import numpy as np
        import re
        from pymatgen.core.lattice import Lattice
        from CRYSTALpytools.geometry import CStructure

        pbc = {0 : (False, False, False),
               1 : (True, False, False),
               2 : (True, True, False),
               3 : (True, True, True)}

        species = []
        coords = []
        latt_mx = np.eye(3)

        empty_line = data[data.map(lambda x: x.strip() == '')].index.to_numpy(dtype=int)
        title_line = data[
            data.str.contains(r'^\s+A\s+B\s+C\s+ALPHA\s+BETA\s+GAMMA\s+$')
        ].index.to_numpy(dtype=int)

        # lattice
        if len(title_line) != 0:
            latt_line = np.array(data.loc[title_line[0]+1].strip().split(), dtype=float)
            ndimen = 3 - len(re.findall('\(ANGSTROM\)', data[title_line[0]+4]))
            latt = Lattice.from_parameters(a=latt_line[0], b=latt_line[1],
                                           c=latt_line[2], alpha=latt_line[3],
                                           beta=latt_line[4], gamma=latt_line[5])
            pbc = pbc[ndimen]
        else: # molecules
            latt = Lattice(np.eye(3)*500)
            ndimen = 0
            pbc = (False, False, False)
        latt_mx = latt.matrix

        # Atom coordinates
        atom_line = data[
            data.str.contains(r'^\s*ATOMS IN THE ASYMMETRIC UNIT')
        ].index.to_numpy(dtype=int)

        bg = atom_line[0]
        ed = empty_line[np.where(empty_line>bg)[0][0]]

        coord_list = data.loc[bg+3:ed-1].map(lambda x: x.strip().split()).tolist()
        if len(coord_list) == 0: raise Exception('Geometry information not found.')

        species = [i[2] for i in coord_list] # use conventional atomic numbers
        coords = [i[4:7] for i in coord_list]
        species = np.array(species, dtype=int)
        coords = np.array(coords, dtype=float)

        if ndimen != 0: # to cartesian coords
            coords[:, 0:ndimen] = coords[:, 0:ndimen] @ latt_mx[0:ndimen, 0:ndimen]

        struc = CStructure(lattice=latt, species=species, coords=coords,
                           coords_are_cartesian=True, pbc=pbc)
        return struc


class SCFBASE():
    """
    A container of basic methods for SCF loop.
    """
    @classmethod
    def get_SCF_blocks(cls, data):
        """
        Get SCF convergence block from output file.

        Args:
            data (DataFrame): Pandas DataFrame of the output.

        Returns:
            nSCF (int): Number of SCF blocks
            SCFrange (array[int, int]): The beginning and ending points of every
                SCF block.
        """
        import numpy as np

        scftitle = data[data.str.contains(r'^\s*T+\s+SDIK\s+TELAPSE')].index.to_numpy(dtype=int)
        scfend = data[data.str.contains(r'^\s*== SCF ENDED')].index.to_numpy(dtype=int)
        # This pattern excludes the initial charge assignment but charge info is not always printed out
        realtitle = data[data.str.contains(r'^\s*T+\s+MOQGAD\s+TELAPSE')].index.to_numpy(dtype=int)
        if len(realtitle) == 0:
            raise Exception('SCF block not found. Does it include SCF results?')

        nSCF = len(scftitle)
        SCFrange = np.zeros([nSCF, 2], dtype=int)

        if len(scfend) < len(scftitle):
            scfend.append(len(data)-1)
        for i in range(nSCF):
            SCFrange[i, 0] = realtitle[np.where(realtitle>scftitle[i])[0][0]]
            SCFrange[i, 1] = scfend[i]

        return nSCF, SCFrange

    @classmethod
    def read_convergence(cls, data):
        """
        Read a SCF convergence block.

        Args:
            data (DataFrame): Pandas DataFrame of the SCF convergence block.

        Returns:
            ncyc (int): Number of cycles
            endflag (str): 'terminated', 'converged', 'too many cycles' and 'unknown'
            e (array): nCYC\*1 array of SCF energy. Unit: eV
            de (array): nCYC\*1 array of SCF energy difference. Unit: eV
            spin (bool): Whether the system is spin-polarised.
            efermi (array): Fermi energy. Unit: eV
            gap (array): Band gap. Unit: eV
        """
        import warnings, copy
        import numpy as np
        from CRYSTALpytools.units import H_to_eV

        # ending lines
        scfend = data[data.str.contains(r'^\s*== SCF ENDED')].index
        if len(scfend) < 1:
            endflag = 'terminated'
            warnings.warn('SCF convergence not achieved or missing.', stacklevel=3)
        elif len(scfend) > 1:
            raise ValueError("The 'read_convergence' method only accepts 1 SCF convergence block. Multiple are found.")
        else:
            if 'CONVERGENCE' in data[scfend[0]]:
                endflag = 'converged'
            elif 'TOO MANY CYCLES' in data[scfend[0]]:
                endflag = 'too many cycles'
                warnings.warn('SCF convergence not achieved or missing.', stacklevel=3)
            else:
                warnings.warn('Unknown termination: {}.'.format(data[scfend[0]]),
                              stacklevel=3)
                endflag = 'unknown'

        stepbg = data[data.str.contains(r'^\s*CHARGE NORMALIZATION FACTOR')].index.to_numpy(dtype=int)
        steped = copy.deepcopy(stepbg[1:])
        ncyc = len(stepbg)

        # energies
        energies = data[data.str.contains(r'\s*CYC\s+[0-9]+\s+ETOT\(AU\)')].index
        e = data[energies].map(lambda x: x.strip().split()[3]).to_numpy(dtype=float)
        de = data[energies].map(lambda x: x.strip().split()[5]).to_numpy(dtype=float)
        ## set the first digit to 0 rather than total energy
        de[0] = 0.
        e = H_to_eV(e); de = H_to_eV(de)

        # spin
        spinflag = data[data.str.contains(r'^\s*SUMMED SPIN DENSITY')].index
        if len(spinflag) > 0: spin = True
        else: spin = False

        # Fermi level and gap
        efermi = [0.,]
        if spin == False: gap = [0.,]
        else: gap = [[0., 0.]]
        for bg, ed in zip(stepbg[:-1], steped): # Step 0 has no Fermi and gap
            datamini = data.loc[bg:ed]
            condu = datamini[datamini.str.contains(r'^\s*POSSIBLY CONDUCTING STATE - EFERMI')].index.tolist()
            insul = datamini[datamini.str.contains(r'^\s*TOP OF VALENCE BANDS')].index.tolist()
            gapline = datamini[datamini.str.contains(r'^\s*.*DIRECT ENERGY BAND GAP')].index.tolist()
            # spinlock
            if len(condu) == 0 and len(insul) == 0:
                efermi.append(0.)
                if spin == False: gap.append(0.)
                else: gap.append([0., 0.])
            # conductor
            elif len(condu) > 0 and len(insul) == 0:
                efermi.append(float(datamini[condu[0]].strip().split()[5]))
                if spin == False: gap.append(0.)
                else: gap.append([0., 0.])
            # insulator
            elif len(insul) > 0 and len(condu) == 0:
                allfermi = datamini[insul].map(lambda x: x.strip().split()[10]).to_numpy(dtype=float)
                allgap = datamini[gapline].map(lambda x: x.strip().split()[4]).to_numpy(dtype=float)
                efermi.append(np.max(allfermi))
                if spin == False:
                    gap.append(np.min(allgap))
                else:
                    gstate = int(len(allgap) / 2)
                    gap.append([np.min(allgap[:gstate]), np.min(allgap[gstate:])])

        efermi = H_to_eV(np.array(efermi, dtype=float))
        gap = np.array(gap, dtype=float)
        return ncyc, endflag, e, de, spin, efermi, gap


class OptBASE():
    """
    A container of basic methods for Opt loop.
    """
    @classmethod
    def get_opt_block(cls, data):
        """
        Get optimization convergence block (every OPT step) from output file.

        Args:
            data (DataFrame): Pandas DataFrame of the output.

        Returns:
            nOPT (int): Number of OPT steps
            OPTrange (array[int, int]): The beginning and ending points of every
                OPT step.
            endflag (str): 'terminated', 'converged', 'failed' and 'unknown'
        """
        import numpy as np
        import warnings

        opttitle = data[data.str.contains(r'^\s*[A-Z]+ OPTIMIZATION - POINT')].index.to_numpy(dtype=int)
        optend = data[data.str.contains(r'^\s*T+ OPTI\s+TELAPSE')].index.to_numpy(dtype=int)
        block_end = data[data.str.contains(r'^\s*\* OPT END')].index.to_numpy(dtype=int)
        # Include initial SCF step
        # opttitle: step 1 to final run
        # optend: initial SCF to last before final run
        if len(opttitle) == 0: raise Exception('Not an optimization output.')
        if len(opttitle) == 1 and len(block_end) == 0: raise Exception('Initial SCF failed. Nothing to substract.')
        # terminated
        if len(block_end) == 0:
            warnings.warn('Job interrupted. Not a complete file.', stacklevel=3)
            block_end = np.array([data.index[-1]], dtype=int)
            endflag = 'terminated'
        # normal
        else:
            if 'CONVERGED' in data.loc[block_end[0]]:
                endflag = 'converged'
            elif 'FAILED' in data.loc[block_end[0]]:
                endflag = 'failed'
                warnings.warn('Convergence not achieved.', stacklevel=3)
            else:
                warnings.warn('Unknown termination: {}.'.format(data.loc[block_end[0]]),
                                  stacklevel=3)
                endflag = 'unknown'
        ## get ranges
        nOPT = len(opttitle)
        OPTrange = np.zeros([nOPT, 2], dtype=int)
        OPTrange[:, 0] = opttitle
        OPTrange[-1, 1] = block_end[0]
        if nOPT > 1:
            OPTrange[:-1, 1] = optend
        return nOPT, OPTrange, endflag

    @classmethod
    def read_opt_block(cls, data):
        """
        Read information of every OPT step from output file.

        Args:
            data (DataFrame): Pandas DataFrame of the output.

        Returns:
            e (float): Final SCF energy with corrections. Unit: eV
            de (float): Final SCF energy difference with last OPT step. Unit: eV
            struc (CStructure): Modified pymatgen structure.
            maxg (float): Max energy gradient convergence. Unit: Hartree / Bohr.
            rmsg (float): RMS energy gradient convergence. Unit: Hartree / Bohr,
            maxd (float): Max displacement convergence. Unit: Bohr.
            rmsd (float): RMS displacement convergence. Unit: Bohr.
        """
        import numpy as np
        from CRYSTALpytools.units import H_to_eV

        eline = data[data.str.contains(r'^\s+TOTAL ENERGY\(DFT\)\(AU\)\(')].index.to_numpy(dtype=int)
        line = data.loc[eline[-1]].strip().split()
        e = H_to_eV(float(line[3]))
        gxline = data[data.str.contains(r'^\s+MAX GRADIENT')].index.to_numpy(dtype=int)
        gmline = data[data.str.contains(r'^\s+RMS GRADIENT')].index.to_numpy(dtype=int)
        maxg = float(data.loc[gxline[-1]].strip().split()[2])
        rmsg = float(data.loc[gmline[-1]].strip().split()[2])

        if 'POINT    1' in data.iloc[0]: # initial step, no structure / displacement
            de = 0.
            struc = None
            maxd = 0.
            rmsd = 0.
        else:
            de = H_to_eV(float(line[6]))
            struc = GeomBASE.read_geom(data)
            dxline = data[data.str.contains(r'^\s+MAX DISPLAC\.')].index.to_numpy(dtype=int)
            dmline = data[data.str.contains(r'^\s+RMS DISPLAC\.')].index.to_numpy(dtype=int)
            maxd = float(data.loc[dxline[-1]].strip().split()[2])
            rmsd = float(data.loc[dmline[-1]].strip().split()[2])
        return e, de, struc, maxg, rmsg, maxd, rmsd


class PhononBASE():
    """
    A container of basic methods for phonon information.
    """
    @classmethod
    def readmode_basic(cls, data, IRREP):
        """
        Read basic frequency information.

        Args:
            data (Series): Pandas series. The block containing frequency info.
            IRREP (list[str]): Irreducible representations in Mulliken symbols.
                Used for phonon dispersion only.

        Returns:
            frequency (array[float]): nmode \* 1
            mode_symm (array): nmode \* 1
            intens (array[float]): nmode \* 1
            IR (array[bool]): nmode \* 1
            Raman (array[bool]): nmode \* 1
        """
        import pandas as pd
        import numpy as np
        import re
        from CRYSTALpytools.units import cm_to_thz

        mode1 = data.map(lambda x: x[0:5]).to_numpy(dtype=int)
        mode2 = data.map(lambda x: x[6:10]).to_numpy(dtype=int)
        # in cm-1 for more decimal places
        freq_tmp = data.map(lambda x: x[24:36]).to_numpy(dtype=float)
        symm_tmp = data.map(lambda x: x[49:52].strip()).tolist()
        # for phonon dispersions, convert indices into symbols
        if len(re.findall(r'[A-Z,a-z]', symm_tmp[0])) == 0:
            if len(IRREP) == 0: symm_tmp = []
            symm_tmp = [IRREP[int(i)-1] for i in symm_tmp]

        # Raman and IR
        if 'I' in data[data.index[0]] or 'A' in data[data.index[0]]: # Raman and IR
            intens_tmp = data.map(lambda x: x[59:68]).to_numpy(dtype=float)
            IR_tmp = data.map(lambda x: x[56]=='A').to_numpy(dtype=bool)
            Raman_tmp = data.map(lambda x: x[73]=='A').to_numpy(dtype=bool)
        else:
            intens_tmp = []; IR_tmp = []; Raman_tmp = []
        # repeat
        frequency = []; mode_symm=[]; intens = []; IR = []; Raman = []
        count = 0
        for m1, m2 in zip(mode1, mode2):
            frequency.append(freq_tmp[count])
            for i in range(m1, m2): frequency.append(freq_tmp[count])
            ## symmetry
            if len(symm_tmp) > 0:
                mode_symm.append(symm_tmp[count])
                for i in range(m1, m2): mode_symm.append(symm_tmp[count])
            ## Raman and IR
            if len(intens_tmp) > 0:
                intens.append(intens_tmp[count]); IR.append(IR_tmp[count]); Raman.append(Raman_tmp[count])
                for i in range(m1, m2): intens.append(intens_tmp[count]); IR.append(IR_tmp[count]); Raman.append(Raman_tmp[count])
            count += 1

        # cm to thz
        frequency = cm_to_thz(np.array(frequency, dtype=float))
        mode_symm = np.array([i for i in mode_symm])
        if len(intens) != 0:
            intens = np.array(intens, dtype=float)
            IR = np.array(IR, dtype=bool)
            Raman = np.array(Raman, dtype=bool)
        return frequency, mode_symm, intens, IR, Raman

    @classmethod
    def readmode_eigenvector(cls, data, nmode):
        """
        Get mode eigenvectors.

        Returns:
            eigvt (array[float]): nmode\*natom\*3 array.
        """
        import pandas as pd
        import numpy as np

        img = data[data.str.contains(r'MODES IN ANTI\-PHASE')].index
        data = data[data.str.contains(r'[X,Y,Z]\s+\-*[0-9]\.')]
        nlines = len(data)

        if len(img) == 1: # complex
            rvt = data[0:int(nlines/2)].map(lambda x: x[13:].strip().split()).tolist()
            ivt = data[int(nlines/2):].map(lambda x: x[13:].strip().split()).tolist()
        else: # real
            rvt = data.map(lambda x: x[13:].strip().split()).tolist()
            ivt = []

        nblock = int(np.ceil(nmode / 6))
        nline = int(len(rvt) / nblock)

        rvttmp = np.array(rvt[0:nline], dtype=float)
        for i in range(1, nblock):
            rvttmp = np.hstack(
                [rvttmp, np.array(rvt[int(nline*i):int(nline*(i+1))], dtype=float)]
            )
        rvttmp = np.reshape(rvttmp.transpose(), [nmode, int(nline/3), 3])

        if len(ivt) != 0:
            ivttmp = np.array(ivt[0:nline], dtype=float)
            for i in range(1, nblock):
                ivttmp = np.hstack(
                    [ivttmp, np.array(ivt[int(nline*i):int(nline*(i+1))], dtype=float)]
                )
            ivttmp = np.reshape(ivttmp.transpose(), [nmode, int(nline/3), 3])
            eigvt = rvttmp + ivttmp * 1j
        else:
            eigvt = rvttmp

        return eigvt

    @classmethod
    def classical_amplitude(cls, struc, freq):
        """
        Get classical amplitude of phonon modes *Under Testing*

        .. math::

            x = \\sqrt{\\frac{\\hbar}{\\mu\\omega}}

        Args:
            struc (Structure): Pymatgen structure
            freq (float|array): Frequency. Unit: THz
        Returns:
            classic_a (array): nfreq\*3natom\*3natom array, or 3natom\*3natom
                if ``freq`` is float. The diagonal matrix of classical amplitude.
        """
        from CRYSTALpytools.units import amu_to_me, thz_to_hartree
        import numpy as np

        if type(freq) == float:
            freq = np.array([freq])

        natom = struc.num_sites
        nmode = int(3*natom)
        mass_rev = np.zeros([nmode, nmode]) # In AU
        for i in range(0, nmode, 3):
            atid = i // 3
            atmass = amu_to_me(float(struc.species[atomid].atomic_mass))
            mass_rev[i, i] = 1 / np.sqrt(atmass)
            mass_rev[i+1, i+1] = 1 / np.sqrt(atmass)
            mass_rev[i+2, i+2] = 1 / np.sqrt(atmass)

        freq_rev = 1 / np.sqrt(thz_to_hartree(freq))

        classic_a = []
        for f in freq_rev:
            classic_a.append(f * mass_rev)
        classic_a = np.array(classic_a)
        if len(freq) == 1:
            classic_a = classic_a[0]

        return classic_a

    @classmethod
    def normalize_eigenvector(cls, eigvt, amplitude=1.):
        """
        Normalize the mode of eigenvectors.

        Args:
            eigvt (array[complex]): nmode\*natom\*3 array.
            amplitude (float): Amplitude of normalization
        Returns:
            eigvt (array[complex]): Normalized eigenvector.
        """
        import numpy as np

        for idx_m, m in enumerate(eigvt):
            eigvt[idx_m] = eigvt[idx_m] / np.linalg.norm(m) * amplitude

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
        if len(overlap) != 0: # overlap
            qp_cord = []; qp_weight = []; nmode = []; frequency = [];
            mode_symm = []; intens = []; IR = []; Raman = []; eigenvector = []
            for idx_q in range(crysout.nqpoint):
                if idx_q in overlap[:, 1]: continue
                qp_cord.append(crysout.qpoint[idx_q][0])
                qp_weight.append(crysout.qpoint[idx_q][1])
                nmode.append(crysout.nmode[idx_q])
                frequency.append(crysout.frequency[idx_q])
                if len(crysout.mode_symm) != 0:
                    mode_symm.append(crysout.mode_symm[idx_q])
                if len(crysout.intens) != 0:
                    intens.append(crysout.intens[idx_q])
                if len(crysout.IR) != 0:
                    IR.append(crysout.IR[idx_q])
                if len(crysout.Raman) != 0:
                    Raman.append(crysout.Raman[idx_q])
                if len(crysout.eigenvector) != 0:
                    eigenvector.append(crysout.eigenvector[idx_q])
            # Update attributes
            crysout.nqpoint = len(qp_cord)
            crysout.qpoint = [[qp_cord[i], qp_weight[i]] for i in range(crysout.nqpoint)]
            crysout.nmode = np.array(nmode, dtype=int)
            crysout.frequency = np.array(frequency, dtype=float)
            if len(crysout.mode_symm) != 0:
                crysout.mode_symm = np.array(mode_symm, dtype=float)
            if len(crysout.intens) != 0:
                crysout.intens = np.array(intens, dtype=float)
            if len(crysout.IR) != 0:
                crysout.IR = np.array(IR, dtype=bool)
            if len(crysout.Raman) != 0:
                crysout.Raman = np.array(Raman, dtype=bool)
            if len(crysout.eigenvector) != 0:
                crysout.eigenvector = np.array(eigenvector, dtype=float)
        return crysout

    @classmethod
    def clean_imaginary(cls, crysout, threshold):
        """
        Set negative frequenceies and related properteis to 0 and print warning
        message. Eigenvectors are kept.

        Args:
            crysout (Crystal_output): :code:`CRYSTALpytools.crystal_io.Crystal_output` object
            threshold (float): The threshold to identify a phonon mode as negative.
        """
        import numpy as np
        import warnings

        for q, freq in enumerate(crysout.frequency):
            neg_rank = np.where(freq <= threshold)[0]
            if len(neg_rank) == 0: continue
            # warnings.warn('Negative frequencies detected.\n  Calculated thermodynamics might be inaccurate. Negative frequencies will be substituted by 0.',
            #               stacklevel=2)

            if len(neg_rank) > 3:
                warnings.warn('MORE THAN 3 IMAGINARY MODES! The structure is highly probable to be unstable.', stacklevel=2)

            crysout.frequency[q, neg_rank] = 0.

            # Eigenvectors are kept.
            if len(crysout.mode_symm) != 0: crysout.mode_symm[q, neg_rank] = ''
            if len(crysout.intens) != 0: crysout.intens[q, neg_rank] = 0.
            if len(crysout.IR) != 0: crysout.IR[q, neg_rank] = False
            if len(crysout.Raman) != 0: crysout.Raman[q, neg_rank] = False

        return crysout


class POutBASE():
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
            raise FileNotFoundError('EXITING: an output file needs to be specified')

    def get_geometry(self):
        """
        Get geometry from properties output calculation.

        Returns:
            struc (CStructure): Modified Pymatgen structure
        """
        import re
        import numpy as np
        from pymatgen.core.lattice import Lattice
        from CRYSTALpytools.geometry import CStructure

        data = self.data
        countline = 0
        lattice = []
        cart_coord = []
        species = []
        ndimen = 0
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s+DIRECT LATTICE VECTOR COMPONENTS', line):
                lattice = [data[countline+1].strip().split(),
                           data[countline+2].strip().split(),
                           data[countline+3].strip().split()]
                countline += 4
                continue
            # get dimension, for 1-3D
            elif re.match(r'^\s+A\s+B\s+C\s+ALPHA\s+BETA\s+GAMMA\s+VOLUME', line):
                countline += 1
                line = data[countline].strip().split()
                [a, b, c, al, be, ga, vol] = [float(i) for i in line]
                s = a * b * np.sin(ga/180*np.pi)
                l = a
                if np.abs(vol-l) < 1e-4:
                    ndimen = 1
                elif np.abs(vol-s) < 1e-4:
                    ndimen = 2
                else:
                    ndimen = 3
            elif re.match(r'^\s+ATOM N\.AT\.\s+SHELL\s+X\(A\)', line):
                countline += 2
                line = data[countline]
                while not re.match(r'^\s*\*+\s*$', line):
                    line_data = line.strip().split()
                    species.append(line_data[1])
                    cart_coord.append(line_data[4:7])
                    countline += 1
                    line = data[countline]

                break
            else:
                countline += 1
                continue

        if len(cart_coord) == 0:
            raise Exception('Valid geometry not found.')

        pbc = {0 : (False, False, False),
               1 : (True, False, False),
               2 : (True, True, False),
               3 : (True, True, True)}
        if ndimen > 0:
            lattice = Lattice(np.array(lattice, dtype=float), pbc=pbc[ndimen])
        else:
            lattice = Lattice(np.eye(3)*500., pbc=(False, False, False))
        species = [int(i) for i in species]
        cart_coord = np.array(cart_coord, dtype=float)
        return CStructure(lattice=lattice, species=species, coords=cart_coord,
                          coords_are_cartesian=True)

    def get_lattice(self):
        """
        Get lattice matrix from properties output calculation. A 3D lattice is
        generated since no dimensionality information is provided.

        Returns:
            matrix (array): 3\*3 lattice matrix
        """
        struc = self.get_geometry()
        return struc.lattice.matrix

    def get_topond_geometry(self):
        """
        Get the cluster geometry and plot plane base (2D only) from TOPOND
        calculation output.

        Returns:
            atomsplt (array): Atomic numbers and coordinates in plotting frame.
            base (array): *Valid for 2D plots only* 3\*3 range of orthogonal
                plotting base x and y. A: (xmin, ymax), B: (xmin, ymin), C:
                (xmax, ymin). Unit: Bohr.
        """
        import re
        import numpy as np
        from CRYSTALpytools.units import angstrom_to_au
        from scipy.spatial.transform import Rotation

        data = self.data
        countline = 0
        istopond = False; atomsplt = []; rotmx = []; xyrange = [];
        while countline < len(data):
            line = data[countline]
            if re.match(r'^\s*\*\s+T O P O N D', line):
                istopond = True
                countline += 1
            elif re.match(r'^\s*\*\*\* ATOMS \(POINTS\)\: AT\. N\. AND TRASFORMED COORD\.\(AU\)',
                          line):
                countline += 2
                line = data[countline]
                while line.strip() != '':
                    atomsplt.append(line.strip().split())
                    countline += 1
                    line = data[countline]
            elif re.match(r'^\s*ROTAT\. MATRIX', line):
                for i in range(3):
                    line = data[countline+i]
                    rotmx.append(line[37:].strip().split()[0:3])
                countline += 3
            elif re.match(r'^\s*ORIGIN AT \( AU\; SYSTEM REF\. FRAME\)', line):
                origin = line[37:].strip().split()[0:3]
                origin = np.array(origin, dtype=float)
                countline += 1
            elif re.match(r'^\s*X AXIS RANGES AND INCREMENTS', line):
                xyrange.append(line[37:].strip().split()[0:3])
                countline += 1
                line = data[countline]
                xyrange.append(line[37:].strip().split()[0:3])
                xyrange = np.array(xyrange, dtype=float)
                if '(ANG)' in line:
                    xyrange = angstrom_to_au(xyrange)
                break
            else:
                countline += 1

        if istopond == False:
            raise Exception("TOPOND output not found. Is it a TOPOND output file?")

        atomsplt = np.array(atomsplt, dtype=float)
        # define rotation
        rotmx = np.array(rotmx, dtype=float)
        rot = Rotation.from_matrix(rotmx)
        originplt = rot.apply(origin)
        # force origin to 0
        baseplt = np.vstack([originplt, originplt, originplt])
        baseplt[0, 0:2] += [xyrange[0, 0], xyrange[1, 1]]
        baseplt[1, 0:2] += [xyrange[0, 0], xyrange[1, 0]]
        baseplt[2, 0:2] += [xyrange[0, 1], xyrange[1, 0]]
        base = rot.inv().apply(baseplt)

        return atomsplt, base

    def get_reciprocal_lattice(self):
        """
        Get reciprocal lattice matrix from properties output calculation. A 3D
        lattice is generated since no dimensionality information is provided.

        Returns:
            matrix (array): 3\*3 reciprocal lattice matrix
        """
        struc = self.get_geometry()
        return struc.lattice.reciprocal_lattice.matrix

    def get_3dkcoord(self):
        """
        BANDS calculation only. Get 3D fractional coordinates of high-symmetry
        and sampled k points from output file.

        Returns:
            tick_pos3d (array): ntick\*3 array of fractional coordinates of
                high symmetry k points
            k_pos3d(array): nkpoint\*3 fractional coordinates of k points
        """
        import re
        import numpy as np

        data = self.data
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

    def get_XRDSPEC(self):
        """
        The keyword 'XRDSPEC' only. Get calculated XRD spectra.
        """
        import pandas as pd
        import numpy as np

        df = pd.DataFrame(self.data)
        title = df[df[0].str.contains(r'^\s+XRD SPECTRUM\s+$')].index
        if len(title) == 0: raise Exception("XRD spectra not found in file.")

        end = df[df[0].str.contains(r'^\s*T+ XRDSPEC\s+TELAPSE')].index
        if len(end) == 0: raise Exception("Abnormal termination. XRD spectra output broken.")

        spec = df[0][title[0]+7 : end[0]].map(lambda x: x.strip().split()).tolist()
        spec = np.array(spec, dtype=float)
        return spec

    def get_Fermi(self):
        """
        Get Fermi energy in eV from the common block.
        """
        import pandas as pd
        from CRYSTALpytools.units import H_to_eV

        df = pd.DataFrame(self.data)
        fline = df[df[0].str.contains(r'^\s*N\. OF SCF CYCLES.+FERMI ENERGY')].index[0]
        return H_to_eV(float(df[0].loc[fline].strip().split()[-1]))


