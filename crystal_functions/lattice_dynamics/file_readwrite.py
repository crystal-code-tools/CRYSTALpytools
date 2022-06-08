#!/usr/bin/env python
# coding: utf-8

from crystal_functions.file_readwrite import Crystal_output
from crystal_functions.convert import cry_out2pmg


class Freq_output(Crystal_output):
    """
    Class Freq_output, inheriated from Crystal_out and with thermodynamic-
    specific attributes, including:
        self.lattice: __init__, Lattice information
        self.edft: get_edft, DFT total energy, with probable corrections.
                   Unit: KJ / mol cell
        self.nqpoint: get_qpoint, Number of q points
        self.qpoint: get_qpoint, Fractional coordinates of qpoints
        self.nmode: get_mode, Number of vibrational modes at all qpoints
        self.frequency: get_mode, Frequencies of all modes at all qpoints.
                        Unit: THz
        self.eigenvector: get_eigenvector, Eigenvectors (classical amplitude)
                          of all atoms, all modes at all qpoints.
                          Unit: Angstrom

        By Spica. Vir., ICL. Jun. 01, 2022.
    """

    def __init__(self, output_name):
        """
        Input:
            The name of '.out' file
        Output:
            self, Crystal_functions output object
            self.lattice, pymatgen structure object, lattice and atom
                          information
        """
        super(Freq_output, self).__init__(output_name)
        self.lattice = cry_out2pmg(self, initial=False, vacuum=500)

    def get_edft(self):
        '''
        Get the DFT total energy of simulation cell. Unit: KJ / mol cell. To
        include probable energy corrections, the value at 'CENTRAL POINT' of
        force constant matrix is adopted.
        Input:
            -
        Output:
            self.edft, float, DFT total energy. Unit: KJ / mol cell
        '''
        import re

        for i, line in enumerate(self.data):
            if re.match(r'\s*CENTRAL POINT', line):
                self.edft = float(line.strip().split()[2]) * 2625.500256
                break

        return self.edft

    def get_qpoint(self):
        """
        Get the qpoints at which the phonon frequency is calculated.
        Input:
            -
        Output:
            self.nqpoint, int, Number of q points where the frequencies are
                          calculated.
            self.qpoint, nq * 3 numpy float array, Fractional coordinates of
                         qpoints.
        """
        import numpy as np
        import re

        self.nqpoint = 0
        self.qpoint = np.array([], dtype=float)

        for i, line in enumerate(self.data):
            if re.search(r'EXPRESSED IN UNITS\s*OF DENOMINATOR', line):
                shrink = int(line.strip().split()[-1])

            if re.match(r'\s*DISPERSION K POINT NUMBER', line):
                coord = np.array(line.strip().split()[7:10], dtype=float)
                self.qpoint = np.append(self.qpoint, coord / shrink)
                self.nqpoint += 1

        self.qpoint = np.reshape(self.qpoint, (-1, 3))
        if self.nqpoint == 0:
            self.nqpoint = 1
            self.qpoint = np.array([0, 0, 0], dtype=float)

        return self.nqpoint, self.qpoint

    def get_mode(self):
        """
        Get corresponding vibrational frequencies and for all modes and
        compute the total number of vibration modes (natoms * 3).

        Input:
            -
        Output:
            self.nmode, nqpoint * 1 numpy int array, Number of vibration modes
                        at each qpoints.
            self.frequency: nqpoint * nmode numpy float array, Harmonic
                            vibrational frequency. Unit: THz
        """
        import numpy as np
        import re

        if not hasattr(self, 'nqpoint'):
            self.get_qpoint()

        self.frequency = np.array([], dtype=float)

        countline = 0
        while countline < len(self.data):
            is_freq = False
            if re.match(r'\s*DISPERSION K POINT NUMBER\s*\d',
                        self.data[countline]):
                countline += 2
                is_freq = True

            if re.match(r'\s*MODES\s*EIGV\s*FREQUENCIES\s*IRREP',
                        self.data[countline]):
                countline += 2
                is_freq = True

            while self.data[countline].strip() and is_freq:
                line_data = re.findall(r'\-*[\d\.]+[E\d\-\+]*',
                                       self.data[countline])
                if line_data:
                    nm_a = int(line_data[0].strip('-'))
                    nm_b = int(line_data[1])
                    freq = float(line_data[4])

                for mode in range(nm_a, nm_b + 1):
                    self.frequency = np.append(self.frequency, freq)

                countline += 1

            countline += 1

        self.frequency = np.reshape(self.frequency, (self.nqpoint, -1))
        self.nmode = np.array([len(i) for i in self.frequency], dtype=float)

        return self.nmode, self.frequency

    def get_eigenvector(self):
        """
        Get corresponding mode eigenvectors for all modes on all
        atoms in the supercell.

        Input:
            -
        Output:
            self.eigenvector, nqpoint * nmode * natom * 3 numpy float array,
                              Eigenvectors expressed in Cartesian coordinate,
                              at all atoms, all modes and all qpoints.
                              Classical amplitude. Unit: Angstrom
        """
        import numpy as np
        import re

        if not hasattr(self, 'nmode'):
            self.get_mode()

        total_mode = np.sum(self.nmode)
        countline = 0
        # Multiple blocks for 1 mode. Maximum 6 columns for 1 block.
        if np.max(self.nmode) >= 6:
            countmode = 6
        else:
            countmode = total_mode

        # Read the eigenvector region as its original shape
        block_label = False
        total_data = []
        while countline < len(self.data) and countmode <= total_mode:
            # Gamma point / phonon dispersion calculation
            if re.match(r'\s*MODES IN PHASE', self.data[countline]) or\
               re.match(r'\s*NORMAL MODES NORMALIZED', self.data[countline]):
                block_label = True
            elif re.match(r'\s*MODES IN ANTI-PHASE', self.data[countline]):
                block_label = False

            # Enter a block
            if re.match(r'\s*FREQ\(CM\*\*\-1\)', self.data[countline]) and\
               block_label:
                countline += 2
                block_data = []
                while self.data[countline].strip():
                    # Trim annotation part (12 characters)
                    line_data = re.findall(r'\-*[\d\.]+[E\d\-\+]*',
                                           self.data[countline][13:])
                    if line_data:
                        block_data.append(line_data)

                    countline += 1

                countmode += len(line_data)
                total_data.append(block_data)

            countline += 1

        total_data = np.array(total_data, dtype=float)

        # Rearrage eigenvectors
        block_per_q = len(total_data) / self.nqpoint
        self.eigenvector = []
        # 1st dimension, nqpoint
        for q in range(self.nqpoint):
            index_bg = int(q * block_per_q)
            index_ed = int((q + 1) * block_per_q)
            q_data = np.hstack([i for i in total_data[index_bg: index_ed]])
        # 2nd dimension, nmode
            q_data = np.transpose(q_data)
        # 3rd dimension, natom
            natom = len(self.lattice.sites)
            q_rearrange = [np.split(m, natom, axis=0) for m in q_data]

            self.eigenvector.append(q_rearrange)

        self.eigenvector = np.array(self.eigenvector) * 0.529177

        return self.eigenvector


def clean_imaginary(self):
    """
    Substitute imaginary modes and corresponding eigenvectors with numpy NaN
    format and print warning message.

    Input:
        -
    Output:
        cleaned attributes.
        self.frequency
        self.eigenvector
    """
    import numpy as np

    for q, freq in enumerate(self.frequency):
        if freq[0] > -1e-4:
            continue

        print('WARNING: Negative frequencies detected - Calculated thermodynamics might be inaccurate.')
        print('WARNING: Negative frequencies will be substituted with NaN.')

        neg_rank = np.where(freq <= -1e-4)[0]
        self.frequency[q, neg_rank] = np.nan

        natom = len(self.lattice.sites)
        nan_eigvt = np.full([natom, 3], np.nan)
        self.eigenvector[q, neg_rank] = nan_eigvt

    return self.nmode, self.frequency, self.eigenvector
