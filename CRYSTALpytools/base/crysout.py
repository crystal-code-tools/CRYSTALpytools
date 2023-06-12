#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods of keywords used in 'crystal' output file.
"""

class PhononBASE():
    """
    A container of basic methods.

    Whether the output is a frequency one is checked during initialization. The
    identifier is :code:`+++ SYMMETRY ADAPTION OF VIBRATIONAL MODES +++`

    Args:
        data (Crystal_output.data): The ``data`` attribute of
            ``CRYSTALpytools.crystal_io.Crystal_output`` object, i.e., a list
            of text.
    """
    def __init__(self, data):
        import re

        is_freq = False

        for line in data:
            if re.match(r'^\s*\+\+\+\sSYMMETRY\sADAPTION\sOF\sVIBRATIONAL\sMODES\s\+\+\+', line):
                is_freq = True
                break
            else:
                continue

        if is_freq == False:
            raise Exception('Not a frequency calculation.')

    def readmode_basic(self, data):
        """
        Read :math:`E_{0}` energies, q points and frequency information..

        Returns:
            self.edft (float): Energy (in kJ/mol) reported in 'CENTERAL POINT'
                line (DFT + corrected energy)
            self.nqpoint (int): Number of q points
            self.qpoint (list[list[array[float], float]]): A nqpoint list of
                2\*1 list whose first element is a 3\*1 array of q point
                fractional coordinates and the second is its weight.
            self.nmode (array[int]): Number of modes at q point.
            self.frequency (array[float]): nqpoint \* nmode array of vibrational
                frequency. Unit: THz
            self.intens (array[float]): nqpoint \* nmode array of harmonic
                intensiy. Unit: km/mol
            self.IR (array[bool]): nqpoint \* nmode array of boolean values
                specifying whether the mode is IR active
            self.Raman (array[bool]): nqpoint \* nmode array of boolean values
                specifying whether the mode is Raman active
        """
        import re
        import numpy as np
        from CRYSTALpytools.units import H_to_kjmol

        self.edft = 0.
        self.nqpoint = 0
        self.qpoint = []
        self.nmode = []
        self.frequency = []
        self.intens = []
        self.IR = []
        self.Raman = []

        countline = 0
        has_spec = False
        while countline < len(data):
            line = data[countline]
            # E_0 with empirical corrections
            if re.match(r'^\s+CENTRAL POINT', line):
                self.edft = H_to_kjmol(float(line.strip().split()[2]))
                countline += 1
                continue

            # Q point info + frequency
            ## Dispersion
            elif re.match(r'^.+EXPRESSED IN UNITS\s+OF DENOMINATOR', line):
                shrink = int(line.strip().split()[-1])
                countline += 1
                continue

            elif re.match(r'\s+DISPERSION K POINT NUMBER', line):
                coord = np.array(line.strip().split()[7:10], dtype=float)
                weight = float(line.strip().split()[-1])
                self.qpoint.append([coord / shrink, weight])
                self.nqpoint += 1
                ### Read phonons
                countline += 2
                while countline < len(data):
                    line = data[countline]
                    if re.match(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+IRREP', line):
                        countline += 2
                        continue
                    elif re.match(r'^\s+\(HARTREE\*\*2\)\s+\(CM\*\*\-1\)', line):
                        countline += 1
                        continue
                    elif line.strip() == '':
                        break
                    else:
                        line_data = line.split()
                        nm_a = int(line_data[0].strip('-'))
                        nm_b = int(line_data[1])
                        freq = float(line_data[4])
                        # IR/Raman analysis, closed by default in dispersion calcs
                        if 'A' in line or 'I' in line:
                            has_spec = True
                            intens = float(line_data[-2].strip(')')) / (nm_b - nm_a + 1)
                            IR = line_data[-4] == 'A'
                            Raman = line_data[-1] == 'A'

                        for mode in range(nm_a, nm_b + 1):
                            self.frequency.append(freq)
                            if has_spec:
                                self.intens.append(intens)
                                self.IR.append(IR)
                                self.Raman.append(Raman)
                        countline += 1

                countline += 1
                continue
            ## Gamma point
            elif re.match(r'^\s+MODES\s+EIGV\s+FREQUENCIES\s+IRREP', line) and self.nqpoint == 0:
                countline += 2
                while countline < len(data):
                    line = data[countline]
                    if line.strip() == '':
                        countline += 1
                        break
                    line_data = line.split()
                    nm_a = int(line_data[0].strip('-'))
                    nm_b = int(line_data[1])
                    freq = float(line_data[4])
                    # IR/Raman analysis, closed by default in dispersion calcs
                    if 'A' in line or 'I' in line:
                        has_spec = True
                        intens = float(line_data[-2].strip(')')) / (nm_b - nm_a + 1)
                        IR = line_data[-4] == 'A'
                        Raman = line_data[-1] == 'A'

                    for mode in range(nm_a, nm_b + 1):
                        self.frequency.append(freq)
                        if has_spec:
                            self.intens.append(intens)
                            self.IR.append(IR)
                            self.Raman.append(Raman)
                    countline += 1

                countline += 1
                continue
            ## Other data
            else:
                countline += 1
                continue

        # HA/QHA Gamma point calculation
        if self.nqpoint == 0:
            self.nqpoint = 1
            self.qpoint = [[np.zeros([3,]), 1.]]
        else:
            for i in range(self.nqpoint): # Update weight
                self.qpoint[i][1] /= self.nqpoint

        self.frequency = np.reshape(np.array(self.frequency), [self.nqpoint, -1])
        self.nmode = np.array([len(i) for i in self.frequency], dtype=int)
        if has_spec:
            self.intens = np.reshape(self.intens, (self.nqpoint, -1))
            self.IR = np.reshape(self.IR, (self.nqpoint, -1))
            self.Raman = np.reshape(self.Raman, (self.nqpoint, -1))

        return self

    def readmode_eigenvector(self, data):
        """
        Get mode eigenvectors.

        Returns:
            self.eigenvector (array[float]): nqpoint\*nmode\*natom\*3 array of
                eigenvectors. Normalized to 1.
        """
        import numpy as np
        import re

        total_mode = np.sum(self.nmode)
        countline = 0
        # Multiple blocks at 1 qpoint. Maximum 6 columns for 1 block.
        if np.max(self.nmode) >= 6:
            countmode = 6
        else:
            countmode = total_mode

        # Read the eigenvector region as its original shape
        block_label = False
        total_data = []
        while countline < len(data) and countmode <= total_mode:
            line = data[countline]
            # Gamma point / phonon dispersion calculation
            if re.match(r'^\s+MODES IN PHASE', line) or re.match(r'^\s+NORMAL MODES NORMALIZED', line):
                block_label = True
                countline += 1
                continue
            elif re.match(r'^\s+MODES IN ANTI-PHASE', line):
                block_label = False
                countline += 1
                continue

            # Found a block
            elif re.match(r'^\s+FREQ\(CM\*\*\-1\)', line) and block_label:
                countline += 2
                block_data = []
                while data[countline].strip() != '':
                    line = data[countline]
                    # Trim annotation part (12 characters)
                    line_data = re.findall(r'\-*[\d\.]+[E\d\-\+]*', line[13:])
                    if line_data:
                        block_data.append(line_data)
                    countline += 1

                countmode += len(line_data)
                total_data.append(block_data)
                countline += 1
                continue
            # Other lines
            else:
                countline += 1
                continue

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
            natom = int(self.nmode[q] / 3)
            q_rearrange = [np.split(m, natom, axis=0) for m in q_data]
            self.eigenvector.append(q_rearrange)

        self.eigenvector = np.array(self.eigenvector)

        # Normalize eigenvectors of each mode to 1
        for idx_q, q in enumerate(self.eigenvector):
            for idx_m, m in enumerate(q):
                self.eigenvector[idx_q, idx_m] = \
                    self.eigenvector[idx_q, idx_m] / np.linalg.norm(m)

        return self


    def clean_q_overlap(self):
        """
        Remove the repeated q points at both ends of line segment when
        dispersion is read. The weight of q points will be updated here.
        """
        import numpy as np
        import warnings

        if self.nqpoint <= 1:
            return self

        overlap = []
        for idx_q1 in range(self.nqpoint - 1):
            qvec1 = self.qpoint[idx_q1][0]
            for idx_q2 in range(idx_q1 + 1, self.nqpoint): # q1 < q2
                qvec2 = self.qpoint[idx_q2][0]
                if np.linalg.norm(qvec1 - qvec2) < 1e-4:
                    warnings.warn('Overlap of q points is detected between q points {:3d} and {:3d}'.format(idx_q1, idx_q2),
                                  stacklevel=2)
                    overlap.append([idx_q1, idx_q2])

        overlap = np.array(overlap, dtype=int)
        weight = np.array([i[1] for i in self.qpoint])
        for idx_q in range(self.nqpoint, 1, -1):
            for idx_o in np.where(overlap[:, 1] == idx_q)[0]:
                self.nqpoint -= 1
                del self.qpoint[overlap[idx_o, 1]]
                weight = np.delete(weight, overlap[idx_o, 1])
                self.nmode = np.delete(self.nmode, overlap[idx_o, 1])
                self.frequency  = np.delete(self.frequency, overlap[idx_o, 1], axis=0)
                if len(self.intens) != 0:
                    self.intens = np.delete(self.intens, overlap[idx_o, 1], axis=0)
                    self.IR = np.delete(self.IR, overlap[idx_o, 1], axis=0)
                    self.Raman = np.delete(self.Raman, overlap[idx_o, 1], axis=0)
                if hasattr(self, 'eigenvector'):
                    if len(self.eigenvector) != 0:
                        self.eigenvector = np.delete(self.eigenvector, overlap[idx_o, 1], axis=0)
                break # Avoid repeatly deleting

        # Update weight
        weight = weight * 1 / np.sum(weight)
        for i in range(self.nqpoint):
            self.qpoint[i][1] = weight[i]

        return self

    def clean_imaginary(self, threshold=-1e-3):
        """
        Substitute imaginary modes and corresponding eigenvectors with numpy
        NaN format and print warning message.

        Args:
            threshold: The threshold to identify a phonon mode as negative.
        """
        import numpy as np
        import warnings

        for q, freq in enumerate(self.frequency):
            if np.isnan(freq[0]) or freq[0] > threshold:
                continue

            warnings.warn('Negative frequencies detected.\n  Calculated thermodynamics might be inaccurate. Negative frequencies will be substituted by NaN.',
                          stacklevel=2)

            neg_rank = np.where(freq <= threshold)[0]
            self.frequency[q, neg_rank] = np.nan

            if hasattr(self, 'eigenvector'):
                if len(self.eigenvector) != 0:
                    natom = int(self.nmode[q] / 3)
                    nan_eigvt = np.full([natom, 3], np.nan)
                    self.eigenvector[q, neg_rank] = nan_eigvt

            if len(self.intens) != 0:
                self.intens[q, neg_rank] = np.nan
                self.IR[q, neg_rank] = False
                self.Raman[q, neg_rank] = False

        return self