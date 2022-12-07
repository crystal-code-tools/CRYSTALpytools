#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A comprehensive module for lattice dynamics based on harmonic and quasiharmonic
approximations.
"""
from crystal_functions.file_readwrite import Crystal_output


class Mode:
    """
    Class Mode - store important information for a given vibrational mode. Can
    be used for a single harmonic phonon calculation and multiple calculations
    on the same system.

    Initialization:
        self.rank, __init__, The rank of the mode object. Start from 1.
        self.ncalc, __init__, The number of harmonic frequency calculations.
        self.frequency, __init__, Frequencies of the mode. Unit: THz
        self.volume, __init__, Cell volumes of harmonic calculations.
                     Unit: Angstrom^3
        self.eigenvector, __init__, Corresponding eigenvectors. Unit: Angstrom

    Limited to ncalc = 1 cases:
        self.zp_energy, get_zp_energy, Zero point energy of the mode.
                        Unit: KJ/mol cell
        self.U_vib, get_U_vib, Vibration contribution to internal energy,
                    including zero-point energy. Unit: KJ/mol cell
        self.entropy, get_entropy, Entropy of the mode. Unit: J/mol cell*K
        self.C_v, get_C_v, Constant volume specific heat. Unit: J/mol cell*K

    Limited to ncalc > 1 cases:
        self.poly_fit, polynomial_fit, dictionary of numpy polynomial objects.
                       Key: orders of power, Value: fitted polynomials
        self.poly_fit_rsqaure, polynomial_fit, dictionary of the goodness o
                               fittings, characterized by R^2.
    """

    def __init__(self, rank=0, frequency=[], volume=[], eigenvector=[]):
        """
        Input:
            rank, int, The rank of the mode object, from 1.
            frequency, ncalc * 1 array / list, Frequencies of the mode.
                       Unit: THz. Note: Not angular frequency, which is
                       frequency * 2pi
            volume, ncalc * 1 array / list, Lattice volumes of harmonic
                    calculations. Unit: Angstrom^3
            eigenvector, ncalc * natom * 3 array / list, Corresponding
                         eigenvectors. Unit: Angstrom
        Output:
            self.rank
            self.ncalc, int, The number of harmonic calculations
            self.frequency
            self.volume
            self.eigenvector
        """
        import numpy as np

        self.rank = rank
        self.ncalc = len(frequency)
        self.frequency = np.array(frequency, dtype=float)
        self.volume = np.array(volume, dtype=float)
        self.eigenvector = np.array(eigenvector, dtype=float)

    def get_zp_energy(self):
        """
        Get the zero-point energy of a single mode. Limited to ncalc = 1 cases.

        Input:
            -
        Output:
            self.zp_energy, float, Zero-point energy of a given mode.
                            Unit: KJ/mol cell
        """
        import numpy as np

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        hbar_freq = self.frequency[0] * 6.022141 * 6.626070E-2
        self.zp_energy = 0.5 * hbar_freq

        return self.zp_energy

    def get_U_vib(self, temperature=298.15):
        """
        Get the vibration contribution to internal energy of a single mode.
        Limited to ncalc = 1 cases. U_vib includes zero-point energy.

        Input:
            temperature: float, the temperature where the thermal contribution
                         is computed.
        Output:
            self.U_vib, float, Vibration contribution to internal energy of a
                        given mode. Unit: KJ/mol cell
        """
        import numpy as np

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        if not hasattr(self, 'zp_energy'):
            self.get_zp_energy()

        if temperature == 0:
            self.U_vib = self.zp_energy
            return self.U_vib

        hbar_freq = self.frequency[0] * 6.022141 * 6.626070E-2
        kb_t = 1.380649E-3 * 6.022141 * temperature
        expon = np.exp(hbar_freq / kb_t)
        self.U_vib = self.zp_energy + hbar_freq / (expon - 1)

        return self.U_vib

    def get_entropy(self, temperature):
        """
        Get the entropy of a single mode. Limited to ncalc = 1 cases.

        Input:
            temperature: float, the temperature where the thermal contribution
                         is computed.
        Output:
            self.entropy, float, The entropy of a given mode.
                          Unit: J/mol cell*K
        """
        import numpy as np

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        if temperature == 0:
            self.entropy = 0
            return self.entropy

        hbar_freq = self.frequency[0] * 6.022141 * 6.626070E-2
        kb_t = 1.380649E-3 * 6.022141 * temperature
        expon = np.exp(hbar_freq / kb_t)
        entS = kb_t * (hbar_freq / kb_t / (expon - 1) - np.log(1 - 1 / expon))
        self.entropy = entS / temperature * 1000

        return self.entropy

    def get_C_v(self, temperature):
        """
        Get the constant volume specific heat of a single mode. Limited to
        ncalc = 1 cases.

        Input:
            temperature: float, the temperature where the thermal contribution
                         is computed.
        Output:
            self.C_v, float, The constant volume specific heat of a given mode.
            Unit: J/mol cell*K
        """
        import numpy as np

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        if temperature == 0:
            self.C_v = 0
            return self.C_v

        hbar_freq = self.frequency[0] * 6.022141 * 6.626070E-2
        kb_t = 1.380649E-3 * 6.022141 * temperature
        expon = np.exp(hbar_freq / kb_t)

        self.C_v = hbar_freq**2 / kb_t / \
            temperature * expon / (expon - 1)**2 * 1000

        return self.C_v

    def polynomial_fit(self, eq_point, order=[2, 3]):
        """
        Fit phonon frequency as the polynomial function of volume with
        perturbation theory. Limited to ncalc > 1 cases.

        Input:
            eq_point, int, The DFT total energy minimum point with equilibrium
                           volume
            order, norder * 1 list, The orders of polynomials to be fitted.
        Output:
            self.poly_fit, norder * 1 dictionary, the dictionary of numpy
                           polynomial objects. Key: orders of power, Value:
                           fitted polynomials
            self.poly_fit_rsquare, norder * 1 dictionary, the dictionary of the
                                   goodness of fittings, characterized by R^2.
        """
        import numpy as np
        import warnings

        if self.ncalc <= 1:
            raise Exception(
                'This modulus is limited to multiple frequency calculations.')

        if max(order) > self.ncalc - 1:
            warnings.warn(
                'Reference data not sufficient for the order of polynomial fitting.')
            warnings.warn('Too high values will be removed.')

        order = list(set(order))
        order = [p for p in order if p <= self.ncalc - 1]

        self.poly_fit = {}
        self.poly_fit_rsqaure = {}

        vol_fit = self.volume - self.volume[eq_point]
        freq_fit = self.frequency - self.frequency[eq_point]

        for i in order:
            func = np.polynomial.polynomial.Polynomial.fit(
                vol_fit, freq_fit, i)
            self.poly_fit.update({i: func})
            if np.all(abs(self.frequency) < 1E-4):
                r_square = 1.
            else:
                ss_res = np.sum((self.frequency - func(self.volume))**2)
                ss_tot = np.sum((self.frequency - np.mean(self.frequency))**2)
                r_square = 1 - ss_res / ss_tot

            self.poly_fit_rsqaure.update({i: r_square})

        return order, self.poly_fit, self.poly_fit_rsqaure


class Harmonic(Crystal_output):
    """
    Class Harmonic, inherited from the Crystal_output class, with thermodynamic
    attributes. Used for harmonic phonon calclulations. Harmonic object can be
    defined by either a harmonic phonon calculation output file or manually set
    the all the information (usually for QHA).

    Inherited attributes:
        self.edft (length should be 1)
        self.nqpoint
        self.qpoint
        self.nmode
        self.frequency
        self.eigenvector

    New attributes:
        self.strucrure, The pymatgen Structure object of the calculation cell.
        self.mode, The list of mode objects.
        self.temperature, The temperature range for HA thermodynamics.
        self.pressure, The pressure range for HA thermodynamics.
        self.helmholtz, Helmholtz free energy of the cell, at HA level.
                        Unit: KJ/mol cell
        self.gibbs, Gibbs free energy of the cell, at HA level.
                    Unit: KJ/mol cell
        * Following attributes are the summed up values of corresponding
          attributes of class Mode.
        self.U_vib, thermodynamics
        self.zp_energy, thermodynamics
        self.entropy, thermodynamics
        self.C_v, thermodynamics
    """

    def __init__(self, temperature=[], pressure=[],
                 write_out=True, filename='HA-thermodynamics.dat'):
        """
        Initialization.

        Input:
            temperature, nTempt * 1 array / list, Temperatures where the
                         thermodynamic properties are computed. Unit: K
            pressure, npress * 1 array / list, Pressures where the
                      thermodyanmic properties are calculated. Unit: GPa
            write_out, bool, Wheter to print out HA thermodynamic properties.
            filename, str, Name of the printed-out file, used only if
                      write_out = True.
        Note: Temperature can also be defined in 'thermodynamics' method, which
              will cover the settings during initialisation.
        Output: 
            -
        """
        import numpy as np

        if len(temperature) > 0:
            self.temperature = np.array(temperature, dtype=float)

        if len(pressure) > 0:
            self.pressure = np.array(pressure, dtype=float)

        self.write_out = write_out
        if self.write_out:
            self.filename = filename
        else:
            self.filename = 'no file'

    def from_file(self, output_name, scelphono=[], read_eigenvector=False,
                  auto_calc=True):
        """
        Generate the Harominc object from a HA file.

        Input:
            output_name, str, Name of the output file.
            scellphono, ndimension * ndimension list, Supercell manipulation
                        matrix used for vibrational properties. Should be the
                        same as 'SCELPHONO' keyword in CRYSTAL input file. 3x3
                        list is also allowed.
            read_eigenvector, bool, Whether reading eigenvectors, i.e.,
                              normalised modes from outputs.
            auto_calc, bool, Whether to automatically launch thermodynamic
                       calculations. Parameters defined during initialization
                       will be used.
        Output:
            self.strucrure, pymatgen Structure object, the calculation cell,
                            reduced by SCELPHONO.
            self.natom, int, The number of atoms in the reduced calculation cell.
            self.volume, float, The volume the reduced calculation cell.
            self.mode, nqpoint * nmode array, List of mode objects at all
                       qpoints.
        """
        import numpy as np
        from crystal_functions.thermodynamics import Mode

        if hasattr(self, "volume"):
            raise Exception(
                "Data exists. Cannot overwriting the existing data.")

        super(Harmonic, self).read_cry_output(output_name)
        super(Harmonic, self).get_mode()
        super(Harmonic, self).clean_imaginary()
        self._generate_structure(scelphono=scelphono)

        if len(self.edft) != 1:
            raise Exception(
                "Only a single frequency calculation is premitted.")
        else:
            self.edft = self.edft[0]

        # Transfer the modes in self.freqency into lists of mode objects
        self.mode = []
        for qpoint in self.frequency:
            qmode = []
            for m, mode in enumerate(qpoint):
                qmode.append(
                    Mode(rank=m + 1, frequency=[mode], volume=[self.volume]))

            self.mode.append(qmode)

        if read_eigenvector:
            super(Harmonic, self).get_eigenvector()

        if auto_calc:
            self.thermodynamics(sumphonon=True)
            self.print_results()

        return self

    def from_data(self, edft, mode, **geometry):
        """
        Generate a Harmonic object by specifying data. Not recommanded to be
        used as a standalone method.

        Input:
            edft: float, Electron total energy
            mode: nqpoint * nmode array, List of mode objects.
            structure: Optional, Pymatgen structure object
            natom: Optional, int, number of atoms
            volume: Optional, float, volume of the simulation cell.
                    Unit Angstrom^3
        Output:
            self.nqpoint, int, Number of qpoints.
            self.nmode, nqpoint * 1 array, Number of modes at each q point.
            self.mode, nqpoint * nmode array, List of mode objects at Gamma.
            self.structure, pymatgen Structure object, the calculation cell,
                            reduced by SCELPHONO.
            self.natom, int, Number of atoms in the reduced calculation cell.
            self.volume, float, The volume the reduced calculation cell.
        """
        import numpy as np

        if hasattr(self, "volume"):
            raise Exception(
                "Data exists. The current command will be ignored.")

        for key, value in geometry.items():
            if key == 'structure':
                self.structure = value
                self.natom = len(value.species)
                self.volume = value.lattice.volume
                break
            elif key == 'natom':
                self.natom = int(value)
            elif key == 'volume':
                self.volume = float(value)

        self.edft = edft
        self.mode = mode
        self.nqpoint = len(mode)
        self.nmode = np.array([len(q) for q in self.mode])

        return self

    def _generate_structure(self, scelphono):
        """
        Eliminate the influences of the keyword 'SCELPHONO' and generate the
        pymatgen 'structure' object of the actual cell used for phonon
        calculation. Not a standalone method.

        Input:
            scellphono, ndimension * ndimension list. 3*3 list is also allowed.
        Output:
            self.structure, pymatgen Structure object. The calculation cell
                            without 'SCELPHONO' expansion.
            self.natom, int, Number of atoms in the cell.
            self.volume, float, pymatgen structure.volume. Cell volume.
                         Unit: Angstrom^3
        """
        from crystal_functions.convert import cry_out2pmg
        from pymatgen.core.structure import Structure
        import numpy as np

        ndimen = self.get_dimensionality()

        if not scelphono or not ndimen:
            self.structure = cry_out2pmg(self, vacuum=100)
            self.natom = len(self.structure.species)
            self.volume = self.structure.volume

            return self

        scell_mx = np.eye(3, dtype=float)
        scell_mx[: ndimen, : ndimen] = np.array(scelphono)[: ndimen, : ndimen]
        shrink_mx = np.linalg.pinv(scell_mx)

        scel = cry_out2pmg(self, vacuum=100)

        pcel_lattice = np.dot(scel.lattice.matrix, shrink_mx)
        all_coord = np.dot(scel.cart_coords,
                           np.linalg.pinv(pcel_lattice)).tolist()
        all_species = scel.species

        pcel_coord = []
        pcel_species = []

        for i, coord in enumerate(all_coord):
            if any(x > 0.5 or x <= -0.5 for x in coord):
                continue
            else:
                pcel_coord.append(coord)
                pcel_species.append(all_species[i])

        pcel_natom = len(pcel_species)
        pcel_charge = int(scel.charge * pcel_natom / len(all_species))

        self.structure = Structure(
            pcel_lattice, pcel_species, pcel_coord, pcel_charge)
        self.natom = len(self.structure.species)
        self.volume = self.structure.volume

        return self

    def _phonon_sumup(self, temperature, calculate_zp):
        """
        Summing up inidival phonon modes at each q point. Translational modes
        with frequencies = 0 are skipped.

        Not a standalone method. For thermodynamics, use 'self.thermodyanmics'
        instead.

        Input:
            temperature, float, The temperature where the U_vib, entropy and
                         C_v are calculated.
            calculate_zp, bool, Calculate zero-point energy or temperature
                          dependent properties.
        Output:
            zp_energy, nqpoint * 1 array, Zero-point energy at a given q point.
                       Returned if temperature = None.
            U_vib, nqpoint * 1 array, Vibrational contribution to internal
                   energy at constant temperature and all q points. Returned
                   if temperature is given.
            entropy, nqpoint * 1 array, Entropy at constant temperature and
                     given q point. Returned if temperature is given.
            C_v, nqpoint * 1 array, Constant volume specific heat at constant
                 temperature and given q point. Returned if temperature is
                 given.
        """
        import numpy as np

        if calculate_zp:
            zp_energy = np.array([], dtype=float)
        else:
            T = temperature
            U_vib = np.array([], dtype=float)
            entropy = np.array([], dtype=float)
            C_v = np.array([], dtype=float)

        for qpoint in self.mode:
            if calculate_zp:
                zp_energy_q = 0.
            else:
                U_vib_q = 0.
                entropy_q = 0.
                C_v_q = 0.
            # Remove the translational modes
            for mode in qpoint:
                if np.isnan(mode.frequency) or mode.frequency <= 1e-5:
                    continue

                if calculate_zp:
                    zp_energy_q += mode.get_zp_energy()
                else:
                    U_vib_q += mode.get_U_vib(temperature=T)
                    entropy_q += mode.get_entropy(temperature=T)
                    C_v_q += mode.get_C_v(temperature=T)

            if calculate_zp:
                zp_energy = np.append(zp_energy, zp_energy_q)
            else:
                U_vib = np.append(U_vib, U_vib_q)
                entropy = np.append(entropy, entropy_q)
                C_v = np.append(C_v, C_v_q)

        if calculate_zp:
            return zp_energy
        else:
            return U_vib, entropy, C_v

    def thermodynamics(self, sumphonon=True, mutewarning=False, **temptpress):
        """
        Calculate the thermodynamic properties (zp_energy, U_vib, entropy, C_v
        and Helmholtz free energy) of the given system, at all qpoints and the
        whole temperature range.

        Input:
            temperature, Optional, nTempt * 1 array / list, Temperatures where
                         the thermodynamic properties are computed. Unit: K
            pressure, Optional, npress * 1 array / list, Pressures where the
                      thermodyanmic properties are calculated. Unit: GPa
            sumphonon, bool, Whether summing up the phonon contributions across
                       the first Brillouin zone and take average.
            mutewarning, bool, Whether print out warning messages of updating
                               temperature and pressure.
        Output:
            self.helmholtz, nqpoint * nTempt numpy array, Helmholtz free
                            energy. Unit: KJ/mol cell
            self.gibbs, nqpoint * nTempt * npress numpy array, Gibbs free 
                        energy. Unit: KJ/mol cell
            self.zp_energy, nqpoint * 1 numpy array, Zero-point energy.
                            Unit: KJ/mol cell
            self.U_vib, nqpoint * nTempt numpy array, Vibrational contribute to
                        internal energy. Unit: KJ/mol cell
            self.entropy, nqpoint * nTempt numpy array, Entropy.
                          Unit: J/mol cell*K
            self.C_v, nqpoint * nTempt numpy array, Constant volume specific
                      heat. Unit: J/mol cell*K

        Generated, not returned output:
            self.temperature, nTempt * 1 array / list, Temperature range.
                              Unit: K
            self.pressure, nPress * 1 array / list, Pressure range.
                              Unit: K

        Note: If sumphonon = True, the thermodyanmic properties of phonon 
              dispersion are summed. In this case, nqpoint = 1 but the
              dimension is kept.
        """
        import warnings
        import numpy as np

        # Generate temperature and pressure series
        if temptpress:
            if 'temperature' in temptpress:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn(
                        'Temperature attribute exists. Input temperatures will be used to update the attribute.')

                self.temperature = np.array(
                    temptpress['temperature'], dtype=float)

            if 'pressure' in temptpress:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn(
                        'Pressure attribute exists. Input pressures will be used to update the attribute.')

                self.pressure = np.array(temptpress['pressure'], dtype=float)
        else:
            if not hasattr(self, 'temperature') or \
               not hasattr(self, 'pressure'):
                raise ValueError(
                    'Temperature and pressure should be specified.')

        self.zp_energy = self._phonon_sumup(temperature=0., calculate_zp=True)
        U_vib = []
        entropy = []
        C_v = []
        helmholtz = []
        gibbs = []

        for T in self.temperature:
            gibbs_t = []
            U_vib_t, entropy_t, C_v_t = self._phonon_sumup(
                temperature=T, calculate_zp=False)
            helm_t = -entropy_t * T / 1000 + U_vib_t + self.edft

            for p in self.pressure:
                gibbs_tp = p * self.volume * 0.602214 + helm_t
                gibbs_t.append(gibbs_tp)

            # nTemp * nqpoint
            U_vib.append(U_vib_t)
            entropy.append(entropy_t)
            C_v.append(C_v_t)
            helmholtz.append(helm_t)
            # nTemp * npress * nqpoint
            gibbs.append(gibbs_t)

        U_vib = np.transpose(np.array(U_vib, dtype=float))
        entropy = np.transpose(np.array(entropy, dtype=float))
        C_v = np.transpose(np.array(C_v, dtype=float))
        helmholtz = np.transpose(np.array(helmholtz, dtype=float))
        gibbs = np.transpose(np.array(gibbs, dtype=float), (2, 0, 1))

        if sumphonon:
            self.U_vib = np.array([np.sum(U_vib, axis=0)]) / self.nqpoint
            self.entropy = np.array([np.sum(entropy, axis=0)]) / self.nqpoint
            self.C_v = np.array([np.sum(C_v, axis=0)]) / self.nqpoint
            self.helmholtz = np.array(
                [np.sum(helmholtz, axis=0)]) / self.nqpoint
            self.gibbs = np.array([np.sum(gibbs, axis=0)]) / self.nqpoint
            self.nqpoint = 1
            self.qpoint = np.array([[0., 0., 0.]], dtype=float)
        else:
            self.U_vib = U_vib
            self.entropy = entropy
            self.C_v = C_v
            self.helmholtz = helmholtz
            self.gibbs = gibbs

        return self.helmholtz, self.gibbs, self.zp_energy, self.U_vib,\
            self.entropy, self.C_v

    def print_results(self):
        """
        Print a single output file for HA thermodynamics. Used if write_out=True.

        Note: Phonon dispersions are forced to be summed to keep the output
        succinct if the automatic scheme (write_out=True) is launched. To get
        verbose outputs, directly use this attribute.

        Input:
            -
        Output:
            filename, text file.
        """
        if not self.write_out:
            print('Harmonic.write_out = False, return to empty.')
            return

        file = open(self.filename, 'w')
        file.write('%-21s%12.4e%-15s%12.4e%-10s\n' %
                   ('# DFT TOTAL ENERGY = ', self.edft / 96.485340,
                    ' eV,         = ', self.edft, ' kJ/mol'))
        file.write('%-21s%12.6f%-15s%12.6f%-10s\n' %
                   ('# CELL VOLUME      = ', self.volume,
                    ' Angstrom^3, = ', self.volume * 0.602214, ' cm^3/mol'))
        file.write('%s\n' % '# LATTICE PARAMETERS (ANGSTROM, DEGREE)')
        file.write('%12s%12s%12s%12s%12s%12s\n' % ('A', 'B', 'C',
                                                   'ALPHA', 'BETA', 'GAMMA'))
        file.write('%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n\n' %
                   (self.structure.lattice.parameters[0:6]))

        for q in range(self.nqpoint):
            file.write('%-40s%5i\n\n' %
                       ('# HARMONIC THERMODYNAMICS AT QPOINT #', q))
            file.write('%-23s%12.6e%8s\n\n' %
                       ('  zero point energy = ', self.zp_energy[q], ' kJ/mol'))
            file.write('%s\n' % '  temperature dependent properties')
            file.write('%8s%18s%18s%18s%18s\n' %
                       ('    T(K)', 'U_vib(kJ/mol)', 'Entropy(J/mol*K)',
                        'C_V(J/mol*K)', 'Helmholtz(kJ/mol)'))
            for t, tempt in enumerate(self.temperature):
                file.write('%8.2f%18.6e%18.6e%18.6e%18.6e\n' %
                           (tempt, self.U_vib[q, t], self.entropy[q, t],
                            self.C_v[q, t], self.helmholtz[q, t]))

            file.write('\n')
            file.write('%s\n' % '  Gibbs free energy')
            file.write('%-30s' % '    rows    : pressure (GPa)  ')
            for p in self.pressure:
                file.write('%8.3f%1s' % (p, ''))

            file.write('\n')
            file.write('%-30s' % '    columns : temperature (K) ')
            for t in self.temperature:
                file.write('%8.3f%1s' % (t, ''))

            for gibbs_t in self.gibbs[q]:
                file.write('\n%4s' % '')
                for gibbs_p in gibbs_t:
                    file.write('%18.6e%2s' % (gibbs_p, ''))

            file.write('\n\n\n')
            file.close()

        return


class Quasi_harmonic:
    """
    Class Quasi_haromic - Generate and arrange harmonic phonons, store the
    fitted, volume dependent QHA phonon information and obtain the QHA
    thermodynamic properties.

    self.ncalc, from_HA_files, The number of phonon calculations.
    self.combined_phonon, from_HA_files, Sampled calculation as Harmonic object
    self.combined_volume, from_HA_files, A list of volumes.
    self.combined_edft, from_HA_files, A list of DFT total energies.
    self.combined_mode, from_HA_files, A list of mode objects.
    self.eos_method, edft_eos_fit, Fitting method of equation of states
    self.eos, edft_eos_fit, Fitted equation of states
    self.fit_order, freq_polynomial_fit, The optimal order of polynomial fit
    self.temerature, thermodynamics, Temperature series. Unit: K
    self.pressure, thermodynamics, Pressure series. Unit: GPa
    self.equilibrium_volume, thermodynamics, V(T, p). Unit: Angstrom^3
    self.helmholtz, thermodynamics, F(T, V). Unit: kJ/mol
    self.gibbs, thermodynamics, G(T, p). Unit: kJ/mol
    self.entropy, thermodynamics, S(T, V). Unit: J/mol*K
    """

    def __init__(self, temperature=[], pressure=[],
                 write_out=True, filename='QHA-Fit.dat'):
        """
        Initialization.

        Input:
            temperature, nTempt * 1 array / list, Temperatures where the
                         thermodynamic properties are computed. Unit: K
            pressure, npress * 1 array / list, Pressures where the
                      thermodyanmic properties are calculated. Unit: GPa
            write_out, bool, Whether to record the key information into a file.
            filename, str, Name of the printed-out file, used only if
                      write_out = True.
        Note: Temperature can also be defined in 'thermodynamics' method, which
              will cover the settings during initialisation.
        Output: 
            -
        """
        import numpy as np

        if len(temperature) > 0:
            self.temperature = np.array(temperature, dtype=float)

        if len(pressure) > 0:
            self.pressure = np.array(pressure, dtype=float)

        self.write_out = write_out
        if self.write_out:
            self.filename = filename
        else:
            self.filename = 'no file'

    def from_HA_files(self, input_files, scelphono=[], overlap=0.4, sort_phonon=True):
        """
        Read data from individual HA calculation outputs.

        Input:
            input_files, ncalc*1 list, List of phonon output filenames.
            scelphono, ndimen*ndimen or 3*3 list / array, Same to the
                       'SCELPHONO' keyword of CRYSTAL17 input.
            overlap, float, The threshold of close mode overlaps
            sort_phonon, bool, Whether to check phonon continuity. 
        Output:
            self.ncalc, int, Number of HA phonon calculations.
            self.combined_phonon, self.combined_volume, self.combined_edft,
            self.combined_mode, refer the method 'combine_data'
        """
        from crystal_functions.thermodynamics import Harmonic
        import warnings

        if hasattr(self, "ncalc"):
            warnings.warn('Data exists. The current command will be ignored.')

            return self

        self.ncalc = len(input_files)
        if self.ncalc == 1:
            warnings.warn(
                'Single frequency calculation detected! QHA is deteriorated to HA.')

        ha_list = [
            Harmonic(write_out=False).from_file(
                file,
                scelphono=scelphono,
                read_eigenvector=True,
                auto_calc=False
            ) for file in input_files
        ]

        self.combined_phonon, self.combined_volume, self.combined_edft, \
            self.combined_mode = self._combine_data(ha_list, overlap=overlap,
                                                    sort_phonon=sort_phonon)

        return self

    def from_QHA_file(self, input_file, scelphono=[], overlap=0.4, sort_phonon=True):
        """
        Read data from a single QHA calculation at Gamma point.

        Input:
            input_file, string or 1*1 list, Name of QHA output file. Only 1
                        file is permitted

        Other Input/output are consistent with 'from_HA_files'
        """
        from crystal_functions.file_readwrite import Crystal_output
        from crystal_functions.thermodynamics import Harmonic
        from crystal_functions.thermodynamics import Mode
        import warnings
        import re
        import numpy as np
        from pymatgen.core import Structure

        if hasattr(self, "ncalc"):
            warnings.warn('Data exists. The current command will be ignored.')
            return self

        if isinstance(input_file, list) and len(input_file) > 1:
            raise Exception("Only a single QHA file is permitted")
        elif isinstance(input_file, list) and len(input_file) == 1:
            input_file = input_file[0]

        file = Crystal_output().read_cry_output(input_file)
        file.get_mode()
        file.get_eigenvector()
        file.clean_imaginary()

        # Get volume/structure/dimensionality. Only to be used with QHA files
        structures = []
        for idx_line, line in enumerate(file.data):
            if re.match(
                r'^\s+GEOMETRY\sFOR\sWAVE\sFUNCTION\s\-\sDIMENSIONALITY', line
            ):
                ndimen = int(line.strip().split()[9])
            elif re.match(
                r'^\s+DIRECT\sLATTICE\sVECTORS\sCARTESIAN\sCOMPONENTS\s\(ANGSTROM\)',
                line
            ):
                idx_line += 2
                vec1 = np.array(file.data[idx_line].strip().split()[
                                0:3], dtype=float)
                vec2 = np.array(
                    file.data[idx_line + 1].strip().split()[0:3], dtype=float)
                vec3 = np.array(
                    file.data[idx_line + 2].strip().split()[0:3], dtype=float)

                idx_line += 9
                all_species = []
                all_coords = np.array([], dtype=float)
                while re.match(
                    r'^\s+[0-9]+\s+[0-9]+\s+[A-Z]+', file.data[idx_line]
                ):
                    line_info = file.data[idx_line].strip().split()
                    all_coords = np.append(
                        all_coords, np.array(line_info[3:], dtype=float))
                    all_species.append(line_info[2].capitalize())
                    idx_line += 1

                all_coords = np.reshape(all_coords, [-1, 3])

                scell_mx = np.eye(3, dtype=float)
                if scelphono:
                    scell_mx[: ndimen, : ndimen] = np.array(
                        scelphono)[: ndimen, : ndimen]
                    shrink_mx = np.linalg.pinv(scell_mx)
                    pcel_lattice = np.dot(
                        np.stack([vec1, vec2, vec3]), shrink_mx)
                    all_coords = np.dot(
                        all_coords, np.linalg.pinv(pcel_lattice)).tolist()

                    pcel_coord = []
                    pcel_species = []
                    for i, coord in enumerate(all_coords):
                        if any(x > 0.5 or x <= -0.5 for x in coord):
                            continue
                        else:
                            pcel_coord.append(coord)
                            pcel_species.append(all_species[i])

                else:
                    pcel_lattice = np.stack([vec1, vec2, vec3])
                    pcel_coord = all_coords
                    pcel_species = all_species

                struc = Structure(lattice=pcel_lattice, species=pcel_species,
                                  coords=pcel_coord, coords_are_cartesian=False)
                structures.append(struc)
            else:
                continue

        self.ncalc = file.nqpoint
        ha_list = []
        for idx_c in range(self.ncalc):
            ha = Harmonic(write_out=False)
            # The first one is pre-opt geom
            ha.structure = structures[idx_c + 1]
            ha.natom = len(ha.structure.species)
            ha.volume = ha.structure.lattice.volume
            ha.edft = file.edft[idx_c]
            ha.nqpoint = 1
            ha.qpoint = [0, 0, 0]
            ha.nmode = np.array([file.nmode[idx_c]])
            ha.frequency = np.array([file.frequency[idx_c]])
            ha.eigenvector = np.array([file.eigenvector[idx_c]])
            ha_list.append(ha)

        self.combined_phonon, self.combined_volume, self.combined_edft, \
            self.combined_mode = self._combine_data(ha_list, overlap=overlap,
                                                    sort_phonon=sort_phonon)

        return self

    def _combine_data(self, ha_list, overlap, sort_phonon):
        """
        Combine the HA calculation data and rearrange it according to modes.
        Not a standalone method.

        NOTE: All the input data will be rearranged in the low-to-high sequence
              according to volumes.

        Input:
            ha_list, ncalc * 1 list, The list of harmonic objects.
            overlap, float, The threshold of close mode overlaps
            sort_phonon, bool, Whether to check phonon continuity. 
        Output:
            combined_phonon, list of Harmonic objects, Sampled calculations
            combined_volume, ncalc * 1 list, A list of volumes.
                             Unit: Angstrom^3
            combined_edft, ncalc * 1 list, A list of DFT total energies.
                           Unit: KJ / mol cell
            combined_mode, nqpoint * nmode list, A list of mode objects. Each
                           mode object stands for a vibrational mode at the
                           given q point and stores ncalc HA values for volume,
                           frequency and eigenvector.
                           mode.volume: ncalc * 1 array
                           mode.frequency: ncalc * 1 array
                           mode.eigenvector: ncalc * natom * 3 array
        """
        import numpy as np
        import warnings
        from crystal_functions.thermodynamics import Mode

        # Sorting data according to volumes
        sorted_vol = np.zeros([self.ncalc, 2])
        nqpoint = ha_list[0].nqpoint
        nmode = ha_list[0].nmode  # nqpoint * 1 array
        natom = ha_list[0].natom  # int
        for index, ha_phonon in enumerate(ha_list):
            sorted_vol[index, :] = [index, ha_phonon.volume]
            # Check whether the numbers of modes and atoms are consistent.
            if (natom - ha_phonon.natom) != 0 or \
                    not np.all((nmode - ha_phonon.nmode) == 0) or \
                    nqpoint - ha_phonon.nqpoint != 0:
                raise Exception(
                    'The number of qpoints, modes or atoms is not consistent across the sampling points'
                )

        sorted_vol = sorted_vol[np.argsort(sorted_vol[:, 1])]
        nmode = nmode[0]

        combined_phonon = []
        # Volume, ncalc * 1 array
        combined_volume = np.zeros(self.ncalc)
        # DFT total energy, ncalc * 1 array
        combined_edft = np.zeros(self.ncalc)
        # Frequency, ncalc * nqpoint * nmode array
        combined_freq = np.zeros([self.ncalc, nqpoint, nmode])
        # Eigenvector, ncalc * nqpoint * nmode * natom * 3 array
        combined_eigvt = np.zeros([self.ncalc, nqpoint, nmode, natom, 3])
        for idx_new, idx_vol in enumerate(sorted_vol):
            ha_phonon = ha_list[int(idx_vol[0])]
            combined_phonon.append(ha_phonon)
            combined_volume[idx_new] = idx_vol[1]
            combined_edft[idx_new] = ha_phonon.edft
            combined_freq[idx_new] = ha_phonon.frequency
            combined_eigvt[idx_new] = ha_phonon.eigenvector

        # ncalc * nqpoint * nmode array to nqpoint * ncalc * nmode array
        combined_freq = np.transpose(combined_freq, axes=[1, 0, 2])
        # ncalc * nqpoint * nmode * natom * 3 array to
        # nqpoint * ncalc * nmode * natom * 3 array
        combined_eigvt = np.transpose(combined_eigvt, axes=[1, 0, 2, 3, 4])
        if sort_phonon:
            close_overlap = np.zeros([nqpoint, self.ncalc, nmode, nmode])
            for idx_q in range(nqpoint):
                combined_freq[idx_q], combined_eigvt[idx_q], close_overlap[idx_q] \
                    = self._phonon_continuity(combined_freq[idx_q],
                                              combined_eigvt[idx_q],
                                              overlap=overlap)
        # nqpoint * ncalc * nmode array to nqpoint * nmode * ncalc array
        combined_freq = np.transpose(combined_freq, axes=[0, 2, 1])
        # nqpoint * ncalc * nmode * natom * 3 array to
        # nqpoint *  nmode * ncalc * natom * 3 array
        combined_eigvt = np.transpose(combined_eigvt, axes=[0, 2, 1, 3, 4])
        if sort_phonon:
            # nqpoint * ncalc * nmode_ref * nmode_sort array to
            # nqpoint * nmode_ref * ncalc * nmode_sort array
            close_overlap = np.transpose(close_overlap, axes=[0, 2, 1, 3])
            for idx_q, qpoint in enumerate(close_overlap):
                overlap_numbers = np.sum(qpoint)
                if overlap_numbers >= 1.:
                    warnings.warn(
                        'Close overlap of phonon modes detected at qpoint: %3i, %6i overlaps out of %6i modes.'
                        % (idx_q, int(overlap_numbers), int(nqpoint * nmode)), stacklevel=2
                    )

        combined_mode = []
        for idx_q in range(nqpoint):
            combined_mode_q = []
            for idx_m in range(nmode):
                combined_mode_q.append(
                    Mode(rank=idx_m + 1,
                         frequency=combined_freq[idx_q, idx_m, :],
                         volume=combined_volume,
                         eigenvector=combined_eigvt[idx_q, idx_m, :])
                )

            combined_mode.append(combined_mode_q)

        if self.write_out:
            file = open(self.filename, 'w')
            file.write('%s\n' % '# COMBINED QHA DATA')
            file.write('%s' % '## SAMPLED VOLUMES(ANGSTROM^3) = ')
            for v in combined_volume:
                file.write('%16.4e' % v)

            file.write('\n')

            file.write('%s' % '## DFT TOTAL ENERGIES(KJ/MOL CELL) = ')
            for e in combined_edft:
                file.write('%16.6e' % e)

            file.write('\n\n')

            file.write('%s\n\n' % '## COMBINED MODES')
            for idx_q, qpoint in enumerate(combined_mode):
                file.write('%-27s%8i\n' %
                           ('### FREQUENCIES AT QPOINT #', idx_q))
                for mode in qpoint:
                    file.write('\n%-8s%22s%22s\n' %
                               ('  Mode #', 'Volume(Angstrom^3)', 'Frequency(THz)'))

                    for i in range(self.ncalc):
                        if i == 0:
                            file.write('%8i' % mode.rank)
                        else:
                            file.write('%8s' % '')

                        file.write('%22.4f%22.4f\n' %
                                   (mode.volume[i], mode.frequency[i]))

                file.write('\n')

            if sort_phonon:
                file.write('%s\n\n' %
                           '## CLOSE OVERLAPS OF PHONON FREQUENCIES')
                for idx_q, qpoint in enumerate(combined_mode):
                    file.write('%-30s%8i\n\n' %
                               ('### CLOSE OVERLAPS AT QPOINT #', idx_q))
                    file.write('%-10s%2s%8s%2s%9s%2s%9s\n' %
                               ('  Calc_Ref', '', 'Mode_Ref', '', 'Calc_Sort',
                                '', 'Mode_Sort'))
                    for idx_mref, mode in enumerate(qpoint):
                        if np.sum(close_overlap[idx_q, idx_mref]) < 1.:
                            continue

                        for idx_csort in range(1, self.ncalc):
                            for idx_msort in range(nmode):
                                if close_overlap[idx_q, idx_mref, idx_csort, idx_msort]:
                                    file.write(
                                        '%10i%2s%8i%2s%9i%2s%9i\n' % (
                                            idx_mref + 1, '', idx_csort - 1, '', idx_csort, '', idx_msort + 1)
                                    )
                                else:
                                    continue

                    file.write('\n')

            file.close()

        return combined_phonon, combined_volume, combined_edft, combined_mode

    def _phonon_continuity(self, freq, eigvt, symm=None, overlap=0.4):
        """
        Rearrange phonon modes by their continuity. If the difference between
        the maximum scalar product of correspondin eigenvectors (normalized to 
        1) and scalar products of other modes is less than 0.4, warning is
        printed due to the potential overlap of modes. Adopted from CRYSTAL17,

        Erba A. J. Chem. Phys., 2014 141 124115.

        Not a standalone method.

        Input:
            freq, ncalc * nmode array, Phonon frequencies.
            eigvt, ncalc * nmode * natom * 3 array, Eigenvectores normalized to
                   1 of corresponding modes
            symm, ncalc * nmode array, Sub-group numbers of corresponding modes
                  (Not supported at the current implementation)
            overlap, float, The threshold of close mode overlaps
        Output:
            freq, ncalc * nmode array, Sorted phonon frequencies
            eigvt, ncalc * nmode * natom * 3 array, Sorted eigenvectores
            close_overlap, ncalc * nmode * nmode boolian array, Whether close
                           overlap is identified at previous calculation (2nd 
                           dimension) and the current calculation (3rd).
        """
        import numpy as np

        # Exclude negative and 0 frequencies
        ncalc = len(freq)
        nmode = len(freq[0])
        ng_mode = 0
        for idx_c, calc in enumerate(freq):
            for idx_f, frequency in enumerate(calc):
                if np.isnan(frequency) or (frequency < 1e-4):
                    ng_mode_c = idx_f
                else:
                    break

            if ng_mode_c > ng_mode:
                ng_mode = ng_mode_c

        # Sort phonon
        products = np.zeros([ncalc, nmode])
        for sort_c in range(1, ncalc):
            ref_c = sort_c - 1
            for ref_m in range(ng_mode + 1, nmode):
                ref_pdt = 0.
                sort_m_save = 0
                for sort_m in range(ng_mode + 1, nmode):
                    if symm and symm[0, ref_m] != symm[sort_c, sort_m]:
                        continue

                    sort_pdt = abs(np.sum(
                        eigvt[ref_c, ref_m] * eigvt[sort_c, sort_m]
                    ))
                    if sort_pdt > ref_pdt:
                        if sort_m < ref_m:
                            check_pdt = abs(np.sum(
                                eigvt[ref_c, sort_m] * eigvt[sort_c, sort_m]
                            ))

                            if check_pdt > sort_pdt:
                                continue

                        ref_pdt = sort_pdt
                        sort_m_save = sort_m

                products[sort_c, ref_m] = ref_pdt
                freq[[sort_c, sort_c], [sort_m_save, ref_m]] \
                    = freq[[sort_c, sort_c], [ref_m, sort_m_save]]
                eigvt[[sort_c, sort_c], [sort_m_save, ref_m]] \
                    = eigvt[[sort_c, sort_c], [ref_m, sort_m_save]]
                if symm:
                    symm[[sort_c, sort_c], [sort_m_save, ref_m]] \
                        = symm[[sort_c, sort_c], [ref_m, sort_m_save]]

        # Look for close overlaps
        close_overlap = np.zeros([ncalc, nmode, nmode])
        for sort_c in range(1, ncalc):
            ref_c = sort_c - 1
            for ref_m in range(ng_mode + 1, nmode):
                ref_pdt = products[sort_c, ref_m]
                for sort_m in range(ng_mode + 1, nmode):
                    if symm and symm[0, ref_m] != symm[sort_c, sort_m]:
                        continue
                    if sort_m == ref_m:
                        continue

                    sort_pdt = abs(np.sum(
                        eigvt[ref_c, ref_m] * eigvt[sort_c, sort_m]
                    ))
                    if ref_pdt - sort_pdt < overlap:
                        close_overlap[ref_c, ref_m, sort_m] = 1

        return freq, eigvt, close_overlap

    def edft_eos_fit(self, method):
        """
        Fit electron total energy according to equation of states. Not a
        standalone method.

        Input:
            method: string, Name of EoS used. Consistent with requirements of
                    pymatgen (https://pymatgen.org/pymatgen.analysis.eos.html).
        Output:
            self.eos_method, string, Equation of State used
            self.eos, pymatgen EOS object, Fitted equation of state.
        """
        from pymatgen.analysis.eos import EOS

        self.eos_method = method
        self.eos = EOS(method).fit(self.combined_volume, self.combined_edft)
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s%s\n' % (
                '# EQUATION OF STATES FITTED FOR ELECTRON TOTAL ENERGY: ', method))
            file.write(
                '%s\n' % '  Electron total energy is fitted as the function of volume, of which the')
            file.write('%s\n\n' %
                       '  formalism is given by equation of states.')

            file.write('%16s%16s%12s%12s\n' %
                       ('E0(kJ/mol)', 'V0(Angstrom^3)', 'B0(GPa)', 'B1'))
            file.write('%16.4f%16.4f%12.4f%12.4f\n' % (self.eos.e0,
                                                       self.eos.v0,
                                                       self.eos.b0 * 1.660539,
                                                       self.eos.b1))
            file.write('\n')
            file.close()

        return self.eos

    def freq_polynomial_fit(self, order):
        """
        Fit phonon frequencies as polynomial functions of volumes. Not a
        standalone method.

        Input:
            order, list/array, List of the highest order of polynomials to be
                   fitted. Default: [2, 3] (quadratic, cubic)
        Output:
            self.fit_order, int, The optimal order of polynomial fit.

        Also see 'self.poly_fit' and 'self.poly_fit_rsquare' attributes of mode
        object
        """
        import numpy as np

        rsquare_tot = np.array([[od, 0] for od in order], dtype=float)
        eq_point = np.argmin(self.combined_edft)

        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# POLYNOMIAL FIT OF MODE FREQUENCY')
            file.write(
                '%s\n' % '  Frequency of each vibrational mode is fitted as the polynomial function of')
            file.write('%s\n' % '  volume, with specified orders of power.')

        for idx_q, mode_q in enumerate(self.combined_mode):
            rsquare_q = {od: 0. for od in order}

            if self.write_out:
                file.write('\n%s%8i\n' %
                           ('## POLYNOMIAL FIT AT QPOINT #', idx_q))

            for mode in mode_q:
                order_new, _, _ = mode.polynomial_fit(eq_point=eq_point,
                                                      order=order)
                for key, value in mode.poly_fit_rsqaure.items():
                    rsquare_q[key] += value / len(mode_q)

                if self.write_out:
                    file.write('%-8s%7s%14s%s\n' %
                               ('  Mode #', 'Order', 'R^2', '  Coeff low to high'))
                    for idx_od, od in enumerate(order_new):
                        if idx_od == 0:
                            file.write('%8i' % mode.rank)
                        else:
                            file.write('%8s' % '')

                        file.write('%7i%2s%12.6f%2s' %
                                   (od, '', mode.poly_fit_rsqaure[od], ''))
                        for c in mode.poly_fit[od].convert().coef:
                            file.write('%12.4e' % c)

                        file.write('\n')

                    file.write('\n')

            rsquare_tot[:, 1] += np.array(
                [rsquare_q[od] / len(self.combined_mode) for od in order]
            )

            if self.write_out:
                file.write('%s%8i\n' %
                           ('## POLYNOMIAL FIT GOODNESS AT QPOINT #', idx_q))
                file.write('%-7s%14s\n' % ('  Order', 'R^2'))
                for od in order:
                    file.write('%7i%2s%12.6f\n' % (od, '', rsquare_q[od]))

        self.fit_order = int(rsquare_tot[np.argmax(rsquare_tot[:, 1]), 0])

        if self.write_out:
            file.write('\n\n')
            file.close()

        return self

    def _get_harmonic_phonon(self, volume):
        """
        Get numerical phonon frequencies from fitted analytical expressions and
        generate harmonic phonon objects. Not a standalone method.

        Input:
            volume, float, The volume of harmonic lattice. Unit: Angstrom^3
        Output:
            ha, Harmonic, Harmonic phonon object with numerical data.
        """
        import numpy as np
        from crystal_functions.thermodynamics import Harmonic
        from crystal_functions.thermodynamics import Mode

        if not hasattr(self, 'fit_order') or not hasattr(self, 'eos'):
            raise Exception('ERROR: Analytical expressions unavailable.')

        eq_point = np.argmin(self.combined_edft)
        num_mode = []
        for mode_q in self.combined_mode:
            num_mode_q = []
            for idx_m, mode in enumerate(mode_q):
                d_vol = volume - self.combined_volume[eq_point]
                freq = mode.poly_fit[self.fit_order](d_vol) \
                    + mode_q[idx_m].frequency[eq_point]
                num_mode_q.append(
                    Mode(rank=idx_m + 1, frequency=[freq], volume=[volume])
                )

            num_mode.append(num_mode_q)

        ha = Harmonic(write_out=False).from_data(self.eos(volume),
                                                 num_mode, volume=volume)

        return ha

    def _minimize_gibbs(self, volume, temperature, pressure):
        """
        Get Gibbs free energy from the Harmonic phonon object. Used only for
        minimizing G(V; T, p) by SciPy. Not a standalone method.

        Input:
            volume, float, The volume of lattice (V). Unit: Angstrom^3
            temperature, float, T, argument. Unit: K
            pressure, float, p, argument. Unit: GPa
        """
        ha = self._get_harmonic_phonon(volume)
        ha.thermodynamics(temperature=[temperature], pressure=[pressure])

        return ha.gibbs[0, 0, 0]

    def thermodynamics(self, eos_method='birch_murnaghan', poly_order=[2, 3],
                       min_method='BFGS', volume_bound=None, mutewarning=False,
                       **temptpress):
        """
        1. Fit E_DFT and frequencies (if that has not been done) according to
        methods specified. 
        2. Calculate the 0 pressure equilibrium volume and pressure-independent
        properties (Helmholtz free energy, Entropy and Constant-volume specific
        heat) at given temperatures.
        3. Calculate pressure-dependent proerties (Gibbs free energy)

        Input:
            eos_method: string, Equation of state used to fit E_DFT. For EOSs
                        supported, refer https://pymatgen.org/pymatgen.analysis.eos.html
            poly_order: list/array, List of the highest order of polynomials to
                        be fitted.
            min_method: string, Minimisation algorithms. Parameterized and
                        tested algos: 
                        * BFGS(no boundary)
                        * L-BFGS-B(with boundary)
            volume_bound: turple-like, Boundary conditions of equilibrium
                          volumes. Unit: Angstrom^3
            mutewarning, bool, Whether print out warning messages.
            temperature: Optional, nTempt*1 list/array, Temperatures. Unit: K
            pressure: Optional, nPress*1 list/array, Pressures. Unit: GPa
        Output:
            self.temperature, nTempt*1 array, List of temperatures. Unit: K
            self.pressure, nPress*1 array, List of pressures. Unit: GPa
            self.equilibrium_volume, nPress*nTempt array, Equilibrium volumes
                                     at given temperature and pressure. Unit:
                                     Angstrom^3
            self.helmholtz, nPress*nTempt array, Helmholtz free energy at given
                            volume. Unit: kJ/mol
            self.gibbs, nPress*nTempt array, Gibbs free energy at given volume.
                        Unit: kJ/mol
            self.entropy, nPress*nTempt array, Entropy at given volume. Unit:
                          J/mol*K

        Optional outputs, see comments in edft_eos_fit, freq_polynomial_fit
        """
        import numpy as np
        import warnings
        from scipy.optimize import minimize

        # Generate temperature and pressure series
        if temptpress:
            if 'temperature' in temptpress:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn(
                        'Temperature attribute exists. Input temperatures will be used to update the attribute.', stacklevel=2)

                self.temperature = np.array(
                    temptpress['temperature'], dtype=float)

            if 'pressure' in temptpress:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn(
                        'Pressure attribute exists. Input pressures will be used to update the attribute.', stacklevel=2)

                self.pressure = np.array(temptpress['pressure'], dtype=float)

        if not hasattr(self, 'temperature') or not hasattr(self, 'pressure'):
            raise ValueError('Temperature and pressure should be specified.')

        # Fit DFT total energy, if not done yet. Otherwise, fitted values will not be covered.
        if hasattr(self, 'eos') and not mutewarning:
            warnings.warn(
                'DFT total energy is already fitted. To keep the consistency, it will not be updated.', stacklevel=2)
        else:
            self.edft_eos_fit(method=eos_method)

        # Fit frequencies, if not done yet. Otherwise, fitted values will not be covered.
        if hasattr(self, 'fit_order') and not mutewarning:
            warnings.warn(
                'Frequency is already fitted to polynomials. To keep the consistency, it will not be updated.', stacklevel=2)
        else:
            self.freq_polynomial_fit(order=poly_order)

        # Define minimization methods
        methods = {
            'BFGS': "vol = minimize(self._minimize_gibbs, v_init, args=(t, p), method='BFGS', jac='3-point')",
            'L-BFGS-B': "vol = minimize(self._minimize_gibbs, v_init, args=(t, p), method='L-BFGS-B', jac='3-point', bounds=volume_bound)",
        }

        # Gibbs(V; T, p) minimization nPress*nTempt list
        self.equilibrium_volume = []
        v_init = np.mean(self.combined_volume)

        for p in self.pressure:
            eq_vol_p = []
            for t in self.temperature:
                params = {'self': self,
                          'minimize': minimize,
                          'v_init': v_init,
                          't': t,
                          'p': p,
                          'volume_bound': volume_bound}
                exec(methods[min_method], params)
                eq_vol_p.append(params['vol'].x[0])

                if (params['vol'].x[0] < min(self.combined_volume)
                    or params['vol'].x[0] > max(self.combined_volume)) \
                   and not mutewarning:
                    warnings.warn(
                        'Optimised volume exceeds the sampled range. Special care should be taken of.',
                        stacklevel=2
                    )
                    warnings.warn(
                        '  Volume: %12.4f, Temperature: %6.2f, Pressure: %6.2f'
                        % (params['vol'].x[0], t, p), stacklevel=2
                    )

            self.equilibrium_volume.append(eq_vol_p)

        self.equilibrium_volume = np.array(self.equilibrium_volume)

        # Calculate other thermodynamic properties
        self.helmholtz = []
        self.gibbs = []
        self.entropy = []
        for idx_p, p in enumerate(self.pressure):
            helmholtz_p = []
            gibbs_p = []
            entropy_p = []
            for idx_t, t in enumerate(self.temperature):
                vol = self.equilibrium_volume[idx_p, idx_t]
                ha = self._get_harmonic_phonon(vol)
                ha.thermodynamics(temperature=[t], pressure=[
                                  p], mutewarning=True)
                helmholtz_p.append(ha.helmholtz[0, 0])
                gibbs_p.append(ha.gibbs[0, 0, 0])
                entropy_p.append(ha.entropy[0, 0])

            self.helmholtz.append(helmholtz_p)
            self.gibbs.append(gibbs_p)
            self.entropy.append(entropy_p)

        self.helmholtz = np.array(self.helmholtz)
        self.gibbs = np.array(self.gibbs)
        self.entropy = np.array(self.entropy)

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA THERMODYNAMIC PROPERTIES')
            file.write('%s\n\n' % '  Thermodynamic properties fitted by QHA.')
            file.write('%s%6i\n' %
                       ('## FREQUENCY POLYNOMIAL ORDER: ', self.fit_order))
            file.write('%s%s\n' %
                       ('## EQUILIBRIUM VOLUME MINIMISATION: ', min_method))
            if volume_bound:
                file.write('%s\n' % (
                    '## CONSTRAINED VOLUME MINIMIZATION LAUNCHED. VOLUME BOUNDARIES (UNIT: ANGSTROM^3):'))
                file.write('%s%8.2f%s%8.2f\n\n' % (
                    '## LOWER: ', volume_bound[0], ' UPPER: ', volume_bound[1]))

            for idx_p, press in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' %
                           ('## THERMODYNAMIC PROPERTIES AT ', press, '  GPa'))
                file.write('%4s%6s%4s%16s%2s%18s%4s%16s%4s%16s\n' %
                           ('', 'T(K)', '', 'Vol(Angstrom^3)', '', 'Helmholtz(kJ/mol)', '', 'Gibbs(kJ/mol)', '', 'Entropy(J/mol*K)'))
                for idx_t, tempt in enumerate(self.temperature):
                    file.write('%4s%6.1f%4s%16.4f%4s%16.8e%4s%16.8e%4s%16.8e\n' %
                               ('', tempt,
                                '', self.equilibrium_volume[idx_p, idx_t],
                                '', self.helmholtz[idx_p, idx_t],
                                '', self.gibbs[idx_p, idx_t],
                                '', self.entropy[idx_p, idx_t]))

                file.write('\n')

            file.write('\n')
            file.close()

        return self

    def eos_thermodynamics(self, eos_method='birch_murnaghan', poly_order=[2, 3],
                           mutewarning=False, **temptpress):
        """
        Obtain thermodynamic properties by fitting equation of states F(V). 
        Exact sorting and fitting of frequency-volume relationship is disabled,
        EoSs are fitted according to the Helmholtz free energy of sampled 
        harmonic phonons. 

        Input:
            temperature: Optional, nTempt*1 list/array, Temperatures. Unit: K
            pressure: Optional, nPress*1 list/array, Pressures. Unit: GPa
            eos_method: string, Equation of state used to fit F. For EOSs 
                        supported, refer https://pymatgen.org/pymatgen.analysis.eos.html
            poly_order: list/array, List of the highest order of polynomials to
                        be fitted (G(T) for entropy).
            mutewarning, bool, Whether print out warning messages.

        Output:
            self.thermo_eos_method: string, Equation of state used to fit F.
            self.thermo_eos: nTempt*1 list, List of adiabatic EoS at given
                             temperatures

        self.temperature, self.pressure, self.equilibrium_volume, 
        self.helmholtz, self.gibbs, self.entropy are consistent with method
        'thermodynamics'.
        """
        import numpy as np
        import warnings
        from crystal_functions.thermodynamics import Harmonic
        from pymatgen.analysis.eos import EOS
        from scipy.optimize import fmin
        from sympy import diff, lambdify, symbols

        # Check the number of calculations
        if self.ncalc < 4:
            raise Exception('Insufficient database. Increase HA phonons')

        # Generate temperature and pressure series
        if temptpress:
            if 'temperature' in temptpress:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn(
                        'Temperature attribute exists. Input temperatures will be used to update the attribute.',
                        stacklevel=2)

                self.temperature = np.array(
                    temptpress['temperature'], dtype=float)

            if 'pressure' in temptpress:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn(
                        'Pressure attribute exists. Input pressures will be used to update the attribute.',
                        stacklevel=2)

                self.pressure = np.array(temptpress['pressure'], dtype=float)

        if not hasattr(self, 'temperature') or not hasattr(self, 'pressure'):
            raise ValueError('Temperature and pressure should be specified.')

        # Data for fitting
        helmholtz = np.zeros([len(self.temperature), self.ncalc], dtype=float)
        for idx_c, calc in enumerate(self.combined_phonon):
            hfe_c, _, _, _, _, _ = calc.thermodynamics(
                sumphonon=True, mutewarning=True, temperature=self.temperature, pressure=[0.]
            )
            helmholtz[:, idx_c] = hfe_c

        # Fit EOS, get p-V-T
        self.thermo_eos_method = eos_method
        self.thermo_eos = []
        self.equilibrium_volume = np.zeros(
            [len(self.pressure), len(self.temperature)], dtype=float)
        self.helmholtz = np.zeros(
            [len(self.pressure), len(self.temperature)], dtype=float)
        self.gibbs = np.zeros(
            [len(self.pressure), len(self.temperature)], dtype=float)
        self.entropy = np.zeros(
            [len(self.pressure), len(self.temperature)], dtype=float)
        v = symbols('v')
        for idx_t, t in enumerate(self.temperature):
            hfe_t = EOS(eos_method).fit(self.combined_volume, helmholtz[idx_t])
            self.thermo_eos.append(hfe_t(v))
            press_t = -diff(hfe_t(v), v)
            for idx_p, p in enumerate(self.pressure):
                pau = p * 0.602214  # GPa --> kJ/mol.Angstrom^3
                lam_p = lambdify(v, (press_t - pau)**2, 'numpy')
                fit = fmin(lam_p, hfe_t.v0, full_output=True, disp=False)
                self.equilibrium_volume[idx_p, idx_t] = fit[0]
                self.helmholtz[idx_p, idx_t] = hfe_t(fit[0])
                self.gibbs[idx_p, idx_t] = hfe_t(fit[0]) + pau * fit[0]

                if (fit[0] < min(self.combined_volume) or fit[0] > max(self.combined_volume)) and not mutewarning:
                    warnings.warn(
                        'Optimised volume exceeds the sampled range. Special care should be taken of.', stacklevel=2)
                    warnings.warn('  Volume: %12.4f, Temperature: %6.2f, Pressure: %6.2f' % (
                        fit[0], t, p), stacklevel=2)

        # Second fit G(T; p), get entropy by S=-(\frac{\partial G}{\partial T})_{p}
        t = symbols('t')
        if max(poly_order) > len(self.temperature) - 1 and not mutewarning:
            warnings.warn(
                'Temperature series not sufficient for the order of polynomial fitting.', stacklevel=2)
            warnings.warn('Too high values will be removed.', stacklevel=2)

        poly_order = list(set(poly_order))
        poly_order = [p for p in poly_order if p <= len(self.temperature) - 1]

        for idx_p, gibbs in enumerate(self.gibbs):
            r_square = np.array([], dtype=float)
            func = []
            for order in poly_order:
                func_order = np.polynomial.polynomial.Polynomial.fit(
                    self.temperature, gibbs, order)
                func.append(func_order)
                ss_res = np.sum((gibbs - func_order(self.temperature))**2)
                ss_tot = np.sum((gibbs - np.mean(gibbs))**2)
                r_square = np.append(r_square, 1 - ss_res / ss_tot)

            fit_func = func[np.argmax(r_square)]
            entropy = -diff(fit_func(t), t) * 1000.
            for idx_t, tempt in enumerate(self.temperature):
                self.entropy[idx_p, idx_t] = entropy.evalf(subs={'t': tempt})

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA THERMODYNAMIC PROPERTIES - EOS FIT')
            file.write(
                '%s\n\n' % '  Thermodynamic properties obtained by overall fitting of equation of states.')
            file.write('%s%s\n' % ('## EQUATION OF STATES: ', eos_method))
            file.write('%s%i\n' % ('## G(T) POLYNOMIAL ORDER: ',
                                   poly_order[np.argmax(r_square)]))
            file.write(
                '%s\n' % '  WARNING: Entropy at low temperature is probably inaccurate due to the poor fitting of G(T) near 0K.')
            for idx_p, press in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' %
                           ('## THERMODYNAMIC PROPERTIES AT ', press, '  GPa'))
                file.write('%4s%6s%4s%16s%2s%18s%4s%16s%4s%16s\n' %
                           ('', 'T(K)', '', 'Vol(Angstrom^3)', '', 'Helmholtz(kJ/mol)', '', 'Gibbs(kJ/mol)', '', 'Entropy(J/mol*K)'))
                for idx_t, tempt in enumerate(self.temperature):
                    file.write('%4s%6.1f%4s%16.4f%4s%16.8e%4s%16.8e%4s%16.8e\n' %
                               ('', tempt,
                                '', self.equilibrium_volume[idx_p, idx_t],
                                '', self.helmholtz[idx_p, idx_t],
                                '', self.gibbs[idx_p, idx_t],
                                '', self.entropy[idx_p, idx_t]))

                file.write('\n')

            file.write('\n')
            file.close()

        return self
