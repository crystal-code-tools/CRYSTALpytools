#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A post-processing module for DFT lattice dynamics by harmonic and quasiharmonic
approximations (HA/QHA).
"""
from CRYSTALpytools.crystal_io import Crystal_output
from CRYSTALpytools import units


class Mode:
    """
    Store important information for a given vibrational mode and do analysis at
    mode-level. Not recommanded to be used individually.

    Args:
        rank (int): The rank of the mode object, from 1.
        frequency (array[float] | list[float]): Frequencies of the mode
            (Ncalc\*1). Unit: THz. Note: **NOT** angular frequency, which is 
            frequency * 2pi.
        volume (array[float] | list[float]): Lattice volumes of harmonic
            calculations (Ncalc\*1). Unit: Angstrom^3
        eigenvector (array[float] | list[float]): Corresponding normalized 
            eigenvectors (Ncalc\*Natom\*3).

    Returns:
        self.rank (int)
        self.ncalc (int): The number of harmonic calculations (typically at
            different volumes)
        self.frequency (array[float])
        self.volume (array[float])
        self.eigenvector (array[float])
    """

    def __init__(self, rank=0, frequency=[], volume=[], eigenvector=[]):
        import numpy as np

        self.rank = rank
        self.ncalc = len(frequency)
        self.frequency = np.array(frequency, dtype=float)
        self.volume = np.array(volume, dtype=float)
        self.eigenvector = np.array(eigenvector, dtype=float)

    def get_zp_energy(self):
        """
        Get the zero-point energy of a single mode with the following equation. 
        Limited to ncalc = 1 cases.

        .. math::

            E^{zp}_{i,\\mathbf{q}}=\\frac{1}{2}\\hbar\\omega_{i,\\mathbf{q}}

        Returns:
            self.zp_energy (float): Zero-point energy. Unit: KJ/mol

        :raise Exception: If ``self.ncalc`` > 1.
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e9
        self.zp_energy = 0.5 * hbar_freq

        return self.zp_energy

    def get_U_vib(self, temperature=298.15):
        """
        Get the vibration contribution to internal energy (including zero-point
        energy) of a single mode with the following equation. Limited to 
        ncalc = 1 cases.

        .. math::

            U^{vib}_{i,\\mathbf{q}}\\left(T\\right)=E^{zp}_{i,\\mathbf{q}}+
            \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{\\exp{\\left(
                \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T}
            \\right)}-1}

        Args:
            temperature (float, optional): Temperature where the quantity is 
                computed. Unit: K

        Returns:
            self.U_vib (float): Vibration contribution to internal energy. 
                Unit: KJ/mol

        :raise Exception: If ``self.ncalc`` > 1.
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        if not hasattr(self, 'zp_energy'):
            self.get_zp_energy()

        if temperature == 0:
            self.U_vib = self.zp_energy
            return self.U_vib

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e9
        kb_t = scst.k * scst.Avogadro * temperature * 1e-3
        expon = np.exp(hbar_freq / kb_t)
        self.U_vib = self.zp_energy + hbar_freq / (expon - 1)

        return self.U_vib

    def get_entropy(self, temperature):
        """
        Get the entropy of a single mode with the following equation. Limited
        to ncalc = 1 cases.

        .. math::

            S_{i,\\mathbf{q}}\\left(T\\right)=k_{B}\\left\\{
                \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T\\left[
                    \\exp{\\left(
                        \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T}
                    \\right)}-1
                \\right]}-\\ln{\\left[
                    1-\\exp{\\left(
                        -\\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T}
                    \\right)}
                \\right]}
            \\right\\}

        Args:
            temperature (float, optional): Unit: K

        Returns:
            self.entropy (float): Entropy. Unit: J/mol\*K

        :raise Exception: If ``self.ncalc`` > 1.
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        if temperature == 0:
            self.entropy = 0
            return self.entropy

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e12
        kb_t = scst.k * scst.Avogadro * temperature
        expon = np.exp(hbar_freq / kb_t)
        entS = kb_t * (hbar_freq / kb_t / (expon - 1) - np.log(1 - 1 / expon))
        self.entropy = entS / temperature

        return self.entropy

    def get_C_v(self, temperature):
        """
        Get the constant volume specific heat of a single mode with the 
        following equation. Limited to ncalc = 1 cases.

        .. math::

            C^{V}_{i,\\mathbf{q}}=
            \\frac{\\left(\\hbar\\omega_{i,\\mathbf{q}}\\right)^{2}}{k_{B}T^{2}}
            \\frac{\\exp{
            \\left(
                \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T}
            \\right)}
            }{\\left[
                \\exp{\\left(
                    \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T}
                \\right)-1}
            \\right]^{2}
            }

        Args:
            temperature (float, optional): Unit: K

        Returns:
            self.C_v (float): Constant volume specific heat. Unit: J/mol\*K

        :raise Exception: If ``self.ncalc`` > 1.
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise Exception(
                'This module is limited to a single frequency calculation.')

        if temperature == 0:
            self.C_v = 0
            return self.C_v

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e12
        kb_t = scst.k * scst.Avogadro * temperature
        expon = np.exp(hbar_freq / kb_t)

        self.C_v = hbar_freq**2 / kb_t / temperature * expon / (expon - 1)**2

        return self.C_v

    def polynomial_fit(self, eq_point, order=[2, 3]):
        """
        Fit phonon frequency as the polynomial function of volume with
        perturbation theory. Limited to ncalc > 1 cases.

        Args:
            eq_point (int): The index (not rank) corresponds to the DFT total
                energy minimum (typically the unstrained geometry).
            order (array[int] | list[int]], optional): Orders of
                polynomials used.

        Returns:
            self.poly_fit (Dict[int, NumPy Polynomial]): Key - orders of power, 
                Value - fitted NumPy polynomials
            self.poly_fit_rsquare (Dict[int, float]): Key - orders of power, 
                Value - goodness of fittings, characterized by R^2.

        :raise Exception: If ``self.ncalc`` <= 1.
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
    Inherited from the Crystal_output class and has thermodynamic-specific 
    attributes. Used for harmonic phonon calclulations. Harmonic object can be
    defined by either a harmonic phonon calculation output file or manually set
    the all the information (usually for QHA).

    Args:
        temperature (array[float] | list[float], optional): Temperatures
            where thermodynamic properties are computed. Unit: K
        pressure (array[float] | list[float], optional): Pressures where
            the thermodyanmic properties are calculated. Unit: GPa
        write_out (bool, optional): Wheter to print out HA thermodynamic 
            properties in a separate text file.
        filename (str, optional): Name of the printed-out file, valid if
            ``write_out`` = True.

    **Note**

    Temperatures and pressures can also be defined by ``self.thermodynamics``,
    whose entries always cover the entries here.

    Usage::

        ha = Harmonic(temperature=[0, 100, 200, 300], pressure=[0.,])
        ha.from_file('harmonic_phonon.out')
    """

    def __init__(self, temperature=[], pressure=[],
                 write_out=True, filename='HA-thermodynamics.dat'):
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
        Generate the Harominc object from a HA output file.

        Args:
            output_name (str): Name of the output file.
            scellphono (array[float] | list[float], optional):
                Supercell manipulation matrix used for vibrational properties. 
                Should be the same as 'SCELPHONO' keyword in CRYSTAL input file.
                By default a 1\*1\*1 'SCELPHONO' is assumed.
            read_eigenvector (bool, optional): Whether to read eigenvectors from
                the output.
            auto_calc (bool, optional): Whether to automatically launch
                thermodynamic calculations. Parameters defined during
                initialization will be used.

        Returns:
            self.strucrure (PyMatGen Structure): Cell reduced by SCELPHONO.
            self.natom (int): Number of atoms in the reduced cell.
            self.volume (float): Volume of the reduced cell. Unit: Angstrom^3
            self.mode (list[Mode]): List of mode objects at all the qpoints.

        :raise Exception: If computational data is stored in the object.
        :raise Exception: If a QHA output file is read.
        """
        import numpy as np
        from CRYSTALpytools.thermodynamics import Mode

        if hasattr(self, "volume"):
            raise Exception(
                "Data exists. Cannot overwrite the existing data.")

        super(Harmonic, self).read_cry_output(output_name)
        super(Harmonic, self).get_mode()
        self._generate_structure(scelphono=scelphono)

        if len(np.unique(self.edft)) != 1:
            raise Exception(
                "Only the frequency calculations at constant volumes are premitted.")
        else:
            self.edft = self.edft[0]

        # Transfer the modes in self.freqency into lists of mode objects
        if read_eigenvector:
            super(Harmonic, self).get_phonon_eigenvector()
        else:
            self.eigenvector = []

        self.from_frequency(self.edft, self.qpoint, self.frequency, self.eigenvector)

        if auto_calc:
            self.thermodynamics(sumphonon=True)
            self.print_results()

        return self

    def from_frequency(self, edft, qpoint, frequency, eigenvector, **kwargs):
        """
        Generate a Harmonic object by specifying frequency and eigenvector and
        clean imaginary frequencies.

        Not recommanded to be used as a standalone method.

        Args:
            edft (float): Electron total energy
            qpoint (list[list[array[float], float]]): Fractional coordinate 
                and weight of qpoint
            frequency (array[float]): Array of frequencies. Unit: THz
            eigenvector (array[float]): Normalized eigenvectors. 
            structure (PyMatGen Structure, optional)
            natom (int, optional)
            volume (float, optional)

        **Note**: The user should define either 'structure' or 'natom' + 'volume'.

        Returns:
            self.nqpoint (int)
            self.nmode (array[int])
            self.mode (list[Mode])
            self.structure (PyMatGen Structure, optional)
            self.natom (int)
            self.volume (float)

        :raise Exception: If computational data is stored in the object.
        :raise Exception: If the 1st dimension (nqpoint) of ``qpoint`` and ``frequency`` are not consistent.
        :raise Exception: If the 2nd dimension (nfreq) of ``frequency`` and ``eigenvector`` are not consistent and ``eigenvector`` is not ``[]``.
        """
        import numpy as np

        if hasattr(self, "mode"):
            raise Exception(
                "Data exists. The current command will be ignored.")

        for key, value in kwargs.items():
            if key == 'structure':
                self.structure = value
                self.natom = len(value.species)
                self.volume = value.lattice.volume
                break
            elif key == 'natom':
                self.natom = int(value)
            elif key == 'volume':
                self.volume = float(value)

        if len(qpoint) != np.size(frequency, 0):
            raise Exception(
                "The 1st dimension (n qpoint) of 'qpoint' and 'frequency' are not consistent.")
        if len(eigenvector) != 0 and np.size(eigenvector, 1) != np.size(frequency, 1):
            raise Exception(
                "The 2nd dimension (n mode) of 'frequency' and 'eigenvector' are not consistent.")

        self.edft = edft
        self.nqpoint = len(qpoint)
        self.qpoint = qpoint
        self.frequency = frequency
        if len(eigenvector) != 0:
            self.eigenvector = eigenvector
        super(Harmonic, self).clean_imaginary()

        # Transfer the modes in self.freqency into lists of mode objects
        self.mode = []
        for q, qpoint_freq in enumerate(frequency):
            qmode = []
            for m, mode_freq in enumerate(qpoint_freq):
                if len(eigenvector) != 0:
                    qmode.append(Mode(rank=m + 1,
                                      frequency=[mode_freq],
                                      volume=[self.volume],
                                      eigenvector=[eigenvector[q, m]])
                                 )
                else:
                    qmode.append(Mode(rank=m + 1,
                                      frequency=[mode_freq],
                                      volume=[self.volume])
                                 )

            self.mode.append(qmode)
        self.nmode = np.array([len(q) for q in self.mode])

        return self

    def _generate_structure(self, scelphono):
        """
        Eliminate the influences of the keyword 'SCELPHONO' and generate the
        PyMatGen Structure object of the actual cell used for phonon 
        calculation. Not a standalone method.

        Args:
            scellphono (list[int] | array[int]): ndimension\*ndimension or 
                3\*3 matrix corresponds to the 'SCELPHONO' keyword.

        Returns:
            self.structure (PyMatGen Structure)
            self.natom (int)
            self.volume (float)
        """
        from CRYSTALpytools.convert import cry_out2pmg
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
        with frequencies = 0 are skipped. Not a standalone method. For 
        thermodynamics, use 'self.thermodyanmics'
        instead.

        Args:
            temperature (float)
            calculate_zp (bool): Calculate zero-point energy or temperature
                dependent properties.

        Returns:
            zp_energy (array[float]): Zero-point energy at a q point. Returned
                if ``calculate_zp = True``.
            U_vib (array[float]): Vibrational contribution to internal energy 
                at constant temperature and a q point. Returned if 
                ``calculate_zp = False``.
            entropy (array[float]): Entropy at constant temperature and a q 
                point. Returned if ``calculate_zp = False``.
            C_v (array[float]): Constant volume specific heat at constant
                 temperature and a q point. Returned if ``calculate_zp = False``.
        """
        import numpy as np

        if calculate_zp:
            zp_energy = []
        else:
            T = temperature
            U_vib = []
            entropy = []
            C_v = []

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
                zp_energy.append(zp_energy_q)
            else:
                U_vib.append(U_vib_q)
                entropy.append(entropy_q)
                C_v.append(C_v_q)

        if calculate_zp:
            zp_energy = np.array(zp_energy, dtype=float)
            return zp_energy
        else:
            U_vib = np.array(U_vib, dtype=float)
            entropy = np.array(entropy, dtype=float)
            C_v = np.array(C_v, dtype=float)
            return U_vib, entropy, C_v

    def thermodynamics(self, sumphonon=True, mutewarning=False, **kwargs):
        """
        Calculate the thermodynamic properties (zp_energy, U_vib, entropy, C_v
        and Helmholtz free energy) of the given system, at all qpoints and the
        whole temperature range.

        The Helmholtz free energy is defined as:

        .. math::

            F(p,V) = E_{DFT} + F_{vib}(T) = E_{DFT} + U_{vib}(T) - TS(T)

        Args:
            temperature (array[float] | list[float], optional): 
                Temperature. Unit: K
            pressure (array[float] | list[float], optional):
                Pressure. Unit: GPa
            sumphonon (bool): Whether to sum up the phonon contributions across
                the sampled q points and take weighted-average.
            mutewarning (bool): Whether print out warning messages of updated
                temperature and pressure (For QHA).

        Returns:
            self.helmholtz (array[float]): Helmholtz free energy 
                (nqpoint\*ntemperature). Unit: KJ/mol
            self.gibbs (array[float]): Gibbs free energy 
                (nqpoint\*ntemperature\*npressure). Unit: KJ/mol
            self.zp_energy (array[float]): Zero-point energy. (nqpoint\*1). 
                Unit: KJ/mol
            self.U_vib (array[float]): Vibrational contribution to internal 
                energy (nqpoint\*ntemperature). Unit: KJ/mol
            self.entropy (array[float]): Entropy (nqpoint\*ntemperature)
                Unit: J/mol\*K
            self.C_v (array[float]): Constant volume specific heat 
                (nqpoint\*ntemperature). Unit: J/mol\*K


        **Note**: If ``sumphonon = True``, nqpoint = 1.

        :raise ValueError: If temperature and pressure are defined neither here nor during initialization
        """
        import warnings
        import numpy as np
        import scipy.constants as scst

        # Generate temperature and pressure series
        if kwargs:
            if 'temperature' in kwargs:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn(
                        'Temperature attribute exists. Input temperatures will be used to update the attribute.')

                self.temperature = np.array(kwargs['temperature'], dtype=float)

            if 'pressure' in kwargs:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn(
                        'Pressure attribute exists. Input pressures will be used to update the attribute.')

                self.pressure = np.array(kwargs['pressure'], dtype=float)
        else:
            if not hasattr(self, 'temperature') or \
               not hasattr(self, 'pressure'):
                raise ValueError(
                    'Temperature and pressure should be specified.')

        zp_energy = self._phonon_sumup(temperature=0., calculate_zp=True)
        U_vib = []
        entropy = []
        C_v = []
        helmholtz = []
        gibbs = []

        for T in self.temperature:
            gibbs_t = []
            U_vib_t, entropy_t, C_v_t = self._phonon_sumup(temperature=T,
                                                           calculate_zp=False)
            helm_t = -entropy_t * T * 1e-3 + U_vib_t + self.edft

            for p in self.pressure:
                gibbs_tp = p * self.volume * scst.Avogadro * 1e-24 + helm_t
                gibbs_t.append(gibbs_tp)

            # nTemp * nqpoint
            U_vib.append(U_vib_t)
            entropy.append(entropy_t)
            C_v.append(C_v_t)
            helmholtz.append(helm_t)
            # nTemp * npress * nqpoint
            gibbs.append(gibbs_t)

        if sumphonon:
            wt = np.array([qp[1] for qp in self.qpoint])
            self.nqpoint = 1
            self.qpoint = [[np.array([0., 0., 0.]), 1.]]
            self.zp_energy = np.array([np.dot(zp_energy, wt)])
            self.U_vib = np.array([np.dot(U_vib, wt)])
            self.entropy = np.array([np.dot(entropy, wt)])
            self.C_v = np.array([np.dot(C_v, wt)])
            self.helmholtz = np.array([np.dot(helmholtz, wt)])
            self.gibbs = np.array([np.dot(gibbs, wt)])
        else:
            self.zp_energy = zp_energy
            self.U_vib = np.transpose(np.array(U_vib, dtype=float))
            self.entropy = np.transpose(np.array(entropy, dtype=float))
            self.C_v = np.transpose(np.array(C_v, dtype=float))
            self.helmholtz = np.transpose(np.array(helmholtz, dtype=float))
            self.gibbs = np.transpose(np.array(gibbs, dtype=float), (2, 0, 1))

        return self.helmholtz, self.gibbs, self.zp_energy, self.U_vib,\
            self.entropy, self.C_v

    def print_results(self):
        """
        Print HA thermodynamic results into an external file. Used if 
        ``write_out = True``.

        **Note**

        Phonon dispersions are forced to be summed to keep the output
        succinct if the automatic scheme (write_out=True) is launched. To get
        verbose outputs, directly call this method.
        """
        import scipy.constants as scst
        from CRYSTALpytools.units import eV_to_H

        if not self.write_out:
            print('Harmonic.write_out = False, return to empty.')
            return

        file = open(self.filename, 'w')
        file.write('%-21s%12.4e%-15s%12.4e%-10s\n' %
                   ('# DFT TOTAL ENERGY = ', eV_to_H(self.edft),
                    ' eV,         = ', self.edft, ' kJ/mol'))
        file.write('%-21s%12.6f%-15s%12.6f%-10s\n' %
                   ('# CELL VOLUME      = ', self.volume,
                    ' Angstrom^3, = ', self.volume * scst.Avogadro * 1e-24, ' cm^3/mol'))
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
    Generate and rearrange harmonic phonons, store the fitted, volume-dependent
    QHA phonon information and obtain the QHA thermodynamic properties.

    Args:
        temperature (array[float] | list[float], optional): Unit: K
        pressure (array[float] | list[float], optional): Unit: GPa
        write_out (bool): Whether to print the key information into a file.
        filename (str): Name of the output file. Valid if 
                ``write_out = True``.

    **Note**

    Temperatures and pressures can also be defined by ``self.thermodynamics``,
    whose entries always cover the entries here.

    Usage::

        qha = Quasi_harmonic()
        qha.from_QHA_file('qha_phonon.out')
        qha.thermo_freq(eos_method='birch_murnaghan', temperature=[0, 100, 200, 300], pressure=[0., 0.5, 1.]):
    """

    def __init__(self, temperature=[], pressure=[],
                 write_out=True, filename='QHA-Fit.dat'):
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

        Args:
            input_files (list[str]): List of phonon output filenames.
            scelphono (array[float] | list[float]], optional): Corresponds
                to the 'SCELPHONO' keyword in CRYSTAL. Either 3\*3 or 
                ndimension\*ndimension. By default a 1\*1\*1 'SCELPHONO' is 
                assumed.
            overlap (float, optional): The threshold of close mode overlaps.
            sort_phonon (bool, optional): Whether to check phonon continuity.

        Returns:
            self.ncalc (int): Number of HA phonon calculations.
            self.combined_phonon (list[Harmonic]): List of Harmonic objects.
            self.combined_volume (list[float]): Volumes. Unit: Angstrom^3
            self.combined_edft (list[float]): DFT total energies. Unit: KJ/mol
            self.combined_mode (list[Mode]): List of mode objects.
        """
        from CRYSTALpytools.thermodynamics import Harmonic
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
        self.nqpoint = ha_list[0].nqpoint
        self.qpoint = ha_list[0].qpoint # consistency of nqpoint is checked, but not qpoint.

        return self

    def from_QHA_file(self, input_file, scelphono=[], overlap=0.4, sort_phonon=True):
        """
        Read data from a single QHA calculation at Gamma point.

        Args:
            input_files (str | list[str]): Only 1 QHA file is permitted.
            scelphono (array[float] | list[float], optional)
            overlap (float, optional)
            sort_phonon (bool, optional)

        Returns:
            self.ncalc (int)
            self.combined_phonon (list[Harmonic])
            self.combined_volume (list[float])
            self.combined_edft (list[float])
            self.combined_mode (list[Mode])

        :raise Exception: If multiple files are defined.
        """
        from CRYSTALpytools.crystal_io import Crystal_output
        from CRYSTALpytools.thermodynamics import Harmonic
        from CRYSTALpytools.thermodynamics import Mode
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
        file.get_phonon_eigenvector()
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
            ha.from_frequency(file.edft[idx_c], np.array([[0, 0, 0]]),
                              np.array([file.frequency[idx_c]]),
                              np.array([file.eigenvector[idx_c]]),
                              structure=structures[idx_c + 1])  # The first one is pre-opt geom
            ha_list.append(ha)

        self.combined_phonon, self.combined_volume, self.combined_edft, \
            self.combined_mode = self._combine_data(ha_list, overlap=overlap,
                                                    sort_phonon=sort_phonon)
        self.nqpoint = 1
        self.qpoint = [[np.array([0., 0., 0.]), 1.]]

        return self

    def _combine_data(self, ha_list, overlap, sort_phonon):
        """
        Combine the HA calculation data and rearrange it in the ascending order 
        of volumes. Not a standalone method.

        Args:
            ha_list (list[Harmonic]): List of harmonic objects.
            overlap (float)
            sort_phonon (bool)

        Returns:
            combined_phonon (list[Harmonic])
            combined_volume (list[float])
            combined_edft (list[float])
            combined_mode (list[Mode])

        :raise Exception: If number of q points, modes or atoms are not consistent across the HA calculations.
        """
        import numpy as np
        import warnings
        from CRYSTALpytools.thermodynamics import Mode

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

    @staticmethod
    def _phonon_continuity(freq, eigvt, symm=None, overlap=0.4):
        """
        Rearrange phonon modes by their continuity. If the difference between
        the maximum scalar product of corresponding eigenvectors (normalized to 
        1) and scalar products of other modes is less than 0.4, warning is
        printed due to the potential overlap of modes. Adopted from CRYSTAL17.

        .. note::

            A. Erba, *J. Chem. Phys.*, 2014, **141**, 124115.

        Not a standalone method.

        Args:
            freq (array[float]): Phonon frequencies. Unit: THz
            eigvt (array[float]): Eigenvectores normalized to 1
            symm (array[float]): Sub-group numbers of corresponding modes.
                *Not implemented*
            overlap (float): The threshold of close mode overlaps.

        Returns:
            freq (array[float]): Sorted phonon frequencies
            eigvt (array[float]): Sorted eigenvectores
            close_overlap (array[bool]):ncalc\*nmode\*nmode. Whether close 
                overlap is identified between the previous calculation (2nd 
                dimension) and the current one (3rd).
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

    def edft_eos_fit(self, method, **kwargs):
        """
        Fit electron total energy according to equation of states. 

        Args:
            method (str): Name of EoS used. Consistent with
                `PyMatGen <https://pymatgen.org/pymatgen.analysis.eos.html>`_.
            order (int): For the DeltaFactor method. *Not implemented*
            min_ndata_factor, max_poly_order_factor, min_poly_order_factor (int):
                For the NumericalEOS method. *Not implemented*

        Returns:
            self.eos_method (string): Name of the fitted equation of state
            self.eos (PyMatGen EOS): The fitted equation of state.
        """
        import re
        from pymatgen.analysis.eos import EOS
        import scipy.constants as scst

        self.eos_method = method
# Commented due to the inhierant problem of pymatgen EOS object
#         if re.findall(r'deltafactor', method, re.I):
#             self.eos = EOS(method).fit(self.combined_volume, self.combined_edft, order=fitargs['order'])

#         elif re.findall(r'numerical_eos', method, re.I):
#             self.eos = EOS(method).fit(self.combined_volume, self.combined_edft,
#                                        min_ndata_factor=fitargs['min_ndata_factor'],
#                                        max_poly_order_factor=fitargs['max_poly_order_factor'],
#                                        min_poly_order_factor=fitargs['min_poly_order_factor'])

#         else:
#             self.eos = EOS(method).fit(self.combined_volume, self.combined_edft)
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
                                                       self.eos.b0 * 10 / scst.Avogadro,
                                                       self.eos.b1))
            file.write('\n')
            file.close()

        return self.eos

    def freq_polynomial_fit(self, order):
        """
        Fit phonon frequencies as polynomial functions of volumes. 

        Args:
            order (list[int] | array[int]): The order of polynomials used.

        Returns:
            self.fit_order (int): The optimal order of polynomial fit.

        Please also refer to ``self.poly_fit`` and ``self.poly_fit_rsquare`` 
        attributes of Mode class.
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

        Args:
            volume (float): Unit: Angstrom^3

        Returns:
            ha (Harmonic): Harmonic phonon object with numerical data.

        :raise Exception: If frequency is not fitted as function of volumes.
        """
        import numpy as np
        from CRYSTALpytools.thermodynamics import Harmonic
        from CRYSTALpytools.thermodynamics import Mode

        if not hasattr(self, 'fit_order') or not hasattr(self, 'eos'):
            raise Exception('ERROR: Analytical expressions unavailable.')

        eq_point = np.argmin(self.combined_edft)
        num_mode = []
        for mode_q in self.combined_mode:
            num_mode_q = []
            for idx_m, mode in enumerate(mode_q):
                d_vol = volume - self.combined_volume[eq_point]
                freq = mode.poly_fit[self.fit_order](d_vol) + mode_q[idx_m].frequency[eq_point]
                num_mode_q.append(freq)
            num_mode.append(num_mode_q)

        num_mode = np.array(num_mode)
        ha = Harmonic(write_out=False).from_frequency(
            self.eos(volume), self.qpoint, num_mode, [], volume=volume
        )

        return ha

    def _minimize_gibbs(self, volume, temperature, pressure):
        """
        Get Gibbs free energy from the Harmonic phonon object. Used only for
        minimizing G(V; T, p) by SciPy. Not a standalone method.

        Args:
            volume (float)
            temperature (float)
            pressure (float)

        Returns:
            ha.gibbs (float): Gibbs free energy. Unit: KJ/mol
        """
        ha = self._get_harmonic_phonon(volume)
        ha.thermodynamics(temperature=[temperature], pressure=[pressure])

        return ha.gibbs[0, 0, 0]

    def thermo_freq(self, eos_method='birch_murnaghan', poly_order=[2, 3],
                    min_method='BFGS', volume_bound=None, mutewarning=False,
                    **kwargs):
        """
        Obtain thermodynamic properties by explicitly fitting phonon 
        frequencies as polynomial functions of volume. DFT total energies are
        fitted as a function of volume by equation of states (EOS).

        Args:
            eos_method (str, optional): EOS used to fit DFT total energies. 
            poly_order (array[int] | list[int], optional): The order of 
                polynomials used to fit frequency as the function of volumes.
            min_method (string, optional): Minimisation algorithms. 
            volume_bound (tuple-like, optional), Boundary conditions of 
                equilibrium volumes. Unit: Angstrom^3
            mutewarning (bool, optional): Whether print out warning messages.
            temperature (array[float], optional): Unit: K
            pressure (array[float], optional): Unit: GPa
            order (int, optional): For DeltaFactor EOS. *Not implemented*
            min_ndata_factor, max_poly_order_factor, min_poly_order_factor (int, optional):
                For Numerical EOS. *Not implemented*

        **Notes**

        1. EOS supported by ``eos_method`` are consistent with `PyMatGen <https://pymatgen.org/pymatgen.analysis.eos.html>`_.
        2. Parameterized and tested algorithms for ``min_method``: 
            * BFGS(no boundary)
            * L-BFGS-B(with boundary)

        Returns:
            self (Quasi_harmonic)

        **New Attributes**

        * ``self.temperature`` in K and ``self.pressure`` in GPa.
        * ``self.equilibrium_volume``, nPressure\*nTemperature. Equilibrium volumes. Unit: Angstrom^3
        * ``self.helmholtz`` and ``self.gibbs``, nPressure\*nTemperature. Helmholtz and Gibbs free energies. Unit: kJ/mol
        * ``self.entropy``, nPressure\*nTemperature, Entropy. Unit: J/mol\*K

        :raise ValueError: If temperature or pressure is defined neither here nor during initialization.
        """
        import numpy as np
        import warnings
        import re
        from scipy.optimize import minimize

        # Generate temperature and pressure series
        if kwargs:
            if 'temperature' in kwargs:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn(
                        'Temperature attribute exists. Input temperatures will be used to update the attribute.', stacklevel=2)

                self.temperature = np.array(kwargs['temperature'], dtype=float)

            if 'pressure' in kwargs:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn(
                        'Pressure attribute exists. Input pressures will be used to update the attribute.', stacklevel=2)

                self.pressure = np.array(kwargs['pressure'], dtype=float)

        if not hasattr(self, 'temperature') or not hasattr(self, 'pressure'):
            raise ValueError('Temperature and pressure should be specified.')

        # Fit DFT total energy, if not done yet. Otherwise, fitted values will not be covered.
        if hasattr(self, 'eos') and not mutewarning:
            warnings.warn(
                'DFT total energy is already fitted. To keep the consistency, it will not be updated.', stacklevel=2)
        else:
            # Commented due to the inhierant problem of pymatgen EOS object
            #             if re.findall(r'deltafactor', eos_method, re.I):
            #                 if 'order' not in kwargs.keys():
            #                     kwargs.update({'order' : 3})

            #                 self.edft_eos_fit(method=eos_method, order=kwargs['order'])

            #             elif re.findall(r'numerical_eos', eos_method, re.I):
            #                 if 'min_ndata_factor' not in kwargs.keys():
            #                     kwargs.update({'min_ndata_factor' : 3})

            #                 if 'max_poly_order_factor' not in kwargs.keys():
            #                     kwargs.update({'max_poly_order_factor' : 5})

            #                 if 'min_poly_order_factor' not in kwargs.keys():
            #                     kwargs.update({'min_poly_order_factor' : 2})

            #                 self.edft_eos_fit(method=eos_method,
            #                                   min_ndata_factor=kwargs['min_ndata_factor'],
            #                                   max_poly_order_factor=kwargs['max_poly_order_factor'],
            #                                   min_poly_order_factor=kwargs['min_poly_order_factor'])

            #             else:
            #                 self.edft_eos_fit(method=eos_method)
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

    def thermo_eos(self, eos_method='birch_murnaghan', poly_order=[2, 3],
                   mutewarning=False, **kwargs):
        """
        Obtain thermodynamic properties by fitting EOS, which is fitted by the 
        Helmholtz free energies of sampled harmonic phonons. The explicit 
        sorting and fitting of frequency-volume relationship is disabled.

        Args:

            eos_method (str, optional)
            poly_order (array[int] | list[int], optional): Order of 
                polynomials used to fit Gibbs free energy to get entropy.
            mutewarning (bool, optional)
            temperature (array[float] | list[float], optional): Unit: K
            pressure (array[float] | list[float], optional): Unit: GPa
            order (int, optional): For DeltaFactor EOS. *Not implemented*
            min_ndata_factor, max_poly_order_factor, min_poly_order_factor (int, optional):
                For Numerical EOS. *Not implemented*

        Returns:
            self (Quasi_harmonic)

        New Attributes are consistent with the ``thermo_freq`` method

        :raise Exception: If the number of HA calculations is less than 4.
        :raise ValueError: If temperature or pressure is defined neither here nor during initialization.
        """
        import numpy as np
        import warnings
        import re
        from CRYSTALpytools.thermodynamics import Harmonic
        from pymatgen.analysis.eos import EOS
        from scipy.optimize import fmin
        from sympy import diff, lambdify, symbols

        # Check the number of calculations
        if self.ncalc < 4:
            raise Exception('Insufficient database. Increase HA phonons')

        # Generate temperature and pressure series
        if kwargs:
            if 'temperature' in kwargs:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn(
                        'Temperature attribute exists. Input temperatures will be used to update the attribute.',
                        stacklevel=2)

                self.temperature = np.array(
                    kwargs['temperature'], dtype=float)

            if 'pressure' in kwargs:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn(
                        'Pressure attribute exists. Input pressures will be used to update the attribute.',
                        stacklevel=2)

                self.pressure = np.array(kwargs['pressure'], dtype=float)

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
            # Commented due to the inhierant problem of pymatgen EOS object
            #             if re.findall(r'deltafactor', eos_method, re.I):
            #                 if 'order' not in kwargs.keys():
            #                     kwargs.update({'order' : 3})

            #                 hfe_t = EOS(eos_method).fit(self.combined_volume, helmholtz[idx_t], order=kwargs['order'])

            #             elif re.findall(r'numerical_eos', eos_method, re.I):
            #                 if 'min_ndata_factor' not in kwargs.keys():
            #                     kwargs.update({'min_ndata_factor' : 3})

            #                 if 'max_poly_order_factor' not in kwargs.keys():
            #                     kwargs.update({'max_poly_order_factor' : 5})

            #                 if 'min_poly_order_factor' not in kwargs.keys():
            #                     kwargs.update({'min_poly_order_factor' : 2})

            #                 hfe_t = EOS(eos_method).fit(self.combined_volume, helmholtz[idx_t],
            #                                             min_ndata_factor=kwargs['min_ndata_factor'],
            #                                             max_poly_order_factor=kwargs['max_poly_order_factor'],
            #                                             min_poly_order_factor=kwargs['min_poly_order_factor'])

            #             else:
            #                 hfe_t = EOS(eos_method).fit(self.combined_volume, helmholtz[idx_t])
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
