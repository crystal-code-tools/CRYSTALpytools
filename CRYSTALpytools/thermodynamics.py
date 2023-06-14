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
    Defines a vibrational mode and do analysis per mode.

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
        Get the zero-point energy of a single mode. *ncalc = 1 only*.

        .. math::

            E^{zp}_{i,\\mathbf{q}}=\\frac{1}{2}\\hbar\\omega_{i,\\mathbf{q}}

        Returns:
            self.zp_energy (float): Zero-point energy. Unit: KJ/mol
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise AttributeError(
                'This module is limited to a single frequency calculation.')

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e9
        self.zp_energy = 0.5 * hbar_freq

        return self.zp_energy

    def get_u_vib(self, temperature=298.15):
        """
        Get the vibration contribution to internal energy (including zero-point
        energy) of a single mode. *ncalc = 1 only*.

        .. math::

            U^{vib}_{i,\\mathbf{q}}\\left(T\\right)=E^{zp}_{i,\\mathbf{q}}+
            \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{\\exp{\\left(
                \\frac{\\hbar\\omega_{i,\\mathbf{q}}}{k_{B}T}
            \\right)}-1}

        Args:
            temperature (float, optional): Temperature where the quantity is
                computed. Unit: K

        Returns:
            self.u_vib (float): Vibration contribution to internal energy.
                Unit: KJ/mol
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise AttributeError(
                'This module is limited to a single frequency calculation.')

        if not hasattr(self, 'zp_energy'):
            self.get_zp_energy()

        if temperature == 0 or self.frequency[0] < 1e-4:
            self.u_vib = self.zp_energy
            return self.u_vib

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e9
        kb_t = scst.k * scst.Avogadro * temperature * 1e-3
        expon = np.exp(hbar_freq / kb_t)
        self.u_vib = self.zp_energy + hbar_freq / (expon - 1)

        return self.u_vib

    def get_entropy(self, temperature):
        """
        Get the entropy of a single mode. *ncalc = 1 only*.

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
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise AttributeError(
                'This module is limited to a single frequency calculation.')

        if temperature == 0 or self.frequency[0] < 1e-4:
            self.entropy = 0
            return self.entropy

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e12
        kb_t = scst.k * scst.Avogadro * temperature
        expon = np.exp(hbar_freq / kb_t)
        entS = kb_t * (hbar_freq / kb_t / (expon - 1) - np.log(1 - 1 / expon))
        self.entropy = entS / temperature

        return self.entropy

    def get_c_v(self, temperature):
        """
        Get the constant volume specific heat of a single mode. *ncalc = 1 only*.

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
            self.c_v (float): Constant volume specific heat. Unit: J/mol\*K
        """
        import numpy as np
        import scipy.constants as scst

        if self.ncalc > 1:
            raise AttributeError(
                'This module is limited to a single frequency calculation.')

        if temperature == 0 or self.frequency[0] < 1e-4:
            self.c_v = 0
            return self.c_v

        hbar_freq = self.frequency[0] * scst.Avogadro * scst.h * 1e12
        kb_t = scst.k * scst.Avogadro * temperature
        expon = np.exp(hbar_freq / kb_t)

        self.c_v = hbar_freq**2 / kb_t / temperature * expon / (expon - 1)**2

        return self.c_v

    def polynomial_fit(self, order=[2, 3]):
        """
        Fit phonon frequency as the polynomial function of volume. *ncalc > 1 only*.

        .. note::

            To improve the accuracy of fittings, :math:`\\Delta\\omega(\\Delta V)`
            is fitted as a polynomial function without the constant term.

            :math:`\\Delta V=V-V_{min}` is used so HA phonons of the most compact
            structure is kept. See 'FIXINDEX' keyword in CRYSTAL manual for
            further information.

        Args:
            order (array[int] | list[int]], optional): Orders of polynomials.

        Returns:
            self.poly_fit (Dict[int, NumPy Polynomial]): Key - orders of power,
                Value - fitted NumPy polynomials
            self.poly_fit_rsquare (Dict[int, float]): Key - orders of power,
                Value - goodness of fittings, characterized by R^2.
        """
        import numpy as np
        import warnings
        from scipy.optimize import least_squares
        from CRYSTALpytools.thermodynamics import Quasi_harmonic

        if self.ncalc <= 1:
            raise AttributeError(
                'This modulus is limited to multiple frequency calculations.')

        if max(order) > self.ncalc - 1:
            warnings.warn(
                'Reference data not sufficient for the order of polynomial fitting.')
            warnings.warn('Too high values will be removed.')

        order = list(set(order))
        order = [p for p in order if p <= self.ncalc - 1]

        self.poly_fit = {}
        self.poly_fit_rsqaure = {}

        if self.rank == 1 or self.rank == 2 or self.rank == 3:
            for i in order:
                self.poly_fit[i] = np.polynomial.polynomial.Polynomial(
                    [0. for i in range(i + 1)])
                self.poly_fit_rsqaure[i] = 1.

            return order, self.poly_fit, self.poly_fit_rsqaure

        qha = Quasi_harmonic(write_out=False)
        idx_vmin = np.argmin(self.volume)
        vmin = self.volume[idx_vmin]
        fmin = self.frequency[idx_vmin]
        dv = self.volume - vmin
        df = self.frequency - fmin
        for i in order:
            opt = least_squares(qha._poly_no_cst,
                                np.array([1. for j in range(i)]),
                                args=(dv, df))
            poly = np.polynomial.polynomial.Polynomial(np.insert(opt.x, 0, 0.))
            self.poly_fit[i] = poly
            self.poly_fit_rsqaure[i] = 1 - np.sum((df - poly(dv))**2) / np.sum((df - np.mean(df))**2)

        return order, self.poly_fit, self.poly_fit_rsqaure

    def get_gruneisen(self, order, volume):
        """
        Return to mode Grüneisen parameter. *ncalc > 1 only*.

        .. math::

            \\gamma = -\\frac{V}{\\omega(V)}\\frac{\\partial\\omega}{\\partial V}

        Args:
            order (int | list[int]): See ``polynomial_fit``
            volume (float | array): Typically the equilibrium volume

        Returns:
            self.gruneisen (dict): Key, order; Value, Gruneisen parameter
        """
        import numpy as np

        if not hasattr(self, 'poly_fit'):
            raise AttributeError(
                'Polynomial fitting is required to get Gruneisen parameters.')

        order = list(set(order))
        self.gruneisen = {}
        if type(volume) == float:
            volume = np.array([volume])
        else:
            volume = np.array(volume)

        if self.rank == 1 or self.rank == 2 or self.rank == 3:
            for i in order:
                self.gruneisen[i] = np.zeros(volume.shape)
            return self.gruneisen

        idx_vmin = np.argmin(self.volume)
        vmin = self.volume[idx_vmin]
        fmin = self.frequency[idx_vmin]
        dv = volume - vmin
        for i in order:
            dpoly = self.poly_fit[i].deriv(1)
            self.gruneisen[i] = -dpoly(dv) * volume / (self.poly_fit[i](dv) + fmin)

        return self.gruneisen


class Harmonic():
    """
    A class for harmonic phonon calclulations. It can be parameterized from a
    CRYSTAL output file, phonopy ouput file or by setting all the information
    (usually for QHA).

    Args:
        temperature (array[float] | list[float], optional): Temperatures
            where thermodynamic properties are computed. Unit: K
        pressure (array[float] | list[float], optional): Pressures where
            the thermodyanmic properties are calculated. Unit: GPa
        write_out (bool, optional): Wheter to print out HA thermodynamic
            properties in a separate text file.
        filename (str, optional): Name of the printed-out file, valid if
            ``write_out`` = True.

    Temperatures and pressures can also be defined by ``self.thermodynamics``,
    whose entries always cover the entries here.

    Usage::

        ha = Harmonic(temperature=[0, 100, 200, 300], pressure=[0.,])
        ha.from_file('harmonic_phonon.out')
    """

    def __init__(self, temperature=[], pressure=[], autocalc=False,
                 write_out=True, filename='HA-thermodynamics.dat'):
        import numpy as np

        if len(temperature) > 0:
            self.temperature = np.array(temperature, dtype=float)

        if len(pressure) > 0:
            self.pressure = np.array(pressure, dtype=float)

        self.autocalc = autocalc
        self.write_out = write_out
        if self.write_out:
            self.filename = filename
        else:
            self.filename = 'no file'

    def from_file(self, output_name, scelphono=[], read_eigvt=False,
                  imaginary_tol=-1e-4, q_overlap_tol=1e-4):
        """
        Generate the Harominc object from a HA output file. Imaginary modes and
        overlapped q points are forced to be cleaned.

        Args:
            output_name (str): Name of the output file.
            scellphono (array[float] | list[float], optional):
                The 'SCELPHONO' keyword in CRYSTAL input file. By default a
                1\*1\*1 'SCELPHONO' is assumed.
            read_eigvt (bool): Whether to read eigenvectors from output.
            imaginary_tol (float): The threshold of negative frequencies.
            q_overlap_tol (float): The threshold of overlapping points, defined
                as the 2nd norm of the difference of fractional q vectors

        Returns:
            self.structure (PyMatGen Structure): Cell reduced by SCELPHONO.
            self.natom (int): Number of atoms in the reduced cell.
            self.volume (float): Volume of the reduced cell. Unit: Angstrom^3
            self.edft (float)
            self.nqpoint (int)
            self.qpoint (list)
            self.nmode (array[int])
            self.mode (list[Mode]): List of mode objects at all the qpoints.

        :raise ValueError: If a QHA output file is read.
        """
        import numpy as np
        from CRYSTALpytools.crystal_io import Crystal_output
        from CRYSTALpytools.thermodynamics import Mode
        import warnings

        if hasattr(self, "volume"):
            warnings.warn("Data exists. Cannot overwrite the existing data.")
            return self

        output = Crystal_output().read_cry_output(output_name)
        output.get_phonon(read_eigvt=read_eigvt, rm_imaginary=False, rm_overlap=False)
        strucs = _restore_pcel(output, scelphono)

        if len(strucs) != 1: # strucs and edft must have only 1 valid entry
            raise ValueError("Only the frequency calculations at constant volumes are premitted.")

        # Transfer the modes in self.freqency into lists of mode objects
        self.from_frequency(output.edft[0], output.qpoint, output.frequency,
                            output.eigenvector, structure=strucs[0],
                            imaginary_tol=imaginary_tol, q_overlap_tol=q_overlap_tol)

        # Autocalc
        if self.autocalc == True:
            self.thermodynamics(sumphonon=True)

        return self

    def from_phonopy(self, phono_yaml, struc_yaml=None, edft=None,
                     imaginary_tol=-1e-4, q_overlap_tol=1e-4, q_id=None, q_coord=None):
        """
        Build a Harmonic object from `Phonopy <https://phonopy.github.io/phonopy/>`_
        'band.yaml' or 'qpoints.yaml' file.

        Args:
            phono_yaml (str): Phonopy band.yaml or qpoint.yaml file
            struc_yaml (str): Phonopy phonopy.yaml or phonopy_disp.yaml file.
                *Needed only if a qpoint.yaml file is read.*
            edft (float): DFT energy
            imaginary_tol (float): The threshold of negative frequencies.
            q_overlap_tol (float): The threshold of overlapping points, defined
                as the 2nd norm of the difference of fractional q vectors
            q_id (list[int]): Specify the id (from 0) of q points to be read.
                nqpoint\*1 list.
            q_coord (list[list]): Specify the coordinates of q points to be
                read. nqpoint\*3 list.

        ``q_id`` and ``q_coord`` should not be set simultaneously. If set,
        ``q_id`` takes priority and ``q_coord`` is ignored. If both are none,
        all the points will be read.

        :raise Exception: If the length unit in yaml file is neither 'au' nor 'angstrom'.
        :raise Exception: If q point is not found.
        """
        import yaml
        import numpy as np
        from CRYSTALpytools.units import au_to_angstrom
        from pymatgen.core.structure import Structure
        import warnings

        if hasattr(self, "volume"):
            warnings.warn("Data exists. Cannot overwrite the existing data.")
            return self
        if edft == None:
            edft = 0.
            warnings.warn('DFT energy is set to 0.')

        file = open(phono_yaml, 'r', errors='ignore')
        data = yaml.safe_load(file)
        file.close()
        if struc_yaml != None:
            struc_file = open(struc_yaml, 'r')
            struc_data = yaml.safe_load(struc_file)
            struc_file.close()
        else:
            struc_data = data

        # Get unit
        try: # band.yaml
            len_unit = struc_data['length_unit']
        except KeyError: # phonopy.yaml
            len_unit = struc_data['physical_unit']['length']

        if len_unit == 'angstrom':
            unit_len = 1.0
        elif len_unit == 'au':
            unit_len = au_to_angstrom(1.0)
        else:
            raise Exception("Unknown length unit. Available options: au, angstrom.")

        # Get structure
        spec = []
        coord = []
        try: # band.yaml
            latt = np.array(struc_data['lattice'], dtype=float) * unit_len
            for idx_a, atom in enumerate(struc_data['points']):
                spec.append(atom['symbol'])
                coord.append(atom['coordinates'])
        except KeyError: # phonopy.yaml
            latt = np.array(struc_data['primitive_cell']['lattice'], dtype=float) * unit_len
            for idx_a, atom in enumerate(struc_data['primitive_cell']['points']):
                spec.append(atom['symbol'])
                coord.append(atom['coordinates'])

        structure = Structure(lattice=latt, species=spec, coords=coord)
        natom = len(spec)

        if q_id == None and q_coord == None:
            nqpoint = data['nqpoint']
            qinfo = np.array(range(nqpoint), dtype=int)
        elif q_id != None:
            qinfo = np.array(q_id, dtype=int)
            nqpoint = len(qinfo)
        elif q_id == None and q_coord != None:
            qinfo = np.array(q_coord, dtype=float)
            nqpoint = len(qinfo)

        qpoint = [[np.zeros([3, 1]), 1 / nqpoint] for i in range(nqpoint)]
        nmode = np.array([3 * natom for i in range(nqpoint)]) # No fragment phonon is assumed.
        frequency = np.zeros([nqpoint, 3 * natom])
        # Read phonon
        real_q = 0
        for idx_p, phonon in enumerate(data['phonon']):
            if real_q == nqpoint:
                break

            if len(qinfo.shape) == 1: # q_id and all q points
                if idx_p == qinfo[real_q]:
                    qpoint[real_q][0] = np.array(phonon['q-position'])
                    frequency[real_q, :] = np.array([i['frequency'] for i in phonon['band']])
                    real_q += 1
                else:
                    continue
            else: # q_coord
                coord = np.array(phonon['q-position'])
                if np.linalg.norm(qinfo[real_q] - coord) < 1e-4:
                    qpoint[real_q][0] = coord
                    frequency[real_q, :] = np.array([i['frequency'] for i in phonon['band']])
                    real_q += 1
                else:
                    continue

        if real_q < nqpoint:
            raise Exception('Some q points are missing from the yaml file.')

        # set object
        self.from_frequency(edft=edft, qpoint=qpoint, frequency=frequency,
                            eigenvector=[], structure=structure,
                            imaginary_tol=imaginary_tol, q_overlap_tol=q_overlap_tol)

        return self

    def from_frequency(self, edft, qpoint, frequency, eigenvector,
                       structure=None, natom=None, volume=None,
                       imaginary_tol=-1e-4, q_overlap_tol=1e-4):
        """
        Generate a Harmonic object by specifying frequency and eigenvector.
        Imaginary modes and overlapped q points are forced to be cleaned.

        Args:
            edft (float): Electron total energy
            qpoint (list[list[array[float], float]]): Fractional coordinate
                and weight of qpoint
            frequency (array[float]): Array of frequencies. Unit: THz
            eigenvector (array[float]): Normalized eigenvectors.
            structure (Pymatgen Structure)
            natom (int)
            volume (float)
            imaginary_tol (float): The threshold of negative frequencies.
            q_overlap_tol (float): The threshold of overlapping points, defined
                as the 2nd norm of the difference of fractional q vectors

        .. note::

            The user should define either ``structure`` or ``natom`` + ``volume``.

        Returns:
            self.structure (PyMatGen Structure): Cell reduced by SCELPHONO.
            self.natom (int): Number of atoms in the reduced cell.
            self.volume (float): Volume of the reduced cell. Unit: Angstrom^3
            self.edft (float)
            self.nqpoint (int)
            self.qpoint (list)
            self.nmode (array[int])
            self.mode (list[Mode]): List of mode objects at all the qpoints.

        :raise AttributeError: If computational data is stored in the object.
        :raise ValueError: If neither of the 2 available options are defined.
        """
        from CRYSTALpytools.base.crysout import PhononBASE
        from CRYSTALpytools.thermodynamics import Mode
        import numpy as np

        if hasattr(self, "mode"):
            raise AttributeError("Data exists. The current command will be ignored.")

        if structure != None:
            self.structure = structure
            self.natom = len(structure.species)
            self.volume = structure.lattice.volume
        elif natom != None and volume != None:
            self.natom = int(natom)
            self.volume = float(volume)
        else:
            raise ValueError('Geometry is not sufficiently defined. Structure or volume + natom are needed.')

        if len(qpoint) != np.size(frequency, 0):
            raise ValueError("The 1st dimension (n qpoint) of 'qpoint' and 'frequency' are not consistent.")
        if len(eigenvector) != 0 and np.size(eigenvector, 1) != np.size(frequency, 1):
            raise ValueError("The 2nd dimension (n mode) of 'frequency' and 'eigenvector' are not consistent.")

        self.edft = edft
        self.nqpoint = len(qpoint)
        self.qpoint = qpoint
        self.nmode = np.array([len(q) for q in frequency])
        self.frequency = frequency
        self.intens = []
        self.IR = []
        self.Raman = []
        self.eigenvector = eigenvector
        ## Note: Harmonic object is not a crystal_output project, but has the
        ## same attributes
        self = PhononBASE.clean_imaginary(self, threshold=imaginary_tol)
        self = PhononBASE.clean_q_overlap(self, threshold=q_overlap_tol)

        # Transfer the modes in self.freqency into lists of mode objects
        self.mode = []
        for q, freq_q in enumerate(self.frequency):
            qmode = []
            for m, freq_m in enumerate(freq_q):
                if len(self.eigenvector) != 0:
                    qmode.append(Mode(rank=m + 1,
                                      frequency=[freq_m],
                                      volume=[self.volume],
                                      eigenvector=[self.eigenvector[q, m]]))
                else:
                    qmode.append(Mode(rank=m + 1,
                                      frequency=[freq_m],
                                      volume=[self.volume]))

            self.mode.append(qmode)

        # Delete useless attribute
        delattr(self, 'intens')
        delattr(self, 'IR')
        delattr(self, 'Raman')

        return self

    def _phonon_sumup(self, temperature, calculate_zp):
        """
        Summing up inidival phonon modes at each q point. Translational modes
        with frequencies = 0 are skipped. For thermodynamics, directly call
        ``self.thermodyanmics()``.

        Args:
            temperature (float)
            calculate_zp (bool): Calculate zero-point energy or temperature
                dependent properties.

        Returns:
            zp_energy (array[float]): Zero-point energy at a q point. Returned
                if ``calculate_zp = True``.
            u_vib (array[float]): Vibrational contribution to internal energy 
                at constant temperature and a q point. Returned if 
                ``calculate_zp = False``.
            entropy (array[float]): Entropy at constant temperature and a q 
                point. Returned if ``calculate_zp = False``.
            c_v (array[float]): Constant volume specific heat at constant
                 temperature and a q point. Returned if ``calculate_zp = False``.
        """
        import numpy as np

        if calculate_zp:
            zp_energy = []
        else:
            T = temperature
            u_vib = []
            entropy = []
            c_v = []

        for qpoint in self.mode:
            if calculate_zp:
                zp_energy_q = 0.
            else:
                u_vib_q = 0.
                entropy_q = 0.
                c_v_q = 0.
            # Remove the translational modes
            for mode in qpoint:
                if np.isnan(mode.frequency) or mode.frequency <= 1e-5:
                    continue

                if calculate_zp:
                    zp_energy_q += mode.get_zp_energy()
                else:
                    u_vib_q += mode.get_u_vib(temperature=T)
                    entropy_q += mode.get_entropy(temperature=T)
                    c_v_q += mode.get_c_v(temperature=T)

            if calculate_zp:
                zp_energy.append(zp_energy_q)
            else:
                u_vib.append(u_vib_q)
                entropy.append(entropy_q)
                c_v.append(c_v_q)

        if calculate_zp:
            zp_energy = np.array(zp_energy, dtype=float)
            return zp_energy
        else:
            u_vib = np.array(u_vib, dtype=float)
            entropy = np.array(entropy, dtype=float)
            c_v = np.array(c_v, dtype=float)
            return u_vib, entropy, c_v

    def thermodynamics(self, sumphonon=True, mutewarning=False, **kwargs):
        """
        Calculate the thermodynamic properties (zp_energy, u_vib, entropy, c_v
        and Gibbs and Helmholtz free energy) of the HA system at all qpoints
        and the whole temperature/pressure range.

        Other parameters are the sum of corresponding attributes of all the
        ``Mode objects``. The Helmholtz and Gibbs free energies are defined as:

        .. math::

            F(p,V) = E_{DFT} + F_{vib}(T) = E_{DFT} + U_{vib}(T) - TS(T)

            G(p, V) = F + pV

        Args:
            temperature (array[float] | list[float], optional): Unit: K
            pressure (array[float] | list[float], optional): Unit: GPa
            sumphonon (bool): Whether to sum up the phonon contributions across
                the sampled q points and take weighted-average.
            mutewarning (bool): Whether print out warning messages of updated
                temperature and pressure (For QHA).

        Returns:
            self.helmholtz (array[float]): nqpoint\*ntemperature. Unit: KJ/mol
            self.gibbs (array[float]): nqpoint\*nPressure\*nTemperature. Unit: KJ/mol
            self.zp_energy (array[float]): Zero-point energy. nqpoint\*1. Unit: KJ/mol
            self.u_vib (array[float]): Vibrational contribution to internal
                energy. nqpoint\*ntemperature. Unit: KJ/mol
            self.entropy (array[float]): nqpoint\*ntemperature. Unit: J/mol\*K
            self.c_v (array[float]): Constant volume specific heat. 
                nqpoint\*ntemperature. Unit: J/mol\*K

        .. note::

            If ``sumphonon = True``, nqpoint = 1.

        :raise AttributeError: If temperature and pressure are defined neither here nor during initialization
        """
        import warnings
        import numpy as np
        import scipy.constants as scst

        # Generate temperature and pressure series
        if kwargs:
            if 'temperature' in kwargs:
                if hasattr(self, 'temperature') and not mutewarning:
                    warnings.warn('Temperature attribute exists. Input temperatures will be used to update the attribute.')
                self.temperature = np.array(kwargs['temperature'], dtype=float)

            if 'pressure' in kwargs:
                if hasattr(self, 'pressure') and not mutewarning:
                    warnings.warn('Pressure attribute exists. Input pressures will be used to update the attribute.')
                self.pressure = np.array(kwargs['pressure'], dtype=float)
        else:
            if not hasattr(self, 'temperature') or not hasattr(self, 'pressure'):
                raise AttributeError('Temperature and pressure should be specified.')

        zp_energy = self._phonon_sumup(temperature=0., calculate_zp=True)
        u_vib = []
        entropy = []
        c_v = []
        helmholtz = []
        gibbs = []

        for T in self.temperature:
            gibbs_t = []
            u_vib_t, entropy_t, c_v_t = self._phonon_sumup(temperature=T,
                                                           calculate_zp=False)
            helm_t = -entropy_t * T * 1e-3 + u_vib_t + self.edft

            for p in self.pressure:
                gibbs_tp = p * self.volume * scst.Avogadro * 1e-24 + helm_t
                gibbs_t.append(gibbs_tp)

            # nTemp * nqpoint
            u_vib.append(u_vib_t)
            entropy.append(entropy_t)
            c_v.append(c_v_t)
            helmholtz.append(helm_t)
            # nTemp * npress * nqpoint
            gibbs.append(gibbs_t)

        if sumphonon:
            wt = np.array([qp[1] for qp in self.qpoint])
            self.nqpoint = 1
            self.qpoint = [[np.array([0., 0., 0.]), 1.]]
            self.zp_energy = np.array([np.dot(zp_energy, wt)])
            self.u_vib = np.array([np.dot(u_vib, wt)])
            self.entropy = np.array([np.dot(entropy, wt)])
            self.c_v = np.array([np.dot(c_v, wt)])
            self.helmholtz = np.array([np.dot(helmholtz, wt)])
            self.gibbs = np.array([np.dot(gibbs, wt)])
            self.gibbs = np.transpose(self.gibbs, (0, 2, 1))
        else:
            self.zp_energy = zp_energy
            self.u_vib = np.transpose(np.array(u_vib, dtype=float))
            self.entropy = np.transpose(np.array(entropy, dtype=float))
            self.c_v = np.transpose(np.array(c_v, dtype=float))
            self.helmholtz = np.transpose(np.array(helmholtz, dtype=float))
            self.gibbs = np.transpose(np.array(gibbs, dtype=float), (2, 1, 0))

        return self.helmholtz, self.gibbs, self.zp_energy, self.u_vib,\
            self.entropy, self.c_v

    def print_results(self):
        """
        Print HA thermodynamic results into an external file. Used if
        ``write_out = True``.

        Phonon dispersions are forced to be summed if the automatic scheme
        (``write_out=True``) is launched. To get verbose outputs, call
        ``self.thermodynamics()`` first and then call ``self.print_results()``.
        """
        import scipy.constants as scst
        from CRYSTALpytools.units import eV_to_H

        if not self.write_out:
            print('Harmonic.write_out = False, return to empty.')
            return

        file = open(self.filename, 'w')
        file.write('%21s%20.9e%15s%20.12e%s\n' %
                   ('# DFT TOTAL ENERGY = ', eV_to_H(self.edft),
                    ' eV,         = ', self.edft, ' kJ/mol'))
        file.write('%21s%20.4f%15s%20.4f%s\n' %
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
            file.write('%s%20.12e%s\n\n' %
                       ('## ZERO POINT ENERGY = ', self.zp_energy[q], ' kJ/mol'))
            file.write('%s\n\n' % '## TEMPERATURE DEPENDENT PROPERTIES')
            file.write('%8s%20s%20s%20s%20s\n' %
                       ('T(K)', 'U_vib(kJ/mol)', 'Entropy(J/mol*K)',
                        'C_V(J/mol*K)', 'Helmholtz(kJ/mol)'))
            for t, tempt in enumerate(self.temperature):
                file.write('%8.2f%20.12e%20.12e%20.12e%20.12e\n' %
                           (tempt, self.u_vib[q, t], self.entropy[q, t],
                            self.c_v[q, t], self.helmholtz[q, t]))

            file.write('\n')
            for idx_p, gibbs_p in enumerate(self.gibbs[q]):
                file.write('%s%8.2f%s\n\n' % ('## GIBBS FREE ENERGY AT', self.pressure[idx_p], ' GPa'))
                file.write('%8s%20s\n' % ('T(K)', 'Gibbs(kJ/mol)'))
                for idx_t, gibbs_t in enumerate(gibbs_p):
                    file.write('%8.2f%20.12e\n' % (self.temperature[idx_t], gibbs_t))
                file.write('\n')
            file.write('\n')

        file.write('\n')
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
        filename (str): Name of the output file. Valid if ``write_out = True``.

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

    def from_HA_files(self, input_files, scelphono=[], imaginary_tol=-1e-4,
                      q_overlap_tol=1e-4, mode_sort_tol=0.4):
        """
        Read data from individual HA calculation outputs. Imaginary modes and
        overlapped q points are forced to be cleaned.

        Args:
            input_files (list[str]): List of phonon output filenames.
            scelphono (array[float] | list[float]): Corresponds to the
                'SCELPHONO' keyword in CRYSTAL. Either 3\*3 or ndimension\*ndimension.
                By default a 1\*1\*1 'SCELPHONO' is assumed.
            imaginary_tol (float): The threshold of negative frequencies.
            q_overlap_tol (float): The threshold of overlapping points, defined
                as the 2nd norm of the difference of fractional q vectors
            mode_sort_tol (float | None): The threshold of close mode
                overlaps. If none, do not sort modes.

        Returns:
            self (Quasi_harmonic)

        **New Attributes**

        * self.ncalc (int): Number of HA phonon calculations.  
        * self.combined_phonon (list[Harmonic]): List of Harmonic objects.  
        * self.combined_volume (list[float]): Volumes. Unit: Angstrom^3  
        * self.combined_edft (list[float]): DFT total energies. Unit: KJ/mol  
        * self.combined_mode (list[Mode]): List of mode objects.  
        """
        from CRYSTALpytools.thermodynamics import Harmonic
        import warnings

        if hasattr(self, "ncalc"):
            warnings.warn('Data exists. The current command will be ignored.')
            return self

        if len(input_files) == 1:
            raise Exception('Only 1 input file! Use Harmonic object or from_QHA_file method.')
        else:
            self.ncalc = len(input_files)

        if mode_sort_tol != None:
            read_eigvt = True
        else:
            read_eigvt = False

        ha_list = [
            Harmonic(write_out=False, autocalc=False).from_file(
                file,
                scelphono=scelphono,
                read_eigvt=read_eigvt,
                imaginary_tol=imaginary_tol,
                q_overlap_tol=q_overlap_tol
            ) for file in input_files
        ]

        self.combined_phonon, self.combined_volume, self.combined_edft, \
        self.combined_mode = self._combine_data(ha_list, mode_sort_tol=mode_sort_tol)
        self.nqpoint = ha_list[0].nqpoint
        self.qpoint = ha_list[0].qpoint # consistency of nqpoint is checked, but not qpoint.

        return self

    def from_QHA_file(self, input_file, scelphono=[], imaginary_tol=-1e-4,
                      q_overlap_tol=1e-4, mode_sort_tol=0.4):
        """
        Read data from a single QHA calculation at Gamma point. Imaginary modes
        and overlapped q points are forced to be cleaned.

        Args:
            input_files (str | list[str]): Only 1 QHA file is permitted.
            scelphono (array[float] | list[float])
            imaginary_tol (float): The threshold of negative frequencies.
            q_overlap_tol (float): The threshold of overlapping points, defined
                as the 2nd norm of the difference of fractional q vectors
            mode_sort_tol (float | None): The threshold of close mode
                overlaps. If none, do not sort modes.

        Returned attributes are consistent with ``Quasi_harmonic.from_HA_files``.

        :raise ValueError: If multiple files are defined.
        """
        from CRYSTALpytools.thermodynamics import Harmonic
        import warnings
        import re
        import numpy as np
        from pymatgen.core import Structure

        if hasattr(self, "ncalc"):
            warnings.warn('Data exists. The current command will be ignored.')
            return self

        if isinstance(input_file, list) and len(input_file) > 1:
            raise ValueError("Only a single QHA file is permitted")
        elif isinstance(input_file, list) and len(input_file) == 1:
            input_file = input_file[0]

        if mode_sort_tol != None:
            read_eigvt = True
        else:
            read_eigvt = False

        output = Crystal_output().read_cry_output(input_file)
        output.get_phonon(read_eigvt=read_eigvt, rm_imaginary=False, rm_overlap=False)
        strucs = _restore_pcel(output, scelphono)

        self.ncalc = output.nqpoint
        self.nqpoint = 1
        self.qpoint = [[np.array([0., 0., 0.]), 1.]]

        ha_list = []
        for idx_c in range(self.ncalc):
            ha = Harmonic(write_out=False, autocalc=False)
            if read_eigvt == True:
                ha.from_frequency(output.edft[idx_c], [[np.zeros([3,]), 1.]],
                                  np.array([output.frequency[idx_c],]),
                                  np.array([output.eigenvector[idx_c],]),
                                  structure=strucs[idx_c])
            else:
                ha.from_frequency(output.edft[idx_c], [[np.zeros([3,]), 1.]],
                                  np.array([output.frequency[idx_c],]),
                                  [], structure=strucs[idx_c])
            ha_list.append(ha)

        self.combined_phonon, self.combined_volume, self.combined_edft, \
            self.combined_mode = self._combine_data(ha_list, mode_sort_tol)

        return self

    def from_phonopy_files(self, phono_yaml, struc_yaml=None, edft=None,
                           imaginary_tol=-1e-4, q_overlap_tol=1e-4,
                           q_id=None, q_coord=None):
        """
        Build a QHA object from `Phonopy <https://phonopy.github.io/phonopy/>`_
        'band.yaml' or 'qpoints.yaml' file.

        Args:
            phono_yaml (list[str]): ncalc\*1 list of Phonopy band.yaml or
                qpoint.yaml files
            struc_yaml (list[str]): ncalc\*1 list of Phonopy phonopy.yaml or
                phonopy_disp.yaml files. *Needed only if a qpoint.yaml file is
                read.*
            edft (list[float]): ncalc\*1 list / array of DFT energies.
            imaginary_tol (float): The threshold of negative frequencies.
            q_overlap_tol (float): The threshold of overlapping points, defined
                as the 2nd norm of the difference of fractional q vectors
            q_id (list[int]): See ``Harmonic.from_phonopy``.
            q_coord (list[list]): See ``Harmonic.from_phonopy``.

        .. note::

            ``q_id`` and ``q_coord`` should be set once and are applicable to
            all the yaml files.

        Returned attributes are consistent with ``Quasi_harmonic.from_HA_files``.
        """
        import numpy as np
        from CRYSTALpytools.thermodynamics import Harmonic
        import warnings

        if hasattr(self, "ncalc"):
            warnings.warn('Data exists. The current command will be ignored.')
            return self

        if len(input_files) == 1:
            raise Exception('Only 1 input file! Use Harmonic object or from_QHA_file method.')
        else:
            self.ncalc = len(input_files)

        if edft == None:
            warnings.warn('DFT energy is set to 0.')
            edft = np.zeros([self.ncalc,])

        if struc_yaml == None:
            struc_yaml = [None for i in range(self.ncalc)]

        ha_list = [
            Harmonic(write_out=False, autocalc=False).from_phonopy(
                phono_yaml=phono_yaml[i],
                struc_yaml=struc_yaml[i],
                edft=edft[i],
                imaginary_tol=imaginary_tol,
                q_overlap_tol=q_overlap_tol,
                q_id=q_id,
                q_coord=q_coord
            ) for i in range(self.ncalc)
        ]

        self.combined_phonon, self.combined_volume, self.combined_edft, \
        self.combined_mode = self._combine_data(ha_list, mode_sort_tol=None) # Eigenvector not available
        self.nqpoint = ha_list[0].nqpoint
        self.qpoint = ha_list[0].qpoint # consistency of nqpoint is checked, but not qpoint.

        return self

    def _combine_data(self, ha_list, mode_sort_tol):
        """
        Combine the HA calculation data and rearrange it in the ascending order
        of volumes.

        Args:
            ha_list (list[Harmonic]): List of harmonic objects.
            mode_sort_tol (float | None)

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
            if (natom - ha_phonon.natom) != 0 or not np.all((nmode - ha_phonon.nmode) == 0) \
            or nqpoint - ha_phonon.nqpoint != 0:
                raise Exception('The number of qpoints, modes or atoms is not consistent across the sampling points')

        sorted_vol = sorted_vol[np.argsort(sorted_vol[:, 1])]
        nmode = nmode[0]
        if ha_list[0].eigenvector == []:
            do_eigvt = False
        else:
            do_eigct = True

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
            if do_eigvt == True:
                combined_eigvt[idx_new] = ha_phonon.eigenvector

        # ncalc * nqpoint * nmode array to nqpoint * ncalc * nmode array
        combined_freq = np.transpose(combined_freq, axes=[1, 0, 2])
        if do_eigvt == True:
            # ncalc * nqpoint * nmode * natom * 3 array to nqpoint * ncalc * nmode * natom * 3 array
            combined_eigvt = np.transpose(combined_eigvt, axes=[1, 0, 2, 3, 4])

        # Sort phonon modes if requested
        if mode_sort_tol != None and do_eigvt == True:
            close_overlap = np.zeros([nqpoint, self.ncalc, nmode, nmode])
            for idx_q in range(nqpoint):
                combined_freq[idx_q], combined_eigvt[idx_q], close_overlap[idx_q] \
                    = self._phonon_continuity(combined_freq[idx_q],
                                              combined_eigvt[idx_q],
                                              mode_sort_tol=mode_sort_tol)
            # nqpoint * ncalc * nmode_ref * nmode_sort array to
            # nqpoint * nmode_ref * ncalc * nmode_sort array
            close_overlap = np.transpose(close_overlap, axes=[0, 2, 1, 3])
            for q, overlap_q in enumerate(close_overlap):
                overlap_numbers = np.sum(overlap_q)
                if overlap_numbers > 0:
                    warnings.warn('Close overlap of phonon modes detected at qpoint: %3i, %6i overlaps out of %6i modes.'
                                  % (q, int(overlap_numbers), int(nqpoint * nmode)),
                                  stacklevel=2)
        elif mode_sort_tol != None and do_eigvt == False:
            warnings.warn('Eigenvectors not read. Mode sorting not available.')

        # nqpoint * ncalc * nmode array to nqpoint * nmode * ncalc array
        combined_freq = np.transpose(combined_freq, axes=[0, 2, 1])
        if do_eigvt == True:
            # nqpoint * ncalc * nmode * natom * 3 array to nqpoint *  nmode * ncalc * natom * 3 array
            combined_eigvt = np.transpose(combined_eigvt, axes=[0, 2, 1, 3, 4])

        combined_mode = []
        for idx_q in range(nqpoint):
            combined_mode_q = []
            for idx_m in range(nmode):
                if do_eigvt == True:
                    combined_mode_q.append(
                        Mode(rank=idx_m + 1,
                             frequency=combined_freq[idx_q, idx_m, :],
                             volume=combined_volume,
                             eigenvector=combined_eigvt[idx_q, idx_m, :])
                    )
                else:
                    combined_mode_q.append(
                        Mode(rank=idx_m + 1,
                             frequency=combined_freq[idx_q, idx_m, :],
                             volume=combined_volume)
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

            if mode_sort_tol != None and do_eigvt == True::
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
                                        '%10i%2s%8i%2s%9i%2s%9i\n' %
                                        (idx_mref + 1, '', idx_csort - 1, '', idx_csort, '', idx_msort + 1)
                                    )
                                else:
                                    continue

                    file.write('\n')

            file.close()

        return combined_phonon, combined_volume, combined_edft, combined_mode

    @staticmethod
    def _phonon_continuity(freq, eigvt, symm=None, mode_sort_tol=0.4):
        """
        Rearrange phonon modes by their continuity. If the difference between
        the maximum scalar product of corresponding eigenvectors (normalized to
        1) and scalar products of other modes is less than 0.4, warning is
        printed due to the potential overlap of modes. Adopted from CRYSTAL17.

        .. note::

            A. Erba, *J. Chem. Phys.*, 2014, **141**, 124115.

        Args:
            freq (array[float]): Phonon frequencies. Unit: THz
            eigvt (array[float]): Eigenvectores normalized to 1
            symm (array[float]): Sub-group numbers of corresponding modes.
                *Not implemented*
            mode_sort_tol (float): The threshold of close mode overlaps.

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

    def eos_fit(self, volume, energy, method, write_out=True, **kwargs):
        """
        Fit energy-volume relationship by equation of states.

        Args:
            volume (array[float]): Unit: Angstrom^3
            energy (array[float]): Unit: kJ/mol
            method (str): Name of EoS used. Consistent with
                `Pymatgen <https://pymatgen.org/pymatgen.analysis.eos.html>`_.
            write_out (bool): Whether to print EOS information.
            order (int): For DeltaFactor / Polynomial methods.
            min_ndata_factor, max_poly_order_factor, min_poly_order_factor (int):
                For the NumericalEOS method.

        Returns:
            eos (Pymatgen EOS): The fitted equation of state.
            eos_method (string): Name of the fitted equation of state
        """
        import re
        from pymatgen.analysis.eos import Murnaghan, Birch, BirchMurnaghan, \
            PourierTarantola, Vinet, DeltaFactor, NumericalEOS, PolynomialEOS
        import scipy.constants as scst


        eos_method = method
        classes = {
            "murnaghan"         : Murnaghan,
            "birch"             : Birch,
            "birch_murnaghan"   : BirchMurnaghan,
            "pourier_tarantola" : PourierTarantola,
            "vinet"             : Vinet,
            "deltafactor"       : DeltaFactor,
            "numerical_eos"     : NumericalEOS,
            "polynomial"        : PolynomialEOS,
        }
        eos = classes[method](volume, energy)
        eos_command = 'eos.fit('
        # Polynomial / Deltafactor / Numerical
        for idx, key in enumerate(kwargs.keys()):
            value = list(kwargs.values())[idx]
            eos_command += ', {}={}'.format(key, value)
        eos_command += ')'
        eval(eos_command)

        if self.write_out and write_out == True:
            file = open(self.filename, 'a+')
            file.write('%s%s\n' % ('# EQUATION OF STATES FITTED FOR ELECTRON TOTAL ENERGY: ', method))
            file.write('%s\n' % '  Electron total energy is fitted as the function of volume, of which the')
            file.write('%s\n\n' % '  formalism is given by equation of states.')
            file.write('%16s%16s%12s%12s\n' % 
                       ('E0(kJ/mol)', 'V0(Angstrom^3)', 'B0(GPa)', 'B1'))
            file.write('%16.4f%16.4f%12.4f%12.4f\n' %
                       (eos.e0, eos.v0, eos.b0 * 1e24 / scst.Avogadro, eos.b1))
            file.write('\n')
            file.close()

        return eos, method

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
                order_new, _, _ = mode.polynomial_fit(order=order)
                for key, value in mode.poly_fit_rsqaure.items():
                    rsquare_q[key] += value / len(mode_q)

                if self.write_out:
                    file.write('%-8s%7s%14s%s\n' %
                               ('  Mode #', 'Order', 'R^2', '  Coeff low to high (Constant term = 0)'))
                    for idx_od, od in enumerate(order_new):
                        if idx_od == 0:
                            file.write('%8i' % mode.rank)
                        else:
                            file.write('%8s' % '')

                        file.write('%7i%2s%12.6f%2s' %
                                   (od, '', mode.poly_fit_rsqaure[od], ''))
                        for idx_c, c in enumerate(mode.poly_fit[od].convert().coef):
                            if idx_c == 0:
                                continue
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

        if not hasattr(self, 'fit_order') or not hasattr(self, 'e0_eos'):
            raise Exception('ERROR: Analytical expressions unavailable.')

        num_freq = []
        for mode_q in self.combined_mode:
            num_freq_q = []
            for mode in mode_q:
                idx_vmin = np.argmin(mode.volume)
                vmin = mode.volume[idx_vmin]
                fmin = mode.frequency[idx_vmin]
                dv = volume - vmin
                num_freq_q.append(mode.poly_fit[self.fit_order](dv) + fmin)
            num_freq.append(num_freq_q)

        num_freq = np.array(num_freq)
        ha = Harmonic(write_out=False).from_frequency(
            self.e0_eos(volume), self.qpoint, num_freq, [], volume=volume)

        return ha

    def _minimize_gibbs(self, volume, temperature, pressure):
        """
        Get Gibbs free energy from the Harmonic phonon object. Used only for
        minimizing :math:`G(V;T, p)` by SciPy.

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

    @staticmethod
    def _poly_no_cst(param, x, y):
        """
        Define a polynomial :math:`\\Delta f(\\Delta x)` without constant term.
        Orders low to high. For SciPy.
        """
        import numpy as np

        express = np.zeros([len(x)])
        for order, p in enumerate(param):
            express += p * x**(order + 1)
        return ((express - y)**2)**0.5

    def _clean_attr(self):
        """
        When temperature / pressure are changed, thermodynamic attributes are
        removed to keep consistency.
        """
        attr_list = ['volume', 'helmholtz', 'gibbs', 'entropy', 'c_v', 'c_p',
                     'k_t', 'k_s', 'fe_eos_method', 'fe_eos', 'gruneisen',
                     'alpha_vgru', 'c_pgru', 'k_sgru', 'alpha_v', 'vol_fit']

        for attr in attr_list:
            if hasattr(self, attr):
                delattr(self,attr)

        return self

    def thermo_freq(self, eos_method='birch_murnaghan', poly_order=[2, 3],
                    min_method='BFGS', volume_bound=None, mutewarning=False,
                    **kwargs):
        """
        Obtain thermodynamic properties by explicitly fitting phonon
        frequencies as polynomial functions of volume. DFT total energies are
        fitted as a function of volume by equation of states (EOS).

        The equilibrium volume is fitted by minimizing Gibbs free energy at
        constant temperature and pressure.

        .. math::

            V(T,p)=\\text{min}[G(V;T,p)]=\\text{min}[E_{0}(V)+F_{vib}(V;T,p)+pV)]

        Args:
            eos_method (str, optional): EOS used to fit DFT total energy and
                Helmholtz free energy (to get bulk modules).
            poly_order (array[int] | list[int], optional): The order of
                polynomials used to fit frequency as the function of volumes.
            min_method (string, optional): Minimisation algorithms.
            volume_bound (tuple-like, optional), Boundary conditions of
                equilibrium volumes. Unit: Angstrom^3
            mutewarning (bool, optional): Whether print out warning messages.
            temperature (array[float], optional): Unit: K
            pressure (array[float], optional): Unit: GPa
            order (int, optional): For DeltaFactor / Polynomial EOSs.
            min_ndata_factor, max_poly_order_factor, min_poly_order_factor (int, optional):
                For Numerical EOS.

        .. note::

            #. Valid entries of ``eos_method`` are consistent with `PyMatGen <https://pymatgen.org/pymatgen.analysis.eos.html>`_.
            #. Parameterized and tested algorithms for ``min_method``:
                * BFGS(no boundary)
                * L-BFGS-B(with boundary)

        Returns:
            self (Quasi_harmonic)

        **New Attributes**

        * ``self.temperature`` in K and ``self.pressure`` in GPa.  
        * ``self.volume``, nPressure\*nTemperature. Equilibrium volumes. Unit: Angstrom^3  
        * ``self.helmholtz`` and ``self.gibbs``, nPressure\*nTemperature. Helmholtz and Gibbs free energies. Unit: kJ/mol  
        * ``self.entropy``, nPressure\*nTemperature, Entropy. Unit: J/mol\*K  
        * ``self.c_v``, nPressure\*nTemperature, Constant volume specific heat. Unit: J/mol\*K  
        * ``self.e0_eos`` and ``self.e0_eos_method`` Pymatgen EOS objects and string. EOS used to fit DFT energy.  

        :raise ValueError: If temperature or pressure is defined neither here nor during initialization.
        """
        import numpy as np
        import warnings
        from scipy.optimize import minimize

        # Generate temperature and pressure series
        if 'temperature' in kwargs:
            if hasattr(self, 'temperature') and not mutewarning:
                warnings.warn('Temperature attribute exists. Input temperatures will be used to update the attribute.',
                              stacklevel=2)
            self.temperature = np.array(kwargs['temperature'], dtype=float)
            self._clean_attr()

        if 'pressure' in kwargs:
            if hasattr(self, 'pressure') and not mutewarning:
                warnings.warn('Pressure attribute exists. Input pressures will be used to update the attribute.',
                              stacklevel=2)
            self.pressure = np.array(kwargs['pressure'], dtype=float)
            self._clean_attr()

        if not hasattr(self, 'temperature') or not hasattr(self, 'pressure'):
            raise ValueError('Temperature and pressure should be specified.')

        # Fit DFT total energy, if not done yet. Otherwise, fitted values will not be covered.
        if hasattr(self, 'e0_eos') and not mutewarning:
            warnings.warn('DFT total energy is already fitted. To keep the consistency, it will not be updated.',
                          stacklevel=2)
        else:
            eos_method = eos_method.casefold()
            eos_command = 'self.eos_fit(self.combined_volume, self.combined_edft, eos_method'
            # Polynomial / Deltafactor / Numerical
            for idx, key in enumerate(kwargs.keys()):
                if key == 'temperature' or key == 'pressure':
                    continue
                value = list(kwargs.values())[idx]
                eos_command += ', {}={}'.format(key, value)
            eos_command += ')'
            self.e0_eos, self.e0_eos_method = eval(eos_command)
        # Fit frequencies, if not done yet. Otherwise, fitted values will not be covered.
        if hasattr(self, 'fit_order') and not mutewarning:
            warnings.warn('Frequency is already fitted to polynomials. To keep the consistency, it will not be updated.',
                          stacklevel=2)
        else:
            self.freq_polynomial_fit(order=poly_order)

        # Define minimization methods
        methods = {
            'BFGS': "vol = minimize(self._minimize_gibbs, v_init, args=(t, p), method='BFGS', jac='3-point')",
            'L-BFGS-B': "vol = minimize(self._minimize_gibbs, v_init, args=(t, p), method='L-BFGS-B', jac='3-point', bounds=volume_bound)",
        }

        # Gibbs(V; T, p) minimization nPress*nTempt list
        self.volume = np.zeros([len(self.pressure), len(self.temperature)])
        v_init = np.mean(self.combined_volume)

        for idx_p, p in enumerate(self.pressure):
            for idx_t, t in enumerate(self.temperature):
                params = {'self': self,
                          'minimize': minimize,
                          'v_init': v_init,
                          't': t,
                          'p': p,
                          'volume_bound': volume_bound}
                exec(methods[min_method], params)
                self.volume[idx_p, idx_t] = params['vol'].x[0]

                if (params['vol'].x[0] < min(self.combined_volume) or params['vol'].x[0] > max(self.combined_volume)) and not mutewarning:
                    warnings.warn('Optimised volume exceeds the sampled range. Special care should be taken of.\n  Volume: %12.4f, Temperature: %6.2f, Pressure: %6.2f\n'
                                  % (params['vol'].x[0], t, p), stacklevel=2)

        # Calculate other thermodynamic properties
        self.helmholtz = np.zeros(self.volume.shape)
        self.gibbs = np.zeros(self.volume.shape)
        self.entropy = np.zeros(self.volume.shape)
        self.c_v = np.zeros(self.volume.shape)
        for idx_p, p in enumerate(self.pressure):
            for idx_t, t in enumerate(self.temperature):
                vol = self.volume[idx_p, idx_t]
                ha = self._get_harmonic_phonon(vol)
                ha.thermodynamics(temperature=[t], pressure=[p], mutewarning=True)
                self.helmholtz[idx_p, idx_t] = ha.helmholtz[0, 0]
                self.gibbs[idx_p, idx_t] = ha.gibbs[0, 0, 0]
                self.entropy[idx_p, idx_t] = ha.entropy[0, 0]
                self.c_v[idx_p, idx_t] = ha.c_v[0, 0]

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA THERMODYNAMIC PROPERTIES - FREQUENCY')
            file.write('%s\n\n' % '  QHA thermodynamic properties by explicitly fitting frequencies.')
            file.write('%s%6i\n' %
                       ('## FREQUENCY POLYNOMIAL ORDER: ', self.fit_order))
            file.write('%s%s\n' %
                       ('## EQUILIBRIUM VOLUME MINIMISATION: ', min_method))
            file.write('%s%s\n' %
                       ('## HELMHOLTZ FREE ENERGY EOS: ', eos_method))
            if volume_bound:
                file.write('%s\n' % (
                    '## CONSTRAINED VOLUME MINIMIZATION LAUNCHED. VOLUME BOUNDARIES (UNIT: ANGSTROM^3):'))
                file.write('%s%8.2f%s%8.2f\n\n' % (
                    '## LOWER: ', volume_bound[0], ' UPPER: ', volume_bound[1]))

            for idx_p, press in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' %
                           ('## THERMODYNAMIC PROPERTIES AT ', press, '  GPa'))
                file.write('%10s%20s%20s%20s%20s%20s\n' %
                           ('T(K)', 'Vol(Angstrom^3)', 'Helmholtz(kJ/mol)',
                            'Gibbs(kJ/mol)', 'Entropy(J/mol*K)', 'C_V(J/mol*K)'))
                for idx_t, tempt in enumerate(self.temperature):
                    file.write('%10.2f%20.4f%20.8e%20.8e%20.8e%20.8e\n' %
                               (tempt, self.volume[idx_p, idx_t],
                                self.helmholtz[idx_p, idx_t],
                                self.gibbs[idx_p, idx_t],
                                self.entropy[idx_p, idx_t],
                                self.c_v[idx_p, idx_t]))

                file.write('\n')

            file.write('\n')
            file.close()

        return self

    def thermo_gruneisen(self, eos_method='birch_murnaghan', min_method='BFGS',
                         volume_bound=None, mutewarning=False, **kwargs):
        """
        Grüneisen parameters and related properties. The macroscopic Grüneisen
        parameter is defined as:

        .. math::

            \\gamma=\\sum_{\\textbf{q}i}\\frac{\\gamma_{\\textbf{q}i}C_{V,\\textbf{q}i}}{C_{V}}

        Thermal expansion coefficient in Grüneisen model:

        .. math::

            \\alpha_{V}^{gru}=\\frac{\\gamma C_{V}}{K_{T}V}

        .. note::

            The Grüneisen model is used to fit frequencies, equivalent to using
            ``self.thermo_freq(poly_order=[1,])``.

        For arguments, see ``self.thermo_freq``.

        Returns:
            self.gamma(array): npressure\*ntemperature array of macroscopic
                Grüneisen parameter. Temperature should > 0.
        """
        import numpy as np
        import scipy.constants as scst
        import warnings

        if hasattr(self, 'fit_order'):
            raise AttributeError('self.gruneisen cannot be used when self.thermo_freq is already used.')

        command = 'self.thermo_freq(eos_method=eos_method, poly_order=[1,], min_method=min_method, volume_bound=volume_bound, mutewarning=mutewarning'

        for idx, key in enumerate(kwargs.keys()):
            value = list(kwargs.values())[idx]
            if type(value) == np.ndarray:
                value = list(value)
            command += ', {}={}'.format(key, value)
        command += ')'
        eval(command)

        # Get mode-specific Grüneisen parameter
        for idx_q, mode_q in enumerate(self.combined_mode):
            for idx_m, mode in enumerate(mode_q):
                self.combined_mode[idx_q][idx_m].get_gruneisen(order=[1,], volume=self.volume)

        # Macroscopic Grüneisen parameter
        sum_gCv = np.zeros(self.volume.shape, dtype=float)
        for idx_q, mode_q in enumerate(self.combined_mode):
            for idx_m, mode in enumerate(mode_q):
                c_v = np.zeros(self.volume.shape, dtype=float)
                if idx_m == 0 or idx_m == 1 or idx_m == 2:
                    continue
                # Get matrix C_v, nTempt*nPress
                idx_vmin = np.argmin(mode.volume)
                vmin = mode.volume[idx_vmin]
                fmin = mode.frequency[idx_vmin]
                dv = self.volume - vmin
                for idx_t, t in enumerate(self.temperature):
                    if t > 1e-4: # > 0K
                        kb_t = scst.k * scst.Avogadro * t
                        hbar_freq = (fmin + mode.poly_fit[1](dv[:, idx_t])) * scst.Avogadro * scst.h * 1e12
                        expon = np.exp(hbar_freq / kb_t)
                        c_v[:, idx_t] = hbar_freq**2 / kb_t / t * expon / (expon - 1)**2
                    else:
                        c_v[:, idx_t] = 0.
                sum_gCv += c_v * mode.gruneisen[1]

        # Get K_T, EOS keywords
        command = 'self.bulk_modulus(adiabatic=False'
        for idx, key in enumerate(kwargs.keys()):
            if key == 'temperature' or key == 'pressure':
                continue
            command += ', {}={}'.format(key, list(kwargs.values())[idx])
        command += ')'
        eval(command)

        self.alpha_vgru = sum_gCv / self.k_t / self.volume / 1e-21 / scst.Avogadro
        self.c_pgru = self.c_v + self.alpha_vgru**2 * self.k_t * self.volume * self.temperature * 1e-21 * scst.Avogadro

        self.gruneisen = np.zeros(self.volume.shape)
        self.k_sgru = np.zeros(self.volume.shape)
        for idx_t, t in enumerate(self.temperature):
            if t < 1e-4: # 0K
                continue
            self.gruneisen[:, idx_t] = sum_gCv[:, idx_t] / self.c_v[:, idx_t]
            self.k_sgru[:, idx_t] = self.k_t[:, idx_t] + \
                self.alpha_vgru[:, idx_t]**2 * self.volume[:, idx_t] * t * self.k_t[:, idx_t]**2 * 1e-21 * scst.Avogadro / self.c_v[:, idx_t]

        # print out options
        if self.write_out == True:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA THERMODYNAMIC PROPERTIES - GRÜENEISEN MODEL')
            file.write('%s\n\n' % '  Linear dependency of frequency with volume is assumed.')
            for idx_p, p in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' % ('## GRÜENEISEN THERMODYNAMICS AT ', p, '  GPa'))
                file.write('%10s%10s%20s%20s%20s%20s%20s\n' %
                           ('T(K)', 'GRÜ PARAM','alpha_VGRÜ(K^-1)',
                            'C_v(J/mol*K)', 'C_pGRÜ(J/mol*K)',
                            'K_T(GPa)', 'K_SGRÜ(GPa)'))
                for idx_t, t in enumerate(self.temperature):
                    file.write('%10.1f%10.4f%20.8e%20.8e%20.8e%20.8e%20.8e\n' %
                               (t, self.gruneisen[idx_p, idx_t], self.alpha_vgru[idx_p, idx_t],
                                self.c_v[idx_p, idx_t], self.c_pgru[idx_p, idx_t],
                                self.k_t[idx_p, idx_t], self.k_sgru[idx_p, idx_t]))
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

        Entropy is obtained by taking the derivation of Gibbs free energy at
        constant pressure.

        .. math::

            S=-\\left(\\frac{\\partial G}{\\partial T}\\right)_{p}

        Constant pressure specific heat is obtained by taking the second
        derivative of :math:`G`.

        .. math::

            C_{p}=-T\\left(\\frac{\\partial^{2}G}{\\partial T^{2}}\\right)_{p}

        .. note::

            ``poly_order`` should >= 2.

        Returns:
            self (Quasi_harmonic)

        **New attributes**

        * ``self.c_p``, nPressure\*nTemperature, Constant pressure specific heat. Unit: J/mol\*K  
        * ``self.fe_eos`` and ``self.fe_eos_method`` nTemperature\*1 list of Pymatgen EOS objects and string. EOSs used to fit HA free energy at constant temperature.

        For arguments and other attributes, see ``Quasi_harmonic.thermo_freq``.

        :raise Exception: If the number of HA calculations is less than 4.
        :raise ValueError: If temperature or pressure is defined neither here nor during initialization.
        """
        import numpy as np
        import warnings
        import re
        from scipy.optimize import fmin, least_squares
        import scipy.constants as scst
        from sympy import diff, lambdify, symbols

        # Check the number of calculations
        if self.ncalc < 4:
            raise Exception('Insufficient database. Increase HA phonons')

        # Generate temperature and pressure series
        if 'temperature' in kwargs:
            if hasattr(self, 'temperature') and not mutewarning:
                warnings.warn('Temperature attribute exists. Input temperatures will be used to update the attribute.',
                              stacklevel=2)
            self.temperature = np.array(kwargs['temperature'], dtype=float)
            self._clean_attr()

        if 'pressure' in kwargs:
            if hasattr(self, 'pressure') and not mutewarning:
                warnings.warn('Pressure attribute exists. Input pressures will be used to update the attribute.',
                              stacklevel=2)
            self.pressure = np.array(kwargs['pressure'], dtype=float)
            self._clean_attr()

        if not hasattr(self, 'temperature') or not hasattr(self, 'pressure'):
            raise ValueError('Temperature and pressure should be specified.')

        # Get data for fitting. Helmholtz: nTempt*nCalc matrix
        helmholtz = np.zeros([len(self.temperature), self.ncalc], dtype=float)
        for idx_c, calc in enumerate(self.combined_phonon):
            hfe_c, _, _, _, _, _ = calc.thermodynamics(
                sumphonon=True, mutewarning=True,
                temperature=self.temperature, pressure=[0.])
            helmholtz[:, idx_c] = hfe_c

        # Fit EOS
        eos_method = eos_method.casefold()
        if hasattr(self, 'fe_eos') and not mutewarning:
            warnings.warn('Harmonic free energy EOS is fitted. To keep the consistency, it will not be updated.',
                          stacklevel=2)
        else:
            self.fe_eos_method = eos_method
            self.fe_eos = []
            for idx_t, t in enumerate(self.temperature):
                eos_command = 'self.eos_fit(self.combined_volume, helmholtz[idx_t, :], eos_method, write_out=False'
                # Polynomial / Deltafactor / Numerical
                for idx, key in enumerate(kwargs.keys()):
                    if key == 'temperature' or key == 'pressure':
                        continue
                    value = list(kwargs.values())[idx]
                    eos_command += ', {}={}'.format(key, value)
                eos_command += ')'
                eos, _ = eval(eos_command)
                self.fe_eos.append(eos)
        # Get thermoproperties
        self.volume = np.zeros([len(self.pressure), len(self.temperature)])
        self.helmholtz = np.zeros(self.volume.shape)
        self.gibbs = np.zeros(self.volume.shape)
        self.entropy = np.zeros(self.volume.shape)
        self.c_p = np.zeros(self.volume.shape)
        v = symbols('v')
        for idx_t, eos in enumerate(self.fe_eos):
            p_eos = -diff(eos(v), v, 1)
            for idx_p, p in enumerate(self.pressure):
                p_kj = p * scst.Avogadro / 1e24  # GPa --> kJ/mol.Angstrom^3
                lam_p = lambdify(v, (p_eos - p_kj)**2, 'numpy')
                fit = fmin(lam_p, eos.v0, full_output=True, disp=False)
                if np.isnan(fit[0]) == True:
                    raise ValueError('EOS fitting failed at %6.2f K, %6.2f GPa. More sampling points needed.' % (self.temperature[idx_t], p))
                if (fit[0] < min(self.combined_volume) or fit[0] > max(self.combined_volume)) and not mutewarning:
                    warnings.warn('Optimised volume exceeds the sampled range. Special care should be taken of.\n  Volume: %12.4f, Temperature: %6.2f, Pressure: %6.2f\n'
                                  % (fit[0], t, p), stacklevel=2)
                self.volume[idx_p, idx_t] = fit[0]
                self.helmholtz[idx_p, idx_t] = eos(fit[0])
                self.gibbs[idx_p, idx_t] = eos(fit[0]) + p_kj * fit[0]

        # Second fit G(T; p), get entropy and C_p
        if max(poly_order) > len(self.temperature) - 1 and not mutewarning:
            warnings.warn('Temperature series not sufficient for the order of polynomial fitting.\n Too high values will be removed.\n',
                          stacklevel=2)

        poly_order = list(set(poly_order))
        poly_order = [p for p in poly_order if p <= len(self.temperature) - 1]

        idx_tmin = np.argmin(self.temperature)
        tmin = self.temperature[idx_tmin]
        dt = self.temperature - tmin
        for idx_p, gibbs in enumerate(self.gibbs):
            r_square = []
            func = []
            for order in poly_order:
                if order < 2:
                    warnings.warn('The minimum order of polynomial is 2. Skip this entry.')
                    continue
                gmin = gibbs[idx_tmin]
                dg = gibbs - gmin
                opt = least_squares(self._poly_no_cst,
                                    np.array([1. for i in range(order)]),
                                    args=(dt, dg))
                poly = np.polynomial.polynomial.Polynomial(np.insert(opt.x, 0, 0.))
                func.append(poly)
                r_square.append(1 - np.sum((dg - poly(dt))**2) / np.sum((dg - np.mean(dg))**2))

            entropy = func[np.argmax(r_square)].deriv(1)
            self.entropy[idx_p, :] = -entropy(dt) * 1000.
            c_p = func[np.argmax(r_square)].deriv(2)
            self.c_p[idx_p, :] = -c_p(dt) * 1000 * self.temperature

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA THERMODYNAMIC PROPERTIES - EOS FIT')
            file.write('%s\n\n' % '  Thermodynamic properties obtained by overall fitting of equation of states.')
            file.write('%s%s\n' % ('## EQUATION OF STATES: ', eos_method))
            file.write('%s%i\n' % ('## G(T) POLYNOMIAL ORDER: ', poly_order[np.argmax(r_square)]))
            file.write('%s\n' % '  WARNING: Entropy at low temperature is probably inaccurate due to the poor fitting of G(T) near 0K.')
            for idx_p, press in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' % ('## THERMODYNAMIC PROPERTIES AT ', press, '  GPa'))
                file.write('%10s%20s%20s%20s%20s%20s\n' %
                           ('T(K)', 'Vol(Angstrom^3)', 'Helmholtz(kJ/mol)', 
                            'Gibbs(kJ/mol)', 'Entropy(J/mol*K)', 'C_p(J/mol*K)'))
                for idx_t, tempt in enumerate(self.temperature):
                    file.write('%10.1f%20.4f%20.8e%20.8e%20.8e%20.8e\n' %
                               (tempt, self.volume[idx_p, idx_t],
                                self.helmholtz[idx_p, idx_t],
                                self.gibbs[idx_p, idx_t],
                                self.entropy[idx_p, idx_t],
                                self.c_p[idx_p, idx_t]))

                file.write('\n')

            file.write('\n')
            file.close()

        return self

    def expansion_vol(self, poly_order=[2, 3], plot=True, fit_fig='expansion_fit.png'):
        """
        Fit the thermal expansion curve and get thermal expansion coefficients
        at equilibrium volumes.

        The volumetric thermal expansion coefficient at constant pressure:

        .. math::

            \\alpha_{V}(T) = \\frac{1}{V(T)}\\left(\\frac{\\partial V(T)}{\\partial T}\\right)_{p}

        Args:
            poly_order (list[int]): *method = 'polynomial'*, order of polynomials.
            plot (bool): Plot V-T curves to examine the goodness of fitting. An
                interactive window will pump out to let user to specify the
                optimial fitting.
            fit_fig (str): File name for fittings. A temperal figure is printed
                to help the user choose the optimal fitting.

        Returns:
            self (Quasi_harmonic)

        **New attributes**

        * ``self.vol_fit`` nPressure\*1 list of Numpy object, the fitted volume V(T)  
        * ``self.alpha_v`` nPressure\*nTemperature array, expansion coefficients at equilibrium volumes
        """
        import numpy as np
        from scipy.optimize import least_squares
        import matplotlib.pyplot as plt
        import warnings

        if not hasattr(self, 'volume'):
            raise AttributeError('Equilibrium volume should be fit first.')

        poly_order = np.array(poly_order)
        # Polynomial fitting
        func = []
        rs = []
        idx_tmin = np.argmin(self.temperature)
        tmin = self.temperature[idx_tmin]
        dt = self.temperature - tmin
        for idx_p, v_p in enumerate(self.volume):
            func_p = []
            rs_p = []
            vmin = v_p[idx_tmin]
            dv = v_p - vmin
            for order in poly_order:
                opt = least_squares(self._poly_no_cst,
                                    np.array([1. for i in range(order)]),
                                    args=(dt, dv))
                poly = np.polynomial.polynomial.Polynomial(np.insert(opt.x, 0, 0.))
                r_square = 1 - np.sum((dv - poly(dt))**2) / np.sum((dv - np.mean(dv))**2)
                func_p.append(poly)
                rs_p.append(r_square)
            func.append(func_p)
            rs.append(rs_p)

        # Get the optimal fit
        rs = np.array(rs) # npress * npolyorder
        rs_mean = np.array([np.mean(rs[:, i]) for i in range(len(poly_order))])
        if plot == False:
            fit_order_idx = np.argmin(rs_mean)
            fit_order = poly_order[fit_order_idx]
        else:
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            cmap = np.vstack([np.linspace([0, 0, 1], [0, 1, 1], 20, endpoint=False),
                              np.linspace([0, 1, 1], [0, 1, 0], 20, endpoint=False),
                              np.linspace([0, 1, 0], [1, 1, 0], 20, endpoint=False),
                              np.linspace([1, 1, 0], [1, 0, 0], 21, endpoint=True)])
            for idx_p, v_p in enumerate(self.volume):
                ax.scatter(self.temperature, v_p, color='k', marker='D', s=40)
                vmin = v_p[idx_tmin]
                dv = v_p - vmin
                for idx_i, i in enumerate(poly_order):
                    t_interp = np.linspace(self.temperature.min(), self.temperature.max(), 1000)
                    c = cmap[int(idx_i / len(poly_order) * 101)]
                    if idx_p == 0:
                        txt = 'Order {:d}, R^2 {:.4f}'.format(i, rs_mean[idx_i])
                        ax.plot(t_interp, vmin + func[idx_p][idx_i](t_interp - tmin), color=c, label=txt)
                    else:
                        ax.plot(t_interp, vmin + func[idx_p][idx_i](t_interp - tmin), color=c)

            ax.legend(loc='lower right')
            fig.savefig(fname=fit_fig, dpi=200)
            # Choose optimal fit
            fit_order = input('Set the optimal fit: ')
            fit_order = int(fit_order)
            for idx, i in enumerate(poly_order):
                if int(i) == fit_order:
                    break
            fit_order_idx = idx

        self.vol_fit = [i[fit_order_idx] for i in func]
        fit_rs = [i[fit_order_idx] for i in rs]

        # Expansion coefficients
        self.alpha_v = np.zeros([len(self.pressure), len(self.temperature)])
        for idx_p, v_p in enumerate(self.volume):
            vmin = v_p[idx_tmin]
            dv = v_p - vmin
            self.alpha_v[idx_p, :] = \
                self.vol_fit[idx_p].deriv(1)(dt) / (self.vol_fit[idx_p](dt) + vmin)

        self.alpha_v[:, idx_tmin] = 0. # Lowest temperature, alpha = 0

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# THERMAL EXPANSION COEFFICIENTS')
            file.write('%s\n\n' % '  To get thermal expansion coefficients, equilibrium volumes are fit as polynomial function of temperature at constant pressure.')
            file.write('%s%i\n' % ('## OPTIMAL ORDER OF POLYNOMIAL: ', fit_order))
            for idx_p, p in enumerate(self.pressure):
                file.write('%s%6.2f%s%6.4f\n\n' %
                           ('## EXPANSIONS AT ', p, '  GPa, R^2 = ', fit_rs[idx_p]))
                file.write('%10s%20s%20s\n' %
                           ('T(K)', 'Vol(Angstrom^3)', 'alpha_V(K^-1)'))
                vmin = self.volume[idx_p, idx_tmin]
                for idx_t, t in enumerate(self.temperature):
                    file.write('%10.1f%20.4f%20.8e\n' %
                               (t, vmin + self.vol_fit[idx_p](t - tmin), self.alpha_v[idx_p, idx_t]))
                file.write('\n')

            file.write('\n')
            file.close()

        return self

    def bulk_modulus(self, adiabatic=True, **kwargs):
        """
        Calculate isothermal and adiabatic bulk moduli at equilibrium volumes.

        The following equations are used:

        .. math::

            K_{T}(p;T) = V(p;T)\\left(\\frac{\\partial^{2}F(V;T)}{\\partial V^{2}}\\right)_{T}

            K_{S} = K_{T} + \\frac{\\alpha^{2}_{V}VTK^{2}_{T}}{C_{V}}

        To get :math:`K_{T}`, Helmholtz free energy is fit as isothermal EOSs.
        For ``self.thermo_eos()``, that means doing nothing; For
        ``self.thermo_freq()``, EOS fitting is required, whose form is the same
        as EOS used for :math:`E_{0}`.

        Args:
            adiabatic (bool): Whether to fit adiabatic bulk modulus. Thermal
                expansion coefficient needed.
            order, min_ndata_factor, max_poly_order_factor, min_poly_order_factor
                (int, optional): To restore EOS.

        Returns:
            self (Quasi_harmonic)

        **New attributes**

        * ``self.k_t`` nPressure\*nTemperature array, isothermal bulk modulus.  
        * ``self.k_s`` nPressure\*nTemperature array, adiabatic bulk modulus.  
        """
        import scipy.constants as scst
        from sympy import diff, lambdify, symbols
        import copy
        import numpy as np

        if adiabatic == True and not hasattr(self, 'alpha_v'):
            raise AttributeError('Expansion coefficient should be fit at first.')

        # Fit EOS
        if not hasattr(self, 'fe_eos'): # thermo_freq
            self.fe_eos = []
            self.fe_eos_method = self.e0_eos_method
            for idx_t, t in enumerate(self.temperature):
                eos_command = 'self.eos_fit(self.volume[:, idx_t], self.helmholtz[:, idx_t], self.e0_eos_method, write_out=False'
                # Polynomial / Deltafactor / Numerical
                for idx, key in enumerate(kwargs.keys()):
                    value = list(kwargs.values())[idx]
                    eos_command += ', {}={}'.format(key, value)
                eos_command += ')'
                eos, _ = eval(eos_command)
                self.fe_eos.append(eos)

        # Get K_T
        self.k_t = np.zeros(self.volume.shape)
        v = symbols('v')
        for idx_t, eos in enumerate(self.fe_eos):
            df = diff(eos(v), v, 2)
            lam_df = lambdify(v, df, 'numpy')
            self.k_t[:, idx_t] = self.volume[:, idx_t] * lam_df(self.volume[:, idx_t]) * 1e24 / scst.Avogadro
        # Get K_S
        if adiabatic == True:
            self.specific_heat()
            for idx_t, t in enumerate(self.temperature):
                if t > 1e-4: #0K
                    self.k_s[:, idx_t] = self.k_t[:, idx_t] + \
                        self.alpha_v[:, idx_t]**2 * self.volume[:, idx_t] * t * self.k_t[:, idx_t]**2 * 1e-21 * scst.Avogadro\
                        / self.c_v[:, idx_t]
                else:
                    self.k_s[:, idx_t] = 0.

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA BULK MODULI')
            file.write('%s\n\n' % '  Isothermal and adiabatic bulk moduli.')
            for idx_p, p in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' % ('## BULK MODULI K_T and K_S AT ', p, '  GPa'))
                if adiabatic == True:
                    file.write('%10s%20s%20s\n' % ('T(K)', 'K_T(GPa)', 'K_S(GPa)'))
                    for idx_t, t in enumerate(self.temperature):
                        file.write('%10.1f%20.8e%20.8e\n' %
                                   (t, self.k_t[idx_p, idx_t], self.k_s[idx_p, idx_t]))
                    file.write('\n')
                else:
                    file.write('%10s%20s\n' % ('T(K)', 'K_T(GPa)'))
                    for idx_t, t in enumerate(self.temperature):
                        file.write('%10.1f%20.8e\n' % (t, self.k_t[idx_p, idx_t]))
                    file.write('\n')

            file.write('\n')
            file.close()

        return self

    def specific_heat(self):
        """
        Calculate constant volume or pressure specific heat at equilibrium
        volumes.

        The following equation is used:

        .. math::

            C_{p} - C_{V} = \\alpha_{V}^{2}K_{T}VT

        Returns:
            self (Quasi_harmonic)

        **New attributes**

        * ``self.c_v`` or ``self.c_p``  Dependents on the fitting method.

        .. note::

            This method fits ``self.c_p`` by ``self.c_v`` when ``thermo_freq``
            and ``thermo_gruneisen`` was used. ``self.c_v`` is obtained by when
            ``thermo_eos`` is used.
        """
        import numpy as np
        import warnings
        import scipy.constants as scst

        if not hasattr(self, 'alpha_v') or not hasattr(self, 'k_t'):
            raise AttributeError(
                'Expansion coefficient and bulk modulus should be fit at first.')

        if not hasattr(self, 'c_p'): # thermo_freq
            self.c_p = self.c_v + self.alpha_v**2 * self.k_t * self.volume * self.temperature * 1e-21 * scst.Avogadro
        elif not hasattr(self, 'c_v'): # thermo_eos
            self.c_v = self.c_p - self.alpha_v**2 * self.k_t * self.volume * self.temperature * 1e-21 * scst.Avogadro
        else:
            warnings.warn("Attributes 'c_v' and 'c_p' both exist. Nothing is updated.")
            return self

        # Print output file
        if self.write_out:
            file = open(self.filename, 'a+')
            file.write('%s\n' % '# QHA SPECIFIC HEAT')
            file.write('%s\n\n' % '  Constant volume and pressure specific heat.')
            for idx_p, p in enumerate(self.pressure):
                file.write('%s%6.2f%s\n\n' % ('## SPECIFIC HEAT C_V and C_P AT ', p, '  GPa'))
                file.write('%10s%20s%20s\n' % ('T(K)', 'C_v(J/mol*K)', 'C_p(J/mol*K)'))
                for idx_t, t in enumerate(self.temperature):
                    file.write('%10.1f%20.8e%20.8e\n' %
                               (t, self.c_v[idx_p, idx_t], self.c_p[idx_p, idx_t]))
                file.write('\n')

            file.write('\n')
            file.close()

        return self


def _restore_pcel(crysout, scelphono):
    """
    Restore the primitive geometry expanded by 'SCELPHONO' and generate the
    Pymatgen Structure of the cell used for phonon calculation.

    Args:
        crysout (Crystal_output): :code:`CRYSTALpytools.io.Crystal_output` object.
        scellphono (list[int] | array[int]): ndimension\*ndimension or 3\*3
            matrix corresponds to the 'SCELPHONO' keyword.

    Returns:
        structures (list[Structure]): A list of Pymatgen Structure objects.
            nCalc\*1. For HA phonons and dispersions, nCalc=1. For QHA,
            nCalc=sampled HA points.
    """
    from pymatgen.core.structure import Structure, Molecule
    from pymatgen.core.lattice import Lattice
    import numpy as np
    import re
    import warnings

    ndimen = crysout.get_dimensionality()
    pbc = {3 : (True, True, True),
           2 : (True, True, False),
           1 : (True, False, False)}
    # Get structure. Address the issue with QHA file
    idx_line = 0
    structures = []
    while idx_line < len(crysout.data):
        if re.match(r'^\s+DIRECT LATTICE VECTORS CARTESIAN COMPONENTS',
                    crysout.data[idx_line]):
            idx_line += 2
            vec1 = np.array(crysout.data[idx_line].strip().split()[0:3], dtype=float)
            vec2 = np.array(crysout.data[idx_line + 1].strip().split()[0:3], dtype=float)
            vec3 = np.array(crysout.data[idx_line + 2].strip().split()[0:3], dtype=float)

            idx_line += 9
            all_species = []
            all_coord = []
            while re.match(r'^\s+[0-9]+\s+[0-9]+\s+[A-Z]+', crysout.data[idx_line]):
                data = crysout.data[idx_line].strip().split()
                all_coord.append(data[3:])
                all_species.append(data[2].capitalize())
                idx_line += 1
            all_coord = np.array(all_coord, dtype=float)
            scel_mx = np.vstack([vec1, vec2, vec3])

            # Molecule 0D
            if ndimen == 0:
                warnings.warn('0D system is used. There is nothing to reduce.')
                structures.append(Molecule(species=all_species, coords=all_coord))
                idx_line += 1
                continue

            # Periodic systems
            if scelphono != []:
                scell_mx = np.eye(3, dtype=float)
                scell_mx[: ndimen, : ndimen] = np.array(scelphono)[: ndimen, : ndimen]
                shrink_mx = np.linalg.pinv(scell_mx)
                pcel_mx = np.dot(scel_latt, shrink_mx)
                pcel_latt = Lattice(pcel_mx, pbc=pbc[ndimen])
                all_coord = np.dot(all_coord, np.linalg.pinv(pcel_mx)).tolist()
                pcel_coord = []
                pcel_species = []
                for i, coord in enumerate(all_coord):
                    if any(x > 0.5 or x <= -0.5 for x in coord[0:ndimen]):
                        continue
                    else:
                        pcel_coord.append(coord)
                        pcel_species.append(all_species[i])
            else:
                pcel_latt = Lattice(scel_mx, pbc=pbc[ndimen])
                pcel_coord = all_coord
                pcel_species = all_species

            struc = Structure(lattice=pcel_latt, species=pcel_species,
                              coords=pcel_coord, coords_are_cartesian=False)
            structures.append(struc)
            idx_line += 1
        else:
            idx_line += 1

    if structures == []:
        raise Exception('Valid structure not found.')
    elif len(structures) > 1: # QHA / HA + PREOPTGEOM, the first entry is pre-optimized geometry
        structures = structures[1:]

    return structures