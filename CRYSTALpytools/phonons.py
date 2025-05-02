#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A post-processing module for visualizing and analyzing phonon-related properties.
"""
from warnings import warn, filterwarnings

import numpy as np
from scipy import constants

from CRYSTALpytools import units


class PhononBand():
    """
    Phonon band object. Frequency unit: THz.

    .. note::

        Even though this class is not directly inherited from the
        ``electronics.ElectronBand`` class, its ``bands`` attribute still has
        the same, nBand\*nKpoint\*nSpin dimentionalities for using the shared
        plotting functions. nSpin is always 1 here.

    Args:
        tick_pos (array): 1\*nTick array of 1D tick coordinates. Unit: Angstrom
        tick_label (list): 1\*nTick of default tick labels
        bands (array): nBand\*nKpoint\*1 array of frequency. Unit: THz
        k_path (array): 1D coordinates of k points. Unit: Angstrom
        geometry (Structure): Pymatgen structure
        reciprocal_latt (array): 3\*3 array of reciprocal lattice matrix. Not
            valid if ``geometry`` is specified.
        tick_pos3d (array): 1\*nTick 3D fractional tick coordinates
        k_path3d (array): nKpoints\*3 3D fractional coordinates of k points
        unit (str): In principle, should always be 'THz': THz-Angstrom.
    """
    def __init__(self, tick_pos, tick_label, bands, k_path, geometry=None,
                 reciprocal_latt=None, tick_pos3d=None, k_path3d=None, unit='THz'):
        self.n_tick = len(tick_pos)
        self.tick_pos = np.array(tick_pos, dtype=float)
        self.tick_label = tick_label
        self.n_bands = len(bands)
        self.bands = np.array(bands, dtype=float)
        self.n_kpoints = len(k_path)
        self.k_path = np.array(k_path, dtype=float)
        self.geometry = geometry
        if np.all(self.geometry!=None):
            self.reciprocal_latt = self.geometry.lattice.reciprocal_lattice.matrix
        else:
            self.reciprocal_latt = reciprocal_latt
        self.tick_pos3d = np.array(tick_pos3d, dtype=float)
        self.k_path3d = np.array(k_path3d, dtype=float)
        self.unit = unit

    @classmethod
    def from_file(cls, output, q_overlap_tol=1e-4):
        """
        Generate an ``PhononBand`` object from the output file of CRYSTAL.

        .. note::

            Currently only the screen output ('.out') file is supported.

        Args:
            output (str): CRYSTAL output file
            q_overlap_tol (float): The threshold for overlapped k points. Only
                used for getting tick positions.
        Returns:
            cls (PhononBand)
        """
        from CRYSTALpytools.crystal_io import Crystal_output

        out = Crystal_output(output).get_phonon_band(q_overlap_tol=q_overlap_tol)
        return out

    def plot(self, **kwargs):
        """
        A wrapper to plot band structure of a single system using matplotlib.
        For input arguments or plotting multiple systems, check
        :ref:`plot.plot_phonon_bands() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``PhononBand`` object). Check documents for
                :ref:`plot.plot_electron_bands() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
        """
        from CRYSTALpytools.plot import plot_phonon_bands

        kwargs['mode'] = 'single'
        fig = plot_phonon_bands(self, **kwargs)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``PhononBand`` object. Internal method.

        Args:
            unit (str): 'THz' or 'cm-1'. Length unit is always in :math:`\\AA`.
        """
        from CRYSTALpytools.units import cm_to_thz, thz_to_cm

        if unit.lower() == self.unit.lower():
            return self

        opt_props = []  # Optional frequency properties
        if unit.lower() == 'thz':
            self.unit = 'THz'
            self.bands = cm_to_thz(self.bands)
            for p in opt_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, cm_to_thz(attrv))
        elif unit.lower() == 'cm-1':
            self.unit = 'cm-1'
            self.bands = thz_to_cm(self.bands)
            for p in opt_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, thz_to_cm(attrv))
        else:
            raise ValueError('Unknown unit.')
        return self


class PhononDOS():
    """
    Phonon DOS object. Frequency unit: THz.

    .. note::

        Even though this class is not directly inherited from the
        ``electronics.ElectronDOS`` class, its ``doss`` attribute still has the
        same, nProj\*nFrequency\*nSpin dimentionalities for using the shared
        plotting functions. nSpin is always 1 here.

    Args:
        doss (array): nProj\*nFrequency\*1 array of DOS.
        frequency (array): Positions of DOS peaks (x axis)
        unit (str): In principle, should always be 'THz': THz-Angstrom.
    """

    def __init__(self, doss, frequency, unit='THz'):
        self.n_proj = np.shape(doss)[0]
        self.doss = np.array(doss, dtype=float)
        self.n_frequency = len(frequency)
        self.frequency = np.array(frequency, dtype=float)
        self.unit = unit

    @classmethod
    def from_file(cls, output, read_INS=False, atom_prj=[], element_prj=[]):
        """
        Generate an ``PhononDOS`` object from the output file of CRYSTAL.

        .. note::

            Currently only the screen output ('.out') file is supported.

        Args:
            output (str): CRYSTAL output file
            read_INS (bool): Read the inelastic neutron scattering spectra.
            atom_prj (list): Read the projections of atoms with specified labels.
            element_prj (list): Read projections of elements with specified
                conventional atomic numbers.

        Returns:
            cls (PhononDOS)
        """
        from CRYSTALpytools.crystal_io import Crystal_output

        out = Crystal_output(output).get_phonon_dos(read_INS=read_INS,
                                                    atom_prj=atom_prj,
                                                    element_prj=element_prj)
        return out

    def plot(self, **kwargs):
        """
        A wrapper to plot density of states of a single system with matplotlib.
        For input arguments or plotting multiple systems, check
        :ref:`plot.plot_phonon_doss() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``PhononDOS`` object). Check documents for :ref:`plot.plot_phonon_doss() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
        """
        from CRYSTALpytools.plot import plot_phonon_doss

        fig = plot_phonon_doss(self, **kwargs)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``PhononDOS`` object.

        Args:
            unit (str): 'THz' or 'cm-1'. Length unit is always in :math:`\\AA`.
        """
        from CRYSTALpytools.units import cm_to_thz, thz_to_cm

        if unit.lower() == self.unit.lower():
            return self

        opt_f_props = [] # Optional frequency properties
        opt_d_props = [] # Optional density properties
        if unit.lower() == 'thz':
            self.unit = 'THz'
            self.frequency = cm_to_thz(self.frequency)
            self.doss = thz_to_cm(self.doss)
            for p in opt_f_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, cm_to_thz(attrv))
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, thz_to_cm(attrv))
        elif unit.lower() == 'cm-1':
            self.unit = 'cm-1'
            self.frequency = thz_to_cm(self.frequency)
            self.doss = cm_to_thz(self.doss)
            for p in opt_f_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, thz_to_cm(attrv))
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, cm_to_thz(attrv))
        else:
            raise ValueError('Unknown unit.')
        return self


class PhononBandDOS():
    """
    Phonon band + dos object. Frequency unit: THz.

    Args:
        band (PhononBand): ``PhononBand`` object
        dos (PhononDOS): ``PhononDOS`` object
    """

    def __init__(self, band, dos):
        self.band = band
        self.dos = dos

    @classmethod
    def from_file(cls, *output, q_overlap_tol=1e-4,
                  read_INS=False, atom_prj=[], element_prj=[]):
        """
        Get PhononBandDOS object from files

        Args:
            *output (str): CRYSTAL screen output file. 2 files, the first one
                is for band the second one is for DOS. Or a single output file
                file with both band and DOS.
            q_overlap_tol (float): The threshold for overlapped k points. Only
                used for getting tick positions.
            read_INS (bool): Read the inelastic neutron scattering spectra.
            atom_prj (list): Read the projections of atoms with specified labels.
            element_prj (list): Read projections of elements with specified
                conventional atomic numbers.
        Returns:
            cls (PhononBandDOS)
        """
        if len(output)==1:
            return cls(PhononBand.from_file(output[0], q_overlap_tol),
                       PhononDOS.from_file(output[0], read_INS, atom_prj, element_prj))
        elif len(output)==2:
            return cls(PhononBand.from_file(output[0], q_overlap_tol),
                       PhononDOS.from_file(output[1], read_INS, atom_prj, element_prj))
        else:
            raise ValueError('Only 1 or 2 entries are permitted.')

    def plot(self, **kwargs):
        """
        A wrapper to plot phonon band structure + density of states of a
        single system with matplotlib. For input arguments, check
        :ref:`plot.plot_phonon_banddos() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``PhononBandDOS`` object). Check documents for
                :ref:`plot.plot_phonon_banddos() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
        """
        from CRYSTALpytools.plot import plot_phonon_banddos

        fig = plot_phonon_banddos(self, **kwargs)
        return fig

    def _set_unit(unit):
        """
        Set units of data of ``PhononBandDOS`` object.

        Args:
            unit (str): 'THz' or 'cm-1'. Length unit is always in :math:`\\AA`.
        """
        self.band._set_unit(unit)
        self.dos._set_unit(unit)


class Phonon():
    """
    The basic class to save, analyze and visualize phonon vibrations. Unit:
    THz, eV and :math:`\\AA`. **Other units are not allowed, except outputs**.

    Args:
        structure (CStructure): Extended Pymatgen structure class. The extended
            structure for phonon dispersion calculations by finite displacement
            method must be reduced accordingly.
        u_0 (float): Internal energy. Unit: eV.
        qpoint (array): nQpoint\*4 array. The first 3 elements are fractional
            coordinates in reciprocal space and the fourth is its weight.
        frequency (array): nQpoint\*nMode array. Phonon frequencies in THz.
        mode_symm (array): nQpoint\*nMode array of irreducible representations.
            In Mulliken symbols.
        eigenvector (array): nQpoint\*nMode\*nAtom\*3 complex array. **Phased
            and mass-weighted** eigenvectors normalized to 1. Unit: :math:`\\AA`.
        fragment_idx (array): Indices of atoms (starting from 1) if a fragement
            calculation is performed.
    Returns:
        self (Phonon): Attribute names same as inputs
    """
    def __init__(self, structure, u_0, qpoint, frequency,
                 mode_symm=[], eigenvector=[], fragment_idx=[]):
        self.structure = structure
        self.u_0 = u_0
        self.qpoint = np.array(qpoint, dtype=float, ndmin=2)
        self.nqpoint = self.qpoint.shape[0]
        self.frequency = np.array(frequency, dtype=float, ndmin=2)
        self.nmode = self.frequency.shape[1]
        if self.frequency.shape[0] != self.nqpoint:
            raise Exception("Inconsistent q points between q coordinate and frequency.")

        self.mode_symm = np.array(mode_symm, dtype=str, ndmin=2)
        if self.mode_symm.shape[-1] > 0:
            if self.nqpoint != self.mode_symm.shape[0]:
                raise Exception("Inconsistent q points between q coordinate and symmetry.")
            if self.nmode != self.mode_symm.shape[1]:
                raise Exception("Inconsistent number of modes between frequency and symmetry.")
        else:
            self.mode_symm = np.array([['' for i in range(self.nmode)] for j in range(self.nqpoint)], dtype=str)

        self.eigenvector = np.array(eigenvector, dtype=complex, ndmin=4)
        if self.eigenvector.shape[-1] > 0:
            if self.nqpoint != self.eigenvector.shape[0]:
                raise Exception("Inconsistent q points between q coordinate and eigenvectors.")
            if self.nmode != self.eigenvector.shape[1]:
                raise Exception("Inconsistent number of modes between frequency and eigenvectors.")
            if self.nmode != self.eigenvector.shape[2]*3:
                raise Exception("Inconsistent number of atoms and number of modes. nAtom*3 must equal to nMode.")
        else:
            self.eigenvector = np.array(
                [[[[] for i in range(int(self.nmode/3))] for j in range(self.nmode)] for k in range(self.nqpoint)]
            )

        if len(fragment_idx) > 0:
            self._isfrag = True
            self._fragidx = np.array(fragment_idx, dtype=int) - 1
            self.natom = self._fragidx.shape[0]
            if self.natom != self.eigenvector.shape[2]:
                raise Exception("Inconsistent number of atoms and eigenvectors.")
        else:
            self._isfrag = False
            self.natom = self.structure.num_sites
            self._fragidx = np.array([i for i in range(self.natom)], dtype=int)
            if self.natom != self.eigenvector.shape[2]:
                raise Exception("Inconsistent number of atoms and eigenvectors.")

    @classmethod
    def from_file(cls, filename, source='crystal', q_id=None, q_coord=None,
                  **kwargs):
        """
        Instantiation from a file. Files generated by the following methods are
        supported (also the accepted entries for ``source``):

        * ``'crystal'``: CRYSTAL harmonic phonon outputs (Gamma, dispersion and band).  
        * ``'crystal-QHA'``: CRYSTAL QHA outputs.  
        * ``'phonopy'``: Phonopy band, qpoints and mesh files.

        .. note::

            In rare cases, the user might want to collect data of a few
            specific q points. ``q_id`` and ``q_coord`` can be set but not
            simultaneously. If set, ``q_id`` takes priority and ``q_coord`` is
            ignored.

        Args:
            filename (str): Name of the output file.
            source (str): The program used for the output file.
            q_id (array): Index (from 0) of q points to be read. 1\*nqpoint.
            q_coord (array): Fractional coordinates of q points to read.
                nqpoint\*3 list.
            \*\*kwargs: Other keywords passed to corresponding I/O functions.
                See below.
            qha_index (int): *crystal-QHA*. The index of calculation to read (from 0).
            struc_yaml (str): *phonopy*. 'qpoints.yaml' only.
            u_0 (float): *phonopy*. Static internal energy in eV.
            irreps_yaml (array[str]): *phonopy*. Read 'irreps.yaml' for mode 
                symmetry symbols. The fractional coordinates in the file are
                used for ``q_coord`` if not specified. Phonopy 'irreps.yaml'
                only contains symmetry info of a q point. Use a list of
                filenames for multiple q points.
        """
        from CRYSTALpytools.crystal_io import Crystal_output
        from CRYSTALpytools.io.phonopy import YAML

        source = source.lower()
        source_list = ['crystal', 'crystal-qha', 'phonopy']
        if source not in source_list:
            raise Exception("Unknown source file: '{}'.".format(source))
        if source == 'crystal':
            out = Crystal_output(filename).get_phonon(imaginary_tol=None,
                                                      q_overlap_tol=None,
                                                      eigvec='original')
            obj = out.phonon
        elif source == 'crystal-qha':
            if 'qha_index' not in kwargs.keys():
                raise Exception("The index of HA calculation in the output must be specified with 'qha_index' keyword.")
            qha_index = kwargs['qha_index']
            out = Crystal_output(filename).get_phonon(imaginary_tol=None,
                                                      q_overlap_tol=None,
                                                      eigvec='original')
            if out.phonon.nqpoint <= qha_index:
                raise Exception("The specified index '{:d}' does not exist in a QHA calculation of {:d} steps.".format(qha_index, out.phonon.nqpoint))
            # Get structure
            qhatitle = out.df[out.df[0].str.contains(r'\s+QUASI\-HARMONIC APPROXIMATION\s*$')].index
            if len(qhatitle) == 0:
                raise Exception("Not a CRYSTAL QHA calculation file.")
            dftitle = out.df[out.df[0].str.contains(r'\s*\*\s+CELL DEFORMATION\s*$')].index.tolist()
            dftitle.append(out.eoo)
            structitle = out.df[out.df[0].str.contains(r'^\s*ATOMS IN THE ASYMMETRIC UNIT')].index
            idx_struc = np.where(structitle < dftitle[qha_index+1])[0][-1]
            struc = out.get_geometry(initial=idx_struc)
            obj = cls(struc,
                      out.phonon.u_0[qha_index],
                      np.array([[0., 0., 0., 1.]]),
                      out.phonon.frequency[qha_index],
                      mode_symm=out.phonon.mode_symm[qha_index],
                      eigenvector=out.phonon.eigenvector[qha_index])
        elif source == 'phonopy':
            if 'u_0' not in kwargs.keys():
                warn("Static internal energy not available. Setting it to 0.", stacklevel=2)
                u_0 = 0.
            else:
                u_0 = kwargs['u_0']

            if 'struc_yaml' in kwargs.keys():
                strucfile = kwargs['struc_yaml']
            else:
                strucfile = filename
            pobj = YAML.read(strucfile, phonon=filename)

            if 'irreps_yaml' in kwargs.keys():
                symmfile = np.array(kwargs['irreps_yaml'], dtype=str, ndmin=1)
                # find the common q points in phonon and symmetry files
                qpts = []; mode_symm = []
                filterwarnings("ignore", "Unknown length unit. 'angstrom' is assumed.")
                for symm in symmfile:
                    sobj = YAML.read(strucfile, phonon=symm)
                    dist = np.linalg.norm(np.subtract(pobj.qpoint, sobj.qpoint[0]), axis=1)
                    idx = np.where(dist<1e-4)[0]
                    if len(idx) < 1:
                        warn("Q point coordinate [{:.2f} {:.2f} {:.2f}] defined in {} not found in {} and is skipped.".format(
                            q[0], q[1], q[2], symm, filename))
                    for i in idx:
                        qpts.append(i)
                        mode_symm.append(sobj.mode_symm[i])
                del sobj
            else:
                qpts = [i for i in range(pobj.nqpoint)]
                mode_symm = pobj.mode_symm

            obj = cls(pobj.structure, u_0, pobj.qpoint[qpts], pobj.frequency[qpts],
                      mode_symm=mode_symm, eigenvector=pobj.eigenvector[qpts])
            del pobj

        # specific q points
        if np.all(q_id==None) and np.all(q_coord==None):
            pass
        elif np.all(q_id!=None):
            qinfo = np.array(q_id, dtype=int)
            obj = cls(obj.structure, obj.u_0, obj.qpoint[q_info],
                      obj.frequency[q_info], mode_symm=obj.mode_symm[q_info],
                      eigenvector=obj.eigenvector[q_info])
        elif np.all(q_id==None) and np.all(q_coord!=None):
            q_coord = np.array(q_coord, dtype=float)
            q_info = []
            for q in q_coord:
                dist = np.linalg.norm(np.subtract(obj.qpoint[:, 0:3], q), axis=1)
                idx = np.where(dist<1e-4)[0]
                if len(idx) < 1:
                    warn("Q point coordinate [{:.2f} {:.2f} {:.2f}] not found and is skipped.".format(q[0], q[1], q[2]))
                    continue
                q_info += idx

            q_info = np.array(q_info, dtype=int)
            obj = cls(obj.structure, obj.u_0, obj.qpoint[q_info],
                      obj.frequency[q_info], mode_symm=obj.mode_symm[q_info],
                      eigenvector=obj.eigenvector[q_info])
        return obj

    # @classmethod
    # def from_GammaHessian(cls, structure, u_0, Hessian):
    #     """
    #     Instantiation from a Hessian matrix, assumed to be at :math:`\\Gamma`.

    #     Args:
    #         structure (CStructure): Extended Pymatgen structure class. The extended
    #             structure for phonon dispersion calculations by finite displacement
    #             method must be reduced accordingly.
    #         u_0 (float): Internal energy. Unit: eV.
    #         Hessian (array): 3nAtom\*3nAtom real Hessian matrix.
    #     Returns:
    #         cls (Vibrations): See references of the class.
    #     """
    # def get_Hessian(self, qpoint=[0., 0., 0.]):
    #     """
    #     Get Hessian on an **existing** q point.

    #     Args:
    #         qpoint (array): 1\*3 array of fractional coordinates in reciprocal space.
    #     Returns:
    #         hessmx (array): 3nAtom\*3nAtom complex dynamic matrix
    #     """


    # def displaced_structure(self, qpoint=[0., 0., 0.]):
    #     """
    #     Get the phased and mass-weighted dynamic matrix at the specified q point.

    #     Args:
    #         qpoint (array): 1\*3 array of fractional coordinates in reciprocal
    #             space.
    #     Returns:
    #         dynmx (array): 3nAtom\*3nAtom complex dynamic matrix
    #     """

    def clean_q_overlap(self, threshold=1e-4):
        """
        Remove the repeated q points and update the weight. Qpoint, frequency,
        symmetry and eigenvector are updated.

        Args:
            threshold (float): The q point overlap threshold. Defined in
                fractional coordinates.
        """
        if self.nqpoint <= 1: return self

        overlap = [[] for i in range(self.nqpoint)]
        for iq1 in range(1, self.nqpoint):
            vq1 = self.qpoint[iq1, 0:3]
            for iq2 in range(iq1): # q2 < q1, only attach q1 to minimum q2
                vq2 = self.qpoint[iq2, 0:3]
                if np.linalg.norm(vq1 - vq2) < threshold:
                    warn('Overlap of q points is detected between q points {:3d} and {:3d}'.format(iq1, iq2),
                         stacklevel=2)
                    overlap[iq2].append(iq1)
                    break

        saved = []; skipped = []
        for i, o in enumerate(overlap):
            if len(o) > 0:
                saved.append(i)
                skipped += o

        if len(saved) > 0:
            saved = np.array(saved, dtype=int)
            skipped = np.array(skipped, dtype=int)

            qpt = []; freq = []; symm = []; eigvt = []
            for iq in range(self.nqpoint):
                if iq in saved:
                    a_qpt = self.qpoint[iq]
                    for jq in overlap[np.where(saved==iq)[0]]:
                        a_qpt[3] += self.qpoint[jq, 3]
                    qpt.append(a_qpt)
                    freq.append(self.frequency[iq])
                    symm.append(self.mode_symm[iq])
                    eigvt.append(self.eigenvector[iq])
                elif iq in skipped:
                    continue
                else:
                    qpt.append(self.qpoint[iq])
                    freq.append(self.frequency[iq])
                    symm.append(self.mode_symm[iq])
                    eigvt.append(self.eigenvector[iq])

            # Update attributes
            self.nqpoint = len(qpt)
            self.qpoint = np.array(qpt, dtype=float)
            self.frequency = np.array(freq, dtype=float)
            self.mode_symm = np.array(symm, dtype=str)
            self.eigenvector = np.array(eigvt, dtype=float)
        return self

    def clean_imaginary(self, threshold=-1e-4):
        """
        Set negative frequenceies and related properteis to 0 and print warning
        message. Only frequency is modified.

        Args:
            threshold (float): The threshold to identify a phonon mode as negative.
        Returns:
            self
        """
        for q, freq in enumerate(self.frequency):
            neg_rank = np.where(freq < threshold)[0]
            if len(neg_rank) == 0:
                continue
            dist = np.round(np.linalg.norm(self.qpoint[q, 0:3]), 4)
            if dist > 1e-4 and len(neg_rank) > 0:
                warn('OFF-GAMMA IMAGINARY MODES! The structure is highly probable to be unstable.', stacklevel=2)
            if dist <= 1e-4 and len(neg_rank) > 3:
                warn('MORE THAN 3 IMAGINARY MODES AT GAMMA! The structure is highly probable to be unstable.', stacklevel=2)
            self.frequency[q, neg_rank] = 0.
        return self

    def scale_frequency(self, scale, scale_range=[]):
        """
        Empirically scaling phonon frequencies by the scale factor. It does not
        update the attribute.

        Args:
            scale (float|array): Scaling factors
            scale_range (array): nScale\*2 array. The frequency range of the
                scaling factor. Including the minimum and excluding the maximum.
                Use empty list for a uniform scaling.
        Returns:
            freq (array): nQpoint\*nMode frequencies in THz.
        """
        scale = np.array(scale, dtype=float, ndmin=1)
        scale_range = np.array(scale_range, dtype=float, ndmin=2)
        if scale_range.shape[-1] == 0:
            scale_range = np.array([[self.frequency.min(), self.frequency.max()+1.]], dtype=float)
        scale_range = np.round(scale_range, 10)
        freq = np.round(self.frequency, 10)
        if scale.shape[0] != scale_range.shape[0]:
            raise Exception("Inconsistent number of scalind factors and scaled frequency range.")

        for s, r in zip(scale, scale_range):
            idx = np.where((freq>=r.min())&(freq<r.max()))
            freq[idx] *= s
        return freq

    def classical_eigvec(self):
        """
        Return to phonon eigenvector with classical amplitude. The attribute is
        not changed.

        Returns:
            eigvec (array): Eigenvector with classical amplitude, in :math:`\\AA`
                and same dimension as ``self.eigenvector``.
        """
        from scipy import constants

        # Hartree units are used for consistency
        amu2me = units.amu_to_me(1.)
        ev2ha = units.eV_to_H(1.)
        ang2br = units.angstrom_to_au(1.)
        br2ang = units.au_to_angstrom(1.)
        thz2au = constants.h / constants.physical_constants['Hartree energy'][0] * 1e12

        allmass = [i.data['Atomic mass']*amu2me for i in self.structure.species]
        mass = np.zeros([1, self.natom*3])
        for i, oldi in enumerate(self._fragidx):
            nrow = int(i * 3)
            mass[0, nrow:nrow+3] = allmass[oldi]

        eigvec = self.unweight_eigvec()
        eigvec = eigvec.reshape([self.nqpoint, self.nmode, self.natom*3], order='C')
        for iq in range(self.nqpoint):
            meff = np.sum(np.multiply((eigvec[iq]*eigvec[iq].conj()).real, mass), axis=1)
            amplitude = np.zeros(meff.shape) + 1
            idx = np.where(self.frequency[iq]>1e-4)[0]
            amplitude[idx] = np.sqrt(1 / (meff[idx]*self.frequency[iq, idx]*thz2au))
            for im in range(self.nmode):
                eigvec[iq, im] = np.multiply(eigvec[iq, im], amplitude[im])
        eigvec *= br2ang
        eigvec = eigvec.reshape([self.nqpoint, self.nmode, self.natom, 3], order='C')
        return eigvec

    def unweight_eigvec(self):
        """
        Return to mass-unweighted phonon eigenvector normalized to 1. The
        attribute is not changed.

        Returns:
            eigvec (array): Mass-unweighted eigenvector normalized to 1.
        """
        allmass = [i.data['Atomic mass'] for i in self.structure.species]
        mass = np.zeros([1, self.natom*3])
        for i, oldi in enumerate(self._fragidx):
            nrow = int(i * 3)
            mass[0, nrow:nrow+3] = allmass[oldi]

        eigvec = np.zeros([self.nqpoint, self.nmode, self.nmode], dtype=complex)
        for iq in range(self.nqpoint):
            eigvec[iq] = self.eigenvector[iq].reshape([self.nmode, self.nmode], order='C') # nMode*3nAtom
            eigvec[iq] = np.multiply(eigvec[iq], 1/np.sqrt(mass)) # 1/\sqrt{m}
            for im in range(self.nmode):
                eigvec[iq, im] /= np.linalg.norm(eigvec[iq, im])

        eigvec = eigvec.reshape([self.nqpoint, self.nmode, self.natom, 3], order='C')
        return eigvec

    @classmethod
    def get_eigvec(cls, struc, freq, eigvec, method, fragment_idx=[]):
        """
        Get eigenvectors suitable for instantiation.

        Args:
            struc (CStructure): Structure
            freq (array): nQpoint\*nMode array of phonon frequencies in THz.
            eigvec (array): nQpoint\*nMode\*nAtom\*3 complex array of phonon
                eigenvectors, must corresponding to the method specified.
            method (str): 'remove classical', 'mass weight'.
            fragment_idx (array): 1\*nAtom array of ints. Atomic indices from 1.
        Returns:
            eigvec (array): Mass-weighted and phased eigenvector normalized to 1.
        """
        from scipy import constants

        method = method.lower()
        if method not in ['remove classical', 'mass weight']:
            raise Exception("Unknown method: '{}'.".format(method))

        if len(fragment_idx) > 0:
            natom = len(fragment_idx)
            fragment_idx = np.array(fragment_idx, dtype=int, ndmin=1) - 1
        else:
            natom = struc.num_sites
            fragment_idx = np.array([i for i in range(natom)], dtype=int, ndmin=1)

        eigvec = np.array(eigvec, dtype=complex, ndmin=4)
        freq = np.array(freq, dtype=float, ndmin=2)
        nmode = freq.shape[1]; nqpoint = freq.shape[0]
        if nmode != natom*3: raise Exception("Inconsistent number of modes and atoms.")

        allmass = [i.data['Atomic mass'] for i in struc.species]
        mass = np.zeros([1, natom*3])
        for i, oldi in enumerate(fragment_idx):
            nrow = int(i * 3)
            mass[0, nrow:nrow+3] = allmass[oldi]

        if method == 'remove classical':
            # Hartree units are used for amplitude and mass. No need to convert back due to normalization
            thz2au = constants.h / constants.physical_constants['Hartree energy'][0] * 1e12
            for iq in range(nqpoint):
                eigvt = eigvec[iq].reshape([nmode, nmode], order='C') # nMode*3nAtom
                eigvt = units.angstrom_to_au(eigvt)
                mass = units.amu_to_me(mass) # me
                meff = np.sum(np.multiply((eigvt*eigvt.conj()).real, mass), axis=1)
                amplitude = np.zeros(meff.shape) + 1
                idx = np.where(freq[iq]>1e-4)[0]
                amplitude[idx] = 1 / np.sqrt(1 / (meff[idx]*freq[iq, idx]*thz2au)) # 1/amplitude
                for im in range(nmode):
                    eigvt[im] = np.multiply(np.multiply(eigvt[im], amplitude[im]), np.sqrt(mass))
                    eigvt[im] /= np.linalg.norm(eigvt[im])
                eigvec[iq] = eigvt.reshape([nmode, natom, 3], order='C')
        elif method == 'mass weight':
            for iq in range(nqpoint):
                eigvt = eigvec[iq].reshape([nmode, nmode], order='C') # nMode*3nAtom
                eigvt = np.multiply(eigvt, np.sqrt(mass))
                for im in range(nmode):
                    eigvt[im] /= np.linalg.norm(eigvt[im])
                eigvec[iq] = eigvt.reshape([nmode, natom, 3], order='C')

        return eigvec
