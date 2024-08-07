#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A post-processing module for phonon properties
"""
from CRYSTALpytools import units
import numpy as np

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


