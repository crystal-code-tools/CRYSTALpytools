#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A post-processing module for electronic properties
"""
from CRYSTALpytools import units
import numpy as np

class ElectronBand():
    """
    Electron band object. Energy unit: eV. E Fermi is aligned to 0.

    Args:
        spin (int): 1, closed shell; 2, open shell
        tick_pos (array): 1\*nTick array of 1D tick coordinates. Unit: Angstrom
        tick_label (list): 1\*nTick of default tick labels
        efermi (float): Fermi energy. Unit: eV.
        bands (array): nBand\*nKpoint\*nSpin array of energy. Unit: eV
        k_path (array): 1D coordinates of k points. Unit: Angstrom
        geometry (Structure): Pymatgen structure
        reciprocal_latt (array): 3\*3 array of reciprocal lattice matrix. Not
            valid if ``geometry`` is specified.
        tick_pos3d (array): 1\*nTick 3D fractional tick coordinates
        k_path3d (array): nKpoints\*3 3D fractional coordinates of k points
        unit (str): In principle, should always be 'eV': eV-Angstrom.
    """

    def __init__(self, spin, tick_pos, tick_label, efermi, bands, k_path,
                 geometry=None, reciprocal_latt=None, tick_pos3d=None,
                 k_path3d=None, unit='eV'):
        import numpy as np

        self.spin = spin
        self.n_tick = len(tick_pos)
        self.tick_pos = np.array(tick_pos, dtype=float)
        self.tick_label = tick_label
        self.efermi = efermi
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
        # old attrs, commensurate with old codes
        self.tick_position = self.tick_pos
        self.k_point_plot = self.k_path
        self.k_point_pos3d = self.k_path3d

    @classmethod
    def from_file(cls, band, output=None):
        """
        Generate an ``ElectronBand`` object from fort.25 / BAND.DAT file.
        Optional 3D space data is read from the output file of 'properties'.

        Args:
            band (str): 'fort.25' or 'BAND.DAT'
            output (str): Properties output file
        Returns:
            cls (ElectronBand)
        """
        from CRYSTALpytools.base.extfmt import CrgraParser, XmgraceParser
        from CRYSTALpytools.crystal_io import Properties_output

        file = open(band)
        flag = file.readline()
        file.close()
        if '-%-' in flag:  # fort.25 file format
            bandout = CrgraParser.band(band)
        else:
            bandout = XmgraceParser.band(band)

        if np.all(output!=None):
            pout = Properties_output(output)
            struc = pout.get_geometry()
            t3d, k3d = pout.get_3dkcoord()
            return cls(spin=bandout[0], tick_pos=bandout[1],
                       tick_label=bandout[2], efermi=bandout[3],
                       bands=bandout[4], k_path=bandout[5], geometry=struc,
                       tick_pos3d=t3d, k_path3d=k3d, unit=bandout[6])
        else:
            return cls(spin=bandout[0], tick_pos=bandout[1],
                       tick_label=bandout[2], efermi=bandout[3],
                       bands=bandout[4], k_path=bandout[5], unit=bandout[6])

    def plot(self, **kwargs):
        """
        A wrapper to plot band structure of a single system using matplotlib.
        For input arguments or plotting multiple systems, check
        :ref:`plot.plot_electron_bands() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``ElectronBand`` object). Check documents for
                :ref:`plot.plot_electron_bands() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
            ax (Axes): Matplotlib axes object
        """
        from CRYSTALpytools.plot import plot_electron_bands

        kwargs['mode'] = 'single'
        fig, ax = plot_electron_bands(self, **kwargs)
        return fig, ax

    @property
    def bandgap(self):
        """
        A shortcut for band gap only.
        """
        return self.get_bandgap()[0]

    def get_bandgap(self):
        """
        Get band gap. For spin-polarized systems, 2\*1 arrays are used for
        :math:`\\alpha` and :math:`\\beta` states. Data is rounded to 6 decimal
        places.

        Returns:
            self.gap (float): Band gap. Default unit: eV
            self.vbm (flat): Valence band maximum, with reference to Fermi
                level. Default unit: eV
            self.cbm (float): Conduction band minimum, with reference to Fermi
                level. Default unit: eV
            self.gap_pos (array): 1D coordinates of vbm (1st element) and cbm
                (2nd element). For spin-polarized cases, ``self.gap_pos[0, :]``
                are vbm and cbm of :math:`\\alpha` state. Default unit: Angstrom
        """
        import numpy as np

        self.gap = np.zeros([2,], dtype=float)
        self.vbm = np.zeros([2,], dtype=float)
        self.cbm = np.zeros([2,], dtype=float)
        self.gap_pos = np.zeros([2, 2], dtype=float)
        for ispin in range(self.spin):
            for nbd, ebd in enumerate(self.bands[:, 0, ispin]):
                if ebd > 0:
                    nvb = nbd - 1
                    ncb = nbd
                    break
                else:
                    continue

            vbm = np.round(np.max(self.bands[nvb, :, ispin]), 6)
            kvbm = self.k_path[np.argmax(self.bands[nvb, :, ispin])]
            cbm = np.round(np.min(self.bands[ncb, :, ispin]), 6)
            kcbm = self.k_path[np.argmin(self.bands[ncb, :, ispin])]
            if vbm > 0. or cbm < 0.:
                gap = 0.
            else:
                gap = cbm - vbm

            self.gap[ispin] = gap
            self.vbm[ispin] = vbm
            self.cbm[ispin] = cbm
            self.gap_pos[ispin, :] = [kvbm, kcbm]

        if self.spin == 1:
            self.gap = self.gap[0]
            self.vbm = self.vbm[0]
            self.cbm = self.cbm[0]
            self.gap_pos = self.gap_pos[0]

        return self.gap, self.vbm, self.cbm, self.gap_pos

    def to_pmg(self, labels=None):
        """
        Get Pymatgen ``BandStructureSymmLine`` object (inherited from ``BandStructure``).
        No projection is available for now.

        .. note::
            3D information for output file is mandatory here.

        Args:
            labels (list[str]): K point labels to display in the band structure.
        Returns:
            BandStructureSymmLine: Pymatgen band structure.
        """
        import warnings

        import numpy as np
        from pymatgen.core.lattice import Lattice
        from pymatgen.electronic_structure.bandstructure import \
            BandStructureSymmLine
        from pymatgen.electronic_structure.core import Spin

        if not hasattr(self, 'tick_pos3d'):
            raise Exception(
                '3D information is unknown: No properties output file was read.')

        # Set unit to eV-Angstrom
        self._set_unit('eV')

        rep_latt = self.reciprocal_latt
        # label dictionary
        labels_dict = {}
        if np.all(labels==None):
            labels = self.tick_label
        else:
            if len(labels) < self.n_tick:
                warnings.warn(
                    '''{:d} ticks available in band object, but {:d} labels are provided.
The default labels will be used for missing ones.'''.format(self.n_tick, len(labels)),
                    stacklevel=2
                )
                for i in range(len(labels), self.n_tick):
                    labels.append(self.tick_label[i])

            elif len(labels) > self.n_tick:
                warnings.warn(
                    '''{:d} ticks available in band object, but {:d} labels are provided.
The redundant labels will be omitted.'''.format(self.n_tick, len(labels)),
                    stacklevel=2
                )
                labels = labels[:self.n_tick]

            else:
                pass

        for i in range(self.n_tick):
            labels_dict[labels[i]] = self.tick_pos3d[i]

        # Energy eigenvalues
        # pymatgen will plot the bands wrt to the Fermi Energy
        band_energy = self.bands + self.efermi
        if self.spin == 1:
            eigenvals = {Spin.up: band_energy[:, :, 0]}
        else:
            eigenvals = {Spin.up: band_energy[:, :, 0],
                         Spin.down: band_energy[:, :, 1]}

        return BandStructureSymmLine(kpoints=self.k_path3d,
                                     eigenvals=eigenvals,
                                     lattice=Lattice(self.reciprocal_latt),
                                     efermi=self.efermi,
                                     labels_dict=labels_dict,
                                     coords_are_cartesian=False)

    def _set_unit(self, unit):
        """
        Set units of data of ``ElectronBand`` object. Internal method.

        Args:
            unit (str): 'eV': Energy unit = eV, Length unit = Angstrom;
                'a.u.': Energy unit = Hartree. Length unit = Bohr
        """
        from CRYSTALpytools.units import (H_to_eV, angstrom_to_au,
                                          au_to_angstrom, eV_to_H)

        if unit.lower() == self.unit.lower():
            return self

        opt_e_props = ['gap', 'vbm', 'cbm']  # Optional energy properties
        opt_d_props = ['gap_pos']  # Optional distance properties
        if unit.lower() == 'ev':
            self.unit = 'eV'
            self.bands = H_to_eV(self.bands)
            self.efermi = H_to_eV(self.efermi)
            self.tick_pos = au_to_angstrom(self.tick_pos)
            self.k_path = au_to_angstrom(self.k_path)
            if np.all(self.reciprocal_latt!=None):
                self.reciprocal_latt = au_to_angstrom(self.reciprocal_latt)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, H_to_eV(attrv))
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, au_to_angstrom(attrv))
        elif unit.lower() == 'a.u.':
            self.unit = 'a.u.'
            self.bands = eV_to_H(self.bands)
            self.efermi = eV_to_H(self.efermi)
            self.tick_pos = angstrom_to_au(self.tick_pos)
            self.k_path = angstrom_to_au(self.k_path)
            if np.all(self.reciprocal_latt!=None):
                self.reciprocal_latt = angstrom_to_au(self.reciprocal_latt)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, angstrom_to_au(attrv))
        else:
            raise ValueError('Unknown unit.')

        return self


class ElectronDOS():
    """
    Electron DOS object. Energy unit: eV. E Fermi is aligned to 0.

    Args:
        spin (int): 1, closed shell; 2, open shell
        efermi (float): Fermi energy. Unit: eV.
        doss (array): n_proj\*n_energy\*spin array of DOS. Positive values
            for both spin up and spin down states
        energy (array): Positions of DOS peaks (x axis)
        unit (str): In principle, should always be 'eV': eV-Angstrom.
    """

    def __init__(self, spin, efermi, doss, energy, unit='eV'):
        import numpy as np

        self.spin = spin
        self.efermi = efermi
        self.n_proj = np.shape(doss)[0]
        self.doss = np.array(doss, dtype=float)
        self.n_energy = len(energy)
        self.energy = np.array(energy, dtype=float)
        self.unit = unit

    @classmethod
    def from_file(cls, dos):
        """
        Generate an ``ElectronDOS`` object from fort.25 / DOSS.DAT file.

        Args:
            band (str): 'fort.25' or 'DOSS.DAT'
        Returns:
            cls (ElectronDOS)
        """
        from CRYSTALpytools.base.extfmt import CrgraParser, XmgraceParser

        file = open(dos)
        flag = file.readline()
        file.close()
        if '-%-' in flag:  # fort.25 file format
            dosout = CrgraParser.dos(dos)
        else:
            dosout = XmgraceParser.dos(dos)

        return cls(spin=dosout[0], efermi=dosout[1], doss=dosout[2],
                   energy=dosout[3], unit=dosout[4])

    def plot(self, **kwargs):
        """
        A wrapper to plot density of states of a single system with matplotlib.
        For input arguments or plotting multiple systems, check
        :ref:`plot.plot_electron_doss() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``ElectronDOS`` object). Check documents for
                :ref:`plot.plot_electron_doss() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
            ax (Axes): Matplotlib axes object
        """
        from CRYSTALpytools.plot import plot_electron_doss

        fig, ax = plot_electron_doss(self, **kwargs)
        return fig, ax

    def _set_unit(self, unit):
        """
        Set units of data of ``ElectronDOS`` object.

        Args:
            unit (str): 'eV': Energy unit = eV;
                'a.u.': Energy unit = Hartree
        """
        from CRYSTALpytools.units import H_to_eV, eV_to_H

        if unit.lower() == self.unit.lower():
            return self

        opt_e_props = []  # Optional energy properties
        opt_d_props = []  # Optional energy inverse properties
        if unit.lower() == 'ev':
            self.unit = 'eV'
            self.efermi = H_to_eV(self.efermi)
            self.energy = H_to_eV(self.energy)
            self.doss = eV_to_H(self.doss)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, H_to_eV(attrv))
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
        elif unit.lower() == 'a.u.':
            self.unit = 'a.u.'
            self.efermi = eV_to_H(self.efermi)
            self.energy = eV_to_H(self.energy)
            self.doss = H_to_eV(self.doss)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, H_to_eV(attrv))
        else:
            raise ValueError('Unknown unit.')

        return self


class ElectronBandDOS():
    """
    Electron band + dos object. Energy unit: eV. E Fermi is aligned to 0.

    Args:
        band (ElectronBand): ``ElectronBand`` object
        dos (ElectronDOS): ``ElectronDOS`` object
    """

    def __init__(self, band, dos):
        self.band = band
        self.dos = dos

    @classmethod
    def from_file(cls, *files, output=None):
        """
        Get ElectronBandDOS object from files

        Args:
            *files (str): 2 files, the first one is for band, 'fort.25' or
                'BAND.DAT'; the second one is for DOS, 'fort.25' or 'DOSS.DAT'.
                Or a single 'fort.25' file with both band and DOS.
            output (str): Property output file
        Returns:
            cls (ElectronBandDOS)
        """
        from CRYSTALpytools.electronics import ElectronBand, ElectronDOS

        if len(files)==1:
            return cls(ElectronBand.from_file(files[0], output),
                       ElectronDOS.from_file(files[0]))
        elif len(files)==2:
            return cls(ElectronBand.from_file(files[0], output),
                       ElectronDOS.from_file(files[1]))
        else:
            raise ValueError('Only 1 or 2 entries are permitted.')

    def plot(self, **kwargs):
        """
        A wrapper to plot electron band structure + density of states of a]
        single system with matplotlib. For input arguments, check
        :ref:`plot.plot_electron_banddos() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``ElectronBandDOS`` object). Check documents for
                :ref:`plot.plot_electron_banddos() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
            ax (Axes): Matplotlib axes object
        """
        from CRYSTALpytools.plot import plot_electron_banddos

        fig, ax = plot_electron_banddos(self, **kwargs)
        return fig, ax

    def _set_unit(unit):
        """
        Set units of data of ``ElectronBandDOS`` object.

        Args:
            unit (str): 'eV': Energy unit = eV, length unit = Angstrom
                'a.u.': Energy unit = Hartree, length unit = Bohr
        """
        self.band._set_unit(unit)
        self.dos._set_unit(unit)


class ChargeDensity():
    """
    Charge (spin) density object. Unit: :math:`e.\\AA^{-3}`. 3D plot under
    developing.

    Args:
        data (array): Plot data. nX\*nY\*nSpin or nX\*nY\*nZ\*nSpin
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC (2D) or 3 base vectors (3D)
        spin (int): 1 or 2.
        dimen (int): Dimensionality of the plot.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """

    def __init__(self, data, base, spin, dimen, struc=None, unit='Angstrom'):
        import numpy as np
        import warnings

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.spin = int(spin)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'ECHG' # Hidden for charge density plot. Useful for TOPOND child class

    @classmethod
    def from_file(cls, *files, method=None):
        """
        Generate a ``ChargeDensity`` object from a single file, or from multiple
        files by substracting values from the first entry. Can be used for
        multiple dimensions (2D only now. 3D under development.)

        Available methods are:

        * 'substact': Substracting data from the first entry based on following
            entries. Multiple entries only.  
        * 'alpha_beta': Save spin-polarized data in :math:`\\alpha` /
            :math:`\\beta` states, rather than charge(:math:`\\alpha+\\beta`)
            / spin(:math:`\\alpha-\\beta`). Single entry only.

        Args:
            \*files (str): Path to the charge density / spin density file(s).
                All the entries must be in the same file format.
        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        file = open(files[0], 'r')
        header = file.readline()
        file.close()
        if '-%-' in header: # 2D plot in fort.25
            cls = Properties_output().read_ECHG(*files, method=method)
        return cls

    def substract(self, *args):
        """
        Substracting data of the same type from the object.

        Args:
            \*args (str|ChargeDensity): File names or ``ChargeDensity`` objects.
                Must be of the same type (check the attribute ``type``).
        Returns:
            self (ChargeDensity) : spin dimension, if there is, is not kept.
                Only charge density difference is substracted.
        """
        from CRYSTALpytools.crystal_io import Properties_output
        import numpy as np

        for i in args:
            if isinstance(i, str):
                obj = Properties_output().read_ECHG(i, method=None)
            elif isinstance(i, ChargeDensity):
                obj = i
            else:
                raise TypeError('Inputs must be file name strings or ChargeDensity objects.')

            # base vector
            if not np.all(np.abs(self.base-obj.base)<1e-6):
                raise ValueError('Inconsistent base vectors between input and object.')
            # dimensionality
            if self.dimension != obj.dimension:
                raise ValueError('Inconsistent dimensionality between input and object.')
            # mesh grid
            for i in range(self.dimension):
                if self.data.shape[i] != obj.data.shape[i]:
                    raise ValueError('Inconsistent mesh grid between input and object.')
            # spin
            if self.spin != obj.spin:
                raise ValueError('Inconsistent spin dimensionalities between input and object.')
            # substract
            self.data = self.data - obj.data

        # spin difference is not kept - meaningless. (if not please remove this block)
        if self.spin == 2:
            oshape = self.data.shape
            chglen = 1
            for s in oshape[:-1]:
                chglen = chglen * s
            chglen = int(chglen)
            self.data = self.data.flatten(order='F')[:chglen]
            self.data = np.reshape(self.data, oshape[:-1])
            self.spin = 1
        return self

    def alpha_beta(self):
        """
        Get the :math:`\\alpha` / :math:`\\beta` state density, rather than
        charge(:math:`\\alpha+\\beta`) / spin(:math:`\\alpha-\\beta`).
        ``spin=2`` only.

        Returns:
            self (ChargeDensity) : The first entry of ``self.data`` is :math:`\\alpha`
                state density and the second is :math:`\\beta` state.
        """
        import numpy as np

        if self.spin != 2:
            raise ValueError('Not a spin-polarized system.')

        # can be used for any dimension
        oldshape = self.data.shape
        lenchg = 1
        for i in oldshape[:-1]:
            lenchg = lenchg * i
        lenchg = int(lenchg)
        alpha = (self.data.flatten(order='F')[:lenchg] + self.data.flatten(order='F')[lenchg:]) / 2
        beta = (self.data.flatten(order='F')[:lenchg] - self.data.flatten(order='F')[lenchg:]) / 2
        self.data = np.reshape(np.hstack([alpha, beta]), oldshape, order='F')
        return self

    def plot_2D(self, unit='Angstrom', option='both', levels=150, lineplot=False,
                linewidth=1.0, isovalues=None, colorplot=True, colormap='jet',
                cbar_label=None, a_range=[], b_range=[], rectangle=False, cellplot=False,
                x_ticks=5, y_ticks=5, add_title=True, figsize=[6.4, 4.8], **kwargs):
        """
        Plot 2D charge/spin density map. A wrapper of ``plot.plot_dens_ECHG``
        and ``plot.plot_spin_ECHG``.

        Available options:

        * 'both' : If spin polarized, plot both charge and spin densities.
            Otherwise plot charge densities.  
        * 'charge': Plot charge density.  
        * 'spin': Plot spin density.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr:math:`^{-3}`.
            option (str): Available options see above.
            levels (int|array): Set levels of contour plot. A number for
                linear scaled plot colors or an array for user-defined levels,
                **must be consistent with ``unit``**. 2\*nLevel can be defined
                when ``option='both'``.
            lineplot (bool): Plot contour lines.
            linewidth (float): Contour linewidth. Useful only if
                ``lineplot=True``. Other properties are not editable. Solid
                black lines for positive values and 0, dotted for negative.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 1\*2 list of colorbar titles can be set for
                spin-polarized systems. 'None' for default.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in
                fractional coordinate.
            rectangle (bool): If :math:`a, b` are non-orthogonal, plot a
                rectangle region and reset :math:`b`. If used together with
                ``b_range``, that refers to the old :math:`b`.
            cellplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range or rectangle
                options).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            add_title (bool): Whether to add property plotted as title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            \*\*kwargs : Other arguments passed to ``axes.contour()`` function
                to set contour lines.

        Returns:
            fig (Figure): Matplotlib figure object
            ax (Axes): Matplotlib axes object
        """
        from CRYSTALpytools.base.plotbase import plot_2Dscalar
        import numpy as np
        import matplotlib.pyplot as plt

        # dimen
        if self.dimension != 2:
            raise Exception('Not a 2D charge density object.')

        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)

        # levels
        if isinstance(levels, int) or isinstance(levels, float):
            if self.spin == 1:
                levels1 = np.linspace(np.min(self.data), np.max(self.data), int(levels))
                levels2 = levels1
            else:
                len_flat = 1
                for i in self.data.shape[:-1]:
                    len_flat = len_flat * i
                chg = self.data.flatten(order='F')[:len_flat]
                spin = self.data.flatten(order='F')[len_flat:]
                levels1 = np.linspace(np.min(chg), np.max(chg), levels)
                levels2 = np.linspace(np.min(spin), np.max(spin), levels)
                del chg, spin
        else:
            if isinstance(levels[0], int) or isinstance(levels[0], float):
                levels1 = np.array(levels, dtype=float)
                levels2 = np.array(levels, dtype=float)
            else:
                levels1 = np.array(levels[0], dtype=float)
                levels2 = np.array(levels[1], dtype=float)
        # color plot
        if colorplot == False:
            colormap = None
        # contour line
        if lineplot == True:
            chgline = [['k', '-', linewidth] for i in levels1]
            spinline = []
            for j in levels2:
                if j >= 0: spinline.append(['k', '-', linewidth])
                else: spinline.append(['k', 'dotted', linewidth])
        else:
            chgline = None; spinline = None
        # cbar_label
        if np.all(cbar_label==None):
            if unit.lower() == 'angstrom':
                cbar_label1 = r'Charge Density ($|e|/\AA^{3}$)'; cbar_label2 = r'Spin Density ($|e|/\AA^{3}$)'
            else:
                cbar_label1 = r'Charge Density ($|e|/Bohr^{3}$)'; cbar_label2 = r'Spin Density ($|e|/Bohr^{3}$)'
        else:
            if isinstance(cbar_label, list):
                cbar_label1 = cbar_label[0]; cbar_label2 = cbar_label[1]; 
            else:
                cbar_label1 = str(cbar_label); cbar_label2 = str(cbar_label)

        # plot
        ## spin
        if self.spin == 1 and (option.lower()=='both' or option.lower()=='spin'):
            warnings.warn("Spin options not available to non spin-polarized cases.",
                          stacklevel=2)
            option = 'charge'

        if option.lower() == 'both':
            fig, ax = plt.subplots(1, 2, figsize=figsize, sharex=True,
                                   sharey=True, layout='tight')
            fig, ax[0] = plot_2Dscalar(
                fig, ax[0], self.data[:, :, 0], self.base, levels1, chgline,
                isovalues, colormap, cbar_label1, a_range, b_range, rectangle,
                cellplot, x_ticks, y_ticks, **kwargs
            )
            fig, ax[1] = plot_2Dscalar(
                fig, ax[1], self.data[:, :, 1], self.base, levels2, spinline,
                isovalues, colormap, cbar_label2, a_range, b_range, rectangle,
                cellplot, x_ticks, y_ticks, **kwargs
            )
        elif option.lower() == 'charge':
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            if self.spin == 1:
                fig, ax = plot_2Dscalar(
                    fig, ax, self.data[:, :], self.base, levels1, chgline,
                    isovalues, colormap, cbar_label1, a_range, b_range, rectangle,
                    cellplot, x_ticks, y_ticks, **kwargs
                )
            else:
                fig, ax = plot_2Dscalar(
                    fig, ax, self.data[:, :, 0], self.base, levels1, chgline,
                    isovalues, colormap, cbar_label1, a_range, b_range, rectangle,
                    cellplot, x_ticks, y_ticks,  **kwargs
                )
        elif option.lower() == 'spin':
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            fig, ax = plot_2Dscalar(
                fig, ax, self.data[:, :, 1], self.base, levels2, spinline,
                isovalues, colormap, cbar_label2, a_range, b_range, rectangle,
                cellplot, x_ticks, y_ticks,  **kwargs
            )
        else:
            raise ValueError("Unknown option: '{}'.".format(option))

        # range and labels
        if option.lower() == 'both':
            if self.unit.lower() == 'angstrom':
                ax[0].set_xlabel(r'$\AA$')
                ax[0].set_ylabel(r'$\AA$')
                ax[1].set_xlabel(r'$\AA$')
            else:
                ax[0].set_xlabel('Bohr')
                ax[0].set_ylabel('Bohr')
                ax[1].set_xlabel('Bohr')
            if add_title == True:
                ax[0].set_title('Charge Density')
                ax[1].set_title('Spin Density')
        else:
            if self.unit.lower() == 'angstrom':
                ax.set_xlabel(r'$\AA$')
                ax.set_ylabel(r'$\AA$')
            else:
                ax.set_xlabel('Bohr')
                ax.set_ylabel('Bohr')
            if add_title == True:
                if option.lower() == 'charge':
                    ax.set_title('Charge Density')
                else:
                    ax.set_title('Spin Density')

        # restore old unit
        self._set_unit(uold)
        return fig, fig.axes

    def _set_unit(self, unit):
        """
        Set units of data of ``ChargeDensity`` object.

        Args:
            unit (str): 'Angstrom', :math:`e.\\AA^{-3}`.
                'a.u.', :math:`e.Bohr^{-3}`.
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom

        if unit.lower() == self.unit.lower():
            return self

        if unit.lower() == 'angstrom':
            self.unit = 'Angstrom'
            cst = au_to_angstrom(1.)
        elif unit.lower() == 'a.u.':
            self.unit = 'a.u.'
            cst = angstrom_to_au(1.)
        else:
            raise ValueError('Unknown unit.')

        self.base = self.base * cst
        self.data = self.data / cst**3
        return self
