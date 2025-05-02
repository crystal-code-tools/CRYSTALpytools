#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A post-processing module for electronic properties
"""
from CRYSTALpytools import units
import numpy as np


class ElectronBand():
    """
    .. _ref-ElectronBand:

    The class for electron band projected along high symmetry lines. For '3D
    band structures', please refer to the :ref:`electronics.FermiSurface <ref-FermiSurface>`
    class. Energy unit: eV. Fermi energy is aligned to 0.

    Args:
        spin (int): 1, closed shell; 2, open shell
        tick_pos (array): 1\*nTick array of 1D tick coordinates. Unit: :math:`\\AA`
        tick_label (list): 1\*nTick of default tick labels.
        efermi (float): Fermi energy. Unit: eV.
        bands (array): nBand\*nKpoint\*nSpin array of energy of band structure.
            Aligned to :math:`E_{F}=0`. Unit: eV.
        k_path (array): 1D coordinates of k points for band structure. Unit: :math:`\\AA`.
        geometry (Structure): Pymatgen structure
        reciprocal_latt (array): 3\*3 array of reciprocal lattice matrix.
        tick_pos3d (array): 1\*nTick 3D fractional tick coordinates.
        k_path3d (array): nKpoints\*3 3D fractional coordinates of k points.
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
            if np.all(reciprocal_latt!=None): self.reciprocal_latt = np.array(reciprocal_latt)
            else: self.reciprocal_latt = reciprocal_latt
        self.tick_pos3d = tick_pos3d
        self.k_path3d = k_path3d
        if np.all(self.tick_pos3d!=None):
            self.tick_pos3d = np.array(tick_pos3d, dtype=float, ndmin=2)
        if np.all(self.k_path3d!=None):
            self.k_path3d = np.array(k_path3d, dtype=float, ndmin=2)
        self.unit = unit

    @classmethod
    def from_file(cls, band, output=None):
        """
        Generate ``ElectronBand`` object from CRYSTAL BAND.DAT / fort.25 file.
        Energy unit: eV. Fermi energy is aligned to 0.

        Args:
            band (str): Name of BAND.DAT / fort.25 file.
            output (str): Properties output file.
        Returns:
            cls (ElectronBand)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_electron_band(band)

    def plot(self, **kwargs):
        """
        For band structure plotting, it is a wrapper of
        :ref:`plot.plot_electron_bands() <ref-plotebands>` with option 'single'
        (i.e., single system). For input arguments or plotting multiple systems,
        check documentations there.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``ElectronBand`` object). Check documents for
                :ref:`plot.plot_electron_bands() <ref-plotebands>`.
        Returns:
            fig (Figure): Matplotlib figure object
        """
        from CRYSTALpytools.plot import plot_electron_bands

        kwargs['mode'] = 'single'
        fig = plot_electron_bands(self, **kwargs)
        return fig

    @property
    def bandgap(self):
        """
        A shortcut for band gap. Unit is consistent with the ``self.unit`` attribute.
        """
        return self.get_bandgap()[0]

    def get_bandgap(self):
        """
        Get band gap. For spin-polarized systems, 2\*1 arrays are used for
        :math:`\\alpha` and :math:`\\beta` states. Data is rounded to 6 decimal
        places. Unit is consistent with the ``self.unit`` attribute.

        Returns:
            self.gap (float): Band gap.
            self.vbm (flat): Valence band maximum, with reference to Fermi
                level.
            self.cbm (float): Conduction band minimum, with reference to Fermi
                level.
            self.gap_pos (array): Coordinates of vbm (1st element) and cbm
                (2nd element) 3D coordinates are returned if ``self.k_path3d``
                is available. For spin-polarized cases, ``self.gap_pos[0]``
                are vbm and cbm of :math:`\\alpha` state.
        """
        import numpy as np

        self.gap = np.zeros([2,], dtype=float)
        self.vbm = np.zeros([2,], dtype=float)
        self.cbm = np.zeros([2,], dtype=float)
        if np.all(self.k_path3d!=None):
            self.gap_pos = np.zeros([2, 2, 3], dtype=float)
        else:
            self.gap_pos = np.zeros([2, 2], dtype=float)

        for ispin in range(self.spin):
            nvb = -1; ncb = -1
            for nbd, bd in enumerate(self.bands[:, :, ispin]):
                if np.all(bd>0):
                    nvb = nbd - 1; ncb = nbd; break
                else:
                    continue
            if nvb < 0 or ncb < 0:
                raise ValueError("Cannot find VB/CB. All the bands are below/over Fermi level.")
            ivbm = np.argmax(self.bands[nvb, :, ispin])
            icbm = np.argmin(self.bands[ncb, :, ispin])
            vbm = np.round(self.bands[nvb, ivbm, ispin], 6)
            cbm = np.round(self.bands[ncb, icbm, ispin], 6)
            if np.all(self.k_path3d!=None):
                kvbm = self.k_path3d[ivbm]
                kcbm = self.k_path3d[icbm]
            else:
                kvbm = self.k_path[ivbm]
                kcbm = self.k_path[icbm]

            if vbm > 0 or cbm < 0: gap = 0.
            else: gap = cbm - vbm
            self.gap[ispin] = gap
            self.vbm[ispin] = vbm
            self.cbm[ispin] = cbm
            self.gap_pos[ispin] = np.array([kvbm, kcbm])

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
        from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
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
            unit (str): 'eV': Energy unit = eV, Length unit = :math:`\\AA^{-1}`;
                'a.u.': Energy unit = Hartree. Length unit = Bohr:math:`^{-1}`
        """
        from CRYSTALpytools.units import (H_to_eV, angstrom_to_au,
                                          au_to_angstrom, eV_to_H)

        if unit.lower() == self.unit.lower():
            return self

        opt_e_props = ['gap', 'vbm', 'cbm']  # Optional energy properties
        opt_d_props = ['gap_pos']  # Optional distance properties, reciprocal
        if unit.lower() == 'ev':
            self.unit = 'eV'
            self.bands = H_to_eV(self.bands)
            self.efermi = H_to_eV(self.efermi)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, H_to_eV(attrv))
            # reciprocal
            self.tick_pos = angstrom_to_au(self.tick_pos)
            self.k_path = angstrom_to_au(self.k_path)
            if np.all(self.reciprocal_latt!=None):
                self.reciprocal_latt = angstrom_to_au(self.reciprocal_latt)
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, angstrom_to_au(attrv))
        elif unit.lower() == 'a.u.':
            self.unit = 'a.u.'
            self.bands = eV_to_H(self.bands)
            self.efermi = eV_to_H(self.efermi)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
            # reciprocal
            self.tick_pos = au_to_angstrom(self.tick_pos)
            self.k_path = au_to_angstrom(self.k_path)
            if np.all(self.reciprocal_latt!=None):
                self.reciprocal_latt = au_to_angstrom(self.reciprocal_latt)
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, au_to_angstrom(attrv))
        else:
            raise ValueError('Unknown unit.')
        return self


class FermiSurface():
    """
    .. _ref-FermiSurface:

    The class for electron band sampled across the first brillouin zone (1BZ),
    aka :math:`E(k)`. For 'normal' band structures, please refer to the
    :ref:`electronics.ElectronBand <ref-ElectronBand>` class. Energy unit: eV.
    Fermi energy is aligned to 0. For 3D and 2D systems only.

    .. note::

        The data grid is defined and saved in reciprocal unit cell, rather than
        1BZ. To be consistent with real-space data grids, the grid is non-periodic,
        i.e., the element at the boundary is repeated.

    Args:
        geometry (array|CStructure): Matrix of reciprocal lattice, or an
            extended pymatgen Structure object.
        bands (array): nBand\*nZ\*nY\*nX\*nSpin of band energy aligned to
            :math:`E_{F}=0`. The non-periodic direction must have the dimension
            of 1. Unit: eV.
        efermi (float): Fermi energy. Unit: eV.
        unit (str): In principle, should always be 'eV': eV-Angstrom.
    Returns:
        self (FermiSurface): Attributes listed below.
        self.rlattice (array): Reciprocal lattice.
        self.dimension (int): 2 or 3
        self.spin (int): 1 or 2
        self.bands (array): Bands. Dimensionality same as input.
        self.efermi (float)
        self.BZ (list): Cartesian coordinates of vertices of 1BZ. Arranged in
            planes (i.e., the element of 4\*3 array represents a plane of 4 vertices).
        self.unit (str): 'eV' or 'a.u.'
    """
    def __init__(self, geometry, bands, efermi=0., unit='eV'):
        import numpy as np
        from scipy.spatial import Voronoi
        from pymatgen.core.structure import Structure

        if isinstance(geometry, Structure):
            self.rlattice = geometry.lattice.reciprocal_lattice
        else:
            self.rlattice = np.array(geometry)

        self.efermi = efermi
        self.bands = np.array(bands, dtype=float)
        if len(self.bands.shape) != 5:
            raise ValueError('Wrong dimensionalities of input band. It must be a nBand*nX*nY*nZ*nSpin array.')
        self.spin = self.bands.shape[-1]
        idx = np.where(np.array(self.bands.shape[1:4])>1)[0]
        self.dimension = len(idx)
        # vertices of 1BZ
        self.BZ = []
        if self.dimension == 2:
            idx = 2 - idx[::-1] # to X,Y,Z
            va = self.rlattice[idx[0], idx]
            vb = self.rlattice[idx[1], idx]
            vor = Voronoi(np.array([
                -va-vb, -vb, va-vb, -va, [0, 0], va, -va+vb, vb, va+vb
            ], dtype=float))
            for i in vor.ridge_dict:
                if i[0] == 4 or i[1] == 4:
                    tmpk = []
                    for j in vor.ridge_dict[i]:
                        tmp = np.zeros([3,])
                        tmp[idx[0]] = vor.vertices[j][0]
                        tmp[idx[1]] = vor.vertices[j][1]
                        tmpk.append(tmp)
                    self.BZ.append(np.array(tmpk, dtype=float))
        elif self.dimension == 3:
            va = self.rlattice[0]
            vb = self.rlattice[1]
            vc = self.rlattice[2]
            vor = Voronoi(np.array([
                -va-vb-vc, -vb-vc, va-vb-vc, -va-vc, -vc, va-vc, -va+vb-vc, vb-vc, va+vb-vc,
                -va-vb, -vb, va-vb, -va, [0, 0, 0], va, -va+vb, vb, va+vb,
                -va-vb+vc, -vb+vc, va-vb+vc, -va+vc, vc, va+vc, -va+vb+vc, vb+vc, va+vb+vc,
            ], dtype=float))
            for i in vor.ridge_dict:
                if i[0] == 13 or i[1] == 13:
                    self.BZ.append([vor.vertices[j] for j in vor.ridge_dict[i]])
        else:
            raise Exception('3D band structure requires 3D/2D grid points. The dimensionality of input grid is {:d} {:d} {:d}'.format(
                k_points[0], k_points[1], k_points[2]))

        self.unit = unit

    @classmethod
    def from_file(cls, band, output=None):
        """
        Generate ``FermiSurface`` object from CRYSTAL fort.35 file. Energy
        unit: eV. Fermi energy is aligned to 0.

        Args:
            band (str): Name of fort.35 files.
            output (str): Properties output file.
        Returns:
            cls (FermiSurface)
        """
        from CRYSTALpytools.crystal_io import Properties_output
        return Properties_output(output).read_Fermi_surface(band)

    @property
    def bandgap(self):
        """
        A shortcut for band gap. Unit is consistent with the ``self.unit`` attribute.
        """
        return self.get_bandgap()[0]

    def get_bandgap(self):
        """
        Get band gap. For spin-polarized systems, 2\*1 arrays are used for
        :math:`\\alpha` and :math:`\\beta` states. Data is rounded to 6 decimal
        places. Unit is consistent with the ``self.unit`` attribute.

        Returns:
            self.gap (float): Band gap.
            self.vbm (flat): Valence band maximum, with reference to Fermi
                level.
            self.cbm (float): Conduction band minimum, with reference to Fermi
                level.
            self.gap_pos (array): Coordinates of vbm (1st element) and cbm
                (2nd element). For spin-polarized cases, ``self.gap_pos[0]``
                are vbm and cbm of :math:`\\alpha` state.
        """
        import numpy as np

        self.gap = np.zeros([2,], dtype=float)
        self.vbm = np.zeros([2,], dtype=float)
        self.cbm = np.zeros([2,], dtype=float)
        self.gap_pos = np.zeros([2, 2, 3], dtype=float)

        nz = self.bands.shape[1]-1 # general grid repeats the last element
        ny = self.bands.shape[2]-1
        nx = self.bands.shape[3]-1
        fz = np.linspace(0, 1, nz, endpoint=False)
        fy = np.linspace(0, 1, ny, endpoint=False)
        fx = np.linspace(0, 1, nx, endpoint=False)
        for ispin in range(self.spin):
            nvb = -1; ncb = -1
            for nbd, bd in enumerate(self.bands[:, :-1, :-1, :-1, ispin]):
                bd = np.round(bd, 4) # lower the accuracy
                if np.all(bd>=0):
                    nvb = nbd - 1; ncb = nbd; break
                else:
                    continue
            if nvb < 0 or ncb < 0:
                raise ValueError("Cannot find VB/CB. All the bands are below/over Fermi level.")
            fvbm = self.bands[nvb, :-1, :-1, :-1, ispin].flatten()
            fcbm = self.bands[ncb, :-1, :-1, :-1, ispin].flatten()
            ivbm = np.argmax(fvbm)
            icbm = np.argmin(fcbm)
            vbm = np.round(fvbm[ivbm], 4)
            cbm = np.round(fcbm[icbm], 4)

            izv = int(ivbm // (ny*nx))
            iyv = int((ivbm - izv*ny*nx) // nx)
            ixv = int(ivbm - izv*ny*nx - iyv*nx)
            kvbm = np.array([fx[ixv], fy[iyv], fz[izv]]) @ self.rlattice

            izc = int(icbm // (ny*nx))
            iyc = int((icbm - izc*ny*nx) // nx)
            ixc = int(icbm - izc*ny*nx - iyc*nx)
            kcbm = np.array([fx[ixc], fy[iyc], fz[izc]]) @ self.rlattice

            gap = np.round(cbm-vbm, 4)
            if gap < 0: gap = 0.
            self.gap[ispin] = gap
            self.vbm[ispin] = vbm
            self.cbm[ispin] = cbm
            self.gap_pos[ispin] = np.array([kvbm, kcbm])

        if self.spin == 1:
            self.gap = self.gap[0]
            self.vbm = self.vbm[0]
            self.cbm = self.cbm[0]
            self.gap_pos = self.gap_pos[0]
        return self.gap, self.vbm, self.cbm, self.gap_pos

    def plot(self,
             band_index='vb',
             isovalue=0.,
             volume_3d=False,
             interp='no interp',
             interp_size=1,
             BZ_plot=True,
             BZ_scale=1.0,
             BZ_color=(0., 0., 0.),
             BZ_linewidth=1.0,
             tick_pos=[],
             tick_label=[],
             show_the_scene=True,
             **kwargs):
        """
        Plot :math:`E(k)` in the first brillouin zone (1BZ).

        .. note ::

            `MayaVi <https://docs.enthought.com/mayavi/mayavi/>`_ is used to
            display :math:`E(k)`, which is not installed by default.

        * For 3D systems, it is displayed as isosurfaces of :math:`E(k)`, or
            its distributions if ``volume_3d=True``.  
        * For 2D systems, the distribution is displayed, so ``isovalue`` is
            disabled.

        For 3D systems, displaying multiple bands at multiple isovalues is
        discouraged. But the user can still visualize the isosurfaces of
        multiple bands and the same ``isovalue`` and ``colormap`` applies to
        all the bands. A barchart is plotted by matplotlib to indicate band
        energy ranges and isovalues.

        If ``volume_3d=True``, only 1 band is allowed. ``isovalue`` is ignored.
        The grid in reciprocal lattice should be, ideally, orthogonal and
        aligned to X Y and Z axes. If not, the grid is linearly interpolated by
        the `scipy.interpolate.griddata() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html>`_
        method. More details are available in the :ref:`base.plotbase.tvtkGrid() <ref-tvtkGrid>`
        function. The interpolated grid has the same size as specified by
        ``interp_size``.

        .. note ::

            Too large ``interp_size`` might be very memory and cpu demanding!
            Similarly, ``volume_3d=True`` might be very demanding in some cases!

        Args:
            band_index (list|str|int): Indices of bands to plot. Starting from
                1. For spin-polarized cases, one can specify '4a' for the
                :math:`\\alpha` state of band 4. Use 'vb' and 'cb' to display
                the highest valance or the lowest conduction band.
            isovalue (float|list): *3D system only* Isovalue of surfaces.
            volume_3d (bool): *3D system only* Plot 3D volume data rather than
                isosurfaces. Only 1 band is permitted.
            interp (str): Interpolate data to smoothen the plot. 'no interp' or
                'linear', 'nearest', 'slinear', 'cubic'. please refer to
                `scipy.interpolate.interpn <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interpn.html>`_
                 The interpolated data is not saved.
            interp_size (list[int]|int): The new size of interpolated data
                (list) or a scaling factor. *Valid only when ``interp`` is not
                'no interp'*.
            BZ_plot (bool): Whether to plot the boundary of 1BZ.
            BZ_scale (float): *Must be between 1 and 2*. Do not truncate data grid
                at the boundary of 1BZ, slightly expanding it by the scaling
                factor. This does not change the 1BZ frame plotted.
            BZ_color (array): Color of the 1BZ plot. See
                `Mayavi colormap <https://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html>`_.
            BZ_linewidth (float): Linewidth of the 1BZ plot.
            tick_pos (array): *Not implemented*
            tick_label (list): *Not implemented*
            show_the_scene (bool):  Display the scene by ``mlab.show()`` and
                return None. Otherwise return the scene object.
            \*\*kwargs: Other keywords passed to MayaVi. See below.
            colormap (str): Colormap of surfaces / isosurfaces. *Not available
                for ``volume_3d=True``*. Default is 'jet'.
            opacity (float):  Opacity from 0 to 1. For ``volume_3d=True``, that
                defines the opacity of the maximum value. The opacity of the
                minimum is half of it.
            transparent (bool): Scalar-dependent opacity. *Not for volume_3d=True*.
            vmax (float): Maximum value of colormap. Default is max of displayed bands.
            vmin (float): Minimum value of colormap. Default is min of displayed bands.
            warp_scale (float): *2D only* The length along energy (z) axis. Default is 1.
            azimuth: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            elevation: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            distance: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
                By default set to 'auto'.
            focalpoint: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            roll: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
        Returns:
            fig: MayaVi scence object, if ``show_the_scene=False``.
        """
        import copy, warnings, re
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.patches import Rectangle
        from CRYSTALpytools.base.plotbase import GridInterpolate
        try:
            from mayavi import mlab
            from tvtk.util.ctf import PiecewiseFunction
        except ModuleNotFoundError:
            raise ModuleNotFoundError('MayaVi is required for this functionality, which is not in the default dependency list of CRYSTALpytools.')
        from CRYSTALpytools.base.plotbase import tvtkGrid

        #---------------------------------------------------------------------#
        #                                NOTE                                 #
        #---------------------------------------------------------------------#
        # For visualization, band has the dimension of nX*nY*nZ. A transpose  #
        # is needed!                                                          #
        #---------------------------------------------------------------------#

        # Input processing and sanity check
        k_points = np.array([self.bands.shape[3], self.bands.shape[2], self.bands.shape[1]], dtype=int)
        iband, ispin = FermiSurface._get_band_index(self.bands, band_index)
        spin_label = {0 : 'alpha', 1 : 'beta'}

        isovalue = np.array(isovalue, dtype=float, ndmin=1)
        if self.dimension == 3 and volume_3d == False:
            if len(iband) > 1 and len(isovalue) > 1:
                    warnings.warn('For 3D systems, the user is strongly recommended to plot only 1 band or 1 isovalue every time. The same isovalue is applied to all the bands.',
                                  stacklevel=2)
        if self.dimension == 3 and volume_3d == True:
            if len(iband) > 1: raise Exception("Only 1 band is permitted when 'volume_3d=True'.")

        # Define and scale 1BZ
        k_path = self.BZ
        nedge = len(k_path)
        k_path = [np.array(i)*BZ_scale for i in k_path]
        norms = np.zeros([nedge, 3]) # norm vector of 1BZ edges
        dists = np.zeros([nedge,]) # projection of (origin-vertex) onto norm vector
        verteices = np.zeros([nedge, 3]) # vertex of vector (origin-vertex)
        for ik, k in enumerate(k_path):
            if self.dimension == 3: norm = np.cross(k[1]-k[0], k[2]-k[0])
            else: norm = np.cross(k[1]-k[0], [0, 0, 1.])
            norm = norm / np.linalg.norm(norm)
            norms[ik] = norm
            dists[ik] = np.dot(-k[0], norm)
            verteices[ik] = k[0]

        # Energy scale
        emax = []; emin = []
        for ibd, isp in zip(iband, ispin):
            emax.append(np.max(self.bands[ibd, :, :, :, isp]))
            emin.append(np.min(self.bands[ibd, :, :, :, isp]))
        allmax = np.max(emax); allmin = np.min(emin)
        ## Add bar chart for 3D isosurfaces
        if self.dimension == 3 and volume_3d == False:
            barfig, barax = plt.subplots(1, 1, layout='tight', figsize=[6, 4.5])
            clist = list(mcolors.TABLEAU_COLORS.keys())
            barpos = np.array([i*2+1 for i in range(len(iband))])
            barlab = []
            for i in range(len(iband)):
                rect = Rectangle((barpos[i], emin[i]), 1, emax[i]-emin[i], color=clist[i%10])
                barax.add_patch(rect)
                barax.text(barpos[i]+0.5, emax[i], '{:.2f}'.format(emax[i]),
                           horizontalalignment='center', verticalalignment='bottom')
                barax.text(barpos[i]+0.5, emin[i], '{:.2f}'.format(emin[i]),
                           horizontalalignment='center', verticalalignment='top')
                barlab.append('Band {:d}\nSpin $\{}$'.format(iband[i]+1, spin_label[ispin[i]]))
            barax.hlines(isovalue, 0, barpos[-1]+2, colors='tab:gray', linestyles='dotted')
            barax.set_xticks(barpos+0.5, labels=barlab)
            barax.set_xlim([0, barpos[-1]+2])
            barax.set_yticks(isovalue)
            barax.set_ylabel(r'$E_{iso}-E_{F}$ (eV)')
            plt.show()

        # 2D, Scale Z axis for 2D, same length as the longer one of x and y
        if self.dimension == 2:
            xrange = np.max(verteices[:, 0])-np.min(verteices[:, 0])
            yrange = np.max(verteices[:, 1])-np.min(verteices[:, 1])
            zscale = np.max([xrange, yrange]) / (allmax - allmin)
        # 2D, indices of periodic and non-periodic dirs
        prdd = np.where(k_points>1)[0]
        isod = np.where(k_points==1)[0]
        if len(isod) == 1: isod = isod[0]
        elif len(isod) == 0: isod = -1
        else: raise Exception("A 2D/3D k mesh must defined. Now only 1 k point is defined along {:d} dimensions.".format(len(isod)))

        # Interpolation
        interp = interp.lower()
        if interp not in ['no interp', 'linear', 'nearest', 'slinear', 'cubic']:
            raise ValueError("Unknown interpolation method : '{}'.".format(interp))
        if volume_3d == True:
            interp = 'no interp' # not using GridInterpolate.

        # Fig scale, 1BZ might be very small for large systems.
        # Get a self-adaptative scaling factor. Applied only on plotted scenes
        fig_scale = 2 / np.max([np.linalg.norm(self.rlattice[0]),
                                np.linalg.norm(self.rlattice[1]),
                                np.linalg.norm(self.rlattice[2])])
        # Figure
        fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))

        # Band data
        for ibd, isp, vmx, vmi in zip(iband, ispin, emax, emin):
            # Range of isovalues
            if volume_3d == False and self.dimension == 3:
                isoplot = isovalue[np.where((isovalue>vmi)&(isovalue<vmx))[0]]
                if isoplot.shape[0] == 0:
                    warnings.warn("Band {:d}, spin {}: No isovalue lies in the energy range.".format(ibd+1, spin_label[isp]),
                                  stacklevel=2)
                    continue
            # Interpolate band
            if interp != 'no interp':
                if self.dimension == 3:
                    base = np.vstack([[0, 0, 0], self.rlattice])
                    intband, _ = GridInterpolate(
                        base,
                        np.transpose(self.bands[ibd, :, :, :, isp], axes=[2,1,0]),
                        interp,
                        interp_size,
                    )
                    intkpts = np.multiply(interp_size, k_points)
                else:
                    base = np.vstack([[0, 0, 0], self.rlattice[prdd]])
                    intband, _ = GridInterpolate(
                        base,
                        np.transpose(self.bands[ibd, :, :, :, isp], axes=[2,1,0]).reshape([k_points[prdd[0]], k_points[prdd[1]]]),
                        interp,
                        interp_size
                    )
                    intkpts = np.multiply(interp_size, k_points)
                    intkpts[isod] = 1
                    intband = intband.reshape(intkpts)
            else:
                intband = np.transpose(self.bands[ibd, :, :, :, isp], axes=[2,1,0])
                intkpts = k_points
            # Expand to a 2x2 supercell, note it's a general, non-periodic grid
            pband = np.zeros(intkpts*2-1)
            pband[:intkpts[0], :intkpts[1], :intkpts[2]] = intband
            pband[intkpts[0]:, :intkpts[1], :intkpts[2]] = intband[1:, :, :]
            pband[:intkpts[0], intkpts[1]:, :intkpts[2]] = intband[:, 1:, :]
            pband[:intkpts[0], :intkpts[1], intkpts[2]:] = intband[:, :, 1:]
            pband[intkpts[0]:, intkpts[1]:, :intkpts[2]] = intband[1:, 1:, :]
            pband[intkpts[0]:, :intkpts[1], intkpts[2]:] = intband[1:, :, 1:]
            pband[:intkpts[0], intkpts[1]:, intkpts[2]:] = intband[:, 1:, 1:]
            pband[intkpts[0]:, intkpts[1]:, intkpts[2]:] = intband[1:, 1:, 1:]
            del intband

            if isod == 0: pband = pband[0]
            elif isod == 1: pband = pband[:, 0, :]
            elif isod == 2: pband = pband[:, :, 0]
            kptnew = np.array(pband.shape, dtype=int)

            if self.dimension == 3: # 3D plot
                fracx = np.linspace(-1, 1, kptnew[0])
                fracy = np.linspace(-1, 1, kptnew[1])
                fracz = np.linspace(-1, 1, kptnew[2])
                cartx = fracx.reshape([-1, 1]) @ self.rlattice[0].reshape([1, 3])
                carty = fracy.reshape([-1, 1]) @ self.rlattice[1].reshape([1, 3])
                cartz = fracz.reshape([-1, 1]) @ self.rlattice[2].reshape([1, 3])

                if volume_3d == False: nullvalue = np.nan
                else: nullvalue = allmin-1.

                for i in range(kptnew[0]):
                    x = cartx[i]
                    for j in range(kptnew[1]):
                        y = carty[j]
                        for k in range(kptnew[2]):
                            z = cartz[k]
                            v =  x + y + z - verteices
                            # if the point on the same side of 1BZ boundary as origin
                            dist = np.array([
                                np.dot(v[l], norms[l]) * dists[l] for l in range(nedge)
                            ])
                            if np.all(dist>-1e-4): continue
                            pband[i, j, k] = nullvalue
                # Plot isosurface
                if volume_3d == False:
                    grid = tvtkGrid([[0., 0., 0.],
                                     self.rlattice[0]*2*fig_scale, # pband defined on 2x2x2 grid
                                     self.rlattice[1]*2*fig_scale,
                                     self.rlattice[2]*2*fig_scale],
                                    pband,
                                    True,
                                    None)

                    keys = ['colormap', 'opacity', 'transparent', 'vmax', 'vmin']
                    keywords = dict(figure=fig,
                                    contours=isoplot.tolist(),
                                    colormap='jet',
                                    vmax=allmax,
                                    vmin=allmin)
                    for k, v in zip(kwargs.keys(), kwargs.values()):
                        if k in keys: keywords[k] = v

                    contour = mlab.pipeline.iso_surface(grid, **keywords)
                # Plot volume
                else:
                    grid = tvtkGrid([[0., 0., 0.],
                                     self.rlattice[0]*2*fig_scale, # pband defined on 2x2x2 grid
                                     self.rlattice[1]*2*fig_scale,
                                     self.rlattice[2]*2*fig_scale],
                                    pband,
                                    True,
                                    interp_size,
                                    fill_value=nullvalue)

                    keys = ['vmax', 'vmin']
                    keywords = dict(figure=fig,
                                    vmax=allmax,
                                    vmin=allmin)
                    for k, v in zip(kwargs.keys(), kwargs.values()):
                        if k in keys: keywords[k] = v

                    vol = mlab.pipeline.volume(grid, **keywords)
                    # transparency
                    otf = PiecewiseFunction()
                    otf.add_point(nullvalue, 0.)
                    if 'opacity' in kwargs.keys():
                        opacity = kwargs['opacity']
                    else:
                        opacity = 1.0
                    for i, o in zip(np.linspace(allmin, allmax, 10),
                                    np.linspace(opacity*0.5, opacity, 10)):
                        otf.add_point(i, o)
                    vol._otf = otf
                    vol._volume_property.set_scalar_opacity(otf)
            else: # 2D plot
                fracx = np.linspace(-1, 1, kptnew[0])
                fracy = np.linspace(-1, 1, kptnew[1])
                cartx = fracx.reshape([-1, 1]) @ self.rlattice[prdd[0]].reshape([1, 3])
                carty = fracy.reshape([-1, 1]) @ self.rlattice[prdd[1]].reshape([1, 3])
                # mask = np.zeros_like(pband, dtype=bool)
                for i in range(kptnew[0]):
                    x = cartx[i]
                    for j in range(kptnew[1]):
                        y = carty[j]
                        v = x + y - verteices
                        # if the point on the same side of 1BZ boundary as origin
                        dist = np.array([
                            np.dot(v[l], norms[l]) * dists[l] for l in range(nedge)
                        ])
                        if np.all(dist>-1e-4): continue
                        pband[i, j] = np.nan
                        # mask[i, j] = True
                # Plot surface, did not find a masking method for pipeline surf.
                keys = ['colormap', 'opacity', 'transparent', 'vmax', 'vmin', 'warp_scale']
                keywords = dict(figure=fig,
                                # mask=mask,
                                colormap='jet',
                                vmax=allmax,
                                vmin=allmin,
                                warp_scale=1)
                for k, v in zip(kwargs.keys(), kwargs.values()):
                    if k in keys: keywords[k] = v

                surf = mlab.surf(pband, **keywords)
                # Non-orthogonal grid, points are shifted by -0.5 for unknown reasons
                polydata = surf.actor.actors[0].mapper.input
                pts = np.array(polydata.points)
                dlatt = np.zeros([2, 2])
                for i, pdir in enumerate(prdd):
                    dlatt[i, :] = 2 * self.rlattice[pdir, :][prdd] / (kptnew[i]-1)
                pts[:, prdd] = (pts[:, prdd]+0.5) @ dlatt
                pts[:, isod] = pts[:, isod] * zscale
                polydata.points = pts * fig_scale

        # Plot 1BZ
        if BZ_plot == True:
            if self.dimension == 3:
                for kp in k_path:
                    kplt = np.vstack([kp, kp[0]]) / BZ_scale * fig_scale
                    mlab.plot3d(kplt[:, 0], kplt[:, 1], kplt[:, 2],
                                figure=fig, color=BZ_color, line_width=BZ_linewidth)
            else:
                for kp in k_path:
                    kplt = kp / BZ_scale * fig_scale
                    mlab.plot3d(kplt[:, 0], kplt[:, 1], kplt[:, 2],
                                figure=fig, color=BZ_color, line_width=BZ_linewidth)

        # colorbar
        mlab.scalarbar(title='E-Ef ({})'.format(self.unit),
                       orientation='vertical',
                       nb_labels=5,
                       label_fmt='%.2f')


        keys = ['azimuth', 'elevation', 'distance', 'focalpoint', 'roll']
        keywords = dict(figure=fig, distance='auto')
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v
        mlab.view(**keywords)

        if show_the_scene == False:
            return fig
        else:
            mlab.gcf().scene.parallel_projection = True
            mlab.show()
            return

    def to_bxsf(self, filename, band_index=[]):
        """
        **3D data only**

        Write data into the `XCrySDen <http://www.xcrysden.org/>`_ BXSF file
        format. Unit: :math:`\\AA^{-1}` and eV. Band labels:

        * \[0\:-2\]: Band index strating from 1;  
        * -2: Useless, 0;  
        * -1: Spin. Alpha / No spin, 1; Beta, 2

        Args:
            filename (str): The output xsf filename.
            band_index (list|str|int): Indices of bands to be saved. Same as
                the ``band_index`` option of the ``plot()`` method. '[]' for all
                bands.
        Returns:
            None
        """
        from CRYSTALpytools.io.xcrysden import BXSF

        if self.dimension != 3:
            raise Exception('BXSF only supports 3D band data.')

        self._set_unit('eV')
        BXSF(self.rlattice, self.bands,
             efermi=self.efermi, band_index=band_index).write(filename)
        return

    @classmethod
    def _get_band_index(cls, band, index):
        """
        The method to process band indices.

        Args:
            band (array): nBand\*nZ\*nY\*nX*nSpin array.
            index (list|str|int): See the ``plot()`` method.
        Returns:
            iband (array): 1\*nIndex indices of band, starting from 0.
            ispin (array): 1\*nIndex indices of spin, 0 or 1.
        """
        import re, warnings
        import numpy as np

        index = np.unique(np.array(index, dtype=str, ndmin=1))
        iband = []; ispin = []; spin = band.shape[-1]
        for i in index:
            j = re.findall(r'[0-9]+', i)
            if len(j) == 0: # cb/vb
                tmp = np.zeros([spin, 4], dtype=float) - 1 # row 0: Spin up, [VB, eVB, CB, eCB]
                for isp in range(spin):
                    for nbd, bd in enumerate(band[:, :, :, :, isp]):
                        if np.all(bd>0):
                            tmp[isp, 0] = nbd - 1
                            tmp[isp, 1] = np.max(band[nbd-1, :, :, :, isp])
                            tmp[isp, 2] = nbd
                            tmp[isp, 3] = np.min(band[nbd, :, :, :, isp])
                            break
                        else:
                            continue
                if np.any(tmp[:,0]<0) or np.any(tmp[:,2]<0):
                    raise Exception('All the bands are over/below Fermi level. VBM/CBM not found.')

                if i.lower() == 'vb':
                    ivb = np.argmax(tmp[:,1])
                    iband.append(tmp[ivb,0]); ispin.append(ivb)
                elif i.lower() == 'cb':
                    icb = np.argmin(tmp[:,3])
                    iband.append(tmp[icb,2]); ispin.append(icb)
                else:
                    raise ValueError("Unknown band specification: '{}'.".format(i))
            else:
                j = int(j[0])
                if 'a' in i.lower():
                    iband.append(j-1); ispin.append(0)
                    if band.shape[-1]==1:
                        warnings.warn("Not a spin-polarized system. Input 'a'/'b' are ignored.",
                                      stacklevel=2)
                # in case of 'ab' as input
                if 'b' in i.lower():
                    if band.shape[-1]==1:
                        warnings.warn("Not a spin-polarized system. Input 'a'/'b' are ignored.",
                                      stacklevel=2)
                        iband.append(j-1); ispin.append(0)
                    else:
                        iband.append(j-1); ispin.append(1)
                if 'a' not in i.lower() and 'b' not in i.lower():
                    if band.shape[-1]==1:
                        iband.append(j-1); ispin.append(0)
                    else:
                        iband.append(j-1); ispin.append(0)
                        iband.append(j-1); ispin.append(1)

        iband = np.array(iband, dtype=int)
        ispin = np.array(ispin, dtype=int)
        return iband, ispin

    def _set_unit(self, unit):
        """
        Set units of data of ``FermiSurface`` object. Internal method.

        Args:
            unit (str): 'eV': Energy unit = eV, Length unit = :math:`\\AA^{-1}`;
                'a.u.': Energy unit = Hartree. Length unit = Bohr:math:`^{-1}`
        """
        from CRYSTALpytools.units import H_to_eV, angstrom_to_au, au_to_angstrom, eV_to_H

        if unit.lower() == self.unit.lower():
            return self
        opt_e_props = ['gap', 'vbm', 'cbm']  # Optional energy properties
        opt_d_props = ['gap_pos']  # Optional distance properties, reciprocal
        if unit.lower() == 'ev':
            self.unit = 'eV'
            self.bands = H_to_eV(self.bands)
            self.efermi = H_to_eV(self.efermi)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, H_to_eV(attrv))
            # reciprocal
            self.BZ = [angstrom_to_au(i) for i in self.BZ]
            self.rlattice = angstrom_to_au(self.rlattice)
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, angstrom_to_au(attrv))
        elif unit.lower() == 'a.u.':
            self.unit = 'a.u.'
            self.bands = eV_to_H(self.bands)
            self.efermi = eV_to_H(self.efermi)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
            # reciprocal
            self.BZ = [au_to_angstrom(i) for i in self.BZ]
            self.rlattice = au_to_angstrom(self.rlattice)
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, au_to_angstrom(attrv))
        else:
            raise ValueError('Unknown unit.')
        return self


class ElectronDOS():
    """
    .. _ref-ElectronDOS:

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
        from CRYSTALpytools.crystal_io import Properties_output
        return Properties_output().read_electron_dos(dos)

    def plot(self, **kwargs):
        """
        A wrapper to plot density of states of a single system with matplotlib.
        For input arguments or plotting multiple systems, check
        :ref:`plot.plot_electron_doss() <ref-plotedoss>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``ElectronDOS`` object). Check documents for
                :ref:`plot.plot_electron_doss() <ref-plotedoss>`.
        Returns:
            fig (Figure): Matplotlib figure object
        """
        from CRYSTALpytools.plot import plot_electron_doss

        fig = plot_electron_doss(self, **kwargs)
        return fig

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
        opt_d_props = []  # Optional energy invers properties
        if unit.lower() == 'ev':
            self.unit = 'eV'
            self.efermi = H_to_eV(self.efermi)
            self.energy = H_to_eV(self.energy)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, H_to_eV(attrv))
            self.doss = eV_to_H(self.doss)
            for p in opt_d_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
        elif unit.lower() == 'a.u.':
            self.unit = 'a.u.'
            self.efermi = eV_to_H(self.efermi)
            self.energy = eV_to_H(self.energy)
            for p in opt_e_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, eV_to_H(attrv))
            self.doss = H_to_eV(self.doss)
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
        A wrapper to plot electron band structure + density of states of a
        single system with matplotlib. For input arguments, check
        :ref:`plot.plot_electron_banddos() <ref-plot>`.

        Args:
            \*\*kwargs: Plot setting parameters (i.e., except the variable for
                ``ElectronBandDOS`` object). Check documents for
                :ref:`plot.plot_electron_banddos() <ref-plot>`.
        Returns:
            fig (Figure): Matplotlib figure object
        """
        from CRYSTALpytools.plot import plot_electron_banddos

        fig = plot_electron_banddos(self, **kwargs)
        return fig

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
    Charge (spin) density object. Unit: :math:`e.\\AA^{-3}`.

    .. note::

        Definition follows the convention of Gaussian CUBE and XCrySDen XSF
        formats, which requires a non-periodic grid defined over \[0, 1\].

    Args:
        data (array): Plot data. nY\*nX\*nSpin (2D) or nZ\*nY\*nX\*nSpin (3D)
        base (array): 3(4)\*3 Cartesian coordinates of the 3(4) points defining
            base vectors BA, BC (2D) or OA, OB, OC (3D). The sequence is (O),
            A, B, C.
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

    @classmethod
    def from_file(cls, *files, output=None, method='normal'):
        """
        Generate a ``ChargeDensity`` object from a single file, or from multiple
        files by subtracting values from the first entry. Can be used for
        multiple dimensions (2D and 3D).

        .. note::

            The standard screen output is required to identify the indices of
            corresponding 2D data maps. Otherwise the code only reads the
            first 1 (2) 2D data maps for spin = 1 (2).

        Available methods are:

        * 'normal': Normal. For 2D fort.25 files, 1 entry. For 3D cube files,
            2 entries for charge and spin if needed.
        * 'subtract': Subtracting data from the first entry based on following
            entries. Multiple entries only.  
        * 'alpha_beta': Save spin-polarized data in :math:`\\alpha` /
            :math:`\\beta` states, rather than charge(:math:`\\alpha+\\beta`)
            / spin(:math:`\\alpha-\\beta`). For 2D fort.25 files, 1 entry. For
            3D cube files, 2 entries for charge and spin if needed.

        Args:
            \*files (str): Path to the charge density / spin density file(s).
                All the entries must be in the same file format.
            output (str): Screen output file.
            method (str): See above.
        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        file = open(files[0], 'r')
        header = file.readline()
        file.close()
        if '-%-' in header: # 2D plot in fort.25
            cls = Properties_output(output).read_ECHG(*files, method=method)
        else: # no identifier for CUBE files
            cls = Properties_output(output).read_ECH3(*files, method=method)
        return cls

    def subtract(self, *args):
        """
        Subtracting data of the same type from the object.

        Args:
            \*args (str|ChargeDensity): File names or ``ChargeDensity`` objects.
        Returns:
            self (ChargeDensity) : spin dimension, if there is, is not kept.
                Only charge density difference is subtracted.
        """
        from CRYSTALpytools.crystal_io import Properties_output
        import numpy as np

        for i in args:
            if isinstance(i, str):
                obj = Properties_output().read_ECHG(i, method='normal')
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
            if self.data.shape != obj.data.shape:
                raise ValueError('Inconsistent mesh grid between input and object.')
            # spin
            if self.spin != obj.spin:
                raise ValueError('Inconsistent spin dimensionalities between input and object.')
            # subtract
            self.data = self.data - obj.data

        # spin difference is not kept - meaningless. (if not please remove this block)
        if self.spin == 2:
            oshape = self.data.shape
            chglen = 1
            for s in oshape[:-1]:
                chglen = chglen * s
            chglen = int(chglen)
            self.data = self.data.flatten(order='F')[:chglen]
            nshape = [i for i in oshape[:-1]]
            nshape.append(1)
            self.data = np.reshape(self.data, nshape, order='F')
            self.spin = 1
        return self

    def substract(self, *args):
        """An old typo"""
        return self.subtract(*args)

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

    def plot_2D(self,
                unit='Angstrom',
                option='both',
                levels=150,
                lineplot=False,
                linewidth=1.0,
                isovalues=None,
                colorplot=True,
                colormap='jet',
                cbar_label='default',
                a_range=[0., 1.],
                b_range=[0., 1.],
                rectangle=False,
                edgeplot=False,
                x_ticks=5,
                y_ticks=5,
                title='default',
                figsize=[6.4, 4.8],
                fig=None,
                ax_index=None,
                **kwargs):
        """
        Plot 2D charge/spin density map. A wrapper of ``plot.plot_dens_ECHG``
        and ``plot.plot_spin_ECHG``.

        3 styles are available:

        1. ``lineplot=True`` and ``colorplot=True``: The color-filled contour
            map with black contour lines. Dotted lines for negative values and
            solid lines for positive values. The solid line twice in width for 0.  
        2. ``lineplot=False`` and ``colorplot=True``: The color-filled contour
            map.  
        3. ``lineplot=True`` and ``colorplot=False``: The color coded contour
            line map. Blue dotted line for negative values and red solid lines
            for positive values. The balck solid line twice in width for 0.

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
                format string, e.g., ``"%.2f``. Used only if ``lineplot=True``.
                None for not adding isovalues
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 1\*2 list of colorbar titles can be set for
                spin-polarized systems. 'default' for unit and symbol. 'None'
                for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in
                fractional coordinate.
            rectangle (bool): If :math:`a, b` are non-orthogonal, plot a
                rectangle region and reset :math:`b`. If used together with
                ``b_range``, that refers to the old :math:`b` (i.e., expansion first).
            edgeplot (bool): Whether to add cell edges represented by the
                original base vectors (not inflenced by a/b range or rectangle
                options).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Titles for both charge and spin densities.
                'default' for default values and 'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (list[int]): *Developer Only*, indices of axes in ``fig.axes``.
            \*\*kwargs: Other arguments passed to ``axes.contour()`` function
                to set contour lines.

        Returns:
            fig (Figure): Matplotlib figure object
            ax (Axes): Matplotlib axes object
        """
        from CRYSTALpytools.base.plotbase import plot_2Dscalar
        import numpy as np
        import matplotlib.pyplot as plt
        import warnings

        # dimen
        if self.dimension != 2:
            raise Exception('Not a 2D charge density object.')

        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)

        # levels
        levels = np.array(levels, dtype=float, ndmin=2)
        if levels.shape[1] == 1:
            if self.spin == 1:
                levels1 = np.linspace(np.min(self.data), np.max(self.data),
                                      int(levels[0, 0]))
                levels2 = levels1
            else:
                levels1 = np.linspace(np.min(self.data[:,:,0]),
                                      np.max(self.data[:,:,0]), int(levels[0, 0]))
                levels2 = np.linspace(np.min(self.data[:,:,1]),
                                      np.max(self.data[:,:,1]), int(levels[0, 0]))
        else:
            if levels.shape[0] == 1:
                levels1 = levels; levels2 = levels
            else:
                levels1 = levels[0]; levels2 = levels[1]

        # color plot
        if colorplot == False:
            colormap = None

        # contour line
        if lineplot == True:
            if colorplot == False: # colored contour lines
                chgline = []
                for j in levels1:
                    if j > 1e-6: chgline.append(['r', '-', linewidth])
                    else: chgline.append(['k', '-', linewidth*2])
                spinline = []
                for j in levels2:
                    if j > 1e-6: spinline.append(['r', '-', linewidth])
                    elif -j > 1e-6: spinline.append(['b', 'dotted', linewidth])
                    else: spinline.append(['k', '-', linewidth*2])
            else: # black contour lines
                chgline = []
                for j in levels1:
                    if j > 1e-6: chgline.append(['k', '-', linewidth])
                    else: chgline.append(['k', '-', linewidth*2])
                spinline = []
                for j in levels2:
                    if j > 1e-6: spinline.append(['k', '-', linewidth])
                    elif -j > 1e-6: spinline.append(['k', 'dotted', linewidth])
                    else: spinline.append(['k', '-', linewidth*2])
        else:
            chgline = None; spinline = None

        # cbar_label
        if np.all(cbar_label==None):
            cbar_label1 = None; cbar_label2 = None
        else:
            cbar_label = np.array(cbar_label, ndmin=1)
            if cbar_label[0].lower() == 'default':
                if unit.lower() == 'angstrom':
                    cbar_label1 = r'$\rho_{\alpha}+\rho_{\beta}$ ($|e|/\AA^{3}$)'
                    cbar_label2 = r'$\rho_{\alpha}-\rho_{\beta}$ ($|e|/\AA^{3}$)'
                else:
                    cbar_label1 = r'$\rho_{\alpha}+\rho_{\beta}$ ($|e|/Bohr^{3}$)'
                    cbar_label2 = r'$\rho_{\alpha}-\rho_{\beta}$ ($|e|/Bohr^{3}$)'
            else:
                if cbar_label.shape[0] > 1:
                    cbar_label1 = cbar_label[0]; cbar_label2 = cbar_label[1];
                else:
                    cbar_label1 = str(cbar_label[0]); cbar_label2 = str(cbar_label[0])

        # plot
        ## spin
        if self.spin == 1 and (option.lower()=='both' or option.lower()=='spin'):
            warnings.warn("Spin options not available to non spin-polarized cases.",
                          stacklevel=2)
            option = 'charge'

        ## get the correct axes
        if np.all(fig!=None):
            if np.all(ax_index==None):
                raise ValueError("Indices of axes must be set when 'fig' is not None.")
            ax_index = np.array(ax_index, dtype=int, ndmin=1)
            if option.lower() == 'both' and len(ax_index) != 2:
                raise ValueError("2 axes needed when option='both'.")
        else:
            if option.lower() == 'both':
                fig, ax = plt.subplots(1, 2, figsize=figsize, sharex=True,
                                       sharey=True, layout='tight')
                ax_index = [0, 1]
            else:
                fig, ax = plt.subplots(1, 1, figsize=figsize, sharex=True,
                                       sharey=True, layout='tight')
                ax_index = [0]

        ax = [fig.axes[i] for i in ax_index]
        ## plot maps
        if option.lower() == 'both':
            fig = plot_2Dscalar(
                fig, ax[0], self.data[:, :, 0], self.base, levels1, chgline,
                isovalues, colormap, cbar_label1, a_range, b_range, rectangle,
                edgeplot, x_ticks, y_ticks, **kwargs
            )
            fig = plot_2Dscalar(
                fig, ax[1], self.data[:, :, 1], self.base, levels2, spinline,
                isovalues, colormap, cbar_label2, a_range, b_range, rectangle,
                edgeplot, x_ticks, y_ticks, **kwargs
            )
        elif option.lower() == 'charge':
            fig = plot_2Dscalar(
                fig, ax[0], self.data[:, :, 0], self.base, levels1, chgline,
                isovalues, colormap, cbar_label1, a_range, b_range, rectangle,
                edgeplot, x_ticks, y_ticks,  **kwargs
            )
        elif option.lower() == 'spin':
            fig = plot_2Dscalar(
                fig, ax[0], self.data[:, :, 1], self.base, levels2, spinline,
                isovalues, colormap, cbar_label2, a_range, b_range, rectangle,
                edgeplot, x_ticks, y_ticks,  **kwargs
            )
        else:
            raise ValueError("Unknown option: '{}'.".format(option))

        # title and labels
        for a in ax:
            if self.unit.lower() == 'angstrom':
                a.set_xlabel(r'$\AA$'); a.set_ylabel(r'$\AA$')
            else:
                a.set_xlabel('Bohr'); a.set_ylabel('Bohr')

        if np.all(title!=None):
            if title.lower() == 'default':
                if option.lower() == 'both':
                    ax[0].set_title('Charge Density')
                    ax[1].set_title('Spin Density')
                elif option.lower() == 'charge':
                    ax[0].set_title('Charge Density')
                else:
                    ax[0].set_title('Spin Density')
            else:
                for a in ax:
                    a.set_title(title)

        # restore old unit
        self._set_unit(uold)
        return fig

    def plot_3D(self,
                unit='Angstrom',
                option='charge',
                isovalue=None,
                volume_3d=False,
                contour_2d=False,
                interp='no interp',
                interp_size=1,
                grid_display_range=[[0,1], [0,1], [0,1]],
                show_the_scene=True,
                **kwargs):
        """
        Visualize **2D or 3D** charge densities with atomic structures using
        `MayaVi <https://docs.enthought.com/mayavi/mayavi/>`_ (*not installed
        by default*).

        * For 2D charge/spin densities, plot 2D heatmap with/without contour lines.
        * For 3D charge/spin densities, plot 3D isosurfaces or volumetric data.

        Args:
            unit (str): 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for Bohr:math:`^{-3}`.
            option (str): 'charge' or 'spin'.
            isovalue (float|array): Isovalues of 3D/2D contour plots. A number
                or an array for user-defined values of isosurfaces, **must be
                consistent with ``unit``**. By default half between max and min
                values.
            volume_3d (bool): *3D only*. Display 3D volumetric data instead of
                isosurfaces. ``isovalue`` is disabled.
            contour_2d (bool): *2D only* Display 2D black contour lines over
                colored contour surfaces.
            interp (str): Interpolate data to smoothen the plot. 'no interp' or
                'linear', 'nearest', 'slinear', 'cubic'. please refer to
                `scipy.interpolate.interpn <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interpn.html>`_
                 The interpolated data is not saved.
            interp_size (list[int]|int): The new size of interpolated data
                (list) or a scaling factor. *Valid only when ``interp`` is not
                'no interp'*.
            grid_display_range (array): 3\*2 array defining the displayed
                region of the data grid. Fractional coordinates a, b, c are
                used but only the periodic directions are applied.
            show_the_scene (bool): Display the scene by ``mlab.show()`` and
                return None. Otherwise return the scene object.
            \*\*kwargs: Optional keywords passed to MayaVi or ``CStructure.visualize()``.
                Allowed keywords are listed below.
            colormap (turple|str): Colormap of isosurfaces/heatmaps. Or a 1\*3
                RGB turple from 0 to 1 to define colors. *Not for volume_3d=True*.
            opacity (float): Opacity from 0 to 1. For ``volume_3d=True``, that
                defines the opacity of the maximum value. The opacity of the
                minimum is half of it.
            transparent (bool): Scalar-dependent opacity. *Not for volume_3d=True*.
            color (turple): Color of contour lines. *'contour_2d=True' only*.
            line_width (float): Width of 2D contour lines. *'contour_2d=True' only*.
            vmax (float): Maximum value of colormap.
            vmin (float): Minimum value of colormap.
            title (str): Colorbar title.
            orientation (str): Orientation of colorbar, 'horizontal' or 'vertical'.
            nb_labels (int): The number of labels to display on the colorbar.
            label_fmt (str): The string formater for the labels, e.g., '%.1f'.
            atom_color (str): Color map of atoms. 'jmol' or 'cpk'.
            bond_color (turple): Color of bonds, in a 1\*3 RGB turple from 0 to 1.
            atom_bond_ratio (str): 'balls', 'large', 'medium', 'small' or
                'sticks'. The relative sizes of balls and sticks.
            cell_display (bool): Display lattice boundaries (at \[0., 0., 0.\] only).
            cell_color (turple): Color of lattice boundaries, in a 1\*3 RGB turple from 0 to 1.
            cell_linewidth (float): Linewidth of plotted lattice boundaries.
            display_range (array): 3\*2 array defining the displayed region of
                the structure. Fractional coordinates a, b, c are used but only
                the periodic directions are applied.
            scale (float): See :ref:`geometry.CStructure.get_bonds() <ref-CStrucGetBonds>`.
            special_bonds (dict): See :ref:`geometry.CStructure.get_bonds() <ref-CStrucGetBonds>`.
            azimuth: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            elevation: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            distance: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
                By default set to 'auto'.
            focalpoint: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
            roll: See `mlab.view() <https://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html#view>`_.
        Returns:
            fig: MayaVi scence object, if ``show_the_scene=False``.
        """
        import numpy as np
        import warnings
        from CRYSTALpytools.base.plotbase import plot_3Dscalar, plot_3Dplane
        try:
            from mayavi import mlab
        except ModuleNotFoundError:
            raise ModuleNotFoundError('MayaVi is required for this functionality, which is not in the default dependency list of CRYSTALpytools.')


        # Input processing and sanity check
        if self.dimension != 2 and self.dimension != 3:
            raise Exception('Not a 3D/2D charge density object.')

        uold = self.unit
        if self.unit.lower() != unit.lower(): self._set_unit(unit)

        if np.all(self.structure==None):
            raise Exception("Geometry structure not available.")

        option = option.lower()
        if option == 'charge': ispin = 0
        elif option == 'spin': ispin = 1
        else: raise ValueError("Unknown option: '{}'.".format(option))

        if np.all(isovalue==None):
            isovalue = (np.max(self.data) - np.min(self.data)) * 0.5 + np.min(self.data)

        interp = interp.lower()
        if interp not in ['no interp', 'linear', 'nearest', 'slinear', 'cubic']:
            raise ValueError("Unknown interpolation method : '{}'.".format(interp))
        if volume_3d == True:
            interp = 'no interp' # not using GridInterpolate.


        # Expansion, check periodicity
        grid_display_range = np.array(grid_display_range, dtype=float)
        if self.dimension == 2 and len(grid_display_range) > 2:
            if grid_display_range[2, 0] != 0 or grid_display_range[2, 1] != 1:
                warnings.warn("For 2D data grid, a 2x2 display range should be defined.",
                              stacklevel=2)
            grid_display_range = grid_display_range[0:2]
        if len(grid_display_range) != self.dimension:
            raise Exception("Grid display range must have the same dimensionality as grid data.")

        if 'display_range' in kwargs.keys():
            display_range = np.array(kwargs['display_range'], dtype=float)
        else:
            display_range = np.array([[0,1], [0,1], [0,1]], dtype=float)

        if np.any(self.structure.pbc==False):
            idir = np.where(self.structure.pbc==False)[0]
            display_range[idir] = [0., 1.]

        idx = np.where(display_range[:,1]-display_range[:,0]<1e-4)[0]
        if len(idx) > 0:
            direct = ['x', 'y', 'z'][idx[0]]
            raise Exception("Structure display range error along {} axis!\n{} min = {:.2f}, {} max = {:.2f}. No data is displayed.".format(
                direct, direct, display_range[idx[0], 0], direct, display_range[idx[0], 1]))
        idx = np.where(grid_display_range[:,1]-grid_display_range[:,0]<1e-4)[0]
        if len(idx) > 0:
            direct = ['x', 'y', 'z'][idx[0]]
            raise Exception("Grid display range error along {} axis!\n{} min = {:.2f}, {} max = {:.2f}. No data is displayed.".format(
                direct, direct, grid_display_range[idx[0], 0], direct, grid_display_range[idx[0], 1]))

        # Plot structure
        keys = ['atom_color', 'bond_color', 'atom_bond_ratio', 'cell_display',
                'cell_color', 'cell_linewidth', 'scale', 'special_bonds']
        keywords = dict(show_the_scene=False, display_range=display_range)
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v

        fig = self.structure.visualize(**keywords)

        # Plot data
        keys = ['colormap', 'opacity', 'transparent', 'line_width', 'color',
                'vmax', 'vmin', 'title', 'orientation', 'nb_labels', 'label_fmt']
        if self.dimension == 3: # 3D isosurfaces
            keywords = dict(fig=fig,
                            base=self.base,
                            data=self.data[:, :, :, ispin],
                            isovalue=isovalue,
                            volume_3d=volume_3d,
                            interp=interp,
                            interp_size=interp_size,
                            display_range=grid_display_range)
            for k, v in zip(kwargs.keys(), kwargs.values()):
                if k in keys: keywords[k] = v

            fig = plot_3Dscalar(**keywords)
        else: # 2D contour plots
            keywords = dict(fig=fig,
                            base=self.base,
                            data=self.data[:, :, ispin],
                            levels=isovalue,
                            contour_2d=contour_2d,
                            interp=interp,
                            interp_size=interp_size,
                            display_range=grid_display_range)
            for k, v in zip(kwargs.keys(), kwargs.values()):
                if k in keys: keywords[k] = v

            fig = plot_3Dplane(**keywords)

        # Final setups
        keys = ['azimuth', 'elevation', 'distance', 'focalpoint', 'roll']
        keywords = dict(figure=fig, distance='auto')
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v
        mlab.view(**keywords)

        if show_the_scene == False:
            return fig
        else:
            mlab.gcf().scene.parallel_projection = True
            mlab.show()
            return

    def to_xsf(self, filename):
        """
        .. _ref-ChgDensToXSF:

        **2D / 3D data only**

        Write data into the `XCrySDen <http://www.xcrysden.org/>`_ XSF file
        format. The code writes into either 2D or 3D data grids depending on the
        ``dimension`` attribute. For spin-polarized cases, 2 maps will be
        written with the title 'alpha+beta' and 'alpha-beta' and into 2 data
        blocks. Objects updated by the ``alpha_beta()`` method will get grids
        with the **same** titles, which needs extra attention.

        .. note::

            Geometry information is mandatory.

        .. note::

            For 3D data grid, if CUBE file(s) are read without output, the XSF
            file output has the 3D 'CRYSTAL' periodicity. As far as the authors
            have been aware of, this only causes Segmentation Fault of XCrySDen
            1.6.2 when dealing with low dimensional systems. Available solution
            includes 1. including output file 2. using other software such
            as `VESTA <https://jp-minerals.org/vesta/en/>`_ 3. Changing the
            keyword manually.

        Args:
            filename (str): The output xsf filename.
        Returns:
            None
        """
        from CRYSTALpytools.io.xcrysden import XSF
        import numpy as np

        if np.all(self.structure==None): raise Exception('Geometry info unavailable.')

        # Definitions of base vector are different in 2D cases
        if self.dimension == 3:
            base = self.base
        elif self.dimension == 2:
            base = np.vstack([self.base[1], self.base[2], self.base[0]])
        else:
            raise Exception('Limited to 2D and 3D charge / spin densities.')

        if self.dimension == 2:
            xsf = XSF(self.structure, grid_base=base, grid_data=self.data[:, :, 0])
            xsf.write(filename, geomonly=False, grid_name='alpha+beta')
            if self.spin == 2:
                xsf = XSF(self.structure, grid_base=base, grid_data=self.data[:, :, 1])
                xsf._write_grid(filename, grid_name='alpha-beta')
        else:
            xsf = XSF(self.structure, grid_base=base, grid_data=self.data[:, :, :, 0])
            xsf.write(filename, geomonly=False, grid_name='alpha+beta')
            if self.spin == 2:
                xsf = XSF(self.structure, grid_base=base, grid_data=self.data[:, :, :, 1])
                xsf._write_grid(filename, grid_name='alpha-beta')
        return

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
            cst = au_to_angstrom(1.)
            self.unit = 'Angstrom'
        elif unit.lower() == 'a.u.':
            cst = angstrom_to_au(1.)
            self.unit = 'a.u.'
        else:
            raise ValueError('Unknown unit.')

        lprops = [] # length units
        dprops = ['data'] # density units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for d in dprops:
            newattr = getattr(self, d) / cst**3
            setattr(self, d, newattr)
        return self
