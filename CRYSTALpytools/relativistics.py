#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for '2c-SCF', i.e., 2-component SCF, and relativistics.
"""
import numpy as np
from CRYSTALpytools import units
from CRYSTALpytools.electronics import ChargeDensity as ChgDens

class ChargeDensity(ChgDens):
    """
    Same as ``electronics.ChargeDensity``, the spin density. But its dimension
    is kept to be commensurate with other keywords of the 'PROPS2COMP' block.

    .. note::

        Its ``type`` attribute is 'ECHG' rather than 'DENSITY'

    """


class VectorField():
    """
    The basic vector field object, containing a nY\*nX\*3 (nZ\*nY\*nX\*3) data
    array for 2D (3D) fields. **3D methods under development**.
    """
    def plot_2D(self, levels=100, quiverplot=True, quiverscale=1.0,
                colorplot=True, colormap='jet', cbar_label=None, a_range=[],
                b_range=[], rectangle=False, edgeplot=False, x_ticks=5,
                y_ticks=5, figsize=[6.4, 4.8], fig=None, ax_index=None, **kwargs):
        """
        Plot 2D vector field.

        3 styles are available:

        1. ``quiverplot=True`` and ``colorplot=True``: The color-filled contour
            illustrates the norm of vectors. The black arrows indicates both
            the directions and norms of in-plane prjections.  
        2. ``quiverplot=True`` and ``colorplot=False``: The arrows are colored
            to indicate the directions and norms of in-plane prjections.  
        3. ``quiverplot=False`` and ``colorplot=True``: The color-filled contour
            illustrates the norm of vectors, similar to the 2D scalar map.

        Args:
            levels (int|array): Set levels of colored contour/quiver plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D.
            quiverplot (bool): Plot 2D field of arrows.
            quiverscale (float): Tune the length of arrows. Useful only if
                ``quiverplot=True``.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if ``colorplot=True``.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in
                fractional coordinate.
            rectangle (bool): If :math:`a, b` are non-orthogonal, plot a
                rectangle region and reset :math:`b`. If used together with
                ``b_range``, that refers to the old :math:`b`.
            edgeplot (bool): Whether to add cell edges represented by the
                original base vectors (not inflenced by a/b range or rectangle
                options).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (list[int]): *Developer Only*, indices of axes in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``axes.quiver()`` function
                to set arrow styles.
        """
        from CRYSTALpytools.base.plotbase import plot_2Dscalar, plot_2Dvector
        import numpy as np
        import matplotlib.pyplot as plt
        import warnings

        # dimen
        if self.dimension != 2:
            raise Exception('Not a 2D charge density object.')

        # levels
        ## get norm
        vnorm = np.linalg.norm(self.data, axis=2)
        if isinstance(levels, int) or isinstance(levels, float):
            levels = np.linspace(np.min(vnorm), np.max(vnorm), int(levels))
        else:
            levels = np.array(levels, dtype=float)
            if levels.ndim > 1: raise ValueError('Levels must be a 1D array.')

        # plot
        if np.all(fig==None):
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            axes = fig.axes
        else:
            if np.all(ax_index==None):
                raise ValueError("Indices of axes must be set when 'fig' is not None.")
            axes = [fig.axes[i] for i in ax_index]

        for ax in axes:
            # plot colormap
            if colorplot == True:
                fig = plot_2Dscalar(
                    fig, ax, vnorm, self.base, levels, None, None, colormap,
                    cbar_label, a_range, b_range, rectangle, edgeplot, x_ticks,
                    y_ticks
                )
            # plot quiver
            if quiverplot == True:
                if colorplot == True:
                    cxlim = ax.get_xlim()
                    cylim = ax.get_ylim() # in case of white edges
                    fig = plot_2Dvector(
                        fig, ax, self.data, self.base, quiverscale, 'k', levels,
                        colormap, cbar_label, a_range, b_range, rectangle, False,
                        x_ticks, y_ticks
                    )
                    ax.set_xlim(cxlim)
                    ax.set_ylim(cylim)
                else:
                    fig = plot_2Dvector(
                        fig, ax, self.data, self.base, quiverscale, 'colored', levels,
                        colormap, cbar_label, a_range, b_range, rectangle, edgeplot,
                        x_ticks, y_ticks
                    )
        return fig


class Magnetization(VectorField):
    """
    The class for magnetization. Unit: 'SI' (length: :math:`\\AA`,
    magnetization: A/m).

    Args:
        data (array): nY\*nX\*3 (nZ\*nY\*nX\*3) array of magnetization vectors
            in 2D (3D) plane.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC (2D) or 3 base vectors (3D)
        dimen (int): Dimensionality of the plot.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'SI' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc, unit='SI'):
        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'MAGNETIZ'

    def plot_2D(self, unit='SI', levels=100, quiverplot=True, quiverscale=1.0,
                colorplot=True, colormap='jet', cbar_label='default', a_range=[],
                b_range=[], rectangle=False, edgeplot=False, x_ticks=5,
                y_ticks=5, title='default', figsize=[6.4, 4.8],
                fig=None, ax_index=None, **kwargs):
        """
        Plot 2D magnetization field.

        3 styles are available:

        1. ``quiverplot=True`` and ``colorplot=True``: The color-filled contour
            illustrates the norm of vectors. The black arrows indicates both
            the directions and norms of in-plane prjections.  
        2. ``quiverplot=True`` and ``colorplot=False``: The arrows are colored
            to indicate the directions and norms of in-plane prjections.  
        3. ``quiverplot=False`` and ``colorplot=True``: The color-filled contour
            illustrates the norm of vectors, similar to the 2D scalar map.

        Args:
            unit (str): Plot unit. 'SI' for :math:`\\AA` and A/m. 'a.u.' for
                Bohr and a.u. magnetization.
            levels (int|array): Set levels of colored contour/quiver plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D.
            quiverplot (bool): Plot 2D field of arrows.
            quiverscale (float): Tune the length of arrows. Useful only if
                ``quiverplot=True``.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if ``colorplot=True``.
            cbar_label (str|None): Label of colorbar. 'default' for unit.
                'None' for no label. Useful only if ``colorplot=True``.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in
                fractional coordinate.
            rectangle (bool): If :math:`a, b` are non-orthogonal, plot a
                rectangle region and reset :math:`b`. If used together with
                ``b_range``, that refers to the old :math:`b`.
            edgeplot (bool): Whether to add cell edges represented by the
                original base vectors (not inflenced by a/b range or rectangle
                options).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for proeprty plotted.
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (list[int]): *Developer Only*, indices of axes in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``axes.quiver()`` function
                to set arrow styles.
        """
        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # cbar_label
        if isinstance(cbar_label, str):
            if cbar_label.lower() == 'default':
                if self.unit.lower() == 'si':
                    cbar_label = 'Unit: A/m'
                else:
                    cbar_label = 'Unit: a.u.'
        else:
            cbar_label = None
        # figure
        if np.all(fig!=None):
            fig = super().plot_2D(
                levels, quiverplot, quiverscale, colorplot, colormap,
                cbar_label, a_range, b_range, rectangle, edgeplot, x_ticks,
                y_ticks, figsize, fig, ax_index, **kwargs
            )
        else:
            fig = super().plot_2D(
                levels, quiverplot, quiverscale, colorplot, colormap,
                cbar_label, a_range, b_range, rectangle, edgeplot, x_ticks,
                y_ticks, figsize, **kwargs
            )
            ax_index = [0]
        # title
        if isinstance(title, str):
            if title.lower() == 'default':
                title = 'Magnetization'
            for iax in ax_index:
                fig.axes[iax].set_title(title)

        # restore old unit
        self._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``Magnetization`` object.

        Args:
            unit (str): 'SI', length: :math:`\\AA`, magnetization: A/m.
                'a.u.', all in a.u..
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, ampere_to_au, au_to_ampere

        if unit.lower() == self.unit.lower():
            return self

        if unit.lower() == 'si':
            self.unit = 'SI'
            mcst = au_to_ampere(1.) * angstrom_to_au(1.)*1e10
            lcst = au_to_angstrom(1.)
        elif unit.lower() == 'a.u.':
            mcst = ampere_to_au(1.) * au_to_angstrom(1.)*1e10
            lcst = angstrom_to_au(1.)
        else:
            raise ValueError('Unknown unit.')

        lprops = ['base'] # length units
        mprops = ['data'] # magnetization units
        for l in lprops:
            newattr = getattr(self, l) * lcst
            setattr(self, l, newattr)
        for m in mprops:
            newattr = getattr(self, m) * mcst
            setattr(self, m, newattr)
        return self


class OrbitalCurrentDensity(VectorField):
    """
    The class for orbital current density. Unit: 'SI' (length: :math:`\\AA`,
    orbital current density: A/m:math:`^{2}`).

    Args:
        data (array): nY\*nX\*3 (nZ\*nY\*nX\*3) array of orbnital current vectors
            in 2D (3D) plane.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC (2D) or 3 base vectors (3D)
        dimen (int): Dimensionality of the plot.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'SI' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc, unit='SI'):
        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'ORBCURDENS'

    def _set_unit(self, unit):
        """
        Set units of data of ``OrbitalCurrentDensity`` object.

        Args:
            unit (str): 'SI', length: :math:`\\AA`, orbital current density:
                A/m:math:`^{2}`. 'a.u.', all in a.u..
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, ampere_to_au, au_to_ampere

        if unit.lower() == self.unit.lower():
            return self

        if unit.lower() == 'si':
            self.unit = 'SI'
            mcst = au_to_ampere(1.) * (angstrom_to_au(1.)*1e10)**2
            lcst = au_to_angstrom(1.)
        elif unit.lower() == 'a.u.':
            mcst = ampere_to_au(1.) * (au_to_angstrom(1.)*1e10)**2
            lcst = angstrom_to_au(1.)
        else:
            raise ValueError('Unknown unit.')

        lprops = ['base'] # length units
        mprops = ['data'] # current density units
        for l in lprops:
            newattr = getattr(self, l) * lcst
            setattr(self, l, newattr)
        for m in mprops:
            newattr = getattr(self, m) * mcst
            setattr(self, m, newattr)
        return self


class SpinCurrentDensity(VectorField):
    """
    The class for spin current density. Unit: 'SI' (length: :math:`\\AA`,
    spin current density: A/m:math:`^{2}`).

    Args:
        data_x (array): nY\*nX\*3 (nZ\*nY\*nX\*3) array of spin current vectors
            :math:`J^{x}` in 2D (3D) plane.
        data_y (array): nY\*nX\*3 (nZ\*nY\*nX\*3) array of spin current vectors
            :math:`J^{y}` in 2D (3D) plane.
        data_z (array): nY\*nX\*3 (nZ\*nY\*nX\*3) array of spin current vectors
            :math:`J^{z}` in 2D (3D) plane.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC (2D) or 3 base vectors (3D)
        dimen (int): Dimensionality of the plot.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'SI' (case insensitive).
    """
    def __init__(self, data_x, data_y, data_z, base, dimen, struc, unit='SI'):
        self.data_x = np.array(data_x, dtype=float)
        self.data_y = np.array(data_y, dtype=float)
        self.data_z = np.array(data_z, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SPICURDENS'

    def _set_unit(self, unit):
        """
        Set units of data of ``SpinCurrentDensity`` object.

        Args:
            unit (str): 'SI', length: :math:`\\AA`, spin current density:
                A/m:math:`^{2}`. 'a.u.', all in a.u..
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, ampere_to_au, au_to_ampere

        if unit.lower() == self.unit.lower():
            return self

        if unit.lower() == 'si':
            self.unit = 'SI'
            mcst = au_to_ampere(1.) * (angstrom_to_au(1.)*1e10)**2
            lcst = au_to_angstrom(1.)
        elif unit.lower() == 'a.u.':
            mcst = ampere_to_au(1.) * (au_to_angstrom(1.)*1e10)**2
            lcst = angstrom_to_au(1.)
        else:
            raise ValueError('Unknown unit.')

        lprops = ['base'] # length units
        mprops = ['data_x', 'data_y', 'data_z'] # current density units
        for l in lprops:
            newattr = getattr(self, l) * lcst
            setattr(self, l, newattr)
        for m in mprops:
            newattr = getattr(self, m) * mcst
            setattr(self, m, newattr)
        return self
