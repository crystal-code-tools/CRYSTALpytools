#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for `TOPOND <https://www.crystal.unito.it/topond.html>`_ topological
analysis of electron density
"""
from CRYSTALpytools import units
from CRYSTALpytools.electronics import ElectronBandDOS
import numpy as np

class Surf(ChargeDensity):
    """
    TOPOND 2D scalar countour plot class. Length unit: :math:`\\AA`;
    Charge unit: e; Energy / field unit: eV.

    Args:
        surfdata (array): 2D plot data.
        surfbase (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        type (list|str): Type of data. Allowed entries are listed above.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """

    @classmethod
    def from_file(cls, file, type='infer', output=None):
        """
        Instantite the class by 2D scalar countour plot files ('SURF*.DAT').


        .. note::

            For the convenience of analysis and plotting, it is important to select
            the correct type for your input file. By default `type='infer'` will
            search for (case insensitive) the following strings:

            ``'SURFRHOO','SURFSPDE','SURFLAPP','SURFLAPM','SURFGRHO','SURFKKIN','SURFGKIN','SURFVIRI','SURFELFB'``

            For their meanings, please refer the `TOPOND manual <https://www.crystal.unito.it/include/manuals/topond.pdf>`_.

        Args:
            file (str): TOPOND formatted 2D plot file
            type (str): 'infer' or specified. Otherwise warning will be given.
            output (str): Standard output of Properties calculation, used to
                get geometry.

        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond2D(file, type)

    def plot(self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
             isovalues=True, colorplot=False, colormap='jet', cbar_label=None,
             x_range=[], y_range=[], x_ticks=7, y_ticks=7, add_title=True,
             figsize=[6.4, 4.8], **kwargs):
        """
        Plot 2D contour lines, color maps or both for the 2D data set.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA` - eV. 'a.u.' for
                Bohr - Hartree.
            levels (array): Set levels of contour plot. 'Default' for built-in,
                property adaptive levels (``unit='Angstrom'``). Otherwise
                entries **must be consistent with ``unit``**.
            lineplot (bool): Plot contour lines.
            linewidth (float): Contour linewidth. Useful only if
                ``lineplot=True``. Other properties are not editable.
            isovalues (bool): Add isovalues to contour lines. Useful only if
                ``lineplot=True``. 2 decimal places is used for 'SURFELFB' and
                3 for others. 4 for difference maps.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``.
            xrange (list): 1\*2 list of x axis range. **Must be consistent with
                ``unit``**.
            yrange (list): 1\*2 list of y axis range. **Must be consistent with
                ``unit``**.
            xticks (int): Number of ticks on x axis.
            yticks (int): Number of ticks on y axis.
            add_title (bool): Whether to add property plotted as title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            \*\*kwargs : Other arguments passed to ``axes.contour()`` function
                to set contour lines.

        Returns:
            fig (Figure): Matplotlib Figure object
            ax (Axes): Matplotlib Axes object
        """
        from CRYSTALpytools.base.plotbase import plot_2Dscalar
        import numpy as np
        import matplotlib.pyplot as plt
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'

        if self.unit.lower() != unit.lower():
            self._set_unit(unit)

        # level
        if np.all(levels=='default'):
            if self.type == 'SURFELFB':
                levels = np.linspace(0, 1, 21)
            elif self.type in ['SURFLAPP', 'SURFLAPM', 'SURFVIRI', 'SURFKKIN']:
                levels = np.array([-8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8], dtype=float)
            elif self.type in ['SURFRHOO', 'SURFGRHO', 'SURFGKIN']:
                levels = np.array([0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20], dtype=float)
            else:
                warnings.warn("Unknown data type: {}. A linear scale is used.".format(self.type))
                levels = np.linspace(np.min(self.surfdata), np.max(self.surfdata), 10)
        else:
            levels = np.array(levels, dtype=float)
        # contour line styles
        if lineplot == True and colorplot == False:
            contourline = []
            for i in levels:
                if i < -1e-6:
                    contourline.append(['b', '--', linewidth])
                elif i > 1e-6:
                    contourline.append(['r', '-', linewidth])
                else:
                    contourline.append(['k', 'dotted', linewidth])
        elif lineplot == True and colorplot == True:
            contourline = []
            for i in levels:
                if i < -1e-6:
                    contourline.append(['k', '--', linewidth])
                elif i > 1e-6:
                    contourline.append(['k', '-', linewidth])
                else:
                    contourline.append(['k', 'dotted', linewidth])
        else:
            contourline = None
        # isovalues
        if isovalues == True:
            if self.type == 'SURFELFB':
                isovalue = '%.2f'
            else:
                isovalue = '%.3f'
        else:
            isovalue = None
        # colormap
        if colorplot == False:
            colormap = None

        # plot
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = plot_2Dscalar(ax, self.surfdata, self.surfbase, contourline, isovalue,
                           colormap, x_ticks, y_ticks, cbar_label, **kwargs)
        if len(x_range) != 0:
            ax.set_xlim(xrange)
        if len(y_range) != 0:
            ax.set_ylim(yrange)
        if unit.lower() == 'angstrom':
            ax.set_xlabel(r'$\AA$')
            ax.set_ylabel(r'$\AA$')
        else:
            ax.set_xlabel(r'$Bohr$')
            ax.set_ylabel(r'$Bohr$')
        if add_title == True:
            ax.set_title(self.type)

        self._set_unit(uold)
        return fig, fig.axes

    def _set_unit(self, unit):
        """
        Set units of data of ``topond.Surf`` object. 2 sets of units allowed:

        1. Energy / field: 'eV', length: ':math:`\\AA`', Charge: 'e'.  
        2. Energy / field: 'Hartree', length: 'Bohr', Charge: 'e'.

        Densities are in :math:`len^{-3}`, gradients are in :math:`len^{-4}`
        and Laplacians are in :math:`len^{5}` etc..

        Args:
            unit (str): 'Angstrom' or 'a.u.'
        """
        import warnings
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, H_to_eV, eV_to_H

        if unit.lower() == self.unit.lower():
            return self

        density = ['SURFRHOO', 'SURFSPDE']
        gradient = ['SURFGRHO']
        laplacian = ['SURFLAPP', 'SURFLAPM']
        e_density = ['SURFKKIN', 'SURFGKIN', 'SURFVIRI']
        normalized = ['SURFELFB']

        if unit.lower() == 'angstrom':
            cst = au_to_angstrom(1.)
            ecst = H_to_eV(1.)
            self.unit = 'Angstrom'
        elif unit.lower() == 'a.u.':
            cst = angstrom_to_au(1.)
            ecst = eV_to_H(1.)
            self.unit = 'a.u.'
        else:
            raise ValueError('Unknown unit.')

        self.surfbase = self.surfbase * cst

        if self.type.upper() in density: # Bohr^-3 <---> AA^-3
            self.surfdata = self.surfdata / cst**3
        elif self.type.upper() in gradient: # Bohr^-4 <---> AA^-4
            if unit.lower() == 'angstrom':
            self.surfdata = self.surfdata / cst**4
        elif self.type.upper() in laplacian: # Bohr^-5 <---> AA^-5
            self.surfdata = self.surfdata / cst**5
        elif self.type.upper() in e_density: # Eh.Bohr^-3 <---> eV.AA^-3
            self.surfdata = self.surfdata / cst**5 * ecst
        elif self.type.upper() in normalized:
            pass
        else:
            warnings.warn('Unknown data type. Using default units.', stacklevel=2)

        return self

