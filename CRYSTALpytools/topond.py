#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for `TOPOND <https://www.crystal.unito.it/topond.html>`_ topological
analysis of electron density
"""
from CRYSTALpytools import units
from CRYSTALpytools.electronics import ChargeDensity
import numpy as np

class Surf(ChargeDensity):
    """
    TOPOND 2D scalar countour plot class. Length unit: :math:`\\AA`;
    Charge unit: e; Energy / field unit: eV.

    Args:
        data (array): 2D Plot data. nX\*nY
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        type (str): See the classmethod ``from_file``.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, struc=None, type='unknown', unit='Angstrom'):
        super().__init__(data, base, 1, 2, struc, unit)
        self.type = type

    @classmethod
    def from_file(cls, file, output=None, type='infer'):
        """
        Instantite the class by 2D scalar countour plot files ('SURF*.DAT').


        .. note::

            For the convenience of analysis and plotting, it is important to select
            the correct type for your input file. By default `type='infer'` will
            search for (case insensitive) the following strings in the filename:

            'SURFRHOO', 'SURFSPDE', 'SURFLAPP', 'SURFLAPM', 'SURFGRHO',
            'SURFKKIN', 'SURFGKIN', 'SURFVIRI', 'SURFELFB'

            For their meanings, please refer the `TOPOND manual <https://www.crystal.unito.it/include/manuals/topond.pdf>`_.

        Args:
            file (str): TOPOND formatted 2D plot file
            output (str): Standard output of Properties calculation, used to
                get geometry.
            type (str): 'infer' or specified. Otherwise warning will be given.

        Returns:
            cls (Surf)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, type)

    def plot(self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
             isovalues='%.4f', colorplot=False, colormap='jet', cbar_label=None,
             a_range=[], b_range=[], x_ticks=5, y_ticks=5, cellplot=False,
             add_title=True, figsize=[6.4, 4.8], **kwargs):
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
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for default.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in
                fractional coordinate.
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
            fig (Figure): Matplotlib Figure object
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

        # levels
        if np.all(levels=='default'):
            if self.type == 'SURFELFB':
                levels = np.linspace(0, 1, 21)
            elif self.type in ['SURFLAPP', 'SURFLAPM', 'SURFVIRI', 'SURFKKIN']:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, 0, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
            elif self.type in ['SURFRHOO', 'SURFGRHO', 'SURFGKIN']:
                levels = np.array([0.02, 0.04, 0.08, 0.2, 0.4, 0.8,
                                   2, 4, 8, 20, 40, 80, 200], dtype=float)
            elif self.type == 'diff' or self.type == 'SURFSPDE': # difference
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
            else:
                warnings.warn("Unknown data type: {}. A linear scale is used.".format(self.type))
                levels = np.linspace(np.min(self.data), np.max(self.data), 10)
        else:
            if isinstance(levels, int) or isinstance(levels, float):
                levels = np.linspace(np.min(self.data), np.max(self.data), levels)
            else:
                levels = np.array(levels, dtype=float)
        # contour line styles
        if self.type != 'SURFELFB':
            blimit = -1e-6
            rlimit = 1e-6
        else:
            blimit = 0.5
            rlimit = 0.5
        if lineplot == True and colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit:
                    contourline.append(['b', '--', linewidth])
                elif i > rlimit:
                    contourline.append(['r', '-', linewidth])
                else:
                    contourline.append(['k', '-', linewidth*2])
        elif lineplot == True and colorplot == True:
            contourline = []
            for i in levels:
                if i < blimit:
                    contourline.append(['k', '--', linewidth])
                elif i > rlimit:
                    contourline.append(['k', '-', linewidth])
                else:
                    contourline.append(['k', '-', linewidth*2])
        else:
            contourline = None
        # colormap
        if colorplot == False:
            colormap = None
        # cbar_label
        if np.all(cbar_label==None):
            if unit.lower() == 'angstrom':
                if self.type in ['SURFRHOO', 'SURFSPDE']: cbar_label=r'$\rho$ ($|e|/\AA^{-3}$)'
                elif self.type in ['SURFGRHO']: cbar_label=r'$\nabla\rho$ ($|e|/\AA^{-4}$)'
                elif self.type in ['SURFLAPP', 'SURFLAPM']: cbar_label=r'$\nabla^{2}\rho$ ($|e|/\AA^{-5}$)'
                elif self.type in ['SURFKKIN', 'SURFGKIN', 'SURFVIRI']: cbar_label=r'$\rho$ ($eV/\AA^{-3}$)'
                else: cbar_label=self.type
            else:
                if self.type in ['SURFRHOO', 'SURFSPDE']: cbar_label=r'$\rho$ ($|e|/Bohr^{-3}$)'
                elif self.type in ['SURFGRHO']: cbar_label=r'$\nabla\rho$ ($|e|/Bohr^{-4}$)'
                elif self.type in ['SURFLAPP', 'SURFLAPM']: cbar_label=r'$\nabla^{2}\rho$ ($|e|/Bohr^{-5}$)'
                elif self.type in ['SURFKKIN', 'SURFGKIN', 'SURFVIRI']: cbar_label=r'$\rho$ ($Hartree/Bohr^{-3}$)'
                else: cbar_label=self.type
        # plot
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        fig, ax = plot_2Dscalar(
            fig, ax, self.data, self.base, levels, contourline, isovalues, colormap,
            cbar_label, a_range, b_range, False, cellplot, x_ticks, y_ticks, **kwargs
        )
        if unit.lower() == 'angstrom':
            ax.set_xlabel(r'$\AA$')
            ax.set_ylabel(r'$\AA$')
        else:
            ax.set_xlabel(r'$Bohr$')
            ax.set_ylabel(r'$Bohr$')
        if add_title == True:
            ax.set_title(self.type)

        self._set_unit(uold)
        return fig

    def substract(self, *args, type='infer'):
        """
        Substracting data of the same type from the object.

        Args:
            \*args (str|Surf): File names or ``Surf`` objects. Must be of the
                same type (check the attribute ``type``).
            type (str): 'infer' or specified. Otherwise warning will be given.
                Useful only when filenames are given. Check the classmethod
                ``from_file()``.
        Returns:
            self (Surf) : Data difference
        """
        from CRYSTALpytools.crystal_io import Properties_output

        objs = []
        for i in args:
            if isinstance(i, str):
                objs.append(Properties_output().read_topond2D(i, type=type))
            elif isinstance(i, Surf):
                objs.append(i)
            else:
                raise TypeError('Inputs must be file name strings or Surf objects.')
        self = super().substract(*objs)
        self.type = 'diff'
        return self

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

        self.base = self.base * cst

        if self.type.upper() in density: # Bohr^-3 <---> AA^-3
            self.data = self.data / cst**3
        elif self.type.upper() in gradient: # Bohr^-4 <---> AA^-4
            self.data = self.data / cst**4
        elif self.type.upper() in laplacian: # Bohr^-5 <---> AA^-5
            self.data = self.data / cst**5
        elif self.type.upper() in e_density: # Eh.Bohr^-3 <---> eV.AA^-3
            self.data = self.data / cst**5 * ecst
        elif self.type.upper() in normalized:
            pass
        else:
            warnings.warn('Unknown data type. Using default units.', stacklevel=2)

        return self


class Traj():
    """
    TOPOND trajectory plot class. Length unit: :math:`\\AA`.

    Args:
        wtraj (list[int]): 1\*nPath, weight of the path, int 0 to 3.
        traj (list[array]): 1\*nPath, list of critical paths. Every array is a
            nPoint\*3 3D ref framework coordinates of points on the path.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        type (str): See the classmethod ``from_file``.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).

    Returns:
        self (Traj): Import attributes: ``trajectory``, nPath\*2 list of
            trajectory and its weight; ``cpoint``, nCpoint\*3 array of critical
            point coordinates in 3D ref framework.
    """
    def __init__(self, wtraj, traj, base, struc, type='unknown', unit='Angstrom'):
        if len(wtraj) != len(traj):
            raise ValueError('Inconsistent lengths of input trajectory and its weight.')

        self.trajectory = [[int(wtraj[i]), np.array(traj[i], dtype=float)]
                           for i in range(len(wtraj))]
        self.base = np.array(base, dtype=float)
        self.structure = struc
        self.type = type
        self.unit = unit
        cpt = []
        for i in self.trajectory:
            if len(i[1]) == 1:
                cpt.append(i[1][0])
        self.cpoint = np.array(cpt)

    @classmethod
    def from_file(cls, file, output, type='infer'):
        """
        Generate a ``topond.Traj`` object from 'TRAJ*.DAT' file and standard
        output of TOPOND.

        .. note::

            For the convenience of analysis and plotting, it is important to select
            the correct type for your input file. By default `type='infer'` will
            search for (case insensitive) the following strings in the filename:

            'TRAJGRAD', 'TRAJMOLG',

            For their meanings, please refer the `TOPOND manual <https://www.crystal.unito.it/include/manuals/topond.pdf>`_.

        Args:
            file (str): TOPOND formatted 2D plot file
            output (str): Standard output of Properties calculation, used to
                get geometry.
            type (str): 'infer' or specified. Otherwise warning will be given.

        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, type=type)

    def plot_2D(self, unit='Angstrom', cpt_marker='o', cpt_c='k', cpt_s=10,
                traj_color='r', traj_linestyle=':', traj_linewidth=0.5,
                x_ticks=5, y_ticks=5, cellplot=False, add_title=True,
                figsize=[6.4, 4.8], overlay_surf=None, **kwargs):
        """
        Get TOPOND trajectory in a 2D plot.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import copy
        from CRYSTALpytools.plotbase import _get_operation

        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)

        # Get bottom surf figure first
        if np.all(overlay_surf!=None) and isinstance(overlay_surf, Surf):
            overlay_surf._set_unit(unit)
            diff_base = np.abs(overlay_surf.base-self.base)
            if np.any(diff_base>1e-3):
                raise Exception("The plotting base of overlayed surface and 2D trajectory are different.")

            kwargs['unit'] = unit; kwargs['figsize'] = figsize
            kwargs['a_range'] = []; kwargs['b_range'] = [] # no periodicity
            kwargs['x_ticks'] = x_ticks; kwargs['y_ticks'] = y_ticks
            kwargs['cellplot'] = cellplot; kwargs['add_title'] = add_title
            fig = overlay_surf.plot(**kwargs)
            ax = fig.axes[0]
        else:
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        # rotate the trajectory to plotting plane
        rot, disp = _get_operation(self.base)
        # plot TRAJ
        baserot = rot.apply(self.base)
        xmx = np.linalg.norm(baserot[2, :]-baserot[1, :])
        ymx = np.linalg.norm(baserot[0, :]-baserot[1, :])
        extra_width = {1 : 0., 2 : 0., 3 : 0.5} # extra linewidth for critical path
        for wt, traj in self.trajectory:
            traj = rot.apply(traj)
            # plot CPT
            if len(traj) == 1:
                v = traj[0] - baserot[1]
                if v[0]>=0 and v[0]<xmx and v[1]>=0 and v[1]<ymx and np.abs(v[2])<=1e-3:
                    ax.scatter(v[0], v[1], marker=cpt_marker, c=cpt_c, s=cpt_s)
            # plot TRAJ
            else:
                plttraj = []
                for v in traj:
                    v = v - baserot[1]
                    if v[0]>=0 and v[0]<xmx and v[1]>=0 and v[1]<ymx and np.abs(v[2])<=1e-3:
                        plttraj.append(v)

                if len(plttraj) == 0:
                    continue
                plttraj = np.array(plttraj)
                if wt != 0:
                    ax.plot(plttraj[:, 0], plttraj[:, 1], color=cpt_c,
                            linestyle='-', linewidth=traj_linewidth+extra_width[wt])
                else:
                    ax.plot(plttraj[:, 0], plttraj[:, 1], color=traj_color,
                            linestyle=traj_linestyle, linewidth=traj_linewidth)

        ax.set_aspect(1.0)
        ax.set_xlim(0, xmx)
        ax.set_ylim(0, ymx)
        self._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``topond.Traj`` object. :math:`\\AA`' or 'Bohr'.

        Args:
            unit (str): 'Angstrom' or 'a.u.'
        """
        import warnings
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom

        if unit.lower() == self.unit.lower():
            return self

        path = ['TRAJGRAD', 'TRAJMOLG']

        if unit.lower() == 'angstrom':
            cst = au_to_angstrom(1.)
            self.unit = 'Angstrom'
        elif unit.lower() == 'a.u.':
            cst = angstrom_to_au(1.)
            self.unit = 'a.u.'
        else:
            raise ValueError('Unknown unit.')

        self.base = self.base * cst
        if self.type.upper() in path: # Bohr <---> AA
            for i in range(len(self.trajectory)):
                self.trajectory[i][1] = self.trajectory[i][1] * cst
            self.cpoint = self.cpoint * cst
        elif self.type.lower() == 'unknown':
            warnings.warn('Unknown data type. Using default units.', stacklevel=2)

        return self
