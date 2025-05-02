#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for `TOPOND <https://www.crystal.unito.it/topond.html>`_ topological
analysis of electron density
"""
from CRYSTALpytools import units
import numpy as np


class ScalarField():
    """
    Basic TOPOND scalar field class, containing a nY\*nX (nZ\*nY\*nX) data
    array for 2D (3D) fields. Call the property-specific child classes below to
    use. **3D methods under development**.
    """
    def plot_2D(
        self, levels=100, lineplot=False, contourline=None, isovalues='%.2f',
        colorplot=False, colormap='jet', cbar_label=None,
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            2D periodicity (``a_range`` and ``b_range``), though available for
            the ``ScalarField`` class, is not suggested as TOPOND plotting
            window does not always commensurate with periodic boundary. The
            ``Trajectory`` class has no 2D periodicity so if ``overlay`` is not
            None, ``a_range``, ``b_range`` and ``edgeplot`` will be disabled.

        3 styles are available:

        1. ``lineplot=True`` and ``colorplot=True``: The color-filled contour
            map with black contour lines. Dotted lines for negative values and
            solid lines for positive values. The solid line twice in width for 0.  
        2. ``lineplot=False`` and ``colorplot=True``: The color-filled contour
            map.  
        3. ``lineplot=True`` and ``colorplot=False``: The color coded contour
            line map. Blue dotted line for negative values and red solid lines
            for positive values. The balck solid line twice in width for 0.

        Args:
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D.
            lineplot (bool): Plot contour lines.
            contourline (list): nLevel\*3 contour line styles. Useful only if
                ``lineplot=True``. For every line, color, line style and line
                width are defined in sequence.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        from CRYSTALpytools.base.plotbase import plot_2Dscalar
        import numpy as np
        import matplotlib.pyplot as plt
        import warnings

        # dimen
        if self.dimension != 2:
            raise Exception('Not a 2D scalar field object.')

        # levels
        levels = np.array(levels, ndmin=1, dtype=float)
        if levels.shape[0] == 1:
            levels = np.linspace(np.min(self.data), np.max(self.data), int(levels[0]))
        if levels.ndim > 1: raise ValueError('Levels must be a 1D array.')

        # contour line styles
        if lineplot == False:
            contourline = None
        # colormap
        if colorplot == False:
            colormap = None

        # overlay
        if np.all(overlay!=None) and isinstance(overlay, Trajectory):
            diff_base = np.abs(overlay.base-self.base)
            if np.any(diff_base>1e-3):
                raise Exception("The plotting base of surface and trajectory are different.")
            a_range = [0., 1.]; b_range=[0., 1.] # no periodicity for Traj

        # plot
        ## layout
        if np.all(fig==None):
            fig, ax = plt.subplots(1, 1, figsize=figsize, layout='tight')
            ax_index = 0
        else:
            if np.all(ax_index==None):
                raise ValueError("Indices of axes must be set when 'fig' is not None.")
            ax_index = int(ax_index)
            ax = fig.axes[ax_index]
        ## surf first
        fig = plot_2Dscalar(
            fig, ax, self.data, self.base, levels, contourline, isovalues, colormap,
            cbar_label, a_range, b_range, False, edgeplot, x_ticks, y_ticks
        )
        ## plot traj
        if np.all(overlay!=None) and isinstance(overlay, Trajectory):
            kwargs['fig'] = fig; kwargs['ax_index'] = ax_index
            kwargs['y_ticks'] = y_ticks; kwargs['x_ticks'] = x_ticks
            kwargs['figsize'] = figsize; kwargs['unit'] = overlay.unit # convert unit in wrappers!!!
            fig = overlay.plot_2D(**kwargs)

        return fig

    def subtract(self, *args):
        """
        Subtracting data of the same type from the object.

        Args:
            \*args (str|ScalarField): File names or ``ScalarField`` objects.
                Must be of the same type (check the attribute ``type``).
        Returns:
            self (ScalarField) : Data difference
        """
        from CRYSTALpytools.crystal_io import Properties_output

        for i in args:
            if isinstance(i, str):
                obj = Properties_output().read_topond(i, type=self.type)
            elif isinstance(i, ScalarField):
                obj = i
            else:
                raise TypeError('Inputs must be file name strings or Surf objects.')

            # type
            if self.type != obj.type:
                raise TypeError('Input is not the same type as object.')
            # base vector
            if not np.all(np.abs(self.base-obj.base)<1e-6):
                raise ValueError('Inconsistent base vectors between input and object.')
            # dimensionality
            if self.dimension != obj.dimension:
                raise ValueError('Inconsistent dimensionality between input and object.')
            # mesh grid
            if self.data.shape != obj.data.shape:
                raise ValueError('Inconsistent mesh grid between input and object.')
            # subtract
            self.data = self.data - obj.data

        self.subtracted = True # Hidden. For plotting.
        return self

    def substract(self, *args):
        """An old typo"""
        return self.subtract(*args)


class Trajectory():
    """
    Basic TOPOND trajectory plot class. Call the property-specific child classes
    below to use.
    """
    def plot_2D(
        self, cpt_marker='o', cpt_color='k', cpt_size=10,
        traj_color='r', traj_linestyle=':', traj_linewidth=0.5, x_ticks=5,
        y_ticks=5, figsize=[6.4, 4.8], overlay=None, fig=None, ax_index=None,
        **kwargs):
        """
        Get TOPOND trajectory in a 2D plot.

        .. note::

            2D periodicity (``a_range`` and ``b_range``) is not available for
            the ``Trajectory`` class. If ``overlay`` is not None, ``a_range``
            and ``b_range`` and ``edgeplot`` will be disabled for the
            ``ScalarField`` object.

        Args:
            cpt_marker (str): Marker of critical point scatter.
            cpt_color (str): Marker color of critical point scatter.
            cpt_size (float|int): Marker size of critical point scatter.
            traj_color (str): Line color of 2D trajectory plot.
            traj_linestyl (str): Line style of 2D trajectory plot.
            traj_linewidth (str): Line width of 2D trajectory plot.
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|ScalarField): Overlapping a 2D plot from the
                ``topond.ScalarField`` object if not None.
            fig (Figure|None): *Developers only*, matplotlib Figure class..
            ax_index (list[int]): *Developer Only*, indices of axes in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``ScalarField.plot_2D()``.
        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import copy
        from CRYSTALpytools.base.plotbase import GridRotation2D

        # overlay
        if np.all(overlay!=None) and isinstance(overlay, ScalarField):
            diff_base = np.abs(overlay.base-self.base)
            if np.any(diff_base>1e-3):
                raise Exception("The plotting base of surface and trajectory are different.")

        # plot
        ## layout
        if np.all(fig==None):
            fig, ax = plt.subplots(1, 1, figsize=figsize, layout='tight')
            ax_index = 0
        else:
            if np.all(ax_index==None):
                raise ValueError("Indices of axes must be set when 'fig' is not None.")
            ax_index = int(ax_index)
            ax = fig.axes[ax_index]

        ## Get bottom surf figure first
        if np.all(overlay!=None) and isinstance(overlay, ScalarField):
            kwargs['a_range'] = [0., 1.]; kwargs['b_range'] = [0., 1.] # no periodicity
            kwargs['edgeplot'] = False; kwargs['figsize'] = figsize
            kwargs['fig'] = fig; kwargs['ax_index'] = ax_index;
            kwargs['x_ticks'] = x_ticks; kwargs['y_ticks'] = y_ticks
            kwargs['unit'] = overlay.unit # convert unit in wrappers!!!

            fig = overlay.plot_2D(**kwargs)
            ax = fig.axes[ax_index]

        ## rotate the trajectory to plotting plane, base defined in OAB
        rot, disp = GridRotation2D(np.vstack([self.base[1], self.base[2], self.base[0]]))
        ## plot TRAJ
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
                    ax.scatter(traj[0,0], traj[0,1], marker=cpt_marker, c=cpt_color, s=cpt_size)
            # plot TRAJ
            else:
                plttraj = [] # traj in plot plane
                for v in traj:
                    v = v - baserot[1]
                    if v[0]>=0 and v[0]<xmx and v[1]>=0 and v[1]<ymx and np.abs(v[2])<=1e-3:
                        plttraj.append(v+baserot[1])

                if len(plttraj) == 0:
                    continue
                plttraj = np.array(plttraj)
                if wt != 0: # Molegraph path
                    ax.plot(plttraj[:, 0], plttraj[:, 1], color=cpt_color,
                            linestyle='-', linewidth=traj_linewidth+extra_width[wt])
                else: # other pathes
                    ax.plot(plttraj[:, 0], plttraj[:, 1], color=traj_color,
                            linestyle=traj_linestyle, linewidth=traj_linewidth)

        ax.set_aspect(1.0)
        ax.set_xlim(baserot[1, 0], baserot[1, 0]+xmx)
        ax.set_ylim(baserot[1, 1], baserot[1, 1]+ymx)
        return fig


class ChargeDensity(ScalarField):
    """
    The charge density object from TOPOND. Unit: 'Angstrom' for
    :math:`\\AA^{-3}` and 'a.u.' for Bohr :math:`^{-3}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFRHOO'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``ChargeDensity`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFRHOO.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFRHOO')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.2f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.array([0.02, 0.04, 0.08, 0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80, 200],
                                  dtype=float)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$|e|/\AA^{-3}$'
            else: ustr = r'$|e|/Bohr^{-3}$'
            if self.subtracted == False: cbar_label=r'$\rho$ ({})'.format(ustr)
            else: cbar_label=r'$\Delta\rho$ ({})'.format(ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Charge Density'
                else: title = 'Charge Density + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``ChargeDensity`` object.

        Args:
            unit (str): ''Angstrom', :math:`\\AA^{-3}`; 'a.u.', Bohr :math:`^{-3}`.
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

        lprops = ['base'] # length units
        dprops = ['data'] # density units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for d in dprops:
            newattr = getattr(self, d) / cst**3
            setattr(self, d, newattr)
        return self


class SpinDensity(ScalarField):
    """
    The spin density object from TOPOND. Unit: 'Angstrom' for :math:`\\AA^{-3}`
    and 'a.u.' for Bohr :math:`^{-3}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFSPDE'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``SpinDensity`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFSPDE.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (SpinDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFSPDE')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.4f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$|e|/\AA^{-3}$'
            else: ustr = r'$|e|/Bohr^{-3}$'
            if self.subtracted == False: cbar_label=r'$\rho$ ({})'.format(ustr)
            else: cbar_label=r'$\Delta\rho$ ({})'.format(ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Spin Density'
                else: title = 'Spin Density + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``SpinDensity`` object.

        Args:
            unit (str): ''Angstrom', :math:`\\AA^{-3}`; 'a.u.', Bohr :math:`^{-3}`.
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

        lprops = ['base'] # length units
        dprops = ['data'] # density units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for d in dprops:
            newattr = getattr(self, d) / cst**3
            setattr(self, d, newattr)
        return self


class Gradient(ScalarField):
    """
    The charge density gradient object from TOPOND. Unit: 'Angstrom' for
    :math:`\\AA^{-4}` and 'a.u.' for Bohr :math:`^{-4}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFGRHO'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``Gradient`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFGRHO.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (Gradient)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFGRHO')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.2f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.array([0.02, 0.04, 0.08, 0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80, 200],
                                  dtype=float)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$|e|/\AA^{-4}$'
            else: ustr = r'$|e|/Bohr^{-4}$'
            if self.subtracted == False: cbar_label=r'$\nabla\rho$ ({})'.format(ustr)
            else: cbar_label=r'$\Delta(\nabla\rho)$ ({})'.format(ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Density Gradient'
                else: title = 'Density Gradient + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``Gradient`` object.

        Args:
            unit (str): ''Angstrom', :math:`\\AA^{-4}`; 'a.u.', Bohr :math:`^{-4}`.
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

        lprops = ['base'] # length units
        gprops = ['data'] # density gradient units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for g in gprops:
            newattr = getattr(self, g) / cst**4
            setattr(self, g, newattr)
        return self


class Laplacian(ScalarField):
    """
    The Laplacian object from TOPOND. Unit: 'Angstrom' for :math:`\\AA^{-5}`
    and 'a.u.' for Bohr :math:`^{-5}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFLAPP'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``Laplacian`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

            It is suggested to name the file with 'LAPP' or 'LAPM'. Otherwise
            it will be read as 'LAPP'. Data is always saved as :math:`\\nabla^{2}\\rho`.

        Args:
            file (str): File name of the 'SURFLAPP.DAT' or 'SURFLAPM.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (Laplacian)
        """
        from CRYSTALpytools.crystal_io import Properties_output
        import warnings

        if 'LAPM' in file.upper(): type = 'SURFLAPM'
        elif 'LAPP' in file.upper(): type = 'SURFLAPP'
        else:
            warnings.warn("Type of data not available from the file name. Using 'SURFLAPP'.",
                          stacklevel=2)
            type = 'SURFLAPP'
        return Properties_output(output).read_topond(file, type)

    def plot_2D(
        self, unit='Angstrom', plot_lapm=False, levels='default', lineplot=True,
        linewidth=1.0, isovalues='%.2f', colorplot=False, colormap='jet',
        cbar_label='default', a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            plot_lapm (bool): Whether to plot :math:`-\\nabla^{2}\\rho`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # lapm
        if plot_lapm == True:
            self.data = -self.data

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, 0, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if plot_lapm == True: pm = '-'
        else: pm= '+'
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$|e|/\AA^{-5}$'
            else: ustr = r'$|e|/Bohr^{-5}$'
            if self.subtracted == False: cbar_label=r'{}$\nabla^2\rho$ ({})'.format(pm, ustr)
            else: cbar_label=r'$\Delta({}\nabla^2\rho)$ ({})'.format(pm, ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Density Laplacian ({})'.format(pm)
                else: title = 'Density Laplacian ({}) + {}'.format(pm, titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        # restore old unit and laplacian
        self._set_unit(uold)
        if plot_lapm == True:
            self.data = -self.data
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``Laplacian`` object.

        Args:
            unit (str): ''Angstrom', :math:`\\AA^{-5}`; 'a.u.', Bohr :math:`^{-5}`.
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

        lprops = ['base'] # length units
        laprops = ['data'] # laplacian units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for la in laprops:
            newattr = getattr(self, la) / cst**5
            setattr(self, la, newattr)
        return self


class HamiltonianKE(ScalarField):
    """
    The Hamiltonian kinetic energy density object from TOPOND. Unit: 'Angstrom'
    for eV.:math:`\\AA^{-3}` and 'a.u.' for Hartree.Bohr :math:`^{-3}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFKKIN'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``HamiltonianKE`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFKKIN.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (HamiltonianKE)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFKKIN')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.2f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, 0, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$eV/\AA^{-3}$'
            else: ustr = r'$Hartree/Bohr^{-3}$'
            if self.subtracted == False: cbar_label=r'$E_k$ ({})'.format(ustr)
            else: cbar_label=r'$\Delta E_k$ ({})'.format(ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Hamiltonian Kinetic Energy'
                else: title = 'Hamiltonian Kinetic Energy + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``HamiltonianKE`` object.

        Args:
            unit (str): ''Angstrom', eV.:math:`\\AA^{-3}`; 'a.u.',
                Hartree.Bohr :math:`^{-3}`.
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, H_to_eV, eV_to_H

        if unit.lower() == self.unit.lower():
            return self

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

        lprops = ['base'] # length units
        eprops = ['data'] # energy density units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for g in gprops:
            newattr = getattr(self, g) / cst**3 * ecst
            setattr(self, g, newattr)
        return self


class LagrangianKE(ScalarField):
    """
    The Lagrangian kinetic energy density object from TOPOND. Unit: 'Angstrom'
    for eV.:math:`\\AA^{-3}` and 'a.u.' for Hartree.Bohr :math:`^{-3}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFGKIN'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``LagrangianKE`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFGKIN.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (LagrangianKE)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFGKIN')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.2f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.array([0.02, 0.04, 0.08, 0.2, 0.4, 0.8,
                                   2, 4, 8, 20, 40, 80, 200], dtype=float)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$eV/\AA^{-3}$'
            else: ustr = r'$Hartree/Bohr^{-3}$'
            if self.subtracted == False: cbar_label=r'$E_k$ ({})'.format(ustr)
            else: cbar_label=r'$\Delta E_k$ ({})'.format(ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Lagrangian Kinetic Energy'
                else: title = 'Lagrangian Kinetic Energy + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``LagrangianKE`` object.

        Args:
            unit (str): ''Angstrom', eV.:math:`\\AA^{-3}`; 'a.u.',
                Hartree.Bohr :math:`^{-3}`.
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, H_to_eV, eV_to_H

        if unit.lower() == self.unit.lower():
            return self

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

        lprops = ['base'] # length units
        eprops = ['data'] # energy density units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for g in gprops:
            newattr = getattr(self, g) / cst**3 * ecst
            setattr(self, g, newattr)
        return self


class VirialField(ScalarField):
    """
    The Virial field density object from TOPOND. Unit: 'Angstrom'
    for eV.:math:`\\AA^{-3}` and 'a.u.' for Hartree.Bohr :math:`^{-3}`.

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFVIRI'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``VirialField`` object from a single file. Can be used for
        multiple dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFVIRI.DAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (VirialField)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFVIRI')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.2f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if np.all(levels=='default') and unit.lower() != 'angstrom':
            warnings.warn("Unit must be 'Angstrom' when the default levels are set. Using Angstrom-eV units.",
                          stacklevel=2)
            unit = 'Angstrom'
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        # overlay's unit
        if np.all(overlay!=None):
            if not isinstance(overlay, Trajectory):
                raise ValueError("The overlaied layer must be a topond.Trajectory class or its child classes.")
            overlay._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.array([0.02, 0.04, 0.08, 0.2, 0.4, 0.8,
                                   2, 4, 8, 20, 40, 80, 200], dtype=float)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = -1e-6; rlimit = 1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if unit.lower() == 'angstrom': ustr = r'$eV/\AA^{-3}$'
            else: ustr = r'$Hartree/Bohr^{-3}$'
            if self.subtracted == False: cbar_label=r'$VF$ ({})'.format(ustr)
            else: cbar_label=r'$\Delta VF$ ({})'.format(ustr)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Virial Field'
                else: title = 'Virial Field + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``VirialField`` object.

        Args:
            unit (str): ''Angstrom', eV.:math:`\\AA^{-3}`; 'a.u.',
                Hartree.Bohr :math:`^{-3}`.
        """
        from CRYSTALpytools.units import angstrom_to_au, au_to_angstrom, H_to_eV, eV_to_H

        if unit.lower() == self.unit.lower():
            return self

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

        lprops = ['base'] # length units
        eprops = ['data'] # energy density units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        for g in gprops:
            newattr = getattr(self, g) / cst**3 * ecst
            setattr(self, g, newattr)
        return self


class ELF(ScalarField):
    """
    The electron localization object from TOPOND. Dimensionless. Unit: 'Angstrom'
    for :math:`\\AA` and 'a.u.' for Bohr (only for plot base vectors).

    Args:
        data (array): 2D (3D) Plot data. (nZ\*)nY\*nX. **3D methods under developing**.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
    """
    def __init__(self, data, base, dimen, struc=None, unit='Angstrom'):
        import numpy as np

        self.data = np.array(data, dtype=float)
        self.base = np.array(base, dtype=float)
        self.dimension = int(dimen)
        self.structure = struc
        self.unit = unit
        self.type = 'SURFELFB'
        self.subtracted = False # Hidden. For plotting.

    @classmethod
    def from_file(cls, file, output=None):
        """
        Generate a ``ELF`` object from a single file. Can be used for multiple
        dimensions (2D only now. 3D under development).

        .. note::

            Though output is not mandatory for plotting proposes, it is highly
            recommended to be added for geometry information and other methods.

        Args:
            file (str): File name of the 'SURFELFBDAT' file
            output (str): Screen output of 'properties' calculation.
        Returns:
            cls (ELF)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'SURFELFB')

    def plot_2D(
        self, unit='Angstrom', levels='default', lineplot=True, linewidth=1.0,
        isovalues='%.2f', colorplot=False, colormap='jet', cbar_label='default',
        a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
        x_ticks=5, y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None,
        fig=None, ax_index=None, **kwargs
    ):
        """
        Plot 2D contour lines, color maps or both for the 2D data set. The user
        can also get the overlapped plot of ``ScalarField`` and ``Trajectory``
        classes.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.ScalarField``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            levels (array|int): Set levels of colored / line contour plot. A
                number for linear scaled plot colors or an array for
                user-defined levels. 1D. 'default' for default levels.
            lineplot (bool): Plot contour lines.
            contourline (list): Width of contour lines. Other styles are pre-
                set and cannot be changed.
            isovalues (str|None): Add isovalues to contour lines and set their
                formats. Useful only if ``lineplot=True``. None for not adding
                isovalues.
            colorplot (bool): Plot color-filled contour plots.
            colormap (str): Matplotlib colormap option. Useful only if
                ``colorplot=True``.
            cbar_label (str): Label of colorbar. Useful only if
                ``colorplot=True``. 'None' for no labels.
            a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in
                fractional coordinate.
            b_range (list): 1\*2 range of :math:`b` axis (y, or BA) in
                fractional coordinate.
            edgeplot (bool): Whether to add cell boundaries represented by the
                original base vectors (not inflenced by a/b range).
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for default values and
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|Trajectory): Overlapping a 2D plot from the
                ``topond.Trajectory`` object if not None.
            fig (Figure): *Developer Only*, matplotlib Figure class.
            ax_index (int): *Developer Only*, index of the axis in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``topond.Trajectory``.

        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)

        # default levels
        if np.all(levels=='default'):
            if self.subtracted == False:
                levels = np.linspace(0, 1, 21)
            else:
                levels = np.array([-80, -40, -20, -8, -4, -2, -0.8, -0.4, -0.2,
                                   -0.08, -0.04, -0.02, -0.008, -0.004, -0.002,
                                   0, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08,
                                   0.2, 0.4, 0.8, 2, 4, 8, 20, 40, 80], dtype=float)
        # contour line styles
        blimit = 0.5-1e-6; rlimit = 0.5+1e-6
        if colorplot == False:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['b', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['r', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        else:
            contourline = []
            for i in levels:
                if i < blimit: contourline.append(['k', 'dotted', linewidth])
                elif i > rlimit: contourline.append(['k', '-', linewidth])
                else: contourline.append(['k', '-', linewidth*2])
        # cbar label
        if cbar_label=='default':
            if self.subtracted == False: cbar_label=r'ELF'
            else: cbar_label=r'$\Delta$ ELF'

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0
        # plot
        fig = super().plot_2D(levels, lineplot, contourline, isovalues,
                              colorplot, colormap, cbar_label, a_range,
                              b_range, edgeplot, x_ticks, y_ticks, figsize,
                              overlay, fig, ax_index, **kwargs)

        # label and title
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'TRAJGRAD' : 'Gradient Trajectory',
                      'TRAJMOLG' : 'Chemical Graph'}
            if title == 'default':
                if np.all(overlay==None): title = 'Electron Localization'
                else: title = 'Electron Localization + {}'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``ELF`` object.

        Args:
            unit (str): ''Angstrom', :math:`\\AA`; 'a.u.', Bohr. Only for plot
                base vectors
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

        lprops = ['base'] # length units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        return self


class GradientTraj(Trajectory):
    """
    The charge density gradient trajectory object from TOPOND. Unit: 'Angstrom'
    for :math:`\\AA` and 'a.u.' for Bohr.

    Args:
        wtraj (list[int]): 1\*nPath, weight of the path, int 0 to 3.
        traj (list[array]): 1\*nPath, list of paths. Every array is a nPoint\*3
            3D ref framework coordinates of points on the path.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).

    Returns:
        self (GradientTraj): Import attributes: ``trajectory``, nPath\*2 list of
            trajectory and its weight; ``cpoint``, nCpoint\*3 array of critical
            point coordinates in 3D ref framework.
    """
    def __init__(self, wtraj, traj, base, struc, unit='Angstrom'):
        if len(wtraj) != len(traj):
            raise ValueError('Inconsistent lengths of input trajectory and its weight.')

        self.trajectory = [[int(wtraj[i]), np.array(traj[i], dtype=float)]
                           for i in range(len(wtraj))]
        self.base = np.array(base, dtype=float)
        self.structure = struc
        self.type = 'TRAJGRAD'
        self.unit = unit
        cpt = []
        for i in self.trajectory:
            if len(i[1]) == 1:
                cpt.append(i[1][0])
        self.cpoint = np.array(cpt)

    @classmethod
    def from_file(cls, file, output):
        """
        Generate a ``GradientTraj`` object from 'TRAJGRAD.DAT' file and the
        standard screen output of CRYSTAL 'properties' calculation.

        .. note::

            Output file is mandatory.

        Args:
            file (str): 'TRAJGRAD.DAT' file by TOPOND.
            output (str): Standard output of Properties calculation, used to
                get geometry and orientation of the plane.
        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'TRAJGRAD')

    def plot_2D(
        self, unit='Angstrom', cpt_marker='o', cpt_color='k', cpt_size=10,
        traj_color='r', traj_linestyle=':', traj_linewidth=0.5, x_ticks=5,
        y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None, fig=None,
        ax_index=None, **kwargs):
        """
        Plot charge density gradient trajectory in a 2D plot.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.Trajectory``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            cpt_marker (str): Marker of critical point scatter.
            cpt_color (str): Marker color of critical point scatter.
            cpt_size (float|int): Marker size of critical point scatter.
            traj_color (str): Line color of 2D trajectory plot.
            traj_linestyl (str): Line style of 2D trajectory plot.
            traj_linewidth (str): Line width of 2D trajectory plot.
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for proeprty plotted.
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|ScalarField): Overlapping a 2D plot from the
                ``topond.ScalarField`` object if not None.
            fig (Figure|None): *Developers only*, matplotlib Figure class..
            ax_index (list[int]): *Developer Only*, indices of axes in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``ScalarField.plot_2D()``.
        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        if np.all(overlay!=None):
            overlay._set_unit(unit)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0

        # plot
        fig = super().plot_2D(cpt_marker, cpt_color, cpt_size, traj_color,
                              traj_linestyle, traj_linewidth, x_ticks, y_ticks,
                              figsize, overlay, fig, ax_index, **kwargs)
        # labels and titles
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'SURFRHOO' : 'Charge Density',
                      'SURFSPDE' : 'Spin Density',
                      'SURFGRHO' : 'Density Gradient',
                      'SURFLAPP' : 'Density Laplacian',
                      'SURFKKIN' : 'Hamiltonian Kinetic Eenergy',
                      'SURFGKIN' : 'Lagrangian Kinetic Eenergy',
                      'SURFVIRI' : 'Virial Field',
                      'SURFELFB' : 'Electron Localization'}
            if title == 'default':
                if np.all(overlay==None): title = 'Gradient Trajectory'
                else: title = '{} + Gradient Trajectory'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``GradientTraj`` object.

        Args:
            unit (str): 'Angstrom' or 'a.u.'
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

        lprops = ['base', 'cpoint'] # length units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        # special: trajectory
        self.trajectory = [[i[0], i[1]*cst] for i in self.trajectory]
        return self


class ChemicalGraph(Trajectory):
    """
    The molecular / crystal graph object from TOPOND. Unit: 'Angstrom' for
    :math:`\\AA` and 'a.u.' for Bohr.

    Args:
        wtraj (list[int]): 1\*nPath, weight of the path, int 0 to 3.
        traj (list[array]): 1\*nPath, list of paths. Every array is a nPoint\*3
            3D ref framework coordinates of points on the path.
        base (array): 3\*3 Cartesian coordinates of the 3 points defining
            vectors BA and BC.
        struc (CStructure): Extended Pymatgen Structure object.
        unit (str): In principle, should always be 'Angstrom' (case insensitive).

    Returns:
        self (ChemicalGraph): Import attributes: ``trajectory``, nPath\*2 list
            of trajectory and its weight; ``cpoint``, nCpoint\*3 array of
            critical point coordinates in 3D ref framework.
    """
    def __init__(self, wtraj, traj, base, struc, unit='Angstrom'):
        if len(wtraj) != len(traj):
            raise ValueError('Inconsistent lengths of input trajectory and its weight.')

        self.trajectory = [[int(wtraj[i]), np.array(traj[i], dtype=float)]
                           for i in range(len(wtraj))]
        self.base = np.array(base, dtype=float)
        self.structure = struc
        self.type = 'TRAJMOLG'
        self.unit = unit
        cpt = []
        for i in self.trajectory:
            if len(i[1]) == 1:
                cpt.append(i[1][0])
        self.cpoint = np.array(cpt)

    @classmethod
    def from_file(cls, file, output):
        """
        Generate a ``ChemicalGraph`` object from 'TRAJMOLG.DAT' file and the
        standard screen output of CRYSTAL 'properties' calculation.

        .. note::

            Output file is mandatory.

        Args:
            file (str): 'TRAJMOLG.DAT' file by TOPOND.
            output (str): Standard output of Properties calculation, used to
                get geometry and orientation of the plane.
        Returns:
            cls (ChargeDensity)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_topond(file, 'TRAJMOLG')

    def plot_2D(
        self, unit='Angstrom', cpt_marker='o', cpt_color='k', cpt_size=10,
        traj_color='r', traj_linestyle=':', traj_linewidth=0.5, x_ticks=5,
        y_ticks=5, title='default', figsize=[6.4, 4.8], overlay=None, fig=None,
        ax_index=None, **kwargs):
        """
        Plot crystal / molecular graph in a 2D plot.

        .. note::

            For more information of plot setups, please refer its parent class,
            ``topond.Trajectory``.

        Args:
            unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
                Bohr :math:`^{-3}`.
            cpt_marker (str): Marker of critical point scatter.
            cpt_color (str): Marker color of critical point scatter.
            cpt_size (float|int): Marker size of critical point scatter.
            traj_color (str): Line color of 2D trajectory plot.
            traj_linestyl (str): Line style of 2D trajectory plot.
            traj_linewidth (str): Line width of 2D trajectory plot.
            x_ticks (int): Number of ticks on x axis.
            y_ticks (int): Number of ticks on y axis.
            title (str|None): Plot title. 'default' for proeprty plotted.
                'None' for no title.
            figsize (list): Matplotlib figure size. Note that axes aspects are
                fixed to be equal.
            overlay (None|ScalarField): Overlapping a 2D plot from the
                ``topond.ScalarField`` object if not None.
            fig (Figure|None): *Developers only*, matplotlib Figure class..
            ax_index (list[int]): *Developer Only*, indices of axes in ``fig.axes``.
            \*\*kwargs : Other arguments passed to ``ScalarField.plot_2D()``.
        Returns:
            fig (Figure): Matplotlib Figure object
        """
        import numpy as np
        import warnings

        # unit
        uold = self.unit
        if self.unit.lower() != unit.lower():
            self._set_unit(unit)
        if np.all(overlay!=None):
            overlay._set_unit(unit)

        # axis index
        if np.all(fig!=None) and np.all(ax_index!=None):
            iax = int(ax_index)
        else:
            iax = 0

        # plot
        fig = super().plot_2D(cpt_marker, cpt_color, cpt_size, traj_color,
                              traj_linestyle, traj_linewidth, x_ticks, y_ticks,
                              figsize, overlay, fig, ax_index, **kwargs)
        # labels and titles
        if unit.lower() == 'angstrom':
            fig.axes[iax].set_xlabel(r'$\AA$')
            fig.axes[iax].set_ylabel(r'$\AA$')
        else:
            fig.axes[iax].set_xlabel(r'$Bohr$')
            fig.axes[iax].set_ylabel(r'$Bohr$')

        if np.all(title!=None):
            titles = {'SURFRHOO' : 'Charge Density',
                      'SURFSPDE' : 'Spin Density',
                      'SURFGRHO' : 'Density Gradient',
                      'SURFLAPP' : 'Density Laplacian',
                      'SURFKKIN' : 'Hamiltonian Kinetic Eenergy',
                      'SURFGKIN' : 'Lagrangian Kinetic Eenergy',
                      'SURFVIRI' : 'Virial Field',
                      'SURFELFB' : 'Electron Localization'}
            if title == 'default':
                if np.all(overlay==None): title = 'Chemical Graph'
                else: title = '{} + Chemical Graph'.format(titles[overlay.type.upper()])
            fig.axes[iax].set_title(title)

        self._set_unit(uold)
        if np.all(overlay!=None): overlay._set_unit(uold)
        return fig

    def _set_unit(self, unit):
        """
        Set units of data of ``ChemicalGraph`` object.

        Args:
            unit (str): 'Angstrom' or 'a.u.'
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

        lprops = ['base', 'cpoint'] # length units
        for l in lprops:
            newattr = getattr(self, l) * cst
            setattr(self, l, newattr)
        # special: trajectory
        self.trajectory = [[i[0], i[1]*cst] for i in self.trajectory]
        return self


