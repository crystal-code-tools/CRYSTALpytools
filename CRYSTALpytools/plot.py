#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to visualize CRYSTAL outputs.
"""
import numpy as np

##############################################################################
#                                                                            #
#                       ELECTRONIC STRUCTURE                                 #
#                                                                            #
##############################################################################

#--------------------------ECHG charge and spin density----------------------#

def plot_ECHG(
    *echg, unit='Angstrom', output=[], option='both', levels=150,
    lineplot=False, linewidth=1.0, isovalues=None, colorplot=True,
    colormap='jet', cbar_label='default', a_range=[0.,1.], b_range=[0.,1.],
    rectangle=False, edgeplot=False, x_ticks=5, y_ticks=5, layout=None,
    title=None, figsize=[6.4, 4.8], sharex=True, sharey=True, fontsize=14,
    **kwargs):
    """
    Read and plot multiple 2D charge density files / objects. The uniform plot
    set-ups are used for comparison.

    3 styles are available:

    1. ``lineplot=True`` and ``colorplot=True``: The color-filled contour map
        with black contour lines. Dotted lines for negative values and solid
        lines for positive values. The solid line twice in width for 0.  
    2. ``lineplot=False`` and ``colorplot=True``: The color-filled contour map.  
    3. ``lineplot=True`` and ``colorplot=False``: The color coded contour line
        map. Blue dotted line for negative values and red solid lines for
        positive values. The balck solid line twice in width for 0.

    Available options:

    * 'both' : If spin polarized, plot both charge and spin densities.
        Otherwise plot charge densities.  
    * 'charge': Plot charge density.  
    * 'spin': Plot spin density.  
    * 'diff': Subtracting charge data from the first entry with the following
        entries. Return to a non spin-polarized object.  

    Args:
        \*echg (electronics.ChargeDensity|str): Extendable. File names or
            ``electronics.ChargeDensity`` objects.
        unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
            Bohr :math:`^{-3}`.
        output (str|list[str]): Output files corresponding to the input data
            file. String for the same output of all the files and list for
            every input file. If not given, read the first 1(2) MAPNET data of
            systems without (with) spin.
        option (str): Available options see above.
        levels (int|array): Set levels of contour plot. A number for linear
            scaled plot colors or an array for user-defined levels, **must be
            consistent with ``unit``**. 2\*nLevel can be defined when
            ``option='both'``.
        lineplot (bool): Plot contour lines.
        linewidth (float): Contour linewidth. Useful only if ``lineplot=True``.
            Other properties are not editable. Solid black lines for positive
            values and 0, dotted for negative.
        isovalues (str|None): Add isovalues to contour lines and set their
            formats. Useful only if ``lineplot=True``. None for not adding isovalues.
        colorplot (bool): Plot color-filled contour plots.
        colormap (str): Matplotlib colormap option. Useful only if ``colorplot=True``.
        cbar_label (str): Label of colorbar. Useful only if ``colorplot=True``.
            1\*2 list of colorbar titles can be set for spin-polarized systems.
            'None' for no labels and 'default' for default label.
        a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in fractional coordinate.
        rectangle (bool): If :math:`a, b` are non-orthogonal, plot a rectangle
            region and reset :math:`b`. If used together with ``b_range``, that
            refers to the old :math:`b`.
        edgeplot (bool): Whether to add cell edges represented by the original
            base vectors (not inflenced by a/b range or rectangle options).
        x_ticks (int): Number of ticks on x axis.
        y_ticks (int): Number of ticks on y axis.
        layout (list|tuple): The layout of subplots, \[nrow, ncol\]. Default is
            2 cols per row.
        title (str|None): The title of the plot. 'None' for no title. The
            default subplot titles are used either way.
        figsize (list): Matplotlib figure size. Note that axes aspects are
            fixed to be equal.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the heightest level of title.
        \*\*kwargs : Other arguments passed to ``axes.contour()`` function to
            set contour lines. Applied to all the subplots.

    Returns:
        figs (list|Figure): Matplotlib Figure object or a list of them.
    """
    from CRYSTALpytools.electronics import ChargeDensity
    import numpy as np
    import matplotlib.pyplot as plt

    # output
    output = np.array(output, ndmin=1)
    if len(output) == 1: # same output for all the inputs
        output = np.repeat(output, len(echg))
    elif len(output) == 0: # no output
        output = [None for i in range(len(echg))]
    # get objects
    obj = []; countstr = 0
    for i in echg:
        if isinstance(i, str):
            if countstr >= len(output):
                raise ValueError("Inconsistent length of input file name and output.")
            obj.append(ChargeDensity.from_file(i, output=output[countstr]))
        elif isinstance(i, ChargeDensity):
            obj.append(i)
        else:
            raise TypeError("Inputs must be either string or electronics.ChargeDensity objects.")

    # subtraction
    if 'diff' in option.lower():
        obj[0].subtract(*[i for i in obj[1:]])
        option = 'charge'
        obj = [obj[0]]
    # set uniform levels
    if isinstance(levels, float) or isinstance(levels, int):
        spin_range = []
        chg_range = []
        for i in obj:
            if i.spin == 1:
                chg_range.append([np.min(i.data), np.max(i.data)])
            else:
                chg_range.append([np.min(i.data[:, :, 0]), np.max(i.data[:, :, 0])])
                spin_range.append([np.min(i.data[:, :, 1]), np.max(i.data[:, :, 1])])
        if spin_range == []:
            levels1 = np.linspace(np.min(chg_range), np.max(chg_range), levels)
            levels2 = levels1
        else:
            levels1 = np.linspace(np.min(chg_range), np.max(chg_range), levels)
            levels2 = np.linspace(np.min(spin_range), np.max(spin_range), levels)
    else:
        levels = np.array(levels, dtype=float, ndmin=2)
        if levels.shape[0] == 1:
            levels1 = levels
            levels2 = levels
        else:
            levels1 = levels[0]
            levels2 = levels[1]
    levels = np.vstack([levels1, levels2])

    # layout of plots
    if option.lower() == 'both':
        nplt = int(len(obj) * 2)
    else:
        nplt = len(obj)

    if np.all(layout!=None):
        if layout[0]*layout[1] < nplt:
            warnings.warn('Layout size smaller than the number of plots. Using default layout.',
                          stacklevel=2)
            layout = None
    if np.all(layout==None):
        if nplt == 1:
            layout = [1, 1]
        else:
            layout = [int(np.ceil(nplt/2)), 2]

    fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                           sharex=sharex, sharey=sharey, layout='tight')

    # plot
    ax_index = 0
    for o in obj:
        if option.lower() == 'both':
            iax = [ax_index, ax_index+1]
            fig = o.plot_2D(unit, 'both', levels, lineplot, linewidth, isovalues,
                            colorplot, colormap, cbar_label, a_range, b_range,
                            rectangle, edgeplot, x_ticks, y_ticks, 'default',
                            figsize, fig, iax, **kwargs)
            ax_index += 2
        else:
            fig = o.plot_2D(unit, option.lower(), levels, lineplot, linewidth,
                            isovalues, colorplot, colormap, cbar_label, a_range,
                            b_range, rectangle, edgeplot, x_ticks, y_ticks,
                            'default', figsize, fig, [ax_index], **kwargs)
            ax_index += 1
    # title
    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)
    return fig


#------------------------------2c-SCF vector field----------------------------#

def plot_relativistics2D(
    *relat, unit='SI', type=[], output=[], direction=['x','y','z'], levels=100,
    quiverplot=True, quiverscale=1.0, colorplot=True, colormap='jet',
    cbar_label='default', a_range=[0., 1.], b_range=[0., 1.], rectangle=False,
    edgeplot=False, x_ticks=5, y_ticks=5, layout=None, title=None, figsize=[6.4, 4.8],
    sharex=True, sharey=True, fontsize=14, **kwargs):
    """
    Plot 2D vector field properties from relativistics (2c-SCF) calculations.

    3 styles are available:

    1. ``quiverplot=True`` and ``colorplot=True``: The color-filled contour
        illustrates the norm of vectors. The black arrows indicates both the
        directions and norms of in-plane prjections.  
    2. ``quiverplot=True`` and ``colorplot=False``: The arrows are colored to
        indicate the directions and norms of in-plane prjections.  
    3. ``quiverplot=False`` and ``colorplot=True``: The color-filled contour
        illustrates the norm of vectors, similar to the 2D scalar map.

    .. note::

        Not for charge density (``relativistics.ChargeDensity``). To visualize
        it, use ``plot_ECHG``.

    Args:
        \*relat (str|Magnetization|OrbitalCurrentDensity|SpinCurrentDensity|):
            extendable input of vector field classes, or input files.
        unit (str): Plot unit. 'SI' for :math:`\\AA` and A/m (A/m :math:`^{2}`).
            'a.u.' for Bohr and a.u. magnetization / current density.
        type (str|list[str]): Properties to plot. Either as a string or a list
            of strings consistent with input filenames. If a list of types and
            a single file is given, read all the types from the input file. In
            other cases error would be given.
        output (str|list[str]): Output files corresponding to the input data
            file. String for the same output of all the files and list for
            every input file.
        direction (list[str]|str): Only for ``SpinCurrentDensity`` classes.
            Direction of spin-current to plot, in 'x', 'y' or 'z'. Applied to
            all the ``SpinCurrentDensity`` objects.
        levels (int|array): Set levels of colored contour/quiver plot. A number
            for linear scaled plot colors or an array for user-defined levels. 1D.
        quiverplot (bool): Plot 2D field of arrows.
        quiverscale (float): Tune the length of arrows. Useful only if ``quiverplot=True``.
        colorplot (bool): Plot color-filled contour plots.
        colormap (str): Matplotlib colormap option. Useful only if ``colorplot=True``.
        cbar_label (str|None): Label of colorbar. 'default' for unit. 'None'
            for no label. Useful only if ``colorplot=True``.
        a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in fractional coordinate.
        rectangle (bool): If :math:`a, b` are non-orthogonal, plot a rectangle
            region and reset :math:`b`. If used together with ``b_range``, that
            refers to the old :math:`b` (i.e., expansion first).
        edgeplot (bool): Whether to add cell edges represented by the original
            base vectors (not inflenced by a/b range or rectangle options).
        x_ticks (int): Number of ticks on x axis.
        y_ticks (int): Number of ticks on y axis.
        layout (list|tuple): The layout of subplots, \[nrow, ncol\]. Default is
            2 cols per row.
        title (str|None): The title of the plot. 'None' for no title. The
            default subplot titles are used either way.
        figsize (list): Matplotlib figure size. Note that axes aspects are
            fixed to be equal.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the heightest level of title.
        \*\*kwargs : Other arguments passed to ``axes.quiver()`` function to
            set arrow styles. Applied to all the subplots.

        Returns:
            fig (Figure): Matplotlib Figure class.
    """
    from CRYSTALpytools.relativistics import ChargeDensity, Magnetization, \
                                             OrbitalCurrentDensity,SpinCurrentDensity
    import numpy as np
    from CRYSTALpytools.crystal_io import Properties_output
    import matplotlib.pyplot as plt
    import warnings

    type = np.array(type, ndmin=1)
    if len(type) == 1: # same type for all the inputs
        type = np.repeat(type, len(relat))
    output = np.array(output, ndmin=1)
    if len(output) == 1: # same output for all the inputs
        output = np.repeat(output, len(relat))

    objs = []
    if len(relat) == 1: # single input, allow for multiple types.
        if isinstance(relat[0], str):
            if len(output) == 0:
                raise ValueError("Outputs must be set for input files. Otherwise use input objects.")
            for t in type:
                objs.append(Properties_output(output[0]).read_relativistics(relat[0], t))
        elif isinstance(relat[0], Magnetization) \
        or isinstance(relat[0], OrbitalCurrentDensity) \
        or isinstance(relat[0], SpinCurrentDensity):
            objs = [relat[0]]
        elif isinstance(relat[0], ChargeDensity):
            raise TypeError("Use 'plot_ECHG' for charge densities from 2c-SCF calulations.")
        else:
            raise TypeError('Unknown input type.')
    else: # multiple inputs.
        countstr = 0
        for r in relat:
            if isinstance(r, str):
                if len(output) == 0:
                    raise ValueError("Outputs must be set for input files. Otherwise use input objects.")
                if countstr >= len(type) or countstr >= len(output):
                    raise ValueError("Inconsistent length of input file name and output / type.")
                objs.append(
                    Properties_output(output[countstr]).read_relativistics(r, type=type[countstr])
                )
                countstr += 1
            elif isinstance(r, Magnetization) \
            or isinstance(r, OrbitalCurrentDensity) \
            or isinstance(r, SpinCurrentDensity):
                objs.append(r)
            elif isinstance(r, ChargeDensity):
                raise TypeError("Use 'plot_ECHG' for charge densities from 2c-SCF calulations.")
            else:
                raise TypeError('Unknown input type.')

    # layout of plots
    nplt = 0; direction = np.array(direction, ndmin=1)
    for o in objs:
        if isinstance(o, SpinCurrentDensity):
            nplt += len(direction)
        else:
            nplt += 1

    if np.all(layout!=None):
        if layout[0]*layout[1] < nplt:
            warnings.warn('Layout size smaller than the number of plots. Using default layout.',
                          stacklevel=2)
            layout = None
    if np.all(layout==None):
        layout = [int(np.ceil(nplt/2)), 2]

    fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                           sharex=sharex, sharey=sharey, layout='tight')
    # plot
    ax_index = 0
    for o in objs:
        if isinstance(o, SpinCurrentDensity):
            iax = [ax_index+i for i in range(len(direction))]
            fig = o.plot_2D(unit, direction, levels, quiverplot, quiverscale,
                            colorplot, colormap, cbar_label, a_range, b_range,
                            rectangle, edgeplot, x_ticks, y_ticks, 'default',
                            figsize, fig, iax, **kwargs)
            ax_index += len(direction)
        else:
            fig = o.plot_2D(unit, levels, quiverplot, quiverscale, colorplot,
                            colormap, cbar_label, a_range, b_range, rectangle,
                            edgeplot, x_ticks, y_ticks, 'default', figsize,
                            fig, [ax_index], **kwargs)
            ax_index += 1
    # title
    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)
    return fig

#--------------------------------BAND STRUCTURES------------------------------#

def plot_electron_bands(*bands, unit='eV', k_label=[], mode='single',
                        not_scaled=False, energy_range=[], k_range=[],
                        band_label=None, band_color=None, band_linestyle=None,
                        band_linewidth=None, fermi_level=0., fermi_color='tab:green',
                        fermi_linestyle='-', fermi_linewidth=1.0, layout=None,
                        title=None, figsize=[6.4, 4.8], legend='lower right',
                        sharex=True, sharey=True, fontsize=14, **kwargs):
    """
    .. _ref-plotebands:

    Plot electron band structures.

    Args:
        \*bands (ElectronBand|str): ``electronics.ElectronBand`` object or
            band structure files of the CRYSTAL properties excutable. Note that
            lattice information is not available if file names are specified.
        unit (str): The unit of energy. Can be 'eV' or 'a.u.'.
        k_label (list): nSystem\*nTick or 1\*nTick list of strings of the label
             for high symmetry points along the path. If a 1D list is given,
             the same labels are used for all the systems. `mathtext <https://matplotlib.org/stable/users/explain/text/mathtext.html>`_
             experssions can also be used  as in matplotlib.
        mode (str): The plotting mode, 'single', 'multi' and 'compare'.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
            Useful with ``mode='compare'``. The multi mode forces scaling.
        energy_range (array): A 2x1 array specifying the energy range.
        k_range (array): A 2x1 array specifying the k-range.
        band_label (str|list): Plot legend. If only one string is given, apply
            it to all plots. 1\*nSystem or nSystem\*2 (spin) plot legend
            otherwise. If spin>1 and 1\*nSystem list is used, they are marked
            with the same label.
        band_color (str|list): Color of band structure. If only one string is
            given, apply it to all plots. For 'single' and 'compare' modes,
            also 1\*2 color list for spin. For the 'multi' mode, 1\*nSystem or
            nSystem\*2 (spin) plot color. If spin>1 and spin dimension is not
            included, spin states are in the same color. 'None' for default
            values ('tab:blue' and other tab series).
        band_linestyle (str|list): Linestyle of band structure. If only one
            string is given, apply it to all plots. For 'single' and 'compare'
            also r 1\*2 linestyle list for spin. For the 'multi' mode,
            1\*nSystem or nSystem\*2 (spin) linestyle string. If spin>1 and
            spin dimension is not included, spin states are in the same style.
            'None' for default values ('-').
        band_linewidth (str|list): Linewidth of band structure. If only one
            number is given, apply it to all plots. For 'single' and 'compare'
            modes, also 1\*2 linewidth list for spin. For the 'multi' mode,
            1\*nSystem or nSystem\*2 (spin) linewidth numbers. If spin>1 and
            spin dimension is not included, spin states are in the same width.
            'None' for default values (1.0).
        fermi_level (float|list|None): Fermi energy in the same unit as input
            band energy. By default the band is aligned to 0. Can be used to
            offset the band. None for not plotting Fermi. For 'compare' mode,
            different offsets can be used.
        fermi_color (str): Color of the Fermi level.
        fermi_linestyle (str): Line style of Fermi level.
        fermi_linewidth(float): Width of the Fermi level.
        layout (list|tuple): For 'compare' mode, the layout of subplots,
            \[nrow, ncol\]. The default is 2 cols per row.
        title (str): The title of the plot.
        figsize (list): The figure size specified as \[width, height\].
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the highest level title and axis labels.
        \*\*kwargs: Other arguments passed to ``Axes.plot()`` of band plots.

    Returns:
        fig (Figure): Matplotlib figure object

    :raise ValueError: If the specified unit is unknown.
    """
    import matplotlib.pyplot as plt
    from CRYSTALpytools.electronics import ElectronBand
    from CRYSTALpytools.base.plotbase import plot_overlap_bands, plot_compare_bands
    from CRYSTALpytools.units import H_to_eV, eV_to_H
    import copy

    # unit
    if unit.lower() == 'ev':
        unit = 'eV'
        is_ev = True
    elif unit.lower() == 'a.u.':
        unit = 'a.u.'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    # instantiation and unit
    bandsplt = []
    for ib, b in enumerate(bands):
        if isinstance(b, str):
            btmp = ElectronBand.from_file(b)
        elif isinstance(b, ElectronBand):
            btmp = copy.deepcopy(b)
        else:
            raise ValueError('Unknown input type for bands.')

        if unit != btmp.unit:
            btmp._set_unit(unit)
            if np.all(fermi_level!=None):
                if unit == 'eV':
                    fermi_level = H_to_eV(fermi_level)
                else:
                    fermi_level = eV_to_H(fermi_level)
        bandsplt.append(btmp)

    # k label
    if np.all(k_label==None):
        k_label = [b.tick_label for b in bandsplt]

    # plot mode
    if mode.lower() == 'single':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = plot_compare_bands(
            [ax], [bandsplt[0].bands], [bandsplt[0].k_path], [bandsplt[0].tick_pos],
            k_label, False, energy_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, fermi_level, fermi_color,
            fermi_linestyle, fermi_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'multi':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        k_xax = [i.k_path for i in bandsplt]
        ax = plot_overlap_bands(
            ax, [b.bands for b in bandsplt], k_xax, [b.tick_pos for b in bandsplt],
            k_label, energy_range, k_range, band_label, band_color, band_linestyle,
            band_linewidth, fermi_level, fermi_color, fermi_linestyle,
            fermi_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'compare':
        if np.all(layout==None):
            layout = [int(np.ceil(len(bandsplt)/2)), 2]

        fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        k_xax = [i.k_path for i in bandsplt]
        _ = plot_compare_bands(
            ax.flat, [b.bands for b in bandsplt], k_xax, [b.tick_pos for b in bandsplt],
            k_label, not_scaled, energy_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, fermi_level, fermi_color,
            fermi_linestyle, fermi_linewidth, legend, **kwargs
        )
    else:
        raise ValueError("Unknown mode input: '{}'.".format(mode))

    # set titles and axes
    if is_ev == True:
        if np.all(fermi_level!=0.0):
            fig.supylabel('Energy (eV)', fontsize=fontsize)
        else:
            fig.supylabel('$E-E_{F}$ (eV)', fontsize=fontsize)
    else:
        if np.all(fermi_level!=0.0):
            fig.supylabel('Energy (a.u.)', fontsize=fontsize)
        else:
            fig.supylabel('$E-E_{F}$ (a.u.)', fontsize=fontsize)

    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)

    return fig


#-------------------------------DENSITY OF STATES-----------------------------#

def plot_electron_doss(*doss, unit='eV', beta='up', overlap=False, prj=[],
                       energy_range=[], dos_range=[], dos_label=None,
                       dos_color=None, dos_linestyle=None, dos_linewidth=None,
                       fermi_level=0., fermi_color='tab:green', fermi_linestyle='-',
                       fermi_linewidth=1.0, title=None, figsize=[6.4, 4.8],
                       legend='lower right', sharex=True, sharey=False,
                       fontsize=14, **kwargs):
    """
    .. _ref-plotedoss:

    Plot electron density of states.

    Args:
        \*doss (ElectronDOS|str): ``electronics.ElectronDOS`` object or DOSS
            files of the CRYSTAL properties excutable. Note that lattice
            information is not available if file names are specified.
        unit (str): 'eV' or 'a.u.'
        beta (str): Plot settings for :math:`\beta` states ('up' or 'down').
        overlap (bool): Plotting multiple projections into the same figure.
            Useful only if a single entry of ``doss`` is plotted. Otherwise
            projections from the same entry will be overlapped into the same
            subplot.
        prj (list): Index of selected projections, consistent with the first
            dimension of the ``doss``, starting from 1. Effective for all the
            subplots.
        energy_range (list): 1\*2 list of energy range
        dos_range (list): 1\*2 list of DOS range
        dos_label (str|list): Plot legend. If only one string is given, apply
            it to all plots. 1\*nPrj or nPrj\*2 (spin) plot legend otherwise.
            If spin>1 and 1\*nSystem list is used, they are marked with the
            same label. Effective for all the subplots.
        dos_color (str|list): Color of DOSS plots. If only one string is given,
            apply it to all plots. 1\*nPrj or nPrj\*2 (spin) plot color. If
            spin>1 and spin dimension is not included, spin states are in the
            same color. 'None' for default values (matplotlib Tableau Palette
            for both spin-up and spin-down states). Effective for all
            the subplots.
        dos_linestyle (str|list): Linestyle of DOSS plot. If only one string is
            given, apply it to all plots. 1\*nPrj or nPrj\*2 (spin) line
            styles. If spin>1 and spin dimension is not included, spin states
            are in the same linestyle. 'None' for default values ('-' for
            spin-ups and '--' for spin-downs). Effective for all the subplots.
        dos_linewidth (str|list): Linewidth of DOSS plot. If only one number is
            given, apply it to all plots. 1\*nPrj or nPrj\*2 (spin) line
            widthes. If spin>1 and spin dimension is not included, spin states
            are in the same linestyle. 'None' for default values (1.0).
            Effective for all the subplots.
        fermi_level (float|list|None): Fermi energy in the same unit as input
            doss energy. By default the doss is aligned to 0. Can be used to
            offset the doss. None for not plotting Fermi.
        fermi_color (str): Color of the Fermi level.
        fermi_linestyle (str): Line style of Fermi level.
        fermi_linewidth(float): Width of the Fermi level.
        title (str): The title of the plot.
        figsize (list): The figure size specified as \[width, height\].
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the highest level title and axis labels.
        \*\*kwargs: Other arguments passed to ``Axes.plot()`` of band plots.

    Returns:
        fig (Figure): Matplotlib figure object
    """
    import matplotlib.pyplot as plt
    from CRYSTALpytools.base.plotbase import plot_doss, _plot_label_preprocess
    from CRYSTALpytools.electronics import ElectronDOS
    from CRYSTALpytools.units import H_to_eV, eV_to_H
    import copy

    # unit
    if unit.lower() == 'ev':
        unit = 'eV'
        is_ev = True
    elif unit.lower() == 'a.u.':
        unit = 'a.u.'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    # instantiation and unit
    dossplt = []
    for id, d in enumerate(doss):
        if isinstance(d, str):
            dtmp = ElectronDOS.from_file(d)
        elif isinstance(d, ElectronDOS):
            dtmp = copy.deepcopy(d)
        else:
            raise ValueError('Unknown input type for doss.')

        if unit != dtmp.unit:
            dtmp._set_unit(unit)
            if unit == 'eV':
                fermi_level = H_to_eV(fermi_level)
            else:
                fermi_level = eV_to_H(fermi_level)
        dossplt.append(dtmp)

    # prj
    prj = np.array(prj, ndmin=1)
    if prj.ndim > 1: raise ValueError('Projections must be 1D.')
    if len(prj) == 0:
        prj = [[int(i+1) for i in range(len(j.doss))] for j in dossplt]
    else:
        prj = [prj for j in dossplt]

    # plot
    ndoss = len(dossplt)
    if ndoss == 1 and overlap == False: # Same system, projecton into different panels

        nprj = len(prj[0])
        fig, ax = plt.subplots(nprj, 1, figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        # prepare doss plot keywords
        ## label, same as overlap
        ## color, new default
        if np.all(dos_color==None):
            dos_color = [['tab:blue', 'tab:blue'] for i in range(nprj)]
        # Dimeonsion issue: dos plot styles must be consistent with length of input dosss
        dossref = [dossplt[0].doss[i-1] for i in prj[0]]
        commands = _plot_label_preprocess(
            dossref, dos_label, dos_color, dos_linestyle, dos_linewidth
        )
        for i in range(4):
            if np.all(commands[i]==None):
                commands[i] = [None for j in range(nprj)]
            else:
                commands[i] = [[commands[i][j]] for j in range(nprj)]
        for i in range(nprj):
            fig.axes[i] = plot_doss(
                fig.axes[i], dossplt[0].doss, dossplt[0].energy, beta, [prj[0][i]],
                energy_range, dos_range, commands[0][i], commands[1][i],
                commands[2][i], commands[3][i], fermi_level, fermi_color,
                fermi_linestyle, fermi_linewidth, legend, False, **kwargs
            )
    else: # Projecton of the same system into the same panel
        fig, ax = plt.subplots(ndoss, 1, figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        for i in range(ndoss):
            fig.axes[i] = plot_doss(
                fig.axes[i], dossplt[i].doss, dossplt[i].energy, beta, prj[i],
                energy_range, dos_range, dos_label, dos_color, dos_linestyle,
                dos_linewidth, fermi_level, fermi_color, fermi_linestyle,
                fermi_linewidth, legend, False, **kwargs
            )

    # set titles and axes
    if is_ev == True:
        fig.supylabel('DOS (states/eV)', fontsize=fontsize)
        fig.supxlabel('Energy (eV)', fontsize=fontsize)
    else:
        fig.supylabel('DOS (a.u.)', fontsize=fontsize)
        fig.supxlabel('Energy (a.u.)', fontsize=fontsize)

    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)

    return fig


#-----------------------------BAND + DENSITY OF STATES------------------------#

def plot_electron_banddos(
    *data, unit='eV', k_label=[], dos_beta='down', dos_overlap=True, dos_prj=[],
    energy_range=[], k_range=[], dos_range=[], band_width=2, band_label=None,
    band_color=None, band_linestyle=None, band_linewidth=None, dos_label=None,
    dos_color=None, dos_linestyle=None, dos_linewidth=None, fermi_level=0.,
    fermi_color='tab:green', fermi_linestyle='-', fermi_linewidth=1.0, title=None,
    figsize=[6.4, 4.8], legend='lower right', fontsize=14, **kwargs):
    """
    Plot electron band structure + dos for a **single** system, i.e., the
    ``bands`` and ``doss`` variables are not extendable.

    Input arguments not in the list are consistent with ``plot_electron_doss`` and
    ``plot_electron_bands``.

    Args:
        \*data: Either 1 or 2 entries. For one enetry, it is fort.25 containing
            both band and DOS, or ``ElectronBandDOS`` object. For 2 entries,
            the first entry is ``bands`` of ``plot_electron_bands`` and the
            second is ``doss`` of ``plot_electron_doss``
        dos_beta (str): ``beta`` of ``plot_electron_doss``.
        dos_overlap (bool): ``overlap`` of ``plot_electron_doss``. The user can
            either plot projections into the same subplot or into separate
            subplots.
        dos_prj (list): ``prj`` of ``plot_electron_doss``.
        band_width (int|float): Relative width of band structure, times of the
            width of a DOS subplot.
    Returns:
        fig (Figure): Matplotlib figure object

    :raise ValueError: If the unit parameter is unknown.
    """
    from CRYSTALpytools.electronics import ElectronBand, ElectronDOS, ElectronBandDOS
    from CRYSTALpytools.base.plotbase import plot_banddos
    from CRYSTALpytools.units import H_to_eV, eV_to_H
    import warnings, copy

    # unit
    if unit.lower() == 'ev':
        unit = 'eV'
        is_ev = True
    elif unit.lower() == 'a.u.':
        unit = 'a.u.'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    # instantiation and unit
    if len(data) == 1:
        if isinstance(data[0], str):
            bands = ElectronBand.from_file(data[0])
            doss = ElectronDOS.from_file(data[0])
        elif isinstance(data[0], ElectronBandDOS):
            bands = copy.deepcopy(data[0].band)
            doss = copy.deepcopy(data[0].dos)
        else:
            raise ValueError('Unknown input data type for the 1st entry.')
    elif len(data) == 2:
        if isinstance(data[0], str):
            bands = ElectronBand.from_file(data[0])
        elif isinstance(data[0], ElectronBand):
            bands = copy.deepcopy(data[0])
        else:
            raise ValueError('Unknown input data type for the 1st entry.')

        if isinstance(data[1], str):
            doss = ElectronDOS.from_file(data[1])
        elif isinstance(data[1], ElectronDOS):
            doss = copy.deepcopy(data[1])
        else:
            raise ValueError('Unknown input data type for the 2nd entry.')
    else:
        raise ValueError('Input parameter length does not meet requirements.')

    if doss.unit != bands.unit:
        warnings.warn("Band and DOS have different units. Fermi energy is assumed to follow the unit of band. If not 0, it might be wrong.",
                      stacklevel=2)

    if unit != doss.unit:
        doss._set_unit(unit)
    if unit != bands.unit:
        bands._set_unit(unit)
        if unit == 'eV':
            fermi_level = H_to_eV(fermi_level)
        else:
            fermi_level = eV_to_H(fermi_level)

    # set projections
    if len(dos_prj) == 0:
        dos_prj = [i+1 for i in range(len(doss.doss))]

    fig = plot_banddos(
        bands, doss, k_label, dos_beta, dos_overlap, dos_prj, energy_range,
        k_range, dos_range, band_width, band_label, band_color, band_linestyle,
        band_linewidth, dos_label, dos_color, dos_linestyle, dos_linewidth,
        fermi_level, fermi_color, fermi_linestyle, fermi_linewidth, figsize,
        legend, **kwargs
    )

    # set titles and axes
    if is_ev == True:
        if np.all(fermi_level!=0.0):
            fig.supylabel('Energy (eV)', fontsize=fontsize)
        else:
            fig.supylabel('$E-E_{F}$ (eV)', fontsize=fontsize)
    else:
        if np.all(fermi_level!=0.0):
            fig.supylabel('Energy (a.u.)', fontsize=fontsize)
        else:
            fig.supylabel('$E-E_{F}$ (a.u.)', fontsize=fontsize)

    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)

    return fig


##############################################################################
#                                                                            #
#                                    TOPOND                                  #
#                                                                            #
##############################################################################

#--------------------------------2D CONTOUR PLOT-----------------------------#

def plot_topond2D(*topond, unit='Angstrom', type='infer', option='normal',
                  levels='default', lineplot=True, linewidth=1.0, isovalues='%.4f',
                  colorplot=False, colormap='jet', cbar_label='default',
                  cpt_marker='o', cpt_color='k', cpt_size=10, traj_color='r',
                  traj_linestyle=':', traj_linewidth=0.5,
                  a_range=[0., 1.], b_range=[0., 1.], edgeplot=False,
                  x_ticks=5, y_ticks=5, layout=None, title=None,
                  figsize=[6.4, 4.8], sharex=True, sharey=True, fontsize=14):
    """
    Read and plot multiple TOPOND 2D plot files / objects. The uniform plot
    set-ups are used for comparison.

     .. note::

        For the convenience of analysis and plotting, it is important to select
        the correct type for your input file. By default `type='infer'` will
        search for (case insensitive) the following strings:

        'SURFRHOO', 'SURFSPDE', 'SURFLAPP', 'SURFLAPM', 'SURFGRHO', 'SURFKKIN',
        'SURFGKIN', 'SURFVIRI', 'SURFELFB', 'TRAJGRAD', 'TRAJMOLG'

        For their meanings, please refer the `TOPOND manual <https://www.crystal.unito.it/include/manuals/topond.pdf>`_.

    Available options:

    * 'normal' : Literally normal.  
    * 'diff' : Subtract data from the first entry using following entries. All
        the entries must have the same ``type`` and must be 'SURF*' types.  
    * 'overlay': Overlapping a 'TRAJ*' object on the 2D 'SURF*' object. Inputs
        must be 1\*2 lists of a ``Surf`` object and a ``Traj`` object. File
        names are not permitted as geometry information is mandatory.

    .. note::

        2D periodicity (``a_range`` and ``b_range``), though available for the
        child classes of the ``topond.ScalarField`` class, is not suggested as
        TOPOND plotting window does not always commensurate with periodic
        boundary. The ``topond.Trajectory`` class has no 2D periodicity so if
        ``option='overlay'``, ``a_range``, ``b_range`` and ``edgeplot`` will be
        disabled.

    .. note::

        The ``plot_lapm`` option is not available for ``Laplacian`` class. Use
        ``Laplacian().plot_2D()``.

    Args:
        \*topond (electronics.ChargeDensity|SpinDensity|Gradient|Laplacian|HamiltonianKE|LagrangianKE|VirialField|ELF|GradientTraj|ChemicalGraph|list|str):
            Extendable. File names, ``topond`` objects or 1\*2 list for
            ``option='overlay'``. Geometry information is not available if file
            names are used - **might lead to errors!**.
        unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
            Bohr :math:`^{-3}`.
        type (str): 'infer' or specified. Otherwise warning will be given.
        option (str): Available options see above.
        levels (array|int): Set levels of contour plot. 'Default' for built-in,
                property adaptive levels (``unit='Angstrom'``). Otherwise as
                array or int for linear scales. Entries **must be consistent
                with ``unit``**.
        lineplot (bool): Plot contour lines.
        linewidth (float): Contour linewidth. Useful only if ``lineplot=True``.
            Other properties are not editable.
        isovalues (str|None): Add isovalues to contour lines and set their
            formats. Useful only if ``lineplot=True``. None for not adding isovalues.
        colorplot (bool): Plot color-filled contour plots.
        colormap (str): Matplotlib colormap option. Useful only if ``colorplot=True``.
        cbar_label (str): Label of colorbar. Useful only if ``colorplot=True``.
            'None' for default.
        cpt_marker (str): Marker of critical point scatter.
        cpt_color (str): Marker color of critical point scatter.
        cpt_size (float|int): Marker size of critical point scatter.
        traj_color (str): Line color of 2D trajectory plot.
        traj_linestyl (str): Line style of 2D trajectory plot.
        traj_linewidth (str): Line width of 2D trajectory plot.
        a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in fractional
            coordinate. 'Surf' plots only.
        b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in fractional
            coordinate. 'Surf' plots only.
        edgeplot (bool): Whether to add cell boundaries represented by the
            original base vectors (not inflenced by a/b range).
        x_ticks (int): Number of ticks on x axis.
        y_ticks (int): Number of ticks on y axis.
        title (str|None): The title of the plot. 'None' for no title. The
            default subplot titles are used either way.
        figsize (list): Matplotlib figure size. Note that axes aspects are
            fixed to be equal.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the heightest level of title.

    Returns:
        fig (Figure): Matplotlib Figure object.
    """
    from CRYSTALpytools.topond import ScalarField, Trajectory
    from CRYSTALpytools.crystal_io import Properties_output
    import matplotlib.pyplot as plt
    import numpy as np
    import warnings

    obj = []
    for i in topond:
        if isinstance(i, str):
            obj.append(Properties_output().read_topond(i, type=type))
        elif isinstance(i, ScalarField) or isinstance(i, Trajectory):
            obj.append(i)
        elif isinstance(i, list) or isinstance(i, tuple):
            if len(i) != 2:
                raise ValueError('Input lists must have 2 elements.')
            if isinstance(i[0], ScalarField) and isinstance(i[1], Trajectory):
                obj.append([i[0], i[1]])
            elif isinstance(i[1], ScalarField) and isinstance(i[0], Trajectory):
                obj.append([i[1], i[0]])
            else:
                raise TypeError("Input type does not follow the requirement of the 'overlay' option.")
        else:
            raise TypeError("Input type does not meet the requirements.")

    # subtraction
    if 'diff' in option.lower():
        if isinstance(obj[0], Trajectory):
            raise TypeError("The 'diff' option is not applicable to 'topond.Trajectory' objects.")
        for i in obj[1:]:
            if isinstance(i, Trajectory): continue
            if obj[0].type != i.type:
                raise TypeError("Different properties are read for input objects / files, 'diff' option not available.")
            obj[0].subtract(i)
        obj = [obj[0]]

    # set uniform levels
    if np.all(levels=='default'):
        levels = 'default'
    else:
        levels = np.array(levels, ndmin=1, dtype=float)
        if levels.shape[0] == 1:
            ranges = []
            for i in obj:
                if isinstance(i, Trajectory): continue
                ranges.append([np.min(i.data), np.max(i.data)])
            levels = np.linspace(np.min(ranges), np.max(ranges), int(levels[0]))
        else:
            if levels.ndim > 1: raise ValueError('Input levels must be a 1D array.')

    # layout of plots
    nplt = len(obj)

    if np.all(layout!=None):
        if layout[0]*layout[1] < nplt:
            warnings.warn('Layout size smaller than the number of plots. Using default layout.',
                          stacklevel=2)
            layout = None
    if np.all(layout==None):
        if nplt == 1:
            layout = [1, 1]
        else:
            layout = [int(np.ceil(nplt/2)), 2]

    fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                           sharex=sharex, sharey=sharey, layout='tight')

    # plot
    ax_index = 0
    for o in obj:
        if isinstance(o, ScalarField):
            fig = o.plot_2D( # name the keywords for compatibility with Laplacian.plot_2D
                unit=unit, levels=levels, lineplot=lineplot, linewidth=linewidth,
                isovalues=isovalues, colorplot=colorplot, colormap=colormap,
                cbar_label=cbar_label, a_range=a_range, b_range=b_range,
                edgeplot=edgeplot, x_ticks=x_ticks, y_ticks=y_ticks,
                title='default', figsize=figsize, overlay=None, fig=fig,
                ax_index=ax_index)
        elif isinstance(o, Trajectory):
            fig = o.plot_2D(
                unit, cpt_marker, cpt_color, cpt_size, traj_color, traj_linestyle,
                traj_linewidth, x_ticks, y_ticks, 'default', figsize, None, fig,
                ax_index)
        else: # overlay plot
            if a_range[0] != 0. or a_range[1] != 1. \
            or b_range[0] != 0. or b_range[1] != 1.:
                warnings.warn("Periodic plotting not available for trajectory objects. Using default ranges.",
                              stacklevel=2)
                a_range_tmp = []; b_range_tmp = []
            fig = o[0].plot_2D( # name the keywords for compatibility with Laplacian.plot_2D
                unit=unit, levels=levels, lineplot=lineplot, linewidth=linewidth,
                isovalues=isovalues, colorplot=colorplot, colormap=colormap,
                cbar_label=cbar_label, a_range=a_range, b_range=b_range,
                edgeplot=edgeplot, x_ticks=x_ticks, y_ticks=y_ticks,
                title='default', figsize=figsize, overlay=o[1], fig=fig,
                ax_index=ax_index, cpt_marker=cpt_marker, cpt_color=cpt_color,
                cpt_size=cpt_size, traj_color=traj_color, traj_linestyle=traj_linestyle,
                traj_linewidth=traj_linewidth)
        ax_index += 1

    # title
    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)
    return fig


##############################################################################
#                                                                            #
#                                 RHOLINE (1D)                               #
#                                                                            #
##############################################################################

#-------------------------------------RHOLINE---------------------------------#


def plot_cry_rholine(rholine_obj):
    """
    Plot the resistivity as a function of distance.

    Args:
        rholine_obj (object): Rholine object containing the data for the resistivity.
        save_to_file (bool, optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
    - Plots the resistivity as a function of distance.
    - Sets the x-axis label as ``\\AA`` and the y-axis label as ``\\rho [\\frac{e}{\\AA^{3}]``.
    - Saves the plot to a file named 'figure_rholine_YYYY-MM-DD_HHMMSS.jpg' in the current directory.
    - If save_to_file is True, saves the plot to a file specified by save_to_file parameter.
    """
    import os
    import time

    import matplotlib.pyplot as plt

    plt.plot(rholine_obj.x, rholine_obj.y)

    plt.xlabel('d  [$\AA$]', fontsize=14)
    plt.ylabel(r'$\rho$  [$\frac{e}{\AA^3}$]', fontsize=16)

    path = os.path.join('./'+'figure_'+'rholine_' +
                        time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(rholine_obj.title, fontsize=15)
    plt.savefig(path, bbox_inches='tight', dpi=600)

    # if save_to_file != False:
    #     save_plot(save_to_file)

    plt.show()

#-----------------------------------LAPLACIAN---------------------------------#


def plot_cry_lapl_profile(lapl_obj):
    """
    Plot the Laplacian profile of a crystal.

    Args:
        lapl_obj (object): Laplacian object containing the data for the Laplacian profile.
        save_to_file (bool, optional): Indicates whether to save the plot to a file. Defaults to False.

    Returns:
        None

    Notes:
        - Plots the Laplacian profile using the data from the Laplacian object.
        - The x-axis represents the distance in angstroms.
        - The y-axis represents the Laplacian in electrons per cubic angstrom to the fifth power (e/A^5).
        - The area under the curve where the Laplacian is negative is filled with a light blue color.
        - The area under the curve where the Laplacian is positive is filled with a light coral color.
        - If save_to_file is set to a file path, the plot is saved to that file.
    """
    import time

    import matplotlib.pyplot as plt

    plt.plot(lapl_obj.datax, lapl_obj.datay)

    plt.fill_between(lapl_obj.datax, lapl_obj.datay, where=(
        lapl_obj.datay < 0), color='lightblue', interpolate=True)
    plt.fill_between(lapl_obj.datax, lapl_obj.datay, where=(
        lapl_obj.datay > 0), color='lightcoral', interpolate=True)

    # plt.xlim(-0.5,0.5)
    # plt.ylim(-200,200)

    plt.xlabel('Distance [A]')
    plt.ylabel('Laplacian [e/A^5]')

    # if save_to_file != False:
    #     save_plot(save_to_file)

    plt.show()

#-----------------------------DENSITY PROFILE---------------------------------#

def plot_cry_density_profile(lapl_obj):
    """
    Plot the density profile of a crystal.

    Args:
        lapl_obj (object): Laplacian object containing the data for the density profile.
        save_to_file (bool, optional): Indicates whether to save the plot to a file. Defaults to False.

    Returns:
        None

    Notes:
        - Plots the density profile using the data from the Laplacian object.
        - The x-axis represents the distance in angstroms.
        - The y-axis represents the density in electrons per cubic angstrom (e/A^3).
        - If save_to_file is set to a file path, the plot is saved to that file.
    """
    import time

    import matplotlib.pyplot as plt

    plt.plot(lapl_obj.datax, lapl_obj.datay)

    plt.xlabel('Distance [A]')
    plt.ylabel('Density [e/A^3]')

    # if save_to_file != False:
    #     save_plot(save_to_file)

    plt.show()

##############################################################################
#                                                                            #
#                             TRANSPORT PROPERTIES                           #
#                                                                            #
##############################################################################

#-----------------------------------TENSOR PROPS------------------------------#
def plot_transport_tensor(
    *boltztra, option='normal', x_axis='potential', x_range=[], direction='xx',
    spin='sum', plot_series=[], plot_label=None, plot_color=None,
    plot_linestyle=None, plot_linewidth=None, zero_color='tab:gray',
    zero_linestyle='-', zero_linewidth=1., layout=None, title=None,
    figsize=[6.4, 4.8], legend='upper left', sharey=True, fontsize=14, **kwargs):
    """
    Plot tensor-like transport properties in multiple ways:

    1. Normal: For 1 entry, same as ``transport.Tensor.plot()``. Otherwise plot
        different transport properties of the same system. A massive plot in
        nRow\*nCol is plotted and each subplot consists of nDir\*1 subplots. The
        ``layout`` option only specifies the layout of property blocks, and the
        ``sharey`` option fixes the y axis scale of 'sub-subplots'.  
    2. Multi: Plot transport properties of different systems together for
        comparison. In this case ``direction`` does not accept list variables,
        and entries of ``plot_series`` will be plotted into subplots.

    With list inputs the user can get thermoelectric power factor
    (:math:`S^{2}\\sigma`) or dimensionless figure of merit
    (:math:`ZT = \\frac{S^{r}\\sigma T}{\\kappa}`).

    * Power Factor: 1\*2 list of any 2 in 'SIGMA', 'SIGMAS' and 'SEEBECK'.  
    * ZT: 1\*3 list of 1 'KAPPA' and any 2 in 'SIGMA', 'SIGMAS' and 'SEEBECK'.

    .. note::

        Make sure all the entries are from the same calculation when using
        'power factor' or 'zt' options. The code only checks the
        dimensionalities of tensors.

    X-axis options:

    1. X_axis: Chemical potential :math:`\\mu`; Plot series: Temperature :math:`T`.  
    2. X_axis: Carrier density :math:`\\rho(\\mu; T)`; Plot series: Temperature :math:`T`.  
    3. X_axis: Temperature :math:`T`; Plot series: Chemical potential :math:`\\mu`.

    Args:
        \*boltztra (str|Tensor): DAT' files by CRYSTAL BOLTZTRA keyword or the
            `transport.Tensor` objects.
        option (str): 'power factor', 'zt' and 'multi'. See above.
        x_axis (str): X axis options, 'potential', 'carrier' or 'temperature'.
            See above.
        x_range (list): Display range of x axis. Y axis is self-adaptive.
        direction (str|list): Depending on the dimensionality of the system,
            including 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy' and 'zz'.
            A list of options will generate nDirect\*1 subplots. The direction
            of each subplot is annotated on the upper left corner. **List entry
            is not allowed for ``option='multi'``**.
        spin (str): Spin-polarized systems only. Electron spins to plot. 'up',
            'down' or 'sum'. Disabled for ``option='spin'``.
        plot_series (list|array|float): **Values** of the plot series. Should
            be commensurate with the choice of ``x_axis``. It will be annotated
            on the upper right corner for ``option='multi'``.
        plot_label (list|str|None): None for default values: ``plot_series``
            and its unit, or sequence number of materials for ``option='multi'``.
            If str, prefixing that entry in front of the default value. If list,
            it should be 1\*nPlotSeries, or 1\*nBoltztra for ``option='multi'``,
            for every plot line.
        plot_color (list|str|None): Similar to ``electronics.ElectronDOS.plot()``.
            If str, use the same color for all the plot lines. If
            1\*nPlotSeries(nBoltztra), use the same color for every plot line.
            If 2\*nPlotSeries(nBoltztra), use different colors for p-type and
            n-type carrier properties.
        plot_linestyle (list|str|None): Similar to ``electronics.ElectronDOS.plot()``.
            See explanations of ``plot_color``.
        plot_linewidth (list|float|None): Similar to ``electronics.ElectronDOS.plot()``.
            See explanations of ``plot_color``.
        zero_color (str): Color of the 0 value line.
        zero_linestyle (str): Linestyle of the 0 value line.
        zero_linewidth (float): Linewidth of the 0 value line.
        layout (list[int]): Layout of subplots. By default it is nDirection\*1
            or nPlot_series\*1 gird of figures.
        title (str): Title
        figsize (list): Figure size.
        legend (str|None): Location of legend. None for not adding legend.
        sharey (bool): Share y axis for multiple subplots. Share x is enforced.
        fontsize (float|int): Font size of the title, subplot capations
            (direction), x and y axis capations.
        \*\*kwargs: Other arguments passed to the matplotlib ``Axes.plot()``
                method. Applied to all the plots.
    Returns:
        fig (Figure): Matplotlib Figure object.
    """
    from CRYSTALpytools.transport import Tensor
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from CRYSTALpytools.base.plotbase import _plot_label_preprocess
    import warnings

    def check_list(b): # check elements of list entry
        for ib in range(len(b)):
            if isinstance(b[ib], str):
                b[ib] = Tensor.from_file(b[ib])
            elif isinstance(b[ib], Tensor):
                pass
            else:
                raise TypeError("List elements must be filename or classes from the 'transport.Tensor' class.")
        return b

    objs = []
    for b in boltztra:
        if isinstance(b, str):
            objs.append(Tensor.from_file(b))
        elif isinstance(b, Tensor):
            objs.append(b)
        elif isinstance(b, list) or isinstance(b, np.ndarray):
            if len(b) == 2:
                b = check_list(b)
                objs.append(Tensor.get_power_factor(b[0], b[1]))
            elif len(b) == 3:
                b = check_list(b)
                objs.append(Tensor.get_zt(b[0], b[1], b[2]))
            else:
                raise ValueError('For list entries, they must be 1D list of either 2 elements for power factor or 3 elements for ZT.')
        else:
            raise TypeError("Inputs must be filename, list or classes from the 'transport.Tensor' class.")

    if option.lower() != 'normal' and option.lower() != 'multi':
        raise ValueError("Unknown options: '{}'.".format(option))

    direction = np.array(direction, ndmin=1)
    if option.lower() == 'multi':
        if len(direction) != 1:
            raise ValueError("For 'multi' option, only one direction in string should be given.")
        direction = str(direction[0])
        for i in range(1, len(objs)):
            if objs[i].type != objs[0].type:
                raise TypeError("For 'multi' option, input entries must have the same type.")

    # normal options
    if option.lower() == 'normal':
        if len(objs) == 1: # use class method
            fig = objs[0].plot(
                x_axis, x_range, direction, spin, plot_series, plot_label,
                plot_color, plot_linestyle, plot_linewidth, zero_color,
                zero_linestyle, zero_linewidth, layout, title, figsize, legend,
                sharey, fontsize, None, None, **kwargs)
            return fig

        # multi properties
        ## layout
        nobjs = len(objs); ndir = len(direction)
        if np.all(layout!=None):
            if layout[0]*layout[1] < nobjs:
                warnings.warn("The specified layout is not sufficient to accommodate subplots. The default value is used.",
                              stacklevel=2)
                layout = None
        if np.all(layout==None):
            layout = [int(np.ceil(nobjs/2)), 2]

        layoutreal = [int(layout[0]*ndir), layout[1]]
        fig, ax = plt.subplots(layoutreal[0], layoutreal[1], sharex=True,
                               sharey=False, figsize=figsize, layout='tight')
        ## plot
        for iobj in range(nobjs):
            obj = objs[iobj]
            objcol = iobj % layout[1]; objrow = iobj // layout[1]
            ax_index = []
            for i in range(ndir):
                ax_index.append(
                    int(objrow * ndir * layout[1] + i * layout[1] + objcol)
                )
            fig = obj.plot(
                x_axis, x_range, direction, spin, plot_series, plot_label,
                plot_color, plot_linestyle, plot_linewidth, zero_color,
                zero_linestyle, zero_linewidth, None, 'default', figsize, legend,
                sharey, fontsize, fig, ax_index, **kwargs)

        if np.all(title!=None):
            fig.suptitle(title, fontsize=fontsize)
        return fig

    # option multi
    # subplots (plot_series)
    ## get common plot series, x range
    def find_common_row_elements(list2d):
        common_elements = list2d[0]
        for i in range(1, len(list2d)):
            common_elements = np.intersect1d(common_elements, list2d[i])
        return common_elements

    if x_axis == 'potential':
        series_tot = find_common_row_elements([i.T for i in objs])
        x_tot = np.array([[np.min(i.mu), np.max(i.mu)] for i in objs])
    elif x_axis == 'carrier':
        series_tot = find_common_row_elements([i.T for i in objs])
        x_tot = np.array([[np.min(np.abs(i.carrier)), np.max(np.abs(i.carrier))] for i in objs])
    elif x_axis == 'temperature':
        series_tot = find_common_row_elements([i.mu for i in objs])
        x_tot = np.array([[np.min(i.T), np.max(i.T)] for i in objs])
    else:
        raise ValueError("Unknown x axis value: '{}'.".format(x_axis))

    ## plot series values
    if plot_series == []:
        series = series_tot
    else:
        plot_series = np.array(plot_series, ndmin=1)
        series = []
        for p in plot_series:
            if len(np.where(series_tot==p)[0]) == 0:
                warnings.warn("The specified plot series value '{:6.2f}' does not appear in all the materials. It will be removed.".format(p),
                              stacklevel=2)
                continue
            else:
                series.append(p)

    nplt = len(series)
    if nplt == 0: raise Exception('Cannot find common plot series.')
    ## captions
    if x_axis == 'potential' or x_axis == 'carrier':
        captions = ['{:>5.0f} K'.format(i) for i in series]
    else:
        captions = ['{:>6.2f} eV'.format(i) for i in series]
    ## x range
    if x_range == []:
        x_range = [np.min(x_tot[:, 0]), np.max(x_tot[:, 1])]
    else:
        x_range = [np.min(x_range), np.max(x_range)]

    # direction
    if objs[0].data.shape[2] == 6:
        indices = {'xx' : 0, 'xy' : 1, 'xz' : 2, 'yx': 1, 'yy' : 3, 'yz' : 4,
                   'zx' : 2, 'zy' : 4, 'zz' : 5}
    else:
        indices = {'xx' : 0, 'xy' : 1, 'yx' : 1, 'yy' : 2}

    for o in objs:
        if objs[0].data.shape[2] != o.data.shape[2]:
            raise ValueError('Inconsistent dimensionalities of input materials. They must be all in 3D, 2D or 1D.')
    dir = indices[direction]

    # projections in the same plot
    nprj = len(objs)
    ## plot commands
    if np.all(plot_label==None):
        plot_label = ['# {:d}'.format(i+1) for i in range(nprj)]
    elif isinstance(plot_label, str):
        plot_label = ['{} {:d}'.format(plot_label, i+1) for i in range(nprj)]
    elif isinstance(plot_label, list) or isinstance(plot_label, np.ndarray):
        nlabel = len(plot_label)
        plot_label = [plot_label[i%nlabel] for i in range(nprj)]
    else:
        raise TypeError('Unknown type of plot label.')
    ## get a pseudo band input
    bands = np.zeros([nprj, 1, 1, 1], dtype=float)
    commands = _plot_label_preprocess(bands, plot_label, plot_color, plot_linestyle, plot_linewidth)

    # plot layout
    if np.all(layout!=None):
        if layout[0]*layout[1] < nplt:
            warnings.warn("The specified layout is not sufficient to accommodate subplots. The default value is used.",
                          stacklevel=2)
            layout = None
    if np.all(layout==None):
        layout = [nplt, 1]

    fig, ax = plt.subplots(layout[0], layout[1], sharex=True, sharey=sharey,
                           figsize=figsize, layout='tight')

    # plot every subplot
    y_range = []
    for iplt, ax in enumerate(fig.axes):
        y_range_plt = []
        ax.hlines(0, x_range[0], x_range[1], colors=zero_color, linestyle=zero_linestyle,
                  linewidth=zero_linewidth)
        # plot every material
        for iprj, obj in enumerate(objs):
            if x_axis == 'potential':
                ## get indices of projection
                iseries = np.where(obj.T==series[iplt])[0][0]
                if obj.spin == 1 or obj.lower() == 'up':
                    y = obj.data[iseries, :, dir, 0]
                    carrier = obj.carrier[iseries, :, 0]
                elif obj.spin == 2 and obj.lower() == 'sum':
                    y = np.sum(obj.data[iseries, :, dir, :], axis=1)
                    carrier = np.sum(obj.carrier[iseries, :, :], axis=1)
                else:
                    y = obj.data[iseries, :, dir, 1]
                    carrier = obj.carrier[iseries, :, 1]
                ## limit plot range
                idx_x = np.where((obj.mu>=x_range[0])&(obj.mu<=x_range[1]))[0]
                x = obj.mu[idx_x]
                y = y[idx_x]
                carrier = carrier[idx_x]
            elif x_axis == 'carrier':
                ax.set_xscale('log')
                ## get indices of projection
                iseries = np.where(obj.T==series[iplt])[0][0]
                if obj.spin == 1 or obj.lower() == 'up':
                    x = np.abs(obj.carrier[iseries, :, 0])
                    y = obj.data[iseries, :, dir, 0]
                    carrier = obj.carrier[iseries, :, 0]
                elif obj.spin == 2 and obj.lower() == 'sum':
                    x = np.abs(np.sum(obj.carrier[iseries, :, :], axis=1))
                    y = np.sum(obj.data[iseries, :, dir, :], axis=1)
                    carrier = np.sum(obj.carrier[iseries, :, :], axis=1)
                else:
                    x = np.abs(obj.carrier[iseries, :, 1])
                    y = obj.data[iseries, :, dir, 1]
                    carrier = obj.carrier[iseries, :, 1]
                ## limit plot range
                idx_x = np.where((x>=x_range[0])&(x<=x_range[1]))[0]
                x = x[idx_x]
                y = y[idx_x]
                carrier = carrier[idx_x]
            else:
                ## get indices of projection
                iseries = np.where(obj.mu==series[iplt])[0][0]
                if obj.spin == 1 or obj.lower() == 'up':
                    y = obj.data[:, iseries, dir, 0]
                    carrier = obj.carrier[:, iseries, 0]
                elif obj.spin == 2 and obj.lower() == 'sum':
                    y = np.sum(obj.data[:, iseries, dir, :], axis=1)
                    carrier = np.sum(obj.carrier[:, iseries, :], axis=1)
                else:
                    y = obj.data[:, iseries, dir, 1]
                    carrier = obj.carrier[:, iseries, 1]
                ## limit plot range
                idx_x = np.where((obj.T>=x_range[0])&(obj.T<=x_range[1]))[0]
                x = obj.T[idx_x]
                y = y[idx_x]
                carrier = carrier[idx_x]

            ## divide the plot by p and n type carriers and plot
            idx_p = np.where(carrier>=0)[0]
            idx_n = np.where(carrier<0)[0]
            if len(idx_p) == 0 and len(idx_n) == 0:
                raise ValueError('Empty data in the specified x range. Check your input.')
            if len(idx_p) > 0:
                ax.plot(x[idx_p], y[idx_p], label=commands[0][iprj][0],
                        color=commands[1][iprj][0], linestyle=commands[2][iprj][0],
                        linewidth=commands[3][iprj][0], **kwargs)
            if len(idx_n) > 0:
                 ax.plot(x[idx_n], y[idx_n], color=commands[1][iprj][1],
                        linestyle=commands[2][iprj][1], linewidth=commands[3][iprj][1],
                        **kwargs)
            ## linearly interpolate the gap, noticed when x_axis = 'potential' only
            if len(idx_p) > 0 and len(idx_n) > 0 and x_axis == 'potential':
                lastppt = [x[idx_p[-1]], y[idx_p[-1]]]
                firstnpt = [x[idx_n[0]], y[idx_n[0]]]
                midpt = [(lastppt[0] + firstnpt[0]) / 2, (lastppt[1] + firstnpt[1]) / 2]
                ax.plot([lastppt[0], midpt[0]], [lastppt[1], midpt[1]],
                        color=commands[1][iprj][0], linestyle=commands[2][iprj][0],
                        linewidth=commands[3][iprj][0], **kwargs)
                ax.plot([midpt[0], firstnpt[0]], [midpt[1], firstnpt[1]],
                        color=commands[1][iprj][1], linestyle=commands[2][iprj][1],
                        linewidth=commands[3][iprj][1], **kwargs)

            y_range_plt.append([np.min(y), np.max(y)])
        y_range.append([np.min(y_range_plt), np.max(y_range_plt)])

    # plot setups
    for iplt, ax in enumerate(fig.axes):
        if np.all(legend!=None) and np.all(commands[0][0][0]!=None) and iplt==0:
            ax.legend(loc=legend) # add legend to the first plot only
        ax.set_xlim(x_range)
        if sharey == False:
            ax.text(x_range[1], y_range[iplt][1], captions[iplt], fontsize=fontsize,
                    horizontalalignment='right', verticalalignment='top')
            ax.set_ylim(y_range[iplt])
        else:
            y_range = [np.min(y_range), np.max(y_range)]
            ax.text(x_range[1], y_range[1], captions[iplt], fontsize=fontsize,
                    horizontalalignment='right', verticalalignment='top')
            ax.set_ylim(y_range)

    fig.supylabel('{} ({})'.format(objs[0].type.capitalize(), objs[0].unit),
                  fontsize=fontsize)
    if x_axis == 'potential':
        fig.supxlabel('Chemical Potential (eV)', fontsize=fontsize)
    elif x_axis == 'carrier':
        fig.supxlabel('Carrier Density (cm$^{-3}$)', fontsize=fontsize)
    else:
        fig.supxlabel('Temperature (K)', fontsize=fontsize)
    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)

    return fig

##############################################################################
#                                                                            #
#                             ELASTIC PROPERTIES                             #
#                                                                            #
##############################################################################
#--------------------------------3D ELASTIC----------------------------------#
def plot_elastics3D(
    *tensor, property, uniform_scale=True, nphi=90, ntheta=90, nchi=180,
    scale_radius=True, range_cbar=None, range_x=None, range_y=None, range_z=None,
    u=None, utext=None, use_cartesian=True, plot_lattice=False, colormap='jet',
    layout=None, title=None, figsize=[6.4, 4.8], fontsize=14, **kwargs):
    """
    A wrapper function of :ref:`Tensor3D <ref-elastics>` objects to plot 3D
    crystal elastic properties. The user can plot multiple properties for
    different systems. The function returns to a matplotlib figure object for
    further processing. The uniform plot set-ups are used for comparison. Only
    matplotlib is used for plotting. Plotly is not available.

    Properties:

    * "young": Young's modulus.
    * "comp": Compressibility.
    * "shear avg": Average shear modulus.
    * "shear min": Minimum shear modulus.
    * "shear max": Maximum shear modulus.
    * "poisson avg": Average Poisson ratio.
    * "poisson min": Minimum Poisson ratio.
    * "poisson max": Maximum Poisson ratio.

    Args:
        \*tensor (str | Tensor3D | numpy.ndarray): Elastic tensor definition.
            Can be CRYSTAL output files, ``Tensor3D`` objects and 6\*6
            **elastic** matrices in Voigt notation, GPa. For files,
            ``conventional_lattice=True``.
        property (str | list[str]): The properties to plot. See above.
        uniform_scale (bool): Use the same color scale for all plots of the
            same property.
        nphi (int): Resolution of azimuth angle :math:`\\phi` on xy plane, in
            radian :math:`[0, 2\\pi)`.
        ntheta (int): Resolution of polar angle :math:`\\theta` in radian
            :math:`[0, 2\\pi)`. In practice only half of the points defined in
            :math:`[0, \\pi)` are used.
        nchi (int): Resolution of auxiliary angle  :math:`\\chi`, in radian
            :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
        scale_radius (bool): To scale the radius by values of the elastic
            property, or plot a sphere with radius = 1.
        range_cbar, range_x, range_y, range_z (list[float,float]): *Not
            suggested* Explicitly specifying the ranges of colorbar, x, y and z
            axes.
        u (numpy.ndarray): 3\*1 or nu\*3 array of vectors to be added into the
            figure. A line segment with doubled radius is plotted to indicate
            the vector. Or 'max' / 'min' / 'bothends', plot vectors
            corresponding to the max and min of elastic properties. For
            'bothends', min first.
        utext (list[str] | str): A string or a list of string. Used to mark the
            vector. If ``None``, the input ``u`` and the corresponding value is
            annotated. If 'value', the value is annotated.
        use_cartesian (bool): Vector is defined as cartesian or fractional
            coordinates. (*Only when lattice information is available*.)
        plot_lattice (bool): Draw the lattice box around the 3D surface.
        colormap (str): Colormap name.
        layout (list|tuple): The layout of subplots, \[nrow, ncol\]. Default is
            nTensor\*nProperty or 1\*nTensor.
        title (str|None): The title of the plot. 'None' for no title. The
            default subplot titles (property) are added either way.
        figsize (list): Matplotlib figure size.
        fontsize (int): Fontsize of the heightest level of title.
        \*\*kwargs: Parameters passed to ``Axes3D.view_init``. Only camera
            position keywords are suggested.

    Returns:
        fig (Figure): Matplotlib figure object.
    """
    from CRYSTALpytools.elastics import Tensor3D, _plot3D_mplib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D, axes3d
    import numpy as np
    import warnings

    # sanity check and preparation
    tensplt = []
    for t in tensor:
        if isinstance(t, str):
            tensplt.append(Tensor3D.from_file(t))
        elif isinstance(t, Tensor3D):
            tensplt.append(t)
        else:
            ttmp = np.array(t, dtype=float)
            if np.shape(ttmp)[0] != 6 or np.shape(ttmp)[1] != 6:
                raise ValueError('Input tensor is not a 6x6 matrix in Voigt notation.')
            tensplt.append(Tensor3D(matrix=ttmp, lattice=None, is_compliance=False))

    property = np.array(property, ndmin=1)
    if len(tensplt) == 1 and len(property) > 1 and uniform_scale == True:
        warnings.warn("'uniform_scale' cannot be used for multiple proeprties of the same system. Using 'uniform_scale = False'.",
                      stacklevel=2)
        uniform_scale = False

    # plot layout
    n_plot = len(tensplt) * len(property)
    if np.all(layout!=None):
        if layout[0]*layout[1] < n_plot :
            warnings.warn('Insufficient layout. Using default values.', stacklevel=2)
            layout = None
    if np.all(layout==None):
        if len(property) == 1: layout = [1, len(tensplt)]
        else: layout = [len(tensplt), len(property)]

    fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                           layout='tight', subplot_kw={'projection' : '3d'})

    if uniform_scale == False: # Non-uniform scale
        ax_index = 0
        for ip, p in enumerate(property):
            for it, t in enumerate(tensplt):
                t.plot_3D(
                    p, nphi, ntheta, nchi, scale_radius, range_cbar, range_x,
                    range_y, range_z, u, utext, use_cartesian, plot_lattice,
                    colormap, None, [6.4, 4.8], None
                )
                fig = _plot3D_mplib(
                    t, range_cbar, range_x, range_y, range_z, colormap,
                    figsize, fig, ax_index, None, **kwargs
                )
                if 'comp' not in p.lower() and 'poisson' not in p.lower():
                    fig.axes[ax_index].set_title('{} (GPa)'.format(p))
                else:
                    fig.axes[ax_index].set_title(p)
                ax_index += 1
    else: # Uniform scale
        ax_index = 0
        for ip, p in enumerate(property):
            tens_all = []; Rrange_all = []
            # get data first
            for it, t in enumerate(tensplt):
                tens_all.append(t.plot_3D(
                    p, nphi, ntheta, nchi, scale_radius, range_cbar, range_x,
                    range_y, range_z, u, utext, use_cartesian, plot_lattice,
                    colormap, None, [6.4, 4.8], None
                ))
                Rrange_all.append([np.min(tens_all[it].R),
                                   np.max(tens_all[it].R)])
            # plot
            Rmax = np.max(Rrange_all); Rmin = np.min(Rrange_all)
            if np.all(range_cbar==None):
                range_cbar = [Rmin, Rmax]
            for it, t in enumerate(tens_all):
                fig = _plot3D_mplib(
                    t, range_cbar, range_x, range_y, range_z, colormap,
                    figsize, fig, ax_index, Rmax, **kwargs
                )
                if 'comp' in p.lower() or 'poisson' in p.lower():
                    fig.axes[ax_index].set_title(p)
                else:
                    fig.axes[ax_index].set_title('{} (GPa)'.format(p))
                ax_index += 1

    if np.all(title!=None):
        fig.suptitle_title(title)
    return fig


#--------------------------------2D ELASTIC----------------------------------#
def plot_elastics2D(
    *tensor, property, plane=[], ntheta=90, nchi=180, plane_definition='miller',
    u=None, utext=None, use_cartesian=True, plot_lattice=True, loop_label=None,
    loop_color=None, loop_linestyle=None, loop_linewidth=None, layout=None,
    title=None, figsize=[6.4, 4.8], legend='upper left', fontsize=14, **kwargs):
    """
    A wrapper function of :ref:`Tensor3D or Tensor2D <ref-elastics>` objects to
    plot 2D crystal elastic properties. The user can plot multiple properties on
    multiple crystal planes (3D only) for multiple systems. The function returns
    to a matplotlib figure object for further processing. The uniform plot
    set-ups (radius) are used for comparison. Base units: GPa, m.

    Properties, depending on the dimensionality of systems:

    * "young": Young's modulus.  
    * "comp": Compressibility.  
    * "shear": Shear modulus between vectors in plane and plane norm (3D) or in-plane vectors (2D).  
    * "shear avg": Average shear modulus (3D).  
    * "shear min": Minimum shear modulus (3D).  
    * "shear max": Maximum shear modulus (3D).  
    * "poisson": Poisson ratio between vectors in plane and plane norm or in-plane vectors (2D).  
    * "poisson avg": Average Poisson ratio (3D).  
    * "poisson min": Minimum Poisson ratio (3D).  
    * "poisson max": Maximum Poisson ratio (3D).

    For 3D systems, the plotting planes have to be defined with any of the
    following methods:

    * 'miller': Miller indices.  
    * 'cartesian': The plane normal vector in Cartesian coordinates.  
    * 'fractional': The plane normal vector in fractional coordinates.

    The default layout is either 1\*nProp or nProp\*nPlane. For multi-system
    plottings, the properties on the same plane are plotted into the same axis
    for comparison. The radius is 'uniform' for the same property, i.e., same
    raidus for plots on various materials and planes.

    .. note::

        For multi-system plotting, the ``loop_*`` inputs do not change with
        the plotting planes, but with systems in order to distinguish them.
        Therefore they should be defined by 1\*nTensor list. In this case,
        ``loop_label`` is used. Otherwise it is not called.

        For multi-system plotting, ``u`` and ``utext`` options are disabled.

    Args:
        \*tensor (str|Tensor3D|numpy.ndarray): Elastic tensor definition. Can be
            CRYSTAL output files, ``Tensor3D`` / ``Tensor2D``objects and 6\*6 /
            3\*3 **elastic** matrices in Voigt notation, GPa. But
            dimensionalities of systems must be consistent. For 3D files,
            ``conventional_lattice=True``.
        property (str | list[str]): The properties to plot. See above.
        plane (numpy.ndarray): *3D only* 3\*1 or nplane\*3 array of planes.
        ntheta (int): Resolution of azimuth angle :math:`\\theta` on plane, in
            radian :math:`[0, 2\\pi)`.
        nchi (int): *3D only* Resolution of auxiliary angle  :math:`\\chi`, in
            radian :math:`[0, 2\\pi)`. Shear modulus and Poisson ratio only.
        plane_definition (str): *3D only* The method used to define the plane.
            See options above.
        u (numpy.ndarray): nDimen\*1 or nu\*nDimen array of vectors to be added
            into the figure. A line segment with doubled radius is plotted to
            indicate the vector. Or 'max' / 'min' / 'bothends', plot vectors
            corresponding to the max and min of elastic properties. For
            'bothends', min first.
        utext (list[str] | str): A string or a list of string. Used to mark
            The vector. If ``None`` or 'value', the value is annotated.
        use_cartesian (bool): Vector is defined as cartesian or fractional
            coordinates. (*Only when lattice information is available*.)
        plot_lattice (bool): Draw lattice base vectors to indicate orientation
            of the plane.
        loop_label (str|list): *Multi-tensor only* Label of loops. None for
            default values: sequence number of materials If str, prefixing that
            entry in front of the default value. If list, it should be in
            1\*nTensor.
        loop_color (str|list): Color of 2D loops. If only one string is given,
            apply it to all loops. 'None' for default values (matplotlib
            Tableau Palette).
        loop_linestyle (str|list): Linestyle of 2D loops. If only one string is
            given, apply it to all plots. 'None' for default values ('-').
        loop_linewidth (str|list): Linewidth of band structure. If only one
            number is given, apply it to all plots. For the 'multi' mode,
            1\*nSystem list of linewidth numbers. 'None' for default values (1.0).
        layout (list|tuple): *2D only* 2\*1 list. The layout of subplots. The
            first element is nrows, the second is ncolumns. By default, 1\*nProp.
        title (str): Title
        figsize (list): The figure size specified as \[width, height\].
        legend (str|None): *Multi-tensor only* Location of legend. None for not
            adding legend.
        fontsize (float|int): Font size of the title.
        \*\*kwargs: Parameters passed to ``PolarAxes.plot`` to customize the
            loops. Applied to all the loops.

    Returns:
        fig (Figure): Matplotlib figure object.
    """
    from CRYSTALpytools.elastics import Tensor3D, Tensor2D, tensor_from_file, _plot2D
    from CRYSTALpytools.base.plotbase import _plot_label_preprocess
    import numpy as np
    import warnings
    import matplotlib.pyplot as plt

    # sanity check and preparation
    tensplt = []
    for t in tensor:
        if isinstance(t, str):
            tensplt.append(tensor_from_file(t))
        elif isinstance(t, Tensor3D) or isinstance(t, Tensor2D):
            tensplt.append(t)
        else:
            ttmp = np.array(t, dtype=float)
            if np.shape(ttmp)[0] == 6 and np.shape(ttmp)[1] == 6:
                tensplt.append(Tensor3D(matrix=ttmp, lattice=None, is_compliance=False))
            elif np.shape(ttmp)[0] == 3 and np.shape(ttmp)[1] == 3:
                tensplt.append(Tensor2D(matrix=ttmp, lattice=None, is_compliance=False))
            else:
                raise ValueError('Input tensor is not a 6x6 or 3x3 matrix in Voigt notation.')

    if isinstance(tensplt[0], Tensor2D):
        ndim = 2
    else:
        ndim = 3
        if plane == []: raise ValueError("For 3D systems, plot plane must be befined.")

    for t in tensplt:
        if (ndim == 2 and isinstance(t, Tensor3D)) \
        or (ndim == 3 and isinstance(t, Tensor2D)):
            raise TypeError('The dimensionalities of input tensors must be consistent.')
    ntensor = len(tensplt)

    property_list_2D = ['young', 'comp', 'shear', 'poisson']
    property_list_3D = ['young', 'comp', 'shear', 'shear avg', 'shear min',
                        'shear max', 'poisson', 'poisson avg', 'poisson min',
                        'poisson max']
    property = np.array(property, ndmin=1)
    nprop = len(property)
    for ip in range(nprop):
        property[ip] = property[ip].lower()
        if (ndim == 2 and property[ip] not in property_list_2D) \
        or (ndim == 3 and property[ip] not in property_list_3D):
                raise ValueError("Unknown property input: '{}'.".format(property[ip]))

    # Single system, use class methods
    if ntensor == 1:
        if ndim == 2:
            fig = tensplt[0].plot_2D(
                property, ntheta, u, utext, use_cartesian, plot_lattice,
                loop_color, loop_linestyle, loop_linewidth, layout, figsize ,
                False, **kwargs)
        else:
            fig = tensplt[0].plot_2D(
                property, plane, ntheta, nchi, plane_definition, u, utext,
                use_cartesian, plot_lattice, loop_color, loop_linestyle,
                loop_linewidth, figsize, False, **kwargs)

        if np.all(title!=None):
            fig.suptitle(title, fontsize=fontsize)
        return fig

    # Multi-system
    if ndim == 2:
        nplane = 1
    else:
        plane = np.array(plane, ndmin=2)
        nplane = plane.shape[0]

    if np.all(u!=None) or np.all(utext!=None):
        warnings.warn("For multi-system plottings, 'u' and 'utext' options are disabled.",
                      stacklevel=2)
    u = None; utext = None
    ## calculate data
    for it, t in enumerate(tensplt):
        if ndim == 2:
            tensplt[it] = t.plot_2D(
                property, ntheta, u, utext, use_cartesian, plot_lattice,
                loop_color, loop_linestyle, loop_linewidth, layout, figsize ,
                True, **kwargs)
        else:
            tensplt[it] = t.plot_2D(
                property, plane, ntheta, nchi, plane_definition, u, utext,
                use_cartesian, plot_lattice, loop_color, loop_linestyle,
                loop_linewidth, figsize, True, **kwargs)

    ## scale of properties
    rmax = np.zeros([nprop, ntensor], dtype=float)
    for it, t in enumerate(tensplt):
        for ip in range(nprop):
            rmax[ip, it] = np.max(t.r[ip])
    rmax = np.max(rmax, axis=1)

    ## plot setups
    if np.all(loop_label==None):
        loop_label = ['# {:d}'.format(i+1) for i in range(ntensor)]
    elif isinstance(loop_label, str):
        loop_label = ['{} {:d}'.format(i+1) for i in range(ntensor)]
    elif isinstance(loop_label, list) or isinstance(loop_label, np.ndarray):
        nlabel = len(loop_label)
        loop_label = [loop_label[i%nlabel] for i in range(ntensor)]
    else:
        raise TypeError('Unknown type of plot label.')

    ## dummy doss obj
    doss = np.zeros([ntensor, 1, 1])
    commands = _plot_label_preprocess(doss, loop_label, loop_color,
                                      loop_linestyle, loop_linewidth)

    ## plot system 1
    fig = _plot2D(
        tensplt[0], rmax, use_cartesian, u, utext, commands[0][0][0],
        commands[1][0][0], commands[2][0][0], commands[3][0][0], layout, figsize,
        None, None, **kwargs)
    ## other systems, no lattice
    ax_index = [i for i in range(len(fig.axes))]
    for it, t in enumerate(tensplt[1:]):
        t.lattice_plot = None
        fig = _plot2D(
            t, rmax, use_cartesian, u, utext, commands[0][it+1][0],
            commands[1][it+1][0], commands[2][it+1][0], commands[3][it+1][0],
            layout, figsize, fig, ax_index, **kwargs
        )
    ## legend
    if np.all(legend!=None): fig.axes[0].legend(loc=legend)

    ## titles
    if nplane == 1 and nprop > 1:
        for ip, p in enumerate(tensplt[0].property):
            if 'shear' in p or 'young' in p: ptitle = '{} (GPa)'.format(p)
            else: ptitle = p
            fig.axes[ip].set_title(ptitle, y=1.05)

        if ndim == 3:
            if plane_definition == 'miller':
                axtitle = '({:<3d}{:^3d}{:>3d})'.format(plane[0,0], plane[0,1], plane[0,2])
            else:
                axtitle = '({:<4.1f}{:^4.1f}{:>4.1f})'.format(plane[0,0], plane[0,1], plane[0,2])

            fig.axes[0].text(
                -0.1, 0.5, axtitle, rotation='vertical', horizontalalignment='right',
                verticalalignment='center', transform=fig.axes[0].transAxes,
                fontsize=fig.axes[0].title._fontproperties._size
            )
    else: # 3D only
        for inorm, norm in enumerate(plane):
            if plane_definition == 'miller':
                axtitle = '({:<3d}{:^3d}{:>3d})'.format(norm[0], norm[1], norm[2])
            else:
                axtitle = '({:<4.1f}{:^4.1f}{:>4.1f})'.format(norm[0], norm[1], norm[2])
            fig.axes[inorm].set_title(axtitle, y=1.05)

        for ip, p in enumerate(tensplt[0].property):
            if 'shear' in p or 'young' in p: ptitle = '{} (GPa)'.format(p)
            else: ptitle = p
            fig.axes[int(ip*nplane)].text(
                -0.1, 0.5, ptitle, rotation='vertical', horizontalalignment='right',
                verticalalignment='center', transform=fig.axes[int(ip*nplane)].transAxes,
                fontsize=fig.axes[0].title._fontproperties._size
            )

    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)
    return fig


##############################################################################
#                                                                            #
#                             VIBRATIONAL PROPERTIES                         #
#                                                                            #
##############################################################################

#-----------------------------PHONON BAND AND DOS----------------------------#

def plot_phonon_bands(*bands, unit='cm-1', q_overlap_tol=1e-4,
                      k_label=None, mode='single', not_scaled=False,
                      frequency_range=[], k_range=[], band_label=None,
                      band_color=None, band_linestyle=None, band_linewidth=None,
                      plot_freq0=True, freq0_color='tab:gray', freq0_linestyle='-',
                      freq0_linewidth=1.0, layout=None, title=None,
                      figsize=[6.4, 4.8], legend='lower right', sharex=True,
                      sharey=True, fontsize=14, **kwargs):
    """
    Plot phonon band structures.

    Args:
        \*bands (PhononBand|str): ``phonons.PhononBand`` object or the standard
            screen output (.out) files of CRYSTAL.
        unit (str): The unit of frequency. Can be 'cm-1' or 'THz'.
        q_overlap_tol (float): The threshold for overlapped k points. Only used
             for getting tick positions. *For filename input only.*
        k_label (list): nSystem\*nTick or 1\*nTick list of strings of the label
             for high symmetry points along the path. If a 1D list is given,
             the same labels are used for all the systems. `mathtext <https://matplotlib.org/stable/users/explain/text/mathtext.html>`_
             experssions can also be used  as in matplotlib.
        mode (str): The plotting mode, 'single', 'multi' and 'compare'.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
            Useful with ``mode='compare'``. The multi mode forces scaling.
        frequency_range (array): A 2x1 array specifying the frequency range.
        k_range (array): A 2x1 array specifying the k-range.
        band_label (str|list): Plot legend. If only one string is given, apply
            it to all plots. 1\*nSystem plot legend otherwise.
        band_color (str|list): Color of band structure. If only one string is
            given, apply it to all plots. For the 'multi' mode, 1\*nSystem list
            of plot colors. 'None' for default values (matplotlib Tableau
            Palette).
        band_linestyle (str|list): Linestyle of band structure. If only one
            string is given, apply it to all plots. For the 'multi' mode,
            1\*nSystem list of linestyle strings. 'None' for default values
            ('-').
        band_linewidth (str|list): Linewidth of band structure. If only one
            number is given, apply it to all plots. For the 'multi' mode,
            1\*nSystem list of linewidth numbers. 'None' for default values
            (1.0).
        plot_freq0 (bool): Whether to plot 0 frequency line.
        freq0_color (str): Color of the 0 frequency line.
        freq0_linestyle (str): Line style of 0 frequency line.
        freq0_linewidth(float): Width of the 0 frequency line.
        layout (list|tuple): For 'compare' mode, the layout of subplots,
            \[nrow, ncol\]. The default is 2 cols per row.
        title (str): The title of the plot.
        figsize (list): The figure size specified as \[width, height\].
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the highest level title and axis labels.
        \*\*kwargs: Other arguments passed to ``Axes.plot()`` of band plots.

    Returns:
        fig (Figure): Matplotlib figure object

    :raise ValueError: If the specified unit is unknown.
    """
    import matplotlib.pyplot as plt
    from CRYSTALpytools.phonons import PhononBand
    from CRYSTALpytools.base.plotbase import plot_overlap_bands, plot_compare_bands
    import copy

    # unit
    if unit.lower() == 'cm-1':
        unit = 'cm-1'
        is_thz = False
    elif unit.lower() == 'thz':
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    # instantiation and unit
    bandsplt = []
    for ib, b in enumerate(bands):
        if isinstance(b, str):
            btmp = PhononBand.from_file(b, q_overlap_tol=q_overlap_tol)
        elif isinstance(b, PhononBand):
            btmp = copy.deepcopy(b)
        else:
            raise ValueError('Unknown input type for bands.')

        if unit != btmp.unit:
            btmp._set_unit(unit)
        bandsplt.append(btmp)

    # k label
    if np.all(k_label==None):
        k_label = [b.tick_label for b in bandsplt]

    # frequency = 0
    if plot_freq0 != True:
        freq0 = None
    else:
        freq0 = 0.

    # plot mode
    if mode.lower() == 'single':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = plot_compare_bands(
            [ax], [bandsplt[0].bands], [bandsplt[0].k_path], [bandsplt[0].tick_pos],
            k_label, False, frequency_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, freq0, freq0_color, freq0_linestyle,
            freq0_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'multi':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        k_xax = [i.k_path for i in bandsplt]
        ax = plot_overlap_bands(
            ax, [b.bands for b in bandsplt], k_xax, [b.tick_pos for b in bandsplt],
            k_label, frequency_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, freq0, freq0_color, freq0_linestyle,
            freq0_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'compare':
        if np.all(layout==None):
            layout = [int(np.ceil(len(bandsplt)/2)), 2]

        fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        k_xax = [i.k_path for i in bandsplt]
        _ = plot_compare_bands(
            ax.flat, [b.bands for b in bandsplt], k_xax, [b.tick_pos for b in bandsplt],
            k_label, not_scaled, frequency_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, freq0, freq0_color, freq0_linestyle,
            freq0_linewidth, legend, **kwargs
        )
    else:
        raise ValueError("Unknown mode input: '{}'.".format(mode))

    # set titles and axes
    if is_thz == True:
        fig.supylabel('Frequency (THz)', fontsize=fontsize)
    else:
        fig.supylabel('Frequency (cm$^{-1}$)', fontsize=fontsize)
    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)

    return fig


def plot_phonon_doss(
    *doss, unit='cm-1', read_INS=False, atom_prj=[], element_prj=[],
    overlap=False, prj=[], gauss=0.0, frequency_range=[], dos_range=[],
    dos_label=None, dos_color=None, dos_linestyle=None, dos_linewidth=None,
    plot_freq0=True, freq0_color='tab:gray', freq0_linestyle='-',
    freq0_linewidth=1.0, title=None, figsize=[6.4, 4.8], legend='lower right',
    sharex=True, sharey=False, fontsize=14, **kwargs):
    """
    Plot phonon density of states.

    Args:
        \*doss (PhononDOS|str): ``phonon.PhononDOS`` object or the standard
            screen output (.out) files of CRYSTAL.
        unit (str): The unit of frequency. Can be 'cm-1' or 'THz'.
        read_INS (bool): Read the inelastic neutron scattering spectra. *For
            filename input only.*
        atom_prj (list): Read the projections of atoms with specified labels.
            *For filename input only.*
        element_prj (list): Read projections of elements with specified
            conventional atomic numbers. *For filename input only.*
        overlap (bool): Plotting multiple projections into the same figure.
            Useful only if a single entry of ``doss`` is plotted. Otherwise
            projections from the same entry will be overlapped into the same
            subplot.
        prj (list): Index of selected projections, consistent with the first
            dimension of the ``doss``, starting from 1. Effective for all the
            subplots.
        gauss (float): The standard deviation of Gaussian broadening, i.e. the
            :math:`\\sigma` of :math:`a\\exp{\\frac{(x-b)^{2}}{2\\sigma^{2}}}`.
            '0' for no Gaussian broadening. The length of data will be tripled.
            Valid only for data on the plot. Data saved in object is unchanged.
        frequency_range (list): 1\*2 list of frequency range
        dos_range (list): 1\*2 list of DOS range
        dos_label (str|list): Plot legend. If only one string is given, apply
            it to all plots. 1\*nPrj plot legend otherwise. Effective for all
            the subplots.
        dos_color (str|list): Color of DOSS plots. If only one string is given,
            apply it to all plots. When ``overlap=True``, 1\*nPrj list of plot
            color. 'None' for default values (matplotlib Tableau Palette).
            Effective for all the subplots.
        dos_linestyle (str|list): Linestyle of DOSS plot. If only one string is
            given, apply it to all plots. When ``overlap=True``, 1\*nPrj list
            of line styles. 'None' for default values ('-'). Effective for all
            the subplots.
        dos_linewidth (str|list): Linewidth of DOSS plot. If only one number is
            given, apply it to all plots. When ``overlap=True``, 1\*nPrj list
            of line widthes. 'None' for default values (1.0). Effective for all
            the subplots.
        plot_freq0 (bool): Whether to plot 0 frequency line.
        freq0_color (str): Color of the 0 frequency line.
        freq0_linestyle (str): Line style of the 0 frequency line.
        freq0_linewidth(float): Width of the 0 frequency line.
        title (str): The title of the plot.
        figsize (list): The figure size specified as \[width, height\].
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        fontsize (int): Fontsize of the highest level title and axis labels.
        \*\*kwargs: Other arguments passed to ``Axes.plot()`` of band plots.

    Returns:
        fig (Figure): Matplotlib figure object
    """
    import matplotlib.pyplot as plt
    from CRYSTALpytools.base.plotbase import plot_doss, _plot_label_preprocess
    from CRYSTALpytools.phonons import PhononDOS
    import copy, warnings

    # unit
    if unit.lower() == 'cm-1':
        unit = 'cm-1'
        is_thz = False
    elif unit.lower() == 'thz':
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    # instantiation and unit
    dossplt = []
    for id, d in enumerate(doss):
        if isinstance(d, str):
            dtmp = PhononDOS.from_file(d, read_INS=read_INS, atom_prj=atom_prj,
                                       element_prj=element_prj)
        elif isinstance(d, PhononDOS):
            dtmp = copy.deepcopy(d)
        else:
            raise ValueError('Unknown input type for doss.')

        if unit != dtmp.unit:
            dtmp._set_unit(unit)
        dossplt.append(dtmp)

    # prj
    prj = np.array(prj, ndmin=1)
    if prj.ndim > 1: raise ValueError('Projections must be 1D.')
    if len(prj) == 0:
        prj = [[int(i+1) for i in range(len(j.doss))] for j in dossplt]
    else:
        prj = [prj for j in dossplt]

    # frequency = 0
    if plot_freq0 != True:
        freq0 = None
    else:
        freq0 = 0.

    # Gaussian broadening
    ndoss = len(dossplt)
    if gauss != 0:
        if gauss < 0:
            warnings.warn("Negative gaussian broadening width! Using its absolute value.", stacklevel=2)
            gauss = -gauss
        for i in range(ndoss):
            new_nfreq = int(len(dossplt[i].frequency) * 3)
            new_freq = np.linspace(dossplt[i].frequency[0],
                                   dossplt[i].frequency[-1], new_nfreq)
            new_dos = np.zeros([len(dossplt[i].doss), new_nfreq, 1])
            for ifreq, freq in enumerate(new_freq):
                for id, d in enumerate(dossplt[i].doss):
                    new_dos[id, ifreq, 0] = np.sum(
                        d[:, 0] * np.exp(-(freq - dossplt[i].frequency)**2 / gauss**2 / 2)
                    )
            dossplt[i].frequency = copy.deepcopy(new_freq)
            dossplt[i].doss = copy.deepcopy(new_dos)

    # plot
    if ndoss == 1 and overlap == False: # Same system, projecton into different panels
        nprj = len(prj[0])
        fig, ax = plt.subplots(nprj, 1, figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        if nprj == 1:
            ax = plot_doss(
                ax, dossplt[0].doss, dossplt[0].frequency, 'up',  [prj[0][0]],
                frequency_range, dos_range, dos_label, dos_color, dos_linestyle,
                dos_linewidth, freq0, freq0_color, freq0_linestyle, freq0_linewidth,
                legend, False, **kwargs
            )
        else:
            # new defaults: all lines in the same color.
            if np.all(dos_color==None):
                dos_color = [['tab:blue', 'tab:blue'] for i in range(nprj)]
            # Dimeonsion issue: dos plot styles must be consistent with length of input dosss
            dossref = [dossplt[0].doss[i-1] for i in prj[0]]
            commands = _plot_label_preprocess(
                dossref, dos_label, dos_color, dos_linestyle, dos_linewidth
            )
            for i in range(4):
                if np.all(commands[i]==None):
                    commands[i] = [None for j in range(nprj)]
                else:
                    commands[i] = [[commands[i][j]] for j in range(nprj)]
            for i in range(nprj):
                fig.axes[i] = plot_doss(
                    fig.axes[i], dossplt[0].doss, dossplt[0].frequency, 'up',
                    [prj[0][i]], frequency_range, dos_range, commands[0][i],
                    commands[1][i], commands[2][i], commands[3][i], freq0,
                    freq0_color, freq0_linestyle, freq0_linewidth, legend,
                    False, **kwargs
                )
    else: # Projecton of the same system into the same panel
        fig, ax = plt.subplots(ndoss, 1, figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        for i in range(ndoss):
            fig.axes[i] = plot_doss(
                fig.axes[i], dossplt[i].doss, dossplt[i].frequency, 'up', prj[i],
                frequency_range, dos_range, dos_label, dos_color, dos_linestyle,
                dos_linewidth, freq0, freq0_color, freq0_linestyle, freq0_linewidth,
                legend, False, **kwargs
            )

    # set titles and axes
    if is_thz == True:
        fig.supylabel('DOS (states/THz)')
        fig.supxlabel('Frequency (THz)')
    else:
        fig.supylabel('DOS (states/cm$^{-1}$)')
        fig.supxlabel('Frequency (cm$^{-1}$)')

    if np.all(title!=None):
        fig.suptitle(title, fontsize=fontsize)
    return fig


def plot_phonon_banddos(
    *data, unit='cm-1', q_overlap_tol=1e-4, read_INS=False, atom_prj=[],
    element_prj=[], k_label=[], dos_overlap=True, dos_prj=[], dos_gauss=0.0,
    frequency_range=[], k_range=[], dos_range=[], band_width=2, band_label=None,
    band_color=None, band_linestyle=None, band_linewidth=None, dos_label=None,
    dos_color=None, dos_linestyle=None, dos_linewidth=None, plot_freq0=True,
    freq0_color='tab:green', freq0_linestyle='-', freq0_linewidth=1.0, title=None,
    figsize=[6.4, 4.8], legend='lower right', fontsize=14, **kwargs):
    """
    Plot phonon band structure + dos for a **single** system, i.e., the
    ``bands`` and ``doss`` variables are not extendable.

    Input arguments not in the list are consistent with ``plot_phonon_doss``
    and ``plot_phonon_bands``.

    Args:
        \*data: Either 1 or 2 entries. For one enetry, it is the standard
            screen output (.out) file with both band and DOS, or ``PhononBandDOS``
            object. For 2 entries, the first entry is ``bands`` of
            ``plot_phonon_bands`` and the second is ``doss`` of ``plot_phonon_doss``.
        dos_overlap (bool): ``overlap`` of ``plot_phonon_doss``. The user can
            either plot projections into the same subplot or into separate
            subplots.
        dos_prj (list): ``prj`` of ``plot_phonon_doss``.
        dos_gauss (float): ``gauss`` of ``plot_phonon_doss``.
        band_width (int|float): Relative width of band structure, times of the
            width of a DOS subplot.
    Returns:
        fig (Figure): Matplotlib figure object

    :raise ValueError: If the unit parameter is unknown.
    """
    from CRYSTALpytools.phonons import PhononBand, PhononDOS, PhononBandDOS
    from CRYSTALpytools.base.plotbase import plot_banddos
    import warnings, copy

    # unit
    if unit.lower() == 'cm-1':
        unit = 'cm-1'
        is_thz = False
    elif unit.lower() == 'thz':
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    # instantiation and unit
    if len(data) == 1:
        if isinstance(data[0], str):
            bands = PhononBand.from_file(data[0], q_overlap_tol=q_overlap_tol)
            doss = PhononDOS.from_file(data[0], read_INS=read_INS,
                                       atom_prj=atom_prj, element_prj=element_prj)
        elif isinstance(data[0], PhononBandDOS):
            bands = copy.deepcopy(data[0].band)
            doss = copy.deepcopy(data[0].dos)
        else:
            raise ValueError('Unknown input data type for the 1st entry.')
    elif len(data) == 2:
        if isinstance(data[0], str):
            bands = PhononBand.from_file(data[0], q_overlap_tol=q_overlap_tol)
        elif isinstance(data[0], PhononBand):
            bands = copy.deepcopy(data[0])
        else:
            raise ValueError('Unknown input data type for the 1st entry.')

        if isinstance(data[1], str):
            doss = PhononDOS.from_file(data[1], read_INS=read_INS,
                                       atom_prj=atom_prj, element_prj=element_prj)
        elif isinstance(data[1], PhononDOS):
            doss = copy.deepcopy(data[1])
        else:
            raise ValueError('Unknown input data type for the 2nd entry.')
    else:
        raise ValueError('Input parameter length does not meet requirements.')

    if unit != doss.unit:
        doss._set_unit(unit)
    if unit != bands.unit:
        bands._set_unit(unit)

    # set projections
    if len(dos_prj) == 0:
        dos_prj = [i+1 for i in range(len(doss.doss))]

    # frequency = 0
    if plot_freq0 != True:
        freq0 = None
    else:
        freq0 = 0.

    # Gaussian broadening
    if dos_gauss != 0:
        if dos_gauss < 0:
            warnings.warn("Negative gaussian broadening width! Using its absolute value.", stacklevel=2)
            dos_gauss = -dos_gauss
        new_nfreq = int(len(doss.frequency) * 3)
        new_freq = np.linspace(doss.frequency[0], doss.frequency[-1], new_nfreq)
        new_dos = np.zeros([len(doss.doss), new_nfreq, 1])
        for ifreq, freq in enumerate(new_freq):
            for id, d in enumerate(doss.doss):
                new_dos[id, ifreq, 0] = np.sum(
                    d[:, 0] * np.exp(-(freq - doss.frequency)**2 / dos_gauss**2 / 2)
                )
        doss.frequency = copy.deepcopy(new_freq)
        doss.doss = copy.deepcopy(new_dos)

    fig = plot_banddos(
        bands, doss, k_label, 'up', dos_overlap, dos_prj, frequency_range,
        k_range, dos_range, band_width, band_label, band_color, band_linestyle,
        band_linewidth, dos_label, dos_color, dos_linestyle, dos_linewidth,
        freq0, freq0_color, freq0_linestyle, freq0_linewidth, figsize,
        legend, **kwargs
    )
    # set titles and axes
    if is_thz == True:
        fig.supylabel('Frequency (THz)')
    else:
        fig.supylabel('Frequency (cm$^{-1}$)')

    return fig

##############################################################################
#                                                                            #
#                                    SPECTRA                                 #
#                                                                            #
##############################################################################

#--------------------------------------XRD------------------------------------#

def plot_XRD(*xrd, option='LP', shift=10, label=None, color=None,
             linestyle=None, linewidth=None, theta_range=[], title=None,
             figsize=[6.4, 4.8], legend='upper left', fontsize=14, **kwargs):
    """
    Plot the XRD spectra of multiple systems into the same plot axes.

    .. note::

        The highest intensity is normalized to 100.

    Args:
        xrd (str|XRD): File name or ``spectra.XRD`` objects. Extendable.
        option (str): *File name inputs only* 'NC' for no correction (The
            'INTENS' col); 'LD' for Lorentz and polarization effects
            ('INTENS-LP') and 'DW' for LD with Debye-Waller thermal factors
            ('INTENS-LP-DW').
        shift (float): If multiple spectra are plotted, shifting them up by the
            given value. Shift length is the value after normalization.
        label (list|str|None): List of plot labels. 'None' for the default
            values ('# \<number\>') and string for prefix the string before
            number. Otherwise should be a 1\*nXRD list.
        color (list|str|None): If str, use the same color for all the plot
            lines. If 1\*nXRD, use the color defined for every plot line.
            'None' for default values (matplotlib Tableau Palette).
        linestyle (list|str|None): See explanations of ``color``.
        linewidth (list|float|None): See explanations of ``color``.
        theta_range (list): 1\*2 list of theta range in degree.
        title (str|None): The title of the plot. 'None' for no title.
        figsize (list): Matplotlib figure size.
        legend (str|None): Location of legend. None for not adding legend.
        fontsize (int): Fontsize of the axis label and title.
        \*\*kwargs: Other parameters passed to matplotlib ``Axes.plot()`` method.

    Returns:
        fig (Figure): Matplotlib figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from CRYSTALpytools.spectra import XRD
    from CRYSTALpytools.base.plotbase import _plot_label_preprocess
    import warnings, copy

    objs = []
    for o in xrd:
        if isinstance(o, str):
            objs.append(XRD.from_file(o, option=option))
        elif isinstance(o, XRD):
            objs.append(o)
        else:
            raise TypeError('Input must be either file name or spectra.XRD class.')
    nplot = len(objs)

    # normalization
    maxintens = [np.max(o.spectra) for o in objs]
    maxintens = np.max(maxintens)

    # labels and other plot settings.
    if np.all(label==None):
        label = ['# {:d}'.format(i+1) for i in range(nplot)]
    else:
        if isinstance(label, str):
            label = ['{} {:d}'.format(label, i+1) for i in range(nplot)]
        else:
            if len(label) != nplot:
                warnings.warn("Inconsistent lengths of number of plots and plot labels. Using default labels.",
                              stacklevel=2)
                label = ['# {:d}'.format(i+1) for i in range(nplot)]

    commands = _plot_label_preprocess(
        np.zeros([nplot, 1, 1]), label, color, linestyle, linewidth) # dummy doss

    # plot
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    for ispec, o in enumerate(objs):
        spectra = copy.deepcopy(o.spectra)
        spectra = spectra / maxintens * 100 + shift * ispec
        ax.plot(o.theta, spectra, label=commands[0][ispec][0],
                color=commands[1][ispec][0], linestyle=commands[2][ispec][0],
                linewidth=commands[3][ispec][0], **kwargs)

    # plot setups
    if len(theta_range) != 0:
        ax.set_xlim([np.min(theta_range), np.max(theta_range)])
    ax.set_xlabel(r'2$\theta^{\circ}$', fontsize=fontsize)
    ax.set_ylabel(r'Intensity (arb. u.)', fontsize=fontsize)
    ax.set_ylim([0, 100 + shift*(nplot-1)])
    _ = ax.get_yaxis().set_ticks([])
    if np.all(title!=None):
        ax.set_title(title, fontsize=fontsize)
    if np.all(legend!=None):
        ax.legend(loc=legend)
    return fig

#--------------------------------------IR-------------------------------------#

def plot_IR(*ir, unit='cm-1', option='LG', shift=0, label=None, color=None,
            linestyle=None, linewidth=None, x_range=[], title=None, figsize=[6.4, 4.8],
            legend='upper left', sharey=True, fontsize=14, **kwargs):
    """
    Plot the IR spectra of multiple systems into the same plot axes. For
    reflectance spectra, nDirection\*1 plots are generated. All the input files
    must have the same symmetry.

    .. note::

        The highest intensity is normalized to 100.

    Args:
        ir (str|IR): IRSPEC.DAT file name or ``spectra.IR`` objects. Extendable.
        option (str): Broadening method. 'LG' for Lorentzian-Gaussian, 'V' for
            Voigt, 'RS' for Rayleigh spherical particles, 'RE' for Rayleigh with
            elipsoid particles, 'REFL' for reflectance spectra with 'LG'.
            *Periodic systems only*.
        shift (float): If multiple spectra are plotted, shifting them up by the
            given value. Shift length is the value after normalization.
        label (list|str|None): List of plot labels. 'None' for the default
            values ('# \<number\>') and string for prefix the string before
            number. Otherwise should be a 1\*nIR list.
        color (list|str|None): If str, use the same color for all the plot
            lines. If 1\*nIR, use the color defined for every plot line. 'None'
            for default values (matplotlib Tableau Palette).
        linestyle (list|str|None): See explanations of ``color``.
        linewidth (list|float|None): See explanations of ``color``.
        x_range (list): 1\*2 list of x axis range.
        title (str|None): The title of the plot. 'None' for no title.
        figsize (list): Matplotlib figure size.
        legend (str|None): Location of legend. None for not adding legend.
        sharey (bool): Whether to share the y-axis among subplots. Share x is
            enforced.
        fontsize (int): Fontsize of the axis label and title.
        \*\*kwargs: Other parameters passed to matplotlib ``Axes.plot()`` method.

    Returns:
        fig (Figure): Matplotlib figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from CRYSTALpytools.spectra import IR
    from CRYSTALpytools.base.plotbase import _plot_label_preprocess
    import warnings, copy

    objs = []
    for o in ir:
        if isinstance(o, str):
            objs.append(IR.from_file(o))
        elif isinstance(o, IR):
            objs.append(copy.deepcopy(o))
        else:
            raise TypeError('Input must be either file name or spectra.IR class.')

    nsystem = len(objs)
    for o in objs:
        if o.type == 'molecule': option='LG'

    valid_option = ['LG', 'V', 'RS', 'RE', 'REFL']
    if option.upper() not in valid_option:
        raise ValueError("Unknown option: '{}'.".format(option))
    option = option.upper()

    # single system
    if nsystem == 1:
        return objs[0].plot(unit, option, True, False, shift, label, color,
                            linestyle, linewidth, x_range, title, figsize,
                            legend, sharey, fontsize, None, **kwargs)

    # plot commands
    if np.all(label==None):
        label = ['# {:d}'.format(i+1) for i in range(nsystem)]
    else:
        if isinstance(label, str):
            label = ['{} {:d}'.format(label, i+1) for i in range(nsystem)]
        else:
            if len(label) != nsystem:
                warnings.warn(
                    "Inconsistent lengths of number of materials and plot labels. Using default labels.",
                    stacklevel=2)
                label = ['# {:d}'.format(i+1) for i in range(nsystem)]

    doss = np.zeros([nsystem, 1, 1]) # pseudo doss
    commands = _plot_label_preprocess(doss, label, color, linestyle, linewidth)
    # repeat labels to walk around default label settings of class plot
    # now commends[0] is a nSystem\*nDir list. Should be called via commands[0][isystem]
    if option != 'REFL':
        commands[0] = [[i[0]] for i in commands[0]]
    else:
        ndir = len(objs[0].reflectance)
        for o in objs:
            if len(o.reflectance) != ndir: raise Exception('Inconsistent symmetry of the entries.')
        commands[0] = [[i[0] for j in range(ndir)] for i in commands[0]]

    # normalize data
    vmax = []
    if option != 'REFL':
        if option == 'LG':
            for io in range(nsystem): vmax.append(np.max(objs[io].absorbance[0]))
        elif option == 'V':
            for io in range(nsystem): vmax.append(np.max(objs[io].absorbance[1]))
        elif option == 'RS':
            for io in range(nsystem): vmax.append(np.max(objs[io].absorbance[2]))
        elif option == 'RE':
            for io in range(nsystem): vmax.append(np.max(objs[io].absorbance[3]))

        vmax = np.max(vmax)
        for io in range(nsystem):
            objs[io].absorbance = objs[io].absorbance / vmax * 100 + shift*io
    else:
        for io in range(nsystem): vmax.append(np.max(objs[io].reflectance))
        vmax = np.max(vmax)
        for io in range(nsystem):
            objs[io].reflectance = objs[io].reflectance / vmax * 100 + shift*io

    # plot first entry
    fig = objs[0].plot(unit, option, False, False, 0, commands[0][0],
                       commands[1][0][0], commands[2][0][0], commands[3][0][0],
                       x_range, title, figsize, legend, sharey, fontsize, None,
                       **kwargs)
    for io in range(1, nsystem):
        fig = objs[io].plot(unit, option, False, False, 0, commands[0][io],
                            commands[1][io][0], commands[2][io][0], commands[3][io][0],
                            x_range, title, figsize, legend, sharey, fontsize, fig,
                            **kwargs)
    # modify legend and annotations
    if np.all(legend!=None) and np.all(label!=None):
        fig.axes[0].legend(loc=legend)

    for iax in range(len(fig.axes)):
        if iax > 0: fig.axes[iax].get_legend().remove()
        if option == 'REFL':
            xpos = fig.axes[iax].get_xlim()[1]
            ypos = fig.axes[iax].get_ylim()[1]
            fig.axes[iax].text(xpos,ypos, 'direction {:d}'.format(iax+1), fontsize=fontsize,
                               horizontalalignment='right', verticalalignment='top')
    return fig

#-------------------------------------Raman-----------------------------------#

def plot_Raman(*raman, option='poly', overlap=True,
               direction=['xx', 'xy', 'xz', 'yy', 'yz', 'zz'], shift=0,
               label=None, color=None, linestyle=None, linewidth=None, x_range=[],
               title=None, figsize=[6.4, 4.8],legend='upper left', sharey=True,
               fontsize=14, **kwargs):
    """
    Plot the Raman spectra of multiple systems into the same plot axes.

    Available options:

    * 'tot': Plot total raman spectra only. Plots with a single panel is generated.  
    * 'poly': Plot total, parallel and perpendicular spectra into 3 subplots.  
    * 'single': Plot single crystal spectra of specified directions into subplots.

    .. note::

        The highest intensity is normalized to 100.

    Args:
        raman (str|Raman): RAMSPEC.DAT file name or ``spectra.Raman`` objects.
            Extendable.
        option (str): 'tot', 'poly' or 'single', see above.
        overlap (bool): If more than 1 inequivalent directions exists, whether
            to plot spectra into the same plot or into subplots.
        direction (list|str): *``option='single'`` only* Specify the directions
            of single crystal spectra to plot.
        shift (float): If multiple spectra are plotted, shifting them up by the
            given value. Shift length is the value after normalization.
        label (list|str|None): List of plot labels. 'None' for the default
            values ('# \<number\>') and string for prefix the string before
            number. Otherwise should be a 1\*nRaman list.
        color (list|str|None): If str, use the same color for all the plot
            lines. If 1\*nRaman, use the color defined for every plot line.
            'None' for default values (matplotlib Tableau Palette).
        linestyle (list|str|None): See explanations of ``color``.
        linewidth (list|float|None): See explanations of ``color``.
        x_range (list): 1\*2 list of x axis range.
        title (str|None): The title of the plot. 'None' for no title.
        figsize (list): Matplotlib figure size.
        legend (str|None): Location of legend. None for not adding legend.
        sharey (bool): Whether to share the y-axis among subplots. Share x is
            enforced.
        fontsize (int): Fontsize of the axis label and title.
        \*\*kwargs: Other parameters passed to matplotlib ``Axes.plot()`` method.

    Returns:
        fig (Figure): Matplotlib figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from CRYSTALpytools.spectra import Raman
    from CRYSTALpytools.base.plotbase import _plot_label_preprocess
    import warnings, copy

    objs = []
    for o in raman:
        if isinstance(o, str):
            objs.append(Raman.from_file(o))
        elif isinstance(o, Raman):
            objs.append(copy.deepcopy(o))
        else:
            raise TypeError('Input must be either file name or spectra.Raman class.')

    nsystem = len(objs)

    valid_option = ['tot', 'single', 'poly']
    if option.lower() not in valid_option:
        raise ValueError("Unknown option: '{}'.".format(option))
    option = option.lower()

    # single system
    if nsystem == 1:
        return objs[0].plot(option, True, False, direction, shift, label, color,
                            linestyle, linewidth, x_range, title, figsize, legend,
                            sharey, fontsize, None, **kwargs)

    # plot commands
    if np.all(label==None):
        label = ['# {:d}'.format(i+1) for i in range(nsystem)]
    else:
        if isinstance(label, str):
            label = ['{} {:d}'.format(label, i+1) for i in range(nsystem)]
        else:
            if len(label) != nsystem:
                warnings.warn(
                    "Inconsistent lengths of number of materials and plot labels. Using default labels.",
                    stacklevel=2)
                label = ['# {:d}'.format(i+1) for i in range(nsystem)]

    doss = np.zeros([nsystem, 1, 1]) # pseudo doss
    commands = _plot_label_preprocess(doss, label, color, linestyle, linewidth)
    # repeat labels to walk around default label settings of class plot
    # now commends[0] is a nSystem\*nDir list. Should be called via commands[0][isystem]
    if option == 'poly':
        commands[0] = [[i[0] for j in range(3)] for i in commands[0]]
    elif option == 'single':
        direction = np.array(direction, ndmin=1).tolist()
        commands[0] = [[i[0] for j in range(len(direction))] for i in commands[0]]
    elif option == 'tot':
        commands[0] = [[i[0]] for i in commands[0]]

    # normalize data
    vmax = []
    dirs = {'xx' : 0, 'xy' : 1, 'xz' : 2, 'yx' : 1, 'yy' : 3, 'yz' : 4, 'zx' : 2, 'zy' : 4, 'zz' : 5,}
    if option == 'tot':
        for io in range(nsystem): vmax.append(np.max(objs[io].poly[0]))
        vmax = np.max(vmax)
        for io in range(nsystem):
            objs[io].poly = objs[io].poly / vmax * 100 + shift*io
    elif option == 'poly':
        for io in range(nsystem): vmax.append(np.max(objs[io].poly))
        vmax = np.max(vmax)
        for io in range(nsystem):
            objs[io].poly = objs[io].poly / vmax * 100 + shift*io
    elif option == 'single':
        try: idxdir = [dirs[i.lower()] for i in direction]
        except KeyError:
            raise ValueError("Unknown direction name: '{}'.".format(i))
        for io in range(nsystem): vmax.append(np.max(objs[io].single[idxdir]))
        vmax = np.max(vmax)
        for io in range(nsystem):
            objs[io].single = objs[io].single / vmax * 100 + shift*io

    # plot first entry
    fig = objs[0].plot(option, False, False, direction, 0, commands[0][0],
                       commands[1][0][0], commands[2][0][0], commands[3][0][0],
                       x_range, title, figsize, legend, sharey, fontsize, None,
                       **kwargs)
    # others
    for io in range(1, nsystem):
        fig = objs[io].plot(option, False, False, direction, 0, commands[0][io],
                            commands[1][io][0], commands[2][io][0], commands[3][io][0],
                            x_range, title, figsize, legend, sharey, fontsize, fig,
                            **kwargs)
    # modify legend and annotations
    if np.all(legend!=None) and np.all(label!=None):
        fig.axes[0].legend(loc=legend)

    poly_text = ['Total', 'Parallel', 'Perpendicular']
    single_text = direction
    for iax in range(len(fig.axes)):
        if iax > 0: fig.axes[iax].get_legend().remove()
        xpos = fig.axes[iax].get_xlim()[1]
        ypos = fig.axes[iax].get_ylim()[1]
        if option == 'poly':
            fig.axes[iax].text(xpos,ypos, poly_text[iax], fontsize=fontsize,
                               horizontalalignment='right', verticalalignment='top')
        elif option == 'single':
            fig.axes[iax].text(xpos,ypos, single_text[iax], fontsize=fontsize,
                               horizontalalignment='right', verticalalignment='top')
    return fig



#-----------------------------------ANHARMONIC--------------------------------#

def plot_cry_spec(transitions, typeS, components=False, bwidth=5, stdev=3, eta=0.5,
                  fmin=None, fmax=None, ylim=None, savefig=False, dpi=300,
                  filetype='png', exp_spec=None, sep=";", show=True,
                  export_csv=False, label=None, xlabel='Wavenumber [cm$^{-1}$]',
                  ylabel='Intensity [arb. u.]', linewidth=2.0, padd=100,
                  fontsize=12, style=None, compstyle=None, nopadding=False,
                  figsize=(16, 6)):
    """
    .. note::

        **This is not for the released feature of CRYSTAL23 v1.0.1**

    This function enables the simulation of vibrational spectra based on a 2D 
    NumPy array containing a list of transition frequencies and the 
    corresponding intensities. The code allows users to model spectral
    broadening according to various profiles (Gaussian, Lorentzian, 
    pseudo-Voigt), or zero broadening (Dirac deltas-like lines). Please, note
    that by turning the optional argument 'component' to `True` you can
    additionally plot contributions arising from each transition.

    Args:
        transitions (float|numpy.ndarray): Array containing transition frequencies
        (axis=0) and corresponding intensities (axis=1).
        typeS (str): String specifying the spectral profile: 'bars',
        'lorentz', 'gauss', 'pvoigt'. 
        components (bool, optional): Whether to plot contributions arising from
        each transition (default is `False`).  
        bwidth (float, optional): Half-width at half-maximum of the Lorentzian 
        profile (default is 5).
        stdev (float, optional): Standard deviation of the Gaussian profile 
        (default is 5).
        eta (float, optional): Fraction of Lorentzian character in pseudo-Voigt
        profile (default is 0.5).
        fmin (float, optional): Minimum frequency.
        fmax(float, optional): Maximum frequency.
        ylim (float, optional): Maximum intensity.
        savefig (bool, optional): Whether to save the figure (default is `False`).
        dpi (float, optional): Dots per inches (default is 300).
        filetype (str, optional): File extension (default is 'png').
        show (bool, optional): Whether to show the figure (default is `True`).
        export_csv (bool, optional): Whether to save plot in csv format (default is 
        `False`).
        xlabel (str, optional): x-axis label (default is 'Wavenumber [cm$^{-1}$]').
        ylabel (str, optional): y-axis label (default is 'Intensity [arb. u.]').
        linewidth (float): Linewidth (default is 2.0).
        padd (float, optional): left- and right- hand side padding expressed in the
        same unit of the quantity reported in x-axis (default is 100).
        fontsize (integer, optional): Fontsize (default is 12).
        style (str, optional): String specifying Matplotlib style. 
        compstyle (str|list, optional): List containing Matplotlib styles to plot
        each component. 
        nopadding (bool, optional): Whether to remove padding (default is `False`).
        figsize (real|list, optional): List of two numbers specifying the aspect
        ratio of the figure (default is [16, 6]).

    Returns:
        None
    """

    import math, warnings
    import time
    from copy import deepcopy

    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import genfromtxt

    warnings.warn('This is not a released feature of CRYSTAL23 v1.0.1, make sure that you know what you are doing.',
                  stacklevel=2)

    if (show):
        plt.figure(figsize=figsize)
    if (ylim is not None):
        plt.ylim(0, ylim)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)

    bars = False
    lorentz = False
    gauss = False
    pseudo_voigt = False

    if typeS == 'bars':
        bars = True

    if typeS == 'lorentz':
        lorentz = True

    if typeS == 'gauss':
        gauss = True

    if typeS == 'pvoigt':
        pseudo_voigt = True

    n = 20000

    if fmin is None:
        fmin = min(transitions[:, 0] - padd)
    if fmax is None:
        fmax = max(transitions[:, 0] + padd)

    x = np.linspace(fmin, fmax, num=n)
    y = np.zeros(n)

    spec_data = np.block([[x], [y]]).T
    sbuff = np.block([[x], [y]]).T

    if bars:
        spec_data = np.concatenate((spec_data, transitions), axis=0)
        spec_data = spec_data[spec_data[:, 0].argsort()]
    elif lorentz:
        iL = 0
        L = []
        for i in range(len(transitions)):
            if transitions[i, 1] == 0:
                continue
            for j, f in enumerate(spec_data[:, 0]):
                lorentz = (1/math.pi)*bwidth / \
                    ((f-transitions[i, 0])**2+bwidth**2)*transitions[i, 1]
                sbuff[j, 1] = lorentz
            L.append(deepcopy(sbuff))
            iL = iL + 1
        if (not components):
            for i in range(len(L)):
                spec_data[:, 1] = spec_data[:, 1] + L[i][:, 1]
        else:
            for i in range(len(L)):
                plt.plot(spec_data[:, 0], L[i][:, 1], linewidth=linewidth)
            for i in range(len(L)):
                spec_data[:, 1] = spec_data[:, 1] + L[i][:, 1]

    elif gauss:
        G = []
        for i in range(len(transitions)):
            if transitions[i, 1] == 0:
                continue
            for j, f in enumerate(spec_data[:, 0]):
                gauss = (1/(stdev*math.sqrt(2*math.pi))) * \
                    math.exp(-((f-transitions[i, 0])**2) /
                             (2*stdev**2))*transitions[i, 1]
                sbuff[j, 1] = gauss
            G.append(deepcopy(sbuff))
        if (not components):
            for i in range(len(G)):
                spec_data[:, 1] = spec_data[:, 1] + G[i][:, 1]
        else:
            for i in range(len(G)):
                plt.plot(spec_data[:, 0], G[i][:, 1], linewidth=linewidth)
            for i in range(len(G)):
                spec_data[:, 1] = spec_data[:, 1] + G[i][:, 1]

    elif pseudo_voigt:
        V = []
        for i in range(len(transitions)):
            if transitions[i, 1] == 0:
                continue
            for j, f in enumerate(spec_data[:, 0]):
                gauss = (1/(stdev*math.sqrt(2*math.pi))) * \
                    math.exp(-((f-transitions[i, 0])**2) /
                             (2*stdev**2))*transitions[i, 1]
                lorentz = (1/math.pi)*bwidth / \
                    ((f-transitions[i, 0])**2+bwidth**2)*transitions[i, 1]
                sbuff[j, 1] = eta*lorentz + (1-eta)*gauss
            V.append(deepcopy(sbuff))
        if (not components):
            for i in range(len(V)):
                spec_data[:, 1] = spec_data[:, 1] + V[i][:, 1]
        else:
            for i in range(len(V)):
                if (compstyle is not None):
                    plt.plot(spec_data[:, 0], V[i][:, 1], compstyle[i],
                             linewidth=linewidth)
                else:
                    plt.plot(spec_data[:, 0], V[i][:, 1], linewidth=linewidth)
            for i in range(len(V)):
                spec_data[:, 1] = spec_data[:, 1] + V[i][:, 1]

    if (exp_spec is not None):
        exp_data = genfromtxt(exp_spec, delimiter=sep)
        area_spec_data = np.trapz(spec_data[:, 1], spec_data[:, 0])
        area_exp_data = np.trapz(exp_data[:, 1], exp_data[:, 0])
        norm_fac = area_spec_data / area_exp_data
        baseline = 0.2
        exp_data[:, 1] = exp_data[:, 1] * norm_fac - baseline  # * 0.5
        plt.plot(exp_data[:, 0], exp_data[:, 1], 'r-', linewidth=linewidth)

    if label is not None:
        plt.plot(spec_data[:, 0], spec_data[:, 1], linewidth=linewidth,
                 label=label)
    elif (style is not None):
        plt.plot(spec_data[:, 0], spec_data[:, 1], style, linewidth=linewidth)
    else:
        plt.plot(spec_data[:, 0], spec_data[:, 1], linewidth=linewidth)

    if (savefig):
        plt.savefig(typeS + time.strftime("%Y-%m-%d_%H%M%S.") + filetype,
                    format=filetype, dpi=dpi)
    if (show):
        plt.show()

    if (export_csv):
        np.savetxt(typeS + time.strftime("%Y-%m-%d_%H%M%S.") + 'csv',
                   spec_data, delimiter=';')


def plot_cry_spec_multi(files, typeS, components=False, bwidth=5, stdev=3,
                        eta=0.5, fmin=None, fmax=None, ylim=None,
                        savefig=False, dpi=300, filetype='png', label=None,
                        xlabel='Wavenumber [cm$^{-1}$]',
                        ylabel='Instensity [arb. u.]', linewidth=2.0, padd=100,
                        fontsize=12, style=None, nopadding=False,
                        figsize=(16, 6)):
    """
    .. note::

        **This is not for the released feature of CRYSTAL23 v1.0.1**

    This function is a wrapper for `plot_spec` function, enablng the simulation 
    of many vibrational spectra coming from a list of NumPy array.  

    Args:
        transitions (float|numpy.ndarray): Array containing transition frequencies
        (axis=0) and corresponding intensities (axis=1).
        typeS (str): String specifying the spectral profile: 'bars',
        'lorentz', 'gauss', 'pvoigt'. 
        components (bool, optional): Whether to plot contributions arising from
        each transition (default is `False`).  
        bwidth (float, optional): Half-width at half-maximum of the Lorentzian 
        profile (default is 5).
        stdev (float, optional): Standard deviation of the Gaussian profile 
        (default is 5).
        eta (float, optional): Fraction of Lorentzian character in pseudo-Voigt
        profile (default is 0.5).
        fmin (float, optional): Minimum frequency.
        fmax(float, optional): Maximum frequency.
        ylim (float, optional): Maximum intensity.
        savefig (bool, optional): Whether to save the figure (default is `False`).
        dpi (float, optional): Dots per inches (default is 300).
        filetype (str, optional): File extension (default is 'png').
        xlabel (str, optional): x-axis label (default is 'Wavenumber [cm$^{-1}$]').
        ylabel (str, optional): y-axis label (default is 'Intensity [arb. u.]').
        linewidth (float): Linewidth (default is 2.0).
        padd (float, optional): left- and right- hand side padding expressed in the
        same unit of the quantity reported in x-axis (default is 100).
        fontsize (integer, optional): Fontsize (default is 12).
        style (str, optional): String specifying Matplotlib style. 
        nopadding (bool, optional): Whether to remove padding (default is `False`).
        figsize (real|list, optional): List of two numbers specifying the aspect
        ratio of the figure (default is [16, 6]).

    Returns:
        None
    """

    import time, warnings

    import matplotlib.pyplot as plt

    warnings.warn('This is not a released feature of CRYSTAL23 v1.0.1, make sure that you know what you are doing.',
                  stacklevel=2)

    plt.figure(figsize=figsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)

    for i, transitions in enumerate(files):
        if (label is not None):
            plot_spec(transitions, typeS, components, bwidth, stdev, eta, fmin,
                      fmax, ylim, show=False, savefig=False, label=label[i],
                      linewidth=linewidth, padd=padd, nopadding=nopadding,
                      fontsize=fontsize, xlabel=xlabel, ylabel=ylabel)
        elif (style is not None):
            plot_spec(transitions, typeS, components, bwidth, stdev, eta, fmin,
                      fmax, ylim, show=False, savefig=False,
                      linewidth=linewidth, padd=padd, nopadding=nopadding,
                      fontsize=fontsize, style=style[i], xlabel=xlabel,
                      ylabel=ylabel)
        else:
            plot_spec(transitions, typeS, components, bwidth, stdev, eta, fmin,
                      fmax, ylim, show=False, savefig=False,
                      linewidth=linewidth, padd=padd, nopadding=nopadding,
                      fontsize=fontsize, xlabel=xlabel, ylabel=ylabel)

    if (label is not None):
        plt.legend(loc='upper left', fontsize=fontsize)

    if (savefig):
        plt.savefig("multi_" + typeS + time.strftime("%Y-%m-%d_%H%M%S.") +
                    filetype, format=filetype, dpi=dpi)

    plt.show()

#------------------------------------------------------------------------------#
#--------------------------------obsolete functions----------------------------#
#------------------------------------------------------------------------------#

def plot_cry_xrd(xrd_obj):
    """
    Deprecated. Use ``plot_XRD``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_XRD' instead.",
                  stacklevel=2)
    return plot_XRD(xrd_obj)

def plot_cry_irspec(irspec, x_unit='cm-1', y_mode='LG', figsize=None, linestyle='-',
                    linewidth=1.5, color='tab:blue', freq_range=None, int_range=None,
                    label=None):
    """
    Deprecated. Use ``plot_IR``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_IR' instead.",
                  stacklevel=2)
    if np.all(figsize==None): figsize = [6.4, 4.8]
    if np.all(freq_range==None): freq_range = []

    fig = plot_IR(irspec, unit=x_unit, option=y_mode, label=label, linewidth=linewidth,
                  color=color, linestyle=linestyle, x_range=freq_range, figsize=figsize)
    ax = fig.axes[0]
    return fig, ax

def plot_cry_ramspec(ramspec,  y_mode='total', figsize=None, linestyle='-',
                     linewidth=1.5, color='tab:blue', freq_range=None, int_range=None,
                     label=None):
    """
    Deprecated. Use ``plot_Raman``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_Raman' instead.",
                  stacklevel=2)
    if np.all(figsize==None): figsize = [6.4, 4.8]
    if np.all(freq_range==None): freq_range = []

    if y_mode.lower() == 'total':
        y_mode = 'tot'
    elif y_mode.lower() in ['parallel', 'perpendicular']:
        y_mode = 'poly'
    elif y_mode.lower() in ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']:
        dirs = y_mode.lower()
        y_mode = 'single'

    fig = plot_Raman(ramspec, option=y_mode, direction=dirs, label=label, linewidth=linewidth,
                     color=color, linestyle=linestyle, x_range=freq_range, figsize=figsize)
    ax = fig.axes[0]
    return fig, ax

def plot_cry_ela(choose, ndeg, *args, dpi=200, filetype=".png", transparency=False):
    """
    Deprecated. Use ``plot_elastics3D``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_elastics3D' instead.",
                  stacklevel=2)

    for nfig, i in enumerate(args):
        fig = plot_elastics3D(i, property=choose, uniform_scale=True, nphi=ndeg,
                              ntheta=ndeg, nchi=ndeg)
        fig.savefig(choose + '{:d}'.format(nfig) + filetype,
                    dpi=dpi, transparent=transparency)
    return

def plot_vecfield2D_m(header, dens, quivscale, name='MAG', levels=150, dpi=400):
    """
    Deprecated and incompatible with new functions. Give error.
    """
    raise DeprecationWarning("The called method is deprecated and incompatible with the new one. Call 'plot_relativistics2D()' instead.")


def plot_vecfield2D_j(header, dens, quivscale, name='SC', levels=150, dpi=400):
    """
    Deprecated and incompatible with new functions. Give error.
    """
    raise DeprecationWarning("The called method is deprecated and incompatible with the new one. Call 'plot_relativistics2D()' instead.")


def plot_vecfield2D_J(header, dens_JX, dens_JY, dens_JZ, quivscale, name='SCD', levels=150, dpi=400):
    """
    Deprecated and incompatible with new functions. Give error.
    """
    raise DeprecationWarning("The called method is deprecated and incompatible with the new one. Call 'plot_relativistics2D()' instead.")


def plot_electron_band(bands, unit='eV', k_labels=None, mode='single',
                       not_scaled=False, energy_range=None, k_range=None,
                       color='blue', labels=None, linestl='-', linewidth=1,
                       fermi='forestgreen', fermiwidth=1.5, fermialpha=1, 
                       title=None, figsize=None, scheme=None, sharex=True,
                       sharey=True, fontsize=12):
    """
    Deprecated. Use ``plot_electron_bands``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_bands' instead.",
                  stacklevel=2)
    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    if np.all(energy_range==None):
        energy_range=[]
    if np.all(k_range==None):
        k_range=[]
    if np.all(figsize==None):
        figsize=[6.4, 4.8]

    fig = plot_electron_bands(
        *bands, unit=unit, k_label=k_labels, mode=mode, not_scaled=not_scaled,
        energy_range=energy_range, k_range=k_range, band_label=labels, band_color=color,
        band_linestyle=linestl, band_linewidth=linewidth, fermi_color=fermi,
        fermi_linewidth=fermiwidth, title=title, figsize=figsize, layout=scheme,
        sharex=sharex, sharey=sharey, fontsize=fontsize
    )
    return fig, fig.axes


def plot_cry_band(bands, k_labels=[], energy_range=[], title=None, not_scaled=True,
                  mode='single', linestl='-', linewidth=1, color='blue',
                  fermi='forestgreen', k_range=[], labels=None, figsize=[6.4, 4.8],
                  scheme=None, sharex=True, sharey=True, fermiwidth=1.5, fermialpha=1,
                  fontsize=12):
    """
    Deprecated. Use ``plot_electron_bands``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_bands' instead.",
                  stacklevel=2)

    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    fig = plot_electron_bands(
        *bands, k_label=k_labels, mode=mode, not_scaled=not_scaled,
        energy_range=energy_range, k_range=k_range, band_label=labels, band_color=color,
        band_linestyle=linestl, band_linewidth=linewidth, fermi_color=fermi,
        fermi_linewidth=fermiwidth, title=title, figsize=figsize, layout=scheme,
        sharex=sharex, sharey=sharey, fontsize=fontsize
    )
    return fig, fig.axes


def plot_electron_dos(doss, unit='eV', beta='up', overlap=False, prj=None,
                      energy_range=None, dos_range=None, color='blue',
                      labels=None, linestl=None, linewidth=1, fermi='forestgreen',
                      title=None, figsize=None):
    """
    Deprecated. Use ``plot_electron_doss``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_doss' instead.",
                  stacklevel=2)
    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    if np.all(energy_range==None):
        energy_range=[]
    if np.all(dos_range==None):
        dos_range=[]
    if np.all(prj==None):
        prj=[]
    if np.all(figsize==None):
        figsize=[6.4, 4.8]

    fig = plot_electron_doss(
        *doss, unit=unit, beta=beta, overlap=overlap, prj=prj,
        energy_range=energy_range, dos_range=dos_range, dos_label=labels,
        dos_color=color, dos_linestyle=linestl, dos_linewidth=linewidth,
        fermi_color=fermi, title=title, figsize=figsize
    )
    return fig, fig.axes


def plot_cry_doss(doss, color='blue', fermi='forestgreen', overlap=False,
                  labels=None, figsize=[6.4, 4.8], linestl=None,
                  linewidth=1.0, title=None, beta='down', energy_range=[],
                  dos_range=[], prj=[]):
    """
    Deprecated. Use ``plot_electron_doss``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_doss' instead.",
                  stacklevel=2)
    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    fig = plot_electron_doss(
        *doss, beta=beta, overlap=overlap, prj=prj,
        energy_range=energy_range, dos_range=dos_range, dos_label=labels,
        dos_color=color, dos_linestyle=linestl, dos_linewidth=linewidth,
        fermi_color=fermi, title=title, figsize=figsize
    )
    return fig, fig.axes


def plot_cry_es(bands, doss, k_labels=[], color_bd='blue', color_doss='blue',
                fermi='forestgreen', energy_range=[], linestl_bd=None,
                linestl_doss=None, linewidth=1.0, prj=[], figsize=[6.4, 4.8],
                labels=None, dos_range=[], title=None, dos_beta='down'):
    """
    Deprecated. Use ``plot_electron_doss``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_banddos' instead.",
                  stacklevel=2)

    fig = plot_electron_banddos(
        bands, doss, k_label=k_labels, dos_beta=dos_beta, dos_prj=prj,
        energy_range=energy_range, dos_range=dos_range, band_color=color_bd,
        band_linestyle=linestl_bd, band_linewidth=linewidth, dos_label=labels,
        dos_color=color_doss, dos_linestyle=linestl_doss, dos_linewidth=linewidth,
        fermi_color=fermi, title=title, figsize=figsize
    )
    return fig, fig.axes

def plot_dens_ECHG(obj_echg, unit='Angstrom',  xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None):
    """
    Deprecated. Use ``plot_ECHG``.
    """
    import warnings
    import numpy as np

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_banddos' instead.",
                  stacklevel=2)
    if np.all(cmap_min!=None) and np.all(cmap_max!=None):
        levels = np.linspace(cmap_min, cmap_max, 150)
    else:
        levels = 150
    fig = plot_ECHG(obj_echg, unit=unit, levels=levels, option='charge',
                    x_ticks=xticks, y_ticks=yticks)
    return fig, fig.axes


def plot_spin_ECHG(obj_echg, unit='Angstrom', levels=150, xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None):
    """
    Deprecated. Use ``plot_ECHG``.
    """
    import warnings
    import numpy as np

    warnings.warn("You are calling a deprecated function. Use 'plot_electron_banddos' instead.",
                  stacklevel=2)
    if np.all(cmap_min!=None) and np.all(cmap_max!=None):
        levels = np.linspace(cmap_min, cmap_max, 150)
    else:
        levels = 150
    fig = plot_ECHG(obj_echg, unit=unit, levels=levels, option='spin',
                    x_ticks=xticks, y_ticks=yticks)
    return fig, fig.axes


def plot_cry_contour(contour_obj):
    """
    Deprecated. Use ``plot_topond2D``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_topond2D' instead.",
                  stacklevel=2)
    return contour_obj.plot_2D()


def plot_cry_contour_differences(contour_obj, contour_obj_ref):
    """
    Deprecated. Use ``plot_topond2D``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_topond2D' instead.",
                  stacklevel=2)
    return plot_topond2D(contour_obj, contour_obj_ref, option='diff')


def plot_cry_seebeck_potential(seebeck_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)
    if seebeck_obj.type != 'SEEBECK':
        raise TypeError('Not a SEEBECK object.')

    dir = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    fig = plot_transport_tensor(seebeck_obj, x_axis='potential', direction=dir)
    return fig


def plot_cry_seebeck_carrier(seebeck_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)
    if seebeck_obj.type != 'SEEBECK':
        raise TypeError('Not a SEEBECK object.')

    dir = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    fig = plot_transport_tensor(seebeck_obj, x_axis='carrier', direction=dir)
    return fig


def plot_cry_multiseebeck(*seebeck):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)

    k = int(input('Insert the index of temperature you want to plot \n(i.e. if your temperature are [T1, T2, T3] indexes are [0, 1, 2])'))
    minpot = float(input('Insert the lower value of chemical potential you want to plot in eV'))
    maxpot = float(input('Inser the higher value of chemical potential you want to plot in eV'))
    dir = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz\n')
    dir = dir[-2:].lower()

    objs = []
    for s in seebeck:
        if s.type == 'SEEBECK':
            objs.append(s)
        else:
            raise TypeError("Input not a SEEBECK object.")

    fig = plot_transport_tensor(*objs, option='multi', x_range=[minpot, maxpot],
                                direction=dir, plot_series=objs[0].T[k])
    return fig


def plot_cry_sigma_potential(sigma_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)

    if sigma_obj.type != 'SIGMA':
        raise TypeError('Not a SIGMA object.')

    dir = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    fig = plot_transport_tensor(sigma_obj, x_axis='potential', direction=dir)
    return fig


def plot_cry_sigma_carrier(sigma_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)

    if sigma_obj.type != 'SIGMA':
        raise TypeError('Not a SIGMA object.')

    dir = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    fig = plot_transport_tensor(sigma_obj, x_axis='carrier', direction=dir)


def plot_cry_multisigma(*sigma):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)

    k = int(input('Insert the index of temperature you want to plot \n(i.e. if your temperature are [T1, T2, T3] indexes are [0, 1, 2])'))
    minpot = float(input('Insert the lower value of chemical potential you want to plot in eV'))
    maxpot = float(input('Inser the higher value of chemical potential you want to plot in eV'))
    dir = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    objs = []
    for s in sigma:
        if s.type == 'SIGMA':
            objs.append(s)
        else:
            raise TypeError("Input not a SIGMA object.")

    fig = plot_transport_tensor(*objs, option='multi', x_range=[minpot, maxpot],
                                direction=dir, plot_series=objs[0].T[k])
    return fig


def plot_cry_powerfactor_potential(seebeck_obj, sigma_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)

    dir = input('Please, choose the direction you want to plot. \nYou can choose among PF_xx, PF_xy, PF_xz, PF_yx, PF_yy, PF_yz, PF_yz, PF_zx, PF_zy, PF_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    fig = plot_transport_tensor([seebeck_obj, sigma_obj], option='normal', direction=dir)
    return fig


def plot_cry_powerfactor_carrier(seebeck_obj, sigma_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)

    dir = input('Please, choose the direction you want to plot. \nYou can choose among PF_xx, PF_xy, PF_xz, PF_yx, PF_yy, PF_yz, PF_yz, PF_zx, PF_zy, PF_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    fig = plot_transport_tensor([seebeck_obj, sigma_obj], option='normal',
                                direction=dir, x_axis='carrier')
    return fig

def plot_cry_zt(seebeck_obj, sigma_obj):
    """
    Deprecated. Use ``plot_transport_tensor``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_transport_tensor' instead.",
                  stacklevel=2)
    warnings.warn("The y axis captions will be displayed in Power Factor. The plot is divided by Kappa.")

    ktot = float(input('Please insert the value of ktot in WK-1m-1'))
    dir = input('Please, choose the direction you want to plot. \nYou can choose among ZT_xx, ZT_xy, ZT_xz, ZT_yx, ZT_yy, ZT_yz, ZT_yz, ZT_zx, ZT_zy, ZT_zz\n')
    dir = dir[-2:].lower()
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    for irow in range(len(sigma_obj.T)):
        sigma_obj.data[irow] = sigma_obj.data[irow] * sigma_obj.T[irow]
    sigma_obj.data = sigma_obj.data / ktot
    fig = plot_transport_tensor([seebeck_obj, sigma_obj], option='normal', direction=dir)
    return fig


