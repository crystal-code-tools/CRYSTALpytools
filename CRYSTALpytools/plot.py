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

def plot_ECHG(*echg, unit='Angstrom', option='both', levels=150, lineplot=False,
              linewidth=1.0, isovalues=None, colorplot=True, colormap='jet',
              cbar_label=None, a_range=[], b_range=[], rectangle=False,
              edgeplot=False, x_ticks=5, y_ticks=5, add_title=True, figsize=[6.4, 4.8],
              **kwargs):
    """
    Read multiple 2D charge density files / objects and return to a list of
    figures and axes. The uniform plot set-ups are used for comparison.

    Available options:

    * 'both' : If spin polarized, plot both charge and spin densities.
        Otherwise plot charge densities.  
    * 'charge': Plot charge density.  
    * 'spin': Plot spin density.  
    * 'diff': Substracting charge data from the first entry with the following
        entries. Return to a non spin-polarized object.  

    .. note::

        If file names are given, the code only reads the first 2D data map in
        fort.25 files.

    Args:
        \*echg (ChargeDensity|str): Extendable. File names or
            ``electronics.ChargeDensity`` objects.
        unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
            Bohr:math:`^{-3}`.
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
            'None' for default.
        a_range (list): 1\*2 range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): 1\*2 range of :math:`b` axis (x, or AB) in fractional coordinate.
        rectangle (bool): If :math:`a, b` are non-orthogonal, plot a rectangle
            region and reset :math:`b`. If used together with ``b_range``, that
            refers to the old :math:`b`(i.e., expansion first).
        edgeplot (bool): Whether to add cell edges represented by the original
            base vectors (not inflenced by a/b range or rectangle options).
        x_ticks (int): Number of ticks on x axis.
        y_ticks (int): Number of ticks on y axis.
        add_title (bool): Whether to add property plotted as title.
        figsize (list): Matplotlib figure size. Note that axes aspects are
            fixed to be equal.
        \*\*kwargs : Other arguments passed to ``axes.contour()`` function
            to set contour lines.
    Returns:
        figs (list|Figure): Matplotlib Figure object or a list of them.
    """
    from CRYSTALpytools.electronics import ChargeDensity
    import numpy as np

    obj = []
    for i in echg:
        if isinstance(i, str):
            obj.append(ChargeDensity.from_file(i, method=None))
        elif isinstance(i, ChargeDensity):
            obj.append(i)
        else:
            raise TypeError("Inputs must be either string or electronics.ChargeDensity objects.")
    # substraction
    if 'diff' in option.lower():
        obj[0].substract(*[i for i in obj[1:]])
        option = 'charge'
        obj = [obj[0]]
    # set uniform levels
    if isinstance(levels, float) or isinstance(levels, int):
        spin_range = []
        chg_range = []
        for i in obj:
            chg_range.append([np.min(i.data[:, :, 0]), np.max(i.data[:, :, 0])])
            if i.spin == 2:
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

    # plot
    figs = []
    for i in obj:
        figs.append(
            i.plot_2D(unit, option, levels, lineplot, linewidth, isovalues,
                      colorplot, colormap, cbar_label, a_range, b_range,
                      rectangle, edgeplot, x_ticks, y_ticks, add_title,
                      figsize, **kwargs)
        )
    if len(obj) == 1:
        figs = figs[0]
    return figs

#----------------------------------SPIN CURRENTS------------------------------#

def plot_relativistics2D(*relat, unit='SI', type=[], output=[], direction=['x','y','z'],
                         levels=100, quiverplot=True, quiverscale=1.0,
                         colorplot=True, colormap='jet', cbar_label='default',
                         a_range=[], b_range=[], rectangle=False, edgeplot=False,
                         x_ticks=5, y_ticks=5, layout=None, title=None,
                         figsize=[6.4, 4.8], sharex=True, sharey=True, **kwargs):
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
        unit (str): Plot unit. 'SI' for :math:`\\AA` and A/m (A/m:math:`^{2}`).
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
            default subplot titles are used.
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
                    raise ValueError("Indices must be set for input files. Otherwise use input objects.")
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
                           sharex=sharex, sharey=sharey, layout='constrained')
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


def plot_vecfield2D_m(header, dens, quivscale, name='MAG', levels=150, dpi=400):
    """
    Plots the 2D magnetization vector field.

    Args:
        header (list): List containing information about the fort.25 header.
        dens (numpy.ndarray): Array containing the vector field data.
        quivscale (float): Scale factor for the quiver plot.
        name (str, optional):  Name used for saving the plots.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions.
        dpi (int, optional): DPI (dots per inch) for the output image. Default is 400.

    Returns:
        None
    """
    import matplotlib.pyplot as plt
    import numpy as np

    for i in range(2, 20):
        if (header[0] % i) == 0:
            nrow_split = int(header[0]/i)

    for i in range(2, 20):
        if (header[1] % i) == 0:
            ncol_split = int(header[1]/i)

    # initializes the arrays
    projx_m = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_m = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_m = np.zeros((header[0], header[1]), dtype=float)

    # Generates the meshgrid
    mesh_x = np.zeros((header[0], header[1]), dtype=float)
    mesh_y = np.zeros((header[0], header[1]), dtype=float)
    mesh_projx = np.zeros((nrow_split, ncol_split), dtype=float)
    mesh_projy = np.zeros((nrow_split, ncol_split), dtype=float)

    T = np.array([[np.sqrt(1 - header[4]**2)*header[2], 0],
                  [header[4]*header[2], header[3]]])  # Change of basis matrix

    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mesh_x[i, j] = np.dot(T, np.array([i, j]))[0]
            mesh_y[i, j] = np.dot(T, np.array([i, j]))[1]

    mesh_x = mesh_x * 0.529177249  # Bohr to Angstrom conversion
    mesh_y = mesh_y * 0.529177249

    r = 0
    s = 0
    step_nrow = int(header[0]/nrow_split)
    step_ncol = int(header[1]/ncol_split)
    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            mesh_projx[i, j] = mesh_x[r, s]
            mesh_projy[i, j] = mesh_y[r, s]
            s += step_ncol
        s = 0
        r += step_nrow

    # Creates the orthogonal vectorial basis for the projections
    CB = header[7] - header[6]
    BA = header[6] - header[5]
    if header[4] == 0:
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB
    else:
        BA = BA - ((np.dot(BA, CB))/(np.dot(BA, BA)))*CB
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB

    # Calculates the modulus of the vectorial field
    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mod_m[i, j] = np.sqrt((dens[i, j, 0]**2) +
                                  (dens[i, j, 1]**2)+(dens[i, j, 2]**2))

    # Calculates the projections of the vectors in the ABC plane
    ABC_normal = np.cross(BA, CB)
    mod_normal = np.sqrt(ABC_normal[0]**2 +
                         ABC_normal[1]**2 + ABC_normal[2]**2)

    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            projx_m[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol), 0], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_m[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol, 0)], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)

    # Plotting

    m = plt.figure()
    m = plt.contourf(mesh_x, mesh_y, mod_m, levels, cmap='cool')
    m = plt.colorbar(mappable=m)
    m = plt.quiver(mesh_projx, mesh_projy, projx_m, projy_m, scale=quivscale)
    m = plt.xlabel('$\AA$')
    m = plt.ylabel('$\AA$')
    m = plt.savefig(name, dpi=dpi)

    plt.show()


def plot_vecfield2D_j(header, dens, quivscale, name='SC', levels=150, dpi=400):
    """
    Plots the 2D vector field of the spin current.

    Args:
        header (list): List containing information about the fort.25 header.
        dens (numpy.ndarray): Array representing the vector field.
        quivscale (float): Scale factor for the quiver plot.
        name (str, optional):  Name used for saving the plots.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions.
        dpi (int, optional): DPI (dots per inch) for the output image. Defaults to 400.

    Returns:
        None
    """
    import matplotlib.pyplot as plt
    import numpy as np

    for i in range(2, 20):
        if (header[0] % i) == 0:
            nrow_split = int(header[0]/i)

    for i in range(2, 20):
        if (header[1] % i) == 0:
            ncol_split = int(header[1]/i)

    # initializes the arrays
    projx_j = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_j = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_j = np.zeros((header[0], header[1]), dtype=float)

    # Generates the meshgrid
    mesh_x = np.zeros((header[0], header[1]), dtype=float)
    mesh_y = np.zeros((header[0], header[1]), dtype=float)
    mesh_projx = np.zeros((nrow_split, ncol_split), dtype=float)
    mesh_projy = np.zeros((nrow_split, ncol_split), dtype=float)

    T = np.array([[np.sqrt(1 - header[4]**2)*header[2], 0],
                  [header[4]*header[2], header[3]]])  # Change of basis matrix

    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mesh_x[i, j] = np.dot(T, np.array([i, j]))[0]
            mesh_y[i, j] = np.dot(T, np.array([i, j]))[1]

    mesh_x = mesh_x * 0.529177249  # Bohr to Angstrom conversion
    mesh_y = mesh_y * 0.529177249

    r = 0
    s = 0
    step_nrow = int(header[0]/nrow_split)
    step_ncol = int(header[1]/ncol_split)
    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            mesh_projx[i, j] = mesh_x[r, s]
            mesh_projy[i, j] = mesh_y[r, s]
            s += step_ncol
        s = 0
        r += step_nrow

    # Creates the orthogonal vectorial basis for the projections
    CB = header[7] - header[6]
    BA = header[6] - header[5]
    if header[4] == 0:
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB
    else:
        BA = BA - ((np.dot(BA, CB))/(np.dot(BA, BA)))*CB
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB

    # Calculates the modulus of the vectorial field
    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mod_j[i, j] = np.sqrt((dens[i, j, 0]**2) +
                                  (dens[i, j, 1]**2)+(dens[i, j, 2]**2))

    # Calculates the projections of the vectors in the ABC plane
    ABC_normal = np.cross(BA, CB)
    mod_normal = np.sqrt(ABC_normal[0]**2 +
                         ABC_normal[1]**2 + ABC_normal[2]**2)

    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            projx_j[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol), 0], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_j[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                                                np.cross(np.array([dens[int(i*step_nrow), int(j*step_ncol), 0], dens[int(i*step_nrow), int(j*step_ncol), 1],
                                                                                   dens[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)

    j = plt.figure()
    j = plt.contourf(mesh_x, mesh_y, mod_j, levels, cmap='winter')
    j = plt.colorbar(mappable=j)
    j = plt.quiver(mesh_projx, mesh_projy, projx_j, projy_j, scale=quivscale)
    j = plt.xlabel('$\AA$')
    j = plt.ylabel('$\AA$')
    j = plt.savefig(name, dpi=dpi)

    plt.show()


def plot_vecfield2D_J(header, dens_JX, dens_JY, dens_JZ, quivscale, name='SCD', levels=150, dpi=400):
    """
    Plots the 2D spin current density vector fields.

    Args:
        header (list): List containing information about the fort.25 header.
        dens_JX (numpy.ndarray): Array representing the X-component of the spin current density.
        dens_JY (numpy.ndarray): Array representing the Y-component of the spin current density.
        dens_JZ (numpy.ndarray): Array representing the Z-component of the spin current density.
        quivscale: Scale factor for the quiver plot.
        name (str, optional): Name used for saving the plots.
        levels (int or array-like, optional): Determines the number and positions of the contour lines/regions.
        dpi (int, optional): DPI (Dots per inch) for saving the plots. Defaults to 400.

    Returns:
        None
    """
    import matplotlib.pyplot as plt
    import numpy as np

    for i in range(2, 20):
        if (header[0] % i) == 0:
            nrow_split = int(header[0]/i)

    for i in range(2, 20):
        if (header[1] % i) == 0:
            ncol_split = int(header[1]/i)

    # initializes the arrays
    projx_JX = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_JX = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_JX = np.zeros((header[0], header[1]), dtype=float)
    projx_JY = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_JY = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_JY = np.zeros((header[0], header[1]), dtype=float)
    projx_JZ = np.zeros((nrow_split, ncol_split), dtype=float)
    projy_JZ = np.zeros((nrow_split, ncol_split), dtype=float)
    mod_JZ = np.zeros((header[0], header[1]), dtype=float)

    # Generates the meshgrid
    mesh_x = np.zeros((header[0], header[1]), dtype=float)
    mesh_y = np.zeros((header[0], header[1]), dtype=float)
    mesh_projx = np.zeros((nrow_split, ncol_split), dtype=float)
    mesh_projy = np.zeros((nrow_split, ncol_split), dtype=float)

    T = np.array([[np.sqrt(1 - header[4]**2)*header[2], 0],
                  [header[4]*header[2], header[3]]])  # Change of basis matrix

    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mesh_x[i, j] = np.dot(T, np.array([i, j]))[0]
            mesh_y[i, j] = np.dot(T, np.array([i, j]))[1]

    mesh_x = mesh_x * 0.529177249  # Bohr to Angstrom conversion
    mesh_y = mesh_y * 0.529177249

    r = 0
    s = 0
    step_nrow = int(header[0]/nrow_split)
    step_ncol = int(header[1]/ncol_split)
    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            mesh_projx[i, j] = mesh_x[r, s]
            mesh_projy[i, j] = mesh_y[r, s]
            s += step_ncol
        s = 0
        r += step_nrow

    # Creates the orthogonal vectorial basis for the projections
    CB = header[7] - header[6]
    BA = header[6] - header[5]
    if header[4] == 0:
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB
    else:
        BA = BA - ((np.dot(BA, CB))/(np.dot(BA, BA)))*CB
        mod_BA = np.sqrt(BA[0]**2 + BA[1]**2 + BA[2]**2)
        mod_CB = np.sqrt(CB[0]**2 + CB[1]**2 + CB[2]**2)
        v1 = BA/mod_BA
        v2 = CB/mod_CB

    # Calculates the modulus of the vectorial field
    for i in range(0, header[0]):
        for j in range(0, header[1]):
            mod_JX[i, j] = np.sqrt(
                (dens_JX[i, j, 0]**2)+(dens_JX[i, j, 1]**2)+(dens_JX[i, j, 2]**2))
            mod_JY[i, j] = np.sqrt(
                (dens_JY[i, j, 0]**2)+(dens_JY[i, j, 1]**2)+(dens_JY[i, j, 2]**2))
            mod_JZ[i, j] = np.sqrt(
                (dens_JZ[i, j, 0]**2)+(dens_JZ[i, j, 1]**2)+(dens_JZ[i, j, 2]**2))

    # Calculates the projections of the vectors in the ABC plane
    ABC_normal = np.cross(BA, CB)
    mod_normal = np.sqrt(ABC_normal[0]**2 +
                         ABC_normal[1]**2 + ABC_normal[2]**2)

    for i in range(0, nrow_split):
        for j in range(0, ncol_split):
            projx_JX[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JX[int(i*step_nrow), int(j*step_ncol), 0], dens_JX[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JX[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_JX[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JX[int(i*step_nrow), int(j*step_ncol), 0], dens_JX[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JX[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)
            projx_JY[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JY[int(i*step_nrow), int(j*step_ncol), 0], dens_JY[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JY[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_JY[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JY[int(i*step_nrow), int(j*step_ncol), 0], dens_JY[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JY[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)
            projx_JZ[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JZ[int(i*step_nrow), int(j*step_ncol), 0], dens_JZ[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JZ[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v2)
            projy_JZ[i, j] = np.dot((1/(mod_normal**2))*np.cross(ABC_normal,
                                    np.cross(np.array([dens_JZ[int(i*step_nrow), int(j*step_ncol), 0], dens_JZ[int(i*step_nrow), int(j*step_ncol), 1],
                                                       dens_JZ[int(i*step_nrow), int(j*step_ncol), 2]]), ABC_normal)), v1)

    # Plotting

    JX = plt.figure()
    JX = plt.contourf(mesh_x, mesh_y, mod_JX, levels, cmap='summer')
    JX = plt.colorbar(mappable=JX)
    JX = plt.quiver(mesh_projx, mesh_projy, projx_JX,
                    projy_JX, scale=quivscale)
    JX = plt.xlabel('$\AA$')
    JX = plt.ylabel('$\AA$')
    JX = plt.savefig(name+'_JX', dpi=dpi)

    JY = plt.figure()
    JY = plt.contourf(mesh_x, mesh_y, mod_JY, levels, cmap='summer')
    JY = plt.colorbar(mappable=JY)
    JY = plt.quiver(mesh_projx, mesh_projy, projx_JY,
                    projy_JY, scale=quivscale)
    JY = plt.xlabel('$\AA$')
    JY = plt.ylabel('$\AA$')
    JY = plt.savefig(name+'_JY', dpi=dpi)

    JZ = plt.figure()
    JZ = plt.contourf(mesh_x, mesh_y, mod_JZ, levels, cmap='summer')
    JZ = plt.colorbar(mappable=JZ)
    JZ = plt.quiver(mesh_projx, mesh_projy, projx_JZ,
                    projy_JZ, scale=quivscale)
    JZ = plt.xlabel('$\AA$')
    JZ = plt.ylabel('$\AA$')
    JZ = plt.savefig(name+'_JZ', dpi=dpi)

    plt.show()


#--------------------------------BAND STRUCTURES------------------------------#

def plot_electron_bands(*bands, unit='eV', k_label=[], mode='single',
                        not_scaled=False, energy_range=[], k_range=[],
                        band_label=None, band_color=None, band_linestyle=None,
                        band_linewidth=None, fermi_level=0., fermi_color='tab:green',
                        fermi_linestyle='-', fermi_linewidth=1.0, layout=None,
                        title=None, figsize=[6.4, 4.8], legend='lower right',
                        sharex=True, sharey=True, fontsize=14, **kwargs):
    """
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
            [ax], [bandsplt[0].bands], [bandsplt[0].tick_pos],
            k_label, False, energy_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, fermi_level, fermi_color,
            fermi_linestyle, fermi_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'multi':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = plot_overlap_bands(
            ax, [b.bands for b in bandsplt], [b.tick_pos for b in bandsplt],
            k_label, energy_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, fermi_level, fermi_color,
            fermi_linestyle, fermi_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'compare':
        if np.all(layout==None):
            layout = [int(np.ceil(len(bandsplt)/2)), 2]

        fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        _ = plot_compare_bands(
            ax.flat, [b.bands for b in bandsplt], [b.tick_pos for b in bandsplt],
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
            same color. 'None' for default values ('tab:blue' and other tab
            series for both spin-up and spin-down states). Effective for all
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
        if ndoss == 1:
            ax = plot_doss(
                ax, dossplt[0].doss, dossplt[0].energy, beta, prj[0],
                energy_range, dos_range, dos_label, dos_color, dos_linestyle,
                dos_linewidth, fermi_level, fermi_color, fermi_linestyle,
                fermi_linewidth, legend, False, **kwargs
            )
        else:
            for i in range(ndoss):
                ax.flat[i] = plot_doss(
                    ax.flat[i], dossplt[i].doss, dossplt[i].energy, beta, prj[i],
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
                  colorplot=False, colormap='jet', cbar_label=None,
                  cpt_marker='o', cpt_color='k', cpt_size=10,
                  traj_color='r', traj_linestyle=':', traj_linewidth=0.5,
                  a_range=[], b_range=[], edgeplot=False, x_ticks=5, y_ticks=5,
                  add_title=True, figsize=[6.4, 4.8]):
    """
    Read multiple TOPOND 2D plot files / objects and return to a list of
    figures and axes. The uniform plot set-ups are used for comparison.

     .. note::

        For the convenience of analysis and plotting, it is important to select
        the correct type for your input file. By default `type='infer'` will
        search for (case insensitive) the following strings:

        'SURFRHOO', 'SURFSPDE', 'SURFLAPP', 'SURFLAPM', 'SURFGRHO', 'SURFKKIN',
        'SURFGKIN', 'SURFVIRI', 'SURFELFB', 'TRAJGRAD', 'TRAJMOLG'

        For their meanings, please refer the `TOPOND manual <https://www.crystal.unito.it/include/manuals/topond.pdf>`_.

    Available options:

    * 'normal' : Literally normal.  
    * 'diff' : Substract data from the first entry using following entries. All
        the entries must have the same ``type`` and must be 'SURF*' types.  
    * 'overlay': Overlapping a 'TRAJ*' object on the 2D 'SURF*' object. Inputs
        must be 1\*2 lists of a ``Surf`` object and a ``Traj`` object. File
        names are not permitted as geometry information is mandatory.

    .. note::

        2D periodicity (``a_range`` and ``b_range``), though available for the
        ``Surf`` class, is not suggested as TOPOND plotting window does not
        always commensurate with periodic boundary. The ``Traj`` class has no
        2D periodicity so if ``option='overlay'``, ``a_range``, ``b_range`` and
        ``edgeplot`` will be disabled.

    Args:
        \*topond (Surf|Traj|str): Extendable. File names, ``topond.Surf`` or
            ``topond.Traj`` objects. Geometry information is not available if
            file names are used - **might lead to errors!**.
        unit (str): Plot unit. 'Angstrom' for :math:`\\AA^{-3}`, 'a.u.' for
            Bohr:math:`^{-3}`.
        type (str): 'infer' or specified. Otherwise warning will be given.
        option (str): Available options see above.
        levels (array): Set levels of contour plot. 'Default' for built-in,
                property adaptive levels (``unit='Angstrom'``). Otherwise
                entries **must be consistent with ``unit``**.
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
        add_title (bool): Whether to add property plotted as title.
        figsize (list): Matplotlib figure size. Note that axes aspects are
            fixed to be equal.
    Returns:
        figs (list|Figure): Matplotlib Figure object or a list of them.
    """
    from CRYSTALpytools.topond import Surf, Traj
    from CRYSTALpytools.crystal_io import Properties_output
    import numpy as np
    import warnings

    obj = []
    for i in topond:
        if isinstance(i, str):
            obj.append(Properties_output().read_topond(i, type=type))
        elif isinstance(i, Surf) or isinstance(i, Traj):
            obj.append(i)
        elif isinstance(i, list) or isinstance(i, tuple):
            if len(i) != 2:
                raise ValueError('Input lists must have 2 elements.')
            if isinstance(i[0], Surf) and isinstance(i[1], Traj):
                obj.append([i[0], i[1]])
            elif isinstance(i[1], Surf) and isinstance(i[0], Traj):
                obj.append([i[1], i[0]])
            else:
                raise TypeError("Input type does not follow the requirement of the 'overlay' option.")
        else:
            raise TypeError("Input type does not meet the requirements.")

    # substraction
    if 'diff' in option.lower():
        if isinstance(obj[0], Traj):
            raise TypeError("The 'diff' option is not applicable to 'topond.Traj' objects.")
        for i in obj[1:]:
            if obj[0].type != i.type:
                raise TypeError("Different properties are read for input objects / files, 'diff' option not available.")
            obj[0].substract(i)
        obj = [obj[0]]
    # set uniform levels
    if np.all(levels=='default'):
        levels = 'default'
    else:
        levels = np.array(levels, dtype=float)

    # plot
    figs = []
    if 'overlay' in option.lower():
        if a_range != [] or b_range != []:
            warnings.warn("Periodic plotting not available for 'topond.Traj' objects. Using default ranges.",
                          stacklevel=2)
            a_range = []; b_range = []
        for i in obj:
            figs.append(i[0].plot(
                unit, levels, lineplot, linewidth, isovalues, colorplot, colormap,
                cbar_label, [], [], False, x_ticks, y_ticks, add_title, figsize,
                overlay=i[1], cpt_marker=cpt_marker, cpt_color=cpt_color,
                cpt_size=cpt_size, traj_color=traj_color,
                traj_linestyle=traj_linestyle, traj_linewidth=traj_linewidth))

    else:
        for i in obj:
            if isinstance(i, Surf):
                figs.append(i.plot(
                    unit, levels, lineplot, linewidth, isovalues, colorplot,
                    colormap, cbar_label, a_range, b_range, edgeplot, x_ticks,
                    y_ticks, add_title, figsize, None))
            else:
                figs.append(i.plot_2D(
                    unit, cpt_marker, cpt_color, cpt_size, traj_color, traj_linestyle,
                    traj_linewidth, x_ticks, y_ticks, add_title, figsize, None, None))

    if len(obj) == 1:
        figs = figs[0]
    return figs

#--------------------------------------XRD------------------------------------#


def plot_cry_xrd(xrd_obj):
    """
    Plot the X-ray diffraction pattern.

    Args:
        xrd_obj (object): XRD object containing the data for the X-ray diffraction pattern.
        save_to_file (bool, optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Plots the X-ray diffraction pattern.
        - Sets the figure size to [16, 9].
        - Sets the x-axis limit to (0, 30).
        - Saves the plot to a file named 'figure_XRD_YYYY-MM-DD_HHMMSS.jpg' in the current directory.
        - If save_to_file is True, saves the plot to a file specified by save_to_file parameter.

    """
    import os
    import time

    import matplotlib.pyplot as plt

    plt.rcParams["figure.figsize"] = [16, 9]

    plt.plot(xrd_obj.x, xrd_obj.y)

    plt.xlim((0, 30))

    path = os.path.join('./'+'figure_'+'XRD_' +
                        time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(xrd_obj.title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=600)

    # if save_to_file != False:
    #     save_plot(save_to_file)

    plt.show()

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
        - Sets the x-axis label as 'd  [$\AA$]' and the y-axis label as r'$\rho$  [$\frac{e}{\AA^3}$]'.
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

    1. Normal: Same as ``transport.Tensor.plot()``. Only 1 entry permitted.  
    2. Power Factor: Power factor :math:`S^{2}\\sigma`. Limited to 2 entries in
        'SIGMA', 'SIGMAS' and 'SEEBECK'.  
    3. ZT: Dimensionless figure of merit :math:`\\frac{S^{r}\\sigma T}{\\kappa}`.
        Limited to 3 entries, 1 of 'KAPPA' and 2 in 'SIGMA', 'SIGMAS' and
        'SEEBECK'.  
    4. Multi: Plot transport properties of different systems together for
        comparison. In this case ``direction`` does not accept list variables,
        and entries of ``plot_series`` will be plotted into subplots.

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
            is not allowed for ``option='multi'``.
        spin (str): Spin-polarized systems only. Electron spins to plot. 'up',
            'down' or 'sum'. Disabled for ``option='spin'``.
        plot_series (list|array|float): **Values** of the plot series. Should
            be commensurate with the choice of ``x_axis``. It will be annotated
            on the upper right corner for ``option='multi'``.
        plot_label (list|str|None): None for default values: ``plot_series``
            and its unit, or sequence number of materials for ``option='multi'``.
            If str, prefixing that entry in front of the default value. If list,
            it should be 1\*nPlotSeries, or 1\*nPlot_series for ``option='multi'``,
            for every plot line.
        plot_color (list|str|None): Similar to ``electronics.ElectronDOS.plot()``.
            If str, use the same color for all the plot lines. If
            1\*nPlotSeries(nPlot_series), use the same color for every plot line.
            If 2\*nPlotSeries(nPlot_series), use different colors for p-type and
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

    objs = []
    for b in boltztra:
        if isinstance(b, str):
            objs.append(Tensor.from_file(b))
        elif isinstance(b, Tensor):
            objs.append(b)
        else:
            raise TypeError("Inputs must be filename or 'transport.Tensor' object.")
    # options
    if option.lower() == 'normal':
        if len(objs) != 1:
            raise ValueError("For 'normal' option, Only 1 entry is needed.")
    elif option.lower() == 'power factor':
        if len(objs) != 2:
            raise ValueError("For 'power factor' option, 2 entries are needed.")
        objs = [Tensor.get_power_factor(objs[0], objs[1])]
    elif option.lower() == 'zt':
        if len(objs) != 3:
            raise ValueError("For 'zt' option, 3 entries are needed.")
        objs = [Tensor.get_zt(objs[0], objs[1], objs[2])]
    elif option.lower() == 'multi':
        if not isinstance(direction, str):
            raise ValueError("For 'multi' option, only one direction in string should be given.")
        for i in range(1, len(objs)):
            if objs[i].type != objs[0].type:
                raise TypeError("For 'multi' option, input entries must have the same type.")
    else:
        raise ValueError("Unknown options: '{}'.".format(option))

    if option.lower() != 'multi':
        fig = objs[0].plot(x_axis, x_range, direction, spin, plot_series,
                           plot_label, plot_color, plot_linestyle, plot_linewidth,
                           zero_color, zero_linestyle, zero_linewidth, layout, False,
                           figsize, legend, sharey, fontsize, **kwargs)
        if np.all(title!=None):
            fig.suptitle(title, fontsize=fontsize)
    else:
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
        if nplt == 0:
            raise Exception('Cannot find common plot series.')
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
        if objs[0].data.shape[2] == 9:
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
            plot_label = ['# {:d}'.format(i) for i in range(nprj)]
        elif isinstance(plot_label, str):
            plot_label = ['{} {:d}'.format(plot_label, i) for i in range(nprj)]
        elif isinstance(plot_label, list) or isinstance(plot_label, array):
            nlabel = len(plot_label)
            plot_label = [plot_label[i%nlabel] for i in range(nprj)]
        else:
            raise TypeError('Unknown type of plot label.')
        ## get a pseudo band input
        bands = np.zeros([nprj, 1, 1, 1], dtype=float)
        commands = _plot_label_preprocess(bands, plot_label, plot_color,
                                          plot_linestyle, plot_linewidth)

        # plot layout
        if np.all(layout==None):
            fig, ax = plt.subplots(nplt, 1, sharex=True, sharey=sharey,
                                   figsize=figsize, layout='constrained')
        else:
            fig, ax = plt.subplots(layout[0], layout[1], sharex=True,
                                   sharey=sharey, figsize=figsize, layout='constrained')

        # plot every subplot
        y_range = []
        for iplt, ax in enumerate(fig.axes):
            y_range_plt = []
            ax.hlines(0, x_range[0], x_range[1], colors=zero_color,
                      linestyle=zero_linestyle, linewidth=zero_linewidth)
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
            if np.all(legend!=None) and np.all(commands[0]!=None):
                ax.legend(loc=legend)
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
def plot_elastics_3D(property, *tensor, uniform_scale=True, add_title=True,
                     **kwargs):
    """
    A wrapper function of :ref:`Tensor3D <ref-elastics>` objects to plot 3D
    crystal elastic properties. The user can plot multiple properties for
    different systems. The function returns to lists of figure objects for
    further processing. Only matplotlib is used for plotting. Plotly is not
    available.

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
        property (str | list[str]): The properties to plot. See above.
        \*tensor (str | Tensor3D | numpy.ndarray): Elastic tensor definition.
            Can be CRYSTAL output files, ``Tensor3D`` objects and 6\*6
            **elastic** matrices in Voigt notation, GPa. For files,
            ``conventional_lattice=True``.
        uniform_scale (bool): Use the same color scale for all plots of the
            same property.
        add_title (bool): Add properties as titles of subplots.
        \*\*kwargs : Plot settings. The settings will be used for all plots.
            Check :ref:`Tensor3D.plot_3D <ref-elastics>`

    Returns:
        figs (Figure): For 1 tensor, 1 property, return to a matplotlib Figure
            object. For 1 tensor and multiple properties or multiple tensors
            and 1 property, n\*1 list of Figure objects. For multiple tensors
            and properties, nTensor\*nProperty list.
    """
    from CRYSTALpytools.elastics import Tensor3D, _plot3D_mplib
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

    if isinstance(property, str):
        property = [property]

    if len(tensplt) == 1 and len(property) > 1 and uniform_scale == True:
        warnings.warn("'uniform_scale' cannot be used for multiple proeprties of the same system. Using 'uniform_scale = False'.",
                      stacklevel=2)
        uniform_scale = False

    if 'plot_lib' in kwargs.keys():
        if kwargs['plot_lib'] == 'plotly':
            warnings.warn("Only matplotlib is available.", stacklevel=2)
    kwargs['plot_lib'] = 'matplotlib'

    # plot
    n_plot = len(tensplt) * len(property)
    figs = [[0 for i in range(len(property))] for j in range(len(tensplt))]

    if uniform_scale == False: # Non-uniform scale
        for ip, p in enumerate(property):
            kwargs['property'] = p
            for it, t in enumerate(tensplt):
                fig = t.plot_3D(**kwargs)
                if add_title == True:
                    fig.axes[1].set_title(p)
                figs[it][ip] = fig
    else: # Uniform scale
        # possible camera position args
        camera_args = {}
        for i in ['elev', 'azim', 'roll', 'vertical_axis', 'share']:
            try:
                camera_args[i] = kwargs[i]
            except KeyError:
                pass

        for ip, p in enumerate(property):
            R_all = []
            X_all = []
            Y_all = []
            Z_all = []
            uplt_all =[]
            utext_all = []
            platt_all = []
            kwargs['property'] = p
            kwargs['plot_lib'] = None
            # Get data first
            for it, t in enumerate(tensplt):
                R, X, Y, Z, scale_radius, uplt, utext, platt = t.plot_3D(**kwargs)
                R_all.append(R)
                X_all.append(X)
                Y_all.append(Y)
                Z_all.append(Z)
                uplt_all.append(uplt)
                utext_all.append(utext)
                platt_all.append(platt)
            # plot
            try: range_cbar = kwargs['range_cbar']
            except KeyError: range_cbar = [np.min(R_all), np.max(R_all)]
            try: range_x = kwargs['range_x']
            except KeyError: range_x = None
            try: range_y = kwargs['range_y']
            except KeyError: range_y = None
            try: range_z = kwargs['range_z']
            except KeyError: range_z = None
            if np.all(platt_all[0]!=None):
                Rref = R_all[
                    np.argmax(np.array(R_all).flatten()) // \
                    (np.shape(R_all[0])[0]* np.shape(R_all[0])[1])
                ]
            else:
                Rref = None

            for it, t in enumerate(tensplt):
                fig = _plot3D_mplib(
                    R_all[it], X_all[it], Y_all[it], Z_all[it], scale_radius,
                    uplt_all[it], utext_all[it], platt_all[it], range_cbar,
                    range_x, range_y, range_z, Rref, **camera_args
                )
                if add_title == True:
                    fig.axes[1].set_title(p)
                figs[it][ip] = fig

    # dimensionality
    if len(tensplt) == 1 and len(property) == 1:
        figs = figs[0][0]
    elif len(tensplt) == 1 and len(property) > 1:
        figs = figs[0]
    elif len(tensplt) > 1 and len(property) == 1:
        figs = [i[0] for i in figs]

    return figs


#--------------------------------2D ELASTIC----------------------------------#
def plot_elastics_2D(property, *tensor, same_fig_2D=True, uniform_scale_2D=True, **kwargs):
    """
    A wrapper function of :ref:`Tensor3D or Tensor2D <ref-elastics>` objects to
    plot 2D crystal elastic properties. The user can plot multiple properties
    for different systems. The function returns to lists of figure objects for
    further processing. Base units: GPa, m.

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

    Args:
        property (str | list[str]): The properties to plot. See above.
        \*tensor (str | Tensor3D | numpy.ndarray): Elastic tensor definition.
            Can be CRYSTAL output files, ``Tensor3D`` / ``Tensor2D``objects
            and 6\*6 / 3\*3 **elastic** matrices in Voigt notation, GPa. But
            dimensionalities of systems must be consistent. For 3D files,
            ``conventional_lattice=True``.
        same_fig_2D (bool): *Valid only for 2D systems*. Plot all the pole
            chart into the same figure. Same dimension as the ``figs`` list.
        uniform_scale_2D (bool): *Valid only for 2D systems*. Use the same
            radial scale for all plots of the same property.  ``uniform_scale``
            can be enabled for different planes of the same system in 3D cases.
        \*\*kwargs : Plot settings. The settings will be used for all plots.
            Check :ref:`Tensor3D.plot_2D and Tensor2D.plot_2D <ref-elastics>`.

    Returns:
        figs (Figure): For 1 tensor, 1 property, return to a matplotlib Figure
            object. For 1 tensor and multiple properties or multiple tensors
            and 1 property, n\*1 list of Figure objects. For multiple tensors
            and properties, nTensor\*nProperty list. If ``same_fig_2D=True``,
            return to a figure.
    """
    from CRYSTALpytools.elastics import Tensor3D, Tensor2D, tensor_from_file, _plot2D_single
    import numpy as np
    import warnings
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

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

    if isinstance(property, str):
        property = [property]

    if len(tensplt) == 1 and len(property) > 1 and uniform_scale_2D == True:
        warnings.warn("'uniform_scale_2D' cannot be used for multiple proeprties of the same system. Using 'uniform_scale_2D = False'.",
                      stacklevel=2)
        uniform_scale = False


    n_plot = len(tensplt) * len(property)
    if isinstance(tensplt[0], Tensor3D) or n_plot == 1.:
        same_fig_2D = False
        uniform_scale_2D = False

    # set subfigures
    if same_fig_2D == True:
        if len(property) == 1 and len(tensplt) > 1:
            figs, axes = plt.subplots(nrows=1, ncols=len(tensplt),
                                      subplot_kw={'projection' : 'polar'}, layout='tight')
        else:
            figs, axes = plt.subplots(nrows=len(tensplt), ncols=len(property),
                                      subplot_kw={'projection' : 'polar'}, layout='tight')
    else:
        figs = [[0 for i in range(len(property))] for j in range(len(tensplt))]

    # set colors, nTensplt*nProperty
    clist = list(mcolors.TABLEAU_COLORS.keys())
    nclist = len(clist)
    if len(property) == 1 and len(tensplt) > 1:
        colors = [[clist[i%nclist]] for i in range(len(tensplt))]
    else:
        colors = [[clist[i%nclist] for i in range(len(property))] for j in range(len(tensplt))]

    # plot
    if uniform_scale_2D == False: # Non-uniform scale
        for ip, p in enumerate(property):
            kwargs['property'] = p
            for it, t in enumerate(tensplt):
                if same_fig_2D == False:
                    fig = t.plot_2D(**kwargs)
                    figs[it][ip] = fig
                else:
                    iplt = int(it*len(property) + ip)
                    kwargs['return_data'] = True
                    chi, r, title, uplt, utextplt, platt = t.plot_2D(**kwargs)
                    axes.flat[iplt] = _plot2D_single(
                        axes.flat[iplt], chi, r, np.array([0, 0, 1], dtype=float),
                        colors[it][ip], title, np.max(r), uplt, utextplt, platt
                    )
    else: # Uniform scale
        for ip, p in enumerate(property):
            chi_all = []
            r_all = []
            title_all = []
            uplt_all =[]
            utext_all = []
            platt_all = []
            kwargs['property'] = p
            kwargs['return_data'] = True
            # Get data first
            for it, t in enumerate(tensplt):
                chi, r, title, uplt, utextplt, platt = t.plot_2D(**kwargs)
                chi_all.append(chi)
                r_all.append(r)
                title_all.append(title)
                uplt_all.append(uplt)
                utext_all.append(utextplt)
                platt_all.append(platt)
            # plot
            rmax = np.max(r_all)
            for it, t in enumerate(tensplt):
                if same_fig_2D == False:
                    fig, ax = plt.subplots(nrows=1, ncols=1,
                                           subplot_kw={'projection' : 'polar'}, layout='tight')
                    ax = _plot2D_single(ax, chi_all[it], r_all[it], np.array([0, 0, 1], dtype=float),
                                        colors[it][ip], title_all[it], rmax,
                                        uplt_all[it], utextplt_all[it], platt_all[it])
                    figs[it][ip] = fig
                else:
                    iplt = int(it*len(property) + ip)
                    axes.flat[iplt] = _plot2D_single(
                        axes.flat[iplt], chi_all[it], r_all[it], np.array([0, 0, 1], dtype=float),
                        colors[it][ip], title_all[it], rmax, uplt_all[it], utext_all[it], platt_all[it]
                    )

    # dimensionality
    if same_fig_2D == False:
        if len(tensplt) == 1 and len(property) == 1:
            figs = figs[0][0]
        elif len(tensplt) == 1 and len(property) > 1:
            figs = figs[0]
        elif len(tensplt) > 1 and len(property) == 1:
            figs = [i[0] for i in figs]

    return figs


##############################################################################
#                                                                            #
#                             VIBRATIONAL PROPERTIES                         #
#                                                                            #
##############################################################################

#-----------------------------PHONON BAND AND DOS----------------------------#

def plot_phonon_bands(*bands, unit='cm-1', k_label=None, mode='single',
                      not_scaled=False, frequency_range=[], k_range=[],
                      band_label=None, band_color=None, band_linestyle=None,
                      band_linewidth=None, plot_freq0=True, freq0_color='tab:gray',
                      freq0_linestyle='-', freq0_linewidth=1.0, layout=None,
                      title=None, figsize=[6.4, 4.8], legend='lower right',
                      sharex=True, sharey=True, fontsize=14, **kwargs):
    """
    Plot phonon band structures.

    Args:
        \*bands (PhononBand|str): ``phonons.PhononBand`` object or band
            structure files of CRYSTAL. Note that lattice information is not
            available if file names are specified.
        unit (str): The unit of frequency. Can be 'cm-1' or 'THz'.
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
            of plot colors. 'None' for default values ('tab:blue' and other tab
            series).
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
            btmp = PhononBand.from_file(b)
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
    if plot_freq0 == True:
        freq0 = None
    else:
        freq0 = 0.

    # plot mode
    if mode.lower() == 'single':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = plot_compare_bands(
            [ax], [bandsplt[0].bands], [bandsplt[0].tick_pos],
            k_label, False, frequency_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, freq0, freq0_color, freq0_linestyle,
            freq0_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'multi':
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax = plot_overlap_bands(
            ax, [b.bands for b in bandsplt], [b.tick_pos for b in bandsplt],
            k_label, frequency_range, k_range, band_label, band_color,
            band_linestyle, band_linewidth, freq0, freq0_color, freq0_linestyle,
            freq0_linewidth, legend, **kwargs
        )
    elif mode.lower() == 'compare':
        if np.all(layout==None):
            layout = [int(np.ceil(len(bandsplt)/2)), 2]

        fig, ax = plt.subplots(layout[0], layout[1], figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        _ = plot_compare_bands(
            ax.flat, [b.bands for b in bandsplt], [b.tick_pos for b in bandsplt],
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
    *doss, unit='cm-1', overlap=False, prj=[], frequency_range=[], dos_range=[],
    dos_label=None, dos_color=None, dos_linestyle=None, dos_linewidth=None,
    plot_freq0=True, freq0_color='tab:gray', freq0_linestyle='-',
    freq0_linewidth=1.0, title=None, figsize=[6.4, 4.8], legend='lower right',
    sharex=True, sharey=False, fontsize=14, **kwargs):
    """
    Plot phonon density of states.

    Args:
        \*doss (PhononDOS|str): ``phonon.PhononDOS`` object or PHONODOSS files
            of CRYSTAL. Note that lattice information is not available if file
            names are specified.
        unit (str): The unit of frequency. Can be 'cm-1' or 'THz'.
        overlap (bool): Plotting multiple projections into the same figure.
            Useful only if a single entry of ``doss`` is plotted. Otherwise
            projections from the same entry will be overlapped into the same
            subplot.
        prj (list): Index of selected projections, consistent with the first
            dimension of the ``doss``, starting from 1. Effective for all the
            subplots.
        frequency_range (list): 1\*2 list of frequency range
        dos_range (list): 1\*2 list of DOS range
        dos_label (str|list): Plot legend. If only one string is given, apply
            it to all plots. 1\*nPrj plot legend otherwise. Effective for all
            the subplots.
        dos_color (str|list): Color of DOSS plots. If only one string is given,
            apply it to all plots. When ``overlap=True``, 1\*nPrj list of plot
            color. 'None' for default values ('tab:blue' and other tab series).
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
    dossplt = []
    for id, d in enumerate(doss):
        if isinstance(d, str):
            dtmp = PhononDOS.from_file(d)
        elif isinstance(d, PhononDOS):
            dtmp = copy.deepcopy(d)
        else:
            raise ValueError('Unknown input type for doss.')

        if unit != dtmp.unit:
            dtmp._set_unit(unit)
        dossplt.append(dtmp)

    # prj
    if len(prj) == 0:
        prj = [[int(i+1) for i in range(len(j.doss))] for j in dossplt]
    else:
        prj = [prj for j in dossplt]

    # frequency = 0
    if plot_freq0 == True:
        freq0 = None
    else:
        freq0 = 0.

    # plot
    ndoss = len(dossplt)
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
                ax.flat[i] = plot_doss(
                    ax.flat[i], dossplt[0].doss, dossplt[0].frequency, 'up',
                    [prj[0][i]], frequency_range, dos_range, commands[0][i],
                    commands[1][i], commands[2][i], commands[3][i], freq0,
                    freq0_color, freq0_linestyle, freq0_linewidth, legend,
                    False, **kwargs
                )
    else: # Projecton of the same system into the same panel
        fig, ax = plt.subplots(ndoss, 1, figsize=figsize,
                               sharex=sharex, sharey=sharey, layout='constrained')
        if ndoss == 1:
            ax = plot_doss(
                ax, dossplt[0].doss, dossplt[0].frequency, 'up', prj[0],
                frequency_range, dos_range, dos_label, dos_color, dos_linestyle,
                dos_linewidth, freq0, freq0_color, freq0_linestyle, freq0_linewidth,
                legend, False, **kwargs
            )
        else:
            for i in range(ndoss):
                ax.flat[i] = plot_doss(
                    ax.flat[i], dossplt[i].doss, dossplt[i].frequency, 'up', prj[i],
                    frequency_range, dos_range, dos_label, dos_color, dos_linestyle,
                    dos_linewidth, freq0, freq0_color, freq0_linestyle,
                    freq0_linewidth, legend, False, **kwargs
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
    *data, unit='cm-1', k_label=[], dos_overlap=True, dos_prj=[], frequency_range=[],
    k_range=[], dos_range=[], band_width=2, band_label=None, band_color=None,
    band_linestyle=None, band_linewidth=None, dos_label=None, dos_color=None,
    dos_linestyle=None, dos_linewidth=None, plot_freq0=True, freq0_color='tab:green',
    freq0_linestyle='-', freq0_linewidth=1.0, title=None, figsize=[6.4, 4.8],
    legend='lower right', fontsize=14, **kwargs):
    """
    Plot phonon band structure + dos for a **single** system, i.e., the
    ``bands`` and ``doss`` variables are not extendable.

    Input arguments not in the list are consistent with ``plot_phonon_doss``
    and ``plot_phonon_bands``.

    Args:
        \*data: Either 1 or 2 entries. For one enetry, it is fort.25 containing
            both band and DOS, or ``PhononBandDOS`` object. For 2 entries, the
            first entry is ``bands`` of ``plot_phonon_bands`` and the second is
            ``doss`` of ``plot_phonon_doss``.
        dos_overlap (bool): ``overlap`` of ``plot_phonon_doss``. The user can
            either plot projections into the same subplot or into separate
            subplots.
        dos_prj (list): ``prj`` of ``plot_phonon_doss``.
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
            bands = PhononBand.from_file(data[0])
            doss = PhononDOS.from_file(data[0])
        elif isinstance(data[0], PhononBandDOS):
            bands = copy.deepcopy(data[0].band)
            doss = copy.deepcopy(data[0].dos)
        else:
            raise ValueError('Unknown input data type for the 1st entry.')
    elif len(data) == 2:
        if isinstance(data[0], str):
            bands = PhononBand.from_file(data[0])
        elif isinstance(data[0], PhononBand):
            bands = copy.deepcopy(data[0])
        else:
            raise ValueError('Unknown input data type for the 1st entry.')

        if isinstance(data[1], str):
            doss = PhononDOS.from_file(data[1])
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
    if plot_freq0 == True:
        freq0 = None
    else:
        freq0 = 0.

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


#------------------------------------SPECTRA----------------------------------#

def plot_cry_irspec(irspec, x_unit='cm-1', y_mode='LG', figsize=None, linestyle='-',
                    linewidth=1.5, color='tab:blue', freq_range=None, int_range=None,
                    label=None):
    """Generates the IR spectra for the IRSPEC.DAT file produced by an IRSPEC calculation

    Args:
        irspec (External_unit object): Object generated by the read_cry_irspec function necessary for the plot
        x_unit (str, optional): Unit measure of the x axes. Avalilable: 'cm-1' and 'nm'. Defaults to 'cm-1'.
        y_mode (str, optional): Peak broadening modality in absorbance and reflectance. 
                                Available: 'LG'(Lorentzian-Gaussian broadening), 'V' (Voight broadening), 'RS' (Rayleigh spherical particles), 'RE' (Rayleigh with elipsoid particles), 'REFL' (Reflectance)
                                Defaults to 'LG'.
        figsize (tuple, optional): Image dimensions correspondig to matplotlib figsize. Defaults to None.
        linestyle (str/list[str], optional): linestyle corresponding to the matplotlib one it can be a list for a multiplot. Defaults to '-'.
        linewidth (float/list[float], optional): linewidth corresponding to the matplotlib one it can be a list for a multiplot. Defaults to 1.5.
        color (str/list[str], optional): Color of the spectra it can accept all matplotlib colors it can be a list for multiplots. Defaults to 'tab:blue'.
        freq_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given frequency window. Defaults to None.
        int_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given intensity window. Defaults to None.
        label (list[str], optional): List of labels for the legend of a multiplot. Defaults to None.
        save_to_file (str, optional): Filename of the spectra to be saved. Defaults to None.
        dpi (int, optional): Resolution of the saved file. Defaults to 300.
        transparency (bool, optional): Enables the transparency of the saved file background. Defaults to False.

    :raise ValueError: The function raises an error when the object to be plotted does not have the required y_mode  
    """

    import sys
    import warnings

    import matplotlib.pyplot as plt
    import numpy as np

    modes = ['single', 'multi']
    accepted_y = ['LG', 'V', 'RS', 'RE', 'REFL']

    if isinstance(irspec, list):
        mode = modes[1]

        if not isinstance(linestyle, list):
            style = linestyle
            linestyle = []
            for i in enumerate(irspec):
                linestyle.append(style)

        if not isinstance(linewidth, list):
            width = linewidth
            linewidth = []
            for i in enumerate(irspec):
                linewidth.append(width)

        if not isinstance(color, list):
            color = ['dimgrey', 'blue', 'indigo', 'slateblue',
                     'thistle', 'purple', 'orchid', 'crimson']

        for file in irspec:
            if (file.calculation == 'molecule') and (y_mode != accepted_y[0]):
                raise ValueError(
                    'This spectra does not contain the y_mode requested: available y_mode'+accepted_y[0])

    else:
        mode = modes[0]

        if (irspec.calculation == 'molecule') and (y_mode != accepted_y[0]):
            raise ValueError(
                'This spectra does not contain the y_mode requested: available y_mode'+accepted_y[0])

    if figsize is not None:
        fig, ax = plt.subplots(figsize=figsize)

    if mode == modes[0]:

        # selection of the x axis unit
        if x_unit == 'cm-1':
            x = irspec.irspec[:, 0]

        elif x_unit == 'nm':
            x = irspec.irspec[:, 1]

        # selection of the intensities mode
        if y_mode == accepted_y[0]:
            y = irspec.irspec[:, 2]

        elif y_mode == accepted_y[1]:
            y = irspec.irspec[:, 5]

        elif y_mode == accepted_y[2]:
            y = irspec.irspec[:, 6]

        elif y_mode == accepted_y[3]:
            y = irspec.irspec[:, 7]

        elif y_mode == accepted_y[4]:
            y = irspec.irspec[:, 8]

        print(x, y)

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)-1
        ymax = max(y)+10

        plt.plot(x, y, linestyle=linestyle, linewidth=linewidth, color=color)

    if mode == modes[1]:

        xmin = []
        xmax = []
        ymin = []
        ymax = []

        for index, file in enumerate(irspec):
            # selection of the x axis unit
            if x_unit == 'cm-1':
                x = file.irspec[:, 0]

            elif x_unit == 'nm':
                x = file.irspec[:, 1]

            # selection of the intensities mode
            if y_mode == accepted_y[0]:
                y = file.irspec[:, 2]

            elif y_mode == accepted_y[1]:
                y = file.irspec[:, 5]

            elif y_mode == accepted_y[2]:
                y = file.irspec[:, 6]

            elif y_mode == accepted_y[3]:
                y = file.irspec[:, 7]

            elif y_mode == accepted_y[4]:
                y = file.irspec[:, 8]

            xmin.append(min(x))
            xmax.append(max(x))
            ymin.append(min(y)-1)
            ymax.append(max(y)+10)

            if label is not None:
                ax.plot(x, y, linestyle=linestyle[index], linewidth=linewidth[index],
                               color=color[index], label=label[index])
                ax.legend()
            else:
                ax.plot(
                    x, y, linestyle=linestyle[index], linewidth=linewidth[index], color=color[index])

        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)

    if freq_range is not None:
        xmin = freq_range[0]
        xmax = freq_range[1]

    if int_range is not None:
        ymin = int_range[0]
        ymax = int_range[1]

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    if x_unit == 'cm-1':
        plt.xlabel('Wavenumber (cm$^{-1}$)')
    elif x_unit == 'nm':
        plt.xlabel('Wavelength (nm)')

    if y_mode != accepted_y[4]:
        plt.ylabel('Absorbance (A.U.)')
    else:
        plt.ylabel('Reflectance (A.U.)')

    return fig, ax


def plot_cry_ramspec(ramspec,  y_mode='total', figsize=None, linestyle='-',
                     linewidth=1.5, color='tab:blue', freq_range=None, int_range=None,
                     label=None):
    """Generates the RAMAN spectra for the RAMSPEC.DAT file produced by an RAMSPEC calculation

    Args:
        irspec (External_unit object): Object generated by the read_cry_ramspec function necessary for the plot
        y_mode (str, optional): Polarization of the spectra for the simulated compound
                                Available: 'total', 'parallel', 'perpendicular' (for powders), 'xx', 'xy', 'xz', 'yy', 'yz', 'zz' (for single crystals)
                                Defaults to 'LG'.
        figsize (tuple, optional): Image dimensions correspondig to matplotlib figsize. Defaults to None.
        linestyle (str/list[str], optional): linestyle corresponding to the matplotlib one it can be a list for a multiplot. Defaults to '-'.
        linewidth (float/list[float], optional): linewidth corresponding to the matplotlib one it can be a list for a multiplot. Defaults to 1.5.
        color (str/list[str], optional): Color of the spectra it can accept all matplotlib colors it can be a list for multiplots. Defaults to 'tab:blue'.
        freq_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given frequency window. Defaults to None.
        int_range (list, optional): Two element list [min, max], that allows to visualize the spectra in a given intensity window. Defaults to None.
        label (list[str], optional): List of labels for the legend of a multiplot. Defaults to None.
        save_to_file (str, optional): Filename of the spectra to be saved. Defaults to None.
        dpi (int, optional): Resolution of the saved file. Defaults to 300.
        transparency (bool, optional): Enables the transparency of the saved file background. Defaults to False.
    """

    import sys
    import warnings

    import matplotlib.pyplot as plt
    import numpy as np

    modes = ['single', 'multi']
    accepted_y = ['total', 'parallel', 'perpendicular',
                  'xx', 'xy', 'xz', 'yy', 'yz', 'zz']

    if isinstance(ramspec, list):
        mode = modes[1]
        if not isinstance(linestyle, list):
            style = linestyle
            linestyle = []
            for i in enumerate(ramspec):
                linestyle.append(style)

        if not isinstance(linewidth, list):
            width = linewidth
            linewidth = []
            for i in enumerate(ramspec):
                linewidth.append(width)

        if not isinstance(color, list):
            color = ['dimgrey', 'blue', 'indigo', 'slateblue',
                     'thistle', 'purple', 'orchid', 'crimson']

    else:
        mode = modes[0]

    if figsize is not None:
        fig, ax = plt.subplots(figsize=figsize)

    if mode == modes[0]:

        x = ramspec.ramspec[:, 0]

        # selection of the intensities mode
        if y_mode == accepted_y[0]:
            y = ramspec.ramspec[:, 1]

        elif y_mode == accepted_y[1]:
            y = ramspec.ramspec[:, 2]

        elif y_mode == accepted_y[2]:
            y = ramspec.ramspec[:, 3]

        elif y_mode == accepted_y[3]:
            y = ramspec.ramspec[:, 4]

        elif y_mode == accepted_y[4]:
            y = ramspec.ramspec[:, 5]

        elif y_mode == accepted_y[5]:
            y = ramspec.ramspec[:, 6]

        elif y_mode == accepted_y[6]:
            y = ramspec.ramspec[:, 7]

        elif y_mode == accepted_y[7]:
            y = ramspec.ramspec[:, 8]

        elif y_mode == accepted_y[8]:
            y = ramspec.ramspec[:, 9]

        xmin = min(x)
        xmax = max(x)
        ymin = min(y)-1
        ymax = max(y)+10

        fig = plt.plot(x, y, linestyle=linestyle,
                       linewidth=linewidth, color=color)

    if mode == modes[1]:
        xmin = []
        xmax = []
        ymin = []
        ymax = []

        for index, file in enumerate(ramspec):
            x = file.ramspec[:, 0]

            # selection of the intensities mode
            if y_mode == accepted_y[0]:
                y = file.ramspec[:, 1]

            elif y_mode == accepted_y[1]:
                y = file.ramspec[:, 2]

            elif y_mode == accepted_y[2]:
                y = file.ramspec[:, 3]

            elif y_mode == accepted_y[3]:
                y = file.ramspec[:, 4]

            elif y_mode == accepted_y[4]:
                y = file.ramspec[:, 5]

            elif y_mode == accepted_y[5]:
                y = file.ramspec[:, 6]

            elif y_mode == accepted_y[6]:
                y = file.ramspec[:, 7]

            elif y_mode == accepted_y[7]:
                y = file.ramspec[:, 8]

            elif y_mode == accepted_y[8]:
                y = file.ramspec[:, 9]

            xmin.append(min(x))
            xmax.append(max(x))
            ymin.append(min(y)-1)
            ymax.append(max(y)+10)

            if label is not None:
                ax.plot(x, y, linestyle=linestyle[index], linewidth=linewidth[index],
                               color=color[index], label=label[index])
                plt.legend()
            else:
                ax.plot(
                    x, y, linestyle=linestyle[index], linewidth=linewidth[index], color=color[index])

        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)

    if freq_range is not None:
        xmin = freq_range[0]
        xmax = freq_range[1]

    if int_range is not None:
        ymin = int_range[0]
        ymax = int_range[1]

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xlabel('Wavenumber (cm$^{-1}$)')

    if y_mode != accepted_y[4]:
        plt.ylabel('Absorbance (A.U.)')
    else:
        plt.ylabel('Reflectance (A.U.)')

    return fig, ax


#-----------------------------------ANHARMONIC--------------------------------#

def plot_cry_spec(transitions, typeS, components=False, bwidth=5, stdev=3, eta=0.5,
                  fmin=None, fmax=None, ylim=None, savefig=False, dpi=300,
                  filetype='png', exp_spec=None, sep=";", show=True,
                  export_csv=False, label=None, xlabel='Wavenumber [cm$^{-1}$]',
                  ylabel='Intensity [arb. u.]', linewidth=2.0, padd=100,
                  fontsize=12, style=None, compstyle=None, nopadding=False,
                  figsize=(16, 6)):
    """
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

    import math
    import time
    from copy import deepcopy

    import matplotlib.pyplot as plt
    import numpy as np
    from numpy import genfromtxt

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

    import time

    import matplotlib.pyplot as plt

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

def plot_cry_ela(choose, ndeg, *args, dpi=200, filetype=".png",
                 transparency=False):
    """
    Deprecated. Use ``plot_elastics_3D``.
    """
    import warnings

    warnings.warn("You are calling a deprecated function. Use 'plot_elastics_3D' instead.",
                  stacklevel=2)
    args = [i for i in args]
    figs = plot_elastics_3D(
        choose, *args, uniform_scale=True, add_title=False, nphi=ndeg,
        ntheta=ndeg, nchi=ndeg
    )

    if not isinstance(figs, list):
        figs = [figs]
    for nfig, fig in enumerate(figs):
        fig.savefig(choose + '{:d}'.format(nfig) + filetype,
                    dpi=dpi, transparent=transparency)
    return


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
    return contour_obj.plot()


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

    fig = plot_transport_tensor(seebeck_obj, sigma_obj, option='power factor', direction=dir)
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

    fig = plot_transport_tensor(seebeck_obj, sigma_obj, option='power factor',
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
    fig = plot_transport_tensor(seebeck_obj, sigma_obj, option='power factor', direction=dir)
    return fig


