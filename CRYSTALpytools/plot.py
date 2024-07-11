#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to visualize CRYSTAL outputs.
"""
import numpy as np

##############################################################################
#                                                                            #
#                       ELECTRONIC STRUCTURE AND PHONONS                     #
#                                                                            #
##############################################################################

#-------------------------------ECHG charge density----------------------------#


def plot_dens_ECHG(obj_echg, unit='Angstrom', levels=150, xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None):
    """
    Plots the 2D ECHG density map from a fort.25 file.

    Args:
        obj_echg (ChargeDensity): Charge/spin density object.
        unit (str): The energy unit for **plotting**. 'Angstrom' for :math:`e.\AA^{-3} or 'a.u.' for :math:`e.Bohr^{-3}`.
        levels (int | array-like): Number and positions of the contour lines/regions.
        xticks (int): Number of ticks in the x direction.
        yticks (int): Number of ticks in the y direction.
        cmap_max(float): Maximum value used for the colormap.
        cmap_min(float): Minimun value used for the colormap.

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import copy

    from CRYSTALpytools.base.plotbase import plot_2Dscalar

    obj = copy.deepcopy(obj_echg)
    obj._set_unit(unit)
    if unit.lower() == 'angstrom':
        cbarlabel = 'Charge Density ($|e|.\AA^{-3}$)'
    else:
        cbarlabel = 'Charge Density ($|e|.Bohr^{-3}$)'

    fig, ax = plot_2Dscalar(obj.chgmap, obj.gridv, levels, xticks, yticks,
                        cmap_max, cmap_min, cbarlabel)

    plt.show()
    return fig, ax

#--------------------------------ECHG spin density-----------------------------#

def plot_spin_ECHG(obj_echg, unit='Angstrom', levels=150, xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None):
    """
    Plots the 2D spin density map from a ECHG output file (fort.25). For charge
    density map please refer to ``plot_dens_ECHG``.

    Args:
        obj_echg (ChargeDensity): Charge/spin density object.
        unit (str): The energy unit for **plotting**. 'Angstrom' for :math:`e.\AA^{-3} or 'a.u.' for :math:`e.Bohr^{-3}`.
        levels (int | array-like): *Optional* Determines the number and positions of the contour lines/regions. Default is 150.
        xticks (int): *Optional* Number of ticks in the x direction. Default is 5.
        yticks (int): *Optional* Number of ticks in the y direction. Default is 5.
        cmap_max(float): *Optional*, Maximum value used for the colormap. Default is None.
        cmap_min(float): *Optional* Minimun value used for the colormap. Default is None.

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import copy

    from CRYSTALpytools.base.plotbase import plot_2Dscalar

    obj = copy.deepcopy(obj_echg)
    obj._set_unit(unit)
    if unit.lower() == 'angstrom':
        cbarlabel = 'Spin Density ($|e|.\AA^{-3}$)'
    else:
        cbarlabel = 'Spin Density ($|e|.Bohr^{-3}$)'

    fig, ax = plot_2Dscalar(obj.spinmap, obj.gridv, levels, xticks, yticks,
                        cmap_max, cmap_min, cbarlabel)

    plt.show()
    return fig, ax

#----------------------------------SPIN CURRENTS------------------------------#

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

def plot_phonon_band(bands, unit='cm-1', k_labels=None, mode='single',
                     not_scaled=False, energy_range=None, k_range=None,
                     color='blue', labels=None, linestl='-', linewidth=1,
                     line_freq0=None, title=None, figsize=None,
                     scheme=None, sharex=True, sharey=True, fontsize=12):
    """
    A wrapper of plot_cry_bands for phonon band structure.

    Args:
        bands (BandsBASE|list): Bands object generated by `CRYSTALpytools.crystal_io.Crystal_output.read_pband` or
            a list of BandsBASE objects.
        unit (str): The unit of frequency. Can be 'cm-1' or 'THz'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be, for example, 'Gamma'.
        mode (str): The plotting mode. Possible values are 'single', 'multi', and 'compare'.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
        energy_range (array): A 2x1 array specifying the energy range.
        k_range (array): A 2x1 array specifying the k-range.
        color (str|list): Color of plot lines. Should be consistent with bands.
        labels (str|list): Plot legend. Should be consistent with bands.
        linestl (str|list): Linestyle string. Should be consistent with bands.
        linewidth (float): The width of the plot lines.
        line_freq0 (str): The color of the frequency=0 line.
        title (str): The title of the plot.
        figsize (list): The figure size specified as [width, height].
        scheme (list|tuple): The layout of subplots.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        save_to_file (str): The file name to save the plot.
        dpi (int): Dots per inch resolution of the saved file.
        fontsize (int): Fontsize of the axis labels.

    Returns:
        None

    :raise ValueError: If the specified unit is unknown.
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALpytools.base.plotbase import plot_cry_bands
    from CRYSTALpytools.units import cm_to_thz, thz_to_cm

    if re.match(r'^cm\-1$', unit, re.IGNORECASE):
        unit = 'cm-1'
        is_thz = False
    elif re.match(r'^THz$', unit, re.IGNORECASE):
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    if np.all(line_freq0==None):
        line_freq0 = (1., 0., 0., 0.)  # Transparent

    for b in bands:
        if unit != b.unit:
            if unit == 'cm-1':
                b.bands[:, :, :] = thz_to_cm(b.bands[:, :, :])
            else:
                b.bands[:, :, :] = cm_to_thz(b.bands[:, :, :])
            b.unit = unit
    if len(bands) == 1:
        bands = bands[0]

    fig, ax = plot_cry_bands(bands, k_labels=k_labels, energy_range=energy_range, title=title,
                         not_scaled=not_scaled, mode=mode, linestl=linestl, linewidth=linewidth,
                         color=color, fermi=line_freq0, k_range=k_range, labels=labels,
                         figsize=figsize, scheme=scheme, sharex=sharex, sharey=sharey, fermialpha=1, fermiwidth=0)
    if is_thz == True:
        fig.supylabel('Frequency (THz)', fontsize=fontsize)
    else:
        fig.supylabel('Frequency (cm$^{-1}$)', fontsize=fontsize)

    return fig, ax


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
        mode (str): The plotting mode, including 'single', 'multi' and
            'compare'.
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
            1\*nSystem or nSystem\*2 (spin) linewidth string. If spin>1 and
            spin dimension is not included, spin states are in the same width.
            'None' for default values (1.0).
        fermi_level (float|list|None): Fermi energy in the same unit as input
            band energy. By default the band is aligned to 0. Can be used to
            offset the band. None for not plotting Fermi. For 'compare' mode,
            different offsets can be used.
        fermi_color (str): Color of the Fermi level.
        fermi_linestyle (str): Line style of Fermi level.
        fermi_linewidth(float): Width of the Fermi level.
        layout (list|tuple): For 'compare' mode, the layout of subplots. The
            default is 2 cols per row.
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
        ax (Axes): Matplotlib axes object

    :raise ValueError: If the specified unit is unknown.
    """
    import matplotlib.pyplot as plt
    from CRYSTALpytools.electronics import ElectronBand
    from CRYSTALpytools.base.plotbase import plot_overlap_bands, plot_compare_bands
    from CRYSTALpytools.units import H_to_eV, eV_to_H

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
            btmp = b
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

    return fig, fig.axes


#-------------------------------DENSITY OF STATES-----------------------------#


def plot_electron_doss(*doss, unit='eV', beta='up', overlap=False, prj=[],
                       energy_range=[], dos_range=[], dos_label=None,
                       dos_color=None, dos_linestyle=None, dos_linewidth=None,
                       fermi_level=0., fermi_color='tab:green', fermi_linestyle='-',
                       fermi_linewidth=1.0, title=None, figsize=[6.4, 4.8],
                       legend='lower right', sharex=True, sharey=True,
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
            apply it to all plots. Also 1\*2 color list for spin. When
            ``overlap=True``, 1\*nPrj or nPrj\*2 (spin) plot color. If spin>1
            and spin dimension is not included, spin states are in the same
            color. 'None' for default values ('tab:blue' and other tab series).
            Effective for all the subplots.
        dos_linestyle (str|list): Linestyle of DOSS plot. If only one string is
            given, apply it to all plots. Also 1\*2 color list for spin. When
            ``overlap=True``, 1\*nPrj or nPrj\*2 (spin) line styles. If spin>1
            and spin dimension is not included, spin states are in the same
            linestyle. 'None' for default values ('-'). Effective for all the
            subplots.
        dos_linewidth (str|list): Linewidth of DOSS plot. If only one number is
            given, apply it to all plots. Also 1\*2 color list for spin. When
            ``overlap=True``, 1\*nPrj or nPrj\*2 (spin) line widthes. If spin>1
            and spin dimension is not included, spin states are in the same
            linestyle. 'None' for default values (1.0). Effective for all the
            subplots.
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
        ax (Axes): Matplotlib axes object
    """
    import matplotlib.pyplot as plt
    from CRYSTALpytools.base.plotbase import plot_doss, _plot_label_preprocess
    from CRYSTALpytools.electronics import ElectronDOS
    from CRYSTALpytools.units import H_to_eV, eV_to_H

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
            dtmp = d
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
        if nprj == 1:
            ax = plot_doss(
                ax, dossplt[0].doss, dossplt[0].energy, beta, [prj[0][0]],
                energy_range, dos_range, dos_label, dos_color, dos_linestyle,
                dos_linewidth, fermi_level, fermi_color, fermi_linestyle,
                fermi_linewidth, legend, False, **kwargs
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
                    ax.flat[i], dossplt[0].doss, dossplt[0].energy, beta, [prj[0][i]],
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

    return fig, fig.axes


def plot_phonon_dos(doss, unit='cm-1', overlap=False, prj=None,
                    freq_range=None, dos_range=None, color='blue',
                    labels=None, linestl=None, linewidth=1, line_freq0=None,
                    title=None, figsize=None):
    """
    A wrapper of plot_cry_doss for electron density of states.

    Args:
        doss (DOSBASE): DOS obect generated by code:`CRYSTALpytools.crystal_io.Crystal_output.read_pdos`.
            Or a list of DOSBASE objects.
        unit (str): 'cm-1' or 'THz'
        overlap (bool): Plotting multiple lines into the same figure
        prj (list): Index of selected projection. Consistent with the
            index of the 2nd dimension of :code:`doss.doss`
        freq_range (list[float]): 2*1 list of frequency range
        dos_range (list[float]): 2*1 list of DOS range
        color (str | list[str]): Color of plot lines. *Should be
            consistent with number of projections.*
        labels (str | list[str]): Plot legend. *Should be consistent with
            number of projections.*
        linestl (str | list[str]): linestyle string. *Should be consistent
            with number of projections.*
        linewidth (float)
        line_freq0 (str): Color of frequency = 0 line.
        title (str)
        figsize (list[float])

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALpytools.base.plotbase import plot_cry_doss
    from CRYSTALpytools.units import cm_to_thz, thz_to_cm

    if re.match(r'^cm\-1$', unit, re.IGNORECASE):
        unit = 'cm-1'
        is_thz = False
    elif re.match(r'^thz$', unit, re.IGNORECASE):
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    for d in doss:
        if unit != d.unit:
            if unit == 'cm-1':
                d.doss[:, 0, :] = thz_to_cm(d.doss[:, 0, :])
                d.doss[:, 1:, :] = cm_to_thz(d.doss[:, 1:, :])
            else:
                d.doss[:, 0, :] = cm_to_thz(d.doss[:, 0, :])
                d.doss[:, 1:, :] = thz_to_cm(d.doss[:, 1:, :])
            d.unit = unit

    if np.all(line_freq0==None):
        line_freq0 = (1., 0., 0., 0.)  # Transparent
    if len(doss) == 1:
        doss = doss[0]

    fig, ax = plot_cry_doss(doss, color=color, fermi=line_freq0, overlap=overlap,
                        labels=labels, figsize=figsize, linestl=linestl,
                        linewidth=linewidth, title=title, beta='up',
                        energy_range=freq_range, dos_range=dos_range, prj=prj)

    if is_thz == True:
        fig.supylabel('DOS (states/THz)')
        fig.supxlabel('Frequency (THz)')
    else:
        fig.supylabel('DOS (states/cm$^{-1}$)')
        fig.supxlabel('Frequency (cm$^{-1}$)')

    return fig, ax


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
        ax (Axes): Matplotlib axes object

    :raise ValueError: If the unit parameter is unknown.
    """
    from CRYSTALpytools.electronics import ElectronBand, ElectronDOS, ElectronBandDOS
    from CRYSTALpytools.base.plotbase import plot_banddos
    from CRYSTALpytools.units import H_to_eV, eV_to_H
    import warnings

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
            bands = data[0].band
            doss = data[0].dos
        else:
            raise ValueError('Unknown input data type for the 1st entry.')
    elif len(data) == 2:
        if isinstance(data[0], str):
            bands = ElectronBand.from_file(data[0])
        elif isinstance(data[0], ElectronBand):
            bands = data[0]
        else:
            raise ValueError('Unknown input data type for the 1st entry.')

        if isinstance(data[1], str):
            doss = ElectronDOS.from_file(data[1])
        elif isinstance(data[1], ElectronDOS):
            doss = data[1]
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

    fig, ax = plot_banddos(
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

    return fig, fig.axes


def plot_phonon_banddos(bands, doss, unit='cm-1', k_labels=None, dos_prj=None,
                        freq_range=None, dos_max_range=None, color_band='blue',
                        color_dos='blue', labels=None, linestl_band='-',
                        linestl_dos=None, linewidth=1, freq0_line=None,
                        title=None, figsize=None):
    """
    A wrapper of plot_cry_es for phonon band structure + dos. Only one pair is permitted.

    Args:
        bands (BandsBASE|list): Bands object generated by CRYSTALpytools.crystal_io.Properties_output.read_bands
            or a list of BandsBASE objects.
        doss (DOSBASE): DOS object generated by CRYSTALpytools.crystal_io.Properties_output.read_doss
            or a list of DOSBASE objects.
        unit (str): Unit of frequency. Valid options are 'cm-1' or 'THz'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be
            represented as strings, for example, 'Gamma'.
        dos_prj (list): Index of selected projection. Consistent with the index of the 2nd dimension
            of doss.doss.
        freq_range (list): A list of two values representing the frequency range to be plotted.
        dos_max_range (float): Maximum DOS range for the y-axis.
        color_band (str): Color of the phonon bands in the plot.
        color_dos (str): Color of the density of states (DOS) in the plot.
        labels (list): A list of labels for the plot legend.
        linestl_band (str): Linestyle of the phonon bands.
        linestl_dos (str): Linestyle of the density of states (DOS).
        linewidth (float): Width of the lines in the plot.
        freq0_line (str): Color of the frequency=0 line.
        title (str): Title of the plot.
        figsize (list[float]): Size of the figure in inches (width, height).

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object

    :raise ValueError: If the unit parameter is unknown.
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALpytools.base.plotbase import plot_cry_es
    from CRYSTALpytools.units import cm_to_thz, thz_to_cm

    if re.match(r'^cm\-1$', unit, re.IGNORECASE):
        unit = 'cm-1'
        is_thz = False
    elif re.match(r'^thz$', unit, re.IGNORECASE):
        unit = 'THz'
        is_thz = True
    else:
        raise ValueError('Unknown unit.')

    if unit != doss.unit:
        if unit == 'cm-1':
            doss.doss[:, 0, :] = thz_to_cm(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = cm_to_thz(doss.doss[:, 1:, :])
        else:
            doss.doss[:, 0, :] = cm_to_thz(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = thz_to_cm(doss.doss[:, 1:, :])
        doss.unit = unit
    if unit != bands.unit:
        if unit == 'cm-1':
            bands.bands[:, :, :] = thz_to_cm(bands.bands[:, :, :])
        else:
            bands.bands[:, :, :] = cm_to_thz(bands.bands[:, :, :])
        bands.unit = unit

    if np.all(line_freq0==None):
        line_freq0 = (1., 0., 0., 0.)  # Transparent

    fig, ax = plot_cry_es(bands=bands, doss=doss, k_labels=k_labels, color_bd=color_band,
                      color_doss=color_dos, fermi=line_freq0, energy_range=energy_range,
                      linestl_bd=linestl_band, linestl_doss=linestl_dos,
                      linewidth=linewidth, prj=dos_prj, figsize=figsize, labels=labels,
                      dos_max_range=dos_max_range, title=title, dos_beta='up')
    if is_thz == True:
        fig.supylabel('Frequency (THz)')
    else:
        fig.supylabel('Frequency (cm$^{-1}$)')

    return fig, ax


##############################################################################
#                                                                            #
#                                     QTAIM                                  #
#                                                                            #
##############################################################################

#----------------------------------CONTOUR PLOT-------------------------------#


def plot_cry_contour(contour_obj):
    """
    Plot a contour plot.

    Args:
        contour_obj (object): Contour object representing the contour plot.
        save_to_file (bool, optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Plots a contour plot based on the data in the contour object.
        - Retrieves the data from the contour object and converts it to a 2D list.
        - Sets the figure size based on x_graph_param and y_graph_param attributes of the contour object.
        - Sets the x-axis and y-axis labels.
        - Creates a meshgrid using the x_points and y_points attributes of the contour object.
        - Defines contour levels, colors, linestyles, and fmt.
        - Plots the contour plot.
        - Saves the plot to a file named 'figure_TIPO_YYYY-MM-DD_HHMMSS.jpg' in the current directory.
        - If save_to_file is True, saves the plot to a file specified by save_to_file parameter.

    """
    import os
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    df = contour_obj.df
    n_punti_x = contour_obj.npx

    for i in range(0, 8):
        df[i] = df[i].astype(float)

    flat_list = [item for sublist in df.values for item in sublist]

    cleaned_list = [x for x in flat_list if ~np.isnan(x)]

    l = [cleaned_list[x:x+n_punti_x]
         for x in range(0, len(cleaned_list), n_punti_x)]

    c = contour_obj.x_graph_param
    d = contour_obj.y_graph_param

    plt.rcParams["figure.figsize"] = [c, d]

    plt.xlabel(r'$\AA$', fontsize=18)
    plt.ylabel(r'$\AA$', fontsize=18)

    X, Y = np.meshgrid(contour_obj.x_points, contour_obj.y_points)

    levels = contour_obj.levels
    colors = contour_obj.colors
    linestyles = contour_obj.linestyles
    fmt = contour_obj.fmt

    # Change here to have or not the isovalues on the plot
    iso = True
    # iso = False

    if (iso == True):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)
        plt.clabel(L, inline=1, fontsize=7, fmt=fmt)
    elif (iso == False):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)

    path = os.path.join('./'+'figure_' + contour_obj.tipo +
                        '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.savefig(path, bbox_inches='tight', dpi=600)
    print('\nThe image has been saved in the current directory')

    # if save_to_file != False:
    #     save_plot(save_to_file)

    plt.show()


def plot_cry_contour_differences(contour_obj, contour_obj_ref):
    """
    Plot the differences between two contour plots.

    Args:
        contour_obj (object): Contour object representing the original contour plot.
        contour_obj_ref (object): Contour object representing the reference contour plot.
        save_to_file (bool, optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Plots the differences between two contour plots.
        - Requires the contour objects to have a tipo attribute with values 'SURFLAPP', 'SURFLAPM', 'SURFRHOO', or 'SURFELFB'.
        - Calculates the difference between the dataframes of the two contour objects.
        - Sets the figure size based on x_graph_param and y_graph_param attributes of the contour object.
        - Sets the x-axis and y-axis labels.
        - Creates a meshgrid using the x_points and y_points attributes of the contour object.
        - Defines contour levels, colors, and linestyles.
        - Plots the contour differences.
        - Saves the plot to a file named 'figure_diff_TIPO_YYYY-MM-DD_HHMMSS.jpg' in the current directory.
        - If save_to_file is True, saves the plot to a file specified by save_to_file parameter.

    """
    import os
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    if (contour_obj.tipo == 'SURFLAPP') or (contour_obj.tipo == 'SURFLAPM') or (contour_obj.tipo == 'SURFRHOO') or (contour_obj.tipo == 'SURFELFB'):
        pass
    else:
        sys.exit(
            'Difference option only allowed for SURFLAPP, SURFLAPM, SURFRHOO and SURFELFB file')

    n_punti_x = contour_obj.npx

    df = contour_obj.df
    for i in range(0, 8):
        df[i] = df[i].astype(float)

    df_ref = contour_obj_ref.df
    for i in range(0, 8):
        df_ref[i] = df_ref[i].astype(float)

    df_diff = df - df_ref

    flat_list = [item for sublist in df_diff.values for item in sublist]

    cleaned_list = [x for x in flat_list if ~np.isnan(x)]

    l = [cleaned_list[x:x+n_punti_x]
         for x in range(0, len(cleaned_list), n_punti_x)]

    c = contour_obj.x_graph_param
    d = contour_obj.y_graph_param

    plt.rcParams["figure.figsize"] = [c, d]

    plt.xlabel(r'$\AA$', fontsize=18)
    plt.ylabel(r'$\AA$', fontsize=18)

    X, Y = np.meshgrid(contour_obj.x_points, contour_obj.y_points)

    ctr1dif = np.array([-8, -4, -2, -0.8, -0.4, -0.2, -0.08, -0.04, -0.02, -0.008, -0.004, -0.002, -0.0008, -0.0004, -0.0002, 0,
                       0.0002, 0.0004, 0.0008, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8, 2, 4, 8])
    colors1dif = ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'k', 'r', 'r', 'r', 'r', 'r', 'r', 'r',
                  'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r']
    ls1dif = ['--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', 'dotted', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

    levels = ctr1dif
    colors = colors1dif
    linestyles = ls1dif
    fmt = '%1.4f'

    # Change here to have or not the isovalues on the plot
    iso = True
    # iso = False

    if (iso == True):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)
        plt.clabel(L, inline=1, fontsize=7, fmt=fmt)
    elif (iso == False):
        L = plt.contour(X, Y, l, levels=levels, colors=colors, linestyles=linestyles, linewidths=0.7,
                        alpha=1)

    path = os.path.join('./'+'figure_diff_' + contour_obj.tipo +
                        '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.savefig(path, bbox_inches='tight', dpi=600)
    print('\nThe image has been saved in the current directory')

    # if save_to_file != False:
    #     save_plot(save_to_file)

    plt.show()

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

#-------------------------------------SEEBACK---------------------------------#


def plot_cry_seebeck_potential(seebeck_obj):
    """
    Plot the Seebeck coefficient as a function of chemical potential.

    Args:
        seebeck_obj (object): Seebeck object containing the data for the Seebeck coefficient.
        save_to_file (bool, optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to choose the direction to plot among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz.
        - Plots the Seebeck coefficient as a function of chemical potential for each temperature.
        - Distinguishes between n-type and p-type conduction with dashed and solid lines, respectively.
        - If save_to_file is True, saves the plot to a file named 'seebeck_potential_at_T_K___YYYY-MM-DD_HHMMSS.jpg' for each temperature, and 'seebeck_potential_different_T_YYYY-MM-DD_HHMMSS.jpg' for all temperatures combined.

    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syx':
        col = 6
    elif case == 'syy':
        col = 7
    elif case == 'syz':
        col = 8
    elif case == 'szx':
        col = 9
    elif case == 'szy':
        col = 10
    elif case == 'szz':
        col = 11

    else:
        sys.exit('please, choose a valid chioce')

    vol = seebeck_obj.volume

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    carrier = []
    for k in range(0, len(seebeck_obj.all_data)):
        carrier.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    y = []
    for k in range(0, len(seebeck_obj.all_data)):
        y.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col])*1000000)))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(carrier[k])):
            if carrier[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    endx = []
    endy = []

    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Seebeck at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.savefig('seebeck_potential_at_' + str(seebeck_obj.temp[k]) + 'K___' + time.strftime(
            "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()

    from matplotlib.pyplot import figure
    figure(figsize=(7, 7))
    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Seebeck at different T')
    plt.savefig('seebeck_potential_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    plt.show()
    # if save_to_file != False:
    #     save_plot(save_to_file)


def plot_cry_seebeck_carrier(seebeck_obj):
    """
    Plot the Seebeck coefficient as a function of charge carrier concentration.

    Args:
        seebeck_obj: Seebeck object containing the data for the Seebeck coefficient.
        save_to_file (optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to choose the direction to plot among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz.
        - Plots the Seebeck coefficient as a function of charge carrier concentration for each temperature, distinguishing between n-type and p-type conduction.
        - If save_to_file is True, saves the plot to a file named 'seebeck_carrier_at_T_K___YYYY-MM-DD_HHMMSS.jpg' for each temperature, and 'seebeck_carrier_different_T_YYYY-MM-DD_HHMMSS.jpg' for all temperatures combined.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syx':
        col = 6
    elif case == 'syy':
        col = 7
    elif case == 'syz':
        col = 8
    elif case == 'szx':
        col = 9
    elif case == 'szy':
        col = 10
    elif case == 'szz':
        col = 11
    else:
        sys.exit('please, choose a valid chioce')

    vol = seebeck_obj.volume

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))
    y = []

    for k in range(0, len(seebeck_obj.all_data)):
        y.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col])*1000000)))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(abs(np.array(xnegfin[k])), ynegfin[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Seebeck at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.xscale('log')
        plt.show()

    from matplotlib.pyplot import figure

    figure(figsize=(7, 7))
    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(abs(np.array(xnegfin[k])), ynegfin[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        plt.xscale('log')
        plt.title('Seebeck at different T')
    plt.savefig('seebeck_carrier_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    plt.show()

    # if save_to_file != False:
    #     save_plot(save_to_file)


def plot_cry_multiseebeck(*seebeck):
    """
    Plot the multiseebeck coefficient for different temperatures.

    Args:
        *seebeck: Variable number of seebeck objects containing the data for the Seebeck coefficient.

    Returns:
        None

    Notes:
        - Prompts the user to input the index of the temperature to plot.
        - Prompts the user to input the lower and higher values of chemical potential to plot in eV.
        - Prompts the user to choose the direction to plot among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz.
        - Plots the multiseebeck coefficient for each seebeck object.
        - Differentiates transport coefficients due to n-type or p-type conduction using dashed and solid lines.
        - Saves the plot to a file named 'multiseebeckYYYY-MM-DD_HHMMSS.jpg', where YYYY-MM-DD_HHMMSS represents the current date and time.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    k = int(input(
        'Insert the index of temperature you want to plot \n(i.e. if your temperature are [T1, T2, T3] indexes are [0, 1, 2])'))
    minpot = float(
        input('Insert the lower value of chemical potential you want to plot in eV'))
    maxpot = float(
        input('Inser the higher value of chemical potential you want to plot in eV'))

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_yz, S_zx, S_zy, S_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syx':
        col = 6
    elif case == 'syy':
        col = 7
    elif case == 'syz':
        col = 8
    elif case == 'szx':
        col = 9
    elif case == 'szy':
        col = 10
    elif case == 'szz':
        col = 11

    else:
        sys.exit('please, choose a valid chioce')

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    i = 0
    for n in seebeck:
        vol = n.volume

        x = []
        for kq in range(0, len(n.all_data)):
            x.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[0]))))

        carrier = []
        for kq in range(0, len(n.all_data)):
            carrier.append(np.array(n.all_data[kq].apply(
                lambda x: (float(x.split()[2])/vol))))

        y = []
        for kq in range(0, len(n.all_data)):
            y.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[col])*1000000)))

        yneg = []
        ypos = []
        xpos = []
        xneg = []
        yposfin = []
        xposfin = []
        ynegfin = []
        xnegfin = []

        for kq in range(0, len(n.all_data)):
            for j in range(0, len(carrier[kq])):
                if carrier[kq][j] >= 0:
                    xpos.append(x[kq][j])
                    ypos.append(y[kq][j])
                else:
                    xneg.append(x[kq][j])
                    yneg.append(y[kq][j])
            yposfin.append(ypos)
            ynegfin.append(yneg)
            xposfin.append(xpos)
            xnegfin.append(xneg)
            xpos = []
            ypos = []
            xneg = []
            yneg = []

        colours = ['royalblue', 'orange', 'green', 'red',
                   'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

        endx = []
        endy = []

        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[i])
        plt.plot(xposfin[k], yposfin[k], color=colours[i], label=str(n.title))
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[i])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)', fontsize=12)
        plt.xlim(minpot, maxpot)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        i = 1+i
    plt.title('MultiSeebeck ' + str(n.temp[k]) + ' K')
    plt.savefig('multiseebeck' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')


#-------------------------------------SIGMA-----------------------------------#

def plot_cry_sigma_potential(sigma_obj):
    """
    Plot the electrical conductivity as a function of chemical potential.

    Args:
        sigma_obj (object): Sigma object containing the data for electrical conductivity.
        save_to_file (bool, optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to choose the direction to plot among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz.
        - Plots the electrical conductivity as a function of chemical potential for each temperature.
        - Distinguishes between n-type and p-type conduction with dashed and solid lines, respectively.
        - If save_to_file is True, saves the plot to a file named 'sigma_potential_different_T_YYYY-MM-DD_HHMMSS.jpg'.

    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syy':
        col = 6
    elif case == 'syz':
        col = 7
    elif case == 'szz':
        col = 8
    else:
        sys.exit('please, choose a valid chioce')

    vol = sigma_obj.volume

    x = []
    for k in range(0, len(sigma_obj.all_data)):
        x.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    carrier = []
    for k in range(0, len(sigma_obj.all_data)):
        carrier.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    y = []
    for k in range(0, len(sigma_obj.all_data)):
        y.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(sigma_obj.all_data)):
        for j in range(0, len(x[k])):
            if carrier[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')
    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']
    endx = []
    endy = []

    for k in range(0, len(sigma_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(sigma_obj.temp[k])+' K')
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Sigma at '+str(sigma_obj.temp[k]) + 'K')
        plt.legend(loc='upper left', fontsize=12)
        plt.savefig('sigma_potential_at_' + str(sigma_obj.temp[k]) + 'K___' + time.strftime(
            "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()

    for k in range(0, len(sigma_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(sigma_obj.temp[k])+' K')
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.title('Sigma at different T')
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
    plt.savefig('sigma_potential_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')

    # if save_to_file != False:
    #     save_plot(save_to_file)


def plot_cry_sigma_carrier(sigma_obj):
    """
    Plot the electrical conductivity as a function of charge carrier concentration.

    Args:
        sigma_obj: Sigma object containing the data for the electrical conductivity.
        save_to_file (optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to choose the direction to plot among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz.
        - Plots the electrical conductivity as a function of charge carrier concentration for each temperature, distinguishing between n-type and p-type conduction.
        - If save_to_file is True, saves the plot to a file named 'sigma_carrier_at_T_K___YYYY-MM-DD_HHMMSS.jpg' for each temperature, and 'sigma_carrier_different_T_YYYY-MM-DD_HHMMSS.jpg' for all temperatures combined.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syy':
        col = 6
    elif case == 'syz':
        col = 7
    elif case == 'szz':
        col = 8
    else:
        sys.exit('please, choose a valid chioce')

    vol = sigma_obj.volume

    x = []
    for k in range(0, len(sigma_obj.all_data)):
        x.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    y = []
    for k in range(0, len(sigma_obj.all_data)):
        y.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    yneg = []
    ypos = []
    xpos = []
    xneg = []
    yposfin = []
    xposfin = []
    ynegfin = []
    xnegfin = []

    for k in range(0, len(sigma_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xpos.append(x[k][j])
                ypos.append(y[k][j])
            else:
                xneg.append(x[k][j])
                yneg.append(y[k][j])
        yposfin.append(ypos)
        ynegfin.append(yneg)
        xposfin.append(xpos)
        xnegfin.append(xneg)
        xpos = []
        ypos = []
        xneg = []
        yneg = []

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    for k in range(0, len(sigma_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(sigma_obj.temp[k])+' K')
        plt.plot(abs(np.array(xnegfin[k])), ynegfin[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Sigma at ' + str(sigma_obj.temp[k]) + 'K')
        plt.legend(loc='upper left', fontsize=12)
        plt.xscale('log')
        plt.savefig('sigma_carrier_at_' + str(sigma_obj.temp[k]) + 'K___' + time.strftime(
            "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()

    for k in range(0, len(sigma_obj.all_data)):
        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xposfin[k], yposfin[k], color=colours[k],
                 label=str(sigma_obj.temp[k])+' K')
        plt.plot(abs(np.array(xnegfin[k])), ynegfin[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.title('Sigma at different T')
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        plt.xscale('log')
    plt.savefig('sigma_carrier_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')

    # if save_to_file != False:
    #     save_plot(save_to_file)


def plot_cry_multisigma(*sigma):
    """
    Plot the multisigma conductivity for different temperatures.

    Args:
        *sigma: Variable number of sigma objects containing the data for the conductivity.

    Returns:
        None

    Notes:
        - Prompts the user to input the index of the temperature to plot.
        - Prompts the user to input the lower and higher values of chemical potential to plot in eV.
        - Prompts the user to choose the direction to plot among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz.
        - Plots the multisigma conductivity for each sigma object.
        - Differentiates transport coefficients due to n-type or p-type conduction using dashed and solid lines.
        - Saves the plot to a file named 'multisigmaYYYY-MM-DD_HHMMSS.jpg', where YYYY-MM-DD_HHMMSS represents the current date and time.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    k = int(input(
        'Insert the index of temperature you want to plot \n(i.e. if your temperature are [T1, T2, T3] indexes are [0, 1, 2])'))
    minpot = float(
        input('Insert the lower value of chemical potential you want to plot in eV'))
    maxpot = float(
        input('Inser the higher value of chemical potential you want to plot in eV'))

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yy, S_yz, S_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syy':
        col = 6
    elif case == 'syz':
        col = 7
    elif case == 'szz':
        col = 8

    else:
        sys.exit('please, choose a valid chioce')

    i = 0
    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')
    for n in sigma:
        vol = n.volume

        x = []
        for kq in range(0, len(n.all_data)):
            x.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[0]))))

        carrier = []
        for kq in range(0, len(n.all_data)):
            carrier.append(np.array(n.all_data[kq].apply(
                lambda x: (float(x.split()[2])/vol))))

        y = []
        for kq in range(0, len(n.all_data)):
            y.append(np.array(n.all_data[kq].apply(
                lambda x: float(x.split()[col]))))

        yneg = []
        ypos = []
        xpos = []
        xneg = []
        yposfin = []
        xposfin = []
        ynegfin = []
        xnegfin = []

        for kq in range(0, len(n.all_data)):
            for j in range(0, len(x[kq])):
                if carrier[kq][j] >= 0:
                    xpos.append(x[kq][j])
                    ypos.append(y[kq][j])
                else:
                    xneg.append(x[kq][j])
                    yneg.append(y[kq][j])
            yposfin.append(ypos)
            ynegfin.append(yneg)
            xposfin.append(xpos)
            xnegfin.append(xneg)
            xpos = []
            ypos = []
            xneg = []
            yneg = []

        colours = []
        colours = ['royalblue', 'orange', 'green', 'red',
                   'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']
        endx = []
        endy = []

        endx = [xposfin[k][-1], xnegfin[k][0]]
        endy = [yposfin[k][-1], ynegfin[k][0]]
        plt.plot(endx, endy, color=colours[i])
        plt.plot(xposfin[k], yposfin[k], color=colours[i], label=str(n.title))
        plt.plot(xnegfin[k], ynegfin[k], '--', color=colours[i])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Electrical Conductivity (S/m)', fontsize=12)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        plt.xlim(minpot, maxpot)
        i = 1+i
    plt.title('MultiSigma ' + str(sigma[0].temp[k]) + ' K')
    plt.savefig('multisigma' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')


#--------------------------------POWERFACTOR----------------------------------#

def plot_cry_powerfactor_potential(seebeck_obj, sigma_obj):
    """
    Plot the power factor for different potentials.

    Args:
        seebeck_obj: Seebeck object containing the data for the Seebeck coefficient.
        sigma_obj: Sigma object containing the data for the electrical conductivity.
        save_to_file (optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to choose the direction to plot among PF_xx, PF_xy, PF_xz, PF_yx, PF_yy, PF_yz, PF_yz, PF_zx, PF_zy, PF_zz.
        - Calculates the power factor using the Seebeck coefficient and electrical conductivity data for each temperature.
        - Plots the power factor for each temperature as a function of the chemical potential, distinguishing between n-type and p-type conduction.
        - If save_to_file is True, saves the plot to a file named 'powerfactor_potential_at_T_K___YYYY-MM-DD_HHMMSS.jpg' for each temperature, and 'powerfactor_potential_different_T_YYYY-MM-DD_HHMMSS.jpg' for all temperatures combined.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among PF_xx, PF_xy, PF_xz, PF_yx, PF_yy, PF_yz, PF_yz, PF_zx, PF_zy, PF_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'pfxx':
        col = 3
    elif case == 'pfxy':
        col = 4
    elif case == 'pfxz':
        col = 5
    elif case == 'pfyx':
        col = 6
    elif case == 'pfyy':
        col = 7
    elif case == 'pfyz':
        col = 8
    elif case == 'pfzx':
        col = 9
    elif case == 'pfzy':
        col = 10
    elif case == 'pfzz':
        col = 11
    else:
        sys.exit('please, choose a valid chioce')

    if case == 'pfxx':
        cols = 3
    elif case == 'pfxy':
        cols = 4
    elif case == 'pfxz':
        cols = 5
    elif case == 'pfyx':
        cols = 4
    elif case == 'pfyy':
        cols = 6
    elif case == 'pfyz':
        cols = 7
    elif case == 'pfzx':
        cols = 5
    elif case == 'pfzy':
        cols = 7
    elif case == 'pfzz':
        cols = 8
    else:
        sys.exit('please, choose a valid chioce')

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    vol = sigma_obj.volume
    carrier = []

    yse = []
    for k in range(0, len(seebeck_obj.all_data)):
        yse.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    ysi = []
    for k in range(0, len(sigma_obj.all_data)):
        ysi.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[cols]))))

    carrier = []
    for k in range(0, len(seebeck_obj.all_data)):
        carrier.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    ysineg = []
    ysipos = []
    xsipos = []
    xsineg = []
    ysiposfin = []
    xsiposfin = []
    ysinegfin = []
    xsinegfin = []

    yseneg = []
    ysepos = []
    xsepos = []
    xseneg = []
    yseposfin = []
    xseposfin = []
    ysenegfin = []
    xsenegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(carrier[k])):
            if carrier[k][j] >= 0:
                xsipos.append(x[k][j])
                ysipos.append(ysi[k][j])
            else:
                xsineg.append(x[k][j])
                ysineg.append(ysi[k][j])
        ysiposfin.append(np.array(ysipos))
        ysinegfin.append(np.array(ysineg))
        xsiposfin.append(np.array(xsipos))
        xsinegfin.append(np.array(xsineg))
        xsipos = []
        ysipos = []
        xsineg = []
        ysineg = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(carrier[k])):
            if carrier[k][j] >= 0:
                xsepos.append(x[k][j])
                ysepos.append(yse[k][j])
            else:
                xseneg.append(x[k][j])
                yseneg.append(yse[k][j])

        yseposfin.append(np.array(ysepos))
        ysenegfin.append(np.array(yseneg))
        xseposfin.append(np.array(xsepos))
        xsenegfin.append(np.array(xseneg))
        xsepos = []
        ysepos = []
        xseneg = []
        yseneg = []

    pf_meta_pos = []
    for i in range(0, len(yseposfin)):
        pf_meta_pos.append(yseposfin[i]*yseposfin[i])

    pf_pos = []
    for i in range(0, len(pf_meta_pos)):
        pf_pos.append(pf_meta_pos[i] * ysiposfin[i])

    pf_meta_neg = []
    for i in range(0, len(ysenegfin)):
        pf_meta_neg.append(ysenegfin[i] * ysenegfin[i])

    pf_neg = []
    for i in range(0, len(pf_meta_neg)):
        pf_neg.append(pf_meta_neg[i] * ysinegfin[i])

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xsiposfin[k][-1], xsinegfin[k][0]]
        endy = [pf_pos[k][-1], pf_neg[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(xsinegfin[k], pf_neg[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Power Factor at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.show()

    from matplotlib.pyplot import figure

    figure(figsize=(7, 7))
    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xsiposfin[k][-1], xsinegfin[k][0]]
        endy = [pf_pos[k][-1], pf_neg[k][0]]
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(xsinegfin[k], pf_neg[k], '--', color=colours[k])
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        plt.title('Power Factor at different T')

    plt.savefig('powerfactor_potential_different_T_' + time.strftime(
        "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    # if save_to_file != False:
    #     save_plot(save_to_file)


def plot_cry_powerfactor_carrier(seebeck_obj, sigma_obj):
    """
    Plot the power factor for different charge carrier concentrations.

    Args:
        seebeck_obj: Seebeck object containing the data for the Seebeck coefficient.
        sigma_obj: Sigma object containing the data for the electrical conductivity.
        save_to_file (optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to choose the direction to plot among PF_xx, PF_xy, PF_xz, PF_yx, PF_yy, PF_yz, PF_yz, PF_zx, PF_zy, PF_zz.
        - Calculates the power factor using the Seebeck coefficient and electrical conductivity data for each temperature.
        - Plots the power factor for each temperature as a function of the charge carrier concentration, distinguishing between n-type and p-type conduction.
        - If save_to_file is True, saves the plot to a file named 'powerfactor_carrier_at_T_K___YYYY-MM-DD_HHMMSS.jpg' for each temperature, and 'powerfactor_carrier_different_T_YYYY-MM-DD_HHMMSS.jpg' for all temperatures combined.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among PF_xx, PF_xy, PF_xz, PF_yx, PF_yy, PF_yz, PF_yz, PF_zx, PF_zy, PF_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'pfxx':
        col = 3
    elif case == 'pfxy':
        col = 4
    elif case == 'pfxz':
        col = 5
    elif case == 'pfyx':
        col = 6
    elif case == 'pfyy':
        col = 7
    elif case == 'pfyz':
        col = 8
    elif case == 'pfzx':
        col = 9
    elif case == 'pfzy':
        col = 10
    elif case == 'pfzz':
        col = 11
    else:
        sys.exit('please, choose a valid chioce')

    if case == 'pfxx':
        cols = 3
    elif case == 'pfxy':
        cols = 4
    elif case == 'pfxz':
        cols = 5
    elif case == 'pfyx':
        cols = 4
    elif case == 'pfyy':
        cols = 6
    elif case == 'pfyz':
        cols = 7
    elif case == 'pfzx':
        cols = 5
    elif case == 'pfzy':
        cols = 7
    elif case == 'pfzz':
        cols = 8
    else:
        sys.exit('please, choose a valid chioce')

    vol = sigma_obj.volume

    x = []
    for k in range(0, len(sigma_obj.all_data)):
        x.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: (float(x.split()[2])/vol))))

    yse = []
    for k in range(0, len(seebeck_obj.all_data)):
        yse.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    ysi = []
    for k in range(0, len(sigma_obj.all_data)):
        ysi.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[cols]))))

    pf_meta = []
    for i in range(0, len(yse)):
        pf_meta.append(yse[i] * yse[i])

    pf = []
    for i in range(0, len(pf_meta)):
        pf.append(pf_meta[i] * ysi[i])

    ysineg = []
    ysipos = []
    xsipos = []
    xsineg = []
    ysiposfin = []
    xsiposfin = []
    ysinegfin = []
    xsinegfin = []

    yseneg = []
    ysepos = []
    xsepos = []
    xseneg = []
    yseposfin = []
    xseposfin = []
    ysenegfin = []
    xsenegfin = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xsipos.append(x[k][j])
                ysipos.append(ysi[k][j])
            else:
                xsineg.append(x[k][j])
                ysineg.append(ysi[k][j])
        ysiposfin.append(np.array(ysipos))
        ysinegfin.append(np.array(ysineg))
        xsiposfin.append(np.array(xsipos))
        xsinegfin.append(np.array(xsineg))
        xsipos = []
        ysipos = []
        xsineg = []
        ysineg = []

    for k in range(0, len(seebeck_obj.all_data)):
        for j in range(0, len(x[k])):
            if x[k][j] >= 0:
                xsepos.append(x[k][j])
                ysepos.append(yse[k][j])
            else:
                xseneg.append(x[k][j])
                yseneg.append(yse[k][j])
        yseposfin.append(np.array(ysepos))
        ysenegfin.append(np.array(yseneg))
        xseposfin.append(np.array(xsepos))
        xsenegfin.append(np.array(xseneg))
        xsepos = []
        ysepos = []
        xseneg = []
        yseneg = []

    pf_meta_pos = []
    for i in range(0, len(yseposfin)):
        pf_meta_pos.append(yseposfin[i]*yseposfin[i])

    pf_pos = []
    for i in range(0, len(pf_meta_pos)):
        pf_pos.append(pf_meta_pos[i] * ysiposfin[i])

    pf_meta_neg = []
    for i in range(0, len(ysenegfin)):
        pf_meta_neg.append(ysenegfin[i] * ysenegfin[i])

    pf_neg = []
    for i in range(0, len(pf_meta_neg)):
        pf_neg.append(pf_meta_neg[i] * ysinegfin[i])

    print('To differentiate transport coefficients due to n-type or p-type conduction (electrons or holes as majority carriers) dashed and solid lines are used, respectively.')

    colours = []
    colours = ['royalblue', 'orange', 'green', 'red',
               'purple', 'brown', 'pink', 'grey', 'olive', 'cyan']

    from matplotlib.pyplot import figure

    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xsiposfin[k][-1], xsinegfin[k][0]]
        endy = [pf_pos[k][-1], pf_neg[k][0]]
        plt.figure()
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(abs(xsinegfin[k]), pf_neg[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('Power Factor at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.xscale('log')
        plt.savefig('powerfactor_carrier_at_' + str(seebeck_obj.temp[k]) + 'K___' + time.strftime(
            "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()

    from matplotlib.pyplot import figure
    figure(figsize=(7, 7))
    for k in range(0, len(seebeck_obj.all_data)):
        endx = [xsiposfin[k][-1], xsinegfin[k][0]]
        endy = [pf_pos[k][-1], pf_neg[k][0]]
        plt.plot(endx, endy, color=colours[k])
        plt.plot(xsiposfin[k], pf_pos[k], color=colours[k],
                 label=str(seebeck_obj.temp[k])+' K')
        plt.plot(abs(xsinegfin[k]), pf_neg[k], '--', color=colours[k])
        plt.xlabel('Charge Carrier Concentration (cm$^{-3}$)', fontsize=12)
        plt.ylabel('Power Factor (10$^{-12}$WK$^{-2}$m$^{-1}$)', fontsize=12)
        plt.title('Power Factor at different T')
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
        plt.xscale('log')
    plt.savefig('powerfactor_carrier_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    # plt.show()

    # if save_to_file != False:
    #     save_plot(save_to_file)

#-------------------------------------ZT--------------------------------------#


def plot_cry_zt(seebeck_obj, sigma_obj):
    """
    Plot the ZT value for different temperatures.

    Args:
        seebeck_obj: Seebeck object containing the data for the Seebeck coefficient.
        sigma_obj: Sigma object containing the data for the electrical conductivity.
        save_to_file (optional): If True, saves the plot to a file. Default is False.

    Returns:
        None

    Notes:
        - Prompts the user to input the value of ktot in W-1K-1m-1.
        - Prompts the user to choose the direction to plot among ZT_xx, ZT_xy, ZT_xz, ZT_yx, ZT_yy, ZT_yz, ZT_yz, ZT_zx, ZT_zy, ZT_zz.
        - Calculates the ZT value using the Seebeck coefficient and electrical conductivity data.
        - Plots the ZT value for each temperature as a function of the chemical potential.
        - If save_to_file is True, saves the plot to a file named 'zt_at_T_K___YYYY-MM-DD_HHMMSS.jpg' for each temperature, and 'zt_different_T_YYYY-MM-DD_HHMMSS.jpg' for all temperatures combined.
    """
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np

    ktot = float(input(
        'Please insert the value of ktot in W-1K-1m-1'))
    case = input(
        'Please, choose the direction you want to plot. \nYou can choose among ZT_xx, ZT_xy, ZT_xz, ZT_yx, ZT_yy, ZT_yz, ZT_yz, ZT_zx, ZT_zy, ZT_zz\n')

    case = case.lower().replace('_', '')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')

    if case == 'ztxx':
        col = 3
    elif case == 'ztxy':
        col = 4
    elif case == 'ztxz':
        col = 5
    elif case == 'ztyx':
        col = 6
    elif case == 'ztyy':
        col = 7
    elif case == 'ztyz':
        col = 8
    elif case == 'ztzx':
        col = 9
    elif case == 'ztzy':
        col = 10
    elif case == 'ztzz':
        col = 11
    else:
        sys.exit('please, choose a valid chioce')

    if case == 'ztxx':
        cols = 3
    elif case == 'ztxy':
        cols = 4
    elif case == 'ztxz':
        cols = 5
    elif case == 'ztyx':
        cols = 4
    elif case == 'ztyy':
        cols = 6
    elif case == 'ztyz':
        cols = 7
    elif case == 'ztzx':
        cols = 5
    elif case == 'ztzy':
        cols = 7
    elif case == 'ztzz':
        cols = 8
    else:
        sys.exit('please, choose a valid chioce')

    x = []
    for k in range(0, len(seebeck_obj.all_data)):
        x.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[0]))))

    yse = []
    for k in range(0, len(seebeck_obj.all_data)):
        yse.append(np.array(seebeck_obj.all_data[k].apply(
            lambda x: float(x.split()[col]))))

    ysi = []
    for k in range(0, len(sigma_obj.all_data)):
        ysi.append(np.array(sigma_obj.all_data[k].apply(
            lambda x: float(x.split()[cols]))))

    pf_meta = []
    for i in range(0, len(yse)):
        pf_meta.append(yse[i] * yse[i])

    pf = []
    for i in range(0, len(pf_meta)):
        pf.append(pf_meta[i] * ysi[i])

    zt = []
    for i in range(0, len(pf_meta)):
        zt.append((pf[i] * seebeck_obj.temp[i])/ktot)

    for k in range(0, len(seebeck_obj.all_data)):
        plt.figure()
        plt.plot(x[k], pf[k], label=str(seebeck_obj.temp[k])+' K')
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('ZT', fontsize=12)
        plt.axhline(0, color='k')
        plt.title('ZT at ' + str(seebeck_obj.temp[k]) + ' K')
        plt.legend(loc='upper left', fontsize=12)
        plt.savefig('zt_at_' + str(seebeck_obj.temp[k]) + 'K___' + time.strftime(
            "%Y-%m-%d_%H%M%S") + '.jpg', format='jpg', dpi=600, bbox_inches='tight')
        plt.show()

    for k in range(0, len(seebeck_obj.all_data)):
        plt.plot(x[k], pf[k], label=str(seebeck_obj.temp[k])+' K')
        plt.xlabel('Chemical Potential (eV)', fontsize=12)
        plt.ylabel('ZT', fontsize=12)
        plt.title('ZT at different T')
        plt.axhline(0, color='k')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
    plt.savefig('zt_different_T_' + time.strftime("%Y-%m-%d_%H%M%S") +
                '.jpg', format='jpg', dpi=100, bbox_inches='tight')
    # if save_to_file != False:
    #     save_plot(save_to_file)


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
    different systems. The function returns to lists of figure and axes objects
    for further processing. Only matplotlib is used for plotting. Plotly is
    not available.

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
        axes (Axes): Matplotlib Axes object. Same dimensionality as ``figs``.
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
    axes = [[0 for i in range(len(property))] for j in range(len(tensplt))]

    if uniform_scale == False: # Non-uniform scale
        for ip, p in enumerate(property):
            kwargs['property'] = p
            for it, t in enumerate(tensplt):
                fig, ax = t.plot_3D(**kwargs)
                if add_title == True:
                    ax[1].set_title(p)
                figs[it][ip] = fig
                axes[it][ip] = ax
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
                fig, ax = _plot3D_mplib(
                    R_all[it], X_all[it], Y_all[it], Z_all[it], scale_radius,
                    uplt_all[it], utext_all[it], platt_all[it], range_cbar,
                    range_x, range_y, range_z, Rref, **camera_args
                )
                if add_title == True:
                    ax[1].set_title(p)
                figs[it][ip] = fig
                axes[it][ip] = ax

    # dimensionality
    if len(tensplt) == 1 and len(property) == 1:
        figs = figs[0][0]
        axes = axes[0][0]
    elif len(tensplt) == 1 and len(property) > 1:
        figs = figs[0]
        axes = axes[0]
    elif len(tensplt) > 1 and len(property) == 1:
        figs = [i[0] for i in figs]
        axes = [i[0] for i in axes]

    return figs, axes


#--------------------------------2D ELASTIC----------------------------------#
def plot_elastics_2D(property, *tensor, same_fig_2D=True, uniform_scale_2D=True, **kwargs):
    """
    A wrapper function of :ref:`Tensor3D or Tensor2D <ref-elastics>` objects to
    plot 2D crystal elastic properties. The user can plot multiple properties
    for different systems. The function returns to lists of figure and axes
    objects for further processing. Base units: GPa, m.

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
        axes (Axes): Matplotlib Axes object. Same dimensionality as ``figs``.
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
        axes = [[0 for i in range(len(property))] for j in range(len(tensplt))]

    # set colors, nTensplt*nProperty
    clist = list(mcolors.TABLEAU_COLORS.keys())
    if len(property) == 1 and len(tensplt) > 1:
        colors = [[clist[i]] for i in range(len(tensplt))]
    else:
        colors = [[clist[i] for i in range(len(property))] for j in range(len(tensplt))]

    # plot
    if uniform_scale_2D == False: # Non-uniform scale
        for ip, p in enumerate(property):
            kwargs['property'] = p
            for it, t in enumerate(tensplt):
                if same_fig_2D == False:
                    fig, ax = t.plot_2D(**kwargs)
                    figs[it][ip] = fig
                    axes[it][ip] = ax
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
                    axes[it][ip] = ax
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
            axes = axes[0][0]
        elif len(tensplt) == 1 and len(property) > 1:
            figs = figs[0]
            axes = axes[0]
        elif len(tensplt) > 1 and len(property) == 1:
            figs = [i[0] for i in figs]
            axes = [i[0] for i in axes]

    return figs, axes


##############################################################################
#                                                                            #
#                             VIBRATIONAL PROPERTIES                         #
#                                                                            #
##############################################################################

#------------------------------------HARMONIC---------------------------------#

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
#--------------------------------obsolute functions----------------------------#
#------------------------------------------------------------------------------#

def plot_cry_ela(choose, ndeg, *args, dpi=200, filetype=".png",
                 transparency=False):
    """
    Obsolute. Use ``plot_elastics_3D``.
    """
    import warnings

    warnings.warn("You are calling an obsolute function. Use 'plot_elastics_3D' instead.",
                  stacklevel=2)
    args = [i for i in args]
    figs, axes = plot_elastics_3D(
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
    Obsolute. Use ``plot_electron_bands``.
    """
    import warnings

    warnings.warn("You are calling an obsolute function. Use 'plot_electron_bands' instead.",
                  stacklevel=2)
    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    if np.all(energy_range==None):
        energy_range=[]
    if np.all(k_range==None):
        k_range=[]
    if np.all(figsize==None):
        figsize=[6.4, 4.8]

    fig, ax = plot_electron_bands(
        *bands, unit=unit, k_label=k_labels, mode=mode, not_scaled=not_scaled,
        energy_range=energy_range, k_range=k_range, band_label=labels, band_color=color,
        band_linestyle=linestl, band_linewidth=linewidth, fermi_color=fermi,
        fermi_linewidth=fermiwidth, title=title, figsize=figsize, layout=scheme,
        sharex=sharex, sharey=sharey, fontsize=fontsize
    )
    return fig, ax


def plot_cry_band(bands, k_labels=[], energy_range=[], title=None, not_scaled=True,
                  mode='single', linestl='-', linewidth=1, color='blue',
                  fermi='forestgreen', k_range=[], labels=None, figsize=[6.4, 4.8],
                  scheme=None, sharex=True, sharey=True, fermiwidth=1.5, fermialpha=1):
    """
    Obsolute. Use ``plot_electron_bands``.
    """
    import warnings

    warnings.warn("You are calling an obsolute function. Use 'plot_electron_bands' instead.",
                  stacklevel=2)

    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    fig, ax = plot_electron_bands(
        *bands, k_label=k_labels, mode=mode, not_scaled=not_scaled,
        energy_range=energy_range, k_range=k_range, band_label=labels, band_color=color,
        band_linestyle=linestl, band_linewidth=linewidth, fermi_color=fermi,
        fermi_linewidth=fermiwidth, title=title, figsize=figsize, layout=scheme,
        sharex=sharex, sharey=sharey, fontsize=fontsize
    )
    return fig, ax


def plot_electron_dos(doss, unit='eV', beta='up', overlap=False, prj=None,
                      energy_range=None, dos_range=None, color='blue',
                      labels=None, linestl=None, linewidth=1, fermi='forestgreen',
                      title=None, figsize=None):
    """
    Obsolute. Use ``plot_electron_doss``.
    """
    import warnings

    warnings.warn("You are calling an obsolute function. Use 'plot_electron_doss' instead.",
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

    fig, ax = plot_electron_doss(
        *doss, unit=unit, beta=beta, overlap=overlap, prj=prj,
        energy_range=energy_range, dos_range=dos_range, dos_label=labels,
        dos_color=color, dos_linestyle=linestl, dos_linewidth=linewidth,
        fermi_color=fermi, title=title, figsize=figsize
    )
    return fig, ax


def plot_cry_doss(doss, color='blue', fermi='forestgreen', overlap=False,
                  labels=None, figsize=[6.4, 4.8], linestl=None,
                  linewidth=1.0, title=None, beta='down', energy_range=[],
                  dos_range=[], prj=[]):
    """
    Obsolute. Use ``plot_electron_doss``.
    """
    import warnings

    warnings.warn("You are calling an obsolute function. Use 'plot_electron_doss' instead.",
                  stacklevel=2)
    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    fig, ax = plot_electron_doss(
        *doss, beta=beta, overlap=overlap, prj=prj,
        energy_range=energy_range, dos_range=dos_range, dos_label=labels,
        dos_color=color, dos_linestyle=linestl, dos_linewidth=linewidth,
        fermi_color=fermi, title=title, figsize=figsize
    )
    return fig, ax


def plot_cry_es(bands, doss, k_labels=[], color_bd='blue', color_doss='blue',
                fermi='forestgreen', energy_range=[], linestl_bd=None,
                linestl_doss=None, linewidth=1.0, prj=[], figsize=[6.4, 4.8],
                labels=None, dos_range=[], title=None, dos_beta='down'):
    """
    Obsolute. Use ``plot_electron_doss``.
    """
    import warnings

    warnings.warn("You are calling an obsolute function. Use 'plot_electron_banddos' instead.",
                  stacklevel=2)

    fig, ax = plot_cry_es(
        bands, doss, k_label=k_labels, dos_beta=dos_beta, dos_prj=prj,
        energy_range=energy_range, dos_range=dos_range, band_color=color_bd,
        band_linestyle=linestl_bd, band_linewidth=linewidth, dos_label=labels,
        dos_color=color_doss, dos_linestyle=linestl_doss, dos_linewidth=linewidth,
        fermi_color=fermi, title=title, figsize=figsize
    )
    return fig, ax





