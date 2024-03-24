#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to visualize CRYSTAL outputs.
"""

##############################################################################
#                                                                            #
#                       ELECTRONIC STRUCTURE AND PHONONS                     #
#                                                                            #
##############################################################################

#-------------------------------ECHG charge density----------------------------#


def plot_dens_ECHG(obj_echg, unit='Angstrom', levels=150, xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None, dpi=400, name=None):
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
        name (str):  Name of the colormap. None for not saving it.
        dpi (int): *Valid if name!=None* Resolution (dots per inch) for the output image.
   Returns:
        None
    """
    import copy

    from CRYSTALpytools.base.plotbase import plot_2Dscalar

    obj = copy.deepcopy(obj_echg)
    obj._set_unit(unit)
    if unit.lower() == 'angstrom':
        cbarlabel = 'Charge Density ($|e|.\AA^{-3}$)'
    else:
        cbarlabel = 'Charge Density ($|e|.Bohr^{-3}$)'

    fig = plot_2Dscalar(obj.chgmap, obj.gridv, levels, xticks, yticks,
                        cmap_max, cmap_min, cbarlabel)

    # if name != None:
    #     save_plot(name, dpi=dpi)

    plt.show()
    return

#--------------------------------ECHG spin density-----------------------------#

def plot_spin_ECHG(obj_echg, unit='Angstrom', levels=150, xticks=5,
                   yticks=5, cmap_max=None, cmap_min=None, dpi=400, name=None):
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
        dpi (int): *Optional* Resolution (dots per inch) for the output image. Default is 400.
        savefig (bool): *Optional* Chose to save the figure or not. Default is False.
        name (str): *Valid if savefig=True* Name for the colormap.

   Returns:
        None
    """
    import copy

    from CRYSTALpytools.base.plotbase import plot_2Dscalar

    obj = copy.deepcopy(obj_echg)
    obj._set_unit(unit)
    if unit.lower() == 'angstrom':
        cbarlabel = 'Spin Density ($|e|.\AA^{-3}$)'
    else:
        cbarlabel = 'Spin Density ($|e|.Bohr^{-3}$)'

    fig = plot_2Dscalar(obj.spinmap, obj.gridv, levels, xticks, yticks,
                        cmap_max, cmap_min, cbarlabel)

    # if name != None:
    #     save_plot(name, dpi=dpi)

    plt.show()
    return

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
        transparency(bool): Background transparency of the saved file,

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

    if line_freq0 == None:
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

    # if save_to_file != None:
    #     save_plot(save_to_file, dpi=dpi, transparency=transparency)

    # plt.show()
    return fig, ax 


def plot_electron_band(bands, unit='eV', k_labels=None, mode='single',
                       not_scaled=False, energy_range=None, k_range=None,
                       color='blue', labels=None, linestl='-', linewidth=1,
                       fermi='forestgreen', fermiwidth=1.5, fermialpha=1, title=None, figsize=None,
                       scheme=None, sharex=True, sharey=True, fontsize=12):
    """
    A wrapper of plot_cry_bands for electron band structure.

    Args:
        bands (BandsBASE|list): Bands object generated by `CRYSTALpytools.crystal_io.Properties_output.read_bands` or
            a list of BandsBASE objects.
        unit (str): The unit of energy. Can be 'eV' or 'Hartree'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be, for example, 'Gamma'.
        mode (str): The plotting mode. Possible values are 'single', 'multi', and 'compare'.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
        energy_range (array): A 2x1 array specifying the energy range.
        k_range (array): A 2x1 array specifying the k-range.
        color (str|list): Color of plot lines. Should be consistent with bands.
        labels (str|list): Plot legend. Should be consistent with bands.
        linestl (str|list): Linestyle string. Should be consistent with bands.
        linewidth (float): The width of the plot lines.
        fermi (str): The color of the Fermi level line.
        fermiwidth (float): The width of the fermi line.
        fermialpha (float): Opacity of the fermi level 0-1.
        title (str): The title of the plot.
        figsize (list): The figure size specified as [width, height].
        scheme (list|tuple): The layout of subplots.
        sharex (bool): Whether to share the x-axis among subplots.
        sharey (bool): Whether to share the y-axis among subplots.
        save_to_file (str): The file name to save the plot.
        dpi (int): Dots per inch resolution of the saved file.
        fontsize (int): Fontsize of the axis labels 
        transparency: Background Transparency of the saved file.

    Returns:
        None

    :raise ValueError: If the specified unit is unknown.
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALpytools.base.plotbase import plot_cry_bands
    from CRYSTALpytools.units import H_to_eV, eV_to_H

    if re.match(r'^eV$', unit, re.IGNORECASE):
        unit = 'eV'
        is_ev = True
    elif re.match(r'^Hartree$', unit, re.IGNORECASE):
        unit = 'Hartree'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(bands, list) or isinstance(bands, tuple)):
        bands = [bands]

    for b in bands:
        if unit != b.unit:
            if unit == 'eV':
                b.bands[:, :, :] = H_to_eV(b.bands[:, :, :])
            else:
                b.bands[:, :, :] = eV_to_H(b.bands[:, :, :])
            b.unit = unit
    if len(bands) == 1:
        bands = bands[0]

    fig, ax = plot_cry_bands(bands, k_labels=k_labels, energy_range=energy_range, title=title,
                         not_scaled=not_scaled, mode=mode, linestl=linestl, linewidth=linewidth,
                         color=color, fermi=fermi, fermiwidth=fermiwidth, fermialpha=fermialpha, k_range=k_range, labels=labels,
                         figsize=figsize, scheme=scheme, sharex=sharex, sharey=sharey)
    if is_ev == True:
        fig.supylabel('$E-E_{F}$ (eV)', fontsize=fontsize)
    else:
        fig.supylabel('$E-E_{F}$ (Hartree)', fontsize=fontsize)

    # if save_to_file != None:
    #     save_plot(save_to_file, dpi=dpi, transparency=transparency)
    #
    # plt.show()
    return fig, ax


#-------------------------------DENSITY OF STATES-----------------------------#


def plot_electron_dos(doss, unit='eV', beta='up', overlap=False, prj=None,
                      energy_range=None, dos_range=None, color='blue',
                      labels=None, linestl=None, linewidth=1, fermi='forestgreen',
                      title=None, figsize=None):
    """
    A wrapper of plot_cry_doss for electron density of states.

    Args:
        doss (DOSBASE): DOS obect generated by code:`CRYSTALpytools.crystal_io.Properties_output.read_doss`.
            Or a list of DOSBASE objects.
        unit (str): 'eV' or 'Hartree'
        beta (str): Plot spin-down state 'up' or 'down'
        overlap (bool): Plotting multiple lines into the same figure
        prj (list): Index of selected projection. Consistent with the
            index of the 2nd dimension of :code:`doss.doss`
        energy_range (list[float]): 2*1 list of energy range
        dos_range (list[float]): 2*1 list of DOS range
        color (str | list[str]): Color of plot lines. *Should be
            consistent with number of projections.*
        labels (str | list[str]): Plot legend. *Should be consistent with
            number of projections.*
        linestl (str | list[str]): linestyle string. *Should be consistent
            with number of projections.*
        linewidth (float)
        fermi (str): Color of Fermi level line.
        title (str)
        figsize (list[float])
        save_to_file (str): File name.

    Returns:
        None
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALpytools.base.plotbase import plot_cry_doss
    from CRYSTALpytools.units import H_to_eV, eV_to_H

    if re.match(r'^ev$', unit, re.IGNORECASE):
        unit = 'eV'
        is_ev = True
    elif re.match(r'^hartree$', unit, re.IGNORECASE):
        unit = 'Hartree'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    if not (isinstance(doss, list) or isinstance(doss, tuple)):
        doss = [doss]

    for d in doss:
        if unit != d.unit:
            if unit == 'eV':
                d.doss[:, 0, :] = H_to_eV(d.doss[:, 0, :])
                d.doss[:, 1:, :] = eV_to_H(d.doss[:, 1:, :])
            else:
                d.doss[:, 0, :] = eV_to_H(d.doss[:, 0, :])
                d.doss[:, 1:, :] = H_to_eV(d.doss[:, 1:, :])
            d.unit = unit
    if len(doss) == 1:
        doss = doss[0]

    fig, ax = plot_cry_doss(doss, color=color, fermi=fermi, overlap=overlap,
                        labels=labels, figsize=figsize, linestl=linestl,
                        linewidth=linewidth, title=title, beta=beta,
                        energy_range=energy_range, dos_range=dos_range, prj=prj)
    if is_ev == True:
        fig.supylabel('DOS (states/eV)')
        fig.supxlabel('Energy (eV)')
    else:
        fig.supylabel('DOS (states/Hartree)')
        fig.supxlabel('Energy (Hartree)')

    # if save_to_file != None:
    #     save_plot(save_to_file)
    #
    # plt.show()
    return fig, ax


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
        save_to_file (str): File name.

    Returns:
        None
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

    if line_freq0 == None:
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

    # if save_to_file != None:
    #     save_plot(save_to_file)
    #
    # plt.show()
    return fig, ax


#-----------------------------BAND + DENSITY OF STATES------------------------#


def plot_electron_banddos(bands, doss, unit='eV', k_labels=None, dos_beta='down',
                          dos_prj=None, energy_range=None, dos_range=None,
                          color_band='blue', color_dos='blue', labels=None, linestl_band='-',
                          linestl_dos=None, linewidth=1, fermi='forestgreen',
                          title=None, figsize=None, legend=False):
    """
    A wrapper of plot_cry_es for electron band structure + dos. For spin-polarized cases, beta state.

    Args:
        bands (BandsBASE|list): Bands object generated by CRYSTALpytools.crystal_io.Properties_output.read_bands
            or a list of BandsBASE objects.
        doss (DOSBASE): DOS object generated by CRYSTALpytools.crystal_io.Properties_output.read_doss
            or a list of DOSBASE objects.
        unit (str): Unit of energy. Valid options are 'eV' or 'Hartree'.
        k_labels (list): A list of high-symmetric k point labels. Greek alphabets should be
            represented as strings, for example, 'Gamma'.
        dos_beta (str): Spin state to plot. Valid options are 'Up' or 'down'. If 'down', the beta
            state will be plotted on the same side as the alpha state, otherwise on the other side.
        dos_prj (list): Index of selected projection. Consistent with the index of the 2nd dimension
            of doss.doss.
        energy_range (list): A list of two values representing the energy range to be plotted.
        dos_range (list): DOS range for the y-axis.
        color_band (str): Color of the electron bands in the plot.
        color_dos (str): Color of the density of states (DOS) in the plot.
        labels (list): A list of labels for the plot legend.
        linestl_band (str): Linestyle of the electron bands.
        linestl_dos (str): Linestyle of the density of states (DOS).
        linewidth (float): Width of the lines in the plot.
        fermi (str): Color of the Fermi level line.
        title (str): Title of the plot.
        figsize (list[float]): Size of the figure in inches (width, height).
        save_to_file (str): File name to save the plot.
        legend (bool): Enables or disables the legend of the density of states (DOS).

    Returns:
        None

    :raise ValueError: If the unit parameter is unknown.
    """
    import re

    import matplotlib.pyplot as plt

    from CRYSTALpytools.base.plotbase import plot_cry_es
    from CRYSTALpytools.units import H_to_eV, eV_to_H

    if re.match(r'^ev$', unit, re.IGNORECASE):
        unit = 'eV'
        is_ev = True
    elif re.match(r'^hartree$', unit, re.IGNORECASE):
        unit = 'Hartree'
        is_ev = False
    else:
        raise ValueError('Unknown unit.')

    if unit != doss.unit:
        if unit == 'eV':
            doss.doss[:, 0, :] = H_to_eV(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = eV_to_H(doss.doss[:, 1:, :])
        else:
            doss.doss[:, 0, :] = eV_to_H(doss.doss[:, 0, :])
            doss.doss[:, 1:, :] = H_to_eV(doss.doss[:, 1:, :])
        doss.unit = unit
    if unit != bands.unit:
        if unit == 'eV':
            bands.bands[:, :, :] = H_to_eV(bands.bands[:, :, :])
        else:
            bands.bands[:, :, :] = eV_to_H(bands.bands[:, :, :])
        bands.unit = unit

    fig, ax = plot_cry_es(bands=bands, doss=doss, k_labels=k_labels, color_bd=color_band,
                      color_doss=color_dos, fermi=fermi, energy_range=energy_range,
                      linestl_bd=linestl_band, linestl_doss=linestl_dos,
                      linewidth=linewidth, prj=dos_prj, figsize=figsize, labels=labels,
                      dos_range=dos_range, title=title, dos_beta=dos_beta, legend=legend)
    if is_ev == True:
        fig.supylabel('Energy (eV)')
    else:
        fig.supylabel('Energy (Hartree)')

    # if save_to_file != None:
    #     save_plot(save_to_file)
    #
    # plt.show()
    return fig, ax


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
        save_to_file (str): File name to save the plot.

    Returns:
        None

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

    if line_freq0 == None:
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

    # if save_to_file != None:
    #     save_plot(save_to_file)
    #
    # plt.show()
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

#--------------------------------YOUNG MODULUS--------------------------------#

def plot_cry_young(theta, phi, S):
    """
    Compute Young's modulus for each direction of the space (i.e., each pair
    of theta and phi angles).

    Args:
        theta (float): Theta value.
        phi (float): Phi value.
        S (numpy.ndarray): Compliance matrix.

    Returns:
        float: Young's modulus values.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    # C2V = Matrix to refer the Cartesian into Voigt's notation
    # Observe that the matrix should be written as is shown below
    # C2V = np.array([[1,6,5],[6,2,4],[5,4,3]])
    # Since python start counting from zero all numbers must be subtracted by 1
    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2]
        ]
    )
    # print("The Matrix to convert Cartesian into Voigs Notation: \n", C2V)
    # creating the 1x3 vector "a"
    a = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )
    # e is a pseudo Young modulus value folowing the relation 1/e
    e = 0.0
    # i,j,k,l are updatable indices refering to cartesian notation that
    # will be converted by C2V into Voigt's
    for i in range(3):
        for j in range(3):
            v = C2V[i, j]
            for k in range(3):
                for l in range(3):
                    u = C2V[k, l]
                    # rf is a factor that must be multipled by the compliance element if
                    # certain conditions are satisfied
                    rf = 1
                    if v >= 3 and u >= 3:
                        rf = 4
                    if v >= 3 and u < 3:
                        rf = 2
                    if u >= 3 and v < 3:
                        rf = 2

                    rtmp = a[i] * a[j] * a[k] * a[l] * (S[v, u] / rf)
                    e = e + rtmp
    E_tmp = 1 / e  # is the Young Modulus of each cycle
    return E_tmp

#----------------------------COMPRESSION PROPERTIES---------------------------#


def plot_cry_comp(theta, phi, S):
    """
    Compute linear compressibility for each direction of the space (i.e., each
    pair of theta and phi angles).

    Args:
        theta (float): Theta value.
        phi (float): Phi value.
        S (numpy.ndarray): Compliance matrix.

    Returns:
        float: Linear compressibility values.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2]
        ]
    )

    a = np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )
    B = 0.0

    for i in range(3):
        for j in range(3):
            v = C2V[i, j]
            for k in range(3):
                u = C2V[k, k]
                rf = 1
                if v >= 3 and u >= 3:
                    rf = 4
                if v >= 3 and u < 3:
                    rf = 2
                if u >= 3 and v < 3:
                    rf = 2

                rtmp = a[i] * a[j] * (S[v, u] / rf)
                B = B + rtmp
    return B


#--------------------------------SHEAR MODULUS--------------------------------#

def plot_cry_shear(theta_1D, phi_1D, S, ndeg, shear_choice):
    """
    For each direction of the space (i.e., for each pair
    of theta and phi angles) the shear modulus is computed for the third angle
    chi and the average, maximum and minimum values are stored.

    Args:
        theta_1D (numpy.ndarray): One-dimensional array of theta values.
        phi_1D (numpy.ndarray): One-dimensional array of phi values.
        S (numpy.ndarray): Compliance matrix.
        ndeg (int): Number of degrees for discretization.
        shear_choice (str): Type of shear property to plot. Options: "avg", "min", "max".

    Returns:
        numpy.ndarray: Shear property array.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2],
        ]
    )
    shear_chi = np.zeros(ndeg)
    shear_min = np.zeros((ndeg, ndeg))
    shear_max = np.zeros((ndeg, ndeg))
    shear_avg = np.zeros((ndeg, ndeg))
    chi_1D = np.linspace(0, 2 * np.pi, ndeg)

    for phi_idx in range(ndeg):
        phi = phi_1D[phi_idx]
        for theta_idx in range(ndeg):
            theta = theta_1D[theta_idx]
            for chi_idx in range(ndeg):
                chi = chi_1D[chi_idx]
                a = np.array(
                    [
                        np.sin(theta) * np.cos(phi),
                        np.sin(theta) * np.sin(phi),
                        np.cos(theta),
                    ]
                )
                b = np.array(
                    [
                        np.cos(theta) * np.cos(phi) * np.cos(chi)
                        - np.sin(phi) * np.sin(chi),
                        np.cos(theta) * np.sin(phi) * np.cos(chi)
                        + np.cos(phi) * np.sin(chi),
                        -np.sin(theta) * np.cos(chi),
                    ]
                )
                shear_tmp = 0
                for i in range(3):
                    for j in range(3):
                        v = C2V[i, j]
                        for k in range(3):
                            for l in range(3):
                                u = C2V[k, l]
                                rf = 1
                                if v >= 3 and u >= 3:
                                    rf = 4
                                if v >= 3 and u < 3:
                                    rf = 2
                                if u >= 3 and v < 3:
                                    rf = 2
                                rtmp = a[i] * b[j] * a[k] * \
                                    b[l] * (S[v, u] / rf)
                                shear_tmp = shear_tmp + rtmp
                shear_chi[chi_idx] = 1 / (4 * shear_tmp)
            shear_min[phi_idx, theta_idx] = np.amin(shear_chi)
            shear_max[phi_idx, theta_idx] = np.amax(shear_chi)
            shear_avg[phi_idx, theta_idx] = np.mean(shear_chi)

    if shear_choice == "avg":
        return shear_avg
    if shear_choice == "min":
        return shear_min
    if shear_choice == "max":
        return shear_max


#------------------------------------POISSON RATIO----------------------------#

def plot_cry_poisson(theta_1D, phi_1D, S, ndeg, poisson_choice):
    """
    For each direction of the space (i.e., for each pair
    of theta and phi angles) the Poisson ratio is computed for the third angle
    chi and the average, maximum and minimum values are stored.

    Args:
        theta_1D (numpy.ndarray): One-dimensional array of theta values.
        phi_1D (numpy.ndarray): One-dimensional array of phi values.
        S (numpy.ndarray): Compliance matrix.
        ndeg (int): Number of degrees for discretization.
        poisson_choice (str): Type of Poisson's ratio to plot. Options: "avg", "min", "max".

    Returns:
        numpy.ndarray: Poisson's ratio array.

    Notes:
        - This function is intended to be called by cry_ela_plot
    """
    import numpy as np

    C2V = np.array(
        [
            [0, 5, 4],
            [5, 1, 3],
            [4, 3, 2],
        ]
    )
    poisson_chi = np.zeros(ndeg)
    poisson_min = np.zeros((ndeg, ndeg))
    poisson_max = np.zeros((ndeg, ndeg))
    poisson_avg = np.zeros((ndeg, ndeg))
    chi_1D = np.linspace(0, 2 * np.pi, ndeg)

    for phi_idx in range(ndeg):
        phi = phi_1D[phi_idx]
        for theta_idx in range(ndeg):
            theta = theta_1D[theta_idx]
            for chi_idx in range(ndeg):
                chi = chi_1D[chi_idx]
                a = np.array(
                    [
                        np.sin(theta) * np.cos(phi),
                        np.sin(theta) * np.sin(phi),
                        np.cos(theta),
                    ]
                )
                b = np.array(
                    [
                        np.cos(theta) * np.cos(phi) * np.cos(chi)
                        - np.sin(phi) * np.sin(chi),
                        np.cos(theta) * np.sin(phi) * np.cos(chi)
                        + np.cos(phi) * np.sin(chi),
                        -np.sin(theta) * np.cos(chi),
                    ]
                )
                poisson_num = 0
                poisson_den = 0
                for i in range(3):
                    for j in range(3):
                        v = C2V[i, j]
                        for k in range(3):
                            for l in range(3):
                                u = C2V[k, l]
                                rf = 1
                                if v >= 3 and u >= 3:
                                    rf = 4
                                if v >= 3 and u < 3:
                                    rf = 2
                                if u >= 3 and v < 3:
                                    rf = 2
                                num = (a[i] * a[j] * b[k] * b[l] * S[v, u])/rf
                                den = (a[i] * a[j] * a[k] * a[l] * S[v, u])/rf
                                poisson_num = poisson_num + num
                                poisson_den = poisson_den + den
                poisson_chi[chi_idx] = - poisson_num / poisson_den
            poisson_min[phi_idx, theta_idx] = np.amin(poisson_chi)
            poisson_max[phi_idx, theta_idx] = np.amax(poisson_chi)
            poisson_avg[phi_idx, theta_idx] = np.mean(poisson_chi)

    if poisson_choice == "avg":
        return poisson_avg
    if poisson_choice == "min":
        return poisson_min
    if poisson_choice == "max":
        return poisson_max


#----------------------------------ELASTIC------------------------------------#

def plot_cry_ela(choose, ndeg, *args, dpi=200, filetype=".png",
                 transparency=False):
    """
    Plot crystal elastic properties on the basis of the elastic tensor. A
    variable number of elastic tensors can be provided in order to get
    multiple plots in one shot, establishing a fixed color scale among them.

    Args:
        choose (str): Property to plot. Options: "young", "comp", "shear avg", 
        "shear min", "shear max", "poisson avg", "poisson min", "poisson max".
        ndeg (int): Number of degrees for discretization.
        *args: Variable number of elastic tensors.
        dpi (int, optional): Dots per inch for saving the plot. Default is 200.
        filetype (str, optional): File format of the output plot. Default is "png".
        transparency (bool, optional): Flag indicating whether to make the plot 
        background transparent. Default is False.

    Returns:
        None
    """
    import math
    import sys
    import time

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import animation, cm, colors
    from mpl_toolkits.mplot3d import Axes3D, axes3d

    i = 0
    R = [None] * len(args)
    tmin = []
    tmax = []

    # Compute elastic properties for each tensor -->
    for C in args:

        # Inverse of the matrix C in GPa (Compliance)
        S = np.linalg.inv(C)

        # One dimentional array of theta from 0 to pi
        theta_1D = np.linspace(0, np.pi, ndeg)
        # One dimentional array of phi from 0 to 2pi
        phi_1D = np.linspace(0, 2 * np.pi, ndeg)
        # Make a 2D array for theta and phi
        theta_2D, phi_2D = np.meshgrid(theta_1D, phi_1D)

        # Call to function
        if choose == "young":
            R[i] = plot_cry_young(theta_2D, phi_2D, S)
        elif choose == "comp":
            R[i] = plot_cry_comp(theta_2D, phi_2D, S)
        elif choose == "shear avg":
            R[i] = plot_cry_shear(theta_1D, phi_1D, S, ndeg, "avg")
        elif choose == "shear min":
            R[i] = plot_cry_shear(theta_1D, phi_1D, S, ndeg, "min")
        elif choose == "shear max":
            R[i] = plot_cry_shear(theta_1D, phi_1D, S, ndeg, "max")
        elif choose == "poisson avg":
            R[i] = plot_cry_poisson(theta_1D, phi_1D, S, ndeg, "avg")
        elif choose == "poisson min":
            R[i] = plot_cry_poisson(theta_1D, phi_1D, S, ndeg, "min")
        elif choose == "poisson max":
            R[i] = plot_cry_poisson(theta_1D, phi_1D, S, ndeg, "max")

        i += 1
    # <--

    # Find highest and lowest values -->
    for k in range(i):
        tmin.append(np.min(R[k]))
        tmax.append(np.max(R[k]))
    vmin = min(tmin)
    vmax = max(tmax)
    # <--

    # Create plot for each tensor -->
    for k in range(i):
        X = R[k] * np.sin(theta_2D) * np.cos(phi_2D)
        Y = R[k] * np.sin(theta_2D) * np.sin(phi_2D)
        Z = R[k] * np.cos(theta_2D)

        norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)
        fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))

        ax.plot_surface(
            X,
            Y,
            Z,
            rstride=1,
            cstride=1,
            facecolors=cm.jet(norm(R[k])),
            antialiased=True,
            alpha=0.75,
        )

        m = cm.ScalarMappable(cmap=cm.jet, norm=norm)
        m.set_array(R[k])
        fig.colorbar(m, ax=ax, shrink=0.7, location="left")

        # Make the planes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Make the grid lines transparent
        #  ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        #  ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        #  ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        # Fixing limits
        ax.set_xlim(-1 * np.max(R), np.max(R))
        ax.set_ylim(-1 * np.max(R), np.max(R))
        ax.set_zlim3d(-1 * np.max(R), np.max(R))
        ax.locator_params(nbins=5)  # tight=True,
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        ax.set_box_aspect(aspect=(1, 1, 1))  # Fix aspect ratio

        plt.show()
        fig.savefig(choose + time.strftime("%Y-%m-%d_%H%M%S.") +
                    filetype, dpi=dpi, transparent=transparency)

        # <--


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

    # if save_to_file != None:
    #     save_plot(save_to_file, dpi, transparency)
    #
    # plt.show()
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

    # if save_to_file != None:
    #     save_plot(save_to_file, dpi, transparency)
    #
    # plt.show()
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


##############################################################################
#                                                                            #
#                             COMMON FUNCTIONS                               #
#                                                                            #
##############################################################################


# def save_plot(path_to_file, dpi, transparency):
#     """
#     Save the plot as a file.
#
#     Args:
#         path_to_file (str): Path to the output file.
#         format (str, optional): File format of the output plot. Default is 'png'.
#
#     Raises:
#         FileNotFoundError: If the specified folder does not exist.
#
#     Returns:
#         None
#     """
#     import warnings
#     from os import path
#
#     import matplotlib.pyplot as plt
#
#     folder = path.split(path_to_file)[0]
#     file = path.split(path_to_file)[1]
#     extension = path.splitext(file)[-1]
#     extension_list = ['.png', '.jpg', '.jpeg', '.tif', '.pdf', '.svg', '.eps']
#
#     if folder == '':
#         folder = '.'
#     if extension != '':
#         if extension in extension_list:
#             format = extension[1:]
#         else:
#             warnings.warn('Unrecognized file format. PNG format is used.',
#                           stacklevel=2)
#
#     if path.exists(folder) == True:
#         plt.savefig('%s/%s.%s' % (folder, file, format),
#                     dpi=dpi, transparent=transparency)
#     else:
#         raise FileNotFoundError('Folder %s does not exist' % path_to_file)
