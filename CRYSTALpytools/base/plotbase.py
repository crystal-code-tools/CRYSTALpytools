#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base functions for plotting 2D and 3D figures
"""

def plot_overlap_bands(ax, bands, k_path, k_label, energy_range, k_range,
                       band_label, band_color, band_linestyle, band_linewidth,
                       fermi, fermi_color, fermi_linestyle, fermi_linewidth, **kwargs):
    """
    The plotting function for overlapped band structures of electron or phonon.
    Also can be used to get a single band.

    .. note::

        You must specify colors in string. List of RGB values are not allowed
        and might lead to unexpected errors.

    Args:
        ax (Axes): Matplotlib Axes object. To be passed from wrapper functions.
        bands (list): Band structure, 1\*nSystem list of nBand\*nKpoint\*nSpin
            numpy arrays.
        k_path (list): Coordinates of high-symmetric k points, 1\*nSystem list
            of 1\*nTick numpy arrays. Unit: :math:`\\AA^{-1}`.
        k_label (list[str] | None): 1\*nTick list of strings of the label for
            high symmetry points along the path. `mathtext <https://matplotlib.org/stable/users/explain/text/mathtext.html>`_
            experssions can also be used as in matplotlib.
        energy_range (list): 1\*2 list of plotting energy range.
        k_range (list): 1\*2 list of plotting k range. Can either be length
            (float) or k label (str). Must be used with ``not_scaled=False``
            and the same set of ``k_label``.
        band_label (list): 1\*nSystem or nSystem\*2 (spin) plot legend. If
            spin>1 and 1\*nSystem list is used, they are marked with the same
            label.
        band_color (list): 1\*nSystem or nSystem\*2 (spin) plot color. If spin
            >1 and 1\*nSystem list is used, they are in the same color.
        band_linestyle (list): 1\*nSystem or nSystem\*2 (spin) linestyle string.
            If spin>1 and 1\*nSystem list is used, they are in the same style.
        band_linewidth (list): 1\*nSystem or nSystem\*2 (spin) width of the plot
            lines. If spin>1 and 1\*nSystem list is used, they are in the same
            width.
        fermi (float | None): Fermi energy. By default the band is aligned to
            0. Can be used to offset the band. None for not plotting Fermi
            level.
        fermi_color (str): Color of the Fermi level.
        fermi_linestyle (str): Line style of Fermi level.
        fermi_linewidth(float): Width of the Fermi level.
        \*\*kwargs: Other commands passed to matplotlib ``axes.plot()`` method
            when plotting bands. Applied to all bands.

    Returns:
        ax (Axes): Matplotlib Axes object.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import copy

    nsys = len(bands)

    if (len(k_label) != 0) and (not isinstance(k_label[0], str)) and (not isinstance(k_label[0], float)):
        raise ValueError("K labels must be string or float number and only 1 set of k labels can be set for 'multi' mode.")

    # preprocessing, always scale k path
    k_path, k_label, energy_range, k_range, commands = _plot_bands_preprocess(
        bands, k_path, k_label, False, energy_range, k_range,
        band_label, band_color, band_linestyle, band_linewidth
    )

    # Start plotting
    ## Fermi level
    if np.all(fermi!=None):
        ax.hlines(fermi, k_range[0], k_range[1], color=fermi_color,
                  linestyle=fermi_linestyle, linewidth=fermi_linewidth)

    ## high symmetry lines
    for k in k_path[0]:
        ax.vlines(k, energy_range[0], energy_range[1], color='k', linewidth=0.5)

    ## bands
    keywords = ['label', 'color', 'linestyle', 'linewidth']

    ilabel = []
    countlabel = 0
    for isys in range(nsys):
        bandsplt = copy.deepcopy(bands[isys])
        if np.all(fermi!=None):
            bandsplt = bandsplt + fermi
            energy_range = energy_range + fermi

        nband, nkpt, nspin = bandsplt.shape
        k_pathplt = np.linspace(np.min(k_path[isys]), np.max(k_path[isys]), nkpt)
        for ispin in range(nspin):
            for icmd in range(4):
                if np.all(commands[icmd]==None):
                    continue
                kwargs[keywords[icmd]] = commands[icmd][isys][ispin]

            ax.plot(k_pathplt, bandsplt[:, :, ispin].transpose(), **kwargs)
            # a label for a set of bands, dimension of bandsplt array might vary
            countlabel = countlabel + nband
        ilabel.append(countlabel-1)

    # a label for a set of bands, dimension of bandsplt array might vary
    if np.all(commands[0]!=None):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend([handles[i] for i in ilabel],
                  [labels[i] for i in ilabel],
                  loc='center right')

    ax.set_xticks(k_path[0], labels=k_label[0])
    ax.set_xlim(k_range)
    ax.set_ylim(energy_range)
    return ax


def plot_compare_bands(ax, bands, k_path, k_label, not_scaled, energy_range, k_range,
                       band_label, band_color, band_linestyle, band_linewidth,
                       fermi, fermi_color, fermi_linestyle, fermi_linewidth,
                       **kwargs):
    """
    The plotting function for band structures of electron or phonon in different
    panels.

    .. note::

        You must specify colors in string. List of RGB values are not allowed
        and might lead to unexpected errors.

    Args:
        ax (Axes): Matplotlib Axes object or a flatted list of them. To be
            passed from wrapper functions.
        bands (list): Band structure, 1\*nSystem list of nBand\*nKpoint\*nSpin
            numpy arrays.
        k_path (list): Coordinates of high-symmetric k points, 1\*nSystem list
            of 1\*nTick numpy arrays. Unit: :math:`\\AA^{-1}`.
        k_label (list): nSystem\*nTick or 1\*nTick list of strings of the label
             for high symmetry points along the path. If a 1D list is given,
             the same labels are used for all the systems. `mathtext <https://matplotlib.org/stable/users/explain/text/mathtext.html>`_
             experssions can also be used  as in matplotlib.
        not_scaled (bool): Whether to scale the x-axis for different volumes.
        energy_range (list): 1\*2 list of plotting energy range.
        k_range (list): 1\*2 list of plotting k range. Can either be length
            (float) or k label (str). Must be used with ``not_scaled=False``
            and the same set of ``k_label``.
        band_label (list): 1\*nSystem or nSystem\*2 (spin) plot legend. If
            spin>1 and 1\*nSystem list is used, they are marked with the same
            label.
        band_color (list|str): A color string or 1\*2 color list for spin. If
            nSpin>1 but a string is given, same color for both spins.
        band_linestyle (list|str): A linestyle string or 1\*2 linestyle list
            for spin.  If nSpin>1 but a string is given, same linestyle for
            both spins.
        band_linewidth (list|float): A linewidth number or 1\*2 linewidth list
            for spin. If nSpin>1 but a string is given, same linewidth for both
            spins.
        fermi (float | None): Fermi energy. By default the band is aligned to
            0. Can be used to offset the band. None for not plotting Fermi
            level.
        fermi_color (str): Color of the Fermi level.
        fermi_linestyle (str): Line style of Fermi level.
        fermi_linewidth(float): Width of the Fermi level.
        \*\*kwargs: Other commands passed to matplotlib ``axes.plot()`` method
            when plotting bands. Applied to all bands.

    Returns:
        ax (Axes): Matplotlib Axes object or a flatted list of them.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import copy

    nsys = len(bands)

    # prepare band plot keywords
    ## label, use overlap plottings
    ## color. New Default
    if np.all(band_color==None):
         band_color = [['tab:blue', 'tab:blue'] for i in range(nsys)]
    ## line style, use overlap plottings
    ## line width, use overlap plottings

    # prepare fermi level
    if np.all(fermi!=None):
        if isinstance(fermi, float) or isinstance(fermi, int):
            fermi = np.array([fermi for i in range(nsys)], dtype=float)
        else:
            if len(fermi) != nsys:
                raise ValueError('Inconsistent numbers of Fermi level and systems')

    # preprocessing
    k_path, k_label, energy_range, k_range, commands  = _plot_bands_preprocess(
        bands, k_path, k_label, not_scaled, energy_range, k_range,
        band_label, band_color, band_linestyle, band_linewidth
    )

    # Start plotting
    keywords = ['label', 'color', 'linestyle', 'linewidth']
    bandsplt = copy.deepcopy(bands)

    for isys in range(nsys):
        bandsplt = copy.deepcopy(bands[isys])
        nband, nkpt, nspin = bandsplt.shape
        ## Fermi level
        if np.all(fermi!=None):
            ax[isys].hlines(fermi[isys], k_range[0], k_range[1], color=fermi_color,
                            linestyle=fermi_linestyle, linewidth=fermi_linewidth)
            bandsplt = bandsplt + fermi[isys]

        ## high symmetry lines
        for k in k_path[isys]:
            ax[isys].vlines(k, energy_range[0], energy_range[1], color='k', linewidth=0.5)

        ## bands
        k_pathplt = np.linspace(np.min(k_path[isys]), np.max(k_path[isys]), nkpt)
        for ispin in range(nspin):
            for icmd in range(4):
                if np.all(commands[icmd]==None): # 4*nsys*2(spin)
                    continue
                kwargs[keywords[icmd]] = commands[icmd][isys][ispin]

            ax[isys].plot(k_pathplt, bandsplt[:, :, ispin].transpose(), **kwargs)

        # a label for a set of bands
        if np.all(commands[0]!=None):
            handles, labels = ax[isys].get_legend_handles_labels()
            ilabel = [int(i*nband) for i in range(nspin)]
            ax[isys].legend([handles[i] for i in ilabel],
                            [labels[i] for i in ilabel],
                            loc='lower right')

        ax[isys].set_xticks(k_path[isys], labels=k_label[isys])
        ax[isys].set_xlim(k_range)
        ax[isys].set_ylim(energy_range)
    return ax


def _plot_bands_preprocess(
    bands, k_path, k_label, not_scaled, energy_range, k_range,
    band_label, band_color, band_linestyle, band_linewidth):
    """
    Do the boring parameters checking jobs for band structures. For the meanings
    of parameters, refer to ``plot_overlap_band()`` (``plot_compare_bands`` has
    less strict requirements).

    ``not_scaled`` is a flag to set whether to set the same length of k pathes
    of different systems.
    """
    import numpy as np
    import matplotlib.colors as mcolors

    nsys = len(bands)

    # For compatibility with old versions
    greek = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta',
             'Theta', 'Iota', 'Kappa', 'Lambda', 'Mu', 'Nu', 'Csi', 'Omicron',
             'Pi', 'Rho', 'Sigma', 'Tau', 'Upsilon', 'Phi', 'Chi', 'Psi',
             'Omega', 'Sigma_1']

    # Prepare k_label
    if len(k_label) != 0:
        if isinstance(k_label[0], str) or isinstance(k_label[0], float): # same definition for all k pathes
            same_klabel = True
            if len(k_label) != np.shape(k_path)[1]:
                raise ValueError('Inconsistent dimensions of k label and k path.')
            k_label = [k_label for i in range(nsys)]
        else:
            same_klabel = False
            for i in range(nsys):
                if len(k_label[i]) != len(k_path[i]):
                    raise ValueError('Inconsistent dimensions of k label and k path.')

        for i, listk in enumerate(k_label):
            for j, k in enumerate(listk):
                k = str(k)
                if k in greek:
                    if k != 'Sigma_1':
                        k_label[i][j] = r'$\{}$'.format(k)
                    else:
                        k_label[i][j] = r'$\Sigma_{1}$'
                else:
                    k_label[i][j] = k
    else:
        k_label = [[str(np.round(j, 2)) for j in i] for i in k_path]

    # scale k path to the longest one
    if not_scaled == False:
        alllen = []
        for i in range(nsys):
            alllen.append(np.max(k_path[i]) - np.min(k_path[i]))
        maxlen = np.max(alllen)
        for i in range(nsys):
            k_path[i] = k_path[i] / alllen[i] * maxlen

    # Prepare energy and k range
    if len(energy_range) != 0:
        energy_range = np.array([np.min(energy_range), np.max(energy_range)], dtype=float)
    else:
        energy_range = np.array([np.min(bands), np.max(bands)], dtype=float)

    if len(k_range) != 0:
        k_range = k_range[0:2]
        if isinstance(k_range[0], str):
            if not_scaled == True or same_klabel != True:
                raise Exception('You must scale k range and use the same set of k labels when restricting them.')
            for i in range(2):
                if k_range[i] in greek:
                    if k_range[i] != 'Sigma_1':
                        k_range[i] = r'$\{}$'.format(k_range[i])
                    else:
                        k_range[i] = r'$\Sigma_{1}$'
                else:
                    k_range[i] = k_range[i]

            # label must be found for every system
            found0 = False
            found1 = False
            for i in range(len(k_label[0])):
                if k_range[0] == k_label[0][i] and found0 != True:
                    found0 = True
                    k_range[0] = k_path[0][i]
                elif k_range[1] == k_label[0][i] and found1 != True:
                    found1 = True
                    k_range[1] = k_path[0][i]

            if not (found0 == True and found1 == True):
                raise Exception('Labelled k range is not found! Check your k_path, k_label and k_range inputs.')
        # for both label and number k ranges
        k_range = np.array([np.min(k_range), np.max(k_range)], dtype=float)
    else:
        k_range = np.array([np.min(k_path), np.max(k_path)], dtype=float)

    # Prepare band plot setups
    if np.all(band_label!=None):
        if isinstance(band_label, str):
            band_label = [[band_label, band_label] for i in range(nsys)]
        else:
            if len(band_label) != nsys:
                raise ValueError('Inconsistent system labels and number of systems.')
            for i in range(nsys):
                nband, nkpt, nspin = bands[i].shape
                ## label, and default setups of band labels
                if not isinstance(band_label[i], list):
                    if nspin == 2:
                        band_label[i] = [r'{} ($\alpha$)'.format(band_label[i]),
                                         r'{} ($\beta$)'.format(band_label[i])]
                    else:
                        band_label[i] = [band_label[i], band_label[i]]
                else:
                    band_label[i] = band_label[i][0:2]

    else:
        band_label = []; any_spin = False
        for i in range(nsys):
            nband, nkpt, nspin = bands[i].shape
            if nspin == 2:
                any_spin = True
                band_label.append([r'$\alpha$', r'$\beta$'])
            else:
                band_label.append(['', ''])

        if any_spin == False:
            band_label = None
    ## color
    if np.all(band_color!=None):
        if isinstance(band_color, str):
            band_color = [[band_color, band_color] for i in range(nsys)]
        else:
            if len(band_color) != nsys:
                raise ValueError('Inconsistent band colors and number of systems.')
            if not isinstance(band_color[0], list):
                band_color = [[i, i] for i in band_color]
            else:
                band_color = [[i[0], i[1]] for i in band_color]
    else: # defalut setups of band color
        clist = list(mcolors.TABLEAU_COLORS.keys())
        band_color = [[clist[i], clist[i]] for i in range(nsys)]
    ## line style
    if np.all(band_linestyle!=None):
        if isinstance(band_linestyle, str):
            band_linestyle = [[band_linestyle, band_linestyle] for i in range(nsys)]
        else:
            if len(band_linestyle) != nsys:
                raise ValueError('Inconsistent band line style and number of systems.')
            if not isinstance(band_linestyle[0], list):
                band_linestyle = [[i, i] for i in band_linestyle]
            else:
                band_linestyle = [[i[0], i[1]] for i in band_linestyle]
    else: # defalut setups of line style
        band_linestyle = [['-', '--'] for i in range(nsys)]
    ## linewidth
    if np.all(band_linewidth!=None):
        if isinstance(band_linewidth, int) or isinstance(band_linewidth, float):
            band_linewidth = [[band_linewidth, band_linewidth] for i in range(nsys)]
        else:
            if len(band_linewidth) != nsys:
                raise ValueError('Inconsistent band line width and number of systems.')
            if not isinstance(band_linewidth[0], list):
                band_linewidth = [[i, i] for i in band_linewidth]
            else:
                band_linewidth = [[i[0], i[1]] for i in band_linewidth]
    else: # defalut setups of linewidth
        band_linewidth = [[1.0, 1.0] for i in range(nsys)]

    commands = [band_label, band_color, band_linestyle, band_linewidth] # ncmd\*nsys\*2(spin)

    return k_path, k_label, energy_range, k_range, commands


def plot_doss(ax, doss, energy, beta, prj, energy_range, dos_range,
              dos_label, dos_color, dos_linestyle, dos_linewidth,
              fermi, fermi_color, fermi_linestyle, fermi_linewidth, **kwargs):
    """
    The base function to plot electron / phonon density of states on one axes.

    Args:
        ax (Axes): Matplotlib Axes object.
        doss (numpy.ndarray): nProj\*nEnergy\*nSpin array of DOS. Positive
            values for both spin up and spin down states.
        energy (numpy.ndarray): 1\*nEnergy array of energy.
        beta (str): Plot settings for :math:`\beta` states ('up' or 'down').
        prj (list): Index of selected projections, consistent with the first
            dimension of the ``doss``, starting from 1.
        energy_range (list): 1\*2 list of energy range.
        dos_range (list): 1\*2 list of DOS range.
        dos_label (list): 1\*nProj or nProj\*2 (spin) plot legend. If spin>1
            and 1\*nProj list is used, they are marked with the same label.
        dos_color (list): 1\*nProj or nProj\*2 (spin) plot color. If spin>1 and
            1\*nProj list is used, they are in the same color.
        dos_linestyle (list): 1\*nProj or nProj\*2 (spin) linestyle string. If
            spin>1 and 1\*nProj list is used, they are in the same style.
        dos_linewidth (list): 1\*nProj or nProj\*2 (spin) width of the plot
            lines. If spin>1 and 1\*nSystem list is used, they are in the same
            width.
        fermi (float|None): Fermi energy in eV. By default the band is aligned
            to 0. Can be used to offset the band. None for not plotting Fermi
            level.
        fermi_color (str|None): Color of the Fermi level.
        fermi_linestyle (str|None): Line style of Fermi level.
        fermi_linewidth(float|None): Width of the Fermi level.
        \*\*kwargs: Other commands passed to matplotlib ``axes.plot()`` method
            when plotting bands. Applied to all bands.

    Returns:
        ax (Axes): Matplotlib axes object
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import warnings, copy
    import matplotlib.colors as mcolors

    # Sanity check
    ## projection
    if len(prj) != 0:
        if len(prj) > len(doss):
            raise ValueError('Specified number of projects are lager than the length of data.')
        dossplt = doss[np.array(prj, dtype=int)-1]
    else:
        dossplt = copy.deepcopy(doss)

    nprj, nenergy,  nspin = dossplt.shape

    ## energy
    if len(energy) != nenergy:
        raise ValueError('Inconsistent DOSS and energy dimensions.')

    ## beta
    if beta.lower() == 'up':
        pass
    elif beta.lower() == 'down':
        if nspin == 2:
            dossplt[:, :, 1] = -dossplt[:, :, 1]
    else:
        raise ValueError("Beta can be defined only as: 'up' or 'down'.")

    # prepare energy and dos ranges
    if len(energy_range) == 0:
        energy_range = np.array([np.min(energy), np.max(energy)])
    else:
        energy_range = np.array([np.min(energy_range), np.max(energy_range)])
    if len(dos_range) == 0:
        dos_range = np.array([np.min(doss), np.max(doss)])
    else:
        dos_range = np.array([np.min(dos_range), np.max(dos_range)])

    # Prepare line label, color style and width
    ## label
    if np.all(dos_label!=None):
        if isinstance(dos_label, str):
            dos_label = [[dos_label, dos_label] for i in range(nprj)]
        else:
            if len(dos_label) != nprj:
                raise ValueError('Inconsistent system labels and number of systems.')
            for i in range(nprj):
                _, _, nspin = doss[i].shape
                ## label, and default setups of dos labels
                if not isinstance(dos_label[i], list):
                    if nspin == 2:
                        dos_label[i] = [r'{} ($\alpha$)'.format(dos_label[i]),
                                        r'{} ($\beta$)'.format(dos_label[i])]
                    else:
                        dos_label[i] = [dos_label[i], dos_label[i]]
                else:
                    dos_label = dos_label[i][0:2]

    else:
        dos_label = []; any_spin = False
        for i in range(nprj):
            _, _, nspin = doss[i].shape
            if nspin == 2:
                any_spin = True
                dos_label.append([r'$\alpha$', r'$\beta$'])
            else:
                dos_label.append(['', ''])

        if any_spin == False:
            dos_label = None
    ## color
    if np.all(dos_color!=None):
        if isinstance(dos_color, str):
            dos_color = [[dos_color, dos_color] for i in range(nprj)]
        else:
            if len(dos_color) != nprj:
                raise ValueError('Inconsistent dos colors and number of systems.')
            if not isinstance(dos_color[0], list):
                dos_color = [[i, i] for i in dos_color]
            else:
                dos_color = [[i[0], i[1]] for i in dos_color]
    else: # defalut setups of dos color
        clist = list(mcolors.TABLEAU_COLORS.keys())
        dos_color = [[clist[i], clist[i]] for i in range(nprj)]
    ## line style
    if np.all(dos_linestyle!=None):
        if isinstance(dos_linestyle, str):
            dos_linestyle = [[dos_linestyle, dos_linestyle] for i in range(nprj)]
        else:
            if len(dos_linestyle) != nprj:
                raise ValueError('Inconsistent dos line style and number of systems.')
            if not isinstance(dos_linestyle[0], list):
                dos_linestyle = [[i, i] for i in dos_linestyle]
            else:
                dos_linestyle = [[i[0], i[1]] for i in dos_linestyle]
    else: # defalut setups of line style
        dos_linestyle = [['-', '--'] for i in range(nprj)]
    ## linewidth
    if np.all(dos_linewidth!=None):
        if isinstance(dos_linewidth, int) or isinstance(dos_linewidth, float):
            dos_linewidth = [[dos_linewidth, dos_linewidth] for i in range(nprj)]
        else:
            if len(dos_linewidth) != nprj:
                raise ValueError('Inconsistent dos line width and number of systems.')
            if not isinstance(dos_linewidth[0], list):
                dos_linewidth = [[i, i] for i in dos_linewidth]
            else:
                dos_linewidth = [[i[0], i[1]] for i in dos_linewidth]
    else: # defalut setups of linewidth
        dos_linewidth = [[1.0, 1.0] for i in range(nprj)]

    commands = [dos_label, dos_color, dos_linestyle, dos_linewidth]
    keywords = ['label', 'color', 'linestyle', 'linewidth']

    # plot
    ## Fermi level
    if np.all(fermi!=None):
        energy = energy + fermi
        ax.vlines(fermi, dos_range[0], dos_range[1], color=fermi_color,
                  linestyle=dos_linestyle, linewidth=dos_linewidth)

    ## DOS=0 line
    if beta.lower() == 'down':
        ax.hlines(0, energy_range[0], energy_range[1], color='k', linewidth=0.5)

    ## DOS
    for iprj in range(nprj):
        for ispin in range(nspin):
            for icmd in range(4):
                if np.all(commands[icmd]==None): # 4*nprj*2(spin)
                    continue
                kwargs[keywords[icmd]] = commands[icmd][isys][ispin]

            ax.plot(energy, dosplt[iprj, :, ispin], **kwargs)

    # a label for a plot
    ax.legend(loc='lower right')
    ax.set_xlim(energy_range)
    ax.set_ylim(dos_range)
    return ax


def plot_cry_es(bands, doss, k_labels, color_bd, color_doss, fermi, energy_range, linestl_bd,
                linestl_doss, linewidth, prj, figsize, labels, dos_range, title, dos_beta):
    """
    The base function to plot electron / phonon band structure + DOS

    Args:
        bands (object): Object containing band structure data
        doss (object): Object containing density of states data
        k_labels (list or None): List of labels for high symmetry points along the path
        color_bd (str): Color for the band structure plot
        color_doss (str or list or tuple): Color(s) for the density of states plot
        fermi (str): Color for the Fermi level lines
        energy_range (list or None): Range of energy values for the y-axis
        linestl_bd (str): Linestyle for the band structure plot
        linestl_doss (str or list or tuple or None): Linestyle(s) for the density of states plot
        linewidth (float): Width of the lines
        prj (list or None): List of projection indices for plotting specific projections
        figsize (tuple): Figure size (width, height)
        labels (str or list or tuple or None): Labels for the density of states plot
        dos_range (list or None): Range of the density of states plot
        title (str or None): Title of the figure
        dos_beta (str): Beta state for the density of states plot ('up' or 'down')

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import warnings
    from os import path

    import matplotlib.pyplot as plt
    import numpy as np

    # Dictionary of greek letters for the band plot in the electronic structure
    greek = {'Alpha': '\u0391', 'Beta': '\u0392', 'Gamma': '\u0393', 'Delta': '\u0394', 'Epsilon': '\u0395', 'Zeta': '\u0396', 'Eta': '\u0397',
             'Theta': '\u0398', 'Iota': '\u0399', 'Kappa': '\u039A', 'Lambda': '\u039B', 'Mu': '\u039C', 'Nu': '\u039D', 'Csi': '\u039E',
             'Omicron': '\u039F', 'Pi': '\u03A0', 'Rho': '\u03A1', 'Sigma': '\u03A3', 'Tau': '\u03A4', 'Upsilon': '\u03A5', 'Phi': '\u03A6',
             'Chi': '\u03A7', 'Psi': '\u03A8', 'Omega': '\u03A9', 'Sigma_1': '\u03A3\u2081'}

    if isinstance(color_doss, list) or isinstance(color_doss, tuple):
        # Error check on linestl_doss length when the prj kwargs is not used
        if (prj is None) and (len(color_doss) != doss.n_proj):
            if len(color_doss) > doss.n_proj:
                warnings.warn('You have a number of linestl_doss elements is greater than the number of projection!',
                              stacklevel=2)
            else:
                raise ValueError(
                    "You don't have enough elements in linestl_doss for the number of projection required")
        # Error check on linestl_doss length when the prj kwargs is used
        elif (prj is not None) and (len(color_doss) != len(prj)):
            if len(color_doss) > len(prj):
                warnings.warn('You have a number of linestl_doss element greater than the number of projection required(prj elements)!',
                              stacklevel=2)
            else:
                raise ValueError(
                    "You don't have enough elements in linestl_doss for the number of projection required(prj elements)")
    else:
        if np.all(prj==None):
            nprj = doss.n_proj
        else:
            nprj = len(prj)
        if nprj > 1:
            warnings.warn(
                'Only one color / no color is given. All DOS projections are in the same color.')
        color_doss = [color_doss for i in range(nprj)]

    if isinstance(labels, list) or isinstance(labels, tuple):
        # Error check on labels length when the prj kwargs is not used
        if (prj is None) and (len(labels) != doss.n_proj):
            if len(labels) > doss.n_proj:
                warnings.warn(
                    'You have a number of labels greater than the number of projection!')
            else:
                raise ValueError(
                    "You don't have enough elements in labels for the numeber of projection required")
        # Error check on labels length when the prj kwargs is used
        elif (prj is not None) and (len(labels) != len(prj)):
            if len(labels) > len(prj):
                warnings.warn('You have a number of linestl_doss element greater than the number of projection required(prj elements)!',
                              stacklevel=2)
            else:
                raise ValueError(
                    "You don't have enough elements in linestl_doss for the number of projection required(prj elements)")
    else:
        if np.all(prj==None):
            nprj = doss.n_proj
        else:
            nprj = len(prj)
        if nprj > 1:
            warnings.warn('Label not given. No label is available for DOS.')
        labels = [None for i in range(nprj)]

    # DOS beta state
    if dos_beta == 'up':
        spin_idx = -1  # DOS in crystal output is distinguished by +/- sign already
        line_0 = False
    elif dos_beta == 'down':
        spin_idx = 1
        line_0 = True
    else:
        raise ValueError("'dos_beta' should be either 'up' or 'down'.")
    if doss.doss.shape[2] == 2:
        doss.doss[:, :, 1] = doss.doss[:, :, 1] * spin_idx

    # Definition and creation of the figure and the axes
    fig, ax = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [2, 1]},
                           sharex=False, sharey=True, figsize=figsize)
    if np.all(title!=None):
        fig.suptitle(title)

    # Definition of the hsp position variables
    hsp = bands.tick_pos

    # Error check on k_labels lenght against the HSP poisitions
    if k_labels is not None:
        if len(hsp) != len(k_labels):
            if len(hsp) > len(k_labels):
                raise ValueError(
                    'You have specified a number of label smaller than the number of High Simmetry Point along the path')
            elif len(hsp) < len(k_labels):
                raise ValueError(
                    'You have more labels than the High Simmetry point along the path')

    # Local variable definition for the band plot
    dx_bd = bands.k_path
    pltband = bands.bands
    no_bands = np.shape(pltband)[0]
    ymin_bd = np.amin(pltband)
    ymax_bd = np.amax(pltband)
    xmin_bd = np.amin(dx_bd)
    xmax_bd = np.amax(dx_bd)
    count1 = 0
    count2 = 0

    # band plot
    for i in range(no_bands):
        if bands.spin == 1:
            ax[0].plot(dx_bd, pltband[i, :], color=color_bd,
                       linestyle=linestl_bd, linewidth=linewidth)

        elif bands.spin == 2:
            if count1 == count2:
                ax[0].plot(dx_bd, pltband[i, :, 0], color=color_bd,
                           linestyle='-', linewidth=linewidth, label='Alpha')
                ax[0].plot(dx_bd, pltband[i, :, 1], color=color_bd,
                           linestyle='--', linewidth=linewidth, label='Beta')
            else:
                ax[0].plot(dx_bd, pltband[i, :, 0], color=color_bd,
                           linestyle='-', linewidth=linewidth)
                ax[0].plot(dx_bd, pltband[i, :, 1], color=color_bd,
                           linestyle='--', linewidth=linewidth)

        count1 += 1

    # Definition of dx for the doss plot
    dx_dos = doss.energy

    # Determination of xmin, xmax, ymin, and ymax
    if np.all(prj!=None):
        argplt = prj
    else:
        argplt = range(1, doss.doss.shape[1])

    # Plot a vertical line at 0 DOS. Only for spin polarized cases
    if line_0 == True and doss.spin == 2:
        xmin_dos = np.amin(doss.doss[:, argplt, 1])
        xmax_dos = np.amax(doss.doss[:, argplt, 0])
        # make the scale symmetric, for comparison
        xmin_dos = -max([abs(xmin_dos), abs(xmax_dos)])
        xmax_dos = -xmin_dos
    else:
        xmin_dos = 0.
        xmax_dos = np.amax(doss.doss[:, argplt, :])

    ymin_dos = np.amin(dx_dos)
    ymax_dos = np.amax(dx_dos)

    # Plot of all projections
    if prj is None:
        for projection in range(doss.n_proj):
            if doss.spin == 1:
                if doss.n_proj > 1:
                    if linestl_doss is None:
                        ax[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[projection],
                                   label=labels[projection], linewidth=linewidth)

                    else:
                        ax[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[projection],
                                   label=labels[projection], linestyle=linestl_doss[projection], linewidth=linewidth)

                else:
                    ax[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[0],
                               linewidth=linewidth)

            elif doss.spin == 2:
                if doss.n_proj > 1:
                    ax[1].plot(doss.doss[:, projection, 0], dx_dos, color=color_doss[projection],
                               label=labels[projection], linestyle='-', linewidth=linewidth)
                    ax[1].plot(doss.doss[:, projection, 1], dx_dos, color=color_doss[projection],
                               label=labels[projection], linestyle='--', linewidth=linewidth)
                else:
                    ax[1].plot(doss.doss[:, projection, 0], dx_dos, color=color_doss[0],
                               linestyle='-', linewidth=linewidth)
                    ax[1].plot(doss.doss[:, projection, 1], dx_dos, color=color_doss[0],
                               linestyle='--', linewidth=linewidth)

    # Plot of a selected number of projections
    else:
        for index, projection in enumerate(prj):
            projection = projection-1
            if doss.spin == 1:
                if doss.n_proj > 1:
                    if linestl_doss is None:
                        ax[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[index],
                                   label=labels[index], linewidth=linewidth)
                    else:
                        ax[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[index],
                                   label=labels[index], linestyle=linestl_doss[index], linewidth=linewidth)
                else:
                    ax[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[0],
                               linewidth=linewidth)
            elif doss.spin == 2:
                if doss.n_proj > 1:
                    ax[1].plot(doss.doss[:, projection, 0], dx_dos, color=color_doss[index],
                               label=labels[index], linestyle='-', linewidth=linewidth)
                    ax[1].plot(doss.doss[:, projection, 1], dx_dos, color=color_doss[index],
                               label=labels[index], linestyle='--', linewidth=linewidth)
                else:
                    ax[1].plot(doss.doss[:, projection, 0], dx_dos, color=color_doss[0],
                               linestyle='-', linewidth=linewidth)
                    ax[1].plot(doss.doss[:, projection, 1], dx_dos, color=color_doss[0],
                               linestyle='--', linewidth=linewidth)

    # Set Y axis, HSP lines (band) and 0 Lines (DOS)
    if energy_range is not None:
        ymin = energy_range[0]
        ymax = energy_range[1]
    else:
        ymin = min([ymin_bd, ymin_dos])
        ymax = max([ymax_bd, ymax_dos])
    ax[0].vlines(hsp, ymin, ymax, color='black', linewidth=0.5)
    ax[0].set_ylim(ymin, ymax)
    hsp_label = []
    if k_labels is not None:
        for n in k_labels:
            if n in greek:
                g = greek.get(n)
                hsp_label.append(g)
            else:
                hsp_label.append(n)
    ax[0].set_xticks(hsp)
    if k_labels is not None:
        ax[0].set_xticklabels(hsp_label)

    if line_0 == True:
        ax[1].vlines(0., ymin, ymax, color='black', linewidth=0.5)
    ax[1].set_ylim(ymin, ymax)

    # Set X axis and fermi lines
    xmax_bd = hsp[len(hsp)-1]
    ax[0].hlines(0., xmin_bd, xmax_bd, color=fermi, linewidth=1.5)
    ax[0].set_xlim(xmin_bd, xmax_bd)

    if dos_range is not None:
        xmin_dos = dos_range[0]
        xmax_dos = dos_range[1]

    # if (prj is None) and (doss.n_proj not in prj):
    #    xmax_dos = np.amax(doss.doss[:, 1:doss.n_proj-1, :])

    ax[1].hlines(0., xmin_dos*1.05, xmax_dos*1.05, color=fermi, linewidth=1.5)
    ax[1].set_xlim(xmin_dos*1.05, xmax_dos*1.05)

    if np.all(labels!=None):
        ax[1].legend()

    return fig, ax


def plot_2Dscalar(datamap, gridv, levels, xticks, yticks, cmap_max, cmap_min, cbar_label):
    """
    Plot 2D scalar field map.

    Args:
        datamap (array): 2D map data.
        gridv (array): 2\*3 base vectors of 2D map.
        levels (int | array-like): Determines the number and positions of the contour lines/regions.
        xticks (int): Number of ticks in the x direction.
        yticks (int): Number of ticks in the y direction.
        cmap_max (float): Maximum value used for the colormap.
        cmap_min (float): Minimun value used for the colormap.
        cbar_label (str): Title of colorbar (typically for quantuity and unit)

    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    vector_ab = gridv[0, :]
    length_ab = np.norm(v1)
    vector_cb = gridv[1, :]
    length_cb = np.norm(v2)
    points_ab = datamap.shape[0]
    points_cb = datamap.shape[1]
    cosxy = np.dot(vector_ab, vector_cb) / \
        np.norm(vector_ab) / np.norm(vector_cb)

    mesh_x = np.zeros((points_ab, points_cb), dtype=float)
    mesh_y = np.zeros((points_ab, points_cb), dtype=float)
    for i in range(0, points_ab):
        for j in range(0, points_cb):
            mesh_y[i, j] = (lenght_ab / points_ab) * i * np.sqrt(1 - cosxy**2)
            mesh_x[i, j] = (lenght_cb / points_cb) * j + \
                (lenght_ab / points_ab) * i * cosxy

    if cmap_max is None:
        max_data = np.amax(datamap)
    else:
        max_data = cmap_max

    if cmap_min is None:
        min_data = np.amin(datamap)
    else:
        min_data = cmap_min

    fig, ax = plt.subplots()
    im = ax.contourf(mesh_x, mesh_y, dens, levels, cmap='gnuplot')
    divider = make_axes_locatable(ax)
    im.set_clim(vmin=min_data, vmax=max_data)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbar_label, rotation=270)
    ax.set_xlabel('$\AA$')
    ax.set_xticks(np.linspace(0, lenght_cb, xticks).tolist())
    ax.set_yticks(np.linspace(0, lenght_ab, yticks).tolist())
    ax.set_ylabel('$\AA$')
    ax.set_aspect(1.0)
    ax.set_xlim(np.amin(mesh_x), np.amax(mesh_x))
    ax.set_ylim(0, np.amax(mesh_y) * np.sqrt(1 - obj_echg.cosxy**2))

    return fig, ax
