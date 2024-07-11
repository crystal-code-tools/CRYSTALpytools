#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base functions for plotting 2D and 3D figures
"""

def plot_overlap_bands(ax, bands, k_path, k_label, energy_range, k_range,
                       band_label, band_color, band_linestyle, band_linewidth,
                       fermi, fermi_color, fermi_linestyle, fermi_linewidth,
                       legend,**kwargs):
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
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
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
    ## Fermi check, must be None, float, int
    if np.all(fermi!=None):
        if not isinstance(fermi, float) and isinstance(fermi, int):
            raise ValueError('Fermi level must be None, float or int.')
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
    if np.all(commands[0]!=None) and np.all(legend!=None):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend([handles[i] for i in ilabel],
                  [labels[i] for i in ilabel],
                  loc=legend)

    ax.set_xticks(k_path[0], labels=k_label[0])
    ax.set_xlim(k_range)
    ax.set_ylim(energy_range)
    return ax


def plot_compare_bands(ax, bands, k_path, k_label, not_scaled, energy_range, k_range,
                       band_label, band_color, band_linestyle, band_linewidth,
                       fermi, fermi_color, fermi_linestyle, fermi_linewidth,
                       legend, **kwargs):
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
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
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
        if np.all(commands[0]!=None) and np.all(legend!=None):
            handles, labels = ax[isys].get_legend_handles_labels()
            ilabel = [int(i*nband) for i in range(nspin)]
            ax[isys].legend([handles[i] for i in ilabel],
                            [labels[i] for i in ilabel],
                            loc=legend)

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
        energy_range = np.array([np.min([np.min(i) for i in bands]),
                                 np.max([np.max(i) for i in bands])], dtype=float)

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

    # Get plot labels colors lines...
    commands = _plot_label_preprocess(bands, band_label, band_color, band_linestyle, band_linewidth)
    return k_path, k_label, energy_range, k_range, commands


def _plot_label_preprocess(bands, band_label, band_color, band_linestyle, band_linewidth):
    """
    Do the boring parameters checking jobs for plots (both band and dos). For
    the meanings of parameters, refer to ``plot_compare_bands``. The same rule
    is applied to DOS.
    """
    import numpy as np
    import matplotlib.colors as mcolors

    nsys = len(bands)
    if np.all(band_label!=None):
        if isinstance(band_label, str):
            band_label = [[band_label, band_label] for i in range(nsys)]
        else:
            if len(band_label) != nsys:
                raise ValueError('Inconsistent system labels and number of systems(band) / projections(DOS).')
            for i in range(nsys):
                nspin = bands[i].shape[-1]
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
            nspin = bands[i].shape[-1]
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
                raise ValueError('Inconsistent band colors and number of systems(band) / projections(DOS).')
            if not isinstance(band_color[0], list):
                band_color = [[i, i] for i in band_color]
            else:
                band_color = [[i[0], i[1]] for i in band_color]
    else: # defalut setups of band color
        clist = list(mcolors.TABLEAU_COLORS.keys())
        nclist = len(clist)
        band_color = [[clist[i%nclist], clist[i%nclist]] for i in range(nsys)]
    ## line style
    if np.all(band_linestyle!=None):
        if isinstance(band_linestyle, str):
            band_linestyle = [[band_linestyle, band_linestyle] for i in range(nsys)]
        else:
            if len(band_linestyle) != nsys:
                raise ValueError('Inconsistent band line style and number of systems(band) / projections(DOS).')
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
                raise ValueError('Inconsistent band line width and number of systems(band) / projections(DOS).')
            if not isinstance(band_linewidth[0], list):
                band_linewidth = [[i, i] for i in band_linewidth]
            else:
                band_linewidth = [[i[0], i[1]] for i in band_linewidth]
    else: # defalut setups of linewidth
        band_linewidth = [[1.0, 1.0] for i in range(nsys)]

    commands = [band_label, band_color, band_linestyle, band_linewidth] # ncmd\*nsys\*2(spin)
    return commands


def plot_doss(ax, doss, energy, beta, prj, energy_range, dos_range,
              dos_label, dos_color, dos_linestyle, dos_linewidth,
              fermi, fermi_color, fermi_linestyle, fermi_linewidth, legend,
              plot_vertical, **kwargs):
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
        fermi_linewidth (float|None): Width of the Fermi level.
        legend (str|None): Loc parameter passed to `axes.legend() <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html>`_
            None for not adding legend.
        plot_vertical (bool): *Developer Only* Get vertical (DOS-Energy) DOS
            plots.
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
        dos_range = np.array([np.min(dossplt), np.max(dossplt)])
    else:
        dos_range = np.array([np.min(dos_range), np.max(dos_range)])

    # Prepare line label, color style and width
    commands = _plot_label_preprocess(dossplt, dos_label, dos_color, dos_linestyle, dos_linewidth)
    keywords = ['label', 'color', 'linestyle', 'linewidth']

    # plot
    ## Fermi level
    ## Fermi check, must be None, float, int
    if np.all(fermi!=None):
        if not isinstance(fermi, float) and isinstance(fermi, int):
            raise ValueError('Fermi level must be None, float or int.')
        energy = energy + fermi
        if plot_vertical == False:
            ax.vlines(fermi, dos_range[0], dos_range[1], color=fermi_color,
                      linestyle=fermi_linestyle, linewidth=fermi_linewidth)
        else:
            ax.hlines(fermi, dos_range[0], dos_range[1], color=fermi_color,
                      linestyle=fermi_linestyle, linewidth=fermi_linewidth)

    ## DOS=0 line
    if beta.lower() == 'down':
        if plot_vertical == False:
            ax.hlines(0, energy_range[0], energy_range[1], color='k', linewidth=0.5)
        else:
            ax.vlines(0, energy_range[0], energy_range[1], color='k', linewidth=0.5)

    ## DOS
    for iprj in range(nprj):
        for ispin in range(nspin):
            for icmd in range(4):
                if np.all(commands[icmd]==None): # 4*nprj*2(spin)
                    continue
                kwargs[keywords[icmd]] = commands[icmd][iprj][ispin]

            if plot_vertical == False:
                ax.plot(energy, dossplt[iprj, :, ispin], **kwargs)
            else:
                ax.plot(dossplt[iprj, :, ispin], energy, **kwargs)

    # a label for a plot
    if np.all(commands[0]!=None) and np.all(legend!=None):
        ax.legend(loc=legend)

    if plot_vertical == False:
        ax.set_xlim(energy_range)
        ax.set_ylim(dos_range)
    else:
        ax.set_xlim(dos_range)
        ax.set_ylim(energy_range)
    return ax


def plot_banddos(bands, doss, k_label, beta, overlap, prj, energy_range, k_range,
                 dos_range, band_width, band_label, band_color, band_linestyle,
                 band_linewidth, dos_label, dos_color, dos_linestyle, dos_linewidth,
                 fermi, fermi_color, fermi_linestyle, fermi_linewidth, figsize,
                 legend, **kwargs):
    """
    The base function to plot electron / phonon band structure + DOS. A single
    system only.

    Input arguments not in the list are consistent with ``plot_doss`` and
    ``plot_compare_bands``.

    Args:
        bands (ElectronBand): A ``electronics.ElectronBand`` object.
        doss (ElectronDOS): A ``electronics.ElectronDOS`` object
        band_width (int|float): Relative width of band structure, times of the
            width of a DOS subplot.
        overlap (bool): Plot DOS projections into the same axes or multiple
            axes.
    Returns:
        fig (Figure): Matplotlib figure object
        ax (Axes): Matplotlib axes object
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Definition and creation of the figure and the axes
    if overlap == False:
        w_ratio = [band_width]
        w_ratio.extend([1 for i in range(len(prj))])
    else:
        w_ratio = [band_width, 1]

    ncol = len(w_ratio)
    fig, ax = plt.subplots(1, ncol, gridspec_kw={'width_ratios': w_ratio},
                           sharex=False, sharey=True, figsize=figsize, layout='constrained')

    # plot band structure
    ax[0] = plot_compare_bands(
        [ax[0]], [bands.bands], [bands.tick_pos], k_label, False, energy_range,
        k_range, band_label, band_color, band_linestyle, band_linewidth,
        fermi, fermi_color, fermi_linestyle, fermi_linewidth, legend, **kwargs
    )

    # plot DOS
    if overlap == False:
        # new defaults: all lines in the same color.
        if np.all(dos_color==None):
            dos_color = [['tab:blue', 'tab:blue'] for i in range(ncol-1)]
        # Dimeonsion issue: dos plot styles must be consistent with length of input dosss
        dossref = [doss.doss[i-1] for i in prj]
        commands = _plot_label_preprocess(
            dossref, dos_label, dos_color, dos_linestyle, dos_linewidth
        )
        for i in range(4):
            if np.all(commands[i]==None):
                commands[i] = [None for j in range(ncol-1)]
            else:
                commands[i] = [[commands[i][j]] for j in range(ncol-1)]
        # subplots
        for i in range(ncol-1):
            ax.flat[i+1] = plot_doss(
                ax.flat[i+1], doss.doss, doss.energy, beta, [prj[i]], energy_range,
                dos_range, commands[0][i], commands[1][i], commands[2][i], commands[3][i],
                fermi, fermi_color, fermi_linestyle, fermi_linewidth, legend, True, **kwargs
            )
    else:
        ax[1] = plot_doss(
            ax[1], doss.doss, doss.energy, beta, prj, energy_range, dos_range,
            dos_label, dos_color, dos_linestyle, dos_linewidth,
            fermi, fermi_color, fermi_linestyle, fermi_linewidth, legend, True, **kwargs
        )

    return fig, fig.axes


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
