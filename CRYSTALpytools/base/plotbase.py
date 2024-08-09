#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base functions for plotting 2D and 3D figures
"""

def plot_overlap_bands(ax, bands, k_xax, k_path, k_label, energy_range, k_range,
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
        k_xax (list): Coordinates of x axis. 1\*nSystem list.
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

    if len(k_label) != 0:
        k_label = np.array(k_label, ndmin=1)
        if k_label.ndim > 1:
            raise ValueError("K labels must be string or float number and only 1 set of k labels can be set for 'multi' mode.")
        k_label = k_label.tolist()

    # preprocessing, always scale k path
    k_path, k_label, energy_range, k_range, commands = _plot_bands_preprocess(
        bands, k_path, k_label, False, energy_range, k_range,
        band_label, band_color, band_linestyle, band_linewidth
    )

    # Start plotting
    ## Fermi level
    ## Fermi check, must be None, float, int
    if np.all(fermi!=None):
        fermi = np.array(fermi, ndmin=1, dtype=float)[0]
        ax.hlines(fermi, k_range[0], k_range[1], color=fermi_color,
                  linestyle=fermi_linestyle, linewidth=fermi_linewidth)

    ## high symmetry lines
    ax.vlines(k_path[0], energy_range[0], energy_range[1], color='k', linewidth=0.5)

    ## bands
    keywords = ['label', 'color', 'linestyle', 'linewidth']
    ilabel = []; countlabel = 0
    idx = np.argmax([i[-1] for i in k_xax])
    k_xax_max = k_xax[idx]
    for isys in range(nsys):
        bandsplt = copy.deepcopy(bands[isys])
        if np.all(fermi!=None):
            bandsplt = bandsplt + fermi
            energy_range = energy_range + fermi

        nband, nkpt, nspin = bandsplt.shape
        k_pathplt = k_xax[isys] / k_xax[isys][-1] * k_xax_max[-1]
        # k_pathplt = np.linspace(np.min(k_path[isys]), np.max(k_path[isys]), nkpt)
        for ispin in range(nspin):
            for icmd in range(4):
                kwargs[keywords[icmd]] = commands[icmd][isys][ispin]
            ax.plot(k_pathplt, bandsplt[:, :, ispin].transpose(), **kwargs)
            # a label for a set of bands, dimension of bandsplt array might vary
            countlabel = countlabel + nband
        ilabel.append(countlabel-1)

    # a label for a set of bands, dimension of bandsplt array might vary
    if np.all(commands[0][0][0]!=None) and np.all(legend!=None):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend([handles[i] for i in ilabel],
                  [labels[i] for i in ilabel],
                  loc=legend)

    ax.set_xticks(k_path[0], labels=k_label[0])
    ax.set_xlim(k_range)
    ax.set_ylim(energy_range)
    return ax


def plot_compare_bands(
    ax, bands, k_xax, k_path, k_label, not_scaled, energy_range, k_range,
    band_label, band_color, band_linestyle, band_linewidth, fermi, fermi_color,
    fermi_linestyle, fermi_linewidth, legend, **kwargs):
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
        k_xax (list): Coordinates of x axis. 1\*nSystem list.
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
    ## label, same as overlap
    ## color
    if np.all(band_color==None):
        band_color = [['tab:blue', 'tab:blue'] for i in range(nsys)]
    else: # New Default
        band_color = np.array(band_color, ndmin=1)
        if band_color.shape[0] > 1:
            band_color = [[band_color[0], band_color[1]] for i in range(nsys)]
        else:
            band_color = str(band_color[0])
    ## line style
    if np.all(band_linestyle!=None):
        band_linestyle = np.array(band_linestyle, ndmin=1)
        if band_linestyle.shape[0] > 1:
            band_linestyle = [[band_linestyle[0], band_linestyle[1]] for i in range(nsys)]
        else:
            band_linestyle = str(band_linestyle[0])
    ## line width
    if np.all(band_linewidth!=None):
        band_linewidth = np.array(band_linewidth, ndmin=1)
        if band_linewidth.shape[0] > 1:
            band_linewidth = [[band_linewidth[0], band_linewidth[1]] for i in range(nsys)]
        else:
            band_linewidth = float(band_linewidth[0])

    # prepare fermi level
    if np.all(fermi!=None):
        fermi = np.array(fermi, ndmin=1, dtype=float)
        if len(fermi) == 1:
            fermi = fermi.repeat(nsys)
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

    # uniform x scale along the longest x axis
    if not_scaled != True:
        idx = np.argmax([i[-1] for i in k_xax])
        k_xax_max = k_xax[idx]
        for i in range(len(k_xax)):
            k_xax[i] = k_xax[i] / k_xax[i][-1] * k_xax_max[-1]

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
        k_pathplt = k_xax[isys]
        # k_pathplt = np.linspace(np.min(k_path[isys]), np.max(k_path[isys]), nkpt)
        for ispin in range(nspin):
            for icmd in range(4): # 4*nsys*2(spin)
                kwargs[keywords[icmd]] = commands[icmd][isys][ispin]
            ax[isys].plot(k_pathplt, bandsplt[:, :, ispin].transpose(), **kwargs)

        # a label for a set of bands
        if np.all(commands[0][0][0]!=None) and np.all(legend!=None):
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
    of parameters, refer to ``plot_overlap_bands()`` (``plot_compare_bands`` has
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
        k_label = np.array(k_label)
        if k_label.ndim == 1: # same definition for all k pathes
            same_klabel = True
            for i in range(nsys):
                if len(k_label) != len(k_path[i]):
                    raise ValueError('Inconsistent dimensions of k label and k path.')
            k_label = [k_label.tolist() for i in range(nsys)]
        else:
            same_klabel = False
            for i in range(nsys):
                if len(k_label[i]) != len(k_path[i]):
                    raise ValueError('Inconsistent dimensions of k label and k path.')
            k_label = k_label.tolist()

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

    1. For None input, generate nsystem\*2 list default values  
    2. For string/float inputs, generate nsystem\*2 list same values  
    3. For nsystem\*1 inputs, generate nsystem\*2 list, spin-up and down states share the same style.  
    4. For nsystem\*2 inputs, do nothing.

    return to a bigger command variable of \[label, color, linestyle, linewidth\].
    """
    import numpy as np
    import matplotlib.colors as mcolors

    nsys = len(bands)
    ## label
    if np.all(band_label!=None):
        if isinstance(band_label, str):
            band_label = [[band_label, band_label] for i in range(nsys)]
        else:
            band_label_ref = np.array(band_label)
            if band_label_ref.shape[0] != nsys:
                raise ValueError('Inconsistent system labels and number of systems(band) / projections(DOS).')
            for i in range(nsys):
                nspin = bands[i].shape[-1]
                if band_label_ref.ndim == 1:
                    if nspin == 2:
                        band_label[i] = [r'{} ($\alpha$)'.format(band_label_ref[i]),
                                         r'{} ($\beta$)'.format(band_label_ref[i])]
                    else:
                        band_label[i] = [band_label[i], band_label[i]]
                else:
                    band_label[i] = band_label[i][0:2]

    else: # defalut setups of band label
        band_label = []; any_spin = False
        for i in range(nsys):
            nspin = bands[i].shape[-1]
            if nspin == 2:
                any_spin = True
                band_label.append([r'$\alpha$', r'$\beta$'])
            else:
                band_label.append(['', ''])
        if any_spin == False:
            band_label = [[None, None] for i in range(nsys)]
    ## color
    if np.all(band_color!=None):
        if isinstance(band_color, str):
            band_color = [[band_color, band_color] for i in range(nsys)]
        else:
            band_color = np.array(band_color)
            if band_color.shape[0] != nsys:
                raise ValueError('Inconsistent band colors and number of systems(band) / projections(DOS).')
            if band_color.ndim == 1:
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
            band_linestyle = np.array(band_linestyle)
            if band_linestyle.shape[0] != nsys:
                raise ValueError('Inconsistent band line style and number of systems(band) / projections(DOS).')
            if band_linestyle.ndim == 1:
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
            band_linewidth = np.array(band_linewidth)
            if band_linewidth.shape[0] != nsys:
                raise ValueError('Inconsistent band line width and number of systems(band) / projections(DOS).')
            if band_linewidth.ndim == 1:
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
            for icmd in range(4): # 4*nprj*2(spin)
                kwargs[keywords[icmd]] = commands[icmd][iprj][ispin]
            if plot_vertical == False:
                ax.plot(energy, dossplt[iprj, :, ispin], **kwargs)
            else:
                ax.plot(dossplt[iprj, :, ispin], energy, **kwargs)

    # a label for a plot
    if np.all(commands[0][0][0]!=None) and np.all(legend!=None):
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
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from CRYSTALpytools.electronics import ElectronDOS
    from CRYSTALpytools.phonons import PhononDOS

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
        [ax[0]], [bands.bands], [bands.k_path], [bands.tick_pos], k_label,
        False, energy_range, k_range, band_label, band_color, band_linestyle,
        band_linewidth, fermi, fermi_color, fermi_linestyle, fermi_linewidth,
        legend, **kwargs
    )

    # plot DOS
    if isinstance(doss, ElectronDOS):
        xaxis = doss.energy
    else:
        xaxis = doss.frequency

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
            commands[i] = [[commands[i][j]] for j in range(ncol-1)]
        # subplots
        for i in range(ncol-1):
            ax.flat[i+1] = plot_doss(
                ax.flat[i+1], doss.doss, xaxis, beta, [prj[i]], energy_range,
                dos_range, commands[0][i], commands[1][i], commands[2][i], commands[3][i],
                fermi, fermi_color, fermi_linestyle, fermi_linewidth, legend, True, **kwargs
            )
    else:
        ax[1] = plot_doss(
            ax[1], doss.doss, xaxis, beta, prj, energy_range, dos_range,
            dos_label, dos_color, dos_linestyle, dos_linewidth,
            fermi, fermi_color, fermi_linestyle, fermi_linewidth, legend, True, **kwargs
        )

    return fig


def plot_2Dscalar(fig, ax, data, base, levels, contourline, isovalue, colormap, cbar_label,
                  a_range, b_range, rectangle, edgeplot, xticks, yticks, **kwargs):
    """
    Plot 2D scalar field map.

    Args:
        fig (Figure): Matplotlib Figure object
        ax (Axes): Matplotlib Axes object
        data (array): 2D map data.
        base (array): 3\*3 Cartesian coordinates of points A, B, C to define a
            2D map. Vectors BA and BC are used.
        levels (array|None): Contour line / color isovalues. It also defines
            the range of data.
        contourline (list|None): If not None, set line styles and colors of
            every contourline. nLevel\*3 list of matplotlib plot color,
            linestyle and linewidth.
        isovalue (str|None): If not None, set the format of isovalues added to
            contourlines. Useful only when ``contourline`` is not None.
        colormap (str|None): If not None, set the colormap of color-filled
            contour plots.
        cbar_label (str): Title of colorbar. Useful only when ``colormap`` is
            not None.
        a_range (list): Range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): Range of :math:`b` axis (x, or AB) in fractional coordinate.
        rectangle (bool): If :math:`a, b` are non-orthogonal, plot a rectangle
            region and reset :math:`b`. If used together with ``b_range``, that
            refers to the old :math:`b`.
        edgeplot (bool): Whether to plot plane edges
        xticks (int): Number of ticks in the x direction.
        yticks (int): Number of ticks in the y direction.
        \*\*kwargs: Other arguments passed to ``axes.contour()`` function to
            set contour lines.

    Returns:
        fig (Figure): Matplotlib Figure object
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from pymatgen.core.lattice import Lattice
    import copy

    # expand and rectangle
    X, Y, data, a_range, b_range = _manipulate_2D_grid(data, base, a_range, b_range, rectangle)

    # plot, put colormap at the back
    if np.all(colormap!=None):
        ax.contourf(X, Y, data, levels, cmap=colormap, vmin=np.min(levels), vmax=np.max(levels))
        norm = colors.Normalize(vmin=np.min(levels), vmax=np.max(levels), clip=False)
        m = cm.ScalarMappable(cmap=colormap, norm=norm)
        m.set_array(levels)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad="5%")
        colorbar = fig.colorbar(m, cax=cax)
        if np.all(cbar_label!=None):
            colorbar.set_label(cbar_label, rotation=270, labelpad=15)

    if np.all(contourline!=None):
        if len(contourline) != len(levels):
            raise ValueError('Inconsistent lengthes of contour line and contour line styles')
        clist = []; stlist = []; wlist = []
        for i in contourline:
            clist.append(i[0]); stlist.append(i[1]); wlist.append(i[2])

        L = ax.contour(X, Y, data, levels, colors=clist, linestyles=stlist,
                       linewidths=wlist, **kwargs)
        if np.all(isovalue!=None):
            ax.clabel(L, inline=1, fmt=isovalue)

    # plot plane edges
    if edgeplot == True:
        ## get shift: always close to the positive side of the plot
        xpath, ypath = _get_2D_base_frame(base, a_range, b_range)
        ax.plot(xpath, ypath,'k-', linewidth=1.0)

    # New ranges due to changes of a b ranges in non-orthogonal axis
    xrange = [np.round(np.min(X), 2), np.round(np.max(X), 2)]
    yrange = [np.round(np.min(Y), 2), np.round(np.max(Y), 2)]
    ax.set_xticks(np.round(np.linspace(xrange[0], xrange[1], xticks), 2))
    ax.set_yticks(np.round(np.linspace(yrange[0], yrange[1], yticks), 2))
    ax.set_aspect(1.0)
    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])
    return fig


def plot_2Dvector(fig, ax, data, base, scale, colorquiver, levels, colormap, cbar_label,
                  a_range, b_range, rectangle, edgeplot, xticks, yticks, **kwargs):
    """
    Plot 2D vector field map.

    Args:
        fig (Figure): Matplotlib Figure object
        ax (Axes): Matplotlib Axes object
        data (array): 2D vector map data, in nY\*nX\*3.
        base (array): 3\*3 Cartesian coordinates of points A, B, C to define a
            2D map. Vectors BA and BC are used.
        scale (float): Tune the length of arrows.
        colorquiver (str): Specify the color of arrows or 'colored' for color-
            coded quiver plots.
        levels (array): Contour color isovalues. It also defines the range of
            data. Useful only if ``colorquiver='colored'``.
        colormap (str|None): Set the colormap of color-filled contour plots.
            Useful only if ``colorquiver='colored'``.
        cbar_label (str): Title of colorbar. Useful only if
            ``colorquiver='colored'``.
        a_range (list): Range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): Range of :math:`b` axis (x, or AB) in fractional coordinate.
        rectangle (bool): If :math:`a, b` are non-orthogonal, plot a rectangle
            region and reset :math:`b`. If used together with ``b_range``, that
            refers to the old :math:`b`.
        edgeplot (bool): Whether to plot plane edges
        xticks (int): Number of ticks in the x direction.
        yticks (int): Number of ticks in the y direction.
        \*\*kwargs: Other arguments passed to ``axes.quiver()`` function to
            set contour lines.

    Returns:
        fig (Figure): Matplotlib Figure object
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from pymatgen.core.lattice import Lattice
    import copy

    # expand and rectangle
    X, Y, data, a_range, b_range = _manipulate_2D_grid(data, base, a_range, b_range, rectangle)

    # get projection and norm of the arrows
    rot, _ = _get_operation(base)
    x3d = rot.inv().apply([1, 0, 0]) # x axis in 3D reference framework
    y3d = rot.inv().apply([0, 1, 0]) # y axis in 3D reference framework
    dataprj = np.dstack([data@x3d, data@y3d])
    ## get norm
    vnorm = np.linalg.norm(data, axis=2)
    del data

    # plot
    if colorquiver == 'colored': # plot colored arrows
        norm = colors.Normalize(vmin=np.min(levels), vmax=np.max(levels), clip=False)
        cmap = plt.get_cmap(colormap)
        ax.quiver(X, Y, dataprj[:,:,0], dataprj[:,:,1], color=cmap(norm(vnorm.flatten())))
        m = cm.ScalarMappable(cmap=colormap, norm=norm)
        m.set_array(levels)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad="5%")
        colorbar = fig.colorbar(m, cax=cax)
        if np.all(cbar_label!=None):
            colorbar.set_label(cbar_label, rotation=270, labelpad=15)
    else: # plot same-color arrows
        ax.quiver(X, Y, dataprj[:,:,0], dataprj[:,:,1], color=colorquiver)

    # plot plane edges
    if edgeplot == True:
        ## get shift: always close to the positive side of the plot
        xpath, ypath = _get_2D_base_frame(base, a_range, b_range)
        ax.plot(xpath, ypath,'k-', linewidth=1.0)
    # New ranges due to changes of a b ranges in non-orthogonal axis
    xrange = [np.round(np.min(X), 2), np.round(np.max(X), 2)]
    yrange = [np.round(np.min(Y), 2), np.round(np.max(Y), 2)]
    ax.set_xticks(np.round(np.linspace(xrange[0], xrange[1], xticks), 2))
    ax.set_yticks(np.round(np.linspace(yrange[0], yrange[1], yticks), 2))
    ax.set_aspect(1.0)
    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])
    return fig


def _manipulate_2D_grid(data, base, a_range, b_range, rectangle):
    """
    Repeat 2D grid data and get rectangle region

    Args:
        data (array): 2D map data.
        base (array): 3\*3 Cartesian coordinates of points A, B, C to define a
            2D map. Vectors BA and BC are used.
        a_range (list): Range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): Range of :math:`b` axis (x, or AB) in fractional coordinate.
        rectangle (bool): If :math:`a, b` are non-orthogonal, plot a rectangle
            region and reset :math:`b`. If used together with ``b_range``, that
            refers to the old :math:`b`.

    Returns:
        X (array): nY\*nX mesh grid.
        Y (array): nY\*nX mesh grid.
        data (array): nY\*nX mesh grid data.
        a_range (array): 1\*2 fractional coordinate of plot range along vector BC.
        b_range (array): 1\*2 fractional coordinate of plot range along vector BA.
    """
    import numpy as np
    import copy

    vx = base[2, :] - base[1, :] # x, BC
    len_vx = np.linalg.norm(vx)
    npt_vx = data.shape[1] # a BA*BC matrix
    vy = base[0, :] - base[1, :] # y, AB
    len_vy = np.linalg.norm(vy)
    npt_vy = data.shape[0]
    cosxy = np.dot(vx, vy) / np.linalg.norm(vx) / np.linalg.norm(vy)
    sinxy = np.linalg.norm(np.cross(vx,vy)) / np.linalg.norm(vx) / np.linalg.norm(vy)

    X, Y = np.meshgrid(np.linspace(0, len_vx, npt_vx),
                       np.linspace(0, len_vy, npt_vy))
    # orthogonality
    X = np.round(cosxy, 3)*Y + X
    Y = np.round(sinxy, 3)*Y

    # periodicity
    ## a, b ranges are in lattice base vector. Can be non-orthogonal
    ulen_vx = len_vx / (npt_vx-1)
    ulen_vy = len_vy / (npt_vy-1)
    extend = False
    if len(a_range) == 2:
        a_range = np.array([np.min(a_range), np.max(a_range)]) * len_vx
        nlen_vx = a_range[1] - a_range[0]
        nnpt_vx = int(np.round(nlen_vx/ulen_vx, 0)) + 1
        extend = True
        if nlen_vx/len_vx - np.round(nlen_vx/len_vx,0) > 1e-3:
            raise ValueError('Use integer supercell sizes.')
    else:
        a_range = np.array([0., 1.]) * len_vx
        nlen_vx = len_vx
        nnpt_vx = npt_vx
    if len(b_range) == 2:
        b_range = np.array([np.min(b_range), np.max(b_range)]) * len_vy
        nlen_vy = b_range[1] - b_range[0]
        nnpt_vy = int(np.round(nlen_vy/ulen_vy, 0)) + 1
        extend = True
        if nlen_vy/len_vy - np.round(nlen_vy/len_vy,0) > 1e-3:
            raise ValueError('Use integer supercell sizes.')
    else:
        b_range = np.array([0., 1.]) * len_vy
        nlen_vy = len_vy
        nnpt_vy = npt_vy
    ## grid
    if extend == True:
        # update base vectors
        X, Y = np.meshgrid(np.linspace(0, nlen_vx, nnpt_vx),
                           np.linspace(0, nlen_vy, nnpt_vy))
        len_vx = nlen_vx
        len_vy = nlen_vy
        X = np.round(cosxy, 3)*Y + X
        Y = np.round(sinxy, 3)*Y

        # represent (x,y) in old (vx,vy) basis.
        shapedata = list(data.shape) # also usable for nY*nX*3 vector data
        shapedata[0] = nnpt_vy
        shapedata[1] = nnpt_vx
        ndata = np.zeros(shapedata, dtype=float)
        for i in range(nnpt_vy):
            for j in range(nnpt_vx):
                idx = j
                idy = i
                while idx >= npt_vx: idx -= npt_vx
                while idx < 0: idx += npt_vx
                while idy >= npt_vy: idy -= npt_vy
                while idy < 0: idy += npt_vy
                ndata[i, j] = data[idy, idx]
        del data
        data = copy.deepcopy(ndata)
        del ndata
        # shift data grid's origin to a_range[0], b_range[0]
        ncol = int(np.round(a_range[0]/ulen_vx, 0))
        nrow = int(np.round(b_range[0]/ulen_vy, 0))
        if ncol < 0:
            tmp = copy.deepcopy(data[:, nnpt_vx+ncol:])
            data[:, -ncol:] = data[:, 0:nnpt_vx+ncol]
            data[:, 0:-ncol] = tmp
        else:
            tmp = copy.deepcopy(data[:, 0:ncol])
            data[:, 0:nnpt_vx-ncol] = data[:, ncol:]
            data[:, nnpt_vx-ncol:] = tmp
        if nrow < 0:
            tmp = copy.deepcopy(data[nnpt_vy+nrow:, :])
            data[-nrow:, :] = data[0:nnpt_vy+nrow, :]
            data[0:-nrow, :] = tmp
        else:
            tmp = copy.deepcopy(data[0:nrow, :])
            data[0:nnpt_vy-nrow, :] = data[nrow:, :]
            data[nnpt_vy-nrow:, :] = tmp

    # get rectangle region
    if rectangle == True and np.abs(cosxy) > 1e-3:
        for i in range(nnpt_vy):
            if cosxy < 0: # case 1, cosxy < 0
                j = 0
                while X[i, j] < 0:
                    tmp = X[i, j]
                    X[i, j:-1] = X[i, j+1:]
                    X[i, -1] = tmp + len_vx
                    tmpd = data[i, j]
                    data[i, j:-1] = data[i, j+1:]
                    data[i, -1] = tmpd
            else: # case 2, cosxy > 0
                j = nnpt_vx-1
                while X[i, j] >= len_vx:
                    tmp = X[i, j]
                    X[i, 1:j+1] = X[i, 0:j]
                    X[i, 0] = tmp - len_vx
                    tmpd = data[i, j]
                    data[i, 1:j+1] = data[i, 0:j]
                    data[i, 0] = tmpd
        len_vy = len_vy * sinxy
        cosxy = 0.
        sinxy = 1.

    return X, Y, data, a_range, b_range


def _get_2D_base_frame(base, a_range, b_range):
    """
    Get the 2D parallelogram plot boundary. Useful when the plot is extended.
    The frame is always shifted to the origin, or positive side of the plot and
    close to the origin.

    Args:
        base (array): 3\*3 Cartesian coordinates of points A, B, C to define a
            2D map. Vectors BA and BC are used.
        a_range (list): Range of :math:`a` axis (x, or BC) in fractional coordinate.
        b_range (list): Range of :math:`b` axis (x, or AB) in fractional coordinate.

    Returns:
        xpath (array): 1\*3 array of x coordinates of parallelogram plot boundary.
        ypath (array): 1\*3 array of y coordinates of parallelogram plot boundary.
    """
    import numpy as np

    vx = base[2, :] - base[1, :]
    len_vx = np.linalg.norm(vx)
    vy = base[0, :] - base[1, :]
    len_vy = np.linalg.norm(vy)
    cosxy = np.dot(vx, vy) / np.linalg.norm(vx) / np.linalg.norm(vy)
    sinxy = np.linalg.norm(np.cross(vx,vy)) / np.linalg.norm(vx) / np.linalg.norm(vy)
    shiftx = (a_range[0]/len_vx - int(a_range[0]/len_vx)) * len_vx
    shifty = (b_range[0]/len_vy - int(b_range[0]/len_vy)) * len_vx
    if shiftx < 0: shiftx += len_vx
    if shifty < 0: shifty += len_vy
    shiftx = shiftx + shifty * cosxy
    shifty = shifty * sinxy
    xpath = np.array([0, len_vx, len_vx+len_vy*cosxy, len_vy*cosxy, 0]) + shiftx
    ypath = np.array([0, 0, len_vy*sinxy, len_vy*sinxy, 0]) + shifty
    return xpath, ypath


def _get_operation(base):
    """
    Get the rotation object and translational movement to align surface norm
    (to z) and BC axis (to x) of 3D reference frame to the plotting frame. The
    translational movement is used to move B (plot origin) to z=0. The plotting
    referance frame is defined by the base vector BC (x) and BA (y).

    Returns:
        rot (Rotation): The Scipy rotation object.
        disp (array): Displacement along x, y, z axes
    """
    from scipy.spatial.transform import Rotation
    import numpy as np

    pltx = base[2, :] - base[1, :]
    plty = base[0, :] - base[1, :]
    pltnorm = np.cross(pltx, plty)
    pltynorm = np.cross(pltnorm, pltx)# Y not necessarily orthogonal to xz plane

    pltx = pltx / np.linalg.norm(pltx)
    pltynorm = pltynorm / np.linalg.norm(pltynorm)
    pltnorm = pltnorm / np.linalg.norm(pltnorm)

    oldv = np.vstack([pltx, pltynorm, pltnorm]).transpose()
    newv = np.eye(3)
    rot = Rotation.from_matrix(newv @ np.linalg.inv(oldv))
    disp = -rot.apply(base[1, :])
    return rot, disp
