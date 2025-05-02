#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic utility functions for plotting 2D and 3D figures.

.. note::

    3D plotting functions are based on `MayaVi <https://docs.enthought.com/mayavi/mayavi/>`_,
    which is not included in the default dependences of CRYSTALpytools.

"""

#--------------------------- 2D plots based on Matplotlib --------------------#

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


#---------------------------- grid general manipulation ----------------------#


def GridCoordinates(base, shape, meshgrid):
    """
    Get cartesian coordinates by equally spaced grid data.

    Args:
        base (array): 4(3)\*3 Cartesian coordinates of points O, A, B(, C) to
            define a 3(2)D map. Vectors OA, OB(, OC) are used.
        shape (array): 1\*nDimen array, nX, nY(, nZ) of data array
        meshgrid (bool): Get nD mesh grids or 1D array of coordinates.
    Returns:
        coords (list): 1\*3 list of x, y(, z) coordinates, either in nA\*nB(\*nC)
            mesh grids or in 1D arrays, depending on ``meshgrid``.
    """
    import numpy as np

    # sanity check
    base = np.array(base)
    shape = np.array(shape, ndmin=1)
    ndim = shape.shape[0]

    if ndim == 2:
        if base.shape[0] != 3 or base.shape[1] != 3:
            raise Exception("2D data grid must be defined on a 3*3 base vector.")
    elif ndim == 3:
        if base.shape[0] != 4 or base.shape[1] != 3:
            raise Exception("3D data grid must be defined on a 4*3 base vector.")
    else:
        raise Exception("Only for 2D/3D data grids.")
    basev = np.vstack([base[i]-base[0] for i in range(1, ndim+1)])

    if meshgrid == False:
        coords = []
        for i in range(ndim):
            frac = np.linspace(0, 1, shape[i], endpoint=False).reshape([-1, 1])
            coords.append(np.add(frac @ [basev[i]], base[0]))
    else:
        fcoords = []
        for i in range(ndim):
            fcoords.append(
                np.linspace(0, 1, shape[i], endpoint=False).reshape([-1, 1])
            )
        coords = np.meshgrid(*fcoords, indexing='ij'); del fcoords
        coords = np.vstack([i.flatten() for i in coords]).T
        coords = np.add((coords @ basev), base[0])  # n*3 cartesian from 0
        coords = coords.T.reshape([3]+[i for i in shape])
    return coords


def GridExpand(base, data, display_range):
    """
    Expand 2D/3D data grid.

    Args:
        base (array): 4(3)\*3 Cartesian coordinates of points O, A, B(, C) to
            define a 3(2)D map. Vectors OA, OB(, OC) are used.
        data (array): nX\*nY(\*nZ) array of data.
        display_range (array): 3(2)\*2 array of fractional display range, in
            \[amin, amax\].
    Returns:
        newbase (array): Expanded base vectors.
        newdata (array): Expanded data grid.
    """
    import numpy as np

    # sanity check
    base = np.array(base)
    data = np.array(data)
    display_range = np.array(display_range)

    if data.ndim == 2:
        if base.shape[0] != 3 or base.shape[1] != 3:
            raise Exception("2D data grid must be defined on a 3*3 base vector.")
    elif data.ndim == 3:
        if base.shape[0] != 4 or base.shape[1] != 3:
            raise Exception("3D data grid must be defined on a 4*3 base vector.")
    else:
        raise Exception("Only for 2D/3D data grids.")
    basev = np.vstack([base[i]-base[0] for i in range(1, data.ndim+1)])

    dispbg = np.round([i[0] for i in display_range], 12)
    disped = np.round([i[1] for i in display_range], 12)
    dist = disped - dispbg

    idx = np.where(dist<1e-12)[0]
    if len(idx) > 0:
        direct = ['a', 'b', 'c'][idx[0]]
        raise Exception("Display range error along {} axis!\n{} min = {:.2f}, {} max = {:.2f}. No data is displayed.".format(
            direct, direct, dispbg[idx[0]], direct, disped[idx[0]]))

    # new data and new base
    newdata = np.zeros(np.array(
        np.round(dist*data.shape, 0), dtype=int
    ))
    newbase = [dispbg @ basev + base[0]]
    for i in range(data.ndim):
        newbase = np.vstack([newbase, basev[i]*dist[i] + newbase[0]])

    # duplicate data
    origin = np.round(np.multiply(dispbg, data.shape), 0)
    origin = np.array(origin, dtype=int)
    end = np.round(np.multiply(disped, data.shape), 0)
    end = np.array(end, dtype=int)
    if data.ndim == 3:
        oldidx1 = np.arange(origin[0], end[0], 1)
        oldidx1 = np.sign(oldidx1) * np.abs(oldidx1) % data.shape[0]
        oldidx2 = np.arange(origin[1], end[1], 1)
        oldidx2 = np.sign(oldidx2) * np.abs(oldidx2) % data.shape[1]
        oldidx3 = np.arange(origin[2], end[2], 1)
        oldidx3 = np.sign(oldidx2) * np.abs(oldidx2) % data.shape[2]
        for i, oi in enumerate(oldidx1):
            for j, oj in enumerate(oldidx2):
                for k, ok in enumerate(oldidx2):
                    newdata[i, j, k] = data[oi, oj, ok]
    else:
        oldidx1 = np.arange(origin[0], end[0], 1)
        oldidx1 = np.sign(oldidx1) * np.abs(oldidx1) % data.shape[0]
        oldidx2 = np.arange(origin[1], end[1], 1)
        oldidx2 = np.sign(oldidx2) * np.abs(oldidx2) % data.shape[1]
        for i, oi in enumerate(oldidx1):
            for j, oj in enumerate(oldidx2):
                newdata[i, j] = data[oi, oj]
    return newbase, newdata


def GridInterpolate(base, data, method, size):
    """
    Interpolate 2D/ 3D grid data by `scipy.interpolate.interpn <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interpn.html>`_.

    Args:
        base (array): 3(4)\*3 Cartesian coordinates of the 3(4) points defining
            base vectors BA, BC (2D) or OA, OB, OC (3D). The sequence is (O),
            A, B, C.
        data (array): nX\*nY(\*nZ) data grid.
        method (str): 'linear', 'nearest', 'slinear', 'cubic'.
        size (int): The new size of interpolated data (list) or a scaling factor.
    Returns:
        DATA (array): nX\*nY(\*nZ) data grid.
        CRDS (array): (nX\*nY\*nZ)\*3 or (nX\*nY)\*3 cartesian coordinates.
            Sequence is consistent with the flattened ``DATA``.
    """
    import numpy as np
    from scipy.interpolate import griddata, interpn

    data = np.array(data, dtype=float)
    base = np.array(base, dtype=float)
    if data.ndim == 3:
        if base.shape[0] != 4 or base.shape[1] != 3:
            raise ValueError('4*3 array of point O, A, B, C must be defined.')
    elif data.ndim == 2:
        if base.shape[0] != 3 or base.shape[1] != 3:
            raise ValueError('3*3 array of point A, B, C must be defined.')
    else:
        raise ValueError('Only for 3D or 2D scalar fields.')

    # method
    method = method.lower()
    if method not in ['linear', 'nearest', 'slinear', 'cubic']:
        raise ValueError("Unknown interpolation method: {}.".format(method))

    # size
    size = np.array(size, ndmin=1, dtype=int)
    if len(size) == 1:
        size = size.repeat(data.ndim)
    if len(size) != data.ndim:
        raise ValueError("Specified dimension of interpolation = {:d}. The system dimension = {:d}.".format(len(size), data.ndim))
    if np.all(size<1):
        raise ValueError("'size' must have values larger than or equal to 1.")

    # Interpolate on coordinates
    newgrid_size = [size[i]*data.shape[i] for i in range(data.ndim)]
    if data.ndim == 3:
        oldgrid = (np.linspace(0, 1, data.shape[0]),
                   np.linspace(0, 1, data.shape[1]),
                   np.linspace(0, 1, data.shape[2]))
        newgrid = np.meshgrid(
            np.linspace(0, 1, newgrid_size[0]),
            np.linspace(0, 1, newgrid_size[1]),
            np.linspace(0, 1, newgrid_size[2]),
            indexing='ij')
    else:
        oldgrid = (np.linspace(0, 1, data.shape[0]),
                   np.linspace(0, 1, data.shape[1]))
        newgrid = np.meshgrid(
            np.linspace(0, 1, newgrid_size[0]),
            np.linspace(0, 1, newgrid_size[1]),
            indexing='ij')

    newgrid = np.array([i.flatten() for i in newgrid]).T
    DATA = interpn(oldgrid, data, newgrid, method=method)
    DATA = DATA.reshape(newgrid_size)
    CRDS = GridCoordinates(base, newgrid_size, meshgrid=True)
    if data.ndim == 3:
        CRDS = np.transpose(CRDS, axes=(1,2,3,0)).reshape([-1, 3])
    else:
        CRDS = np.transpose(CRDS, axes=(1,2,0)).reshape([-1, 3])
    del data

    return DATA, CRDS


def GridRectangle2D(base, data):
    """
    Get 2D grid data on a rectangle regular grid by manipulating the data
    matrix and orient it to xOy plane.

    Args:
        base (array): 3\*3 Cartesian coordinates of points O, A, B to define a
            2D map. Vectors OA, OB are used.
        data (array): Grid data in nX\*nY.
    Returns:
        basenew (array)
        data (array)
        X (array): nX\*nY mesh grid of x coordinates
        Y (array): nX\*nY mesh grid of y coordinates
    """
    import numpy as np
    import copy

    va = base[1] - base[0]
    lena = np.linalg.norm(va)
    vb = base[2] - base[0]
    lenb = np.linalg.norm(vb)
    theta = np.arccos(np.clip(np.dot(va/lena, vb/lenb), -1.0, 1.0))
    cosab = np.cos(theta)
    sinab = np.sin(theta)

    X, Y = np.meshgrid(np.linspace(0, lena, data.shape[0], endpoint=False),
                       np.linspace(0, lenb, data.shape[1], endpoint=False),
                       indexing='ij')
    X = np.round(cosab*Y + X, 12)
    Y = np.round(sinab*Y, 12)

    if abs(cosab) > 1e-4:
        if cosab > 0: # triangle from end to bg
            for j in range(data.shape[1]):
                idx = np.where(X[:, j]>=lena)[0]
                if idx.shape[0] < 1: continue
                i = np.min(idx)
                tmp = copy.deepcopy(data[i:, j])
                data[-i:, j] = data[:i, j]
                data[:-i, j] = tmp
                tmp = copy.deepcopy(X[i:, j])
                X[-i:, j] = X[:i, j]
                X[:-i, j] = tmp - lena
        else: # triangle from bg to end
            for j in range(data.shape[1]):
                idx = np.where(X[:, j]<0)[0]
                if idx.shape[0] < 1: continue
                i = np.max(idx)
                tmp = copy.deepcopy(data[:i+1, j])
                data[:-i-1, j] = data[i+1:, j]
                data[-i-1:, j] = tmp
                tmp = copy.deepcopy(X[:i+1, j])
                X[:-i-1, j] = X[i+1:, j]
                X[-i-1:, j] = tmp + lena

        # Bug fix: new X and Y mesh grid cannot be generated by basenew.
        # Otherwise tiny shifts might occur
        newb = np.cross(np.cross(va, vb)/sinab/lena/lenb, va/lena) * sinab * lenb + base[0]
        basenew = np.vstack([base[0], base[1], newb])
    else:
        basenew = base

    rot, _ = GridRotation2D(basenew)
    base2d = rot.apply(basenew)
    X += base2d[0, 0]
    Y += base2d[0, 1]
    return basenew, data, X, Y


def GridRotation2D(base):
    """
    Get the rotation object and translational movement to align surface norm
    (to z) and OA axis (to x) of 3D reference frame to the plotting frame. The
    translational movement is used to move O (plot origin) to z=0. The plotting
    referance frame is defined by the base vector OA (x) and OB (y).

    Returns:
        rot (Rotation): The Scipy rotation object.
        disp (array): Displacement along x, y, z axes
    """
    from scipy.spatial.transform import Rotation
    import numpy as np

    pltx = base[1] - base[0]
    plty = base[2] - base[0]
    pltnorm = np.cross(pltx, plty)
    pltynorm = np.cross(pltnorm, pltx)# Y not necessarily orthogonal to xz plane

    pltx = pltx / np.linalg.norm(pltx)
    pltynorm = pltynorm / np.linalg.norm(pltynorm)
    pltnorm = pltnorm / np.linalg.norm(pltnorm)

    oldv = np.vstack([pltx, pltynorm, pltnorm]).T
    newv = np.eye(3)
    rot = Rotation.from_matrix(newv @ np.linalg.inv(oldv))
    disp = -rot.apply(base[0])
    return rot, disp


def tvtkGrid(base, data, CenterOrigin, InterpGridSize, **kwargs):
    """
    .. _ref-tvtkGrid:

    Define a 3D/2D tvtk Grid for `MayaVi <https://docs.enthought.com/mayavi/mayavi/data.html>`_.
    Only for representing **3D/2D scalar field** defined in the periodic,
    uniform mesh.

    * For orthogonal data grid aligned to x, y and z axes, return to the
        ``ImageData`` class.  
    * For non-orthogonal data grid or not aligned data grid, return to the
        ``StructuredGrid`` classes, or interpolated ``ImageData`` class by
        the `scipy.interpolate.griddata() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html>`_
        method and linear interpolation. Might be **very time consuming**.

    .. note::

        For generality, the input data grid is a non-periodic one.

    Args:
        base (array): Base vectors, 3(4)\*3 array of origin and point A, B(, C)
        data (array): Scalar field data in column-major order, i.e., nA\*nB\*nC.
        CenterOrigin (bool): Put origin of base vectors in the center. Usually
            for reciprocal space visualization.
        InterpGridSize (array|int|None): Interpolate non-orthogonal data into
            orthogonal grid. 'None' for no interpolation. Integer input for
            interpolation sizes. For volume data representation.
        \*\*kwargs: Passed to `scipy.interpolate.griddata <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html>`_.
    Returns:
        grid (ImageData|StructuredGrid): vtk grid classes.
    """
    import numpy as np
    import copy
    from scipy.interpolate import griddata
    try:
        from tvtk.api import tvtk
    except ModuleNotFoundError:
        raise ModuleNotFoundError('MayaVi is required for this functionality, which is not in the default dependency list of CRYSTALpytools.')

    data = np.array(data, dtype=float)
    base = np.array(base, dtype=float)
    ndim = data.ndim
    if ndim == 3:
        if base.shape[0] != 4 or base.shape[1] != 3:
            raise ValueError('4*3 array of point O, A, B, C must be defined.')
    elif ndim == 2:
        if base.shape[0] != 3 or base.shape[1] != 3:
            raise ValueError('3*3 array of point O, A, B must be defined.')
        data = data.reshape([data.shape[0], data.shape[1], 1])
    else:
        raise ValueError('For 3D/2D data grid only.')

    # alignment
    bvec = np.array([base[i]-base[0] for i in range(1, ndim+1)])
    if ndim == 2: bvec = np.vstack([bvec, np.cross(bvec[0], bvec[1])])
    bvnorm = np.linalg.norm(bvec, axis=1)
    align = bvec @ np.eye(3)

    # Regular interpolation size
    if np.all(InterpGridSize==None):
        InterpGrid = False
    else:
        InterpGrid = True
        size = np.array(InterpGridSize, ndmin=1, dtype=int)
        if len(size) == 1:
            size = size.repeat(ndim)
        if len(size) != ndim:
            raise ValueError("Specified dimension of interpolation = {:d}, not commensurate with grid dimensionality.".format(len(size)))
        if np.all(size<1):
            raise ValueError("'InterpGridSize' must have values larger than or equal to 1.")

    # ImageData
    isOrthogonal = True
    for i in range(ndim):
        if abs(align[i,i]-bvnorm[i]) >= 1e-4: isOrthogonal = False; break
    if isOrthogonal == True:
        if CenterOrigin == False:
            origin = (base[0,0], base[0,1], base[0,2])
        else:
            origin = base[0] - np.sum(bvec, axis=0)*0.5
            origin = (origin[0], origin[0], origin[0])

        if ndim == 3:
            spacing = (bvnorm[0]/data.shape[0],
                       bvnorm[1]/data.shape[1],
                       bvnorm[2]/data.shape[2])
        else:
            spacing = (bvnorm[0]/data.shape[0],
                       bvnorm[1]/data.shape[1], 1.)
        grid = tvtk.ImageData(spacing=spacing, origin=origin)
        data = np.transpose(data) # To z, y, x as required by vtk
        grid.point_data.scalars = data.flatten()
        grid.point_data.scalars.name = 'scalars'
        grid.dimensions = data.shape

    # StructuredGrid
    elif InterpGrid == False:
        if CenterOrigin == False:
            fbg = 0.; fed = 1.
        else:
            fbg = -0.5; fed = 0.5
        pts = np.meshgrid(
            np.linspace(fbg, fed, data.shape[0]),
            np.linspace(fbg, fed, data.shape[1]),
            np.linspace(fbg, fed, data.shape[2]),
            indexing='ij')
        pts = np.array([i.flatten() for i in pts]).T @ bvec # (nx,ny,nz) * 3
        # Use x, y, z, consistent with coordinates
        grid = tvtk.StructuredGrid(dimensions=(data.shape[2], data.shape[1], data.shape[0]))
        grid.points = pts + base[0]
        grid.point_data.scalars = data.flatten()
        grid.point_data.scalars.name = 'scalars'

    # Interpolated ImageData grid.
    else:
        # Old grid
        pts = np.meshgrid(
            np.linspace(0, 1, data.shape[0]),
            np.linspace(0, 1, data.shape[1]),
            np.linspace(0, 1, data.shape[2]),
            indexing='ij'
        )
        pts = np.array([i.flatten() for i in pts]).T @ bvec
        # New grid
        x = np.linspace(np.min(pts[:,0]), np.max(pts[:,0]), size[0]*data.shape[0])
        y = np.linspace(np.min(pts[:,1]), np.max(pts[:,1]), size[1]*data.shape[1])
        z = np.linspace(np.min(pts[:,2]), np.max(pts[:,2]), size[2]*data.shape[2])
        if CenterOrigin == False:
            origin = base[0] + np.array([x[0], y[0], z[0]])
        else:
            origin = base[0] + np.array([x[0], y[0], z[0]]) - np.sum(bvec, axis=0)*0.5
        origin = (origin[0], origin[1], origin[2])
        # interp on new grid
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        DATA = griddata(pts, data.flatten(), (X, Y, Z), **kwargs)
        del pts, data, X, Y, Z

        DATA = np.transpose(DATA) # To z, y, x as required by vtk
        grid = tvtk.ImageData(spacing=(x[1]-x[0], y[1]-y[0], z[1]-z[0]),
                              origin=origin)
        grid.point_data.scalars = DATA.flatten()
        grid.point_data.scalars.name = 'scalars'
        grid.dimensions = DATA.shape
    # ------------------------------------------------------------------------#
    # Note: The following block is for vtk UnstructuredGrid class. Hexadegra  #
    # in volume plotting are represented by discrete triangles. A tetrahedron #
    # method is needed to divide a hexadegron into 6 detrahedra, which hugely #
    # increases the cost and memory requirement. The code below bulids        #
    # hexadegronal mesh and rarely tested                                     #
    # ------------------------------------------------------------------------#
    # pts = np.meshgrid(
    #     np.linspace(0, 1, data.shape[0]),
    #     np.linspace(0, 1, data.shape[1]),
    #     np.linspace(0, 1, data.shape[2]),
    #     indexing='ij'
    # )
    # pts = np.array([i.flatten() for i in pts]).T @ bvec
    # # list all the hexahedrons
    # ncell = int((data.shape[0]-1)*(data.shape[1]-1)*(data.shape[2]-1))
    # cells = np.zeros([ncell, 9], dtype=int)
    # ptsnew = np.zeros([8*ncell, 3], dtype=float)

    # nslice = int(data.shape[1] * data.shape[2])
    # ncol = data.shape[2]
    # for i in range(data.shape[0]-1):
    #     countx = int(i * (data.shape[1]-1) * (data.shape[2]-1))
    #     for j in range(data.shape[1]-1):
    #         county = int(j * (data.shape[2]-1))
    #         for k in range(data.shape[2]-1):
    #             countz = countx + county + k
    #             nextz = countz + 1
    #             # vertices
    #             countv = int(8*countz)
    #             cells[countz] = [8, countv, countv+1, countv+2, countv+3,
    #                              countv+4, countv+5, countv+6, countv+7]
    #             # hexahedron vertices: O, A, A+B, B, C, A+C, A+B+C, B+C
    #             ptsnew[countv:countv+8] = [
    #                 pts[countz], pts[countz+nslice], pts[countz+nslice+ncol], pts[countz+ncol],
    #                 pts[nextz],  pts[nextz+nslice],  pts[nextz+nslice+ncol],  pts[nextz+ncol]
    #             ]
    # del pts
    # cells = cells.flatten()
    # offset = np.array([i*9 for i in range(ncell)], dtype=int)
    # cell_types = np.array([tvtk.Hexahedron().cell_type for i in range(ncell)])
    # cell_array = tvtk.CellArray()
    # cell_array.set_cells(ncell, cells)
    # grid = tvtk.UnstructuredGrid(points=ptsnew)
    # grid.set_cells(cell_types, offset, cell_array)
    # grid.point_data.scalars = data.flatten()
    # grid.point_data.scalars.name = 'scalars'
    # ------------------------------------------------------------------------#
    return grid


#--------------------------- 2D fields based on Matplotlib -------------------#


def plot_2Dscalar(fig, ax, data, base, levels, contourline, isovalue, colormap, cbar_label,
                  a_range, b_range, rectangle, edgeplot, xticks, yticks, **kwargs):
    """
    Plot 2D scalar field map.

    Args:
        fig (Figure): Matplotlib Figure object
        ax (Axes): Matplotlib Axes object
        data (array): 2D map data, in nY\*nX,
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

    # General grid manipulations are written in nX*nY and O, A, B
    basenew = np.vstack([base[1], base[2], base[0]])
    data = data.T

    basenew, data = GridExpand(basenew, data, [a_range, b_range])
    if rectangle == True:
        basenew, data, X, Y = GridRectangle2D(basenew, data)
        rot,_ = GridRotation2D(basenew) # used later
        X = X.T; Y = Y.T; data = data.T
    else:
        rot,_ = GridRotation2D(basenew)
        basenew = rot.apply(basenew)
        CRDS = GridCoordinates(basenew, data.shape, meshgrid=True)
        X = CRDS[0].T; Y = CRDS[1].T; data = data.T

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

    # plot plane edges at the origin
    if edgeplot == True:
        va = rot.apply((base[2]-base[1]))
        vb = rot.apply((base[0]-base[1]))
        path = np.vstack([[0.,0.,0.], va, va+vb, vb, [0.,0.,0.]])
        ax.plot(path[:,0], path[:,1],'k-', linewidth=1.0)

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

    # General grid manipulations are written in nX*nY and O, A, B
    basenew = np.vstack([base[1], base[2], base[0]])
    data = data.transpose([1, 0, 2])

    baseold = copy.deepcopy(basenew)
    basenew, datax = GridExpand(baseold, data[:, :, 0], [a_range, b_range])
    _, datay = GridExpand(baseold, data[:, :, 1], [a_range, b_range])
    _, dataz = GridExpand(baseold, data[:, :, 2], [a_range, b_range])
    del data
    if rectangle == True:
        baseold = copy.deepcopy(basenew)
        data = np.dstack([datax, datay, dataz])
        basenew, datax, X, Y = GridRectangle2D(baseold, data[:, :, 0])
        _, datay, _, _ = GridRectangle2D(baseold, data[:, :, 1])
        _, dataz, _, _ = GridRectangle2D(baseold, data[:, :, 2])
        rot, _ = GridRotation2D(baseold)
        X = X.T; Y = Y.T;
        data = np.dstack([datax, datay, dataz]).transpose([1, 0, 2])
        del datax, datay, dataz
    else:
        rot, _ = GridRotation2D(basenew)
        basenew = rot.apply(basenew)
        CRDS = GridCoordinates(basenew, datax.shape, meshgrid=True)
        X = CRDS[0].T; Y = CRDS[1].T;
        data = np.dstack([datax, datay, dataz]).transpose([1, 0, 2])
        del datax, datay, dataz

    # get projection and norm of the arrows
    xyz3d = rot.apply(np.eye(3)) # xyz axes in 3D reference framework to plot coord
    dataprj = data @ xyz3d
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
        va = rot.apply((base[2]-base[1]))
        vb = rot.apply((base[0]-base[1]))
        path = np.vstack([[0.,0.,0.], va, va+vb, vb, [0.,0.,0.]])
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


#--------------------------- 2D/3D fields based on MayaVi---------------------#


def plot_3Dscalar(fig, base, data, isovalue, volume_3d, interp, interp_size,
                  display_range, **kwargs):
    """
    Plot 3D scalar field.

    Args:
        fig: MayaVi scence object
        base (array): 4\*3 array of base vectors defining O, A, B, C.
        data (array): nZ\*nY\*nX array of plot data.
        isovalue (float|array): Isovalues of 3D/2D contour plots. A number
                or an array for user-defined values of isosurfaces, **must be
                consistent with ``unit``**.
        volume_3d (bool): Display 3D volumetric data instead of isosurfaces.
        interp (str): Interpolation method. 'no interp', 'linear', 'nearest',
            'slinear', 'cubic'.
        interp_size (list[int]|int): The new size of interpolated data (list)
            or a scaling factor.
        display_range (array): 3\*2 array defining the displayed region.
            Fractional coordinates a, b, c are used.
        \*\*kwargs: Optional keywords passed to MayaVi, listed below.
        colormap (turple|str): Colormap of isosurfaces/heatmaps. Or a 1\*3
            RGB turple from 0 to 1 to define colors. *Not for volume_3d=True*.
        opacity (float): Opacity from 0 to 1. For ``volume_3d=True``, that
            defines the opacity of the maximum value. The opacity of the
            minimum is half of it.
        transparent (bool): Scalar-dependent opacity. *Not for volume_3d=True*.
        vmax (float): Maximum value of colormap.
        vmin (float): Minimum value of colormap.
        title (str): Colorbar title.
        orientation (str): Orientation of colorbar, 'horizontal' or 'vertical'.
        nb_labels (int): The number of labels to display on the colorbar.
        label_fmt (str): The string formater for the labels, e.g., '%.1f'.
    Returns:
        fig: MayaVi scence object
    """
    import copy, warnings, re
    import numpy as np
    try:
        from mayavi import mlab
        from tvtk.util.ctf import PiecewiseFunction
    except ModuleNotFoundError:
        raise ModuleNotFoundError('MayaVi is required for this functionality, which is not in the default dependency list of CRYSTALpytools.')

    #---------------------------------------------------------------------#
    #                                NOTE                                 #
    #---------------------------------------------------------------------#
    # For visualization, data has the dimension of nX*nY*nZ. A transpose  #
    # is needed!                                                          #
    #---------------------------------------------------------------------#

    # Input processing and sanity check
    if data.ndim != 3 or base.shape[0] != 4:
        raise Exception("A nZ*nY*nX array is needed for input data and a 4*3 array is needed for base vectors.")
    data = data.T

    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = np.min(data)
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = np.max(data)

    if volume_3d == False:
        isovalue = np.round(np.array(isovalue, ndmin=1), 12)
        if np.any(isovalue<vmin) or np.any(isovalue>vmax):
            warings.warn("Some of the isovalues are not within the visualized range, vmin = {:.4f}, vmax = {:.4f}.".format(vmin, vmax),
                         stacklevel=2)
        isovalue = isovalue[np.where((isovalue>=vmin)&(isovalue<=vmax))]
        if len(isovalue) < 1:
            raise Exception("No isovalue exists in the visulized range.")
        elif len(isovalue) > 1:
            if 'vmin' not in kwargs.keys():
                vmin = np.min(isovalue)
            if 'vmax' not in kwargs.keys():
                vmax = np.max(isovalue)

    # Interpolation
    if interp != 'no interp':
        data, _ = GridInterpolate(base, data, interp, interp_size)

    # Expansion
    base, data = GridExpand(base, data, display_range)

    # Visualization
    if volume_3d == False:
        grid = tvtkGrid(base, data, CenterOrigin=False, InterpGridSize=None)
        keys = ['colormap', 'opacity', 'transparent', 'vmax', 'vmin']
        keywords = dict(figure=fig,
                        contours=isovalue.tolist(),
                        colormap='jet',
                        vmax=vmax,
                        vmin=vmin)
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v

        plot = mlab.pipeline.iso_surface(grid, **keywords)
    else:
        if vmin*vmax >= 0: nullvalue = vmin - (vmax-vmin)
        else: nullvalue = 0.

        grid = tvtkGrid(base, data, CenterOrigin=False,
                        InterpGridSize=interp_size, fill_value=nullvalue)
        keys = ['vmax', 'vmin']
        keywords = dict(figure=fig, vmax=vmax, vmin=vmin)
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v

        plot = mlab.pipeline.volume(grid, **keywords)
        # transparency
        if 'opacity' in kwargs.keys():
            opacity = kwargs['opacity']
        else:
            opacity = 1.0

        if vmin*vmax >= 0: # vmin has the lowest opacity
            vals = np.linspace(vmin, vmax, 10)
            opci = np.linspace(opacity*0.5, opacity, 10)
        else: # 0 has the lowest opacity
            vals = np.hstack([np.linspace(vmin, 0, 5, endpoint=False),
                              np.linspace(vmax, 0, 5, endpoint=False)[::-1]])
            opci = np.hstack([np.linspace(opacity, opacity*0.5, 5, endpoint=False),
                              np.linspace(opacity*0.5, opacity, 5, endpoint=False)])

        otf = PiecewiseFunction()
        otf.add_point(nullvalue, 0.)
        for i, o in zip(vals, opci):
            otf.add_point(i, o)
        plot._otf = otf
        plot._volume_property.set_scalar_opacity(otf)

    keys = ['title', 'orientation', 'nb_labels', 'label_fmt']
    keywords = dict(object=plot)
    for k, v in zip(kwargs.keys(), kwargs.values()):
        if k in keys: keywords[k] = v
    mlab.scalarbar(**keywords)

    return fig


def plot_3Dplane(fig, base, data, levels, contour_2d, interp, interp_size,
                 display_range, **kwargs):

    """
    Plot oriented 2D scalar fields in 3D structure.

    Args:
        fig: MayaVi scence object
        base (array): 3\*3 array of base vectors defining A, B, C.
        data (array): nY\*nX array of plot data.
        levels (float|array): Number of Isovalues of 2D contour plots, equally
            spaced between ``vmin`` and ``vmax``, or an array for user-define
            values of contour lines, **'contour_2d=True' only**.
        contour_2d (bool): Display black contour lines over the colored surface.
        interp (str): Interpolation method. 'no interp', 'linear', 'nearest',
            'slinear', 'cubic'.
        interp_size (list[int]|int): The new size of interpolated data (list)
            or a scaling factor.
        display_range (array): 2\*2 array defining the displayed region.
            Fractional coordinates a, b are used.
        \*\*kwargs: Optional keywords passed to MayaVi, listed below.
        colormap (turple|str): Colormap of heatmaps. Or a 1\*3 RGB turple from
            0 to 1 to define colors.
        color (turple): Color of contour lines. *'contour_2d=True' only*.
        line_width (float): Line width of contour plots. *'contour_2d=True' only*.
        opacity (float): Opacity from 0 to 1. For ``volume_3d=True``, that
            defines the opacity of the maximum value. The opacity of the
            minimum is half of it.
        transparent (bool): Scalar-dependent opacity. *Not for volume_3d=True*.
        vmax (float): Maximum value of colormap.
        vmin (float): Minimum value of colormap.
        title (str): Colorbar title.
        orientation (str): Orientation of colorbar, 'horizontal' or 'vertical'.
        nb_labels (int): The number of labels to display on the colorbar.
        label_fmt (str): The string formater for the labels, e.g., '%.1f'.
    Returns:
        fig: MayaVi scence object
    """
    import copy, warnings, re
    import numpy as np
    try:
        from mayavi import mlab
    except ModuleNotFoundError:
        raise ModuleNotFoundError('MayaVi is required for this functionality, which is not in the default dependency list of CRYSTALpytools.')

    #---------------------------------------------------------------------#
    #                                NOTE                                 #
    #---------------------------------------------------------------------#
    # For visualization, data has the dimension of nX*nY. A transpose is  #
    # needed!                                                             #
    #---------------------------------------------------------------------#

    # Input processing and sanity check
    if data.ndim != 2 or base.shape[0] != 3:
        raise Exception("A nY*nX array is needed for input data and a 3*3 array is needed for base vectors.")
    base = np.vstack([base[1], base[2], base[0]])
    data = data.T

    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = np.min(data)
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = np.max(data)

    if contour_2d == True:
        levels = np.round(np.array(levels, ndmin=1), 12)
        if levels.shape[0] >= 1:
            if np.any(levels<vmin) or np.any(levels>vmax):
                warnings.warn("Some of the contours are not within the visualized range, vmin = {:.4f}, vmax = {:.4f}.".format(vmin, vmax),
                              stacklevel=2)
            levels = levels[np.where((levels>=vmin)&(levels<=vmax))]

        if levels.shape[0] < 1:
            raise Exception("No contour exists in the visulized range. vmin = {:.4f}, vmax = {:.4f}.".format(vmin, vmax))

    # Interpolation
    if interp != 'no interp':
        data, _ = GridInterpolate(base, data, interp, interp_size)

    # Expansion
    base, data = GridExpand(base, data, display_range)

    # Visualization
    grid = tvtkGrid(base, data, CenterOrigin=False, InterpGridSize=None)
    keys = ['colormap', 'opacity', 'transparent', 'vmax', 'vmin']
    keywords = dict(figure=fig,
                    colormap='jet',
                    vmax=vmax,
                    vmin=vmin)

    for k, v in zip(kwargs.keys(), kwargs.values()):
        if k in keys: keywords[k] = v

    surf = mlab.pipeline.surface(grid, **keywords)

    if contour_2d == True:
        keys = ['color', 'line_width', 'opacity', 'transparent', 'vmax', 'vmin']
        keywords = dict(figure=fig,
                        color=(0,0,0),
                        contours=levels.tolist(),
                        vmax=vmax,
                        vmin=vmin)
        for k, v in zip(kwargs.keys(), kwargs.values()):
            if k in keys: keywords[k] = v

        lines = mlab.pipeline.contour_surface(grid, **keywords)

    keys = ['title', 'orientation', 'nb_labels', 'label_fmt']
    keywords = dict(object=surf)
    for k, v in zip(kwargs.keys(), kwargs.values()):
        if k in keys: keywords[k] = v
    mlab.scalarbar(**keywords)

    return fig


#-----------------------------------------------------------------------------#
# Note: The following function passed tests of Mayavi examples, but not the   #
# more complicated practical examples. Also, the volume rendering is not      #
# available for x3d.                                                          #
#-----------------------------------------------------------------------------#
# def MayaViRenderHTML(x3d, html):
#     """
#     Convert X3D (XML based) scene saved by MayaVi into HTML based X3DOM file
#     for embeddings and savings.

#     Args:
#         x3d (str): X3D filename.
#         html (str): HTML filename.
#     Returns:
#         None
#     """
#     import re

#     file = open(x3d, 'r')
#     data = file.readlines()
#     file.close()

#     # Find the scene and match the keywords
#     scenebg = 0; sceneed = 0
#     keysave = ''
#     for i in range(len(data)):
#         if re.match(r'^\s*<Scene>\s*$', data[i], re.IGNORECASE):
#             scenebg = i
#         elif re.match(r'^\s*<\/Scene>\s*$', data[i], re.IGNORECASE):
#             sceneed = i+1
#         elif re.match(r'^\s*<[A-Z,a-z]+.*\/>\s*$', data[i], re.IGNORECASE):
#             indent = ''
#             for j in range(len(data[i])):
#                 if data[i][j] != ' ': break
#                 indent += data[i][j]
#             line = data[i].strip()
#             keyword = line.split()[0][1:]
#             data[i] = indent + line[:-2] + '></{}>\n'.format(keyword)
#         elif re.match(r'^\s*<[A-Z,a-z]+((?!>).)*$', data[i], re.IGNORECASE): # Not closed bracket
#             keysave = data[i].strip().split()[0][1:]
#         elif re.match(r'^((?!<).)*\/>\s*$', data[i], re.IGNORECASE):
#             indent = ''
#             for j in range(len(data[i])):
#                 if data[i][j] != ' ': break
#                 indent += data[i][j]
#             line = data[i].strip()
#             data[i] = indent + line[:-2] + '></{}>\n'.format(keysave)

#     # write into HTML
#     htmldata = """\
# <html>
#   <head>
#     <title>
#       CRYSTALpytools 3D visualization
#     </title>
#     <link rel="stylesheet" type="text/css" href="https://www.x3dom.org/release/x3dom.css">
#     </link>
#     <script type="text/javascript" src="https://www.x3dom.org/release/x3dom.js">
#     </script>
#   </head>
# <body>
# <X3D>
# <!--Inserting Generated X3D Scene-->
# {}
# <!--End of Inserted Scene-->
# </X3D>
# </body>
# </html>
# """.format(''.join([i for i in data[scenebg:sceneed]]))
#     file = open(html, 'w')
#     file.write("%s" % htmldata)
#     file.close()
#     return

