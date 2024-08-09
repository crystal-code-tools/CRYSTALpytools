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