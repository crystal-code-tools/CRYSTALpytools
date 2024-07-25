#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module for electron transport property analysis.
"""
class Tensor():
    """
    The class for electron transport properties described by tensors, including
    'KAPPA', 'SIGMA', 'SIGMAS' and 'SEEBECK'. Units must be consistent with the
    output files.

    Args:
        temperature (array): Temperature in K
        potential (array): Chemical potential in eV.
        carrier_density (array): nT\*nPot\*nSpin array of carrier density in
            cm:math:`^{-3}`.
        tensor (array): nT\*nPot\*nDimen\*nSpin array of flattened tensor
            elements. nDimen = 6 for 3D systems, 3 for 2D and 1 for 1D.
        struc (CStructure): Extended Pymatgen Structure object.
        type (str): 'KAPPA', 'SIGMA', 'SIGMAS' or 'SEEBECK'.
        unit (str): Same as output file.

    Returns:
        self (Tensor): Attributes: 'T', 'mu', 'carrier', 'data', 'type',
            'struc', 'unit' and 'spin'
    """
    def __init__(self, temperature, potential, carrier_density, tensor, struc,
                 type, unit):
        import numpy as np

        self.T = np.array(temperature, dtype=float)
        self.mu = np.array(potential, dtype=float)
        self.carrier = np.array(carrier_density, dtype=float)
        self.data = np.array(tensor, dtype=float)
        self.type = type
        self.struc = struc
        self.unit = unit
        if len(self.T) != self.carrier.shape[0] or len(self.T) != self.data.shape[0]:
            raise ValueError("Inconsistent dimensionalities between temperature and input data.")
        if len(self.mu) != self.carrier.shape[1] or len(self.mu) != self.data.shape[1]:
            raise ValueError("Inconsistent dimensionalities between chemical potential and input data.")
        self.spin = self.data.shape[-1]

    @classmethod
    def from_file(cls, *file, output=None, method='normal'):
        """
        Read electron transport properties by the BOLTZTRA keyword, including
        'KAPPA', 'SIGMA', 'SIGMAS', and 'SEEBECK'. Though currently the
        geometry information is not required, it is saved if the standard
        output file is given.

        Also allow for generating the following properites not printed by
        CRYSTAL:

        * Power Factor: :math:`S^{2}\\sigma`. Any 2 from 'SEEBECK', 'SIGMA' and
            'SIGMAS' are required.  
        * ZT: :math:`\\frac{S^{r}\\sigma T}{\\kappa}`. 'KAPPA' and any 2
            from 'SEEBECK', 'SIGMA' and 'SIGMAS' are required.

        .. note::

            Make sure all the entries are from the same calculation. The code
            only checks the dimensionalities of tensors.

        Args:
            file (str): 'DAT' files by CRYSTAL BOLTZTRA keyword. Extendable
                only when ``method='power factor'`` or ``method='zt'``.
            output (str): Properties output file.
            method (str): 'normal' to read 1 file, or 'power factor' or 'zt'.
                For meanings, see above.
        Returns:
            cls (Tensor)
        """
        from CRYSTALpytools.crystal_io import Properties_output
        import copy
        import numpy as np

        if method.lower() == 'normal':
            obj = Properties_output(output).read_transport(file[0])
        elif method.lower() == 'power factor':
            if len(file) != 2:
                raise Exception("For power factor you need 2 files.")

            obj1 = Properties_output(output).read_transport(file[0])
            obj2 = Properties_output(output).read_transport(file[1])
            obj = cls.get_power_factor(obj1, obj2)
        elif method.lower() == 'zt':
            if len(file) != 3:
                raise Exception("To get power factor you need 3 files.")

            obj1 = Properties_output(output).read_transport(file[0])
            obj2 = Properties_output(output).read_transport(file[1])
            obj3 = Properties_output(output).read_transport(file[2])
            obj = cls.get_zt(obj1, obj2, obj3)
        else:
            raise ValueError("Unknown method: '{}'.".format(method))

        return obj

    @classmethod
    def get_power_factor(cls, obj1, obj2):
        """
        Get thermoelectric power factor object from any 2 objects of 'SEEBECK',
        'SIGMA' or 'SIGMAS' types.

        .. math::

            \\text{Power Factor} = S^{2}\\sigma

        .. note::

            Make sure all the entries are from the same calculation. The code
            only checks the dimensionalities of tensors. The geometry follows
            the first entry.

        Args:
            obj1 (Tensor): Tensor object 1 of 'SEEBECK', 'SIGMA' or 'SIGMAS' types.
            obj2 (Tensor): Tensor object 2 of 'SEEBECK', 'SIGMA' or 'SIGMAS' types.
        Returns:
            cls (Tensor): 'POWERFACTOR' type of object, in 'W/m/K^2'.
        """
        if obj1.shape != obj2.shape:
            raise Exception("Inconsistent shapes for input objects.")

        if obj1.type == 'SEEBECK' and obj2.type == 'SIGMA':
            tensnew = obj1.data**2 * obj2.data
        elif obj1.type == 'SIGMA' and obj2.type == 'SEEBECK':
            tensnew = obj2.data**2 * obj1.data
        elif obj1.type == 'SEEBECK' and obj2.type == 'SIGMAS':
            tensnew = obj1.data**2 * obj2.data/obj1.data
        elif obj1.type == 'SIGMAS' and obj2.type == 'SEEBECK':
            tensnew = obj2.data**2 * obj1.data/obj2.data
        elif obj1.type == 'SIGMA' and obj2.type == 'SIGMAS':
            tensnew = (obj2.data/obj1.data)**2 * obj1.data
        elif obj1.type == 'SIGMAS' and obj2.type == 'SIGMA':
            tensnew = (obj1.data/obj2.data)**2 * obj2.data
        else:
            TypeError("For power factor you need 2 from 'SEEBECK', 'SIGMA' or 'SIGMAS'.")

        return cls(obj1.T, obj1.mu, obj1.carrier, tensnew, obj1.struc, 'POWERFACTOR', 'W/m/K^2')

    @classmethod
    def get_zt(cls, obj1, obj2, obj3):
        """
        Get thermoelectric dimensionless figure of merit (ZT) object from a
        'KAPPA' object and any 2 from 'SEEBECK', 'SIGMA' and 'SIGMAS'.

        .. math::

            ZT = \\frac{S^{r}\\sigma T}{\\kappa}

        .. note::

            Make sure all the entries are from the same calculation. The code
            only checks the dimensionalities of tensors. The geometry follows
            the 'KAPPA' entry.

        Args:
            obj1 (Tensor): Tensor object 1 of 'KAPPA', 'SEEBECK', 'SIGMA' or 'SIGMAS' types.
            obj2 (Tensor): Tensor object 2 of 'KAPPA', 'SEEBECK', 'SIGMA' or 'SIGMAS' types.
            obj3 (Tensor): Tensor object 3 of 'KAPPA', 'SEEBECK', 'SIGMA' or 'SIGMAS' types.
        Returns:
            cls (Tensor): 'ZT' type of object, in 'dimensionless'.
        """
        import numpy as np

        objs = [obj1, obj2, obj3]
        types = np.array([obj1.type, obj2.type, obj3.type])
        if (obj1.shape != obj2.shape) or (obj1.shape != obj3.shape) \
        or (obj2.shape != obj3.shape):
            raise Exception("Inconsistent shapes for input objects.")

        if len(np.where(types=='KAPPA')[0]) != 1:
            raise TypeError("For ZT, you need 1 and only one 'KAPPA' file.")
        else:
            objk = objs[np.where(types=='KAPPA')[0][0]]
            obja = objs[np.where(types!='KAPPA')[0][0]]
            objb = objs[np.where(types!='KAPPA')[0][1]]

        # S^2 \sigma / \kappa
        obj = cls.get_power_factor(obja, objb)
        tensnew = copy.deepcopy(obj.data); del obj, obja, objb
        tensnew = tensnew / objk.data
        # S^2 \sigma T / \kappa
        for irow in range(len(tensnew)):
            tensnew[i] = tensnew[i] * objk.T[i]
        return cls(objk.T, objk.mu, objk.carrier, tensnew, objk.struc, 'ZT', 'dimensionless')

    def plot(self, x_axis='potential', x_range=[], direction='xx', spin='sum',
             plot_series=[], plot_label=None, plot_color=None, plot_linestyle=None,
             plot_linewidth=None, zero_color='tab:gray', zero_linestyle='-',
             zero_linewidth=1., layout=None, add_title=True, figsize=[6.4, 4.8],
             legend='upper left', sharex=True, sharey=True, fontsize=14,
             fig=None, **kwargs):
        """
        Plot tensor-like transport properties in multiple ways:

        1. X_axis: Chemical potential :math:`\\mu`; Plot series: Temperature :math:`T`.  
        2. X_axis: Carrier density :math:`\\rho(\\mu; T)`; Plot series: Temperature :math:`T`.  
        3. X_axis: Temperature :math:`T`; Plot series: Chemical potential :math:`\\mu`.

        Args:
            x_axis (str): X axis options, 'potential', 'carrier' or
                'temperature'. For meanings, see above.
            x_range (list): Display range of x axis. Y axis is self-adaptive.
            direction (str|list): Depending on the dimensionality of the system,
                including 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy' and
                'zz'. A list of options will generate nDirect\*1 subplots. The
                direction of each subplot is annotated on the upper left corner.
            spin (str): Spin-polarized systems only. Electron spins to plot.
                'up', 'down' or 'sum'.
            plot_series (list|array|float): **Values** of the plot series.
                Should be commensurate with the choice of ``x_axis``.
            plot_label (list|str|None): None for default values: ``plot_series``
                and its unit. If str, prefixing that entry in front of the
                default value. If list, it should be 1\*nPlotSeries for every
                plot line.
            plot_color (list|str|None): Similar to ``electronics.ElectronDOS.plot()``.
                If str, use the same color for all the plot lines. If
                1\*nPlotSeries, use the same color for every plot line. If
                2\*nPlotSeries, use different colors for p-type and n-type
                carrier properties.
            plot_linestyle (list|str|None): Similar to ``electronics.ElectronDOS.plot()``.
                See explanations of ``plot_color``.
            plot_linewidth (list|float|None): Similar to ``electronics.ElectronDOS.plot()``.
                See explanations of ``plot_color``.
            zero_color (str): Color of the 0 value line.
            zero_linestyle (str): Linestyle of the 0 value line.
            zero_linewidth (float): Linewidth of the 0 value line.
            layout (list[int]): Layout of subplots. By default it is nDirection\*1
                gird of figures.
            add_title (bool): Whether to add the plotted property as title.
            figsize (list): Figure size.
            legend (str|None): Location of legend. None for not adding legend.
            sharex (bool): Share x axis for multiple subplots.
            sharey (bool): Share y axis for multiple subplots.
            fontsize (float|int): Font size of the title, subplot capations
                (direction), x and y axis capations.
            fig (Figure): *Developer Only*
            \*\*kwargs: Other arguments passed to the matplotlib ``Axes.plot()``
                method. Applied to all the plots.
        Returns:
            fig (Figure): Matplotlib Figure object.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from CRYSTALpytools.base.plotbase import _plot_label_preprocess

        # x axis
        x_axis = x_axis.lower()
        if x_axis != 'potential' and x_axis != 'carrier' and x_axis != 'temperature':
            raise ValueError("Unknown x-axis style: '{}'.".format(x_axis))
        if x_range == []:
            if x_axis == 'potential':
                x_range = [np.min(self.mu), np.max(self.mu)]
            elif x_axis == 'carrier':
                x_range = [np.min(np.abs(self.carrier)), np.max(np.abs(self.carrier))]
            else:
                x_range = [np.min(self.T), np.max(self.T)]
        x_range = [np.min(x_range), np.max(x_range)]

        # projections in the same plot
        if plot_series == []:
            if x_axis == 'potential' or x_axis == 'carrier':
                prj = range(len(self.T))
            else:
                prj = range(len(self.mu))
        else:
            plot_series = np.array(plot_series, ndmin=1)
            prj = []
            if x_axis == 'potential' or x_axis == 'carrier':
                for i in plot_series:
                    lst = np.where(self.T==i)[0]
                    if len(lst) == 0:
                        raise ValueError("Input plot_series {:.2f} does not exist in object's temperature series.".format(i))
                    prj.append(lst[0])
            else:
                for i in plot_series:
                    lst = np.where(self.mu==i)[0]
                    if len(lst) == 0:
                        raise ValueError("Input plot_series {:.2f} does not exist in object's potential series.".format(i))
                    prj.append(lst[0])
        nprj = len(prj)

        # directions (number of subplots)
        if self.data.shape[2] == 9:
            indices = {'xx' : 0, 'xy' : 1, 'xz' : 2, 'yx' : 3, 'yy' : 4, 'yz' : 5,
                       'zx' : 6, 'zy' : 7, 'zz' : 8}
        else:
            indices = {'xx' : 0, 'xy' : 1, 'yx' : 2, 'yy' : 3}

        direction = np.array(direction, ndmin=1)
        dir = [indices[i] for i in direction]
        nplt = len(dir)

        # plot setups
        if x_axis == 'potential' or x_axis == 'carrier':
            dflabel = ['{:>5.0f} K'.format(self.T[i]) for i in prj]
        else:
            dflabel = ['{:>6.2f} eV'.format(self.mu[i]) for i in prj]

        if np.all(plot_label==None):
            plot_label = dflabel
        elif isinstance(plot_label, str):
            plot_label = ['{} {}'.format(plot_label, i) for i in deflabel]
        elif isinstance(plot_label, list) or isinstance(plot_label, array):
            nlabel = len(plot_label)
            plot_label = [plot_label[i%nlabel] for i in range(nprj)]
        else:
            raise TypeError('Unknown type of plot label.')
        ## get a pseudo band input
        bands = np.zeros([nprj, 1, 1, 1], dtype=float)
        commands = _plot_label_preprocess(bands, plot_label, plot_color,
                                          plot_linestyle, plot_linewidth)

        # plotting
        if np.all(fig==None):
            added_plot = False
            if np.all(layout==None):
                fig, ax = plt.subplots(len(direction), 1, sharex=sharex,
                                       sharey=sharey, figsize=figsize,
                                       layout='constrained')
            else:
                fig, ax = plt.subplots(layout[0], layout[1], sharex=sharex,
                                       sharey=sharey, figsize=figsize,
                                       layout='constrained')
        else:
            added_plot = True

        y_range = []
        for idir, ax in enumerate(fig.axes):
            ax.hlines(0, x_range[0], x_range[1], colors=zero_color,
                      linestyle=zero_linestyle, linewidth=zero_linewidth)
            ## spins
            if x_axis == 'potential':
                if self.spin == 1 or spin.lower() == 'up':
                    y = self.data[prj, :, dir[idir], 0]
                    carrier = self.carrier[prj, :, 0]
                elif self.spin == 2 and spin.lower() == 'sum':
                    y = np.sum(self.data[prj, :, dir[idir], :], axis=2)
                    carrier = np.sum(self.carrier[prj, :, :], axis=2)
                else:
                    y = self.data[prj, :, dir[idir], 1]
                    carrier = self.carrier[prj, :, 1]
                ## limit plot range
                idx_x = np.where((self.mu>=x_range[0])&(self.mu<=x_range[1]))[0]
                x = self.mu[idx_x]
                y = y[:, idx_x]
                carrier = carrier[:, idx_x]
                ## divide the plot by p and n type carriers and plot
                for p in range(nprj):
                    idx_p = np.where(carrier[p]>=0)[0]
                    idx_n = np.where(carrier[p]<0)[0]
                    if len(idx_p) > 0:
                        ax.plot(x[idx_p], y[p, idx_p].flatten(), label=commands[0][p][0],
                                color=commands[1][p][0], linestyle=commands[2][p][0],
                                linewidth=commands[3][p][0], **kwargs)
                    if len(idx_n) > 0:
                        ax.plot(x[idx_n], y[p, idx_n].flatten(), color=commands[1][p][1],
                                linestyle=commands[2][p][1], linewidth=commands[3][p][1], **kwargs)
                    # linearly interpolate the gap, noticed when x_axis = 'potential' only
                    if len(idx_p) > 0 and len(idx_n) > 0:
                        lastppt = [x[idx_p[-1]], y[p, idx_p[-1]]]
                        firstnpt = [x[idx_n[0]], y[p, idx_n[0]]]
                        midpt = [(lastppt[0] + firstnpt[0]) / 2, (lastppt[1] + firstnpt[1]) / 2]
                        ax.plot([lastppt[0], midpt[0]], [lastppt[1], midpt[1]],
                                color=commands[1][p][0], linestyle=commands[2][p][0],
                                linewidth=commands[3][p][0], **kwargs)
                        ax.plot([midpt[0], firstnpt[0]], [midpt[1], firstnpt[1]],
                                color=commands[1][p][1], linestyle=commands[2][p][1],
                                linewidth=commands[3][p][1], **kwargs)
            elif x_axis == 'carrier':
                if self.spin == 1 or spin.lower() == 'up':
                    x = np.abs(self.carrier[prj, :, 0])
                    y = self.data[prj, :, dir[idir], 0]
                    carrier = self.carrier[prj, :, 0]
                elif self.spin == 2 and spin.lower() == 'sum':
                    x = np.abs(np.sum(self.carrier[prj, :, :], axis=2))
                    y = np.sum(self.data[prj, :, dir[idir], :], axis=2)
                    carrier = np.sum(self.carrier[prj, :, :], axis=2)
                else:
                    x = np.abs(self.carrier[prj, :, 1])
                    y = self.data[prj, :, dir[idir], 1]
                    carrier = self.carrier[prj, :, 1]
                newx = []; newy = []
                for p in range(nprj):
                    ## limit plot range: dimens of x, y elements might not be uniform
                    idx_x = np.where((x[p, :]>=x_range[0])&(x[p, :]<=x_range[1]))[0]
                    newx = x[p, idx_x]
                    newy = y[p, idx_x]
                    newcarrier = carrier[p, idx_x]
                    ## divide the plot by p and n type carriers and plot
                    idx_p = np.where(newcarrier>=0)[0]
                    idx_n = np.where(newcarrier<0)[0]
                    if len(idx_p) > 0:
                        ax.plot(newx[idx_p], newy[idx_p], label=commands[0][p][0],
                                color=commands[1][p][0], linestyle=commands[2][p][0],
                                linewidth=commands[3][p][0], **kwargs)
                    if len(idx_n) > 0:
                        ax.plot(newx[idx_n], newy[idx_n], color=commands[1][p][1],
                                linestyle=commands[2][p][1], linewidth=commands[3][p][1],
                                **kwargs)
                    ax.set_xscale('log')
            else:
                if self.spin == 1 or spin.lower() == 'up':
                    y = self.data[:, prj, dir[idir], 0]
                    carrier = self.carrier[:, prj, 0]
                elif self.spin == 2 and spin.lower() == 'sum':
                    y = np.sum(self.data[:, prj, dir[idir], :], axis=2)
                    carrier = np.sum(self.carrier[:, prj, :], axis=2)
                else:
                    y = self.data[:, prj, dir[idir], 1]
                    carrier = self.carrier[:, prj, 1]
                ## limit plot range
                idx_x = np.where((self.T>=x_range[0])&(self.T<=x_range[1]))[0]
                x = self.T[idx_x]
                y = y[idx_x, :]
                carrier = carrier[idx_x, :]
                ## divide the plot by p and n type carriers and plot
                for p in range(nprj):
                    idx_p = np.where(carrier[:, p]>=0)[0]
                    idx_n = np.where(carrier[:, p]<0)[0]
                    if len(idx_p) > 0:
                        ax.plot(x[idx_p], y[idx_p, p].flatten(), label=commands[0][p][0],
                                color=commands[1][p][0], linestyle=commands[2][p][0],
                                linewidth=commands[3][p][0], **kwargs)
                    if len(idx_n) > 0:
                        ax.plot(x[idx_n], y[idx_n, p].flatten(), color=commands[1][p][1],
                                linestyle=commands[2][p][1], linewidth=commands[3][p][1],
                                **kwargs)

            y_range.append([np.min(y), np.max(y)])

        # plot setups
        for idir, ax in enumerate(fig.axes):
            if np.all(legend!=None) and np.all(commands[0]!=None):
                ax.legend(loc=legend)
            ax.set_xlim(x_range)
            if sharey == False:
                ax.text(x_range[1], y_range[idir][1], direction[idir], fontsize=fontsize,
                        horizontalalignment='right', verticalalignment='top')
                ax.set_ylim(y_range[idir])
            else:
                y_range = [np.min(y_range), np.max(y_range)]
                ax.text(x_range[1], y_range[1], direction[idir], fontsize=fontsize,
                        horizontalalignment='right', verticalalignment='top')
                ax.set_ylim(y_range)

        if added_plot == False:
            fig.supylabel('{} ({})'.format(self.type.capitalize(), self.unit),
                          fontsize=fontsize)
            if x_axis == 'potential':
                fig.supxlabel('Chemical Potential (eV)', fontsize=fontsize)
            elif x_axis == 'carrier':
                fig.supxlabel('Carrier Density (cm$^{-3}$)', fontsize=fontsize)
            else:
                fig.supxlabel('Temperature (K)', fontsize=fontsize)
            if add_title == True:
                fig.suptitle(self.type)

        return fig


class Distribution():
    """
    The class for electron transport properties described by distribution,
    currently transport distribution function (TDF) only. Units must be
    consistent with the output files.

    Args:
        energy (array): Energy in eV.
        distr (array): nEnergy\*nDimen\*nSpin Distribution function
        carrier_density (array): nT\*nPot\*nSpin array of carrier density in
            cm:math:`^{-3}`. nDimen = 6 for 3D systems, 3 for 2D and 1 for 1D.
        struc (CStructure): Extended Pymatgen Structure object.
        type (str): 'TDF'.
        unit (str): Same as output file.

    Returns:
        self (Distribution): Attributes: 'energy', 'function', 'type', 'struc',
            'unit' and 'spin'
    """
    def __init__(self, energy, distr, struc, type, unit):
        import numpy as np

        self.energy = np.array(energy, dtype=float)
        self.function = np.array(distr, dtype=float)
        self.type = type
        self.struc = struc
        self.unit = unit
        if len(self.energy) != self.function.shape[0]:
            raise ValueError("Inconsistent dimensionalities between energy and distribution function.")
        self.spin = self.function.shape[-1]

    @classmethod
    def from_file(cls, boltztra_out, output=None):
        """
        Read electron transport distribution functions ('TDF.DAT') by the
        BOLTZTRA keyword. Though currently the geometry information is not
        required, it is saved if the standard output file is given.

        Args:
            boltztra_out (str): 'DAT' files by CRYSTAL BOLTZTRA keyword.

        Returns:
            cls (Distribution)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_transport(file)

