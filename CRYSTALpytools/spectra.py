#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and methods for spectra.
"""
class XRD():
    """
    The class for X-ray diffraction spectra.

    Args:
        theta (array): **2:math:`\\theta`** values in degree.
        spectra (array): Spectra intensities.
    """
    def __init__(self, theta, spectra):
        import numpy as np

        self.theta = np.array(theta, dtype=float)
        self.spectra = np.array(spectra, dtype=float)

    @classmethod
    def from_file(cls, output, option='LP'):
        """
        Read XRD spectra from the standard screen output of properties
        calculation.

        Args:
            output (str): Output filename.
            option (str): 'NC' for no correction (The 'INTENS' col); 'LP' for
                Lorentz and polarization effects ('INTENS-LP') and 'DW' for LP
                with Debye-Waller thermal factors ('INTENS-LP-DW').
        Returns:
            cls (XRD)
        """
        from CRYSTALpytools.crystal_io import Properties_output
        return Properties_output(output).read_XRDspec(option=option)

    def plot(self, theta_range=[], normalize=True, title=None,
             figsize=[6.4, 4.8], fontsize=14, **kwargs):
        """
        Plot XRD spectra.

        Args:
            theta_range (list): 1\*2 list of theta range in degree.
            normalize (bool): Normalize the maximum intensity to 100.
            title (str|None): The title of the plot. 'None' for no title.
            figsize (list): Matplotlib figure size.
            fontsize (int): Fontsize of the axis label and title.
            \*\*kwargs: Other parameters passed to matplotlib ``Axes.plot()``
                method.

        Returns:
            fig (Figure): Matplotlib figure.
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import copy

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        spectra = copy.deepcopy(self.spectra)
        if normalize == True:
            spectra = spectra / np.max(spectra) * 100; maxv = 100
        else:
            maxv = np.max(spectra)

        ax.plot(self.theta, self.spectra, **kwargs)

        if len(theta_range) == 0:
            theta_range = [np.min(self.theta), np.max(self.theta)]

        ax.set_xlim([np.min(theta_range), np.max(theta_range)])
        ax.set_xlabel(r'2$\theta^{\circ}$', fontsize=fontsize)
        ax.set_ylabel(r'Intensity (arb. u.)', fontsize=fontsize)
        ax.set_ylim([0, maxv])
        _ = ax.get_yaxis().set_ticks([])
        if np.all(title!=None): ax.set_title(title, fontsize=fontsize)
        return fig


class IR():
    """
    The class for infrared spectrum. Unit: in principle, should always in
    'cm :math:`^{-1}`'.

    Args:
        freq (array): Frequencies. Must be commensurate with unit.
        absorb (array): Absorbance spectrum, nType\*nFreq array.
        reflec (array): Reflectance spectrum along inequivalent polarization
            directions (periodic systems only), nDir\*nFreq array.
        type (str): 'molecule' or 'crystal'
        unit (str): 'cm-1' or 'THz'.
    """
    def __init__(self, freq, absorb, reflec, type, unit='cm-1'):
        import numpy as np

        self.frequency = np.array(freq, dtype=float)
        self.absorbance = np.array(absorb, dtype=float)
        self.reflectance = np.array(reflec, dtype=float)
        self.type = type.lower()
        self.unit = unit
        if self.type == 'molecule':
            self.reflectance = np.array([])

    @classmethod
    def from_file(self, specfile, output=None):
        """
        Generate an ``IR`` object from output and 'IRSPEC.DAT' files of CRYSTAL.

        Args:
            specfile (str): The 'IRSPEC.DAT' file.

        Returns:
            cls (IR)
        """
        from CRYSTALpytools.crystal_io import Crystal_output
        return Crystal_output(output).get_spectra(specfile, 'IRSPEC')

    def plot(self, unit='cm-1', option='LG', normalize=True, REFL_overlap=True,
             shift=0, label=None, color=None, linestyle=None, linewidth=None,
             x_range=[], title=None, figsize=[6.4, 4.8], legend='upper left',
             sharey=True, fontsize=14, fig=None, **kwargs):
        """
        Plot IR spectra into the same axis.

        Args:
            unit (str): X axis unit. 'cm :math:`^{-1}`' or 'nm'.
            option (str): Broadening method. 'LG' for Lorentzian-Gaussian, 'V'
                for Voigt, 'RS' for Rayleigh spherical particles, 'RE' for
                Rayleigh with elipsoid particles, 'REFL' for reflectance spectra
                with 'LG'. *Periodic systems only*.
            normalize (bool): Normalize the maximum intensity to 100.
            REFL_overlap (bool): *For ``option='REFL'`` only* If more than 1
                inequivalent directions exists, whether to plot REFL data into
                the same plot or into subplots.
            shift (float): *For ``option='REFL'`` and ``REFL_overlap=False``
                only* If multiple spectra are plotted, shift them up by the
                given value. Shift length is the value after normalization if
                ``normalize=True``.
            label (list|str|None): *For ``option='REFL'`` only* List of
                plot labels. 'None' for the default values ('# \<number\>') and
                string for prefix the string before number. Otherwise should be
                a 1\*nPlot list.
            color (list|str|None): *For ``option='REFL'`` only*  If str,
                use the same color for all the plot lines. If 1\*nPlot, use the
                color defined for every plot line. 'None' for default values
                (matplotlib Tableau Palette).
            linestyle (list|str|None): *For ``option='REFL'`` only* See
                explanations of ``color``.
            linewidth (list|float|None): *For ``option='REFL'`` only* See
                explanations of ``color``.
            x_range (list): 1\*2 list of x axis range.
            title (str|None): The title of the plot. 'None' for no title.
            figsize (list): Matplotlib figure size.
            legend (str|None): Location of legend. None for not adding legend.
            sharey (bool): Whether to share the y-axis among subplots. Share x
                is enforced.
            fontsize (int): Fontsize of the axis label and title.
            fig (Figure): *Developer Only* Matplotlib figure object.
            \*\*kwargs: Other parameters passed to matplotlib ``Axes.plot()``
                method.

        Returns:
            fig (Figure): Matplotlib figure.
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import copy
        from CRYSTALpytools.base.plotbase import _plot_label_preprocess

        ounit = self.unit
        self._set_unit('cm-1')
        xaxis = copy.deepcopy(self.frequency)
        if unit.lower() == 'nm': xaxis = 1 / xaxis * 1e7
        if unit.lower() not in ['nm', 'cm-1']:
            raise ValueError("Unknown unit: '{}'.".format(unit))

        # options
        valid_option = ['LG', 'V', 'RS', 'RE', 'REFL']
        if option.upper() not in valid_option:
            raise ValueError("Unknown option: '{}'.".format(option))
        option = option.upper()
        if self.type == 'molecule': option = 'LG'

        if option == 'LG':
            if self.absorbance.shape[0] < 1: raise Exception('Not a standard CRYSTAL IRSPEC.DAT data. Option not valid')
            yaxis = copy.deepcopy(self.absorbance[0])
        elif option == 'V':
            if self.absorbance.shape[0] < 2: raise Exception('Not a standard CRYSTAL IRSPEC.DAT data. Option not valid')
            yaxis = copy.deepcopy(self.absorbance[1])
        elif option == 'RS':
            if self.absorbance.shape[0] < 3: raise Exception('Not a standard CRYSTAL IRSPEC.DAT data. Option not valid')
            yaxis = copy.deepcopy(self.absorbance[2])
        elif option == 'RE':
            if self.absorbance.shape[0] < 4: raise Exception('Not a standard CRYSTAL IRSPEC.DAT data. Option not valid')
            yaxis = copy.deepcopy(self.absorbance[3])
        else:
            if len(self.reflectance) < 1: raise Exception('Reflectance spectra not found.')
            yaxis = copy.deepcopy(self.reflectance)

        yaxis = np.array(yaxis, ndmin=2)
        # normalize and shift
        if normalize == True:
            yaxis = yaxis / np.max(yaxis) * 100

        if REFL_overlap == True and option == 'REFL':
            for i in range(yaxis.shape[0]):
                yaxis[i] = yaxis[i] + shift*i

        if sharey == True:
            maxv = [np.max(yaxis) for i in range(yaxis.shape[0])]
        else:
            maxv = np.max(yaxis, axis=1)

        # plot
        if np.all(fig==None):
            if REFL_overlap != True and option == 'REFL':
                fig, ax = plt.subplots(yaxis.shape[0], 1, sharex=True,
                                       sharey=sharey, figsize=figsize, layout='tight')
            else:
                fig, ax = plt.subplots(1, 1, figsize=figsize)
            called = False
        else:
            called = True
        # note: Subplot is generated here rather than in 'plot.plot_IR', so
        # ax_index always start from 0
        ax_index = [i for i in range(len(fig.axes))]

        # plot commands
        if option == 'REFL':
            if np.all(label==None):
                label = ['# {:d}'.format(i+1) for i in range(yaxis.shape[0])]
            else:
                if isinstance(label, str):
                    label = ['{} {:d}'.format(label, i+1) for i in range(yaxis.shape[0])]
                else:
                    if len(label) != yaxis.shape[0]:
                        warnings.warn(
                            "Inconsistent lengths of number of plots and plot labels. Using default labels.",
                             stacklevel=2)
                        label = ['# {:d}'.format(i+1) for i in range(yaxis.shape[0])]

        doss = np.zeros([yaxis.shape[0], 1, 1]) # pseudo doss
        commands = _plot_label_preprocess(doss, label, color, linestyle, linewidth)

        if REFL_overlap == False and option == 'REFL':
            for iax in ax_index:
                ax = fig.axes[iax]
                ax.plot(xaxis, yaxis[iax], label=commands[0][iax][0], color=commands[1][iax][0],
                        linestyle=commands[2][iax][0], linewidth=commands[3][iax][0], **kwargs)
        else:
            for iy, y in enumerate(yaxis):
                ax = fig.axes[ax_index[0]]
                ax.plot(xaxis, y, label=commands[0][iy][0], color=commands[1][iy][0],
                        linestyle=commands[2][iy][0], linewidth=commands[3][iy][0], **kwargs)

        # range and labels
        if len(x_range) == 0: x_range = [np.min(xaxis), np.max(xaxis)]
        else: x_range = [np.min(x_range), np.max(x_range)]
        fig.axes[-1].set_xlim(x_range)
        if unit.lower() == 'cm-1': fig.axes[-1].xaxis.set_inverted(True)
        else: fig.axes[-1].set_xscale('log')

        for iiax, iax in enumerate(ax_index):
            fig.axes[iax].set_ylim([0, maxv[iiax]])
            _ = fig.axes[iax].get_yaxis().set_ticks([])

        if called == False:
            if unit.lower() == 'nm':
                fig.axes[-1].set_xlabel('Wavelength (nm)', fontsize=fontsize)
            else:
                fig.axes[-1].set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=fontsize)

            if option != 'REFL':
                fig.supylabel(r'Absorbance (arb. u.)', fontsize=fontsize)
            else:
                fig.supylabel(r'Reflectance (arb. u.)', fontsize=fontsize)

            if np.all(title!=None): fig.suptitle(title, fontsize=fontsize)
            if np.all(legend!=None) and np.all(commands[0][0][0]!=None):
                for i in ax_index: fig.axes[i].legend(loc=legend)

        self._set_unit(ounit)
        return fig

    def _set_unit(self, unit):
        """
        Set unit of frequency.

        Args:
            unit (str): 'cm-1' or 'THz'
        """
        from CRYSTALpytools.units import cm_to_thz, thz_to_cm

        if unit.lower() == self.unit.lower():
            return self

        opt_f_props = [] # Optional frequency properties
        if unit.lower() == 'thz':
            self.unit = 'THz'
            self.frequency = cm_to_thz(self.frequency)
            for p in opt_f_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, cm_to_thz(attrv))
        elif unit.lower() == 'cm-1':
            self.unit = 'cm-1'
            self.frequency = thz_to_cm(self.frequency)
            for p in opt_f_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, thz_to_cm(attrv))
        else:
            raise ValueError('Unknown unit.')
        return self


class Raman():
    """
    The class for Raman spectrum. Unit: in principle, should always in
    'cm :math:`^{-1}`'.

    Args:
        freq (array): Frequencies. Must be commensurate with unit.
        poly (array): Polycrystalline isotropic spectrum, 3\*nFreq array for
            total, parallel and perpendicular directions.
        single (array): Single crystal spectrum, 6\*nFreq array for xx, xy, xz,
            yy, yz, zz directions.
        unit (str): 'cm-1' or 'THz'.
    """
    def __init__(self, freq, poly, single, unit='cm-1'):
        import numpy as np

        self.frequency = np.array(freq, dtype=float)
        self.poly = np.array(poly, dtype=float)
        self.single = np.array(single, dtype=float)
        self.unit = unit

    @classmethod
    def from_file(self, specfile, output=None):
        """
        Generate an ``Raman`` object from output and 'RAMSPEC.DAT' files of CRYSTAL.

        Args:
            specfile (str): The 'RAMSPEC.DAT' file.

        Returns:
            cls (IR)
        """
        from CRYSTALpytools.crystal_io import Crystal_output
        return Crystal_output(output).get_spectra(specfile, 'RAMSPEC')

    def plot(self, option='poly', normalize=True, overlap=True,
             direction=['xx', 'xy', 'xz', 'yy', 'yz', 'zz'], shift=0,
             label=None, color=None, linestyle=None, linewidth=None, x_range=[],
             title=None, figsize=[6.4, 4.8], legend='upper left', sharey=True,
             fontsize=14, fig=None, ax_index=None, **kwargs):
        """
        Plot Raman spectra into the same axis.

        Args:
            option (str): 'tot', 'poly' or 'single'. 'tot' plots the total raman
                spectrum of polycrystalline material. 'poly' plots total, parallel
                and perpendicular spectra. 'single' plots spectra along xx, xy,
                xz, yy, yz, zz directions.
            normalize (bool): Normalize the maximum intensity to 100.
            overlap (bool): If more than 1 inequivalent directions exists,
                whether to plot spectra into the same plot or into subplots.
            direction (list|str): *``option='single'`` only* Specify the
                directions of single crystal spectra to plot.
            shift (float): If multiple spectra are plotted, shifting them up by
                the given value.
            label (list|str|None): List of plot labels. 'None' for the default
                values ('total', 'parallel' series or 'xx' 'yy' series) and
                string for prefix the string before default values. Otherwise
                should be a 1\*nPlot list.
            color (list|str|None): If str, use the same color for all the plot
                lines. If 1\*nPlot, use the color defined for every plot line.
                'None' for default values (matplotlib Tableau Palette).
            linestyle (list|str|None): See explanations of ``color``.
            linewidth (list|float|None): See explanations of ``color``.
            x_range (list): 1\*2 list of x axis range.
            title (str|None): The title of the plot. 'None' for no title.
            figsize (list): Matplotlib figure size.
            legend (str|None): Location of legend. None for not adding legend.
            sharey (bool): Whether to share the y-axis among subplots. Share x
                is enforced.
            fontsize (int): Fontsize of the axis label and title.
            fig (Figure): *Developer Only* Matplotlib figure object.
            ax_index (int): *Developer Only* Index of the Axes object in fig.
            \*\*kwargs: Other parameters passed to matplotlib ``Axes.plot()``
                method.

        Returns:
            fig (Figure): Matplotlib figure.
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import copy
        from CRYSTALpytools.base.plotbase import _plot_label_preprocess

        ounit = self.unit
        self._set_unit('cm-1')
        xaxis = copy.deepcopy(self.frequency)

        # options
        valid_option = ['tot', 'poly', 'single']
        if option.lower() not in valid_option:
            raise ValueError("Unknown option: '{}'.".format(option))
        option = option.lower()

        dirs = {'xx' : 0, 'xy' : 1, 'xz' : 2,
                'yx' : 1, 'yy' : 3, 'yz' : 4,
                'zx' : 2, 'zy' : 4, 'zz' : 5,}

        if option == 'tot':
            yaxis = copy.deepcopy(self.poly[0])
        elif option == 'poly':
            yaxis = copy.deepcopy(self.poly)
        elif option == 'single':
            direction = np.array(direction, ndmin=1)
            try: idxdir = [dirs[i.lower()] for i in direction]
            except KeyError:
                raise ValueError("Unknown direction name: '{}'.".format(i))
            yaxis = copy.deepcopy(self.single[idxdir])

        yaxis = np.array(yaxis, ndmin=2)

        # normalize and shift
        if normalize == True:
            yaxis = yaxis / np.max(yaxis) * 100
        for i in range(yaxis.shape[0]):
            yaxis[i] = yaxis[i] + shift*i

        if sharey == True:
            maxv = [np.max(yaxis) for i in range(yaxis.shape[0])]
        else:
            maxv = np.max(yaxis, axis=1)

        # plot
        if np.all(fig==None):
            if overlap == True:
                fig, ax = plt.subplots(1, 1, figsize=figsize)
            else:
                fig, ax = plt.subplots(yaxis.shape[0], 1, sharex=True, sharey=sharey,
                                       figsize=figsize, layout='tight')
            called = False
        else:
            called = True
        # note: Subplot is generated here rather than in 'plot.plot_Raman', so
        # ax_index always start from 0
        ax_index = [i for i in range(len(fig.axes))]

        # plot commands
        if option == 'poly':
            if np.all(label==None):
                label = ['total', 'parallel', 'perpendicular']
            else:
                if isinstance(label, str):
                    labeltmp = ['total', 'parallel', 'perpendicular']
                    label = ['{} {}'.format(label, labeltmp[i]) for i in range(3)]
                else:
                    if len(label) != 3:
                        warnings.warn("Inconsistent lengths of number of plots and plot labels. Using default labels.",
                                      stacklevel=2)
                        label = ['total', 'parallel', 'perpendicular']
        elif option == 'single':
            if np.all(label==None):
                label = direction.tolist()
            else:
                if isinstance(label, str):
                    labeltmp = direction.tolist()
                    label = ['{} {}'.format(label, labeltmp[i]) for i in range(6)]
                else:
                    if len(label) != len(direction):
                        warnings.warn("Inconsistent lengths of number of plots and plot labels. Using default labels.",
                                      stacklevel=2)
                        label = direction.tolist()
        doss = np.zeros([yaxis.shape[0], 1, 1]) # pseudo doss
        commands = _plot_label_preprocess(doss, label, color, linestyle, linewidth)

        if overlap == False:
            for iax in ax_index:
                ax = fig.axes[iax]
                ax.plot(xaxis, yaxis[iax], label=commands[0][iax][0], color=commands[1][iax][0],
                        linestyle=commands[2][iax][0], linewidth=commands[3][iax][0], **kwargs)
        else:
            for iy, y in enumerate(yaxis):
                ax = fig.axes[ax_index[0]]
                ax.plot(xaxis, y, label=commands[0][iy][0], color=commands[1][iy][0],
                        linestyle=commands[2][iy][0], linewidth=commands[3][iy][0], **kwargs)

        # range and labels
        if len(x_range) == 0: x_range = [np.min(xaxis), np.max(xaxis)]
        else: x_range = [np.min(x_range), np.max(x_range)]
        fig.axes[-1].set_xlim(x_range)
        fig.axes[-1].xaxis.set_inverted(True)

        for iiax, iax in enumerate(ax_index):
            fig.axes[iax].set_ylim([0, maxv[iiax]])
            _ = fig.axes[iax].get_yaxis().set_ticks([])

        if called == False:
            fig.axes[-1].set_xlabel(r'Wavenumber (cm$^{-1}$)', fontsize=fontsize)
            fig.supylabel(r'Intensity (arb. u.)', fontsize=fontsize)

            if np.all(title!=None): fig.suptitle(title, fontsize=fontsize)
            if np.all(legend!=None) and np.all(commands[0][0][0]!=None):
                for i in ax_index: fig.axes[i].legend(loc=legend)

        self._set_unit(ounit)
        return fig

    def _set_unit(self, unit):
        """
        Set unit of frequency.

        Args:
            unit (str): 'cm-1' or 'THz'
        """
        from CRYSTALpytools.units import cm_to_thz, thz_to_cm

        if unit.lower() == self.unit.lower():
            return self

        opt_f_props = [] # Optional frequency properties
        if unit.lower() == 'thz':
            self.unit = 'THz'
            self.frequency = cm_to_thz(self.frequency)
            for p in opt_f_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, cm_to_thz(attrv))
        elif unit.lower() == 'cm-1':
            self.unit = 'cm-1'
            self.frequency = thz_to_cm(self.frequency)
            for p in opt_f_props:
                if hasattr(self, p):
                    attrv = getattr(self, p)
                    setattr(self, p, thz_to_cm(attrv))
        else:
            raise ValueError('Unknown unit.')
        return self

