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

        if len(theta_range) != 0:
            ax.set_xlim([np.min(theta_range), np.max(theta_range)])
        ax.set_xlabel(r'2$\theta^{\circ}$', fontsize=fontsize)
        ax.set_ylabel(r'Intensity (arb. u.)', fontsize=fontsize)
        ax.set_ylim([0, maxv])
        _ = ax.get_yaxis().set_ticks([])
        if np.all(title!=None):
            ax.set_title(title, fontsize=fontsize)
        return fig

