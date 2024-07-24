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
        self.carrier = np.array(carrier, dtype=float)
        self.data = np.array(tensor, dtype=float)
        self.type = type
        self.struc = struc
        self.unit = unit
        if len(self.T) != self.carrier.shape[0] or len(self.T) != self.data.shape[0]:
            raise ValueError("Inconsistent dimensionalities between temperature and input data.")
        if len(self.mu) != self.carrier.shape[1] or len(self.T) != self.data.shape[1]:
            raise ValueError("Inconsistent dimensionalities between chemical potential and input data.")
        self.spin = self.data.shape[-1]

    @classmethod
    def from_file(cls, boltztra_out, output=None):
        """
        Read electron transport properties by the BOLTZTRA keyword, including
        'KAPPA', 'SIGMA', 'SIGMAS', and 'SEEBECK'. Though currently the
        geometry information is not required, it is saved if the standard
        output file is given.

        Args:
            boltztra_out (str): 'DAT' files by CRYSTAL BOLTZTRA keyword.

        Returns:
            cls (Tensor)
        """
        from CRYSTALpytools.crystal_io import Properties_output

        return Properties_output(output).read_transport(file)

    def plot(self, x_axis='potential', x_range=[], temperature=[], direction='xx',
             data_label=None, data_color=None, data_linestyle=None, data_linewidth=None,
             zero_color='tab:gray', zero_linestyle='-', zero_linewidth=1.,
             add_title=True, figsize=[6.4, 4.8], legend='upper left', sharex=True,
             sharey=True, fontsize=14, fig=None, **kwargs):
        """
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from CRYSTALpytools.base.plotbase import _plot_label_preprocess


        # x axis
        if x_axis.lower() == 'potential':
            usepot = True
            if x_range == []:
                x_range = [np.min(self.mu), np.max(self.mu)]
        elif x_axis.lower() == 'carrier':
            usepot = False
            if x_range == []:
                x_range = [np.min(np.abs(self.carrier)), np.max(np.abs(self.carrier))]
        else:
            raise ValueError("Unknown x-axis style: '{}'.".format(x_axis))
        x_range = [np.min(x_range), np.max(x_range)]

        # temperature (projection)
        if temperature == []:
            temperature = self.T
        else:
            temperature = np.array(temperature, dtype=float, ndmin=1)
        nprj = len(temperature)

        # directions (subplots)
        if is isinstance(direction, str):
            direction = [direction]
        nplt = len(direction)

        # plot setups
        if np.all(label==None):
            label = ['{:>5.0f} K'.format(i) for i in temperature]
        elif isinstance(label, str):
            label = ['{} {:>5.0f} K'.format(label, i) for i in temperature]
        else: # length
            nlabel = len(label)
            label = [label[i%label] for i in range(nlabel)]
        ## get a pseudo band input
        bands = np.zeros([nplt, nprj, 1, 2], dtype=float)
        commands = _plot_label_preprocess(
            bands, data_label, data_color, data_linestyle, data_linewidth)

        # plotting
        if np.all(fig==None):
            fig, ax = plt.subfigures(len(direction), 1, sharex=True,
                                     figsize=figsize, layout='constrained')
        for idir, ax in enumerate(fig.axes):
            ax.hlines(0, x_range[0], x_range[1], colors=zero_color,
                      linestyle=zero_linestyle, linewidth=zero_linewidth)

        


        # assign plotting commands
        ## label
        
        ## color
        if np.all(color==None):
            clist = list(mcolors.TABLEAU_COLORS.keys())
            nclist = len(clist)
            color = [clist[i%nclist] for i in range(ntemp)]
        elif isinstance(color, str):
            color = [color for i in range(ntemp)]
        else: # length
            ncolor = len(color)
            color = [color[i%ncolor] for i in range(ncolor)]
        ## linestyle
        if isinstance(linestyle, str):
            linestyle = [linestyle, linestyle]
        else: # length
            linestyle = [linestyle[i%len(linestyle)] for i in range(2)]
            
            
        


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

