#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/03/2022
"""


from matplotlib.pyplot import ylim
from pyparsing import counted_array


def plot_cry_bands(bands, k_labels=None, energy_range=None, title=False, not_scaled=False, mode='single', linestl='-',
                   linewidth=1, color='blue', fermi='forestgreen', k_range=None, labels=None, figsize=None, scheme=None,
                   sharex=True, sharey=True):

    import matplotlib.pyplot as plt
    # import matplotlib.lines as mlines
    import numpy as np
    import sys

    greek = {'Alpha': '\u0391', 'Beta': '\u0392', 'Gamma': '\u0393', 'Delta': '\u0394', 'Epsilon': '\u0395', 'Zeta': '\u0396', 'Eta': '\u0397',
             'Theta': '\u0398', 'Iota': '\u0399', 'Kappa': '\u039A', 'Lambda': '\u039B', 'Mu': '\u039C', 'Nu': '\u039D', 'Csi': '\u039E',
             'Omicron': '\u039F', 'Pi': '\u03A0', 'Rho': '\u03A1', 'Sigma': '\u03A3', 'Tau': '\u03A4', 'Upsilon': '\u03A5', 'Phi': '\u03A6',
             'Chi': '\u03A7', 'Psi': '\u03A8', 'Omega': '\u03A9', 'Sigma_1': '\u03A3\u2081'}

    # Error check on the mode flag
    modes = ['single', 'multi', 'compare', 'surface']
    if mode not in modes:
        print('Error,The selected mode '+mode+' is not among the possible ones: ' +
              modes[0]+', ' + modes[1] + ', '+modes[2] + ', or '+modes[3])
        sys.exit(1)

    # Error chenk on k_label
    if k_labels is not None:

        if isinstance(k_labels, list):
            for element in k_labels:
                if not isinstance(element, str):
                    print('Error, k_label must be a list of strings:' +
                          str(element)+' is not a string')
                    sys.exit(1)

        else:
            print('Error, k_labels must be a list of strings')
            sys.exit(1)

    # Error check on energy range
    if energy_range is not None:

        if isinstance(energy_range, list):

            if len(energy_range) > 2:
                print('Error, energy_range must be a list of two int or float (min,max)')
                sys.exit(1)

            for element in energy_range:
                if (not isinstance(element, int)) and (not isinstance(element, float)):
                    print('Error, energy_range must be a list of two int or float (min,max): ' +
                          str(element)+', is neither a float or an int')
                    sys.exit(1)

        else:
            print('Error, energy_range must be a list of two int or float (min,max)')
            sys.exit(1)

    # Error check on k_range
    if k_range is not None:

        if isinstance(k_range, list):
            if len(k_range) > 2:
                print('Error, k_range must be a list of two strings')
                sys.exit(1)

            for element in k_range:
                if not isinstance(element, str):
                    print('Error, k_label must be a list of two strings:' +
                          str(element)+' is not a string')
                    sys.exit(1)

        else:
            print('Error, k_range must be a list of two strings')
            sys.exit(1)

    # Error check on title
    if title is not False:
        if not isinstance(title, str):
            print('Error, title needs to be a string')
            sys.exit(1)

    if (mode == modes[0]) or (mode == modes[1]) or (mode == modes[3]):

        # plotting of a single band object
        if mode == modes[0]:

            dx = bands.k_point_plot

            pltband = bands.bands
            no_bands = np.shape(pltband)[0]
            ymin = np.amin(pltband)
            ymax = np.amax(pltband)
            xmin = min(dx)
            xmax = max(dx)
            count1 = 0
            count2 = 0

            # band plot
            for i in range(no_bands):

                if figsize is not None:
                    plt.figure()
                else:
                    plt.figure(figsize=figsize)

                if bands.spin == 1:
                    plt.plot(dx, pltband[i, :], color=color,
                             linestyle=linestl, linewidth=linewidth)

                elif bands.spin == 2:
                    if count1 == count2:
                        plt.plot(dx, pltband[i, :, 0], color='red',
                                 linestyle=linestl, linewidth=linewidth, label='Alpha')
                        plt.plot(dx, pltband[i, :, 1], color='black',
                                 linestyle=linestl, linewidth=linewidth, label='Beta')
                    else:
                        plt.plot(dx, pltband[i, :, 0], color='red',
                                 linestyle=linestl, linewidth=linewidth)
                        plt.plot(dx, pltband[i, :, 1], color='black',
                                 linestyle=linestl, linewidth=linewidth)
                    count1 += 1

        # plot of multiple band objects on a single plot
        elif mode == modes[1]:

            # Error check on the band on the 'multi' mode flag
            if not isinstance(bands, list):
                print('Error, When you choose a ' +
                      modes[1]+' plot bands needs to be a list of band objects')
                sys.exit(1)

            # Error check on color for the 'multi' mode flag
            if isinstance(color, list):
                if len(color) > len(bands):
                    print(
                        'Error, The number of colors is greater than the number of objects you want to plot')
                    sys.exit(1)

            else:
                color = ['dimgrey', 'blue', 'indigo', 'slateblue',
                         'thistle', 'purple', 'orchid', 'crimson']

            # Warning comparison with band.spin==2
            for m in bands:
                if m.spin == 2:
                    print(
                        "Warning: the 'multi' plot is not fully implemented at the moment for file with NSPIN = 2")

            # scaling that enables the comparison of band structure calculated at different pressures
            if not_scaled is False:
                reference = xmax = np.amax(bands[0].k_point_plot)
                xmin = np.amin(bands[0].k_point_plot)

            else:
                xmax = []
                xmin = []

            ymin = []
            ymax = []

            if figsize is not None:
                plt.figure(figsize=figsize)

            # plot of all the bands obj present in the list
            for index, data in enumerate(bands):
                # scaling that enables the comparison of band structure calculated at different pressures
                if not_scaled is False:
                    k_max = np.amax(data.k_point_plot)
                    dx = (data.k_point_plot/k_max)*reference
                else:
                    dx = data.k_point_plot
                    xmin.append(np.amin(dx))
                    xmax.append(np.amax(dx))

                pltband = data.bands
                no_bands = np.shape(pltband)[0]
                ymin.append(np.amin(pltband))
                ymax.append(np.amax(pltband))

                count1 = 0
                count2 = 0

                # Effective plot
                for j in range(no_bands):

                    if (count1 == count2) and (labels is not None):
                        if data.spin == 1:
                            if isinstance(linestl, list):
                                plt.plot(dx, pltband[j, :], color=color[index],
                                         linestyle=linestl[index], linewidth=linewidth, label=labels[index])
                            else:
                                plt.plot(dx, pltband[j, :], color=color[index],
                                         linestyle=linestl, linewidth=linewidth, label=labels[index])
                        elif data.spin == 2:
                            plt.plot(dx, pltband[j, :, 0], color=color[index],
                                     linestyle=linestl, linewidth=linewidth, label=labels[index]+' Alpha')
                            plt.plot(dx, pltband[j, :, 1], color=color[index],
                                     linestyle='--', linewidth=linewidth, label=labels[index]+' Beta')

                    else:
                        if data.spin == 1:
                            if isinstance(linestl, list):
                                plt.plot(dx, pltband[j, :], color=color[index],
                                         linestyle=linestl[index], linewidth=linewidth)
                            else:
                                plt.plot(
                                    dx, pltband[j, :], color=color[index], linestyle=linestl, linewidth=linewidth)
                        elif data.spin == 2:
                            plt.plot(
                                dx, pltband[j, :, 0], color=color[index], linestyle=linestl, linewidth=linewidth)
                            plt.plot(
                                dx, pltband[j, :, 1], color=color[index], linestyle='--', linewidth=linewidth)

                    count1 += 1

            if (isinstance(xmin, list)) or (isinstance(xmax, list)):
                xmin = min(xmin)
                xmax = max(xmax)

            ymin = min(ymin)
            ymax = max(ymax)

        if mode == modes[3]:
            print('Warning, the surface bands is not ready yet')
            sys.exit(0)

        # HSP line plot
        if isinstance(bands, list):
            hsp = bands[0].tick_position
        else:
            hsp = bands.tick_position

        if k_labels is not None:
            if len(hsp) != len(k_labels):
                if len(hsp) > len(k_labels):
                    print(
                        'Error, you have specified a number of label smalle than the number of High Simmetry Point along the path')
                elif len(hsp) < len(k_labels):
                    print(
                        'Error, you have more labels than the High Simmetry point along the path')
                sys.exit(1)

        hsp[len(hsp)-1] = xmax

        y_band = np.linspace(ymin-3, ymax+3, 2)

        for j in hsp:
            x_band = np.ones(2)*j
            plt.plot(x_band, y_band, color='black', linewidth=0.5)

        # plot of the fermi level
        x = np.linspace(xmin, xmax, 2)
        y = np.zeros(2)
        plt.plot(x, y, color=fermi, linewidth=2.5)

        # definition of the plot title
        if title is not False:
            plt.title(title)

        if not isinstance(bands, list):
            if bands.spin == 2:
                plt.legend()
        else:
            plt.legend()

        plt.ylabel('$E-E_F$ (eV)')

    # compare mode plot
    elif mode == modes[2]:

        # Error check on type(bands)
        if not isinstance(bands, list):
            print('Error, When you choose a ' +
                  modes[2]+' plot bands needs to be a list of band objects')
            sys.exit(1)

        # Error check on sharex and sharey option
        if (not isinstance(sharex, bool)) or (not isinstance(sharey, bool)):
            accepted_str = ['row', 'col']

            if (isinstance(sharex, str)) or (isinstance(sharey, str)):
                if sharex not in accepted_str:
                    print('Error, sharex can only be equal to row or col')
                    sys.exit(1)
                elif sharey not in accepted_str:
                    print('Error, sharey can only be equal to row or col')
                    sys.exit(1)

            else:
                print('Error, sharex and sharey have to be boolean')
                sys.exit(1)

        # Error check and definition of scheme
        if scheme is None:
            n_rows = 1
            n_col = len(bands)
        else:
            if (not isinstance(scheme, tuple)) and (not isinstance(scheme, list)):
                print('Error, scheme needs to be a tuple or a list')
                sys.exit(1)
            elif len(scheme) > 2:
                print(
                    'Error, scheme needs to be a tuple or a list of two elements (nrow,ncol)')
                sys.exit(1)
            else:
                n_rows = scheme[0]
                n_col = scheme[1]
        # Creation of the subplots
        if figsize is None:
            fig, axs = plt.subplots(nrows=n_rows, ncols=n_col,
                                    sharex=sharex, sharey=sharey)
        else:
            fig, axs = plt.subplots(nrows=n_rows, ncols=n_col,
                                    sharex=sharex, sharey=sharey)
        # Scaling with different size of the same brillouin zone
        if not_scaled is False:
            reference = xmax = np.amax(bands[0].k_point_plot)
            xmin = np.amin(bands[0].k_point_plot)

        else:
            xmax = []
            xmin = []

        ymin = []
        ymax = []
        count3 = 0

        # Plot of the different band structure into the subplots
        for col in range(n_col):
            for row in range(n_rows):
                data = bands[count3]
                if count3 == 0:
                    hsp = data.tick_position
                pltband = data.bands
                no_bands = np.shape(pltband)[0]
                if not_scaled is False:
                    k_max = np.amax(data.k_point_plot)
                    dx = (data.k_point_plot/k_max)*reference
                else:
                    dx = data.k_point_plot
                    xmin.append(np.amin(dx))
                    xmax.append(np.amax(dx))
                ymin.append(np.amin(pltband))
                ymax.append(np.amax(pltband))
                # Effective plotting action
                for j in range(no_bands):
                    if data.spin == 1:
                        if n_rows == 1:
                            axs[col].plot(
                                dx, pltband[j, :], color=color, linestyle=linestl, linewidth=linewidth)
                        else:
                            axs[row, col].plot(
                                dx, pltband[j, :], color=color, linestyle=linestl, linewidth=linewidth)
                    elif data.spin == 2:
                        if n_rows == 1:
                            axs[col].plot(dx, pltband[j, :, 0], color=color,
                                          linestyle='-', linewidth=linewidth, label='Alpha')
                            axs[col].plot(dx, pltband[j, :, 0], color=color,
                                          linestyle='--', linewidth=linewidth, label='Beta')
                        else:
                            axs[row, col].plot(
                                dx, pltband[j, :, 0], color=color, linestyle='-', linewidth=linewidth, label='Alpha')
                            axs[row, col].plot(
                                dx, pltband[j, :, 0], color=color, linestyle='--', linewidth=linewidth, label='Beta')

                # Plot of the HSPs lines
                yhsp = np.linspace(np.amin(pltband)+5, np.amax(pltband)+5, 2)
                for j in hsp:
                    xhsp = np.ones(2)*j
                    if n_rows == 1:
                        axs[col].plot(
                            xhsp, yhsp, color='black', linewidth=0.5)
                    else:
                        axs[row, col].plot(
                            xhsp, yhsp, color='black', linewidth=0.5)

                # Fermi level line plot
                xfermi = np.linspace(np.amin(pltband), np.amax(pltband), 2)
                yfermi = np.zeros(2)
                if n_rows == 1:
                    axs[col].plot(
                        xfermi, yfermi, color=fermi, linewidth=2.5)
                else:
                    axs[row, col].plot(
                        xfermi, yfermi, color=fermi, linewidth=2.5)

                # Definition of x and y limits
                if n_rows == 1:
                    if sharex is not True:
                        """hsp_label = []
                        for element in k_labels:
                            if element in k_labels:
                                g = greek.get(element)
                                hsp_label.append(g)
                        axs[col].set_xticks(hsp)
                        if k_labels is not None:
                            axs[col].set_xlabels(hsp_label)"""
                        print(
                            'Warning, the sharex = False option has not been developed yet')
                    axs[col].set_xlim([np.amin(dx), np.amax(dx)])
                    if (sharey is not True) and (energy_range is not None):
                        axs[col].set_ylim([energy_range[0], energy_range[1]])
                    else:
                        axs[col].set_ylim([np.amin(pltband), np.amax(pltband)])
                else:
                    if sharex is not True:
                        """hsp_label = []
                        for element in k_labels:
                            if element in k_labels:
                                g = greek.get(element)
                                hsp_label.append(g)
                        axs[row, col].set_xticks(hsp)
                        if k_labels is not None:
                            axs[row, col].set_xlabels(hsp_label)"""
                        print(
                            'Warning, the sharex = False option has not been developed yet')
                    axs[row, col].set_xlim([np.amin(dx), np.amax(dx)])
                    if (sharey is not True) and (energy_range is not False):
                        axs[row, col].set_ylim(
                            [energy_range[0], energy_range[1]])
                    else:
                        axs[row, col].set_ylim(
                            [np.amin(pltband), np.amax(pltband)])

                count3 += 1

        fig.text(.06, 0.5, '$E-E_F$ (eV)', ha='center',
                 va='center', rotation='vertical')

        if (isinstance(ylim, list)) or (isinstance(ymax, list)):
            ymin = min(ymin)
            ymax = max(ymax)

    hsp_label = []
    high_sym_point = []
    if k_labels is not None:
        for n in k_labels:
            if n in greek:
                g = greek.get(n)
                hsp_label.append(g)
                high_sym_point.append(n)
            else:
                hsp_label.append(n)
                high_sym_point.append(n)

    # give the possibility through the k_range to select a shorter path than the one calculated
    high_sym_point2 = high_sym_point
    count = 0
    for i in high_sym_point2:
        repeat = high_sym_point2.count(i)
        if repeat != 1:
            for p in range(0, len(high_sym_point2)):
                if p != count:
                    repeat_count = 1
                    q = high_sym_point2[p]
                    r = high_sym_point[p]
                    if (q == i) & (q == r):
                        high_sym_point2[p] = i+str(repeat_count)
                        repeat_count += 1
                        if repeat_count > repeat:
                            repeat_count = 0
        count += 1

    path_dict = dict(zip(high_sym_point2, hsp))

    # definition of the ylim
    if (energy_range is not None) or (sharey is True):
        ymin = energy_range[0]
        ymax = energy_range[1]

    # definition of the xlim
    if k_range is not None:
        xmin = path_dict[k_range[0]]
        xmax = path_dict[k_range[1]]

    if k_labels is not None:
        plt.xticks(hsp, hsp_label)
    else:
        plt.xticks(hsp)

    plt.ylim(ymin, ymax)
    if (mode == modes[0]) or (mode == modes[1]):
        plt.xlim(xmin, xmax)
    elif (mode == modes[2]) and (k_range is not None):
        print('Warning, the k_range is not available yet for the compare mode')

    plt.show()

    """kpoints = bands.tick_position 
    efermi_band = bands.efermi    
                
    x1 = bands.bands[:,0,:]
    print(x1)
    y1 = bands.bands[:,1:,:]
    
    fig, ax1 = plt.subplots(1,1)
    
    if len(bands.bands[0,0,:]) > 1:  
        
        ax1.plot(x1[:,1], y1[:,:,1],'-',  c='red' )
        ax1.plot(x1[:,0], y1[:,:,0],'-' , c='black')        
        # ax1.legend( bbox_to_anchor=(0.9, .75))
        spin_band = [mlines.Line2D([], [], color='black', label='Alpha'),
                     mlines.Line2D([], [], color='red', label='Beta')]
        ax1.legend(spin_band,['Alpha electrons', 'Beta electrons'],facecolor='white', framealpha=1,bbox_to_anchor=(.83,.90))
    
    else:
        ax1.plot(x1[:,0], y1[:,:,0],'-', c='black' )  
    
    # Display E Fermi on band structure
    if not_scaled == True:
        fermi_line = [[0, np.amax(bands.bands[:,1:,0])+1],[efermi_band,efermi_band]] 
    else:
        fermi_line = [[0, np.amax(bands.bands[:,1:,0])+1],[0,0]] 
        
    ax1.plot(fermi_line[0],fermi_line[1], '--',color='grey')
    ax1.set_title('Band structure', size = 18)
    ax1.set_xlabel('k point', size =12)
    ax1.set_ylabel('E-E Fermi (eV)', size =12)
    ax1.set_xticks(kpoints)
    if k_labels is not None:
        ax1.set_xticklabels(k_labels, size=12)
    ax1.set_xlim([0, bands.bands[-1,0,0]])
    ax1.set_ylim(energy_range)
    ax1.grid()
    
        
    fig.set_size_inches(8, 8)
    if title != False:
        fig.suptitle(title, size=22)
        plt.subplots_adjust(wspace=0.2, top=0.88)"""


def plot_cry_multibands(bands_list, k_labels=None, energy_range=None, title=False, not_scaled=False):
    # Filippo's function
    # bands list is a Crystal_band object
    pass


def plot_cry_bands_compare():
    pass


def plot_cry_doss():
    pass


def plot_cry_es():
    pass


def plot_contour(contour_obj, diff=False):
    # Ale C's function
    pass
