#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/03/2022
"""



def plot_cry_bands(bands, k_labels=None, energy_range=None, title=False, not_scaled=False, save_to_file = False, mode='single', linestl='-',
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
            if labels is not None:
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
                                    sharex=sharex, sharey=sharey, figsize=figsize)
        else:
            fig, axs = plt.subplots(nrows=n_rows, ncols=n_col,
                                    sharex=sharex, sharey=sharey, figsize=figsize)
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

    if save_to_file != False:
        save_plot(save_to_file)
    
    plt.show()


def plot_cry_doss(doss, color='blue', fermi: str = 'forestgreen', save_to_file = False, overlap: bool = False, labels=None, figsize=None, linestl=None,
                  linewidth=1, title: str = None, beta: str = 'up', energy_range=None, dos_range=None, prj: list = None):

    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    from os import path

    # Error check on beta
    accepted_beta = ['up', 'down']
    if beta not in accepted_beta:
        print('Error!, beta can be defined only as: ' +
              accepted_beta[0] + ' or ' + accepted_beta[1])
        sys.exit(1)

    # Plot for overlap == True
    if overlap == True:

        # Error check on color for color in overlap == true and default list of color
        if (not isinstance(color, list)) and (doss.n_proj > 1):
            print('Warning!, When overlap is true color has to be a list!')
            color = ['blue', 'indigo', 'midgrey', 'crimson']
            if (len(color) < doss.n_proj) or (len(color) < len(prj)):
                print('Error!, When overlap is true color has to be a list!')
                sys.exit(1)

        # Error check on labels for overlap == true
        if labels is not None:
            if not isinstance(labels, list) and (doss.n_proj > 1):
                print('Error, When overlap is true labels has to be a list')
                sys.exit(1)

        # Default for figsize
        if figsize is None:
            figsize = (12, 5)

        # Error check on color, labels, and linestl
        if (isinstance(color, list)) or (isinstance(labels, list)) or (isinstance(linestl, list)):
            if color is not None:
                if prj is None:
                    if (len(color) > doss.n_proj) or (len(color) < doss.n_proj):
                        if len(color) > doss.n_proj:
                            print(
                                'Error!, the number of colors is greater than the number of projections!')
                            sys.exit(1)
                        elif len(color) < doss.n_proj:
                            print(
                                'Error!, the number of colors is greater than the number of projections!')
                            sys.exit(1)
                else:
                    if (len(color) > len(prj)) or (len(color) < len(prj)):
                        if len(color) > len(prj):
                            print(
                                'Error!, the number of colors is greater than the number of projections required!')
                            sys.exit(1)
                        elif len(color) < len(prj):
                            print(
                                'Error!, the number of colors is greater than the number of projections required!')
                            sys.exit(1)

            if labels is not None:
                if prj is None:
                    if (len(labels) > doss.n_proj) or (len(labels) < doss.n_proj):
                        if len(labels) > doss.n_proj:
                            print(
                                'Error!, the number of labels is greater than the number of projections!')
                            sys.exit(1)
                        elif len(labels) < doss.n_proj:
                            print(
                                'Error!, the number of labels is greater than the number of projections!')
                            sys.exit(1)
                else:
                    if (len(labels) > len(prj)) or (len(labels) < len(prj)):
                        if len(labels) > len(prj):
                            print(
                                'Error!, the number of labels is greater than the number of projections required!')
                            sys.exit(1)
                        elif len(color) < len(prj):
                            print(
                                'Error!, the number of labels is greater than the number of projections required!')
                            sys.exit(1)

            if linestl is not None:
                if prj is None:
                    if (len(linestl) > doss.n_proj) or (len(linestl) < doss.n_proj):
                        if len(linestl) > doss.n_proj:
                            print(
                                'Error!, the number of linestl is greater than the number of projections!')
                            sys.exit(1)
                        elif len(linestl) < doss.n_proj:
                            print(
                                'Error!, the number of linestl is greater than the number of projections!')
                            sys.exit(1)
                else:
                    if (len(linestl) > len(prj)) or (len(linestl) < len(prj)):
                        if len(linestl) > len(prj):
                            print(
                                'Error!, the number of linestl is greater than the number of projections required!')
                            sys.exit(1)
                        elif len(color) < len(prj):
                            print(
                                'Error!, the number of linestl is greater than the number of projections required!')
                            sys.exit(1)

        # Creation of the figure
        plt.figure(figsize=figsize)

        # Creation of dx for the plot
        if doss.spin == 1:
            dx = doss.doss[:, 0]
        elif doss.spin == 2:
            dx = doss.doss[:, 0, :]
            dx_alpha = doss.doss[:, 0, 0]
            dx_beta = doss.doss[:, 0, 1]

        # Determination of xmin, xmax, ymin, and ymax
        xmin = np.amin(dx)
        xmax = np.amax(dx)
        ymin = np.amin(doss.doss[1:, :, :])
        ymax = np.amax(doss.doss[1:, :, :])

        # Plot of all projections
        if prj is None:
            for projection in range(1, doss.n_proj+1):
                # for projection in range(1, 2):
                if doss.spin == 1:
                    if doss.n_proj > 1:
                        if linestl is None:
                            plt.plot(dx, doss.doss[:, projection], color=color[projection-1],
                                     label=labels[projection-1], linewidth=linewidth)
                        else:
                            plt.plot(dx, doss.doss[:, projection], color=color[projection-1],
                                     label=labels[projection-1], linestyle=linestl[projection-1], linewidth=linewidth)
                    else:
                        plt.plot(
                            dx, doss.doss[:, projection], color=color, linewidth=linewidth)
                elif doss.spin == 2:
                    if beta == accepted_beta[0]:
                        if doss.n_proj > 1:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0], color=color[projection-1],
                                     label=labels[projection-1], linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, -doss.doss[:, projection, 1], color=color[projection-1],
                                     label=labels[projection-1], linestyle='--', linewidth=linewidth)
                        else:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0],
                                     linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, -doss.doss[:, projection, 1],
                                     linestyle='--', linewidth=linewidth)
                    elif beta == accepted_beta[1]:
                        if doss.n_proj > 1:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0], color=color[projection-1],
                                     label=labels[projection-1], linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, doss.doss[:, projection, 1], color=color[projection-1],
                                     label=labels[projection-1], linestyle='--', linewidth=linewidth)

                        else:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0], color=color,
                                     linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, doss.doss[:, projection, 1], color=color,
                                     linestyle='--', linewidth=linewidth)

        # Plot of selected projections
        else:
            for index, projection in enumerate(prj):
                if doss.spin == 1:
                    if doss.n_proj > 1:
                        if linestl is None:
                            plt.plot(dx, doss.doss[:, projection], color=color[index],
                                     label=labels[index], linewidth=linewidth)
                        else:
                            plt.plot(dx, doss.doss[:, projection], color=color[index],
                                     label=labels[index], linestyle=linestl[index], linewidth=linewidth)
                    else:
                        plt.plot(
                            dx, doss.doss[:, projection], color=color, linewidth=linewidth)
                elif doss.spin == 2:
                    if beta == accepted_beta[0]:
                        if doss.n_proj > 1:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0], color=color[index],
                                     label=labels[index], linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, -doss.doss[:, projection, 1], color=color[index],
                                     label=labels[index], linestyle='--', linewidth=linewidth)
                        else:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0],
                                     linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, -doss.doss[:, projection, 1],
                                     linestyle='--', linewidth=linewidth)
                    elif beta == accepted_beta[1]:
                        if doss.n_proj > 1:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0], color=color[index],
                                     label=labels[index], linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, doss.doss[:, projection, 1], color=color[index],
                                     label=labels[index], linestyle='--', linewidth=linewidth)
                        else:
                            plt.plot(dx_alpha, doss.doss[:, projection, 0], color=color,
                                     linestyle='-', linewidth=linewidth)
                            plt.plot(dx_beta, doss.doss[:, projection, 1], color=color,
                                     linestyle='--', linewidth=linewidth)

        yfermi_level = np.linspace(ymin-6, ymax+6, 2)
        xfermi_level = np.zeros(2)
        if beta == accepted_beta[0]:
            plt.plot(xfermi_level, yfermi_level,
                     color=fermi, linewidth=1.5)
        else:
            yfermi_level = np.linspace(-ymax-5, ymax+5, 2)
            plt.plot(xfermi_level, yfermi_level,
                     color=fermi, linewidth=1.5)

        y_zero = np.zeros(2)
        x_zero = np.linspace(xmin, xmax, 2)
        plt.plot(x_zero, y_zero, color='black', linewidth=0.4)

        if dos_range is not None:
            ymin = dos_range[0]
            ymax = dos_range[1]
            plt.ylim(ymin, ymax)
        else:
            if doss.spin == 1:
                ymin = 0
                plt.ylim(ymin, ymax+5)
            elif doss.spin == 2:
                if beta == accepted_beta[0]:
                    ymin = 0
                    plt.ylim(ymin, ymax+5)
                elif beta == accepted_beta[1]:
                    plt.ylim(ymin-5, ymax+5)

        plt.ylabel('DOS (a.u)')

        if title is not None:
            plt.title(title)

        if labels is not None:
            plt.legend()

    # Plot for overlap == False
    else:

        # Error checks on colors, labels, and linestl
        if isinstance(color, list):
            print('Warning!, When overlap is false color should be a string!')
            color = 'blue'

        if isinstance(labels, list):
            print('Warning!, When overlap is false labels should be a string!')
            labels = None

        if isinstance(linestl, list):
            print('Warning!, When overlap is false color should be a string!')
            linestl = '-'

        # Creation of dx for the plot
        if doss.spin == 1:
            dx = doss.doss[:, 0]
        elif doss.spin == 2:
            dx = doss.doss[:, 0, :]
            dx_alpha = doss.doss[:, 0, 0]
            dx_beta = doss.doss[:, 0, 1]

        # Determination of xmin and xmax
        xmin = np.amin(dx)
        xmax = np.amax(dx)

        # Creation xfermi, x_zero, and y_zero
        xfermi = np.zeros(2)
        x_zero = np.linspace(xmin, xmax, 2)
        y_zero = np.zeros(2)

        # Creation of subplots for all projections
        if prj is None:
            fig, axs = plt.subplots(
                nrows=doss.n_proj, ncols=1, sharex=True, figsize=figsize)

        # Creation of subplots for selected projections
        else:
            fig, axs = plt.subplots(
                nrows=len(prj), ncols=1, sharex=True, figsize=figsize)

        # Plot for all projections
        if prj is None:
            for projection in range(doss.n_proj):
                if doss.spin == 1:
                    ymin = 0
                    ymax = np.amax(doss.doss[:, projection+1])
                    yfermi = np.linspace(ymin, ymax, 2)
                    axs[projection].plot(dx, doss.doss[:, projection+1],
                                         color=color, linestyle=linestl, linewidth=linewidth)
                    axs[projection].plot(
                        xfermi, yfermi, color=fermi, linewidth=1.5)
                    if dos_range is not None:
                        ymin = dos_range[0]
                        ymax = dos_range[1]
                    axs[projection].set_ylim(ymin, ymax)
                elif doss.spin == 2:
                    if beta == accepted_beta[0]:
                        ymin = 0
                        ymax = np.amax(doss.doss[:, projection+1, :])
                        yfermi = np.linspace(ymin, ymax, 2)
                        axs[projection].plot(dx_alpha, doss.doss[:, projection+1, 0], color=color,
                                             linestyle='-', linewidth=linewidth, label='Alpha')
                        axs[projection].plot(dx_beta, -doss.doss[:, projection+1, 1], color=color,
                                             linestyle='--', linewidth=linewidth, label='Beta')
                        axs[projection].plot(
                            xfermi, yfermi, color=fermi, linewidth=1.5)
                        if dos_range is not None:
                            ymin = dos_range[0]
                            ymax = dos_range[1]
                        axs[projection].set_ylim(ymin, ymax)
                    elif beta == accepted_beta[1]:
                        ymin = -np.amax(doss.doss[:, projection+1, :])
                        ymax = np.amax(doss.doss[:, projection+1, :])
                        yfermi = np.linspace(ymin, ymax, 2)
                        axs[projection].plot(dx_alpha, doss.doss[:, projection+1, 0], color=color,
                                             linestyle='-', linewidth=linewidth, label='Alpha')
                        axs[projection].plot(dx_beta, doss.doss[:, projection+1, 1], color=color,
                                             linestyle='--', linewidth=linewidth, label='Beta')
                        axs[projection].plot(x_zero, y_zero, linewidth=0.4)
                        axs[projection].plot(
                            xfermi, yfermi, color=fermi, linewidth=2.5)
                        if dos_range is not None:
                            ymin = dos_range[0]
                            ymax = dos_range[1]
                        axs[projection].set_ylim(ymin, ymax)
        # Plot for selected projections
        else:
            for index, projection in enumerate(prj):
                if doss.spin == 1:
                    ymin = 0
                    ymax = np.amax(doss.doss[:, projection])
                    yfermi = np.linspace(ymin, ymax, 2)
                    axs[index].plot(dx, doss.doss[:, projection],
                                    color=color, linestyle=linestl, linewidth=linewidth)
                    axs[index].plot(
                        xfermi, yfermi, color=fermi, linewidth=1.5)
                    if dos_range is not None:
                        ymin = dos_range[0]
                        ymax = dos_range[1]
                    axs[index].set_ylim(ymin, ymax)
                elif doss.spin == 2:
                    if beta == accepted_beta[0]:
                        ymin = 0
                        ymax = np.amax(doss.doss[:, projection, :])
                        yfermi = np.linspace(ymin, ymax, 2)
                        axs[index].plot(dx_alpha, doss.doss[:, projection, 0], color=color,
                                        linestyle='-', linewidth=linewidth, label='Alpha')
                        axs[index].plot(dx_beta, -doss.doss[:, projection, 1], color=color,
                                        linestyle='--', linewidth=linewidth, label='Beta')
                        axs[index].plot(
                            xfermi, yfermi, color=fermi, linewidth=1.5)
                        if dos_range is not None:
                            ymin = dos_range[0]
                            ymax = dos_range[1]
                        axs[index].set_ylim(ymin, ymax)
                    elif beta == accepted_beta[1]:
                        ymin = -np.amax(doss.doss[:, projection, :])
                        ymax = np.amax(doss.doss[:, projection, :])
                        yfermi = np.linspace(ymin, ymax, 2)
                        axs[index].plot(dx_alpha, doss.doss[:, projection, 0], color=color,
                                        linestyle='-', linewidth=linewidth, label='Alpha')
                        axs[index].plot(dx_beta, doss.doss[:, projection, 1], color=color,
                                        linestyle='--', linewidth=linewidth, label='Beta')
                        axs[index].plot(x_zero, y_zero, linewidth=0.4)
                        axs[index].plot(
                            xfermi, yfermi, color=fermi, linewidth=1.5)
                        if dos_range is not None:
                            ymin = dos_range[0]
                            ymax = dos_range[1]
                        axs[index].set_ylim(ymin, ymax)

        if (doss.spin == 2) and (beta == accepted_beta[1]):
            plt.legend()

        fig.text(0.06, 0.5, 'DOS (a.u.)', ha='center',
                 va='center', rotation='vertical')

    if energy_range is not None:
        xmin = energy_range[0]
        xmax = energy_range[1]
    plt.xlim(xmin, xmax)

    plt.xlabel('Energy (eV)')

    if save_to_file != False:       
        folder = path.split(save_to_file)[0]
        if path.exists(folder) == True:
            plt.savefig('%s.png'%save_to_file)
        else:
            print('Error: folder %s does not exist'%save_to_file)
            sys.exit(1)
    
    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()




def plot_cry_es(bands, doss, k_labels: list = None, save_to_file = False, color_bd='blue', color_doss='blue', fermi='forestgreen', energy_range: list = None, linestl_bd='-',
                linestl_doss=None, linewidth=1, prj: list = None, figsize=None, labels: list = None, dos_max_range: float = None):

    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    from os import path

    greek = {'Alpha': '\u0391', 'Beta': '\u0392', 'Gamma': '\u0393', 'Delta': '\u0394', 'Epsilon': '\u0395', 'Zeta': '\u0396', 'Eta': '\u0397',
             'Theta': '\u0398', 'Iota': '\u0399', 'Kappa': '\u039A', 'Lambda': '\u039B', 'Mu': '\u039C', 'Nu': '\u039D', 'Csi': '\u039E',
             'Omicron': '\u039F', 'Pi': '\u03A0', 'Rho': '\u03A1', 'Sigma': '\u03A3', 'Tau': '\u03A4', 'Upsilon': '\u03A5', 'Phi': '\u03A6',
             'Chi': '\u03A7', 'Psi': '\u03A8', 'Omega': '\u03A9', 'Sigma_1': '\u03A3\u2081'}

    fig, axs = plt.subplots(nrows=1, ncols=2,
                            gridspec_kw={'width_ratios': [2, 1]},
                            sharex=False, sharey=True, figsize=figsize,
                            )
    hsp = bands.tick_position

    if k_labels is not None:
        if len(hsp) != len(k_labels):
            if len(hsp) > len(k_labels):
                print(
                    'Error!, you have specified a number of label smalle than the number of High Simmetry Point along the path')
            elif len(hsp) < len(k_labels):
                print(
                    'Error!, you have more labels than the High Simmetry point along the path')
            sys.exit(1)

    dx_bd = bands.k_point_plot

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
            axs[0].plot(dx_bd, pltband[i, :], color=color_bd,
                        linestyle=linestl_bd, linewidth=linewidth)

        elif bands.spin == 2:
            if count1 == count2:
                axs[0].plot(dx_bd, pltband[i, :, 0], color=color_bd,
                            linestyle='-', linewidth=linewidth, label='Alpha')
                axs[0].plot(dx_bd, pltband[i, :, 1], color=color_bd,
                            linestyle='--', linewidth=linewidth, label='Beta')
            else:
                axs[0].plot(dx_bd, pltband[i, :, 0], color=color_bd,
                            linestyle='-', linewidth=linewidth)
                axs[0].plot(dx_bd, pltband[i, :, 1], color=color_bd,
                            linestyle='--', linewidth=linewidth)

        count1 += 1

    if doss.spin == 1:
        dx_dos = doss.doss[:, 0]
    elif doss.spin == 2:
        dx_dos = doss.doss[:, 0, :]
        dx_alpha = doss.doss[:, 0, 0]
        dx_beta = doss.doss[:, 0, 1]

     # Determination of xmin, xmax, ymin, and ymax
    xmin_dos = np.amin(doss.doss[:, 1:, :])
    xmax_dos = np.amax(doss.doss[:, 1:, :])
    ymin_dos = np.amin(dx_dos)
    ymax_dos = np.amax(dx_dos)

    # Plot of all projections
    if prj is None:
        for projection in range(1, doss.n_proj+1):
            if doss.spin == 1:

                if doss.n_proj > 1:

                    if linestl_doss is None:
                        axs[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[projection-1],
                                    label=labels[projection-1], linewidth=linewidth)

                    else:
                        axs[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[projection-1],
                                    label=labels[projection-1], linestyle=linestl_doss[projection-1], linewidth=linewidth)

                else:
                    axs[1].plot(doss.doss[:, projection], dx_dos,
                                color=color_doss, linewidth=linewidth)

            elif doss.spin == 2:
                if doss.n_proj > 1:
                    axs[1].plot(doss.doss[:, projection, 0], dx_alpha, color=color_doss[projection-1],
                                label=labels[projection-1], linestyle='-', linewidth=linewidth)
                    axs[1].plot(-doss.doss[:, projection, 1], dx_beta, color=color_doss[projection-1],
                                label=labels[projection-1], linestyle='--', linewidth=linewidth)
                else:
                    axs[1].plot(doss.doss[:, projection, 0], dx_alpha, color=color_doss,
                                linestyle='-', linewidth=linewidth)
                    axs[1].plot(-doss.doss[:, projection, 1], dx_beta, color=color_doss,
                                linestyle='--', linewidth=linewidth)

    else:
        for index, projection in enumerate(prj):
            if doss.spin == 1:
                if doss.n_proj > 1:
                    if linestl_doss is None:
                        axs[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[index],
                                    label=labels[index], linewidth=linewidth)
                    else:
                        axs[1].plot(doss.doss[:, projection], dx_dos, color=color_doss[index],
                                    label=labels[index], linestyle=linestl_doss[index], linewidth=linewidth)
                else:
                    axs[1].plot(doss.doss[:, projection], dx_dos,
                                color=color_doss, linewidth=linewidth)
            elif doss.spin == 2:
                if doss.n_proj > 1:
                    axs[1].plot(doss.doss[:, projection, 0], dx_alpha, color=color_doss[index],
                                label=labels[index], linestyle='-', linewidth=linewidth)
                    axs[1].plot(-doss.doss[:, projection, 1], dx_beta, color=color_doss[index],
                                label=labels[index], linestyle='--', linewidth=linewidth)
                else:
                    axs[1].plot(doss.doss[:, projection, 0], dx_alpha, color=color_doss,
                                linestyle='-', linewidth=linewidth)
                    axs[1].plot(-doss.doss[:, projection, 1], dx_beta, color=color_doss,
                                linestyle='--', linewidth=linewidth)

    if ymin_bd > ymin_dos:
        ymin = ymin_dos
    elif ymin_bd <= ymin_dos:
        ymin = ymin_bd

    if ymax_bd >= ymax_dos:
        ymax = ymax_bd
    elif ymax_bd < ymax_dos:
        ymax = ymax_dos

    xmax_bd = hsp[len(hsp)-1]
    xmin_dos = 0

    yhsp = np.linspace(ymin-5, ymax+5, 2)
    for j in hsp:
        xhsp = np.ones(2)*j
        axs[0].plot(xhsp, yhsp, color='black', linewidth=0.5)

    hsp_label = []
    if k_labels is not None:
        for n in k_labels:
            if n in greek:
                g = greek.get(n)
                hsp_label.append(g)
            else:
                hsp_label.append(n)

    axs[0].set_xticks(hsp)
    if k_labels is not None:
        axs[0].set_xticklabels(hsp_label)

    xfermi_bd = np.linspace(xmin_bd, xmax_bd, 2)
    xfermi_dos = np.linspace(xmin_dos, xmax_dos, 2)
    yfermi = np.zeros(2)

    axs[0].plot(xfermi_bd, yfermi, color=fermi, linewidth=1.5)
    axs[1].plot(xfermi_dos, yfermi, color=fermi, linewidth=1.5)

    axs[0].set_xlim(xmin_bd, xmax_bd)

    if dos_max_range is not None:
        xmax_dos = dos_max_range

    if (prj is None) and (doss.n_proj not in prj):
        xmax_dos = np.amax(doss.doss[:, 1:doss.n_proj-1, :])

    axs[1].set_xlim(xmin_dos, xmax_dos+5)

    if energy_range is not None:
        ymin = energy_range[0]
        ymax = energy_range[1]

    plt.ylim(ymin, ymax)

    fig.text(.06, 0.5, '$E-E_F$ (eV)', ha='center',
             va='center', rotation='vertical')

    plt.legend()

    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()


def plot_cry_contour(contour_obj, diff=False, save_to_file = False):

    import matplotlib.pyplot as plt
    import os
    import numpy as np
    import time
    
    df = contour_obj.df
    n_punti_x = contour_obj.npx
    
    for i in range(0,8):
        df[i] = df[i].astype(float) 

    flat_list = [item for sublist in df.values for item in sublist]

    cleaned_list = [x for x in flat_list if ~np.isnan(x)]
    
    l = [cleaned_list[x:x+n_punti_x] for x in range(0, len(cleaned_list),n_punti_x)]
    
    c = contour_obj.x_graph_param
    d = contour_obj.y_graph_param
    
    plt.rcParams["figure.figsize"] = [c,d]

    plt.xlabel(r'$\AA$',fontsize=18)
    plt.ylabel(r'$\AA$',fontsize=18)

    X,Y = np.meshgrid(contour_obj.x_points,contour_obj.y_points) 
    
    
    levels = contour_obj.levels
    colors = contour_obj.colors
    linestyles = contour_obj.linestyles
    fmt = contour_obj.fmt

    #Change here to have or not the isovalues on the plot
    iso = True
    #iso = False

    if (iso == True):
        L = plt.contour(X, Y, l, levels = levels, colors = colors, linestyles = linestyles, linewidths = 0.7,
                    alpha = 1)
        plt.clabel(L, inline = 1, fontsize = 7, fmt = fmt)
    elif (iso == False):
        L = plt.contour(X, Y, l, levels = levels, colors = colors, linestyles = linestyles, linewidths = 0.7,
                    alpha = 1)

    
    path = os.path.join('./'+'figure_' + contour_obj.tipo + '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.savefig(path, bbox_inches = 'tight',dpi=600)
    print('\nThe image has been saved in the current directory')

    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()


def plot_cry_contour_differences(contour_obj, contour_obj_ref, save_to_file = False):
    
    import matplotlib.pyplot as plt
    import os
    import numpy as np
    import time
    
    n_punti_x = contour_obj.npx
    
    df = contour_obj.df  
    for i in range(0,8):
        df[i] = df[i].astype(float) 
        
    df_ref = contour_obj_ref.df
    for i in range(0,8):
        df_ref[i] = df_ref[i].astype(float)   
        
    df_diff = df - df_ref
        
    flat_list = [item for sublist in df_diff.values for item in sublist]

    cleaned_list = [x for x in flat_list if ~np.isnan(x)]
    
    l = [cleaned_list[x:x+n_punti_x] for x in range(0, len(cleaned_list),n_punti_x)]
    
    c = contour_obj.x_graph_param
    d = contour_obj.y_graph_param
    
    plt.rcParams["figure.figsize"] = [c,d]

    plt.xlabel(r'$\AA$',fontsize=18)
    plt.ylabel(r'$\AA$',fontsize=18)

    X,Y = np.meshgrid(contour_obj.x_points,contour_obj.y_points) 
    
    
    ctr1dif = np.array([-8,-4,-2,-0.8,-0.4,-0.2,-0.08,-0.04,-0.02,-0.008,-0.004,-0.002,-0.0008,-0.0004,-0.0002,0,
                       0.0002,0.0004,0.0008,0.002,0.004,0.008,0.02,0.04,0.08,0.2,0.4,0.8,2,4,8])
    colors1dif = ['b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','k','r','r','r','r','r','r','r',
                  'r','r','r','r','r','r','r','r']
    ls1dif = ['--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','dotted','-','-','-','-',
              '-','-','-','-','-','-','-','-','-','-','-']

    levels = ctr1dif
    colors = colors1dif
    linestyles = ls1dif
    fmt = '%1.4f'

    #Change here to have or not the isovalues on the plot
    iso = True
    #iso = False

    if (iso == True):
        L = plt.contour(X, Y, l, levels = levels, colors = colors, linestyles = linestyles, linewidths = 0.7,
                    alpha = 1)
        plt.clabel(L, inline = 1, fontsize = 7, fmt = fmt)
    elif (iso == False):
        L = plt.contour(X, Y, l, levels = levels, colors = colors, linestyles = linestyles, linewidths = 0.7,
                    alpha = 1)
    
    path = os.path.join('./'+'figure_diff_' + contour_obj.tipo + '_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.savefig(path, bbox_inches = 'tight',dpi=600)
    print('\nThe image has been saved in the current directory')

    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()


def plot_cry_xrd(xrd_obj, save_to_file = False):

    import matplotlib.pyplot as plt
    import os
    import time
    
    plt.rcParams["figure.figsize"] = [16,9]

    plt.plot(xrd_obj.x,xrd_obj.y)

    plt.xlim((0, 30))
        
    path = os.path.join('./'+'figure_'+'XRD_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(xrd_obj.title, fontsize=20)
    plt.savefig(path, bbox_inches = 'tight',dpi=600)
    
    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()

   
def plot_cry_rholine(rholine_obj, save_to_file = False):
    
    import matplotlib.pyplot as plt
    import os
    import time

    plt.plot(rholine_obj.x,rholine_obj.y)

    plt.xlabel('d  [$\AA$]',fontsize=14)
    plt.ylabel(r'$\rho$  [$\frac{e}{\AA^3}$]',fontsize=16)
        
    path = os.path.join('./'+'figure_'+'rholine_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(rholine_obj.title, fontsize=15)
    plt.savefig(path, bbox_inches = 'tight',dpi=600)

    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()


def plot_cry_seebeck(seebeck_obj, save_to_file = False):

    import sys
    import matplotlib.pyplot as plt
    import numpy as np
    import time
    
    case = input('Please, choose the direction you want to plot. \nYou can choose among S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_zx, S_zy, S_zz\n')

    case = case.lower().replace('_','')

    if case.isalpha() == True:
        pass
    else:
        sys.exit('Please, select a valid chioce')
        
    if case == 'sxx':
        col = 3
    elif case == 'sxy':
        col = 4
    elif case == 'sxz':
        col = 5
    elif case == 'syx':
        col = 6
    elif case == 'syy':
        col = 7
    elif case == 'syz':
        col = 8
    elif case == 'szx':
        col = 9
    elif case == 'szy':
        col = 10
    elif case == 'szz':
        col = 11
    else:
        sys.exit('please, choose a valid chioce')
        
    x = [] # qui metto i potenziali alle diverse T (che saranno sempre uguali)
    for k in range(0,len(seebeck_obj.all_data)):
            x.append(np.array(seebeck_obj.all_data[k].apply(lambda x: float(x.split()[0]))))
            
    y = [] # qui metto i valori di y che vuole plottare l'utente
    for k in range(0,len(seebeck_obj.all_data)):
            y.append(np.array(seebeck_obj.all_data[k].apply(lambda x: float(x.split()[col])*1000000)))
            
    for k in range(0,len(seebeck_obj.all_data)):
        plt.figure()
        plt.plot(x[k],y[k],label=str(seebeck_obj.temp[k])+' K')
        plt.xlabel('Chemical Potential (eV)',fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)',fontsize=12)
        plt.axhline(0, color='k')
        plt.title(seebeck_obj.title)
        plt.legend(loc='upper left',fontsize=12)
        plt.savefig('seebeck_at_' + str(seebeck_obj.temp[k]) + 'T___' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg',format='jpg',dpi=600,bbox_inches='tight')
        plt.show()
        
        
    for k in range(0,len(seebeck_obj.all_data)):
        plt.plot(x[k],y[k],label=str(seebeck_obj.temp[k])+' K')
        plt.xlabel('Chemical Potential (eV)',fontsize=12)
        plt.ylabel('Seebeck Coefficient ($\mu$V/K)',fontsize=12)
        plt.title(seebeck_obj.title)
        plt.axhline(0, color='k')
        plt.legend(loc='upper left',fontsize=12)
    
    if save_to_file != False:
        save_plot(save_to_file)


def plot_cry_lapl_profile(lapl_obj, save_to_file = False):

    import matplotlib.pyplot as plt
    import time

    plt.plot(lapl_obj.datax,lapl_obj.datay)

    plt.fill_between(lapl_obj.datax,lapl_obj.datay,where=(lapl_obj.datay < 0),color='lightblue', interpolate=True)
    plt.fill_between(lapl_obj.datax,lapl_obj.datay,where=(lapl_obj.datay > 0),color='lightcoral', interpolate=True)

    #plt.xlim(-0.5,0.5)
    #plt.ylim(-200,200)

    plt.xlabel('Distance [A]') 
    plt.ylabel('Laplacian [e/A^5]')
        
    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()

def plot_cry_density_profile(lapl_obj, save_to_file = False):
    
    import matplotlib.pyplot as plt
    import time


    plt.plot(lapl_obj.datax,lapl_obj.datay) 

    plt.xlabel('Distance [A]') 
    plt.ylabel('Density [e/A^3]')
        
    plt.savefig('Density_profile' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg',format='jpg',dpi=600,bbox_inches='tight')

    if save_to_file != False:
        save_plot(save_to_file)

    plt.show()

def save_plot(path_to_file):

    from os import path
    import sys
    import matplotlib as plt
       
    folder = path.split(path_to_file)[0]
    if path.exists(folder) == True:
        plt.savefig('%s.png'%path_to_file)
    else:
        print('Error: folder %s does not exist'%path_to_file)
        sys.exit(1)
    