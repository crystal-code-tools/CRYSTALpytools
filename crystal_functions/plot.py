#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/03/2022
"""


from pyparsing import counted_array


def plot_cry_bands(bands, k_labels=None, energy_range=None, title=False, not_scaled=False, mode='single', linestl='-', linewidth=1, color='blue', fermi='forestgreen', k_range=None, labels=None):

    import matplotlib.pyplot as plt
    # import matplotlib.lines as mlines
    import numpy as np
    import sys

    greek = {'Alpha': '\u0391', 'Beta': '\u0392', 'Gamma': '\u0393', 'Delta': '\u0394', 'Epsilon': '\u0395', 'Zeta': '\u0396', 'Eta': '\u0397', 'Theta': '\u0398', 'Iota': '\u0399',
             'Kappa': '\u039A', 'Lambda': '\u039B', 'Mu': '\u039C', 'Nu': '\u039D', 'Csi': '\u039E', 'Omicron': '\u039F', 'Pi': '\u03A0', 'Rho': '\u03A1', 'Sigma': '\u03A3', 'Tau': '\u03A4',
             'Upsilon': '\u03A5', 'Phi': '\u03A6', 'Chi': '\u03A7', 'Psi': '\u03A8', 'Omega': '\u03A9', 'Sigma_1': '\u03A3\u2081'}

    # Error check on the mode flag
    modes = ['single', 'multi', 'compare', 'surface']
    if mode not in modes:
        print('The selected mode '+mode+' is not among the possible ones: ' +
              modes[0]+', ' + modes[1] + ', '+modes[2] + ', or '+modes[3])
        sys.exit()

    # Error chenk on k_label
    if k_labels != None:

        if type(k_labels) != list:
            print('k_labels must be a list of strings')
            sys.exit()

        elif type(k_labels) == list:
            for element in k_labels:
                if type(element) != str:
                    print('k_label must be a list of strings:' +
                          str(element)+' is not a string')
                    sys.exit()

    # Error check on energy range
    if energy_range != None:

        if type(energy_range) != list:
            print('energy_range must be a list of two int or float (min,max)')
            sys.exit()

        elif type(energy_range) == list:

            if len(energy_range) > 2:
                print('energy_range must be a list of two int or float (min,max)')
                sys.exit()

            for element in energy_range:
                if (type(element) != int) or (type(element) != float):
                    print('energy_range must be a list of two int or float (min,max): ' +
                          str(element)+', is neither a float or an int')
                    sys.exit()

    # Error check on k_range
    if k_range != None:

        if type(k_range) != list:
            print('k_range must be a list of two strings')
            sys.exit()

        elif type(k_range) == list:
            if len(k_range) > 2:
                print('k_range must be a list of two strings')
                sys.exit()

            for element in k_range:
                if type(element) != str:
                    print('k_label must be a list of two strings:' +
                          str(element)+' is not a string')
                    sys.exit()

    # Error check on title
    if title != False:
        if type(title) != str:
            print('title needs to be a string')
            sys.exit

    # plotting of a single band object
    if mode == modes[0]:

        dx = bands.k_point_plot

        pltband = bands.bands
        no_bands = np.shape(pltband)[0]
        ymin = np.amin(pltband)
        ymax = np.amax(pltband)
        xmin = min(dx)
        xmax = max(dx)

        # band plot
        for i in range(no_bands):

            if bands.spin == 1:
                plt.plot(dx, pltband[i, :], color=color,
                         linestyle=linestl, linewidth=linewidth)

            elif bands.spin == 2:
                plt.plot(dx[:, 0], pltband[i, :, 0], color='red',
                         linestyle=linestl, linewidth=linewidth, label='Alpha')
                plt.plot(dx[:, 1], pltband[i, :, 1], color='black',
                         linestyle=linestl, linewidth=linewidth, label='Beta')

    # plot of multiple band objects on a single plot
    elif mode == modes[1]:

        # Error check on the band on the 'multi' mode flag
        if type(bands) != list:
            print('When you choose a ' +
                  modes[1]+' plot bands needs to be a list of band objects')
            sys.exit()

        # Error check on color for the 'multi' mode flag
        if type(color) != list:
            print('When you choose a ' +
                  modes[1]+' plot color needs to be a list')
            sys.exit()

        elif type(color) == list:
            if len(color) > len(bands):
                print(
                    'The number of colors is greater than the number of objects you want to plot')
                sys.exit()

        # Warning comparison with band.spin==2
        for m in bands:
            if m.spin == 2:
                print(
                    "Warning: the 'multi' plot is not available at the moment for file with NSPIN = 2")

        # scaling that enables the comparison of band structure calculated at different pressures
        if not_scaled == False:
            reference = xmax = np.amax(bands[0].k_point_plot)
            xmin = np.amin(bands[0].k_point_plot)

        else:
            xmax = []
            xmin = []

        ymin = []
        ymax = []

        # plot of all the bands obj present in the list
        for index, data in enumerate(bands)
        # scaling that enables the comparison of band structure calculated at different pressures
           if not_scaled == False:
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

            for j in range(no_bands):
                
                if count1 == count2:
                    if type(linestl)==list:
                        plt.plot(dx, pltband[j, :], color=color[index],linestyle=linestl[index], linewidth=linewidth, label=labels[index])
                    else:
                        plt.plot(dx, pltband[j, :], color=color[index],linestyle=linestl, linewidth=linewidth, label=labels[index])
                        
                else:
                    if type(linestl)==list:
                        plt.plot(dx, pltband[j, :], color=color[index],linestyle=linestl[index], linewidth=linewidth)
                    else:
                        plt.plot(dx, pltband[j, :], color=color[index],linestyle=linestl, linewidth=linewidth)
                count1+=1
                
    # HSP line plot
    if type(bands) == list:
        hsp = bands[0].tick_position
    else:
        hsp = bands.tick_position
    y_band = np.linspace(ymin-3, ymax+3, 2)
    high_sym_point = []
    hsp_label = []
    for j in hsp:
        x_band = np.ones(2)*j
        plt.plot(x_band, y_band, color='black', linewidth=0.5)

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

    # plot of the fermi level
    x = np.linspace(xmin, xmax, 2)
    y = np.zeros(2)
    plt.plot(x, y, color=fermi, linewidth=2.5)

    # definition of the ylim
    if energy_range != None:
        ymin = energy_range[0]
        ymax = energy_range[1]

    # definition of the xlim
    if k_range != None:
        for i in range(0, len(k_range)):
            j = k_range[i]
            if j in greek:
                g = greek.get(j)
                k_range[i] = g
        xmin = path_dict[k_range[0]]
        xmax = path_dict[k_range[1]]

    # definition of the plot title
    if title != False:
        plt.title(title)

    if bands.spin == 2:
        plt.legend()

    plt.xticks(hsp, hsp_label)
    plt.ylabel('$E-E_F$ (eV)')
    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax) 

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
