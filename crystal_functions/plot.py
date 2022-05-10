#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import sys
import ase
from ase import Atoms 
from ase.visualize import view 
from ase.visualize.plot import plot_atoms


# In[2]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/03/2022
"""

def plot_cry_bands(bands,k_labels=None,energy_range=None,title=False,not_scaled=False):
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import numpy as np
    
    kpoints = bands.tick_position 
    efermi_band = bands.efermi    
                
    x1 = bands.bands[:,0,:]
    print(x1)
    y1 = bands.bands[:,1:,:]
    
    fig, ax1 = plt.subplots(1,1)
    
    if len(bands.bands[0,0,:]) > 1:  
        
        ax1.plot(x1[:,1], y1[:,:,1],'-',  c='red' )
        ax1.plot(x1[:,0], y1[:,:,0],'-' , c='black')        
        #ax1.legend( bbox_to_anchor=(0.9, .75))
        spin_band = [mlines.Line2D([], [], color='black', label='Alpha'),
                     mlines.Line2D([], [], color='red', label='Beta')]
        ax1.legend(spin_band,['Alpha electrons', 'Beta electrons'],facecolor='white', framealpha=1,bbox_to_anchor=(.83,.90))
    
    else:
        ax1.plot(x1[:,0], y1[:,:,0],'-', c='black' )  
    
    #Display E Fermi on band structure
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
        plt.subplots_adjust(wspace=0.2, top=0.88)

def plot_cry_multibands(bands_list,k_labels=None,energy_range=None,title=False,not_scaled=False):
    #Filippo's function
    #bands list is a Crystal_band object
    pass

def plot_cry_bands_compare():
    pass

def plot_cry_doss():
    pass

def plot_cry_es():
    pass


#coss
def plot_contour(contour_obj):
    
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

    plt.show()
    
#coss
def plot_contour_differences(contour_obj, contour_obj_ref):
    
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

    plt.show()
    

#coss
def plot_XRD(xrd_obj):
    
    plt.rcParams["figure.figsize"] = [16,9]

    plt.plot(xrd_obj.x,xrd_obj.y)

    plt.xlim((0, 30))
        
    path = os.path.join('./'+'figure_'+'XRD_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(xrd_obj.title, fontsize=20)
    plt.savefig(path, bbox_inches = 'tight',dpi=600)
    
    plt.show()

    
#coss    
def plot_rholine(rholine_obj):
    
    plt.plot(rholine_obj.x,rholine_obj.y)

    plt.xlabel('d  [$\AA$]',fontsize=14)
    plt.ylabel(r'$\rho$  [$\frac{e}{\AA^3}$]',fontsize=16)
        
    path = os.path.join('./'+'figure_'+'rholine_' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg')
    plt.title(rholine_obj.title, fontsize=15)
    plt.savefig(path, bbox_inches = 'tight',dpi=600)

    plt.show()
    
    
#coss
def plot_out_molecule(mol_obj):
    return view(mol_obj.atoms,viewer='ngl')

#coss
def plot_out_opt_molecule(mol_opt_obj):
    return view(mol_opt_obj.atoms,viewer='ngl')

#coss
def plot_out_crystal(cry_obj):
    return view(cry_obj.cell,viewer='ngl')

#coss
def plot_out_opt_crystal(cry_opt_obj):
    return view(cry_opt_obj.cell,viewer='ngl')

#coss
def plot_seebeck(seebeck_obj):
    
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
    plt.savefig('seebeck_different_T_at_' + str(seebeck_obj.temp[k]) + 'T___' + time.strftime("%Y-%m-%d_%H%M%S") + '.jpg',format='jpg',dpi=600,bbox_inches='tight')
        #Do NOT put here plt.show()
        


# In[ ]:




