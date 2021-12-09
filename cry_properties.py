#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:35:48 2021

@author: brunocamino
"""

def newk(shrink1,shrink2,Fermi=1,print_option=0,title=None):
    
    newk_block = ['NEWK\n','%s %s\n' %(shrink1,shrink2),
                  '%s %s\n'%(Fermi,print_option)]
    
    return newk_block
    
def bands(k_path,n_kpoints,first_band,last_band,print_eig=0,print_option=1,
          title='Bands'):
    import numpy as np
    
    bands_block = []
    
    k_path = np.array(k_path)
    k_unique = np.unique(k_path)
    
    #Find the shrinking factor
    k_unique = np.array(np.around(k_unique,4)*10000,dtype=int)
    if len(k_unique) > 2:
        gcd = np.gcd.reduce(k_unique) 
    else:
        gcd = np.gcd(k_unique[0],k_unique[1])
    k_path = np.array((k_path/gcd)*10000,dtype=int)
    shrink = int(10000/gcd)
    
    bands_block.append('BAND\n')
    if title == None:
        bands_block.append('BAND STRUCTURE CALCULATION\n')
    else:
        bands_block.append(title+'\n')
    
    bands_block.append(str(len(k_path)-1)+' '+str(shrink)+' '+str(n_kpoints)+
                       ' '+str(first_band)+' '+str(last_band)+' '+
                       str(print_eig)+' '+str(print_option)+'\n')
    
    
    #Add the symmetry lines
    for i in range(len(k_path[:-1])):
        bands_block.append(' '.join([str(x) for x in k_path[i]])+'  '+
                           ' '.join([str(x) for x in k_path[i+1]])+'\n')
    
    bands_block.append('END\n')
    
    print(bands_block)
    
    

###TESTING
k_path = [[0,0,0],[0.5,0,0],[0.5,0.5,0.5],[0.25,0,0.5]]

bands(k_path,200,1,26)
    