#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:35:48 2021

@author: brunocamino
"""

def cry_newk(shrink1,shrink2,Fermi=1,print_option=0,title=None):
    newk_block = ['NEWK\n','%s %s\n' %(shrink1,shrink2),
                  '%s %s\n'%(Fermi,print_option)]
    return newk_block
    
def cry_bands(k_path,n_kpoints,first_band,last_band,print_eig=0,print_option=1,
          title='BAND STRUCTURE CALCULATION'):
    #k_path can be:
        ##list of list
        ##pymatgen HighSymmKpath object
    
    import numpy as np
    import sys
    
    bands_block = []
    
    if 'HighSymmKpath' in str(type(k_path)):
        k_path_flat = [item for sublist in k_path.kpath['path'] for item in sublist]
        k_path_pmg = []
        for i in k_path_flat:
            #This is a pmg HighSymmKpath object
            k_path_pmg.append(k_path.kpath['kpoints'][i].tolist())  
        k_path = np.array(k_path_pmg)
    
    elif type(k_path[0]) == list:
        #This is a list of lists
        k_path = np.array(k_path)           
    
    else:    
        print('EXITING: k_path type must be a list of list (k coordinates) or\
              a pymatgen HighSymmKpath object. %s selected' %type(k_path) )
        sys.exit(1)
    
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
    bands_block.append(title+'\n')
    
    bands_block.append(str(len(k_path)-1)+' '+str(shrink)+' '+str(n_kpoints)+
                       ' '+str(first_band)+' '+str(last_band)+' '+
                       str(print_option)+' '+str(print_eig)+'\n')
    
    
    #Add the symmetry lines
    for i in range(len(k_path[:-1])):
        bands_block.append(' '.join([str(x) for x in k_path[i]])+'  '+
                           ' '.join([str(x) for x in k_path[i+1]])+'\n')
    
    bands_block.append('END\n')
    
    return bands_block

'''###TESTING
k_path = [[0,0,0],[0.5,0,0],[0.5,0.5,0.5],[0.25,0,0.5]]
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer'''



def cry_doss(n_points=200,band_range=None,e_range=None,plotting_option=2,
             poly=12,print_option=1):
     #e_range in eV
    import sys
    
    doss_block = []
    if band_range == None and e_range == None:
        print('EXITING: please specify either band_range or e_range. None selected')
        sys.exit(1)
    elif band_range != None and e_range != None:
        print('EXITING: please specify either band_range or e_range. Both selected')
        sys.exit(1)
    elif type(band_range) == list and len(band_range) == 2:
        doss_range = band_range
    elif type(e_range) == list and len(e_range) == 2:
        doss_range = [-1,-1]
        
    else:
        print('EXITING: either the band_range argument or the e_range argument\
              do not match the required format (2 item list)')
        sys.exit(1)
    
    doss_block.append('DOSS\n')
    doss_block.append(str(0)+' '+str(n_points)+' '+str(doss_range[0])+' '+
                      str(doss_range[1])+' '+str(plotting_option)+' '+str(poly)+' '+
                      str(print_option)+'\n')
    
    if doss_range == [-1,-1]:
        doss_block.append(str(e_range[0]/27.2114)+' '+str(e_range[1]/27.2114)+'\n')
    
    doss_block.append('END\n')
    
    return doss_block
    
'''###TESTING
print(cry_doss(e_range=[-5,10]))'''

def cry_pdoss(projections,proj_type='atom',output_file=None,n_points=200,band_range=None,
              e_range=None,plotting_option=2,poly=12,print_option=1):

    #projections is a list of lists
    import sys
    from crystal_functions.file_readwrite import Crystal_output
    pdoss_block = []
    if band_range == None and e_range == None:
        print('EXITING: please specify either band_range or e_range. None selected')
        sys.exit(1)
    elif band_range != None and e_range != None:
        print('EXITING: please specify either band_range or e_range. Both selected')
        sys.exit(1)
    elif type(band_range) == list and len(band_range) == 2:
        pdoss_range = band_range
    elif type(e_range) == list and len(e_range) == 2:
        pdoss_range = [-1,-1]
        
    else:
        print('EXITING: either the band_range argument or the e_range argument\
              do not match the required format (2 item list)')
        sys.exit(1)
    
    pdoss_block.append('DOSS\n')
    pdoss_block.append(str(len(projections))+' '+str(n_points)+' '+str(pdoss_range[0])+' '+
                      str(pdoss_range[1])+' '+str(plotting_option)+' '+str(poly)+' '+
                      str(print_option)+'\n')
    
    if pdoss_range == [-1,-1]:
        pdoss_block.append(str(e_range[0]/27.2114)+' '+str(e_range[1]/27.2114)+'\n')
    
    flat_proj = [x for sublist in projections for x in sublist]
    if all(isinstance(x, int) for x in flat_proj):        
        if proj_type == 'atom':
            for proj in projections:                
                pdoss_block.append(str(-len(proj))+' '+' '.join([str(x) for x in proj])+'\n')
        if proj_type == 'ao':
            for proj in projections:  
                pdoss_block.append(str(len(proj))+' '+' '.join([str(x) for x in proj])+'\n') 
        elif proj_type != 'atom' and proj_type != 'ao':
            print('EXITING: please specify either atom or ao projection. %s selected' % proj_type)
            sys.exit(1)
    elif all(isinstance(x, str) for x in flat_proj):
       if output_file == None:
           print('EXITING: please specify an outut file to use the atoms projection.')
           sys.exit(1)
       else:
           output = Crystal_output(output_file)
           output.extract_last_geom()
           atoms_symbols = output.atom_symbols
           atoms_symbols.insert(0, 0)
           
           for proj in projections:
               atom_positions_list = []
               for element in proj:
                   index = [i for i,ele in enumerate(atoms_symbols) if ele==element.upper()]
                   atom_positions_list.append([str(x) for x in index])
               pdoss_block.append(str(-len(index))+' '+' '.join([str(x) for x in index])+'\n')        
   
    pdoss_block.append('END\n')
    
    return pdoss_block

'''###TESTING
projections = [[1,2],[3,4]]
projections = [['Mg'],['O']]
output_file = 'examples/data/mgo.out'
print(cry_pdoss(projections,proj_type='atom',output_file=output_file,e_range=[-5,10]))'''
    