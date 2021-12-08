#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:35:48 2021

@author: brunocamino
"""

def newk(shrink1,shrink2,Fermi=1,print_option=0):
    
    newk_block = ['NEWK\n','%s %s\n' %(shrink1,shrink2),
                  '%s %s\n'%(Fermi,print_option)]
    
    return newk_block
    
def bands(k_path,n_kpoints,first_band,last_band,print_eig=0,print_option=1):
    pass
    