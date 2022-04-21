#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29/03/2022

"""

def help(folder='./'):
    
    import os.path

    if os.path.isfile(os.path.join(folder,'help.txt')):
        print ("Display the txt file")
    elif os.path.isfile(os.path.join(folder,'help.md')) or os.path.isfile(os.path.join(folder,'help.html')):
        print ("Display the md text")
    elif os.path.isfile(os.path.join(folder,'help.jpg')) or os.path.isfile(os.path.join(folder,'help.png')):
        print ("Display the image")


def view_pmg(pmg_structure):
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.visualize import view
    
    return view(AseAtomsAdaptor().get_atoms(pmg_structure))
