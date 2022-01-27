#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:29:16 2021

@author: brunocamino
"""

def cry_out2pmg(output,initial=False,dimensionality=3,vacuum=10):
    #output_file is a crystal output object
    
    from pymatgen.core import Structure
    
    import numpy as np
   
    output.extract_last_geom(write_gui_file=False,print_cart=False)
    output.primitive_lattice(initial=initial)
    
    if dimensionality == 3:
        structure = Structure(output.primitive_vectors, output.atom_numbers, 
                              output.atom_positions_cart, coords_are_cartesian=True)
        
    elif dimensionality == 2:
        thickness = np.amax(np.array(output.atom_positions_cart)[:,2])-\
                    np.amin(np.array(output.atom_positions_cart)[:,2])
        
        vectors = output.primitive_vectors
        vectors[2,2] = thickness + vacuum
        
        structure = Structure(vectors, output.atom_numbers, 
                              output.atom_positions_cart, coords_are_cartesian=True)
    
    return structure
    
    
'''###TESTING
cry_out2pmg('examples/data/mgo.out')'''

def cry_bands2pmg(output,bands,labels=None):
    #WORK IN PROGRESS
    #Function to transform a crystal band file into a pmg band object
    #Format of the pmg object:
    #classBandStructure(kpoints, eigenvals, lattice, efermi, labels_dict=None, coords_are_cartesian=False, structure=None, projections=None)
    
    import numpy as np
    from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine #,Kpoint
    
    from pymatgen.core.lattice import Lattice
    from pymatgen.electronic_structure.core import Spin
    

    output.reciprocal_lattice()
    labels_dict = {}
    
    if labels != None:        
        for i,j in enumerate(bands.n_points):
            labels_dict[labels[i]] = bands.k_point_coordinates[j-1]
    
    '''#This defines the Kpoint objects. Not needed at the moment, but might be useful in the future
    k_points = []
    for i, coord in enumerate(bands.k_point_coordinates):
        k_points.append(Kpoint(np.array(coord), 
                               Lattice(output.reciprocal_vectors)))
        if len(bands.n_points) > 1:
            if i+1 in bands.n_points[1:-1]:
                print(i+1)
                k_points.append(Kpoint(np.array(coord), 
                               Lattice(output.reciprocal_vectors)))'''
    
    #List of k points coordinates as symmetry lines
    band_energy = bands.bands 
    
    #pymatgen will plot the bands wrt to the Fermi Energy
    band_energy[1:,:,:] = band_energy[1:,:,:] + output.fermi_energy()
    k_points_coordinates = []
    
    for i, coord in enumerate(bands.k_point_coordinates):
        k_points_coordinates.append(bands.k_point_coordinates)
        if len(bands.n_points) > 1:
            if i+1 in bands.n_points[1:-1]:
                k_points_coordinates.append(bands.k_point_coordinates)
                
    k_points_coordinates = bands.k_point_coordinates
    if len(bands.n_points) > 1:
        for i,point in enumerate(bands.n_points[1:-1]):
            k_points_coordinates.insert(point-1+i,k_points_coordinates[point-1+i])
            band_energy = np.insert(band_energy,point-1+i,band_energy[:,point-1+i,:],axis=1)
    eigenvals = {}
    eigenvals[Spin.up] = band_energy[:,:,0]
    
    if len(bands.bands[0,0,:]) > 1:
        eigenvals[Spin.down] = band_energy[:,:,1]
    
    return BandStructureSymmLine(k_points_coordinates, eigenvals, 
                                 Lattice(output.reciprocal_vectors), 
                                 bands.efermi,labels_dict,
                                 coords_are_cartesian=False)
    
'''###TESTING
from cry_file_readwrite import Crystal_bands, Crystal_output
from pymatgen.core.lattice import Lattice

mgo_bands = Crystal_bands('examples/data/mgo_SPIN_BAND_dat.BAND')
mgo_out = Crystal_output('examples/data/mgo_SPIN.out')
#print(mgo_bands.n_points,labels=['A','B','C','D','E'])
#mgo_out.reciprocal_lattice()


#cry_bands2pmg(bands,output_name)
#mgo_output_name = 'data/mgo_BAND.outp'
cry_bands2pmg(mgo_out,mgo_bands,labels=['A','B','C','D','E'])'''

def cry_gui2pmg(gui_file):
    
    pass



