#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:28:28 2021
"""

class Crystal_input:
    #This creates a crystal_input object
    
    def __init__(self,input_name):
        import sys
        
        self.name = input_name
        try: 
            if input_name[-3:] != 'd12':
                input_name = input_name+'.d12'
            file = open(input_name, 'r')
            data = file.readlines()
            file.close()
        except:
            print('EXITING: a .d12 file needs to be specified')
            sys.exit(1)
        
        if 'BASISSET\n' in data:
            end_index = []
            if 'OPTGEOM\n' in data:
                end_index = [i for i, s in enumerate(data) if 'END' in s]
                end_index.insert(1,data.index('BASISSET\n')+1)
            else:
                end_index.append(data.index('BASISSET\n')-1)
                end_index.append(data.index('BASISSET\n')+1)
                end_index.extend([i for i, s in enumerate(data[end_index[1]:]) if 'END' in s])
            self.is_basisset = True
        else:        
            end_index = [i for i, s in enumerate(data) if 'END' in s]
            self.is_basisset = False
        self.geom_block = []
        #self.optgeom_block = []
        self.bs_block = []
        self.func_block = []
        self.scf_block = []
        

        if len(end_index) == 4:
            self.geom_block = data[:end_index[0]+1]
            #self.optgeom_block = []
            self.bs_block = data[end_index[0]+1:end_index[1]+1]
            self.func_block = data[end_index[1]+1:end_index[2]+1]
            #The following loop ensures that keyword+values (such as TOLINTEG 7 7 7 7 14
            #are written as a list)
            for i in range(end_index[2]+1,end_index[-1]):
                if data[i+1][0].isnumeric():
                    self.scf_block.append([data[i],data[i+1]])
                else:
                    if data[i][0].isalpha():
                        self.scf_block.append(data[i])
                    else:
                        pass
            self.scf_block.append('ENDSCF\n') #The loop cannot go over the last element
            #This is the old one, remove if not needed: self.scf_block = data[end_index[2]+1:]
        elif len(end_index) == 5:
            self.geom_block = data[:end_index[1]+1]
            #self.optgeom_block = data[end_index[0]+1:end_index[1]+1]
            self.bs_block = data[end_index[1]+1:end_index[2]+1]
            self.func_block = data[end_index[2]+1:end_index[3]+1]
            #The following loop ensures that keyword+values (such as TOLINTEG 7 7 7 7 14
            #are written as a list)
            for i in range(end_index[3]+1,end_index[-1]):
                if data[i+1][0].isnumeric():
                    self.scf_block.append([data[i],data[i+1]])
                else:
                    if data[i][0].isalpha():
                        self.scf_block.append(data[i])
                    else:
                        pass
            self.scf_block.append('ENDSCF\n') #The loop cannot go over the last element
            
    def add_ghost(self,ghost_atoms):
        if self.is_basisset == True:
            self.bs_block.append('GHOSTS\n')
            self.bs_block.append('%s\n' % len(ghost_atoms))  
            self.bs_block.append(' '.join([str(x) for x in ghost_atoms])+'\n' )
        else:
            self.bs_block.insert(-1, 'GHOSTS\n')
            self.bs_block.insert(-1,'%s\n' % len(ghost_atoms))  
            self.bs_block.insert(-1,' '.join([str(x) for x in ghost_atoms])+'\n' )   

    def opt_to_sp(self):   
        
        if 'OPTGEOM\n' in self.geom_block:
            init = self.geom_block.index('OPTGEOM\n')
            final = self.geom_block.index('END\n')
            del self.geom_block[init:final+1]
    
    def sp_to_opt(self):           
        if 'OPTGEOM\n' not in self.geom_block:
            if self.is_basisset == True:
                self.geom_block.append('OPTGEOM\n')
                self.geom_block.append('END\n')
            else:
                self.geom_block.insert(-1, 'OPTGEOM\n')
                self.geom_block.insert(-1, 'END\n')                

       
'''###TESTING
mgo = Crystal_input('examples/data/mgo.d12')
mgo.add_ghost([1,2,3])     
#print(mgo.name)
#print(mgo.scf_block)
#mgo.scf_block.remove('DIIS\n')
#print(mgo.scf_block)'''

class Crystal_output:
#This class reads a CRYSTAL output and generates an object

    def __init__(self,output_name):

        import sys
        import re

        self.name = output_name
        #Check if the file exists
        try: 
            if output_name[-3:] != 'out' and  output_name[-4:] != 'outp':
                output_name = output_name+'.out'
            file = open(output_name, 'r')
            self.data = file.readlines()
            file.close()
        except:
            print('EXITING: a .out file needs to be specified')
            sys.exit(1)

        #Check the calculation converged
        self.converged = False

        for i,line in enumerate(self.data[::-1]):
            if re.match(r'^ EEEEEEEEEE TERMINATION',line):
                self.converged = True
                #This is the end of output
                self.eoo = len(self.data)-1-i
                break

        
        if self.converged == False:
            self.eoo = len(self.data)
            #print('WARNING: the calculation did not converge. Proceed with care!')

    def get_dimensionality(self):

        import re
        
        for line in self.data:
            if re.match(r'\sGEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM',line) != None:
                self.dimensionality = int(line.split()[9])    
                return self.dimensionality
        
    def get_final_energy(self): 
        
        import re
        
        self.final_energy = None
        for line in self.data[self.eoo::-1]:
            if re.match(r'\s\W OPT END - CONVERGED',line) != None:
                self.final_energy = float(line.split()[7])*27.2114
            elif re.match(r'^ == SCF ENDED',line) != None:
                self.final_energy =  float(line.split()[8])*27.2114 
        
        if self.final_energy == None:
            print('WARNING: no final energy found in the output file. energy = None')
        
        return self.final_energy
    
    def get_scf_convergence(self,all_cycles=False):
        
        import re
        import numpy as np
        
        self.scf_energy = []
        self.scf_deltae = []
        
        scf_energy = []
        scf_deltae = []                    
        
        for line in self.data:
            
            
            if re.match(r'^ CYC ',line):
                scf_energy.append(float(line.split()[3]))
                scf_deltae.append(float(line.split()[5]))
            
            if re.match(r'^ == SCF ENDED - CONVERGENCE ON ENERGY',line):
                if all_cycles == False:
                    self.scf_energy = np.array(scf_energy)*27.2114
                    self.scf_deltae = np.array(scf_deltae)*27.2114
                    
                    return self.scf_energy, self.scf_deltae
                
                elif all_cycles == True:
                    self.scf_energy.append(scf_energy)
                    self.scf_deltae.append(scf_deltae)
                    scf_energy = []
                    scf_deltae = []    

            self.scf_convergence = [self.scf_energy, self.scf_deltae]
        return self.scf_convergence


    def get_num_cycles(self):

        import re

        for line in self.data[::-1]:
            if re.match(r'^ CYC ', line):
                self.num_cycles = int(line.split()[1])
                return self.num_cycles
        return None

    def get_fermi_energy(self):

        import re
        
        self.fermi_energy = None
        
        for i,line in enumerate(self.data[len(self.data)::-1]):     
            #This is in case the .out is from a DOSS calculation
            if re.match(r'^ TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT BAND',self.data[len(self.data)-(i+4)]) != None:
                for j,line1 in enumerate(self.data[len(self.data)-i::-1]):
                    if re.match(r'^ ENERGY RANGE ',line1):
                        self.fermi_energy = float(line1.split()[7])*27.2114  
                        #Define from what type of calcualtion the Fermi energy was exctracted
                        self.efermi_from = 'band'
                        break
            #This is in case the .out is from a DOSS calculation  
            if re.match(r'^ TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT DOSS',self.data[len(self.data)-(i+4)]) != None:
                for j,line1 in enumerate(self.data[len(self.data)-i::-1]):
                    if re.match(r'^ N. OF SCF CYCLES ',line1):
                        self.fermi_energy = float(line1.split()[7])*27.2114  
                        #Define from what type of calcualtion the Fermi energy was exctracted
                        self.efermi_from = 'doss'
                        break
            #This is in case the .out is from a sp/optgeom calculation
            #For non metals think about top valence band
            else:      
                for j,line1 in enumerate(self.data[:i:-1]):
                    if re.match(r'^   FERMI ENERGY:',line1) != None:
                        self.fermi_energy = float(line1.split()[2])*27.2114
                        self.efermi_from = 'scf'
                        break
                    if re.match(r'^ POSSIBLY CONDUCTING STATE - EFERMI',line1) != None:
                        self.fermi_energy = float(line1.split()[5]) * 27.2114
                        self.efermi_from = 'scf'
                        break
                if self.fermi_energy == None:
                    for j,line1 in enumerate(self.data[:i:-1]):
                        if re.match(r'^ TOP OF VALENCE BANDS',line1) != None:
                            self.fermi_energy = float(line1.split()[10])*27.2114
                            self.efermi_from = 'scf_top_valence'
                            break
        
        if self.fermi_energy == None:
            print('WARNING: no Fermi energy found in the output file. efermi = None')

        return self.fermi_energy
            
    def get_primitive_lattice(self,initial=True):
        #Initial = False reads the last lattice vectors. Useful in case of optgeom
        import re
        import numpy as np
        
        lattice = []
        self.primitive_lattice = None
        if initial == True:
            for i,line in enumerate(self.data):
                if re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN',line):
                    for j in range(i+2,i+5):
                        lattice_line = [float(n) for n in self.data[j].split()]
                        lattice.append(lattice_line)
                    self.primitive_lattice = np.array(lattice)
                    break
        elif initial == False:
            for i,line in enumerate(self.data[::-1]):                
                if re.match(r'^ DIRECT LATTICE VECTORS CARTESIAN',line):
                    for j in range(len(self.data)-i+1,len(self.data)-i+4):
                        lattice_line = [float(n) for n in self.data[j].split()]
                        lattice.append(lattice_line)
                    self.primitive_lattice = np.array(lattice)
                    break

                
        
        if lattice == []:
            print('WARNING: no lattice vectors found in the output file. lattice = []')
        
        return self.primitive_lattice
    
    def get_reciprocal_lattice(self,initial=True):
        import re
        import numpy as np
        
        lattice = []
        if initial == True:
            for i,line in enumerate(self.data):
                if re.match(r'^ DIRECT LATTICE VECTORS COMPON. \(A.U.\)',line):
                    for j in range(i+2,i+5):
                        lattice_line = [float(n)/0.52917721067121 for n in self.data[j].split()[3:]]
                        lattice.append(lattice_line)
                    self.reciprocal_lattice = np.array(lattice)
                    return self.reciprocal_lattice
        elif initial == False:
            for i,line in enumerate(self.data[::-1]):
                if re.match(r'^ DIRECT LATTICE VECTORS COMPON. \(A.U.\)',line):
                    for j in range(len(self.data)-i+1,len(self.data)-i+4):
                        lattice_line = [float(n)/0.52917721067121 for n in self.data[j].split()[3:]]
                        lattice.append(lattice_line)
                    self.reciprocal_lattice = np.array(lattice)
                    return self.reciprocal_lattice
        
        return None
            

    def get_band_gap(self):#,spin_pol=False):
        import re
        import numpy as np
        
        self.spin_pol = False
        for line in self.data:
            if re.match(r'^ SPIN POLARIZED',line):
                self.spin_pol = True
                break
                
  
        for i, line in enumerate(self.data[len(self.data)::-1]):
            if self.spin_pol == False:
                if re.match(r'^\s\w+\s\w+ BAND GAP',line):
                    self.band_gap = float(line.split()[4])
                    return self.band_gap
                elif re.match(r'^\s\w+ ENERGY BAND GAP',line):
                    self.band_gap = float(line.split()[4])
                    return self.band_gap
                elif re.match(r'^ POSSIBLY CONDUCTING STATE',line):
                    self.band_gap = False
                    return self.band_gap 
            else:
                #This might need some more work
                band_gap_spin = []
                if re.match(r'\s+ BETA \s+ ELECTRONS',line):
                    band_gap_spin.append(float(self.data[len(self.data)-i-3].split()[4]))
                    band_gap_spin.append(float(self.data[len(self.data)-i+3].split()[4]))
                    self.band_gap = np.array(band_gap_spin)
                    return self.band_gap
        if band_gap_spin == []:
            print('DEV WARNING: check this output and the band gap function in file_readwrite')
                #elif re.match(r'^\s\w+ ENERGY BAND GAP',line1) != None:
                    #band_gap = [float(data[len(data)-i-j-7].split()[4]),float(line1.split()[4])]

    def get_last_geom(self,write_gui_file=True,symm_info='pymatgen'):
        import re
        from mendeleev import element
        import numpy as np
        import sys
        from pymatgen.core.structure import Structure
        
        
        self.get_primitive_lattice(initial=False)
        
        self.opt_converged = False
        for line in self.data:
            if re.match(r'^  FINAL OPTIMIZED GEOMETRY',line):
                self.opt_converged = True
                break
        
        for i,line in enumerate(self.data):
            if re.match(r' TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL',line):
                trans_matrix_flat = [float(x) for x in self.data[i+1].split()]
                self.trans_matrix = []
                for i in range(0,len(trans_matrix_flat),3):
                    self.trans_matrix.append(trans_matrix_flat[i:i+3])
                self.trans_matrix = np.array(self.trans_matrix)
        
        
        for i,line in enumerate(self.data[len(self.data)::-1]):
            if re.match(r'^ T = ATOM BELONGING TO THE ASYMMETRIC UNIT',line):
                self.n_atoms = int(self.data[len(self.data)-i-3].split()[0])
                self.atom_positions = [] 
                self.atom_symbols = []
                self.atom_numbers = []
                for j in range(self.n_atoms):
                    atom_line = self.data[len(self.data)-i-2-int(self.n_atoms)+j].split()[3:]
                    self.atom_symbols.append(str(atom_line[0]) )
                    self.atom_positions.append([float(x) for x in atom_line[1:]]) #These are fractional
                a,b,c,alpha,beta,gamma = self.data[len(self.data)-i-2-int(self.n_atoms)-5].split()
                #DELout2cif(file_name,a,b,c,alpha,beta,gamma,atom_positions)
                #DELout_name = str(file_name[:-4]+'.cif')
                for atom in self.atom_symbols:    
                    self.atom_numbers.append(element(atom.capitalize()).atomic_number)
                
                    
                self.atom_positions_cart = np.matmul(np.array(self.atom_positions),self.primitive_lattice)              
                self.cart_coords = []
                for i in range(len(self.atom_numbers)):
                    self.cart_coords.append([self.atom_numbers[i], self.atom_positions_cart[i][0],self.atom_positions_cart[i][1],self.atom_positions_cart[i][2]])
                self.cart_coords = np.array(self.cart_coords)
                
                
                if write_gui_file == True:                                    
                    #Write the gui file 
                    #This is a duplication from write_gui, but the input is different
                    #It requires both the output and gui files with the same name and in the same directory
                    if symm_info == 'pymatgen':
                        if self.name[-3:] == 'out':
                            gui_file = self.name[:-4]+'.gui'
                
                        elif self.name[-4:] == 'outp':
                            gui_file = self.name[:-5]+'.gui'
                        else:
                            gui_file = self.name+'.gui'
                        
                        structure = Structure(self.get_primitive_lattice(initial=False), self.atom_numbers, 
                              self.atom_positions_cart, coords_are_cartesian=True)
                        write_cry_gui(gui_file,structure)
                    else:
                        gui_file = symm_info
                        try:   
                            file = open(gui_file, 'r')
                            gui_data = file.readlines()
                            file.close()
                        except:
                            print('EXITING: a .gui file with the same name as the input need to be present in the directory.')
                            sys.exit(1)
                            
                        #Replace the lattice vectors with the optimised ones
                        for i,vector in enumerate(self.get_primitive_lattice(initial=False).tolist()):
                            gui_data[i+1]= ' '.join([str(x) for x in vector])+'\n'
        
                        n_symmops = int(gui_data[4])
                        for i in range(len(self.atom_numbers)):
                            gui_data[i+n_symmops*4+6] = '{} {}\n'.format(self.atom_numbers[i],' '.join(str(x) for x in self.atom_positions_cart[i][:]))
                        
                        with open(gui_file[:-4]+'_last.gui','w') as file:
                            for line in gui_data:
                                file.writelines(line)
        self.last_geom = [self.primitive_lattice.tolist(), self.atom_numbers, self.atom_positions_cart.tolist()]
        return self.last_geom
        
    def get_symm_ops(self):
        import re
        import numpy as np
        
        symmops = []

        for i,line in enumerate(self.data):
            if re.match(r'^ \*\*\*\*   \d+ SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS',line):
                self.n_symm_ops = int(line.split()[1])
                for j in range(0,self.n_symm_ops):
                    symmops.append(self.data[i+3+j].split()[2:])
                self.symm_ops = np.array(symmops)
                
                return self.symm_ops
    
    def get_forces(self,initial=False,grad=False):
        import re
        import numpy as np
        
        self.forces_atoms = []
        self.forces_cell = []
        #Number of atoms
        for i,line in enumerate(self.data[len(self.data)::-1]):
            if re.match(r'^ T = ATOM BELONGING TO THE ASYMMETRIC UNIT',line):
                self.n_atoms = int(self.data[len(self.data)-i-3].split()[0])
                break
        
        if grad == True:
            self.grad = []
            self.rms_grad = []
            self.disp = []
            self.rms_disp = []
            for i,line in enumerate(self.data):
                if re.match(r'^ MAX GRADIENT',line):
                    self.grad.append(line.split()[2])
                if re.match(r'^ RMS GRADIENT',line):
                    self.rms_grad.append(line.split()[2])
                if re.match(r'^ MAX DISPLAC.',line):
                    self.disp.append(line.split()[2])
                if re.match(r'^ RMS DISPLAC.',line):
                    self.rms_disp.append(line.split()[2])
                 
        if initial == True:
            for i,line in enumerate(self.data):
                if re.match(r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)',line):
                    for j in range(i+2,i+2+self.n_atoms):
                        self.forces_atoms.append([float(x) for x in self.data[j].split()[2:]])
                    self.forces_atoms = np.array(self.forces_atoms)
                if re.match(r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR',line):
                    for j in range(i+4,i+7):
                        self.forces_cell.append([float(x) for x in self.data[j].split()])
                    self.forces_cell = np.array(self.forces_cell)
                    self.forces = [self.forces_cell,self.forces_atoms]
                    return self.forces
                
        elif initial == False:
            for i,line in enumerate(self.data[::-1]):
                if re.match(r'^ GRADIENT WITH RESPECT TO THE CELL PARAMETER IN HARTREE/BOHR',line):
                    for j in range(len(self.data)-i+3,len(self.data)-i+6):
                        self.forces_cell.append([float(x) for x in self.data[j].split()])
                    self.forces_cell = np.array(self.forces_cell)
                
                if re.match(r'^ CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\)',line):
                    for j in range(len(self.data)-i+1,len(self.data)-i+1+self.n_atoms):
                        self.forces_atoms.append([float(x) for x in self.data[j].split()[2:]])
                    self.forces_atoms = np.array(self.forces_atoms)
                    self.forces = [self.forces_cell,self.forces_atoms]
                    return self.forces
        

    def get_mulliken_charges(self):
        
        import re

        self.mulliken_charges = []
        for i,line in enumerate(self.data):
                if re.match(r'^ MULLIKEN POPULATION ANALYSIS',line):
                    for j in range(len(self.data[i:])):
                        line1 = self.data[i+4+j].split()
                        if line1 == []:
                            return self.mulliken_charges
                        elif line1[0].isdigit() == True:
                            self.mulliken_charges.append(float(line1[3]))
                              

    def get_config_analysis(self):
        import re
        import numpy as np

        try:
            begin = self.data.index('                             CONFIGURATION ANALYSIS\n')
        except:
            return "WARNING: this is not a CONFCNT analysis."

        for i,line in enumerate(self.data[begin:]):
            if re.match(r'^ COMPOSITION', line):
                self.n_classes = line.split()[9]
                original_atom = str(line.split()[2])
                begin = begin+i

        non_subs = []
        subs = []
        class_index = []
        config_list = []

        for line in self.data[begin:]:
            if not re.match(r'^   WARNING', line):
                config_list.extend(line.split())
        config_list = np.array(config_list)
        warning = np.where(config_list == 'WARNING')
        config_list = np.delete(config_list,warning)
        atom1_begin = np.where(config_list == original_atom)[0]
        atom1_end = np.where(config_list == '--------------------------------------------------------------------------')[0]
        atom2_begin = np.where(config_list == 'XX')[0]
        atom2_end = np.where(config_list == '<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>')[0]
        end = np.where(config_list == '===============================================================================')[0][-1]
        atom2_end = np.append(atom2_end,end)
        atom_type1 = []
        atom_type2 = []
        config_list = config_list.tolist()
        for i in range(len(atom1_end)):
            atom_type1.append([int(x) for x in config_list[atom1_begin[i+1]+1:atom1_end[i]]])
            atom_type2.append([int(x) for x in config_list[atom2_begin[i]+1:atom2_end[i]]])

        self.atom_type1 = atom_type1
        self.atom_type2 = atom_type2
        return [self.atom_type1, self.atom_type2]
        


###TESTING
#a = Crystal_output('../examples/data/mgo_optgeom.out')
#a = Crystal_output('../examples/data/LTS_CONFCNT_ONLY')
#print(a.config_analysis())
#print('final_energy\n',a.final_energy())
#print('fermi\n',a.fermi_energy())
#print('primitive\n',a.primitive_lattice(initial=False))
#print('band_gap\n',a.band_gap())
#print('spin\n',a.spin_pol)
#print('reciprocal\n',a.reciprocal_lattice())
#print('last geom\n',a.extract_last_geom(print_cart=False))
#print('symmops\n',a.symm_ops())
#print('forces\n',a.forces(gradient=True))
#print('grad\n',a.grad)
#print('scf convergence\n',a.scf_convergence(all_cycles=True))'''
#print(a.get_mulliken_charges())

class Properties_input:
    #This creates a crystal_input object
    
    def __init__(self,input_name=None):
        import sys
        
        self.name = input_name
        if input_name is not None:
            try: 
                if input_name[-3:] != 'd12':
                    input_name = input_name+'.d12'
                file = open(input_name, 'r')
                self.data = file.readlines()
                file.close()
            except:
                print('EXITING: a .d3 file needs to be specified')
                sys.exit(1)
        
            if 'NEWK\n' in self.data:
                self.is_newk = True
                self.newk_block = self.data[0:2]
                self.property_block = self.data[2:]
            else:
                self.is_newk = False
                self.property_block = self.data

    def make_newk_block(self,shrink1,shrink2,Fermi=1,print_option=0,title=None):
        self.newk_block = ['NEWK\n','%s %s\n' %(shrink1,shrink2),
                    '%s %s\n'%(Fermi,print_option)]

        return self.newk_block


    def make_bands_block(self,k_path,n_kpoints,first_band,last_band,print_eig=0,print_option=1,
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
        
        self.bands_block = bands_block

        return self.bands_block      

    def make_doss_block(self,n_points=200,band_range=None,e_range=None,plotting_option=2,
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
        
        self.doss_block = doss_block

        return self.doss_block
    
    def make_pdoss_block(self,projections,proj_type='atom',output_file=None,n_points=200,band_range=None,
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
        
        self.pdoss_block = pdoss_block
        return self.pdoss_block
     

class Properties_output:

    def __init__(self,properties_output):
       
        import sys

        self.file_name = properties_output

        try: 
            file = open(self.file_name, 'r')
            self.data = file.readlines()
            file.close()
        except:
            print('EXITING: a CRYSTAL properties file needs to be specified')
            sys.exit(1)

    def read_cry_bands(self):
        #This class contains the bands objects created from reading the 
        #band files created by different electronic structure codes
        #Returns an array where the band energy is expressed in eV

        import re
        import numpy as np

        data = self.data
            
        #Read the information about the file
        self.n_kpoints = int(data[0].split()[2])
        self.n_bands = int(data[0].split()[4])
        self.spin = int(data[0].split()[6])
        self.n_tick = int(data[1].split()[2])+1
        self.k_point_inp_coordinates = []
        self.n_points = []
        for i in range(self.n_tick):
            self.n_points.append(int(data[2+i].split()[1]))
            coord = []
            for j in range(3):
                l = re.findall('\d+',data[2+i].split()[2])
                coord.append(float(l[j])/float(l[3]))
            self.k_point_inp_coordinates.append(coord)  
        self.k_point_inp_coordinates = np.array(self.k_point_inp_coordinates)
        self.k_point_coordinates = [self.k_point_inp_coordinates[0]]
        for i in range(1,self.n_tick):
            step = (self.k_point_inp_coordinates[i]-self.k_point_inp_coordinates[i-1])/float(self.n_points[i]-self.n_points[i-1])
            for j in range(self.n_points[i]-self.n_points[i-1]):
                self.k_point_coordinates.append((self.k_point_inp_coordinates[i-1]+step*float(j+1)).tolist())
        self.tick_position = []
        self.tick_label = []
        for i in range(self.n_tick):            
            self.tick_position.append(float(data[16+self.n_tick+i*2].split()[4]))
            self.tick_label.append(str(data[17+self.n_tick+i*2].split()[3][2:]))
        self.efermi = float(data[-1].split()[3])*27.2114
        
        #Allocate the bands as np arrays
        self.bands = np.zeros((self.n_bands,self.n_kpoints,self.spin),dtype=float)        
        
        
        #line where the first band is. Written this way to help identify
        #where the error might be if there are different file lenghts
        first_k = 2 + self.n_tick + 14 + 2*self.n_tick + 2   
        
        #Read the bands and store them into a numpy array
        for i,line in enumerate(data[first_k:first_k+self.n_kpoints]):
            self.bands[:self.n_bands+1,i,0] = np.array([float(n) for n in line.split()[1:]])
        
        if self.spin == 2:
            #line where the first beta band is. Written this way to help identify
            first_k_beta = first_k + self.n_kpoints + 15 + 2*self.n_tick + 2
            for i,line in enumerate(data[first_k_beta:-1]):
                self.bands[:self.n_bands+1,i,1] = np.array([float(n) for n in line.split()[1:]])
        
        #Convert all the energy to eV    
        self.bands[1:,:,:] = self.bands[1:,:,:]*27.2114
        return self  
    
    '''###TESTING
    mgo_bands = Bands('data/mgo_BAND_dat.BAND') 
    mgo_file = mgo_bands.read_cry_band()
    print(mgo_file.bands[-1,0,0])
    print(mgo_file.tick_position[-1])'''

    def read_cry_doss(self):
        
        import re
        import numpy as np

        data = self.data

        #Read the information about the file
        self.n_energy = int(data[0].split()[2])
        self.n_proj = int(data[0].split()[4])
        self.spin = int(data[0].split()[6])
        self.efermi = float(data[-1].split()[3])*27.2114
        
        first_energy = 4
        
        #Allocate the doss as np arrays
        self.doss = np.zeros((self.n_energy,self.n_proj+1,self.spin),dtype=float)        
        
        #Read the doss and store them into a numpy array
        for i,line in enumerate(data[first_energy:first_energy+self.n_energy]):
            self.doss[i,:self.n_proj+1,0] = np.array([float(n) for n in line.split()])
            
        if self.spin == 2:
            #line where the first beta energy is. Written this way to help identify
            first_energy_beta = first_energy + self.n_energy + 3
            for i,line in enumerate(data[first_energy_beta:-1]):
                self.doss[i,:self.n_proj+1,1] = np.array([float(n) for n in line.split()])
        
        #Convert all the energy to eV    
        self.doss[0,:,:] = self.doss[0,:,:]*27.2114
        
        return self 
    
    def read_cry_contour(self):
        #Functions that extract useful info and stores it in sel.attributes
        pass
    

        
###TESTING
#contour_obj = Crystal_properties('../examples/data/SURFRHOO.DAT').read_cry_contour()
#a = contour_obj.read_cry_contour()
#bands = Crystal_properties('../examples/data/mgo_BAND_dat.BAND').read_cry_bands()
#print(bands.bands)



def write_crystal_input(input_name,crystal_input=None,crystal_blocks=None,external_obj=None,comment=None):
    #input_name is the name of the imput that is going to be written (.d12)
    #crystal_input is an object belonging to the crystal_input Class.
    #external_obj is the ASE or pymatgen object that is going to be written in the fort.34
    
    import itertools
    import sys
    import re
    from ase.io.crystal import write_crystal
    from pymatgen.io.ase import AseAtomsAdaptor
    
    if (crystal_input == None and crystal_blocks == None) or (crystal_input == True and crystal_blocks == True):
        print('EXITING: please specify either a CRYSTAL input or CRYSTAL input blocks')
        sys.exit(1)
    if crystal_blocks != None:
        if type(crystal_blocks) == list and len(crystal_blocks) == 4:
            geom_block = crystal_blocks[0]
            bs_block = crystal_blocks[1]
            func_block = crystal_blocks[2]
            scf_block = crystal_blocks[3]
        else:
            print('EXITING: the CRYSTAL blocks are not in the correct format.')
            sys.exit(1)
    elif crystal_input != None:
        geom_block = crystal_input.geom_block
        bs_block = crystal_input.bs_block
        func_block = crystal_input.func_block
        scf_block = crystal_input.scf_block
    #print(geom_block,bs_block,func_block,scf_block)    

        
    
    #if there is an external object, we want to have the EXTERNAL
    #keyword in the geom_block. If it's not present, this means
    #adding it and keeping the rest of the input
    if external_obj != None:    
        if 'EXTERNAL\n' not in geom_block:
            new_geom_block = []
            if comment == None:
                new_geom_block.append(geom_block[0])
            else:
                new_geom_block.append(comment+'\n')
            new_geom_block.append('EXTERNAL\n')
            for i,line in enumerate(geom_block[2:]):
                if line.split()[0].replace('.','',1).isdigit() == False:
                    for line1 in geom_block[i+2:]:
                        new_geom_block.append(line1)
                    break      
            geom_block = new_geom_block
        if 'ase.atoms' in  str(type(external_obj)):
            write_crystal(input_name[:-4]+'.gui',external_obj)
        elif 'pymatgen.core' in str(type(external_obj)):
            gui_file_name = input_name[:-4]+'.gui'
            write_cry_gui(gui_file_name, external_obj)
        else:
            print('EXITING: external object format not recognised, please specfy an ASE or pymatgen object')
            sys.exit(1)
        
    with open(input_name, 'w') as file:
        cry_input = list(itertools.chain(geom_block, bs_block,
                                         func_block,scf_block))
        for line in cry_input:
            file.writelines(line)

###TESTING
'''from pymatgen.core import Structure, Lattice             
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
substrate = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.61491), ["Cu"], [[0, 0, 0]])
substrate_conv = SpacegroupAnalyzer().get_conventional_standard_structure() 

#mgo = Crystal_input('examples/data/mgo.d12') 
#write_cry_input('examples/data/mgo_TEST.d12',crystal_input = mgo,external_obj=substrate_conv,comment='YES')

mgo = Crystal_input('examples/data/mgo.d12') 
print(mgo.geom_block)
write_cry_input('examples/data/mgo_TEST.d12',crystal_blocks= [mgo.geom_block,mgo.bs_block,mgo.func_block,mgo.scf_block],external_obj=substrate_conv,comment='YES')'''

def write_properties_input(input_name,property_block,newk=False):
    
    import sys
    import itertools
    
    if newk == False:
        property_input = property_block
    if newk != False and type(newk) != list:
        print('EXITING: newk must be a newk_block list')
        sys.exit(1)
    elif type(newk) == list:
        property_input = list(itertools.chain(property_block,newk))
            
        
    with open(input_name, 'w') as file:        
        for line in property_input:
            file.writelines(line)
            
def write_cry_gui(gui_file,atoms,dimensionality=3,symm=True):
    #gui_file is the name of the gui that is going to be written (including .gui)
    #atoms is the structure object from ase or pymatgen
    
    from ase.io.crystal import write_crystal
    
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import center_slab
    
    import numpy as np

    #ASE object
    if 'ase.atoms' in  str(type(atoms)):
        write_crystal(gui_file,atoms)
        
    #pymatgen object
    elif 'pymatgen' in str(type(atoms)):
        if 'Slab' in str(type(atoms)):
            dimensionality = 2

    #Save the coordinates before symmetry transformations
    atomic_numbers = []
    atom_coords = []
    for i in range(atoms.num_sites):
        atomic_numbers.append(atoms.atomic_numbers[i])
        atom_coords.append(' '.join(str(np.around(n, 5)) for n in atoms.cart_coords[i]))

    with open(gui_file, 'w') as file:
        if symm == True:
            try:
                atoms = SpacegroupAnalyzer(atoms).get_primitive_standard_structure()
                atomic_numbers = []
                atom_coords = []
                for i in range(atoms.num_sites):
                    atomic_numbers.append(atoms.atomic_numbers[i])
                    atom_coords.append(' '.join(str(np.around(n, 5)) for n in atoms.cart_coords[i]))
                #Is the structure symmetrysed?
                if 'SymmetrizedStructure' not in str(type(atoms)):
                    atoms_symm = SpacegroupAnalyzer(atoms).get_symmetrized_structure()
                symmetry = True
            except:
                atoms_symm = atoms
                symmetry = False
        elif symm == False:
            atoms_symm = atoms
            symmetry = False

        if dimensionality == 3:

            #First line
            file.writelines('3   1   1\n')

            #Cell vectors
            for vector in atoms.lattice.matrix:
                file.writelines(' '.join(str(n) for n in vector)+'\n')

            #N symm ops
            if symmetry == True:
                n_symmops = len(SpacegroupAnalyzer(atoms_symm).get_space_group_operations())

                file.writelines('{}\n'.format(str(n_symmops)))

                #symm ops
                for symmops in SpacegroupAnalyzer(atoms_symm).get_symmetry_operations(cartesian=True):
                    file.writelines('{}\n'.format(' '.join(str(np.around(n,8)) for n in symmops.rotation_matrix[0])))
                    file.writelines('{}\n'.format(' '.join(str(np.around(n,8)) for n in symmops.rotation_matrix[1])))
                    file.writelines('{}\n'.format(' '.join(str(np.around(n,8)) for n in symmops.rotation_matrix[2])))
                    file.writelines('{}\n'.format(' '.join(str(np.around(n,8)) for n in symmops.translation_vector)))
            else:
                n_symmops = 1
                # symm ops
                file.writelines('{}\n'.format(str(n_symmops)))
                identity = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,0.]])
                file.writelines('{}\n'.format(' '.join(str(np.around(n,8)) for n in identity[0])))
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[1])))
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[2])))
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[3])))

        elif dimensionality == 2:
            file.writelines('2   1   1\n')
            #Cell vectors
            z_component = ['0.0000', '0.00000', '500.00000']
            for i,vector in enumerate(atoms.lattice.matrix[0:3,0:2]):
                file.writelines(' '.join(str(n) for n in vector)+' '+z_component[i]+'\n')
            #Center the slab
            #First center at z = 0.5
            atoms = center_slab(atoms)

            #Then center at z=0.0
            translation = np.array([0.0, 0.0, -0.5])
            atoms.translate_sites(list(range(atoms.num_sites)),
                                         translation, to_unit_cell=False)

            if symmetry == True:
                #Remove symmops with z component
                sg = SpacegroupAnalyzer(atoms)
                ops = sg.get_symmetry_operations(cartesian=True)
                symmops = []
                for op in ops:
                    if op.translation_vector[2] == 0.:
                        symmops.extend(op.rotation_matrix.tolist())
                        symmops.extend([op.translation_vector.tolist()])

                #N symm ops
                n_symmops = int(len(symmops)/4)
                file.writelines('{}\n'.format(n_symmops))

                #symm ops
                for symmop in symmops:
                    file.writelines('{}\n'.format(' '.join(str(np.around(n,8)) for n in symmop)))
            else:
                n_symmops = 1
                # symm ops
                identity = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 0., 0.]])
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[0])))
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[1])))
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[2])))
                file.writelines('{}\n'.format(' '.join(str(np.around(n, 8)) for n in identity[3])))

        if symmetry == True:
            # N atoms
            file.writelines('{}\n'.format(len(atoms_symm.equivalent_indices)))

            # atom number + coordinates cart
            for i in range(len(atoms_symm.equivalent_indices)):
                file.writelines('{} {}\n'.format(atomic_numbers[atoms_symm.equivalent_indices[i][0]],
                                                 atom_coords[atoms_symm.equivalent_indices[i][0]]))
        else:
            #N atoms
            file.writelines('{}\n'.format(atoms.num_sites))

            #atom number + coordinates cart
            for i in range(atoms.num_sites):
                file.writelines('{} {}\n'.format(atomic_numbers[i],atom_coords[i]))

        #space group + n symm ops
        if symmetry == True:
            file.writelines('{} {}'.format(SpacegroupAnalyzer(atoms).get_space_group_number(),len(SpacegroupAnalyzer(atoms).get_space_group_operations())))
        else:
            file.writelines('1 1')
            
'''###TESTING
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator


#substrate = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.597), ["Cu"], [[0, 0, 0]])
substrate = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.217), ["Mg",'O'], [[0, 0, 0],[0.5,0.5,0.5]])
substrate = SlabGenerator(substrate, (1,0,0), 5., 10., center_slab=False).get_slab()


write_cry_gui('examples/data/cu_100_TEST.gui',substrate,dimensionality=2)'''
    

class Crystal_density:
    
    def __init__(self, fort98_unit):
        self.file_name = fort98_unit
    
    #def cry_read_density(self):
            
        import sys
        import numpy as np
        import re
    
        try: 
            file = open(self.file_name, 'r')
            data = file.readlines()
            file.close()
        except:
            print('EXITING: a CRYSTAL .f98 file needs to be specified')
            sys.exit(1)

        self.all_file = data

        #the keyword BASATO appears twice, this is a check to see which basato
        #is being read
        basato1 = False
        
        for i,line in enumerate(data):
            
            if re.match(r'^LIMINF LIMTOL LIMPAR', line):
                  inf_vec_len, tol_vec_len, par_vec_len = [int(x) for x in data[i+1].split()]
                
            elif re.match(r'^INF', line):
                self.inf_vec = []
                inf_n_lines = int(np.ceil(inf_vec_len/8))
                for j in range(inf_n_lines):
                    self.inf_vec.extend([int(x) for x in data[i+1+j].split()])   
                n_symmops = self.inf_vec[0]
                n_atoms = self.inf_vec[23]
                n_shells = self.inf_vec[19]
                '''if self.inf_vec[26] == 0:
                    self.spin_pol == False
                elif self.inf_vec[26] == 1:
                    self.spin_pol == True'''
                n_prim_gto = self.inf_vec[74]
                f_irr_len = (self.inf_vec[63]+1)*self.inf_vec[227]
                p_irr_len = (self.inf_vec[63]+1)*self.inf_vec[18]
                nnnc_len = self.inf_vec[190]
                la3_len = self.inf_vec[55]
                #n_symmops_noinv = inf_vec[1]
            
            elif re.match(r'^TOL', line):
                self.tol_vec = []
                tol_n_lines = int(np.ceil(tol_vec_len/8))
                for j in range(tol_n_lines):
                    self.tol_vec.extend([int(x) for x in data[i+1+j].split()])
            
            elif re.match(r'^PAR', line):
                self.par_vec = []
                par_n_lines = int(np.ceil(par_vec_len/4))
                for j in range(par_n_lines):
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #The line below fixes that issue
                    for item in range(0,int(len(data[i+1+j])/20)):
                        self.par_vec.append(float(data[i+1+j][(item)*20:(item+1)*20])) 
                        
            
            elif re.match(r'^XYVGVE', line):
                #This vector contains the rotations, translation,
                #lattice vectors and transformation matrix from primitive to 
                #crystallographic cell 
                #Read all of it first and separate later
                xyvgve_n_lines = int(np.ceil((n_symmops*12+18)/4))
                xyvgve_vec = []
                for j in range(xyvgve_n_lines):
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #The line below fixes that issue
                    for item in range(0,int(len(data[i+1+j])/20)):
                        xyvgve_vec.append(float(data[i+1+j][(item)*20:(item+1)*20]))  
                #Now let's split the xyvgve_vec
                self.rotations_vec = xyvgve_vec[0:n_symmops*9]
                self.translations_vec = xyvgve_vec[n_symmops*9:n_symmops*9+n_symmops*3]
                self.direct_lattice_vec = xyvgve_vec[n_symmops*12:n_symmops*12+9]
                self.transf_matrix = xyvgve_vec[-9:]
                
            elif re.match(r'^BASATO', line):
                if basato1 == False:
                    basato_n_lines = int(np.ceil((n_atoms*4+n_shells*5+n_prim_gto*7)/4))
                    basato_vec = []
                    for j in range(basato_n_lines):
                        #The negative elements appear connected to the previous one
                        #eg:  0.0000000000000E+00-1.0000000000000E+00
                        #The line below fixes that issue
                        for item in range(0,int(len(data[i+1+j])/20)):
                            basato_vec.append(float(data[i+1+j][(item)*20:(item+1)*20]))  
                    #Extract the iformation we need from basato
                    
                    #Atom coordinates
                    self.atom_coord = []
                    for j in range(0,3*n_atoms,3):
                        self.atom_coord.append(basato_vec[(n_atoms+j):(n_atoms+j+3)])
                    #self.atom_coord = np.array(self.atom_coord)
                    
                    #Assign the shell to the atom
                    self.shell_coord = []
                    #The beginning of the part of BASATO I need here
                    init = 4*n_atoms + 2*n_shells 
                    for j in range(0,3*n_shells,3):
                        self.shell_coord.append(basato_vec[(init+j):(init+j+3)])
                    #self.shell_coord = np.array(self.shell_coord)
                    
                    #Array that defines which atom a shell belongs to
                    self.shell_to_atom = []
                    for coord in self.shell_coord:
                        self.shell_to_atom.append(self.atom_coord.index(coord))
                    basato1 = True

                      
                
                elif basato1 == True:
                    self.basato2 = []
                    j = i + 1
                    while 'SPINOR' not in data[j].split()[0]:
                        # The negative elements appear connected to the previous one
                        # eg:  0.0000000000000E+00-1.0000000000000E+00
                        # As opposite to the loops above where the float read was 20
                        # characters long, this ones are 21
                        # The line below fixes that issue
                        self.basato2.extend([int(x) for x in data[j].split()])
                        j += 1
                    self.atom_shell = self.basato2[-n_shells:]
                                  
            elif re.match(r'^SPINOR', line):
                self.f_irr = []
                self.charges = []
                #self.spin = [0]*(2*n_atoms) #tmp
                self.spin = []
                self.ghost = []
                n_ghost = 0
                n_spin_lines = int(np.ceil((n_atoms*2)/8))
                n_basold = 0
                if 'BASOLD' in data[i+n_spin_lines+1]:
                    n_basold = 9+ 3*self.inf_vec[1]+ n_shells +3*n_atoms+3*n_shells+self.inf_vec[4]+1+3*self.inf_vec[78]
                n_basold_lines = int(np.ceil((n_basold)/4))
                n_charge_lines = int(np.ceil(n_atoms/4))
                skip = n_spin_lines + n_charge_lines +1 + n_basold_lines
                for j in range(n_spin_lines):
                    self.spin.extend([int(x) for x in data[i + j + 1].split()])
                if 'IGHOST' in data[i+n_spin_lines+1]:
                    n_ghost = int(np.ceil((n_atoms)/8)) +1
                    skip = skip + n_ghost
                    for j in range(n_ghost-1):
                        self.ghost.extend([float(x) for x in data[i + j + n_spin_lines + 2].split()])
                f_irr_n_lines = int(np.ceil(f_irr_len/4))
                for j in range(n_charge_lines):
                    self.charges.extend([float(x) for x in data[i+j+n_spin_lines+n_basold+n_ghost+1].split()])
                    '''for item in range(0,int(len(data[i+skip+j])/21)):
                        print(data[i+j+n_spin_lines+1])
                        self.charges.append(float(data[i+j+n_spin_lines+1][(item)*21:(item+1)*21]))'''
                for j in range(f_irr_n_lines):
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #As opposite to the loops above where the float read was 20
                    #characters long, this ones are 21
                    #The line below fixes that issue
                    for item in range(0,int(len(data[i+skip+j])/21)):
                        self.f_irr.append(float(data[i+skip+j][(item)*21:(item+1)*21]))    
                self.p_irr = []
                p_irr_n_lines = int(np.ceil(p_irr_len/4))
                skip += 1
                for k in range(i+skip+j,i+skip+j+p_irr_n_lines):
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #As opposite to the loops above where the float read was 20
                    #characters long, this ones are 21
                    #The line below fixes that issue
                    for item in range(0,int(len(data[k])/21)):
                        self.p_irr.append(float(data[k][(item)*21:(item+1)*21]))  
                        
            elif re.match(r'^   NCF', line):
                #The ncf vector contains the pointers to the symmetry irerducible
                #shell couples la3, la4
                self.ncf = []
                j = i+1
                while 'NSTATG' not in data[j].split()[0]:
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #As opposite to the loops above where the float read was 20
                    #characters long, this ones are 21
                    #The line below fixes that issue
                    self.ncf.extend([int(x) for x in data[j].split()])
                    j += 1

            elif re.match(r'^NSTATG', line):
                #The nstatg vector contains the pointers to the starting point
                #of each couple set in P_irr and F_irr
                self.nstatg = []
                j = i+1
                while 'NSTAFG' not in data[j].split()[0]:
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #As opposite to the loops above where the float read was 20
                    #characters long, this ones are 21
                    #The line below fixes that issue
                    self.nstatg.extend([int(x) for x in data[j].split()])
                    j += 1
                    
            elif re.match(r'^  NNNC', line):
                #The nnnc points the starting position in the P matrix for 
                #each couple, and its size corresponds to the total number 
                #of shell couple in the shell couple sets
                self.nnnc = []
                nnnc_n_lines = int(np.ceil(nnnc_len/8))
                for j in range(nnnc_n_lines):
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #As opposite to the loops above where the float read was 20
                    #characters long, this ones are 21
                    #The line below fixes that issue
                    for item in range(0,int(len(data[i+1+j])/10)):
                        self.nnnc.append(int(data[i+1+j][(item)*10:(item+1)*10]))  
            elif re.match(r'^   LA3', line):
               #The nnnc points the starting position in the P matrix for 
               #each couple, and its size corresponds to the total number 
               #of shell couple in the shell couple sets
               self.la3 = []
               #nnnc_n_lines = int(np.ceil(nnnc_len/8))
               j = i+1 
               while 'LA4' not in data[j].split()[0]:
                   #The negative elements appear connected to the previous one
                   #eg:  0.0000000000000E+00-1.0000000000000E+00
                   #As opposite to the loops above where the float read was 20
                   #characters long, this ones are 21
                   #The line below fixes that issue
                   self.la3.extend([int(x) for x in data[j].split()])
                   j += 1   
            elif re.match(r'^   LA4', line):
                #The nnnc points the starting position in the P matrix for 
                #each couple, and its size corresponds to the total number 
                #of shell couple in the shell couple sets
                self.la4 = []
                #nnnc_n_lines = int(np.ceil(nnnc_len/8))
                j = i+1
                while 'IROF' not in data[j].split()[0]:
                    #The negative elements appear connected to the previous one
                    #eg:  0.0000000000000E+00-1.0000000000000E+00
                    #As opposite to the loops above where the float read was 20
                    #characters long, this ones are 21
                    #The line below fixes that issue
                    self.la4.extend([int(x) for x in data[j].split()])
                    j += 1           
        return self
                                     
            #elif re.match('', line):
            #elif re.match('', line):
        
                
def cry_combine_density(density1,density2,density3,new_density='new_density.f98',spin_pol=False):
    import sys
    import numpy as np
    
    try:
        density1_data = Density(density1).cry_read_density() #substrate
        density2_data = Density(density2).cry_read_density()
        ###density3_data_obj = Density(density3).cry_read_density()
        density3_data_obj = Density(density1).cry_read_density()
        file = open(density3, 'r')
        density3_data = file.readlines()
        file.close()
    except:
        print('EXITING: a CRYSTAL .f98 file needs to be specified')
        sys.exit(1)

    #Find P_irr <-> atom correspondence
    fragment_1 = []
    fragment_2 = []

    for i,j in enumerate(density1_data.ncf):
        #density1_data.la3[j] is the shell number
        #density1_data.atom_shell[density1_data.la3[j]] is the atom position number (1-6)
        #density1_data.ghost[density1_data.atom_shell[density1_data.la3[j]]] is either 0 or atomic number depending on ghost or not
        #This tells me if the shell belongs to this fragment
        n_elements = density1_data.nstatg[i] - density1_data.nstatg[i - 1]

        #print(i, j, len(density1_data.la3), len(density1_data.la4), len(density1_data.atom_shell),
              #len(density1_data.ghost))
        '''print(i, j,
              density1_data.la3[j-1],
              density1_data.la4[j-1],
              density1_data.atom_shell[density1_data.la3[j-1]-1],
              density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1])'''
        if density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1] == 0 and \
            density1_data.ghost[density1_data.atom_shell[density1_data.la4[j-1]-1]-1] == 0:
            fragment_1.extend([True]*n_elements)
        else:
            fragment_1.extend([False]*n_elements)
        if density1_data.ghost[density1_data.atom_shell[density1_data.la3[j-1]-1]-1] != 0 and \
            density1_data.ghost[density1_data.atom_shell[density1_data.la4[j-1]-1]-1] != 0:
            fragment_2.extend([True]*n_elements)
        else:
            fragment_2.extend([False]*n_elements)

    if spin_pol == True:
        spin_p1 = fragment_1.copy()
        spin_p2 = fragment_2.copy()
        fragment_1.extend(spin_p1)
        fragment_2.extend(spin_p2)
    '''for i in range(len(fragment_2)):
        if fragment_1[i] == True and fragment_2==True:
            print(fragment_1[i],fragment_2[i])'''
    beginning = density3_data.index('SPINOR\n')
    end = density3_data.index('   NCF\n')
    sum_density = np.array(density1_data.p_irr)+np.array(density2_data.p_irr)
    sum_fock = np.array(density1_data.f_irr)+np.array(density2_data.f_irr)
    sum_charges = np.array(density1_data.charges)+np.array(density2_data.charges)
    #sum_charges = [12.,12.,8.,8.,8,8.]
    spinor = ['SPINOR\n']
    charges = []
    fock = []
    density = []
    new_fock = sum_fock #TMP
    new_fock = [0] * len(density3_data_obj.f_irr)
    new_p = []

    for i in range(len(fragment_1)):
        if fragment_1[i] == True and fragment_2[i] == False:
            #new_fock.append(density1_data.f_irr[i])
            new_p.append(density1_data.p_irr[i])
        elif fragment_1[i] == False and fragment_2[i] == True:
            ##new_fock.append(density2_data.f_irr[i])
            new_p.append(density2_data.p_irr[i])
        elif fragment_1[i] == False and fragment_2[i] == False:
            #new_fock.append(0.)
            new_p.append(sum_density[i])
            #new_p.append(0.)
            #new_p.append(density3_data_obj.p_irr[i])
    #print('fock',len(density3_data_obj.f_irr),len(new_fock))
    #print('density',len(density3_data_obj.p_irr),len(new_p))
    for i in range(0,len(density3_data_obj.spin),8):
        spinor.append(' '.join([str(x) for x in density3_data_obj.spin[i:i+8]])+'\n')
    for i in range(0,len(sum_charges),4):
        charges.append(' '.join(["{:.13e}".format(x) for x in sum_charges[i:i+4]])+'\n')
    for i in range(0,len(new_fock),4):
        fock.append(' '.join(["{:.13e}".format(x) for x in new_fock[i:i+4]])+'\n')
    for i in range(0,len(new_p),4):
        density.append(' '.join(["{:.13e}".format(x) for x in new_p[i:i+4]])+'\n')
    #print(len(fock),len(density))
    final_fort98 = density3_data[0:beginning]+spinor+charges+fock+density+density3_data[end:] 
    with open(new_density, 'w') as file:
        for line in final_fort98:
            file.writelines(line)

def write_cry_density(fort98_name,new_p,new_fort98):

    import numpy as np

    file = open(fort98_name, 'r')
    data = file.readlines()
    file.close()

    density = Density(fort98_name).cry_read_density()

    n_spin_lines = int(np.ceil((density.inf_vec[23] * 2) / 8))
    n_charges_lines = int(np.ceil((density.inf_vec[23] ) / 4))
    beginning = data.index('SPINOR\n') + n_spin_lines + n_charges_lines + 1
    end = data.index('   NCF\n')



    new_fock_vect = [0] * len(density.f_irr)

    new_fock = []
    for i in range(0, len(new_fock_vect), 4):
        new_fock.append(' '.join(["{:.13e}".format(x) for x in new_fock_vect[i:i + 4]]) + '\n')

    new_density = []
    for i in range(0,len(new_p),4):
        new_density.append(' '.join(["{:.13e}".format(x) for x in new_p[i:i+4]])+'\n')


    final_fort98 = data[0:beginning]+new_fock+new_density+data[end:]
    with open(new_fort98, 'w') as file:
        for line in final_fort98:
            file.writelines(line)
            
###TESTING
#H_density =  Density('examples/data/h_bulk.f98').cry_read_density()
#H_density =  Density('examples/data/mgo.f98').cry_read_density()     
#ghost_density = Density('../examples/data/Mg2O2_O1_100_1.f98').cry_read_density()
#ghost_density = Density('../examples/data/Mg2O2_O1_100_1_BSSE_ads.f98').cry_read_density()
#print(ghost_density.charges)
#density_basold = Density('/Users/brunocamino/Desktop/Imperial/cmsg_icl/solid-solutions/data/classification/LTS_2_B3LYP_Ahl_b.f98').cry_read_density()
'''
density1 = '../examples/data/density/substrate.f98'
density2 = '../examples/data/density/adsorbate.f98'
density3 = '../examples/data/density/system.f98'
cry_combine_density(density1,density2,density3,'../examples/data/density/new_density.f98')'''

    