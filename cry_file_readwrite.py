#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:28:28 2021

@authors: brunocamino
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
        
        end_index = [i for i, s in enumerate(data) if 'END' in s]
        
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
            
            

       
'''TESTING
mgo = crystal_input('data/mgo.d12')     
#print(mgo.name)
print(mgo.scf_block)
#mgo.scf_block.remove('DIIS\n')
#print(mgo.scf_block)'''