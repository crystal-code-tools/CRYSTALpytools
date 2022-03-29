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
