#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 19:49:16 2021

@author: comecattin

lecture de fichier .xyz



"""

import numpy as np

def read_xyz(file_name):
    """
    Load .xyz file
    """
    
    file = open(file_name,'r')
    data = np.array([0,0,0,0])
    
    for i, line in enumerate(file):    
        
        if i == 0 :
            N_atoms = np.int(line)
        
        elif line.split() != [] :
            data = np.vstack((data,np.array(line.split())))
            
    file.close()
    
    data = data[1:]
    atoms = data[:,0]
    coord = data[:,1:].astype(np.float)
    return N_atoms , atoms , coord

def charge():
    """
    Return dictionnary charge of each ion
    """
    dico_charge = {'H' : 1,'He' : 2, 'Li' : 3 ,'Be' : 4, 'B' : 5, 'C' : 6,
                   'N' : 7, 'O' : 8, 'F' : 9, 'Ne' : 10}
    
    return dico_charge

def read_basis(file_name):
    """
    read basis gausian file
    only with STO-3G in this version
    """
    file = open(file_name,'r')
    
    data = []
    atoms = []
    exp_coeff_S = []
    exp_coeff_SP = []
    
    for i, line in enumerate(file):
        
        data.append(line.replace('*','').replace('D','e').split())
        
        if line[0:2] == 'S ':
            
            atoms.append(data[i-1][0])
            
        if line[0] == '*':
            if data[i-4][0] == 'S':
                exp_coeff_S.append(data[i-3:i])
            
            if data[i-4][0] == 'SP':
                exp_coeff_SP.append(data[i-3:i])
                exp_coeff_S.append(data[i-7:i-4])
                                
    file.close()
        
    data = data[12:]
    
    return data , atoms, np.array(exp_coeff_S).astype(np.float), np.array(exp_coeff_SP).astype(np.float)
    










if __name__ == '__main__':
    file_name = 'luminol.xyz'
    N_atoms, atoms , coord = read_xyz(file_name)
    dico_charge = charge()
    data , atoms, exp_coeff_S, exp_coeff_SP = read_basis('sto-3g.1.gbs')