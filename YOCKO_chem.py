#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:03:16 2021

@author: comecattin

The chemistry part of YOCKO
"""

import numpy as np
import YOCKO_tools
import YOCKO_math

def get_integrals(file_name):
    """
    Compute the integrals of the system
    """
    
    size = YOCKO_tools.basis_size(file_name)
    N_atoms, atoms_list , coord = YOCKO_tools.read_xyz(file_name)
    data , atoms_basis, exp_coeff_S, exp_coeff_SP = YOCKO_tools.read_basis('sto-3g.1.gbs')
    
    #Initialisation
    S = np.zeros((size,size)) #Overlap
    T = np.zeros((size,size)) #Kinetic
    V = np.zeros((size,size)) #Potential
    multi_elec = np.zeros((size,size,size,size)) #Multi electron tensor
    
    #Atoms
    for i , atoms in enumerate(atoms_list):
        
        R = coord[i]
        Z = YOCKO_tools.charge()[atoms]
        
        #Quantum numbers
        for j in range(YOCKO_tools.quantum_number()[atoms]):
            
            #Pas terminé....
            #get alpha
            alpha_list = 1
            
            print(j)
        
        
    

    


if __name__ == '__main__':
    file_name = 'H2.xyz'
    get_integrals(file_name)