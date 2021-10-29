#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 19:49:16 2021

@author: comecattin

read needed file for YOCKO



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

def quantum_number():
    """
    Return dictionnary max quantum number of each atom
    """
    dico_quantum_number = {'H' : 1,'He' : 1, 'Li' : 2 ,'Be' : 2, 'B' : 2, 'C' : 2,
                   'N' : 2, 'O' : 2, 'F' : 2, 'Ne' : 2}
    
    return dico_quantum_number

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
    
    exp_coeff_S = np.array(exp_coeff_S).astype(np.float)
    exp_coeff_SP = np.array(exp_coeff_SP).astype(np.float)
    
    exp_coeff_S = np.vstack((exp_coeff_S,exp_coeff_SP[:,:,0::2]))
    exp_coeff_P = exp_coeff_SP[:,:,0:2]
    
    
    
    
    return data , atoms, exp_coeff_S, exp_coeff_P
    

def basis_size(file_name_xyz):
    """
    Return the size of the basis considered
    """
    atoms = read_xyz(file_name_xyz)[1]
    size = 0
    dico_quantum_number = quantum_number()
    
    for i in atoms:
        size = size + dico_quantum_number[i]
        
    return size

# def pro_gauss(gauss_A, gauss_B):
# 	"""
# 	Product of two gaussian gauss_A and gauss_B
# 	"""
# 	a , center_A = gauss_A.alpha , gauss_A.center
# 	b , center_B = gauss_B.alpha , gauss_B.center
# 	p= a+b
# 	diff=np.linalg.norm( center_A - center_B )**2
# 	Norm= (4*a*b/(np.pi**2))*0.75
# 	New_prefac= Norm*np.exp(-a*b/p*diff)
# 	New_center=(a*center_A+b*center_B)/p
# 	
# 	   
# 	   #Input needs to be the parameters (tuple type) of each gaussian function
# 	   #(a, center_A) where (1/(sig*sqrt(2pi)))*exp(-(x-mu)**2/(2sig**2)) is gauss_A 
# 	   #as a function of x
# 	   #a=1/(2sig**2)
# 	   #mu=center_A
# 	   #Na prefactor factor=(1/(sig*sqrt(2pi))))
# 	#Parameters of the new gaussian
# 	
# 		
# 	return p, diff,New_prefac,New_center
# 	
# def Gauss_overlap(gauss_A,gauss_B):
# 	"""
# 	Calculate the overlap of two gaussian function gauss_A and gauss_B
# 	"""
# 	p , diff , New_prefac , New_center = pro_gauss(gauss_A,gauss_B)
# 	
# 	A = (np.pi/p) ** 1.5
# 	val = A * New_prefac
# 	
# 	return val

        
    
if __name__ == '__main__':
    
    file_name_xyz = 'H2.xyz'
    file_name_basis = 'sto-3g.1.gbs'
    
    N_atoms, atoms , coord = read_xyz(file_name_xyz)
    data , atoms_basis, exp_coeff_S, exp_coeff_P = read_basis(file_name_basis)
    
    dico_charge = charge()
    dico_quantum_number = quantum_number()
    
    basis_size = basis_size(file_name_xyz)
    
	