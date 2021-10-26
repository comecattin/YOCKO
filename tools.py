#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 19:49:16 2021

@author: comecattin

lecture de fichier .xyz



"""

import numpy as np
from scipy import random

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




def pro_gauss(gauss_A, gauss_B):
	a,center_A=gauss_A
	b, center_B= gauss_B
	p= a+b
	diff=np.linalg.norm( center_A - center_B )**2
	Norm= (4*a*b/(math.pi**2))*0.75
	New_prefac= N*exp(-a*b/p*diff)
	New_center=(a*center_A+b*center_B)/p
	
	   
	   #Input needs to be the parameters (tuple type) of each gaussian function
	   #(a, center_A) where (1/(sig*sqrt(2pi)))*exp(-(x-mu)**2/(2sig**2)) is gauss_A 
	   #as a function of x
	   #a=1/(2sig**2)
	   #mu=center_A
	   #Na prefactor factor=(1/(sig*sqrt(2pi))))
	#Parameters of the new gaussian
	
		
	return (p, diff,New_prefac,New_center)
	
def Gauss_overlap(A,B):
        
    
if __name__ == '__main__':
    file_name = 'luminol.xyz'
    N_atoms, atoms , coord = read_xyz(file_name)
    dico_charge = charge()
	