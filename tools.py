#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 19:49:16 2021

@author: comecattin

lecture de fichier .xyz



"""

import numpy as np
from scipy import random


class gaussian :
	"""
	Definition of gaussian function as a class
	"""
	def __init__(self,alpha,mean) :
		self.alpha = alpha
		self.center = mean
		


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
    

def pro_gauss(gauss_A, gauss_B):
	"""
	Product of two gaussian gauss_A and gauss_B
	"""
	a , center_A = gauss_A.alpha , gauss_A.center
	b , center_B = gauss_B.alpha , gauss_B.center
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
	
		
	return p, diff,New_prefac,New_center
	
def Gauss_overlap(gauss_A,gauss_B):
	"""
	Calculate the overlap of two gaussian function gauss_A and gauss_B
	"""
	p , diff , New_prefac , New_center = pro_gauss(gauss_A,gauss_B)
	
	A = (np.pi/p) ** 1.5
	val = A * New_prefac
	
	return val

        
    
if __name__ == '__main__':
    file_name = 'luminol.xyz'
    N_atoms, atoms , coord = read_xyz(file_name)
    dico_charge = charge()
    data , atoms, exp_coeff_S, exp_coeff_SP = read_basis('sto-3g.1.gbs')
	