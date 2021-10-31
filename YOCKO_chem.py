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
    data , atoms_basis, exp_coeff_S, exp_coeff_P = YOCKO_tools.read_basis('sto-3g.1.gbs')
    
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
            
            if j == 0:
                #S type
                #get alpha
                alpha_list = exp_coeff_S[i,:,0]
                #get coefficient
                coeff_list = exp_coeff_S[i,:,1]
                
            if j == 1:
                #P type
                alpha_list = exp_coeff_P[i,:,0]
                #get coefficient
                coeff_list = exp_coeff_P[i,:,1]
                
            #Gaussian
            for k in range(len(alpha_list)):
                gauss_A = YOCKO_math.gaussian(alpha_list[k],coeff_list[k])
                print(gauss_A.alpha , gauss_A.center)
                
#B is the basis set size, and the dimention of the S T V martices
#The difference bewteen tthe two most recent successive density matrices 
def Diff_Succ_DM(P_1,P_2):
	x=0
	size = YOCKO_tools.basis_size(file_name)
	for i in range(size):
		for j in range(size):
			x+=B**-2*(P_1[i,j]-P_2[i,j])**2
	return x**0.5

#The initial guess for P will be the zero matrice



                
                

if __name__ == '__main__':
	
    file_name = 'H2.xyz'
    get_integrals(file_name)
size = YOCKO_tools.basis_size(file_name)
	
P=np.zeros((size,size))
P_prev=np.zeros((size,size))
P_list=[]
Hcore= T+V	
threshold=100

while threshold>10**-4:
	G=np.zeros((size,size))
	for i in range(size):
		for j in range(size):
			for x in range(size):
				for y in range(size):
				   G[i,j]=+=P[x,y]*(multi_elec_tensor[i,j,y,x]-0.5*multi_elec_tensor[i,x,y,j])
	Fock = Hcore + G
	
	
Fock_orth=dot(X.T,dor(Fock,X))	
evalFock_orth,C_orth=np.linalg.eig(Fock_orth)

index=evalFock_orth.argsort()
evalFock_orth=evalFock_orth[index]
C_orth=C_orth[:,index]
C=dot(X,C_orth)
N= #number of electrons, je ne trouve pas dans le code le nom de cette variable x)	
	
#Form new density matrix P (sum over electron pairs, not over the entire basis set)
for i in range (size):
	for j in range (size):
		for a in range(int(N/2)):
			P[i,j]=2*C[i,a]*C[j,a]
P_list.append(P)
treshold=Diff_Succ_DM(P_prev, P)
P_prev=P.copy()
	
	
	