#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:03:16 2021

@author: comecattin

The chemistry part of YOCKO
"""

import numpy as np
import matplotlib.pyplot as plt
import YOCKO_tools
import YOCKO_math

def get_integrals(file_name_xyz,file_name_basis):
    """
    Compute the integrals of the system
    """
    
    size = YOCKO_tools.basis_size(file_name_xyz)
    N_atoms, atoms_list , coord = YOCKO_tools.read_xyz(file_name_xyz)
    data , atoms_basis, exp_coeff_S, exp_coeff_P = YOCKO_tools.read_basis(file_name_basis)
    #Initialisation
    S = np.zeros((size,size)) #Overlap
    T = np.zeros((size,size)) #Kinetic
    V = np.zeros((size,size)) #Potential
    multi_elec = np.zeros((size,size,size,size)) #Multi electron tensor
  



    #Atoms
    for i , atoms_A in enumerate(atoms_list):
        
        R_A = coord[i]
        Z_A = YOCKO_tools.charge()[atoms_A]
        
        #Orbital
        for j in range(YOCKO_tools.quantum_number()[atoms_A]):
            
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
                
                
                
                
                #Other atom
                for i_prime , atoms_B in enumerate(atoms_list):
                    R_B = coord[i_prime]
                    Z_B = YOCKO_tools.charge()[atoms_B]
                    
                    #Other Orbital
                    for j_prime in range(YOCKO_tools.quantum_number()[atoms_B]):
                        
                        if j_prime == 0:
                            #S type
                            #get alpha
                            alpha_list_B = exp_coeff_S[i_prime,:,0]
                            #get coefficient
                            coeff_list_B = exp_coeff_S[i_prime,:,1]
                            
                            
                        if j_prime == 1:
                            #P type
                            alpha_list_B = exp_coeff_P[i_prime,:,0]
                            #get coefficient
                            coeff_list_B = exp_coeff_P[i_prime,:,1]
                            
                        #Other Gaussian
                        for k_prime in range(len(alpha_list_B)):
                            gauss_B = YOCKO_math.gaussian(alpha_list_B[k_prime],
                                                          coeff_list_B[k_prime])
                            
                            a = (i + 1) * (j+1) - 1
                            b = (i_prime + 1) * (j_prime +1) - 1
                            
                            #Overlap & kinetic matrix
                            S[a,b] += gauss_A.center * gauss_B.center * gauss_A.overlap(gauss_B)
                            T[a,b] += gauss_A.center * gauss_B.center * gauss_A.kinetic(gauss_B)
                            for l in range(N_atoms):
                                V += gauss_A.center * gauss_B.center * YOCKO_math.potential(gauss_A, 
                                                                                            gauss_B, 
                                                                                            l, 
                                                                                            file_name_xyz)
                                
                            
                            #Multi-electron tensor
                            #Atoms
                            for i_pp , atoms_C in enumerate(atoms_list):
                                R_C = coord[i_pp]
                                Z_C = YOCKO_tools.charge()[atoms_C]
                                #Orbital
                                for j_pp in range(YOCKO_tools.quantum_number()[atoms_C]):
                                    
                                    if j_pp == 0:
                                        #S type
                                        #get alpha
                                        alpha_list_C = exp_coeff_S[i_pp,:,0]
                                        #get coefficient
                                        coeff_list_C = exp_coeff_S[i_pp,:,1]
                                        
                                        
                                    if j_pp == 1:
                                        #P type
                                        alpha_list_C = exp_coeff_P[i_pp,:,0]
                                        #get coefficient
                                        coeff_list_C = exp_coeff_P[i_pp,:,1]
                                        
                                    #Gaussian
                                    for k_pp in range(len(alpha_list_C)):
                                        gauss_C = YOCKO_math.gaussian(alpha_list_C[k_pp],
                                                                      coeff_list_C[k_pp])
                                        
                                        
                                        #Other Atoms
                                        for i_ppp , atoms_D in enumerate(atoms_list):
                                            R_D = coord[i_ppp]
                                            Z_D = YOCKO_tools.charge()[atoms_D]
                                            #Orbital
                                            for j_ppp in range(YOCKO_tools.quantum_number()[atoms_D]):
                                                
                                                if j_ppp == 0:
                                                    #S type
                                                    #get alpha
                                                    alpha_list_D = exp_coeff_S[i_ppp,:,0]
                                                    #get coefficient
                                                    coeff_list_D = exp_coeff_S[i_ppp,:,1]
                                                    
                                                    
                                                if j_ppp == 1:
                                                    #P type
                                                    alpha_list_D = exp_coeff_P[i_ppp,:,0]
                                                    #get coefficient
                                                    coeff_list_D = exp_coeff_P[i_ppp,:,1]
                                                    
                                                #Gaussian
                                                for k_ppp in range(len(alpha_list_D)):
                                                    gauss_D = YOCKO_math.gaussian(alpha_list_D[k_ppp],
                                                                                  coeff_list_C[k_ppp])
                                                    
                                                    c = (i_pp + 1) * (j_pp + 1) - 1
                                                    d = (i_ppp + 1) * (j_ppp +1) - 1
                                                    
                                                    multi_elec[a,b,c,d] = gauss_A.center * gauss_B.center * gauss_C.center * gauss_D.center * YOCKO_math.multi_electron(gauss_A,
                                                                                                                                                                        gauss_B,
                                                                                                                                                                        gauss_C,
                                                                                                                                                                        gauss_D)
                                                    
                                
    H_core = T + V
                                
    return S , T , V , multi_elec, H_core
                
                

def orthogonalisation(S):
    """
    Symetric orthogonalisation of the basis
    """
    evalS , U = np.linalg.eig(S)
    diagS = np.dot(U.T,np.dot(S,U))
    diagS_sqrt = np.diag(np.diagonal(diagS)**-0.5)
    X = np.dot(U,np.dot(diagS_sqrt,U.T))
    
    return X

                
                
                
#The difference bewteen tthe two most recent successive density matrices 
def Diff_Succ_DM(P_1,P_2,file_name):
     x=0
     size = YOCKO_tools.basis_size(file_name)
     for i in range(size):
         for j in range(size):
             x += size**-2*(P_1[i,j]-P_2[i,j])**2
     return x**0.5
 
 
 
def Algo(N,file_name_xyz,file_name_basis):
    """
    SCF Convergence algorithm
    """
    
    S , T , V , multi_elec, H_core = get_integrals(file_name_xyz, file_name_basis)
    size = YOCKO_tools.basis_size(file_name_xyz)
    
    P=np.zeros((size,size))	
    P_prev=np.zeros((size,size))
    P_list=[]
    
    threshold=100

    while threshold>10**-4:
        G=np.zeros((size,size))
        for i in range(size):
            for j in range(size):
                for x in range(size):
                    for y in range(size):
                        G[i,j]+=P[x,y]*(multi_elec[i,j,y,x]-0.5*multi_elec[i,x,y,j])
        X = orthogonalisation(S)
        Fock = H_core + G	
        Fock_orth = np.dot(X.T,np.dot(Fock,X))
        evalFock_orth, C_orth= np.linalg.eig(Fock_orth)	

        index=evalFock_orth.argsort()
        evalFock_orth=evalFock_orth[index]
        C_orth=C_orth[:,index]
        C = np.dot(X,C_orth)	

        #Form new density matrix P (sum over electron pairs, not over the entire basis set)
        for i in range (size):
            for j in range (size):
                for a in range(int(N/2)):
                    P[i,j]=2*C[i,a]*C[j,a]
        P_list.append(P)
        threshold=Diff_Succ_DM(P_prev, P, file_name_xyz)
        P_prev=P.copy()
    print('\n')
    print("YOCKO took {} step to converge".format(len(P_list)))
    print('\n')
    print("The orbital energies are {} and {} Ha".format(evalFock_orth[0],evalFock_orth[1]))
    print('\n')
    print(f'The orbital matrix is: \n\n{C}')
    print('\n')
    print(f'The density/bound order matrix is \n\n{P}')
    return (P_prev, evalFock_orth)

def geometry_opt(file_name_xyz,file_name_basis,N):
    """
    Optimise the geometry of the molecule
    """
    
    end = range(100)
    x = []
    y = []
    
    for i in end :
        #read
        N_atoms , atoms , coord = YOCKO_tools.read_xyz(file_name_xyz)
        x_1 , x_2 = coord[:,0]
        y_1 , y_2 = coord[:,1]
        z_1 , z_2 = coord[:,2]
        
        #Distance
        d = np.sqrt((x_1 - x_2)**2 + (y_1 - y_2)**2 + (z_1 - z_2)**2)
        x.append(d)
        
        #Result
        P_prev, evalFock_orth = Algo(N,file_name_xyz, file_name_basis)
        y.append(evalFock_orth[0])
        
        #Update
        coord[:,2] = coord[:,0] - 1/20
        
        YOCKO_tools.change_xyz(file_name_xyz,coord)
        file_name_xyz = 'new.xyz'
        
        
        print('SCF Done')
    return x,y
        

        
def geometry_opt_plot(x,y):
    """
    Plot the energy versus the inter-atomic distance 
    """
    
    fig = plt.figure()
    plt.plot(x,y,'*')
    plt.show()
    
    
    




if __name__ == '__main__' :
    
    file_name_xyz = 'H2.xyz'
    file_name_basis = 'sto-3g.1.gbs'
    N = 2
    S , T , V, multi_elec, H_core = get_integrals(file_name_xyz,file_name_basis)
    X = orthogonalisation(S)
    #Algo(N,file_name_xyz,file_name_basis)
    x , y = geometry_opt(file_name_xyz,file_name_basis,N)
    geometry_opt_plot(x,y)
	
	
	