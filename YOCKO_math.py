#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 18:51:16 2021

@author: comecattin

All the math needed in YOCKO

"""
import numpy as np
import math
import YOCKO_tools

class gaussian:
    """
    Definition of gaussian function as a class
    """
    def __init__(self,alpha,mean):
        self.alpha = alpha
        self.center = mean

    def product(self,gauss_B):
        """
        Product between 2 gaussians
        """
        a , center_A = self.alpha , self.center
        b , center_B = gauss_B.alpha , gauss_B.center
        
        p = a+b
        
        diff=np.linalg.norm( center_A - center_B )**2
        Norm= (4*a*b/(np.pi**2))**0.75
        New_prefac= Norm*np.exp(-a*b/p*diff)
        New_center=(a*center_A+b*center_B)/p

        return p, diff,New_prefac,New_center
    
    def overlap(self,gauss_B):
        """
        Calculate the overlap of two gaussian function gauss_A and gauss_B
        """
        p , diff , New_prefac , New_center = self.product(gauss_B)

        A = (np.pi/p) ** 1.5
        val = A * New_prefac
        
        return val
    
    def kinetic(self,gauss_B):
        """
        Calculate the kinetic integral 
        """
        p , diff , New_prefac , New_center = self.product(gauss_B)
        
        A = (np.pi/p) ** 1.5
        r_exp = self.alpha * gauss_B.alpha / p
        val = r_exp * (3 - 2*r_exp*diff) * A * New_prefac
        
        return val

def F_0(alpha):
    """
    Function needed to calculate the potential and the electron electronrepulsion integrals
    """
    if alpha == 0:
        val = 1
    else :
        val = 0.5 * (np.pi/alpha)**0.5 * math.erf(alpha**0.5)
        
    return val


def potential(gauss_A,gauss_B,index_atom,file_name):
    """
    Return the nuclear/electron potential of the atom indexed by index_atom in the file 'file_name.xyz'
    """
    p , diff , New_prefac , New_center = gauss_A.product(gauss_B)
    N_atoms, atoms , coord = YOCKO_tools.read_xyz(file_name)
    dico_charge = YOCKO_tools.charge()
    
    #Position of the nuclei
    R = coord[index_atom]
    #Charge of the nuclei
    Z = dico_charge[atoms[index_atom]]
    
    alpha = p * np.linalg.norm(New_center - R)**2
    
    val = (-2*np.pi*Z/p)*New_prefac*F_0(alpha)
    
    return val

def multi_electron(gauss_A,gauss_B,gauss_C,gauss_D):
    """
    Return the multi-electron integral
    """
    p_ab , diff_ab , New_prefac_ab , New_center_ab = gauss_A.product(gauss_B)
    p_cd , diff_cd , New_prefac_cd , New_center_cd = gauss_C.product(gauss_D)
    alpha = p_ab * p_cd / (p_ab + p_cd) * np.linalg.norm(New_center_ab-New_center_cd)**2
    prefactor = 2 * np.pi**2.5 * (p_ab*p_cd*(p_ab+p_cd)**0.5)** -1
    
    val = prefactor * New_prefac_ab * New_prefac_cd * F_0(alpha)

    return val






if __name__ == '__main__':
    file_name = 'H2.xyz'
    data , atoms_basis, exp_coeff_S, exp_coeff_P = YOCKO_tools.read_basis('sto-3g.1.gbs')
    H_1S_1 = gaussian(exp_coeff_S[0,0,0],exp_coeff_S[0,0,1])
    H_1S_2 = gaussian(exp_coeff_S[0,1,0],exp_coeff_S[0,1,1])
    H_1S_3 = gaussian(exp_coeff_S[0,2,0],exp_coeff_S[0,2,1])
    H_1S = [H_1S_1,H_1S_2,H_1S_3]
    prod_test = H_1S_1.product(H_1S_2)
    overlap_test = H_1S_2.overlap(H_1S_2)
    kinetic_test = H_1S_1.kinetic(H_1S_1)
    potential_test = potential(H_1S_1,H_1S_2,1,file_name)
    multi_electron_test = multi_electron(H_1S_1,H_1S_2,H_1S_3,H_1S_1)

