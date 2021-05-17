# -*- coding: utf-8 -*-
"""
Created on Mon May 17 15:24:02 2021

@author: eliam
"""


import numpy as np
import scipy.linalg as spl
from sistemilineari import Usolve

def metodoQR(x, y, n):
    """ metodoQR(x, y, n)
    
    Desc
        
    Parameters
    ----------
    x : vector
        Vettore colonna con le ascisse dei punti
    y : vector
        Vettore colonna con le ordinate dei punti
    n : number
        Grado del polinomio
        
    Returns
    -------
    vector
        a : Vettore colonna contenente i coefficienti incogniti
    
    """
    
    H = np.vander(x, n+1)
    Q, R = spl.qr(H)
    y1 = np.dot(Q.T, y)
    a, flag = Usolve(R[0:n+1, :], y1[0:n+1])
    return a
