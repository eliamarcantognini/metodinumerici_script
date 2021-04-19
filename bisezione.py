# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:12:32 2021

@author: eliam
"""

import numpy as np
import math
from utils import sign

def bisez(fname, a, b, tol):
    """ Metodo di bisezione per la ricerca di zeri di funzione.
    
    Descrizione della funzione
    
    Parameters
    ----------
    fname : lambdified function
        Funzione principale
    a : number
        Limite sinistro dell'intervallo
    b : number
        Limite destro dell'intervallo
    tol : number
        tolleranza per il test d'arresto sull'incremento
        
    Returns
    -------
    list
        x : soluzione
        it : numero di iterazioni
        xk : iterabili che convergono allo zero della funzione
    
    """

    eps = np.spacing(1) # np.spacing(x) Restituisce la distanza tra x e il numero adiacente più vicino.
                        # np.spacing(1)  restituisce quindi l' eps di macchina.
    fa = fname(a)
    fb = fname(b)
    if sign(fa) == sign(fb):
        print("Intervallo non corretto. Algoritmo non applicabile.")
        return [], 0, []
    else:
        maxit = math.log((fb - fa)/tol, 10)/math.log(2, 10)
        print("Numero massimo di iterazioni: ", maxit)
        it = 0
        xk = []
        while it < maxit and abs(b - a) >= tol + eps*max([abs(a), abs(b)]):
            m = a + (b - a)/2 # punto medio fra a e b
            xk.append(m) # salvo in xk i vari valori di m
            fm = fname(m) # valore della funzione nell'ascissa m 
            if sign(fm) == 0: 
                # se fm è 0, sono sull'ascissa quindi ho trovato x
                break
            elif sign(fm) == sign(fa): 
                # se fm ha lo stesso segno di fa vuol dire che devo prendere
                # m come nuovo limite sinistro dell'intervallo
                a = m # il punto medio diventa il nuovo limite sinistro
                fa = fm
            elif sign(fm) == sign(fb):
                # se fm ha lo stesso segno di fb vuol dire che devo prendere
                # m come nuovo limite destro dell'intervallo
                b = m # il punto medio diventa il nuovo limite destro
                fb = fm
            
            it += 1
        x = m 
        return x, it, xk
                
                
                
            
            
        