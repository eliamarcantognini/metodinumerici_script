# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:43:19 2021

@author: eliam
"""

import numpy as np
import math
from utils import sign

def newton(fname, fpname, x0, tolx, tolf, nmax):
    """ Metodo di Newton per la ricerca di zeri di funzione.
    
    Descrizione della funzione
    
    Parameters
    ----------
    fname : lambdified function
        Funzione principale
    fpname : lambdified function
        Derivata prima della funzione principale
    x0 : number
        Iterata iniziale
    tolx : number
        tolleranza per il test d'arresto sull'incremento
    tolf: number
        tolleranza per il test del residuo
    nmax: int
        numero massimo di iterazioni possibili prima di fermare l'algoritmo
        
    Returns
    -------
    list
        x : soluzione
        it : numero di iterazioni
        xk : iterabili che convergono allo zero della funzione
    
    """
    
    xk = []
    fx0 = fname(x0)
    # il valore di m nel metodo di newton è la derivata della funzione in x0
    m = fpname(x0)
    # se la derivata non è nulla proseguo con l'algoritmo
    if abs(m) > np.spacing(1):
        # geometricamente si prende come nuova approssimazione l'intersezione dell'asse
        # delle ascisse con la retta tangente a f in (x0, f(x0))
        x1 = x0 - fx0 / m # calcolo xi +1
        fx1 = fname(x1) # calcolo il valore di f in xi+1
        xk.append(x1)
        it = 1 # un'iterazione l'ho già fatta, calcolando xi+1
    else:
        print("Derivata nulla in x0.")
        return [], 0, xk
        
    while it < nmax and abs(fx1) >= tolf and abs(fx0 / m)>=tolx*abs(x1):
        # procedo con l'iterazione successiva, dando come valore precedente x1,
        # come valore precedente al precedente(xm1) x0 e calcolando il successivo xi+1
        x0 = x1
        fx0 = fname(x0)
        # calcolo m_i 
        m = fpname(x0)
        if abs(m) > np.spacing(1): 
            fx0 = fx1 
            x1 = x0 - fx0 / m
            fx1 = fname(fx1)
            xk.append(x1)
            it += 1
        else:
            print("Derivata nulla in x0.")
            return x1, it, xk
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk