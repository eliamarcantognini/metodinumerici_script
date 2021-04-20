# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:35:48 2021

@author: eliam
"""

import numpy as np
import math
from utils import sign

def corde(fname, fpname, x0, tolx, tolf, nmax):
    """ Metodo delle corde per la ricerca di zeri di funzione.
        
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
    # il valore di m nel metodo delle corde è costante ed è convenzione prendere
    # come valore costante la valutazione in x0 della derivata prima di f
    # che è il coefficiente angolare della tangente in x0
    m = fpname(x0) 
    # x1 è l'ascissa del punto di intersezione tra la retta che passa per il punto
    # (xi, f(xi)) con ha pendenza uguale a m e l'asse x
    x1 = x0 - fx0 / m # calcolo xi +1
    fx1 = fname(x1) # calcolo il valore di f in xi+1
    it = 1 # un'iterazione l'ho già fatta, calcolando xi+1
    xk.append(x1)
    while it < nmax and abs(fx1) >= tolf and abs(fx0 / m)>=tolx*abs(x1):
        # procedo con l'iterazione successiva, dando come valore precedente x1
        # e calcolando il successivo xi+1
        fx0 = fx1 
        x1 = x0 - fx0 / m
        fx1 = fname(fx1)
        xk.append(x1)
        it += 1
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk
            