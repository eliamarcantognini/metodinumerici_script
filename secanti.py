# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:43:26 2021

@author: eliam
"""

import numpy as np
import math
from utils import sign

def secanti(fname, xm1, x0, tolx, tolf, nmax):
    """ Metodo delle secanti per la ricerca di zeri di funzione.
    
    Descrizione della funzione
    
    Parameters
    ----------
    fname : lambdified function
        Funzione principale
    xm1 : number
        Iterata (x-1)
    x0 : number
        Iterata iniziale
    a : number
        Limite sinistro dell'intervallo
    b : number
        Limite destro dell'intervallo
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
    fxm1 = fname(xm1)
    # il valore di m nel metodo delle secanti è un approssimazione della derivata
    # prima di f come rapporto incrementale tra le ultime due valutazioni di
    # funzione f(xm1) e f(x0)
    m = (fx0 - fxm1) / (x0 - xm1)
    # il calcolo di x1 è lo stesso del metodo delle corde e di newton, a cambiare è m
    # geometricamente, ad ogni iterazione si approssima il grafico della funzione f
    # con la retta che passa per i punti (xm1, f(xm1)) e (x0, f(x0))
    x1 = x0 - fx0 / m # calcolo xi +1
    fx1 = fname(x1) # calcolo il valore di f in xi+1
    xk.append(x1)
    it = 1 # un'iterazione l'ho già fatta, calcolando xi+1
    while it < nmax and abs(fx1) >= tolf and abs(fx0 / m)>=tolx*abs(x1):
        # procedo con l'iterazione successiva, dando come valore precedente x1,
        # come valore precedente al precedente(xm1) x0 e calcolando il successivo xi+1
        xm1 = x0
        x0 = x1
        fx0 = fname(x0)
        fxm1 = fname(xm1)
        # calcolo m_i 
        m = (fx0 - fxm1) / (x0 - xm1)    
        fx0 = fx1 
        x1 = x0 - fx0 / m
        fx1 = fname(fx1)
        xk.append(x1)
        it += 1
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk