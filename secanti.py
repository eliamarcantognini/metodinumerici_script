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
    
    return 1