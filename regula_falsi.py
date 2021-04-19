# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:29:16 2021

@author: eliam
"""

import numpy as np
import math
from utils import sign

def regula_falsi(fname, a, b, tol, nmax):
    """ Metodo di falsa posizione per la ricerca di zeri di funzione.
    
    Descrizione della funzione
    
    Parameters
    ----------
    fname : lambdified function
        Funzione principale
    a : number
        Limite sinistro dell'intervallo
    b : number
        Limite destro dell'intervallo
    tol: number
        tolleranza per il test del residuo
    nmax: int
        numero massimo di iterazioni possibili prima di fermare l'algoritmo
        
    Returns
    -------
    list
        x1 : soluzione
        it : numero di iterazioni
        xk : iterabili che convergono allo zero della funzione
    
    """

    return 1    
