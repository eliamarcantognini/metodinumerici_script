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
    
    eps = np.spacing(1) # np.spacing(x) Restituisce la distanza tra x e il numero adiacente più vicino.
                        # np.spacing(1)  restituisce quindi l' eps di macchina.
    
    