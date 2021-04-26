# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:29:16 2021

@author: eliam
"""

import numpy as np
from utils import sign

def regula_falsi(fname, a, b, tol, nmax):
    """ Metodo di falsa posizione per la ricerca di zeri di funzione.
    
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
    eps = np.spacing(1) # np.spacing(x) Restituisce la distanza tra x e il numero adiacente più vicino.
                        # np.spacing(1)  restituisce quindi l' eps di macchina.  
    fa = fname(a)
    fb = fname(b)
    it = 0
    xk = []
    if sign(fa) == sign(fb):
        print("Intervallo non corretto. Algoritmo non applicabile.")
        return [], it, xk
    else:
        fx1 = fname(a) # necessario per usare la condizione d'arresto generica su fx1
        while (
                it < nmax and
                abs(b-a) >= tol + eps*max(abs(a), abs(b)) and
                abs(fx1) >= tol
            ):
            x1 = a - fa * (b-a) / (fb - fa) # trovo x1 che è l'ascissa del punto (x1, fname(x1))
            xk.append(x1)
            fx1 = fname(x1) # il valore della funzione in x1
            if fx1 == 0: 
                # se fx1 è 0, ho trovato x (alpha)
                break
            elif sign(fx1) == sign(fa):
                # se fx1 ha lo stesso di fa significa che devo prendere x1
                # come nuovo limite sinistro dell'intervallo
                a = x1
                fa = fx1
            elif sign(fx1) == sign(fb):
                # se fx1 ha lo stesso di fb significa che devo prendere x1
                # come nuovo limite destro dell'intervallo
                b = x1
                fb = fx1
            it += 1 # proseguo con l'iterazione successiva
        
        if it >= nmax:
            # se ho raggiunto il numero massimo di iterazioni possibili
            # stampo l'informazione così che l'utente sappia che non
            # è stata raggiunta la convergenza entro nmax iterazioni
            print("Raggiunto numero massimo di iterazioni possibili.")
        
        return x1, it, xk
