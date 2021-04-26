# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:32:29 2021

@author: eliam
"""

def iterazione(gname, x0, tolx, nmax):
    """ Metodo generico per la ricerca di zeri di funzione.
        
    Riscriviamo f(x) = 0 nella forma x = g(x) per una g opportuna, così che la 
    ricerca degli zeri possa essere ricondotta allo studio di g, ossia
    all'approssimazione di alpha t.c.
    f(alpha) = 0 se e solo se g(alpha) = alpha
    
    Parameters
    ----------
    gname : lambdified function
        Funzione g(x) t.c. 
        si ha la funzione f(x) = 0 riscritta nella forma x = g(x) 
    x0 : number
        Iterata iniziale
    tolx : number
        tolleranza per il test d'arresto sull'incremento
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
    xk.append(x0)
    x1 = gname(x0)
    xk.append(x1)
    it = 1
    while it < nmax and abs(x1 - x0) >= tolx * abs(x1):
        x0 = x1
        x1 = gname(x0)
        xk.append(x1)
        it += 1
    
    if it >= nmax:
        print("Numero massimo di iterazioni possibile raggiunto.")
    
    return x1, it, xk