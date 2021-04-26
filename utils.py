# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:48:27 2021

@author: eliam
"""
import math
import numpy as np

def stima_ordine(xk,iterazioni):
    """ Stima numerica dell'ordine di convergenza
    
    Parameters
    ----------
    xk : list
        Lista delle iterate
    iterazioni: int
        Numero di iterazioni
    
    Returns
    -------
    int
        Ordine di convergenza
    """
    p=[]

    for k in range(iterazioni-3):
        p.append(np.log(abs(xk[k+2]-xk[k+3])/abs(xk[k+1]-xk[k+2]))/np.log(abs(xk[k+1]-xk[k+2])/abs(xk[k]-xk[k+1])));
    ordine=p[-1]
    return ordine
  

# =============================================================================
#     Il core Python non possiede la funzione sign.
#     La funzione copysign(a,b)  del modulo math restituisce un valore numerico che ha il valore assoluto di
#     a e segno di b.
#     Per avere il segno di un valore numerico b si può usare math.copysign(1,b)
#     che resistuisce 1 se b>0, -1 se b<0, 0 se b è zero
# =============================================================================
def sign(x):
    """Funzione sign(x)     
    
    Parameters
    ----------
    x : number
        Il numero di cui si vuol sapere il segno
    
    Returns
    -------
    int
        1 se b>0, -1 se b<0, 0 se b=0
    """
    return math.copysign(1, x)
