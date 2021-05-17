# -*- coding: utf-8 -*-
"""
Created on Mon May 17 17:27:13 2021

@author: eliam
"""
import numpy as np
# =============================================================================
# 
# =============================================================================
def plagr(xnodi, k):
    """ plagr(xnodi, k)
    
    Restituisce i coefficienti del k-esimo polinomio di Lagrange associato
    ai punti del vettore xnodi
        
    Parameters
    ----------
    xnodi : vector
        Vettore con i nodi
    k : number
        numero del nodo
        
    Returns
    -------
    list
        p : coefficienti del polinomio di Lagrange
    
    """
    
    xzeri = np.zeros_like(xnodi)
    n = xnodi.size
    if k == 0:
        xzeri = xnodi[1:n]
    else:
        # elimo il k-esimo nodo
        xzeri = np.append(xnodi[0:k], xnodi[k+1:n])
    # poly prende un vettore e mi restituisce un vettore con i coefficienti
    # di quel polinomio.
    # Ad esempio se faccio np.poly(np.array[2, 3]) mi restituisce (x-2)(x-3)
    # cioè x^2 -5x +6 sotto forma di vettore = [1, -5, 6].
    # i coefficienti sono il numeratore di p
    num = np.poly(xzeri)
    # polyval valuta un polinomio (num) nei punti dati (xnodi[k])
    # il denominatore di p è dato dal numeratore valutato in xi che sarebbe
    # l'indice del polinomio di base che sto valutando (l'argomento k)
    den = np.polyval(num, xnodi[k])
    p = num / den
    return p
# =============================================================================
# 
# =============================================================================
def InterpL(x, f, xx):
    """ InterpL(x, f, xx)
    
    Funzione che determina in un insieme di punti il valore del polinomio
    interpolante ottenuto dalla formula di Lagrange
        
    Parameters
    ----------
    x : vector
        Vettore con i nodi dell'interpolazione
    f : vector
        Vettore con i valori dei nodi
    xx : vector
        Vettore con i punti in cui si vuole calcolare il polinomio
        
    Returns
    -------
    vector
        y : vettore contenente i valori assunti dal polinomio interpolante
    
    """
    n = x.size
    m = xx.size
    L = np.zeros((n, m))
    for k in range(n):
        p = plagr(x, k)
        L[k, :] = np.polyval(p, xx)
    return np.dot(f, L)
    