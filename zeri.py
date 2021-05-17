# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:17:33 2021

@author: eliam
"""

import numpy as np
import math
from utils import sign
# =============================================================================
# 
# =============================================================================
def bisez(fname, a, b, tol):
    """ Metodo di bisezione per la ricerca di zeri di funzione.

    Parameters
    ----------
    fname : lambdified function
        Funzione principale
    a : number
        Limite sinistro dell'intervallo
    b : number
        Limite destro dell'intervallo
    tol : number
        tolleranza per il test d'arresto sull'incremento
        
    Returns
    -------
    list
        x : soluzione
        it : numero di iterazioni
        xk : iterabili che convergono allo zero della funzione
    
    """

    eps = np.spacing(1) # np.spacing(x) Restituisce la distanza tra x e il numero adiacente più vicino.
                        # np.spacing(1)  restituisce quindi l' eps di macchina.
    fa = fname(a)
    fb = fname(b)
    if sign(fa) == sign(fb):
        print("Intervallo non corretto. Algoritmo non applicabile.")
        return [], 0, []
    else:
        maxit = math.ceil(math.log((b - a)/tol)/math.log(2))
        print("Numero massimo di iterazioni: ", maxit)
        it = 0
        xk = []
        while it < maxit and abs(b - a) >= tol + eps*max([abs(a), abs(b)]):
            m = a + (b - a)/2 # punto medio fra a e b
            xk.append(m) # salvo in xk i vari valori di m
            fm = fname(m) # valore della funzione nell'ascissa m 
            if fm == 0: 
                # se fm è 0, sono sull'ascissa quindi ho trovato x
                break
            elif sign(fm) == sign(fa): 
                # se fm ha lo stesso segno di fa vuol dire che devo prendere
                # m come nuovo limite sinistro dell'intervallo
                a = m # il punto medio diventa il nuovo limite sinistro
                fa = fm
            elif sign(fm) == sign(fb):
                # se fm ha lo stesso segno di fb vuol dire che devo prendere
                # m come nuovo limite destro dell'intervallo
                b = m # il punto medio diventa il nuovo limite destro
                fb = fm
            
            it += 1
        x = m 
        return x, it, xk
# =============================================================================
#     
# =============================================================================
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
        x0 = x1
        fx0 = fx1 
        x1 = x0 - fx0 / m
        fx1 = fname(x1)
        xk.append(x1)
        it += 1
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk
# =============================================================================
#             
# =============================================================================
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
        it = 1 
    else:
        print("Derivata nulla in x0.")
        return [], 0, []
        
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
            fx1 = fname(x1)
            xk.append(x1)
            it += 1
        else:
            print("Derivata nulla in x0.")
            return x1, it, xk
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk
# =============================================================================
# 
# =============================================================================
#Newton Modificato
def newton_m(fname, fpname, x0, m, tolx, tolf, nmax):
    """ Metodo di Newton modificato per la ricerca di zeri di funzione.
    
    Descrizione della funzione
    
    Parameters
    ----------
    fname : lambdified function
        Funzione principale
    fpname : lambdified function
        Derivata prima della funzione principale
    x0 : number
        Iterata iniziale
    m : number
        Molteplicità della radice
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
    xk.append(x0)
    fx0 = fname(x0)
    # il valore di m nel metodo di newton è la derivata della funzione in x0
    dfx0 = fpname(x0)
    # se la derivata non è nulla proseguo con l'algoritmo
    if abs(dfx0) > np.spacing(1):
        # geometricamente si prende come nuova approssimazione l'intersezione dell'asse
        # delle ascisse con la retta tangente a f in (x0, f(x0))
        x1 = x0 - m * fx0 / dfx0 # calcolo xi +1
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
        dfx0 = fpname(x0)
        if abs(dfx0) > np.spacing(1): 
            fx0 = fx1 
            x1 = x0 - m * fx0 / dfx0
            fx1 = fname(x1)
            xk.append(x1)
            it += 1
        else:
            print("Derivata nulla in x0.")
            return x1, it, xk
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk
# =============================================================================
# 
# =============================================================================
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
# =============================================================================
# 
# =============================================================================
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
        fx1 = fname(x1)
        xk.append(x1)
        it += 1
    if it >= nmax:
        print("Raggiunto numero massimo di iterazioni possibili.")
    return x1, it, xk