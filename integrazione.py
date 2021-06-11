# -*- coding: utf-8 -*-
"""
Created on Tue May 25 18:01:35 2021

@author: eliam
"""

import numpy as np

# Di seguito sono implementati gli algoritmi Trapezio e Simpson compositi
# utili al calcolo dell'integrazione numerica
# I(f; a, b) con I integrale, f : [a, b] -> R polinomio interpolante,
# a e b limiti dell'intervallo
# L'espressione generale delle formule di quadratura è data dalla
# sommatoria da i=1 a n di ( c_i * f(x_i) ) con x_i nodi appartenenti ad [a, b],
# c_i coefficienti fissati (detti pesi)
#
# Nella formula del trapezio si approssima la f con un polinomio di grado 1 (retta)
# graficando di fatto un trapezio formato dai 4 punti (a,0);(a,f(a));(b,f(b));(b,0) 
# La formula del trapezio è data da:
#   I(f) = Integrale di P1(x) con
#       P1(x) il polinomio di grado 1 che interpola (a, f(a)), (b, f(b))
#       P1(x) = [(f(b) - f(a))/(b - a)]*x + (b*f(a) - a*f(b))/(b - a)
#   => I(f) = (b - a)/2 * (f(a) + f(b)) con (b-a)/2 l'altezza del trapezio
# L'errore (resto) è dato da:
#   R(f) = -1/12 * (b - a)^3 * f"(csi)  
#
# Nella formula di Simpson si approssima f con un polinomio di grado 2 (parabola)
# graficando di fatto un """rettangolo""" con il """tetto""" formato da una parabola
# formato sempre dai 4 punti: (a,0);(a,f(a));(b,f(b));(b,0) a cui si aggiunge
# il vertice della parabola in ((a+b)/2, f((a+b)/2))
# La formula di Simpson è data da
#   I(f) = Integrale di P2(x) con
#       P2(x) = 1/(b-a)^2 * [(a+b)*b*f(a) -4*a*b*f((a+b)/2) +(a+b)*a*f(b)]
#   => I(f) = (b - a)/6 * (f(a) + 4*f((a + b)/2) + f(b))
# L'errore (resto) è dato da:
#   R(f) = -1/90 * [(b-a)/2]^5 * f""(csi)
#
# Per migliorare l'accuratezza del calcolo d'integrazione numerica
# si possono utilizzare le formule su scritte in N sottointervalli 
# dell'intervallo [a, b]. Queste nuove formule sono dette composite.
# => Si suddivide l'intervallo di integrazione [a, b] in N sottointervalli [z_k, z_k+1]
# con k = 0...N-1 tutti con la stessa ampiezza. In ciascuno dei sotto intervalli
# si applica la forma del trapezio o di Simpson. 
# In generale, si può scrivere che la formula composita è data dalla sommatoria
# da k=0 a N-1 della formula d'integrazione calcolata nell'intervallo k che è
# (I_n)'k con k intervallo e n formula di quadratura (n=1 trapezio, n=2 Simpson)
# Il resto, similmente, è dato dalla sommatoria da k=0 a N-1 del resto calcolato
# in ogni sottointervallo, che è
# (R_n)'k con k intervallo e n formula di quadratura (n=1 trapezio, n=2 Simpson)
# 
# => La formula dei trapezi composita è data da:
#   I(f) = [(b - a)*2*N] * [f(a) + 2*SUM_k=1->N-1[f(z_k)] + f(b)]
# e il suo resto è quindi:
#   R(f) = -1/12 * [(b-a)^3]/N^2 * f"(csi)
#
# => La formula di Simpson composita è data da:
#   I(f) = [(b - a)*6*N] * [f(a) + 2*SUM_k=1->N-1[f(z_k)] + 4*SUM_k=1->N-1{f[(z_k+z_k+1)/2]} + f(b)]
# e il suo resto è quindi:
#   R(f) = -1/2880 * [(b-a)^5]/N^4 * f""(csi)


def TrapComp(fname, a, b, n):
    h = (b - a)/n
    nodi = np.arange(a, b+h, h)
    f = fname(nodi)
    I = (f[0] + 2*np.sum(f[1:n]) + f[n]) * h/2
    return I


def SimpComp(fname, a, b, n):
    h = (b - a)/(2*n)
    nodi = np.arange(a, b+h, h)
    f = fname(nodi)
    I = (f[0] + 2*np.sum(f[2 : 2*n : 2]) + 4*np.sum(f[1 : 2*n : 2]) + f[2*n]) * h/3
    return I

def traptoll(fun, a, b, tol):
    nmax = 2048
    err = 1
    N = 1
    IN = TrapComp(fun, a, b, N)
    while N<=nmax and err>tol:
        N = 2*N
        I2N = TrapComp(fun, a, b, N)
        err = abs(IN - I2N)/3
        IN = I2N
    if N > nmax:
        print("traptoll: raggiunto nmax di intervalli")
        N = 0
        IN = []
    return IN, N

def simptoll(fun, a, b, tol):
    nmax = 2048
    err = 1
    N = 1
    IN = SimpComp(fun, a, b, N)
    while N<=nmax and err>tol:
        N = 2*N
        I2N = SimpComp(fun, a, b, N)
        err = abs(IN - I2N)/15
        IN = I2N
    if N > nmax:
        print("simptoll: raggiunto nmax di intervalli")
        N = 0
        IN = []
    return IN, N