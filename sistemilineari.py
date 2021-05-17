# -*- coding: utf-8 -*-
"""
Created on Mon May 17 15:28:14 2021

@author: eliam
"""

import numpy as np
# =============================================================================
# METODI RISOLUTIVI
# =============================================================================
def Lsolve(L, b):
    """ Lsolve(L, b)
    
    Risoluzione con procedura forward di Lx=b con L triangolare inferiore
    
    Parameters
    ----------
    L : matrix
        Matrice triangolare inferiore
    b : vector
        Vettore termine noto
        
    Returns
    -------
    list
        x : soluzione del sistema lineare
        flag : 0 se sono soddisfatti i criteri di applicabilità, 1 altrimenti
    
    """
    flag = 0
    # test dimensione
    m, n = L.shape
    if n != m:
        print("Lsolve: Matrice non quadrata")
        flag = 1
        return [], flag
    # test singolarità
    if np.all(np.diag(L)) != True:
        print("Lsolve: Matrice triangolare inferiore")
        flag = 1
        return [], flag
    # preallocazione vettore soluzione
    x = np.zeros((n, 1))
    for i in range(n):
        # scalare = vettore riga * vettore colonna
        # sarebbe la sommatoria aij*xj
        s = np.dot(L[i, :i], x[:i])
        # xi = (bi - risultato della sommatoria)/aii
        x[i] = (b[i] - s) / L[i, i]
    
    return x, flag
# =============================================================================
# 
# =============================================================================
def Usolve(U, b):
    """ Usolve(U, b)
    
    Risoluzione con procedura backward di Ux=b con U triangolare superiore
    
    Parameters
    ----------
    U : matrix
        Matrice triangolare superiore
    b : vector
        Vettore termine noto
        
    Returns
    -------
    list
        x : soluzione del sistema lineare
        flag : 0 se sono soddisfatti i criteri di applicabilità, 1 altrimenti
    
    """
    flag = 0
    # test dimensione
    m, n = U.shape
    if m != n:
        print("Usolve: Matrice non quadrata")
        flag = 1
        return [], flag
    # test singolarità
    if np.all(np.diag(U)) != True:
        print("Usolve: Elemento diagonale nullo")
        flag = 1
        return [], flag
    # preallocazione vettore soluzione
    x = np.zeros((n, 1))
    
    for i in range(n-1, -1, -1):
        # scalare = vettore riga * vettore colonna
        s = np.dot(U[i, i+1:n], x[i+1:n])
        # xi = (bi - risultato della sommatoria)/aii
        x[i] = (b[i] - s) / U[i, i]
        
    return x, flag
# =============================================================================
# Dopo essersi costruiti LU in PA = LU tramite fattorizzazione LU ora
# Il nostro obiettivo è risolvere Ax = b -> PAx = Pb -> LUx = Pb
# Poniamo {Ux = y} => Ly = Pb
# Quindi Ax = b equivale a risolvere due sistemi lineari
# y = LSolve(L, Pb)
# x = Usolve(U, y)
# =============================================================================
def LUsolve(L, U, P, b):
    """ LUsolve(L, U, P, b)
    
    Risoluzione a partire da PA = LU assegnata
    
    Parameters
    ----------
    L : matrix
        Matrice triangolare inferiore
    U : matrix
        Matrice triangolare superiore
    P : matrix
        Matrice identità
    b : vector
        Vettore termine noto
        
    Returns
    -------
    list
        x : soluzione del sistema lineare
        flag : 0 se sono soddisfatti i criteri di applicabilità, 1 altrimenti
    
    """
    # LUx = Pb pongo {Ux = y} => Ly = Pb 
    Pb = np.dot(P, b)
    # Trovo y
    y, flag = Lsolve(L, Pb)
    if flag == 0:
        # Ora posso risolvere Ux = y
        x, flag = Usolve(U, y)
    else:
        return [], flag
    return x, flag
# =============================================================================
# 
# =============================================================================
def LU_nopivot(A):
    """ LU_nopivot(A)
    
    Fattorizzazione A = PA = LU
    
    Parameters
    ----------
    A : matrix
        Matrice quadrata
        
    Returns
    -------
    list
        L : matrice triangolare inferiore
        U : matrice triangolare superiore
        P : matrice identità
        flag : 0 se sono soddisfatti i criteri di applicabilità, 1 altrimenti
    
    """
    flag = 0
    # test dimensione
    m, n = A.shape[0]
    if n != m:
        print("LU_nopivot: matrice non quadrata")
        flag = 1
        return [], [], [], flag
    # creo la matrice identità P e conservo A
    P = np.eye(n)
    # creo una copia di A per rendere U un oggetto distinto da A
    U = A.copy()
    # Fattorizzazione
    for k in range(n-1): # da 0 a n-2 perché n-1 è escluso
        # test pivot
        if U[k, k] == 0:
            print("LU_nopivot: elemento diagonale nullo")
            flag = 1
            return [], [], [], flag
        # memorizza i moltiplicatori
        U[k+1:n, k] /= U[k, k] 
        # eliminazione gaussiana sulla matrice:    
        # np.outer è il prodotto esterno, ossia il prodotto della colonna k
        # per la riga k che dà una matrice come risultato
        # il prodotto scalare è detto anche prodotto interno
        # moltiplica la riga k-esima per la riga dei moltiplicatori
        U[k+1:n, k+1:n] -= np.outer(U[k+1:n, k], U[k, k+1:n])
    # estrae i moltiplicatori e aggiunge la diagonale identità
    L = np.tril(U, -1) + np.eye(n)
    # estrae la parte triangolare superiore
    U = np.triu(U)
    return P, L, U, flag    
# =============================================================================
# METODI DI FATTORIZZAZIONE
# =============================================================================
# =============================================================================
# 
#     |a00 a01 a02|
#     |a10 a11 a12|
#     |a20 a21 a22|
#     
#     k=0
#         m10 = a10 / a00 
#         m10 * a1 -> annullo a10 -> a10 = 0
#         salvo m10 in a10 per risparmiare spazio
#         così nella triangolare inferiore (esclusa la diagonale) ho salvato
#         i moltiplicatori di lagrange.
#         tanto non avrei vantaggi nel tenere la triangolare inferiore piena
#         di zeri
# 
#       |\  |
#       | \U|   (quadrata)
#       |L \|
#
#     np.tril(U, -1) estrae la triangolare inferiore esclusa la diagonale
#     +np.eye(n) aggiunge l'identità per avere 1 sulla diagonale di L
#     |1 0 0|
#     |  1 0|
#     |L   1|
# =============================================================================
def LU_pivot(A):
    """ LU_pivot(A)
    
    Fattorizzazione A = PA = LU con pivoting parziale
    
    Parameters
    ----------
    A : matrix
        Matrice quadrata
        
    Returns
    -------
    list
        L : matrice triangolare inferiore
        U : matrice triangolare superiore
        P : matrice identità
        flag : 0 se sono soddisfatti i criteri di applicabilità, 1 altrimenti
    
    """
    flag = 0
    # test dimensione
    m, n = A.shape
    if m != n:
        print("LU_pivot: matrice non quadrata")
        flag = 1
        return [], [], [], flag
    # creo la matrice identità P
    P = np.eye(n)
    # creo una copia di A per rendere U un oggetto distinto da A
    U = A.copy()
    for k in range(n-1):
        # test pivot
        if U[k, k] == 0:
            print("LU_pivot: elemento diagonale nullo")
            flag = 1
            return [], [], [], 1

        # il k ci va perché argmax restituisce l'indice del sottovettore estratto
        # mentre nella funzione l'indice all'interno della matrice corrisponde al
        # risultato + k
        p = np.argmax(abs(U[k:n,k])) + k 
        if k != p:
            # U[[k,p],:] = U[[p,k],:] scambio la riga k e p
            # visto che scambiamo le righe, dobbiamo scambiare anche le stesse righe
            # nella matrice dei termini noti e ne teniamo traccia nella matrice identità
            U[[k, p], :] = U[[p, k], :]
            P[[k, p], :] = P[[p, k], :]
        # eliminazione gaussiana
        for i in range(k+1, n):
            # memorizza i moltiplicatori
            U[i, k] /= U[k, k]
            # eliminazione gaussiana sulla matrice
            U[i, k+1:n] -= U[i, k] * U[k, k+1:n]

        # estrae i moltiplicatori e aggiungere la diagonale identità
        L = np.tril(U, -1) + np.eye(n)
        # estra la triangolare superiore
        U = np.triu(U)
        return P, L, U, flag
        
        
        
        