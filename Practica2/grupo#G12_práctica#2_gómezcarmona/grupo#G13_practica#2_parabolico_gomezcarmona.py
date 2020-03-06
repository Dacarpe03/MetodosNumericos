# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 11:29:42 2019

@author: Daniel y Lucas
"""
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import mpmath as mp
import math as ma
import matplotlib.pyplot as plt
import numpy as np

def phi(t):
    if(t <= 0.5):
        return 2*t
    else:
        return 2*(1-t)

def calcular_N_o_M(intervalo, paso):
    #Calculamos el número de puntos y lo devolvemos
    N = int(1 + (intervalo[1]-intervalo[0])/paso)
    return N
#fin calcular_N_o_M
    
def calcular_soporte(intervalo, paso, N):
    soporte = list()
    #Calculamos desde t_1 hasta t_N
    for i in range(1, N+1):
        t_i = intervalo[0] + (i-1)*paso
        soporte.append(t_i)
    #Devolvemos el soporte
    return soporte
#fin calcular_soporte
    
def crear_soporte_bidimensional(soporte_espacio, soporte_tiempo):
    #Vamos a crear una matriz con los puntos del soporte haciendo todas las combinaciones de los puntos de los soportes
    #Obtendremos una matriz NxM
    soporte_bidimensional = list()
    for e in soporte_espacio: #m
        fila = list()
        for t in soporte_tiempo: #n 
            #En la posicion n m guardamos los valores del soporte_espacio n y soporte_tiempo m
            soportenm = [e,t]
            fila.append(soportenm)
        soporte_bidimensional.append(fila)
    return soporte_bidimensional
#fin crear_soporte_bidimensional
    
def calcular_r(h, tau):
    r = (tau/h**2)
    return r
#fin calcular_r
    
def resolver_edp_parabolica(M, N, r, soporte_espacio, soporte_tiempo, u0, u1, phi):
    s = list()
    #Calculamos la primera columna para n=1
    columna1 = list()
    for m in range(M):
        um1 = phi(soporte_espacio[m])
        columna1.append(um1)
    s.append(columna1)
    
    #Vamos a por las demás columna
    for n in range(1, N):
        columnan = list()
        for m in range(M):
            #Si la fila es 1 entonces u1n = u0
            if (m==0):
                umn = u0
            #Si la fila es M entonces uMn = u1
            elif (m==M-1):
                umn = u1
            else:
                a = s[n-1][m-1] #u m-1, n-1
                b = s[n-1][m]   #u m,   n-1
                c = s[n-1][m+1] #u m+1, n-1
                #umn = r*u m-1,n-1 + (1-2r)*u m,n-1 +  r*u m+1, n-1
                umn = r*a + (1-2*r)*b + r*c
            columnan.append(umn)
        s.append(columnan)
    return s
    
def mostrar_matriz(A):
    print('   M   N-->  1 %s' % (' '.join('%8s' % str(j+1) for j in range(1, len(A[0])))))
    i = 1
    for row in A:
        print ('{:4} |%s|'.format(i) % (' '.join('%8s' % round(l,4) for l in row)))
        i += 1
#fin mostrar_matriz
        
        
def tras(matriz):
    '''
    Encuentra la transpuesta de un matriz
    assumiendo hileras de la misma longitud
    '''
    hileras = len (matriz)
    columnas = len(matriz[0])
    t = [[0 for x in range(hileras)] for y in range(columnas)] 

    for i in range(hileras):
        for j in range(columnas):
            t[j][i] = matriz[i][j]
    return t
#fin tras
    
def pintar_grafica3d(vector_a_pintar):
    espacio  = vector_a_pintar[0]
    tiempo   = vector_a_pintar[1]
    solucion = vector_a_pintar[2]
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    ax.scatter3D(espacio, tiempo, solucion, c=solucion);
    ax.set_xlabel('Espacio')
    ax.set_ylabel('Tiempo')
#fin pintar_grafica3d
    
def mostrar_salida(lista):
    Tabla = """\
+---------------------------------------------------------------------------------------------------------------+
| m               n               xm                        tn                        u_mn                      |
|---------------------------------------------------------------------------------------------------------------|
{}
+---------------------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<15} {:<15} {:<25} {:<25} {:<25} |".format(*fila)
    for fila in lista)))
    print(Tabla)
#fin mostrar_salida
    
def crear_tabla(sol_buena, soporte_espacio, soporte_tiempo, M, N):
    vector_tabla = list()
    vector_espacio = list()
    vector_tiempo = list()
    vector_soluciones = list()
    
    for m in range(M):
        for n in range(N):
            vector_tabla.append([m, n, soporte_espacio[m], soporte_tiempo[n], sol_buena[m][n]])
            vector_espacio.append(soporte_espacio[m])
            vector_tiempo.append(soporte_tiempo[n])
            vector_soluciones.append(sol_buena[m][n])
            
    mostrar_salida(vector_tabla)
    
    vector_grafica = [vector_espacio, vector_tiempo, vector_soluciones]
    
    return vector_grafica
#fin crear_tabla
    
if __name__ == "__main__":
    intervalo_espacio = [0,1]
    h = 0.02
    
    intervalo_tiempo = [0,1]
    tau = 0.01
    
    M = calcular_N_o_M(intervalo_espacio, h)
    N = calcular_N_o_M(intervalo_tiempo, tau)
    
    u0 = 0.0;
    u1 = 0.0;
    
    soporte_espacio = calcular_soporte(intervalo_espacio, h, M)
    soporte_tiempo = calcular_soporte(intervalo_tiempo, h, N)
    
    r = calcular_r(h, tau)
    
    print("Nuestro intervalo de espacio es ({},{}) que discretizando con un paso h={} nos da el soporte:".format(intervalo_espacio[0], intervalo_espacio[1], h))
    print(soporte_espacio)
    print("Nuestro intervalo de tiempo es ({},{}) que discretizando con un paso \u03C4={} nos da el soporte:".format(intervalo_tiempo[0], intervalo_tiempo[1], tau))
    print(soporte_tiempo)
    print("Tenemos una r={}\n".format(round(r,4)))
     
    soluciones = resolver_edp_parabolica(M, N, r, soporte_espacio, soporte_tiempo, u0, u1, phi)
 
    sol_buena = tras(soluciones)
    mostrar_matriz(sol_buena)
    
    vectores_para_pintar = crear_tabla(sol_buena, soporte_espacio, soporte_tiempo, M, N)
    pintar_grafica3d(vectores_para_pintar)