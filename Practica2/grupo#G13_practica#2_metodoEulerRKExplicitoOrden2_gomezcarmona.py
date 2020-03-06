# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 11:52:12 2019

@author: Daniel y Lucas
"""
import sympy as sp
import mpmath as mp
import math as m
import matplotlib.pyplot as plt
import numpy as np
    
def f_real(t):
    x = (t+1)**2 - 0.5 * np.exp(t)
    return x
#fin f_real
    
def f_u(v):
    return v
#fin_f_u
    
def f_v(u, v, B, k, M):
    v1 = ((-B * np.abs(v) * v) - (k * u))/M
    return v1
#fin_f_v
    
def calcular_N(intervalo, h):
    #Calculamos el número de puntos y lo devolvemos
    N = int(1 + (intervalo[1]-intervalo[0])/h)
    return N
#fin calcular_N
    
def calcular_soporte(intervalo, h, N):
    soporte = list()
    #Calculamos desde t_1 hasta t_N
    for i in range(1, N+1):
        t_i = intervalo[0] + (i-1)*h
        soporte.append(t_i)
    #Devolvemos el soporte
    return soporte
#fin calcular_soporte

def calcular_error(t, y):
    #Calculamos el error absoluto
    error = abs(y - f_real(t))
    error_relativo = error/f_real(t)
    return error_relativo
#fin calcular_error

def metodo_Euler_Muelle(f_u, f_v, intervalo, h, u_inicial, v_inicial, M, B, k):
    #Calculamos N
    N = calcular_N(intervalo, h)
    #Creamos la solución
    soluciones = list()
    soporte = calcular_soporte(intervalo, h, N)
    #Añadimos el punto inicial y_1
    soluciones.append([soporte[0], u_inicial, v_inicial])
    #Calculamos de y_2 a y_N
    for i in range(1, N):
        solucion = list()
        #Cogemos la u y la v anterior
        u_anterior = soluciones[i-1][1]
        v_anterior = soluciones[i-1][2]
        #Calculamos la u y la v actual
        u = u_anterior + h*f_u(v_anterior)
        v = v_anterior + h*f_v(u_anterior, v_anterior, B, k, M)
        #Calculamos el error
        #Guardamos la solución de esta iteración
        solucion = [soporte[i], u, v]
        #La añadimos al conjunto de iteraciones
        soluciones.append(solucion)
        
    return soluciones

def mostrar_salida(lista):
    Tabla = """\
+------------------------------------------------------------------------------+
|t_n                        u_n                       v_n                      |
|------------------------------------------------------------------------------|
{}
+------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<25} {:<25} {:<25}|".format(*fila)
    for fila in lista)))
    print(Tabla)
#fin mostrar_salida
    
def pintar_funciones(soluciones, soporte, M, B, k):
    u = list()
    uprima = list()
    for s in soluciones:
        u.append(s[1])
        uprima.append(s[2])
    
    #Definimos el tamaño de los puntos de la función y las representamos
    tamañoPunto = 3
    plt.figure(figsize=(9,6))
    plt.plot()
    plt.plot(soporte ,u ,'o', markersize=tamañoPunto, label='u  (m)')
    plt.plot(soporte, uprima, 'ro', markersize=tamañoPunto, label='u\' (m/s)')
    
    plt.title('Oscilación muelle')

    # Add X and y Label
    plt.xlabel('Tiempo (segundos)')
    plt.ylabel('Altura (metros)')

    # Add a grid
    plt.grid(alpha=.4,linestyle='--')
    plt.text(6.1, 1, 'M={} kg'.format(M))
    plt.text(6.1, .9, 'B={} Ns\u00b2/m\u00b2'.format(B))
    plt.text(6.1, .8, 'k={} N/m'.format(k))
    plt.legend()
    plt.show()
#fin pintar_funciones
    
if __name__ == "__main__":
    T = 10 #segundos
    M = 10 #kg
    B = 50 #Ns2/m2
    k = 200#N/m
    
    intervalo = [0,T]
    h = 0.01
    u_inicial = 0
    v_inicial = 1
    
    N = calcular_N(intervalo, h)
    soporte = calcular_soporte(intervalo, h, N)
    
    print("La aproximación de problema para los valores\nM={} kg, B={} Ns\u00b2/m\u00b2, k={} N/m en el intervalo de [0, {}] segundos con h={}".format(M, B, k, T, h))
    soluciones = metodo_Euler_Muelle(f_u, f_v, intervalo, h, u_inicial, v_inicial, M, B, k)
    mostrar_salida(soluciones)
    pintar_funciones(soluciones,soporte, M, B, k)