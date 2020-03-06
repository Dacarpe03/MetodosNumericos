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

def f(t,y):
    x = y - t**2 + 1
    return x
#fin f

def f_real(t):
    x = (t+1)**2 - 0.5 * np.exp(t)
    return x
#fin f_real

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

def metodo_Euler(f, f_real, intervalo, h, valor_inicial):
    #Calculamos N
    N = calcular_N(intervalo, h)
    #Calculamos el soporte
    soporte = calcular_soporte(intervalo, h, N)
    #Creamos la solución
    soluciones = list()
    #Añadimos el punto inicial y_1
    soluciones.append([1, soporte[0], valor_inicial, valor_inicial, 0])
    #Calculamos de y_2 a y_N
    for i in range(1, N):
        solucion = list()
        #El t_n actual
        t_i = soporte[i]
        #La solución anterior
        y_i_anterior = soluciones[i-1][3]
        #Calculamos la solución actual
        #Esta fórmula la hemos sacado despejando la fórmula de Euler implícita
        y_i = (y_i_anterior -  (h *(t_i)**2)+ h)/(1-h)
        #Calculamos el error
        error = calcular_error(t_i, y_i)
        #Guardamos la solución de esta iteración
        solucion = [i+1, t_i, f_real(t_i), y_i, error]
        #La añadimos al conjunto de iteraciones
        soluciones.append(solucion)

    return soluciones

def mostrar_salida(lista):
    Tabla = """\
+----------------------------------------------------------------------------------------------------+
| n               t_n                  f_real               y_n                 error  
|----------------------------------------------------------------------------------------------------|
{}
+----------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<15} {:<20} {:<20} {:<20} {:<20}|".format(*fila)
    for fila in lista)))
    print(Tabla)
#fin mostrar_salida

def pintar_funciones(soluciones, soporte ,f_real):
    y = list()
    for s in soluciones:
        y.append(s[3])

    t = np.arange(soporte[0], soporte[len(soporte)-1], 0.01)
    plt.plot()
    plt.plot(soporte ,y , 'bs', label='Método Euler')
    plt.plot(t, f_real(t), label='Función real')

    plt.title('Comparación')

    # Add X and y Label
    plt.xlabel('Eje X')
    plt.ylabel('Eje Y')

    # Add a grid
    plt.grid(alpha=.4,linestyle='--')

    plt.legend()
    plt.show()
#fin pintar_funciones



if __name__ == "__main__":

    intervalo = [0,1]
    h = 0.01
    valor_inicial = 0.5
    solucion = metodo_Euler(f, f_real, intervalo, h, valor_inicial)
    mostrar_salida(solucion)

    N = calcular_N(intervalo, h)
    soporte = calcular_soporte(intervalo, h, N)
    pintar_funciones(solucion, soporte, f_real)
