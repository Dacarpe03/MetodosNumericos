# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 17:27:54 2019

@authors: Lucas Gómez y Daniel Carmona
"""

import sympy as sp
x = sp.symbols('x')

def funcion(t):
    return t**3 - 3*t + 2
#fin_funcion
    
def calcular_k(en, en_anterior, k):
    return abs(en)/(abs(en_anterior**k))
#fin calcular_k
    
#Función para calcular la raíz mediante el metodo newton raphson
def metodo_newton_raphson(pto_inicial, tolerancia, funcion, raiz_real):
    #Declaramos una lista para guardar la información de las iteraciones
    info_iteraciones = []
    
    #Guardamos la derivada de la funcion
    derivada = sp.diff(funcion(x),x)
    
    #Variable para el error
    error = raiz_real - pto_inicial
    info_iteraciones.append([0, pto_inicial, error, 0, 0, 0])
    #Variables para las x_i de las iteraciones
    x_i = pto_inicial;
    x_i1 = 0;
    
    iteracion = 1
    
    #Ahora empezamos el bucle para buscar la raíz
    while(abs(error) > tolerancia):
        #Calculamos valores para esta iteración
        valorxi = sp.expand(funcion(x_i)).evalf()
        valorderivadaxi= sp.expand(derivada.subs(x,x_i)).evalf()
        #Calculamos x_i+1
        x_i1 = x_i - valorxi/valorderivadaxi
        #Calculamos el error en esta iteración
        error = x_i1 - x_i
        #Calculamos En en la iteración
        e_n = x_i1 - raiz_real
        #Calculamos los en de la iteración anterior
        k1 = calcular_k(e_n, info_iteraciones[iteracion-1][2], 1)
        k2 = calcular_k(e_n, info_iteraciones[iteracion-1][2], 2)
        k3 = calcular_k(e_n, info_iteraciones[iteracion-1][2], 3)
        #Guardamos los en de la iteración anterior
        info_iteraciones[iteracion-1][3] = k1
        info_iteraciones[iteracion-1][4] = k2
        info_iteraciones[iteracion-1][5] = k3
        #Guardamos la información de esta iteración
        info_it_i = [iteracion, x_i1, e_n, 0, 0, 0]
        info_iteraciones.append(info_it_i)
        #Siguiente iteración
        iteracion+=1
        x_i = x_i1
        
    return info_iteraciones

def mostrar_salida(lista):
    Tabla = """\
+-----------------------------------------------------------------------------------------------------------------------------------------------+
| n          x_i                       En                        k=1                       k=2                       k=3                        |
|-----------------------------------------------------------------------------------------------------------------------------------------------|
{}
+-----------------------------------------------------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<10} {:<40} {:<40} {:<40} {:<40} {:<40}  |".format(*fila)
    for fila in lista)))
    print(Tabla)
    
if __name__ == "__main__":
    pto_inicial = -3
    raiz_real = -2
    tolerancia = 1/1000000
    print("Comprobamos la convergencia en la función f(x)={} a partir del punto x={} con una tolerancia de {}".format(funcion(x), pto_inicial, tolerancia))
    lista_iteraciones = metodo_newton_raphson(pto_inicial, tolerancia, funcion, raiz_real)
    for l in lista_iteraciones:
        print("Iteración {}:   Aproximación={}   En={}   k=1: {}   k=2: {}   k=3: {}".format(l[0], l[1], l[2], l[3], l[4], l[5]))
    #mostrar_salida(lista_iteraciones)