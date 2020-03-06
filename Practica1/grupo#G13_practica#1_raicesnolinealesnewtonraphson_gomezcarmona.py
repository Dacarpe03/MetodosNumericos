# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 18:10:38 2019

@authors: Lucas Gómez y Daniel Carmona
"""

import sympy as sp
import mpmath as mp
import math as m
x = sp.symbols('x')

def funcion(t):
    return t*sp.sin((t**2)/2) + sp.exp(-t)
#fin_funcion
    
#Función para calcular la raíz mediante el metodo newton raphson
def metodo_newton_raphson(pto_inicial, tolerancia, funcion):
    #Declaramos una lista para guardar la información de las iteraciones
    info_iteraciones = []
    
    #Guardamos la derivada de la funcion
    derivada = sp.diff(funcion(x),x)
    
    #Variable para el error
    error = abs(pto_inicial)
    
    #Variables para las x_i de las iteraciones
    x_i = pto_inicial;
    x_i1 = 0;
    
    iteracion = 1
    
    #Ahora empezamos el bucle para buscar la raíz
    while(error > tolerancia):
        
        #Calculamos valores para esta iteración
        valorxi = sp.expand(funcion(x_i)).evalf()
        valorderivadaxi= sp.expand(derivada.subs(x,x_i)).evalf()
        
        #Calculamos x_i+1
        x_i1 = x_i - valorxi/valorderivadaxi
        #Calculamos el error en esta iteración
        error = abs(x_i1-x_i)
        #Guardamos la información de esta iteración
        info_it_i = [iteracion, x_i, x_i1, error, sp.expand(funcion(x_i1)).evalf()]
        info_iteraciones.append(info_it_i)
        
        #Siguiente iteración
        iteracion+=1
        x_i = x_i1
        
    return info_iteraciones

def mostrar_salida(lista):
    Tabla = """\
+-----------------------------------------------------------------------------------------------------------+
| Iteracion       x_i-1                x_i                  error                          f(x_i)                                 |
|-----------------------------------------------------------------------------------------------------------|
{}
+-----------------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<15} {:<20} {:<20} {:<30} {:<20}  |".format(*fila)
    for fila in lista)))
    print(Tabla)
#fin mostrar_salida
    

    
if __name__ == "__main__":
    pto_inicial = 4.5
    tolerancia = 1/100000
    print("Con el método Newton-Raphson amos a buscar la raíz de la función f(x)={} a partir del punto x={} con una tolerancia de {}".format(funcion(x), pto_inicial, tolerancia))
    lista_iteraciones = metodo_newton_raphson(pto_inicial, tolerancia, funcion)
    mostrar_salida(lista_iteraciones)