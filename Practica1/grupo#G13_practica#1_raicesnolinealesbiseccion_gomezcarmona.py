# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 11:09:07 2019

@authors: Lucas Gómez y Daniel Carmona
"""
import sympy as sp
import mpmath as mp
import math as m
x = sp.symbols('x')

def funcion(t):
    return t*sp.sin((t**2)/2) + sp.exp(-t)
#fin_funcion
    
def calcular_pto_medio(a, b):
    return (a+b)/2
#fin calcular_pto_medio
    
#Función para calcular el número de iteraciones mínimas
def calcular_iteraciones(a, b, tolerancia):
    #Calculamos la n con la fórmula vista en clase
    n = m.ceil(mp.log((b-a)/tolerancia)/mp.log(2))
    return n
#fin calcular_iteraciones
    
def tienen_mismo_signo(x, y):
    return ((x<0 and y<0) or (x>0 and y>0))

def metodo_biseccion(extremo_inf, extremo_sup, tolerancia, funcion):
    pto_medio = extremo_sup
    #Declaramos una lista para guardar la información de las iteraciones
    info_iteraciones = []
    
    #Primero calculamos el número de iteraciones
    n = calcular_iteraciones(extremo_inf, extremo_sup, tolerancia)
    
    #Variable para el error
    error = extremo_sup - extremo_inf;
        
    #Ahora empezamos el bucle para buscar la raíz
    for i in range(n):
        #Calculamos el error de esta iteración (pto_medio de la anterior iteracion - pto medio de esta iteracion)
        error = abs(pto_medio - calcular_pto_medio(extremo_inf, extremo_sup))
        print(error)
        #Calculamos el punto medio de la iteración i
        pto_medio = calcular_pto_medio(extremo_inf, extremo_sup)
        
        #Calculamos los valores de la funcion en los extremos y en el punto medio
        valor_sup = sp.expand(funcion(extremo_sup)).evalf()
        valor_inf = sp.expand(funcion(extremo_inf)).evalf()
        valor_medio = sp.expand(funcion(pto_medio)).evalf()
        
        #Guardamos la información de esta iteracion
        info_it_i = [i+1, extremo_inf, extremo_sup, pto_medio, valor_medio]
        info_iteraciones.append(info_it_i)
        
        #Ahora hacemos más pequeño el intervalo de búsqueda
        if(valor_inf == 0 or error <= tolerancia):#Si la función en el pto_medio es 0 es que es raíz, si el error es menor que la tolerancia paramos
            return info_iteraciones
        #Si f(extremo superior) tiene el mismo signo que f(punto medio) el nuevo intervalo será [extremo inferior, punto medio]
        elif(tienen_mismo_signo(valor_sup, valor_medio)):
            extremo_sup = pto_medio
        #Si f(extremo inferior) tiene el mismo signo que f(punto medio) el nuevo intervalo será [punto medio, extremo superior]
        else:
            extremo_inf = pto_medio
        
        
    return info_iteraciones

def mostrar_salida(lista):
    Tabla = """\
+-----------------------------------------------------------------------------------------------------------+
| Iteracion           Extremo_Izq           Extremo_Der          Pto_Medio           f(Pto.Medio)       |
|-----------------------------------------------------------------------------------------------------------|
{}
+-----------------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<20} {:<20} {:<20} {:<20} {:<30}  |".format(*fila)
    for fila in lista)))
    print(Tabla)
    
if __name__ == "__main__":
    punto_inicial_inferior = 3;
    punto_inicial_superior = 4;
    tolerancia = 1/100000;
    print("Vamos a calcular una raíz de f(x)={} entre los puntos x={} y x={} con una tolerancia e={}".format(funcion(x), punto_inicial_inferior, punto_inicial_superior, tolerancia))
    info = metodo_biseccion(punto_inicial_inferior, punto_inicial_superior, tolerancia, funcion)
    mostrar_salida(info)