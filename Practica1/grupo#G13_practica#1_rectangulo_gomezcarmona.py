# -*- coding: utf-8 -*-
"""
@authors: Lucas Gómez y Daniel Carmona
"""
import sympy as sp


x = sp.symbols('x')
soporte = [0,2*sp.pi]




def calcular_error_relativo(aproximacion, integral):
    return (abs(abs(aproximacion-integral))/abs(integral))
#fin de calcular_error

def formulaRectangulo(funcion,soporte):

    #Empieza el programa
    print("Nuestra función es: f(x) =", funcion(x))
    print("Queremos integrarla en:", soporte)
    #Calculamos el valor real
    I_r = (sp.integrate(funcion(x), (x, soporte[0], soporte[1]))).evalf()
    print("La integral exacta es I_r =", I_r)

    #Calculamos las aproximaciones de la integral para cada lado y sus errores relativos correspondientes
    extremoIzquierdo = (funcion(soporte[0])*(soporte[1]-soporte[0])).evalf()
    print("En el extremo izquierdo  la integral es ", extremoIzquierdo)
    puntoMedio = (funcion((soporte[0]+soporte[1])/(2))*(soporte[1]-soporte[0])).evalf()
    print("En el punto medio la integral es ",puntoMedio)
    extremoDerecho = (funcion(soporte[1])*(soporte[1]-soporte[0])).evalf()
    print("En el extremo derecho  la integral es ", extremoDerecho)

    #Calculamos y mostramos el error relativo
    error_relativo = calcular_error_relativo(extremoIzquierdo, I_r)
    print("El error relativo en el extremo izquierdo es", error_relativo.evalf())
    error_relativo = calcular_error_relativo(puntoMedio, I_r)
    print("El error relativo en el punto medio es", error_relativo.evalf())
    error_relativo = calcular_error_relativo(extremoDerecho, I_r)
    print("El error relativo en el extremo derecho es", error_relativo.evalf())

def funcion (t):
    return sp.cos((t**2)-1)





formulaRectangulo(funcion,soporte)


