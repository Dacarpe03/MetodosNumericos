# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:13:00 2019

@authors: Lucas Gómez y Daniel Carmona
"""
import sympy as sp
import fractions as f
x = sp.symbols('x')

#Definimos la función a integrar
def funcion_uno(t):
    return sp.cos((t**2)-1)
#fin de funcion_uno

#Definimos la función de la parte opcional
def funcion_opcional(t):
    return (sp.cos(t) - t*sp.sin(t))
#fin de funcion_opcional

#Función para calcular los x_i
def calcular_puntos_del_soporte(n, soporte):
    #Calculo la distancia entre los nuevos puntos del soporte
    num_puntos = n+1;
    distancia = (soporte[1]-soporte[0])/n;
    #Creo un nuevo soporte vacío
    nuevo_soporte = list()
    #Bucle para calcular los puntos del soporte
    for i in range(num_puntos):
        nuevo_punto = soporte[0] + i*distancia
        nuevo_soporte.append(nuevo_punto)
    return nuevo_soporte
#fin de calcular_puntos_del_soporte

def calcular_h(a, b, n):
    return (a-b)/n
#fin de calcular_h

def formula_Newton_Cotes_cerrada(diccionario, soporteab, funcion):
    #Cojo los extremos del soporte
    a = soporteab[0]
    b = soporteab[1]
    #Cojo los valores de la tabla
    n = diccionario["n"]
    c_i = diccionario["c_i"]
    D = diccionario["D"]
    #Calculo los puntos del nuevo soporte según n
    sop = calcular_puntos_del_soporte(n, soporteab)
    #I es para guardar el resultado de la aproximación
    I=0
    #Calculo el sumatorio
    for i in range(n+1):#Llega hasta n
        I+=c_i[i]*funcion(sop[i])
    #Ahora terminamos la aproximación multiplicando por (b-a)/D el sumatorio obtenido
    I*=(b-a)/D

    return I
#fin de formula_Newton_Cotes_cerrada

def calcular_error_relativo(aproximacion, integral):
    return (abs(abs(aproximacion-integral))/abs(integral))
#fin de calcular_error

def formula_compuesta_N(diccionario, soporteab, funcion, N):
    I = 0;
    #Vamos a dividir el intervalo subintervalos de la misma longitud.
    #Para ello calcularemos puntos equidistantes que cogeremos como extremos de los subintervalos
    subintervalo = calcular_puntos_del_soporte(N, soporteab)
    puntero_a = 0;
    puntero_b = 1;
    for i in range(N):
        #Calculamos los extremos del nuevo subintervalo
        a = subintervalo[puntero_a]
        b = subintervalo[puntero_b]

        #Hacemos la formula para el subintervalo [a,b]
        I+=formula_Newton_Cotes_cerrada(diccionario, [a,b], funcion)

        #Pasamos al siguiente subintervalo
        puntero_a += 1
        puntero_b += 1

    return I
#fin formula_compuesta

#Defino los valores de los coeficientes de la tabla de cerrados en diccionarios
diccionario_trapecio = {"Nombre": "Trapecio", "n":1, "c_i":[1,1], "D": 2, "K": f.Decimal(1/12)}
diccionario_simpson_1_3 = {"Nombre": "Simpson 1/3", "n":2, "c_i":[1,4,1], "D": 6, "K": f.Decimal(1/90)}
diccionario_simpson_3_8 = {"Nombre": "Simpson 3/8", "n":3, "c_i":[1, 3, 3, 1], "D": 8, "K": f.Decimal(3/80)}
diccionario_milne = {"Nombre": "Milne", "n":4, "c_i":[7, 32, 12, 32, 7], "D": 90, "K": f.Decimal(8/945)}
diccionario_sinnombre = {"Nombre": "Sin nombre", "n":5, "c_i":[19, 75, 50, 50, 75, 19], "D": 288, "K": f.Decimal(275/12096)}
diccionario_weddle = {"Nombre": "Weddle", "n":6, "c_i":[41, 216, 27, 272, 27, 216, 41], "D": 840, "K": f.Decimal(9/1400)}

tabla_cerradas = (diccionario_trapecio, diccionario_simpson_1_3, diccionario_simpson_3_8, diccionario_milne, diccionario_sinnombre, diccionario_weddle)

#Defino el soporte
soporte = [0,2*sp.pi]

#Empieza el programa
print("Nuestra función es: f(x) =", funcion_uno(x))
print("Queremos integrarla en:", soporte)
print("Elige con qué método:")
for i in range(6):
    print("{}.".format(i+1) ,tabla_cerradas[i]["Nombre"])

elegido = int(input("Introduce un número del 1 al 6 y pulsa intro:\n"))
print("Has elegido el método de", tabla_cerradas[elegido-1]["Nombre"])

#Calculamos la aproximación
aproximacion = formula_Newton_Cotes_cerrada(tabla_cerradas[elegido-1], soporte, funcion_uno)

#Calculamos el valor real
I_r = (sp.integrate(funcion_uno(x), (x, soporte[0], soporte[1]))).evalf()
print("La integral exacta es I_r =", I_r)

#Mostramos la aproximación
print("La aproximación de la integral de", funcion_uno(x), "en", soporte, "usando el método de", tabla_cerradas[elegido-1]["Nombre"], "es I =", aproximacion.evalf())

#Calculamos y mostramos el error
error_relativo = calcular_error_relativo(aproximacion, I_r)
print("El error relativo es", error_relativo.evalf())

#Parte opcional
print("\n\n----- PARTE OPCIONAL -----")
soporte = [1,3]
N = int(input("Introduce la N para la fórmula de Milne compuesta:\n"))
print("Ahora vamos a hacer la fórmula de Milne compuesta con N =", N, "para la función f(x) =", funcion_opcional(x))

#Calculamos el valor real
I_r = (sp.integrate(funcion_opcional(x), (x, soporte[0], soporte[1]))).evalf()
print("La integral exacta es I_r =", I_r)

#Calculamos la aproximación
I = formula_compuesta_N(tabla_cerradas[3], soporte, funcion_opcional, N)
print("La aproximación de f(x) con la fórmula de Milne compuesta es I =", I)

#Calculamos y mostramos el
error_relativo = calcular_error_relativo(I, I_r)
print("El error relativo es", error_relativo.evalf())
#Fin de Parte opcional
