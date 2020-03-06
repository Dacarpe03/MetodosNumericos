# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:42:16 2019

@author: Lucas y Daniel
"""
import numpy as np
import matplotlib.pyplot as plt

def f(t, y):
    return -2*t**3 + 12*t**2 - 20*t + 8.5

def f_real(t):
    return -0.5*t**4 + 4*t**3 - 10*t**2 + 8.5*t + 1

def calcular_error_relativo(t, y):
    #Calculamos el error relativo
    error = abs(y - f_real(t))
    error_relativo = 100*error/f_real(t)
    return error_relativo
#fin calcular_error
    
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

def RK2(diccionario, intervalo, h, solucion_inicial):
    #Calculamos el número de puntos del soporte
    N = calcular_N(intervalo, h)
    #Discretizamos el intervalo
    soporte = calcular_soporte(intervalo, h, N)
    #Ahora cogemos los coeficientes del diccionario o tabla
    a1 = diccionario["a1"]
    a2 = diccionario["a2"]
    p1 = diccionario["p1"]
    q11 = diccionario["q11"]
    #Creamos una lista de soluciones para cada iteracion
    soluciones = list()
    #Ahora hacemos la aproximación
    soluciones.append([1, soporte[0], solucion_inicial, solucion_inicial, 0])
    for i in range(1,N):
        #Cogemos los datos necesarios para calcular y_n+1
        t_n = soluciones[i-1][1]
        y_n = soluciones[i-1][3]
        #Calculamos las k
        k1 = f(t_n, y_n)
        k2 = f(t_n + p1*h, y_n + q11*k1*h)
        #Calculamos y_n+1
        y_n1 = y_n + (a1*k1 + a2*k2)*h
        #Calculamos el valor real
        y_x = f_real(soporte[i])
        #Calculamos el error
        error_rel = calcular_error_relativo(soporte[i], y_n1)
        #Guardamos la solución de la iteración
        soluciones.append([i+1, soporte[i], y_x, y_n1, error_rel])
    return soluciones

def RK3(diccionario, intervalo, h, solucion_inicial):
    #Calculamos el número de puntos del soporte
    N = calcular_N(intervalo, h)
    #Discretizamos el intervalo
    soporte = calcular_soporte(intervalo, h, N)
    #Ahora cogemos los coeficientes del diccionario o tabla
    a1 = diccionario["a1"]
    a2 = diccionario["a2"]
    a3 = diccionario["a3"]
    p1 = diccionario["p1"]
    p2 = diccionario["p2"]
    q11 = diccionario["q11"]
    q21 = diccionario["q21"]
    q22 = diccionario["q22"]
    #Creamos una lista de soluciones para cada iteracion
    soluciones = list()
    #Ahora hacemos la aproximación
    soluciones.append([1, soporte[0], solucion_inicial, solucion_inicial, 0])
    for i in range(1,N):
        #Cogemos los datos necesarios para calcular y_n+1
        t_n = soluciones[i-1][1]
        y_n = soluciones[i-1][3]
        #Calculamos las k
        k1 = f(t_n, y_n)
        k2 = f(t_n + p1*h, y_n + q11*k1*h)
        k3 = f(t_n + p2*h, y_n + q21*k1*h + q22*k2*h)
        #Calculamos y_n+1
        y_n1 = y_n + (a1*k1 + a2*k2 + a3*k3)*h
        #Calculamos el valor real
        y_x = f_real(soporte[i])
        #Calculamos el error
        error_rel = calcular_error_relativo(soporte[i], y_n1)
        #Guardamos la solución de la iteración
        soluciones.append([i+1, soporte[i], y_x, y_n1, error_rel])
    return soluciones
#fin_RK3
    
def RK4(diccionario, intervalo, h, solucion_inicial):
    #Calculamos el número de puntos del soporte
    N = calcular_N(intervalo, h)
    #Discretizamos el intervalo
    soporte = calcular_soporte(intervalo, h, N)
    #Ahora cogemos los coeficientes del diccionario o tabla
    a1 = diccionario["a1"]
    a2 = diccionario["a2"]
    a3 = diccionario["a3"]
    a4 = diccionario["a4"]
    p1 = diccionario["p1"]
    p2 = diccionario["p2"]
    p3 = diccionario["p3"]
    q11 = diccionario["q11"]
    q21 = diccionario["q21"]
    q31 = diccionario["q31"]
    q22 = diccionario["q22"]
    q32 = diccionario["q32"]
    q33 = diccionario["q33"]
    #Creamos una lista de soluciones para cada iteracion
    soluciones = list()
    #Ahora hacemos la aproximación
    soluciones.append([1, soporte[0], solucion_inicial, solucion_inicial, 0])
    for i in range(1,N):
        #Cogemos los datos necesarios para calcular y_n+1
        t_n = soluciones[i-1][1]
        y_n = soluciones[i-1][3]
        #Calculamos las k
        k1 = f(t_n, y_n)
        k2 = f(t_n + p1*h, y_n + q11*k1*h)
        k3 = f(t_n + p2*h, y_n + q21*k1*h + q22*k2*h)
        k4 = f(t_n + p3*h, y_n + q31*k1*h + q32*k2*h + q33*k3*h)
        #Calculamos y_n+1
        y_n1 = y_n + (a1*k1 + a2*k2 + a3*k3 + a4*k4)*h
        #Calculamos el valor real
        y_x = f_real(soporte[i])
        #Calculamos el error
        error_rel = calcular_error_relativo(soporte[i], y_n1)
        #Guardamos la solución de la iteración
        soluciones.append([i+1, soporte[i], y_x, y_n1, error_rel])
    return soluciones
#fin_RK4
    
def mostrar_salida(lista):
    Tabla = """\
+----------------------------------------------------------------------------------------------------+
| n               t_n                  f_real               y_n                  error relativo (%) 
|----------------------------------------------------------------------------------------------------|
{}
+----------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<15} {:<20} {:<20} {:<20} {:<20}|".format(*fila)
    for fila in lista)))
    print(Tabla)
#fin mostrar_salida
    
def pintar_funciones(soluciones, soporte ,f_real, orden, nombre):
        y = list()
        for s in soluciones:
            y.append(s[3])

        t = np.arange(soporte[0], soporte[len(soporte)-1], 0.01)
        plt.plot()
        plt.plot(soporte ,y , 'bs', label='Runge-Kutta orden {} con método {}'.format(orden, nombre))
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
    
    #Defino los distintos métodos para Runge-Kutta de orden 2 fijando a2 y obteniendo
    #a1=1-a2; p1 = 1/2a2; q11=1/2a2
    diccionario_RK2_heun =     {"Nombre": "Heun",        "a1": 1/2, "a2": 1/2, "p1": 1,   "q11": 1}
    diccionario_RK2_ptomedio = {"Nombre": "Punto medio", "a1": 0,   "a2": 1,   "p1": 1/2, "q11": 1/2}
    diccionario_RK2_ralston =  {"Nombre": "Ralston",     "a1": 1/3, "a2": 2/3, "p1": 3/4, "q11": 3/4}
    
    #Defino los coeficientes para RK3
    diccionario_RK3 = {"Nombre": "común", "a1": 1/6, "a2": 4/6, "a3": 1/6, "p1": 1/2, "p2":1, "q11": 1/2, "q21":-1, "q22":2}
    
    #Defino los coeficientes para RK4
    diccionario_RK4 = {"Nombre": "común", "a1": 1/6, "a2": 1/3, "a3": 1/3, "a4":1/6 ,"p1": 1/2, "p2":1/2, "p3":1, "q11": 1/2, "q21":0, "q31": 0, "q22":1/2, "q32":0, "q33":1}
    
    tabla_orden2 = [diccionario_RK2_heun, diccionario_RK2_ptomedio, diccionario_RK2_ralston]
    #Parámetros de discretización y solución inicial
    intervalo = [0,1]
    pasoh = 0.5
    solucion_inicial = 1
    
    N=calcular_N(intervalo, pasoh)
    soporte = calcular_soporte(intervalo, pasoh, N)
    
    print("Intervalo {} con paso {}".format(intervalo, pasoh))
    
    for i in range(3):
        print("Runge-Kutta de orden 2 con el método", tabla_orden2[i]["Nombre"])
        soluciones = RK2(tabla_orden2[i], intervalo, pasoh, solucion_inicial)
        mostrar_salida(soluciones)
        pintar_funciones(soluciones, soporte, f_real, 2, tabla_orden2[i]["Nombre"])
    
    print("Runge-Kutta de orden 3 con el método común")
    soluciones = RK3(diccionario_RK3, intervalo, pasoh, solucion_inicial)
    mostrar_salida(soluciones)
    pintar_funciones(soluciones, soporte, f_real, 3, diccionario_RK3["Nombre"])
    
    print("Runge-Kutta de orden 4 con el método común")
    soluciones = RK4(diccionario_RK4, intervalo, pasoh, solucion_inicial)
    mostrar_salida(soluciones)
    pintar_funciones(soluciones, soporte, f_real, 4, diccionario_RK4["Nombre"])
    
    