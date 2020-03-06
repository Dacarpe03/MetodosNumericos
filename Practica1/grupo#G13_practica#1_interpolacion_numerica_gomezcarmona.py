# -*- coding: utf-8 -*-
"""
Editor de Spyder

@authors: Lucas Gómez y Daniel Carmona
"""
import sympy as sp
x = sp.symbols('x')

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

def polinomios_de_base_de_Lagrange(soporte):
    num_puntos = len(soporte)#Calculamos el número de puntos del soporte
    polinomios_base = [1]*num_puntos #Para almacenar todos los Li creamos una lista vacía
    #Empezamos el for para calcular los Li
    for i in range(num_puntos): #De 0 a numPuntos-1
        polinomio_i = 1; #El polinomio Li
        for j in range(num_puntos): #Pitorio para calcular Li
            if(i != j): #Si i es distinto de j
                polinomio_i *= (x-soporte[j])*sp.Rational(1, (soporte[i]-soporte[j])) #Fórmula del pitorio
        #Guardamos el polinomio Li
        polinomios_base[i] = polinomio_i
    
    return polinomios_base #Devolvemos los polinomios
#####FIN polinomios_de_base_de_Lagrange
    
def polinomio_interpolador_de_Lagrange(soporte, funcion):
    #Primero calculamos los polinomios base de los puntos del soporte
    polinomios_base_Lagrange = polinomios_de_base_de_Lagrange(soporte) 
    #Mostrar los polinomios de Lagrange
    for i in range(len(polinomios_base_Lagrange)):
        print("El polinomio base L{} es:".format(i), sp.expand(polinomios_base_Lagrange[i]).evalf())
    #Calcular el polinomio interpolador de Lagrange
    polinomio_Lagrange = 0
    for i in range(len(polinomios_base_Lagrange)):#Fórmula del sumatorio
        polinomio_Lagrange += polinomios_base_Lagrange[i] * funcion(soporte[i])
    return polinomio_Lagrange
#####FIN polinomio_interpolador_de_Lagrange
    
def funcion_a_aproximar(t):
    return sp.exp(-t)+sp.cos(4*t/sp.pi)
#####FIN funcion_a_aproximar


print("Vamos a aproximar la funcion f(x)={}".format(funcion_a_aproximar(x)))
soporte = [0,2]
puntos_soporte = 6;
soporte_con_mas_puntos = calcular_puntos_del_soporte(puntos_soporte-1,soporte)
print("Puntos del soporte:",soporte_con_mas_puntos)
#Calculamos solución
polinomio_solucion = polinomio_interpolador_de_Lagrange(soporte_con_mas_puntos,funcion_a_aproximar)
#Mostramos solución
print("El polinomio interpolador de Lagrange con {} puntos es:".format(len(soporte_con_mas_puntos)),sp.expand(polinomio_solucion.evalf()))
