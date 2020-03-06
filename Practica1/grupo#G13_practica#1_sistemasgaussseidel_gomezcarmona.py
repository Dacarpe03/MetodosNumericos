# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 20:09:59 2019

@authors: Lucas Gómez y Daniel Carmona
"""
import numpy as np

def determinante(A):
    return np.linalg.det(A)

def determinante_distinto_de_cero(A):
    return determinante(A)!=0

def es_diagonalmente_dominante(A):
    for i in range(len(A)):
        suma_linea = 0;
        valor_elem_diagonal = abs(A[i][i])
        #Sumamos los valores absolutos de los elementos de la línea menos el elemento de la diagonal
        for j in range(len(A)):
            if(i!=j):
                suma_linea += abs(A[i][j])
        #Si el valor abs del elemento Aii es menor que la suma entonces no es diag. dominante
        if(valor_elem_diagonal <= suma_linea):
            print("No es diagonalmente dominante, puede que la solución no converja\n")
    #Si llegamos aquí es que la matriz es diag. dominante
    print("El sistema es diagonalmente dominante.\n")    


def calcular_error_max(vector_actual, vector_anterior):
    error_max = 0
    for i in range(len(vector_actual)):
        error_i = abs(vector_actual[i]-vector_anterior[i])
        if(error_i > error_max):
            error_max = error_i
    return error_max
    
def metodo_Gauss_Seidel(A, B, vector_solucion, tolerancia):
    #Si el determinante es distinto de cero
    if(determinante_distinto_de_cero(A)):
        es_diagonalmente_dominante(A)
        #Inicializo el vector k-1
        vector_aux = vector_solucion.copy()
        #Error máximo en la iteración k
        error_max_k = calcular_error_max(vector_solucion, vector_aux)
        #Iteración k
        iteracion_k = 0;
        while(abs(error_max_k) > tolerancia or iteracion_k == 0):
            #Actualizamos el vector anterior copiando la solución de la iteración anterior
            vector_aux = vector_solucion.copy()
            #Estamos en la siguiente iteración
            iteracion_k += 1
            #Calculamos el vector de la iteración k
            for i in range(len(vector_solucion)):
                #Primero calculamos los sumatorios
                sumatorio_ant = 0;
                sumatorio_act = 0;
                numerador = 0;
                #Del vector actual
                for j in range(i):
                    sumatorio_act += A[i][j]*vector_solucion[j]
                #Del vector anterior
                for j in range(i+1, len(vector_solucion)):
                    sumatorio_ant += A[i][j]*vector_aux[j]
                #Ahora calculamos el numerador
                numerador = B[i] - sumatorio_act - sumatorio_ant
                #Finalmente calculamos x_i_k
                vector_solucion[i] = numerador/A[i][i]
            #Por último calculamos el error de la iteración
            error_max_k = calcular_error_max(vector_solucion, vector_aux)
            print("Iteración {}: vector_k = {}   Error máximo={}".format(iteracion_k, vector_solucion, error_max_k))
    return vector_solucion;


if __name__ == "__main__":
    #Tolerancia
    tolerancia = 0.01
    #Matriz del sistema
    A = [[3, -0.1, -0.2],
         [0.1, 7, -0.3], 
         [0.3, -0.2, 10]]
    #Vector B del sistema
    B = [7.85, -19.3, 71.4]
    #Vector inicial de soluciones
    v = [0, 0, 0]
    #Calculamos la solución
    v_sol = metodo_Gauss_Seidel(A, B, v, tolerancia)
    