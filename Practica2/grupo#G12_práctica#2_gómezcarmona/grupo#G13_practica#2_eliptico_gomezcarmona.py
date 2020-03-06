# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 10:10:08 2019

@author: Daniel y Lucas
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def f(xm, yn):
    return -2
#fin f
    
def funcion_contorno(x, y):
    return 0
#fin funcion_contorno

def calcular_N_o_M(intervalo, paso):
    #Calculamos el número de puntos y lo devolvemos
    N = int(1 + (intervalo[1]-intervalo[0])/paso)
    return N
#fin calcular_N_o_M
    
def calcular_soporte(intervalo, paso, N):
    soporte = list()
    #Calculamos desde t_1 hasta t_N
    for i in range(1, N+1):
        t_i = intervalo[0] + (i-1)*paso
        soporte.append(t_i)
    #Devolvemos el soporte
    return soporte
#fin calcular_soporte

def determinar_matriz_coeficientes(M, N, h, k):
    A = []
        
    for i in range(N*M): #Filas
        #print()
        fila = [0]*M*N
        if(i<M or i>((M-1)*(N))-1):#Si estamos en una condición de contorno de y
            fila[i] = 1 #Asignamos 1 en la diagonal
        elif((i+1)%(M)==0 or (i+1)%(M)==1):#Si estamos en una condición de contorno de x
            fila[i]=1 #asignamos 1 en la diagonal
        else: #Si no es una condición de contorno entonces
            #n
            fila[i]   = -2*(k**2+h**2)  #m,n
            fila[i-1] = k**2            #m-1,n
            fila[i+1] = k**2            #m+1,n
            #n-1
            fila[i-M] = h**2            #m,n-1
            #n+1
            fila[i+M] = h**2            #m,n+1
        A.append(fila)
    return A
#fin determinar_matriz_coeficientes
    
def determinar_matriz_B(soporte_espacio, soporte_tiempo, M, N, f, funcion_contorno):
    print("La matriz B:")
    B=[]
    for j in range(N*M):
        if(j<M): #j=1
            v = funcion_contorno(soporte_espacio[j], soporte_tiempo[0])
            print("| {} |".format(v))
            B.append(v)
        elif(j>((M-1)*(N))-1): #j=N
            m = j%(M-1)
            v = funcion_contorno(soporte_espacio[m], soporte_tiempo[N-1])
            print("| {} |".format(v))
            B.append(v)
        elif((j+1)%M==0): #i=M
            n = j%(N-1)
            v = funcion_contorno(soporte_espacio[0], soporte_tiempo[n])
            print("| {} |".format(v))
            B.append(v)
        elif((j+1)%M==1): #i=1
            n = j%(N-1)
            v = funcion_contorno(soporte_espacio[M-1], soporte_tiempo[n])
            print("| {} |".format(v))
            B.append(v)
        else: #Los demás intermedios
            n = j%(N-1)
            m = j%(M-1)
            v = (h**2)*(k**2)*f(soporte_espacio[m], soporte_tiempo[n])
            print("| {} |".format(v))
            B.append(v)
    print()
    return B
#fin determinar_matriz_B
    
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
            print("\nA no es diagonalmente dominante, puede que la solución no converja.\n")
            return
    #Si llegamos aquí es que la matriz es diag. dominante
    print("El sistema es diagonalmente dominante.\n")
    
def determinante(A):
    return np.linalg.det(np.array(A))

def determinante_distinto_de_cero(A):
    return determinante(A)!=0

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
                sumatorio_ant = 0
                sumatorio_act = 0
                numerador = 0
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
            #print("Iteración {}: vector_k = {}   Error máximo={}".format(iteracion_k, vector_solucion, error_max_k))
    return vector_solucion

def mostrar_matriz(A):
    print("La matriz A:")
    for row in A:
        print ('|%s|' % (' '.join('%5s' % round(i,4) for i in row)))
#fin mostrar_matriz      
        
def mostrar_salida(lista):
    Tabla = """\
+---------------------------------------------------------------------------------------------------------------+
| n               m               tn                        xm                        u_mn                      |
|---------------------------------------------------------------------------------------------------------------|
{}
+---------------------------------------------------------------------------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<15} {:<15} {:<25} {:<25} {:<25} |".format(*fila)
    for fila in lista)))
    print(Tabla)
#fin mostrar_salida
    
def crear_tabla(sol_buena, soporte_espacio, soporte_tiempo, M, N):
    vector_tabla = list()
    vector_espacio = list()
    vector_tiempo = list()
    vector_soluciones = list()
    
    i=0
    for n in range(N):
        for m in range(M):
            vector_tabla.append([n, m, soporte_tiempo[n], soporte_espacio[m], sol_buena[i]])
            vector_espacio.append(soporte_espacio[m])
            vector_tiempo.append(soporte_tiempo[n])
            vector_soluciones.append(sol_buena[i])
            i+=1     
    mostrar_salida(vector_tabla)
    
    vector_grafica = [vector_espacio, vector_tiempo, vector_soluciones]
    
    return vector_grafica

def pintar_grafica3d(vector_a_pintar):
    espacio  = vector_a_pintar[0]
    tiempo   = vector_a_pintar[1]
    solucion = vector_a_pintar[2]
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    ax.scatter3D(espacio, tiempo, solucion, c=solucion);
    ax.set_xlabel('Espacio')
    ax.set_ylabel('Tiempo')
#fin pintar_grafica3d
    
if __name__ == "__main__":
    #x
    intervalo_espacio = [-1,1]
    h = 0.4
    
    #y
    intervalo_tiempo = [-1,1]
    k = 0.1
    
    tolerancia = 0.00000001
    M = calcular_N_o_M(intervalo_espacio, h)
    N = calcular_N_o_M(intervalo_tiempo, k)
    soporte_x = calcular_soporte(intervalo_espacio, h, M)
    soporte_y = calcular_soporte(intervalo_tiempo, k, N)
    
    A = determinar_matriz_coeficientes(M, N, h, k)
    B = determinar_matriz_B(soporte_x, soporte_y, M, N, f, funcion_contorno)
    mostrar_matriz(A)
    
    vector_inicial = [0]*N*M
    v_solucion = metodo_Gauss_Seidel(A, B, vector_inicial, tolerancia)
    
    v_pintar = crear_tabla(v_solucion, soporte_x, soporte_y, M, N)
    pintar_grafica3d(v_pintar)