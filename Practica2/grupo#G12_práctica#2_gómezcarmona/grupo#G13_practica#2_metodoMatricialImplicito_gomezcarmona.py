# -*- coding: utf-8 -*-
"""


@author: Daniel y Lucas
"""
import sympy as sp
import mpmath as mp
import math as m
import matplotlib.pyplot as plt
import numpy as np

def string_f_real():
    return 'e**t-1/e-1'

def f_real(t):
    u = (np.exp(t)-1)/(np.exp(1)-1)
    return u
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
    error_relativo = 100*error/f_real(t)
    return error_relativo
#fin calcular_error

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


def calcular_alfa(D,h,V):
    x = (-D/h**2)-(V/(2*h))
    return x
#fin calcular_alfa

def calcular_beta(D,h,q):
    x = ((2*D)/(h**2))+q
    return x

#fin calcular_beta

def calcular_gamma(D,h,V):
    x = (-D/h**2)+(V/(2*h))
    return x

#fin gamma

def determinar_A(alfa,beta,gamma,N):
    A = []
    for i in range(N):#recorremos la matriz desde 0 a N-1
        Fila = [0]*N #matriz de ceros para cada posicion de las filas
        A.append(Fila)
    for i in range(N):

        #print()
        for j in range(N):

            if(i==0):#primera fila matriz
                if(j==0): #primera columna matriz

                    A[i][j]=1 #asignamos uno
                else:
                    A[i][j]=0 #asignamos cero
            elif(i == N-1):#ultima fila matriz
                if(j == N-1):#ultima columna matriz
                    A[i][j]= 1 #asignamos un 1
                else:
                    A[i][j] = 0 #asignamos un 0
            else:
                if(i==j):#diagonal
                    A[i][j] = beta #asignamos beta
                elif(j==i-1): #antes de la diagonal
                    A[i][j] = alfa #asignamos alfa
                elif (j == i+1): #después de la diagonal
                    A[i][j]= gamma #asignamos gamma
                else:#resto de casos
                    A[i][j]=0 #asignamos cero
            #print(A[i][j],end='  ')

    return A
#FinDeterminarA

def determinar_B(u0, uL, N, Ta, k):
    B = []
    for i in range(N):
        Fila = [0]*N #matriz de ceros para cada posicion de las filas
        B.append(Fila)
    for i in range(N):
        if(i== 0):#primera fila matriz
            B[i] = u0 #asignamos u0
        elif(i == N-1):#ultima fila matriz
            B[i] = uL #asignamos uL
        else:
            B[i] = k*Ta #asignamos 0
        print("| {} |".format(B[i]))

    return B
#FinDeterminarB

def pintar_funciones(soluciones, soporte , h, tolerancia, L, Ta, k):
    t = np.arange(soporte[0], soporte[len(soporte)-1], 0.01)
    plt.plot()
    plt.plot(soporte ,soluciones , 'bs', label='Método Implícito Matricial')

    plt.title('Comparación')

    # Add X and y Label
    plt.xlabel('Metros')
    plt.ylabel('Temperatura Cº')

    # Add a grid
    plt.grid(alpha=.4, linestyle='--')

    plt.text(0, 180, 'h={}'.format(h))
    plt.text(0, 170, 'tolerancia Gauss-Seidel={}'.format(tolerancia))
    plt.text(0, 160, 'L={} m'.format(L))
    plt.text(0, 150, 'Ta={} Cº'.format(Ta))
    plt.text(0, 140, 'k={} m\u00b2'.format(k))




    plt.legend()
    plt.show()
#fin pintar_funciones


def crear_tabla(N, soporte, v_solucion):
    tabla = list()
    for i in range(N):
        linea = [i+1, soporte[i], v_solucion[i]]
        tabla.append(linea)
    return tabla

def mostrar_matriz(A):
    for row in A:
        print ('|%s|' % (' '.join('%7s' % round(i,4) for i in row)))

def mostrar_salida(lista):
    Tabla = """\
+---------------------------------------------------+
| n               t_n        y_n                  
|---------------------------------------------------|
{}
+------------------------------------------+\
"""
    Tabla = (Tabla.format('\n'.join("| {:<15} {:<10} {:<25}|".format(*fila)
                                    for fila in lista)))
    print(Tabla)
#fin mostrar_salida

if __name__ == "__main__":
    
    toleranciaGS = 0.0000001
    L = 10
    intervalo = [0,L]
    h = 0.1
    
    u_0 = 40 #T0
    u_L = 200 #TL
    
    k = 1
    Ta = 20
    D = 1
    V = 0
    q = k
    
    print("La función real es: ", string_f_real())
    N = calcular_N(intervalo, h)
    soporte = calcular_soporte(intervalo, h, N)

    
    alfa = calcular_alfa(D,h,V)
    beta = calcular_beta(D,h,q)
    gamma = calcular_gamma(D,h,V)
    print("Primero calculamos los coeficientes:")
    print("{} = {} \u2248 {}".format(chr(945), alfa, round(alfa,4)))
    print("{} = {} \u2248 {}".format(chr(946), beta, round(beta,4)))
    print("{} = {} \u2248 {}".format(chr(947), gamma, round(gamma,4)))

    print("\nAhora calculamos la Matriz A:")
    A = determinar_A(alfa,beta,gamma,N)
    mostrar_matriz(A)
    print("\nAhora calculamos la Matriz B:")
    B = determinar_B(u_0,u_L,N,Ta, k)
    
    v_inicial = [0]*N
    v_solucion = metodo_Gauss_Seidel(A,B,v_inicial, toleranciaGS)
    print("\nAhora aplicamos Gauss-Seidel y obtenemos el vector solución {}".format(v_solucion))
    
    tabla = crear_tabla(N, soporte, v_solucion)
    mostrar_salida(tabla)
    pintar_funciones(v_solucion, soporte, h, toleranciaGS, L, Ta, k)



