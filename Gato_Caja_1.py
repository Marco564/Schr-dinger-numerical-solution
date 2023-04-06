# El objetivo de este programa es resolver la ecuación de Schrodinger para
# cualquier potencial.



import numpy as np
import time #Esto solo para ver el tiempo de ejecución
import sympy as sp

#Bibliotecas hechas por mi
from Funciones_Dominio import gauss, Dominio
from Graficacion import graficar

#_________________________________________________________________________



#_________________________________________________________________________________
# Básciamente esta función obtiene un arreglo con los índices de los puntos que se
# encuentran en la vecindad del punto i.
def Add_Condition_W(W):
    for i in range( len(W) ):
        W[i].append(0)
    
    # esto es para agregar una diagonal con 1 a la matriz w
    aux = [0 for i in range( len(W[0]) )]
    aux[ len(aux) - 1 ] = 1

    W.append(aux)

    return W
def Indices_Vecindad(i,D,a,sigma):
        V_indices = [] # V de vecindad, o sea los índices de los puntos que están
                    # en la vecindad de i
        for j in range(len(D)):
            peso = gauss(a, sigma, D[i], D[j])
            if round(peso,5) != 0:
                V_indices.append(j)
        
        n = len(V_indices) # numero de puntos en vecindad de i
        return V_indices,n
#_________________________________________________________________________________

#_________________________________________________________________________________
def Build_M(condicion,D,i,V_indices,n,N):
    M = []
    for j in range(n):
        dx = D[ V_indices[j] ][0] - D[i][0]
        dy = D[ V_indices[j] ][1] - D[i][1]
        M.append( [1 , dx, dy, dx**2/2, dx*dy, dy**2/2] )
    M.append( [0, 0, 0, 1, 0, 1] )

    if condicion == 1:
        M =  np.array(M)
        return M
    elif condicion == 2:
        M.append( [1, 0, 0, 0, 0, 0] )
        M =  np.array(M)
        return M
    else:
        M.append( [0, N[i][2], N[i][3], 0, 0, 0] )
        M =  np.array(M)
        return M

    
#_________________________________________________________________________________





#_________________________________________________________________________________
def Build_W(condicion,D,i,a,sigma,V_indices,n):
    W = []
    
    for j in range(n):
        aux = [0 for j in range(n)]
        aux[j] = gauss(a,sigma,D[i],D[ V_indices[j] ])
        W.append(aux)
    W = Add_Condition_W(W)
    
    if condicion == 1:
        W =  np.array(W)
        return W
    else:
        W = Add_Condition_W(W)
        W =  np.array(W)
        return W
    



#_________________________________________________________________________________
def Build_b(condicion,indices,n,f,g,Z,D,i):
    
    b = []
    for j in range(n):
        b.append( Z[ indices[j] ] )
    b.append( float(f.subs({'x': D[i][0],'y':D[i][1]}).evalf())  )
    
    if condicion != 1:
        b.append( float(g.subs({'x': D[i][0],'y':D[i][1]}).evalf()) )

    b = np.array( [b] )
    b = np.transpose(b)
    
    return b
#_________________________________________________________________________________




#_________________________________________________________________________________
def Aproximacion(condicion,D,i,indices,n,N,a,sigma,f,g,Z):
    
    M = Build_M(condicion,D,i,indices,n,N)
    W = Build_W(condicion,D,i,a,sigma,indices,n)
    b = Build_b(condicion,indices,n,f,g,Z,D,i)

    aux_1 = np.matmul( np.transpose(M), W )
    aux_2 = np.linalg.inv(np.matmul( aux_1, M ))
    aux_3 = np.matmul( aux_2, aux_1  )
    
    A = np.matmul( aux_3, b )

    return A[0,0]
#_________________________________________________________________________________




# Este es como el programa principal.
def main():

    FE,a,sigma,NE,D = Dominio(1,25)
    # El dominio es cuadrado y con esquina en el origen

    
    #_________________________________________________________________________
    # En esta parte ya obtengo un resultado que son los z que le corresponde a
    # cada punto en el plano. Contruyo matrices M y W en cada iteración y con
    # las condiciones del problema ya saco un resultado.

    # O sea puede que a algunos puntos en la frontera aplique dirichlet y a otros
    # la neumann. O puede que a todos en la frontera se les aplique o puede que 
    # estén combinados.

    #Haré primero la aprox del paper donde se inicializa en cero todo
    Z = [0 for j in range(len(D))] # Esto irá cambiando conforme se itere el programa.
    # Supondré que al punto Z[i] le corresponde el punto D[i]

    #Puedo agregarle una tercer componente a todos los D[i], que sería un valor 1,2 o 3
    # 1 si es punto interior
    # 2 si es condición dirichlet
    # 3 si es condición neumann
    
    #Pongo condiciones a cada punto
    for j in range( len(D) ):
        if j >= len(FE):
            D[j].append(1)
        else:
            D[j].append(2) # Dirichlet
            #D[j].append(3) #Neumann
    
    '''
    # Aquí es por si combino condiciones, sino quiero combinar solo debo comentar este for
    for j in range( len(D) ):
        if j < len(FE):
            if D[j][0] == 1 or D[j][1]==1:
                D[j][2] = 2
            else:
                D[j][2] = 3
    '''


    # Esto para poder tener, f y g como cualquier función real de x,y
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    
    #Inicia la aproximación
    f = sp.sympify("-2") #Con esto se supone ya son funciones de x,y. -cos(pi*x)
    g = sp.sympify("0") #Siempre se debe elegir o g o s o los dos sino no se podría resolver numéricamente
    s = sp.sympify("0")
    
    sigue = True
    while sigue:
        Zi = Z.copy() # Copio como inicia para al final sacar el error
        for i in range(len(D)):

            indices,n = Indices_Vecindad(i,D,a,sigma)
            
            if D[i][2] == 1:
                Z[i] = Aproximacion(1,D,i,indices,n,'NaN',a,sigma,f,g,Z)

            elif D[i][2] == 2:
                Z[i] = Aproximacion(2,D,i,indices,n,'NaN',a,sigma,f,g,Z)
                    
            else:
                Z[i] = Aproximacion(3,D,i,indices,n,NE,a,sigma,f,s,Z)
        
        # en esta parte checo el error
        sum1 = 0
        sum2 = 0
        for i in range(len(Z)):
            sum1 += abs(Z[i] - Zi[i])
            sum2 += abs(Z[i])
        crit = sum1/sum2
        if crit < 0.1:
            sigue = False
        else:
            sigue = True
                
    
    # En este punto ya tengo la lista de matrices M y la lista de matrices W.
    #_________________________________________________________________________

    #Esto solo para que no cuente el tiempo que yo veo la gráfica
    # o que se toma en graficar.
    end_time = time.time() 

    graficar(D,Z,NE)

    return end_time

if __name__ == '__main__':
    start_time = time.time()
    
    end_time = main()

    print("\nTiempo de ejecución:", end_time - start_time, "segundos")
    print("Tiempo de ejecución:", (end_time - start_time)/60, "minutos")
    print("Tiempo de ejecución:", (end_time - start_time)/3600, "horas\n")
#_________________________________________________________________________________