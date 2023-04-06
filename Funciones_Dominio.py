import shapely as sh
import math as mh
from shapely import Point, Polygon
import random as rd
import numpy as np

def NormalesProm(F,where_n):
    # En esta función saco las normales pero por el método de promedios
    # para este programa lo usaré más que nada para mejorar lo de las 
    # normales por optimización.
    poligono = Polygon(F)
    
    # los vectores normales constarán de dos pares de puntos, digamos entonces que la estructura de los vectores normales a las 
    # rectas que se encuentran entre cada par de puntos es n = [[xi,yi,xf,yf],[xi,yi,xf,yf],...,[xi,yi,xf,yf]]. Otra cosa importante
    # a notar es que este proceso solo incolucra a las fronteras interiores y exteriores, de momento solo lo haré para la exterior
    # pero se podrá generalizar para las interiores también tal vez con unas modificaciones.
    n = []

    for i in range(len(F)):
        aux = [0,0,0,0]
        
        f = i + 1
        if f == len(F):
            f = 0
        
        aux[0] = (F[i][0] + F[f][0])/2
        aux[1] = (F[i][1] + F[f][1])/2

        dx = F[f][0] - F[i][0]
        dy = F[f][1] - F[i][1]

        if dy == 0:
            th = mh.pi/2
        else:
            m = -dx/dy
            th = mh.atan(m)
            
        aux[2] = aux[0] + mh.cos( th )
        aux[3] = aux[1] + mh.sin( th )

        #Aquí lo divido entre mil para asegurar que el punto esté o adentro o afuera y no curce la figura
        if sh.contains(poligono, Point(aux[0] + mh.cos( th )/1000, aux[1] + mh.sin( th )/1000)) == True and (where_n=='out'):
            aux[2] = aux[0] + mh.cos( th + mh.pi )
            aux[3] = aux[1] + mh.sin( th + mh.pi )

        elif sh.contains(poligono, Point(aux[0] + mh.cos( th )/1000, aux[1] + mh.sin( th )/1000)) == False and (where_n=='in'):
            aux[2] = aux[0] + mh.cos( th + mh.pi )
            aux[3] = aux[1] + mh.sin( th + mh.pi )

        n.append(aux)

    # Ahora voy a sacar los buenos, estos son los que salen de los puntos y son normales a la superficie, esto solo es sacar los
    # promedios de cada par de vectores que tengo.

    aux = [] # Aqui guardare los buenos, cada dato en la frontera debe tener uno, y el formato será el mismo[[x1,y1,x2,y2]...]
        # O sea donde empieza y donde termina el vector.
    for i in range(len(n)):
        
        k = i + 1
        if k >= len(n):
            k = 0

        promx = (n[i][2] + n[k][2])/2
        promy = (n[i][3] + n[k][3])/2
        
        aux.append([F[k][0], F[k][1], promx - F[k][0], promy - F[k][1]])

    #Lo único es que está recorrido N o sea empieza en el punto 1 y no en el 0 en referencia al FE, entonces solo los recorro
    # para que esté de acuerdo a lo más conveniente.

    N = [] # Aquí guardaré ya el producto final de todo este rollo, osea los vectores normales de interés
    for i in range(len(aux)):
        k = i - 1
        if k == -1:
            k = len(aux) - 1
        N.append(aux[k])

    return N

# Función para sacar la distancia entre dos puntos, es para no repetir tanto la fórmula
def dist(P1,P2):
    # Recibo dos arreglos tipo [x1,y1] y [x2,y2]
    d = ( (P2[0] - P1[0])**2 + (P2[1] - P1[1])**2 )**0.5
    return d
# Esta función básicamente solo es para generar un punto y revisar que si cumpla que este a una distancia mayor a r de
# todos los demás puntos.

def generar(x,y, FE, poligono, r, PM):

    #Ahora a este punto le genero otro punto a una distancia r, pero a un ángulo random y tiene que cumplir que no puede
    # atravesar el radio de otro punto, si es así no pongo nada y comienzo de nuevo.
    ran = rd.random()
    Pg = [x + r*mh.cos(2*mh.pi*ran), y + r*mh.sin(2*mh.pi*ran)]
    
    # Pasa me dirá si es válido o no el punto
    pasa = True
    for i in range(len(FE)):
        if dist(Pg,FE[i]) <= r:
            pasa = False
            break

    if pasa == True:
        for i in range(len(PM)):
            if dist(Pg,PM[i]) <= r:
                pasa = False
                break
    
    if pasa == True:
        if sh.contains(poligono, Point(Pg[0],Pg[1])) == False:
            pasa = False
    
    # Ya terminé de checarlo con todos los demás puntos ahora si lo apendeo si es válido sino que empiece de nuevo el ciclo
    if pasa == True:
        PM.append(Pg)
        return pasa

    elif pasa == False:
        return pasa


def gauss(a,h,r1,r2):
    #h = sigma de la distribución
    if dist(r2,r1)/h <= 1:
        f = mh.exp(-a * ( dist(r2,r1)**2 )/h**2  )
    else:
        f = 0
    return f


#Normales debe recibir F = [[x1,y1],[x2,y2],...,[xn,yn]] osea un arreglo que
# que tiene los puntos de un contorno de esa forma. Y where_n='in' o 'out'
# un string que me dicta si quiero los vectores apuntando hacia adentro
# o hacia afuera, in para contornos interiores y out para el exterior.

# Y tiene que regresar N = [[xi1,yi1,xf1,yf1],[xi2,yi2,xf2,yf2],...,[xin,yin,xfn,yfn]]
# Osea un arreglo que tiene el punto donde inicia el vector y donde termina. En general
# inicia donde están los puntos del contorno y debe apuntar normal al contorno.
def NormalesOp(F,where_n, r, a, sigma):

    N_P = NormalesProm(F,where_n)

    n = 3 #esta es para tomar n puntos a la derchea y a la izquierda de i.

    N = []
    # Comencemos con i = 0.
    # La función gaussiana es f(x) = -a*exp[(x-b)^2/(2*c^2)], b me dicta alrededor de que punto
    # a es algo arbitrario es como cuanto se eleva la campana, y c es la desviación standard.
    for i in range(len(F)):
        M = [[0,0],[0,0]]  #Esta matriz es como el cero en la sumatoria, aquí se irán acumulando los valores

        for j in range(len(F)):
            #Con esta condición tomo n puntos a la dercha y n a la izquierda de i y que j sea diferente de i no lo
            # comparo consigo mismo.
            if j != i and j<=i+n and j>=i-n:
                rji = np.array( [[F[j][0] - F[i][0]], [F[j][1] - F[i][1]]] )
                rji_T = rji.transpose()
                k = gauss(a, sigma, F[i], F[j])/dist(F[j], F[i])**2 #Digamos un a = 1 y sigma = 1, si elijo sigma=r se rompe.

                M += np.dot(rji , rji_T) * k
        
        #Una vez terminada la sumatoria de matrices ya se puede plantear el problema de eigenvalores y eigenvectores
        B,V = np.linalg.eig(M) #B: lista de eigenvalores. V: lista de eigenvectores normalizados
        
        # En esta parte solo renombro los eigenvectores
        EV1 = np.array( [V[0,0], V[1,0]] )
        EV2 = np.array( [V[0,1], V[1,1]] )

        # Aquí tomo el Vector Comparación que viene de las normales promediadas, y voy a tomarla como
        # una referencia
        VC = np.array( [N_P[i][2], N_P[i][3]] )

        # EV1, EV2 y VC se supone ya están normalizados, entonces el producto punto sería el coseno del
        # ángulo entre ellos, por tanto, d1 o d2 = 1 - abs(cos), al compararlas la que sea más pequeña
        # será la que más cercana estará al promedio y ese eigenvector es el que tomo.
        d1 = 1 - abs( np.dot(EV1,VC) )
        d2 = 1 - abs( np.dot(EV2,VC) )

        if d1 < d2:
            aux = [F[i][0], F[i][1], V[0,0], V[1,0]]
        else: 
            aux = [F[i][0], F[i][1], V[0,1], V[1,1]]
        
        #Voy a poner una condición para específicamente el cuadrado que sería que si me encuentro en
        # una esquina que tome la promedio.
        V1 = np.array( [F[i-1][0] - F[i][0] , F[i-1][1] - F[i][1]] )
        k = i + 1
        if k == len(F):
            k = 0
        V2 = np.array( [F[k][0] - F[i][0] , F[k][1] - F[i][1]] )
        
        if np.dot(V1,V2) == 0:
            aux = [F[i][0], F[i][1], N_P[i][2], N_P[i][3]]


        # pero debo checar aun así si queda adentro, sino lo debo girar 180°.
        #Aquí lo divido entre mil para asegurar que el punto esté o adentro o afuera y no cruce la figura
        poligono = Polygon(F)
        k = 1000
        if where_n == 'out':
            if sh.contains(poligono, Point(aux[0] + aux[2]/k, aux[1] + aux[3]/k)) == True:
                aux[2] = -aux[2]
                aux[3] = -aux[3]

        elif where_n=='in':
            if sh.contains(poligono, Point(aux[0] + aux[2]/k, aux[1] + aux[3]/k)) == False:
                aux[2] = -aux[2]
                aux[3] = -aux[3]
        
        N.append( aux )
    return N

def ContornoCuadradoE(l, d):
    
    #l: de que tamaño quiero el cuadro
    #d: en cuanto divido los segmentos, osea habría l*d puntos en cada lado del cuadrado

    FE = [] #Frontera exterior supondré que tiene el formato [[x1,y1],[x2,y2],...,[xn,yn]], y que van en orden,
        #osea podría yo dibujarlo si sigo el orden de los puntos.
    
    h = l/d

    #Este ciclo es para la base del cuadrado
    i = 0
    while i < l:
        FE.append([i,0])
        i += h
    
    #Este ciclo es para el lado derecho del cuadrado
    i = 0
    while i < l:
        FE.append([l,i])
        i += h
    
    #Este es para la tapa del cuadrado
    i = 0
    while i < l:
        FE.append([l-i,l])
        i+=h

    i = 0
    #Este es para el lado derecho del cuadrado
    while i < l:
        FE.append([0,l-i])
        i+=h

    return FE


def Dominio(a,b):
        #_________________________________________________________________________
    # En esta parte voy a construir el contorno. Para propósitos prácticos en
    # este programa consturiré un contorno exterior cuadrado.
    FE = ContornoCuadradoE(a,b)#a: de que tamaño quiero el caudrado. b: en cuanto divido sus lados
    
    #_________________________________________________________________________
    #Aquí saco la distancia máxima de separación en base a la distancia entre
    # los puntos de las fronteras interiores y exteriores. 
    # Esto será como una distancia de máxima proximidad.
    r = 0
    for i in range(len(FE)):
        r += dist(FE[i],FE[i - 1])
    r /= len(FE)

    # En sí saqué el promedio

    # r va a ser mi distancia máxima entre puntos no pueden estar más cerca que esto.
    #_________________________________________________________________________

    import matplotlib.pyplot as plt

    #_________________________________________________________________________
    PM = []#Puntos de la malla

    poligono = Polygon(FE)

    intentos = 50

    for i in range(len(FE)):
        for j in range(intentos):
            pasa = generar(FE[i][0],FE[i][1], FE, poligono, r, PM)
    i = 0
    while i < len(PM):
        for j in range(intentos):
            pasa = generar(PM[i][0],PM[i][1], FE, poligono, r, PM)
        i += 1
    #_________________________________________________________________________
    

    
    #_________________________________________________________________________
    # Ahora voy a hacer los vectores normales.
    a = 1
    sigma = r*2**2

    NE = NormalesOp(FE,'out', r, a, sigma)
    #_________________________________________________________________________




    #_________________________________________________________________________
    # Aquí solo junto todos los puntos de mi dominio en un solo arreglo
    D = [] #Aquí voy a guardar cada punto del dominio D[i][0] = componente x del 
    # punto i, el 1 es para la y.

    for i in range(len(FE)): #Notar que los primer de D son la frontera
        D.append(FE[i])
    for i in range(len(PM)):
        D.append(PM[i])

    
    return FE,a,sigma,NE,D