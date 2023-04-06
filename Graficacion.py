import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Función para darle el formato para graficar con matplotlib y también servirá para la 
# del círculo y en general todo lo que venga con el formato [[x1,y1],[x2,y2],...,[xn,yn]]
def form_graph(D):
    X = []
    Y = []
    Z = []
    for i in range(len(D)):
        X.append(D[i][0])
        Y.append(D[i][1])
        Z.append(D[i][2])
    return(X,Y,Z)
def graficar(D,Z,NE):
    #_________________________________________________________________________
    # Aquí solo le doy forma a las gráficas y las ploteo y muestro
    P1 = [] #puntos interiores
    P2 = [] #puntos dirichlet
    P3 = [] #puntos neumann
    for i in range(len(D)):
        if D[i][2] == 1:
            P1.append( [D[i][0],D[i][1],Z[i]] )
        elif D[i][2] == 2:
            P2.append( [D[i][0],D[i][1],Z[i]] )
        else:
            P3.append( [D[i][0],D[i][1],Z[i]] )

    X1,Y1,Z1 = form_graph(P1)
    X2,Y2,Z2 = form_graph(P2)
    X3,Y3,Z3 = form_graph(P3)

    #Creo la figura
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    #grafico normales
    for i in range(len(NE)):
        ax.quiver(NE[i][0], NE[i][1], Z[i], NE[i][2], NE[i][3], 0, length = 0.1, normalize=True)

    # Graficar los puntos en 3D
    ax.scatter(X1, Y1, Z1, c='black', marker='o', label = 'Interiores')
    ax.scatter(X2, Y2, Z2, c='green', marker='o', label = 'Dirichlet')
    ax.scatter(X3, Y3, Z3, c='blue', marker='o', label = 'Neumann')


    # Añadir etiquetas de los ejes
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')

    ax.legend()

    plt.show()
    #_________________________________________________________________________

    return 