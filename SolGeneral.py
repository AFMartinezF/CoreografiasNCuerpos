import numpy as np 


# funciones EDO segundo orden en poscion de cada masa

def DV(n,M,R): # posicion de un cuerpo respecto a otro 

    DV = np.zeros((n, 3)) # dv(numero de masas, cordenadas xyz) 
    G = 1

    for i in range (n):
        for j in range (n):
            if i != j:
                rij = R[j] - R[i] #posicion de un cuerpo respecto a otro
                DV[i] += G *( M[j] * rij / np.linalg.norm(rij)**3 ) 
    
    return DV




#metodo de euler
 
def Euler(NPasos, tFinal, R, V): # posicion de un cuerpo respecto  matriz 3d

    dt = tFinal / NPasos
    RNew = np.zeros((NPasos+1, n, 3)) # RNew (tienpo ,numero de masas, cordenadas xyz) nueva posicion de cada cuerpo
    RNew[0] = R
    VNew = V

    for i in range(NPasos):
        RNew[i+1] = RNew[i] + VNew * dt 
        VNew = VNew + DV(n, M, RNew[i]) * dt 

    #RNew[:, 0, : ]) #posiciones xyz cuerpo 1
    #RNew[:, 1, : ]) #posiciones xyz  cuerpo 2
    #RNew[:, 2, : ]) #posiciones xyz  cuerpo 3

    #Rsol organiza los datos de los pociones xyz para cada masa horizontalmente 
    RSol = RNew[:, 0, : ]
    for i in range(1,n):
        RSol = np.concatenate((RSol, RNew[:, i, : ]), axis=1)

    return RSol




#metodo de Runge-Kutta

def RK4(NPasos, tFinal, R, V): # posicion de un cuerpo respecto  matriz 3d

    dt = tFinal / NPasos
    RNew = np.zeros((NPasos+1, n, 3)) # RNew (tienpo ,numero de masas, cordenadas xyz) nueva posicion de cada cuerpo
    RNew[0] = R
    VNew = V
   
    for i in range(NPasos):
        
        k1= VNew 
        k1v = DV(n, M, RNew[i]) 
        k2= VNew + k1v * dt /2
        k2v = DV(n, M, RNew[i]+ k1 * dt/2)
        k3 = VNew  + k2v * dt /2
        k3v = DV(n, M, RNew[i]+ k2 * dt/2)
        k4 = VNew + k3v * dt 
        k4v = DV(n, M, RNew[i]+ k3 *dt)

        RNew[i+1] = RNew[i] + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        VNew =  VNew + (k1v + 2*k2v + 2*k3v + k4v) * dt / 6
    
    #RNew[:, 0, : ]) #posiciones xyz cuerpo 1
    #RNew[:, 1, : ]) #posiciones xyz  cuerpo 2
    #RNew[:, 2, : ]) #posiciones xyz  cuerpo 3

    #Rsol organiza los datos de los pociones xyz para cada masa horizontalmente 
    RSol = RNew[:, 0, : ]
    for i in range(1,n):
        RSol = np.concatenate((RSol, RNew[:, i, : ]), axis=1)

    return RSol




#metodo de Leap-Frog

def LF(NPasos, tFinal, R, V): # posicion de un cuerpo respecto  matriz 3d

    dt = tFinal / NPasos
    RNew = np.zeros((NPasos+1, n, 3)) # RNew (tienpo ,numero de masas, cordenadas xyz) nueva posicion de cada cuerpo
    RNew[0] = R
    VNew = V
    VOld = V - dt/2 * DV(n, M, R) # este es V n-1/2

    for i in range(NPasos):
        
        VNew =  VOld +  dt * DV(n, M, RNew[i]) 
        RNew[i+1] = RNew[i] + 2* dt * VNew
        VOld = VNew

    #RNew[:, 0, : ]) #posiciones xyz cuerpo 1
    #RNew[:, 1, : ]) #posiciones xyz  cuerpo 2
    #RNew[:, 2, : ]) #posiciones xyz  cuerpo 3

    #Rsol organiza los datos de los pociones xyz para cada masa horizontalmente 
    RSol = RNew[:, 0, : ]
    for i in range(1,n):
        RSol = np.concatenate((RSol, RNew[:, i, : ]), axis=1)

    return RSol




#metodo de Verlet

def Verlet(NPasos, tFinal, R, V): # posicion de un cuerpo respecto  matriz 3d

    dt = tFinal / NPasos
    RNew = np.zeros((NPasos+1, n, 3)) # RNew (tienpo ,numero de masas, cordenadas xyz) nueva posicion de cada cuerpo
    RNew[0] = R
    ROld = R + V * -dt  # este es R n-1

    for i in range(NPasos):
        
        RNew[i+1] = 2 * RNew[i] - ROld + dt**2 * DV(n, M, RNew[i]) 
    
    #RNew[:, 0, : ]) #posiciones xyz cuerpo 1
    #RNew[:, 1, : ]) #posiciones xyz  cuerpo 2
    #RNew[:, 2, : ]) #posiciones xyz  cuerpo 3

    #Rsol organiza los datos de los pociones xyz para cada masa horizontalmente 
    RSol = RNew[:, 0, : ]
    for i in range(1,n):
        RSol = np.concatenate((RSol, RNew[:, i, : ]), axis=1)

    return RSol







# condiciones iniciales 

n = 3 # numero de masas (numero entero)
M = [1,2,3] # masa de cada cuerpo

R = np.zeros((n, 3)) #posicion de cada cuerpo (xyz en la horizontal)
R[0, 2] = 1
R[1, 1] = 1
R[2, 0] =1

V = np.zeros((n, 3)) #velocidad de cada cuerpo (xyz en la horizontal)
V[0, 0] = 1
V[1, 1] = 1
V[2, 2] = 1


REuler = Euler(100,3,R,V)
np.savetxt("MetodoEuler.cvs ",REuler ,delimiter=",")
