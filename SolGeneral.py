import numpy as np 


# condiciones iniciales 

n = 3 # numero de masas (numero entero)
M = [1,2,3] # masa de cada cuerpo

R = np.zeros((n, 3)) #posicion de cada cuerpo
R[0, 2] = 1
R[1, 1] = 1
R[2, 0] =1

V = np.zeros((n, 3)) #velocidad de cada cuerpo
V[0, 0] = 1
V[1, 1] = 1
V[2, 2] = 1



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
        RNew[i+1] = VNew * dt + RNew[i]
        VNew = DV(n, M, RNew[i]) * dt + VNew

    return RNew



REuler = Euler(10,1,R,V)
print(REuler[:, 0, : ]) #posiciones xyz cuerpo 1
print(REuler[:, 1, : ]) #posiciones xyz  cuerpo 2
print(REuler[:, 2, : ]) #posiciones xyz  cuerpo 3



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
        k2v = DV(n, M, RNew[i]+ k2/2)
        k3 = VNew  + k2v * dt /2
        k3v = DV(n, M, RNew[i]+ k3/2)
        k4 = VNew + k3v * dt 
        k4v = DV(n, M, RNew[i]+ k3)

        RNew[i+1] = RNew[i] + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
        VNew =  VNew + dt * (k1v + 2*k2v + 2*k3v + k4v) / 6
    
    return RNew

RRK4 = RK4(10,1,R,V)
print(RRK4[:, 0, : ]) #posiciones xyz  cuerpo 1
print(RRK4[:, 1, : ]) #posiciones xyz  cuerpo 2
print(RRK4[:, 2, : ]) #posiciones xyz  cuerpo 3

