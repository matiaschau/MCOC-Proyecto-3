# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 22:26:49 2018

@author: Tomás
"""

import scipy as sp

def AlturaNormal(B,ss,yini,C0,S0,n):
    yn = yini
    q = 30 #m**3/s
    tol = 10**(-6)+1
    while tol>(10**(-6)):
        A = B*yn + ss*(yn**2)
        P = B + 2*yn*sp.sqrt(ss**2+1)
        dAdyn = B + 2*ss*yn
        dPdyn = 2*sp.sqrt(ss**2+1)
        Q = (C0*(A**(5./3.)*(P**(-2./3.))*sp.sqrt(S0))/n)

        fyn = (A**(5./3.))*(P**(-2./3.))-(n*q)/(C0*sp.sqrt(S0))
        dfyn = (5./3.)*(A**(2./3.))*(P**(-2./3.))*dAdyn - (2./3.)*(A**(5./3.))*(P**(-5./3.))*dPdyn
        yn2 = yn - (fyn/dfyn)

        tol = abs(yn2-yn)

        yn = yn2
    return yn


def geom(y, b, ss):
    A = b * y + ss * (y**2)
    dAdy = b + 2*y*ss
    l = y * sqrt((ss**2) + 1)
    P = b + 2*l
    dPdy = 2 * sqrt((ss**2) + 1)
    R = A/P
    dRdy = (1/P) * dPdy - (A/(P**2)) * dPdy
    
    valores = [A, dAdy, P, dPdy, R, dRdy]
    
    return valores
    
    
    
def NR_FGV(yini, b, Q, n, So, ss, Co, dx, y1, z1):
    y2 = yini
    
    if Co == 1:
        g = 9.81
    else:
        g = 32.2
    
    tol = 10
    iter = 0
    
    while tol > 0.00001 and iter < 100:
        iter = iter + 1
        valores1 = geom(y1, b, ss)
        valores2 = geom(y2, b, ss)
        
        z2 = z1 - So*dx
        h1 = z1 + y1 + ((Q**2)/(2*g*(valores1[0]**2)))
        h1 = z2 + y2 + ((Q**2)/(2*g*(valores2[0]**2)))
        sf1 = ((Q*n)**2) / ( (Co * valores1[0] * ((valores1[4])**(2/3))) **2 )
        sf2 = ((Q*n)**2) / ( (Co * valores2[0] * ((valores2[4])**(2/3))) **2 )    

        fy2 = ((5/3)*(valores2[0] ** (2/3)) * (valores2[2] ** (- 2/3)) * valores2[1]) - ((2/3)*(valores2[0] ** (5/3)) * (valores2[2] ** (- 5/3)) * valores2[3])
        dfy2dy2 = 1 - ((Q**2)/(g*(valores2[0] ** 3)) * valores2[1]) - dx*(((sf2/valores2[0]) * valores2[1]) + (2/3)*(sf2/valores2[4])*valores2[5])
        
        yf = y2 - (fy2 / dfy2dy2)
        tol = abs(yf - y2)
        y2 = yf
    
    valores = [y2, z2]    
    return valores
        

def FGV(Q, b, So, ss, Co, x1, y1, z1, dx):
    y_ini = 1.
    x = [x1]
    y = [y1]
    z = [z1]
    
    if Co == 1:
        g = 9.81
    else:
        g = 32.2
        
    yn = AlturaNormal(b,ss,yini,C0,S0,n)
    
    diff = 0.001
    dist = 0
    i = -1 #valor de entrada
    while y[i+1] > (yn - diff):
        i = i + 1
        dist = dist + dx
        x.append(x[i]+dx)
        hola = NR_FGV(y_ini, b, Q, n, So, ss, Co, dx, y[i], z[i])
        y.append(hola[0])
        z.append(hola[1])
    
    elsa = hola[0] + hola[1]     #elevacion superficie de agua = y + z
    
    valores = [x, y, z, elsa, dist]
    
    return valores


# Datos Fisicos:
yini = 1.                                       # Altura inicial dada por uno para comenzar a iterar
B    = 10.                                      # Ancho de fondo del canal
ss   = 2./1                                     # Relacion entre H/V
yn   = yini                                     # Altura normal del agua 
So   = 0.001                                    # Pendiente de fondo
n    = 0.013                                    # n Manning
Co   = 1.                                       # Sistema internacional
A    = b*yn+ss*yn**2                            # Area
Pm   = b + 2*yn*sp.sqrt(ss**2 + 1)
Q    = Co/n*A*sp.sqrt(So)*(A/Pm)**(2./3)       # Caudal
#Q    = 30.
tol  = 10.                                      # Toleracia
[x1, y1, z1] = [0., 3., 0.]                     # Punto inicial conocido
dx = 5.                                         # delta x


# Desarrollo:
[x, y, z, elev_sa, dist] = FGV(Q, B, So, ss, Co, x1, y1, z1, dx)
