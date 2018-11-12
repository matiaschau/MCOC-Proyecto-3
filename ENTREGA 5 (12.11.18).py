# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 22:26:49 2018

@author: Tomás
"""

import scipy as sp

#Datos Fisicos
yini = 1. #Altura inicial dada por uno para comenzar a iterar

b = 10. #Ancho de fondo del canal
ss = 1./2 #Relacion entre H/V
yn = yini #Altura normal del agua 
So = 0.001 #Pendiente de fondo
n = 0.013
Co = 1. #Depende del materia ###CAMBIAR
A = b*yn+ss*yn**2 #Area
Pm = b + 2*yn*sp.sqrt(ss**2 + 1)
Q = Co/n*A*sp.sqrt(So)*(A/Pm)**(2./3) #Caudal

tol = 10.


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
    y_ini = yini
    x1 = x1
    y1 = y1
    z1 = z1
    dx = dx
    
    if Co == 1:
        g = 9.81
    else:
        g = 32.2
        
    i = 0 #valor de entrada
    while tol > 1./10**6.:
	A = b*yn+ss*yn**2 #Area
	Ap = b+2.*ss*yn #Derivada del area
	Pm = b + 2.*yn*sp.sqrt(ss**2 +1) #Perimetro mojado
	Pmp = 2.*sp.sqrt(ss**2 +1) #Derivada del perimetro mojado c/r a yn

	fyn = A**(5./3)*Pm**(-2./3) - (n*Q)/(Co*sp.sqrt(So))
	fynp = 5./3 * Pm**(-2./3) * A**(2./3) * Ap - 2./3 *A**(5./3) * Pm**(-5./3) *Pmp
	yn2 = yn - fyn/fynp

	tol = abs(yn2-yn)
	yn = yn2
	if tol < 1./10**6:
		print yn
    yn = yn
    
    diff = 0.001
    
    while y[i+1] > (yn - diff):
        i = i + 1
        dist = dist + dx
        hola = NR_FGV(y_ini, b, Q, n, So, ss, Co, dx, y[i], z[i])
    
    elsa = hola[0] + hola[1]     #elevacion superficie de agua = y + z
    
    valores = [x, y, z, elsa, dist]
    
    return valores
    
