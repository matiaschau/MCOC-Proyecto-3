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