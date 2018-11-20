# -*- coding: utf-8 -*-

import sys
reload(sys)
sys.setdefaultencoding('utf8')

import scipy as sp
import matplotlib.pyplot as plt

def yNormal(yini,B,Q,ss,n,S0,C0):
    yn = yini
    q = 30 #m**3/s
    tol = 10
    while tol>(10**(-6)):
        A = B*yn + ss*(yn**2)
        P = B + 2*yn*sp.sqrt(ss**2+1)
        dAdyn = B + 2*ss*yn
        dPdyn = 2*sp.sqrt(ss**2+1)

        fyn = (A**(5./3.))*(P**(-2./3.))-(n*q)/(C0*sp.sqrt(S0))
        dfyn = (5./3.)*(A**(2./3.))*(P**(-2./3.))*dAdyn - (2./3.)*(A**(5./3.))*(P**(-5./3.))*dPdyn
        yn2 = yn - (fyn/dfyn)

        tol = abs(yn2-yn)

        yn = yn2
    return yn

def NRalfa(Q,B,ss,n,Co,S0):
	y_a = 5.
	beta = 3./5

	tol = 10
	iteracion = 0
	while tol > (10**(-6)) and iteracion < 100:
		iteracion += 1
        A = B*y_a + ss*y_a**2
        dAdy = B + 2*ss*y_a
        l = y_a * sp.sqrt((ss**2)+1)
        P = B + 2*l
        dPdy = 2 * sp.sqrt((ss**2)+1)
        
        alpha = ((n*P**(2./3.)) / (Co*sp.sqrt(S0)))**beta
        fy = A - alpha*(Q**beta)
        dfdy = dAdy - (2./5.) * (((n*Q)/(Co*S0))**beta) * dPdy * P**(-beta)
        yf = y_a - fy/dfdy
        
        tol = abs(yf-y_a)
        
        y_a = yf

	return [A, yf]


# Datos
N    = 1000 		# Numero de nodos (divisiones en la malla 1D)
L    = 150000 		# longitud del tramo del rio a investigar
B    = 100 			# ancho de la base del canal
S    = 0.001 		# pendiente longitudinal del canal
n    = 0.045 		# coeficiente de Manning
ss   = 0 			# pendiente lateral (si ss=0, canal rectangular)
Co   = 1.485 		# factor de conversionde unidades
NC   = 1 			# numero de Courant, V_onda/V_numerica < 1. cond. de estabilidad
Tfin = 600 			# tiempo final de simulacion en minutos
Tfin = Tfin*60 		# tiempo final en segundos
g    = 32.2;


# SECCION 2
# Calculo de delta x y asignar valores al vector x, conociendo dx
dx = L/(N-1)
x = sp.zeros((1,N))
for i in range(N):
	x[0,i]=(i-1)*dx

# SECCION 3
# Condiciones iniciales (k=1, t=0)
Q = sp.zeros((N,N))
Q[0,:] = 250
y = sp.zeros((N,N))
y[0,:] = yNormal(1,B,Q[0,0],ss,n,S,Co)
A = sp.zeros((N,N))
A[0,:] = B*y[0,0]
V = sp.zeros((N,N))
V[0,:] = (Q[0,0]/A[0,0])

# SECCION  4
# Loop desde t = 0  a  t = Tfin
t = 0
k = 0
tiempo = [0]
while t < Tfin:
	# dt = NC*dx/(g*mean(y(k,1:N)))^.5
	# print V[k]
	dt = NC*dx/(sp.mean(V[k]))
	t += dt
	tiempo.append(t)
	tshow = t/60
# Condiciones  de  borde
	if t <= 150*60:
		Q[k+1 ,0] = 250 + 750/sp.pi * (1-sp.cos(sp.pi*t/(60*75))) # hidrograma de la descarga de la represa
	else:
		Q[k+1 ,0] = 250 # caudal constante luego de la descarga (volviendo a las condiciones originales)
		
	[A[k+1 ,0], y[k+1 ,0]] = NRalfa(Q[k+1 ,0],B,ss,n,Co,S) # NR para encontar valor de A e Y iteracionrando
	V[k+1 ,0] = Q[k+1 ,0]/A[k+1 ,0]
	# Moviendo estencil en el tiempo k, desde nodo 1 al nodo N
	for i in range(N-1):
		Q[k+1,i+1] = Q[k+1,i]-dt/dx*(A[k+1,i]-A[k,i]) # ecuacion de continuidad
		[A[k+1,i+1], y[k+1,i+1]] = NRalfa(Q[k+1,i+1],B,ss,n,Co,S) # ecuacion de momentum
		V[k+1,i+1] = Q[k+1,i+1]/A[k+1,i+1]

	# graficando...
	plt.figure(1)
	plt.draw() 	# este comando obliga a graficar para cada ciclo del loop
	plt.subplot(2,1,1)
	plt.plot(x[0,:],y[k,:], linewidth=2)  # altura de agua vs X
	plt.xlabel('x (m)')
	plt.ylabel('altura de agua (m)')
	ejes1 = [0, 150000, 0, 8]
	plt.axis(ejes1)

	plt.draw()
	plt.subplot(2,1,2)
	plt.plot(x[0,:],Q[k,:], linewidth=2)  # caudal de agua vs X
	plt.xlabel('x (m)')
	plt.ylabel('Caudal (ft^3/s)')
	ejes2 = [0, 150000, 0, 1600]
	plt.axis(ejes2)
	plt.title(['tiempo: '+str(tshow)])
	plt.show()

	k += 1
