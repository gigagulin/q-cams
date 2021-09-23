# q-frequency CA102  by Gigagulin 2021 Aug.

from blueqat import Circuit
import numpy as np
import qcam
from fractions import Fraction
import time

# -------- Setting ------------------------------------------------------------------

n_counter,n_unitary=10,7
N=n_counter+n_unitary						# number of qubit
initial_a=np.array([0.5]*n_counter+[0,0,0,0,0,0,1],dtype='float')	# initial set
xprob=np.array([0]*2**n_counter,dtype='float')			# Probability of counter

# --------  inversed QFT ------------------------------------------------------------

def IQFT():

	for i in range(n_counter):
		for j in range(i):
			c.cu1(-np.pi/(2**(i-j)))[i,j]
		c.h[i]

	return

# --------  CA102  --------------------------------------------------------------

def qproc(bf):

	for i in range(n_counter,N):	
		c.ccx[bf,1+i, i]

	return

# -------- Main Body --------------------------------------------------------------


c=Circuit(N)		
qcam.propinit(N,c,initial_a)


t0=time.time()

for m in range(n_counter):
	power=2**m
	for k in range(power):
		qproc(m)

IQFT()

master_a=np.array(c.run())
t1=time.time()

num,vector_a,prob_a=qcam.qvextract(N,1,1,master_a)

print(' counter phase probability guess')

for j in range(num):
	xdeci=0
	for i in range(n_counter):
		xdeci=int(vector_a[j,i])*2**(N-n_unitary-1-i)+xdeci
	xprob[xdeci]=prob_a[j]+xprob[xdeci]

for i in range(2**n_counter):
	if xprob[i]>=0.00005:
		phase=i/(2**n_counter)
		guess=str(Fraction(phase).limit_denominator(100))
		print('{:>5}'.format(i),'   ','{:>5,.4f}'.format(phase),end='  ')
		print('   ','{:>5,.5f}'.format(xprob[i]),end=' ')
		spa=' '*(6-len(guess))
		print(spa,':'+str(guess))

print(' Quantum Process Time =>',t1-t0)