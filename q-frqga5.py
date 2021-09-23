﻿# q-frequency grover 4bit by Gigagulin 2021 Aug.

from blueqat import Circuit
import numpy as np
import qcam
import qmcn as q
from fractions import Fraction

# -------- Setting ------------------------------------------------------------------

n_counter,n_unitary=10,5
N=n_counter+n_unitary						# number of qubit
initial_a=np.array([0.5]*n_counter+[0.5,0.5,0.5,0.5,0.5],dtype='float')	# initial set
xprob=np.array([0]*2**n_counter,dtype='float')			# Probability of counter
ga0=n_counter
ga1=n_counter+1
ga2=n_counter+2
ga3=n_counter+3
ga4=n_counter+4

# --------  inversed QFT ------------------------------------------------------------------

def IQFT():

	for i in range(n_counter):
		for j in range(i):
			c.cu1(-np.pi/(2**(i-j)))[i,j]
		c.h[i]

	return
# --------  q-process  ------------------------------------------------------------------

def qproc(bf):

	c.h[ga4]
	q.cccccx(c,bf,ga0,ga1,ga2,ga3,ga4)
	c.h[ga4]

	c.cx[bf,ga0].ry(-np.pi/4)[ga0].cx[bf,ga0].ry(np.pi/4)[ga0].cx[bf,ga0]
	c.cx[bf,ga1].ry(-np.pi/4)[ga1].cx[bf,ga1].ry(np.pi/4)[ga1].cx[bf,ga1]
	c.cx[bf,ga2].ry(-np.pi/4)[ga2].cx[bf,ga2].ry(np.pi/4)[ga2].cx[bf,ga2]
	c.cx[bf,ga3].ry(-np.pi/4)[ga3].cx[bf,ga3].ry(np.pi/4)[ga3].cx[bf,ga3]
	c.cx[bf,ga4].ry(-np.pi/4)[ga4].cx[bf,ga4].ry(np.pi/4)[ga4].cx[bf,ga4]
	c.cx[bf,ga0].cx[bf,ga1].cx[bf,ga2].cx[bf,ga3].cx[bf,ga4]


	c.h[ga4]
	q.cccccx(c,bf,ga0,ga1,ga2,ga3,ga4)
	c.h[ga4]

	c.cx[bf,ga0].cx[bf,ga1].cx[bf,ga2].cx[bf,ga3].cx[bf,ga4]
	c.cx[bf,ga4].ry(-np.pi/4)[ga4].cx[bf,ga4].ry(np.pi/4)[ga4].cx[bf,ga4]
	c.cx[bf,ga3].ry(-np.pi/4)[ga3].cx[bf,ga3].ry(np.pi/4)[ga3].cx[bf,ga3]
	c.cx[bf,ga2].ry(-np.pi/4)[ga2].cx[bf,ga2].ry(np.pi/4)[ga2].cx[bf,ga2]
	c.cx[bf,ga1].ry(-np.pi/4)[ga1].cx[bf,ga1].ry(np.pi/4)[ga1].cx[bf,ga1]
	c.cx[bf,ga0].ry(-np.pi/4)[ga0].cx[bf,ga0].ry(np.pi/4)[ga0].cx[bf,ga0]

	return
# -------- Main Body --------------------------------------------------------------


c=Circuit(N)		
qcam.propinit(N,c,initial_a)

for m in range(n_counter):
	power=2**m
	for k in range(power):
		qproc(m)
IQFT()

master_a=np.array(c.run())

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
