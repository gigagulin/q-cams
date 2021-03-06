﻿# q-cam rule110 by Gigagulin 15th Aug 2020 

from blueqat import Circuit
import numpy as np
import qcam

# -------- Setting ------------------------------------------------------------------

N=6								# number of cells
R=4								# number of registers

initial_a=np.array([0,0,0,0,0,1],dtype='float')			# initial probability distribution
max_step=200							# maximum steps
stepdist_a=np.array([[0]*N]*max_step,dtype='float')		# [step-number,cell-number]
probability_a=np.array([0]*N,dtype='float')
stepdist_a[0,0:N]=initial_a[0:N]
stdev_a=np.array([0]*max_step,dtype='float')

Reg0sb=0							# start qubit number of Reg0
Reg1sb=N							# start qubit number of Reg1
Reg2sb=N*2							# start qubit number of Reg2
Reg3sb=N*3							# start qubit number of Reg3
Reg0lb=N-1							# last qubit number of Reg0
Reg1lb=N*2-1							# last qubit number of Reg1
Reg2lb=N*3-1							# last qubit number of Reg2
Reg3lb=N*4-1							# last qubit number of Reg3

# --------  Q-Process -----------------------------------------------------------------------

def qproc():
	

	c.ccx[Reg0sb,Reg0lb,Reg1sb]
	for i in range(1,N):
		c.ccx[i-1, i, Reg1sb+i]
	c.ccx[Reg0lb,Reg0sb+1,Reg2sb]
	for i in range(1,N-1):
		c.ccx[i-1, i+1, Reg2sb+i]
	c.ccx[Reg0sb,Reg0lb-1,Reg2lb]
	for i in range(N-1):
		c.ccx[i,i+1, Reg3sb+i]
	c.ccx[Reg0sb,Reg0lb,Reg3lb]
	for i in range(N):
		c.ccx[Reg1sb+i,Reg2sb+i,Reg3sb+i]
	for i in range(N):
		c.cx[i+1,i]
	for i in range(N):
		c.cx[Reg3sb+i,i]				# error correction
			
	return
		
# -------- Main Body --------------------------------------------------------------

pstep=0
ret='y'

while ret=='y':

	pstep+=1

	c=Circuit(N*R)

	pinitial_a=probability_a	
	if pstep==1:
		pinitial_a=initial_a
			
	qcam.propinit(N,c,pinitial_a)
	qproc()
	master_a=np.array(c.run())
	
	num,vector_a,prob_a=qcam.qvextract(N,R,1,master_a)
	stdev,probability_a=qcam.qcalcd(N, num,vector_a,prob_a)
	qcam.qcresultout(N, pstep, num,pinitial_a,vector_a,prob_a,stdev,probability_a)
	
	stepdist_a[pstep,0:N]=probability_a[0:N]
	stdev_a[pstep]=stdev

	ret=input(' NEXT(Y/N)?')

qcam.qcfinal(N, pstep, initial_a ,stepdist_a, stdev_a )
qcam.qcamplot(N, pstep, stepdist_a)

