# q-cam rule30 by Gigagulin

from blueqat import Circuit
import numpy as np
import qcam

# -------- Setting ------------------------------------------------------------------

N=8								# number of cells
R=3								# number of registries

initial_a=np.array([0,0,0,0,1,0,0,0],dtype='float')		# initial probability distribution
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
	
	for i in range(N-1):
		c.ccx[i, 1+i, Reg1sb+i]
	c.ccx[Reg0lb,Reg0sb,Reg1lb]
	for i in range(N):
		c.cx[i, Reg2sb+i]
		
	for i in range(N-1):
		c.cx[Reg2sb+i,i+1]
	c.cx[Reg2lb,Reg0sb]
	for i in range(N-1):
		c.cx[Reg2sb+i+1, i]
	c.cx[Reg2sb,Reg0lb]

	for i in range(N):
		c.cx[Reg1sb+i,i]				# error correction
			
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

