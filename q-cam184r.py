# q-cam rule184 by Gigagulin

from blueqat import Circuit
import numpy as np
import qcam
import time

# -------- Setting ------------------------------------------------------------------

N=5								# number of cells
R=4								# number of registries
initial_a=np.array([1,0,1,0,1],dtype='float')			# initial probability distribution
max_step=200							# maximum steps

stepdist_a=np.array([[0]*N]*max_step,dtype='float')		# [step-number,cell-number]
probability_a=np.array([0]*N,dtype='float')
stepdist_a[0,0:N]=initial_a[0:N]
fr_a=np.array([0]*max_step,dtype='float')
stdev_a=np.array([0]*max_step,dtype='float')

Reg0sb=0							# start qubit number of Reg0
Reg1sb=N							# start qubit number of Reg1
Reg2sb=N*2							# start qubit number of Reg2
Reg3sb=N*3							# start qubit number of Reg3
Reg0lb=N-1							# last qubit number of Reg0
Reg1lb=N*2-1							# last qubit number of Reg1
Reg2lb=N*3-1							# last qubit number of Reg2
Reg3lb=N*4-1							# last qubit number of Reg3

# --------  Q-Process   ----------------------------------------------------------

def qproc():
	
	for i in range(N-1):
		c.ccx[i, 1+i, Reg1sb+i].x[Reg1sb+i]		#indicating congestion cells
	c.ccx[Reg0lb,Reg0sb,Reg1lb].x[Reg1lb]

	for i in range(N):
		c.cx[i, Reg2sb+i].x[Reg2sb+i]
	c.ccx[Reg2lb,Reg2sb,Reg3sb]

	for i in range(1,N):
		c.ccx[Reg2sb+i, Reg2sb-1+i, Reg3sb+i]		#indicating hopping correction cells  

	for i in range(N):
		c.cx[Reg1sb+i,Reg0sb+i]				# hopping with error
		c.cx[Reg3sb+i,Reg0sb+i]				# error correction
	
	return
		
# -------- Main Body --------------------------------------------------------------

pstep=0
ret='y'

while ret=='y':

	pstep+=1
	cst=time.time()

	c=Circuit(N*4)

	pinitial_a=probability_a	
	if pstep==1:
		pinitial_a=initial_a
			
	qcam.propinit(N,c,pinitial_a)

	t0=time.time()
	qproc()
	master_a=c.run()

	t1=time.time()		
	print(' ')
	
	num,vector0_a,prob_a=qcam.qvextract(N,R,1,master_a)
	num,vector1_a,prob_a=qcam.qvextract(N,R,2,master_a)
	stdev,probability_a=qcam.qcalcd(N, num,vector0_a,prob_a)
	flowrate,flowcell_a=qcam.calcflow(N,num,vector0_a,vector1_a,prob_a)
	qcam.jamout(N,pstep,num,pinitial_a,vector0_a,prob_a,stdev,probability_a,flowrate,flowcell_a)
	
	cet = time.time()

	print(' Quantum Process Time =>',t1-t0)
	print(' Total Process Time   =>',cet-cst)

	stepdist_a[pstep,0:N]=probability_a[0:N]
	fr_a[pstep]=flowrate
	stdev_a[pstep]=stdev

	ret=input(' NEXT(Y/N)?')

qcam.jamfinal(N, pstep, initial_a ,stepdist_a, stdev_a, fr_a)
qcam.qcamplot(N, pstep, stepdist_a)

