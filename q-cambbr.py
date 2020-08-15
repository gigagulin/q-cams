# q-cambb (q-cam periodical Ball & Box / carrior model ) by Gigagulin  15th Aug 2020

from blueqat import Circuit
import numpy as np
import qcam
import qmcn as q
import time

# -------- Setting ------------------------------------------------------------------

N=7								# number of box (qubit) 
R=3								# number of registers

initial_a=np.array([0,1,0,0,0,1,1],dtype='float')		# initial probability distribution 

prob=np.array([0]*100,dtype='float')
vector_a=np.array([[0]*N]*(2^N),dtype='int')
csum_a=np.array([[0]*N]*(2^N),dtype='float')
stepdist_a=np.array([[0]*N]*100,dtype='float')			# [step-number,cell-number]
probability_a=np.array([0]*N,dtype='float')
stepdist_a[0,0:N]=initial_a[0:N]
stdev_a=np.array([0]*100,dtype='float')

R0=list(range(N))
R0[N:N]=list(range(N))
R1=list(range(N,2*N))
R1[N:N]=list(range(N,2*N))
 
# to reduce number of qubits
R2=[14,15,16,17,14,15,16]
CA0=18								# ball carrier
CA1=19
CA2=20
XC0=21								# carrier empty or occupied flag
XC1=22
XC2=23
							

#-------  Q-Process Ball and Box (Soliton)  ----------------------------------

def qproc():
	
	for i in range(N):
		c.cx[R0[i],R1[i]]				# copy to R1

	c.x[XC0,XC1,XC2]					#  carrier flag sets when the carrier is empty, XCn=1

	for i in range(N):					# first round

		q.cccx(c,R0[i],R1[i],XC0,CA0)			# if R0=R1=1 and CA0(XC0=1) is empty, loading to CA0
		c.ccx[R1[i],CA0,XC0]				# if CA0 is loaded, XC0 turns
		
		q.ccccx(c,R0[i],R1[i],XC0,XC1,CA1)		# if R0=R1=1=1 and CA0 is occupied and CA1(XC1=1) is empty, loading to CA1
		c.ccx[R1[i],CA1,XC1].ccx[R1[i],CA1,XC0]		# if CA1 is loaded, XC1 turns
	
		q.cccccx(c,R0[i],R1[i],XC0,XC1,XC2,CA2)		# if R0=R1=1=1and CA0 and CA1 are occupied and CA2(XC2=1) is empty, loading to CA2
		c.ccx[R1[i],CA2,XC2]				# if CA2 is loaded, XC2 turns

		c.x[R0[i],CA2,XC2].ccx[R0[i],XC2,CA2]		# if R0 and R1 are empty and CA2(XC2=0) is occupied, unloding from CA2			 
		c.ccx[XC2,CA2,R0[i]].x[R0[i],R1[i]]		# if CA2 is changed, R0 turns		
		q.cccx(c,R0[i],R1[i],CA2,XC2)			# if CA2 is changed, R0[i] are ocuippied and R1[i] is empty, XC2 turns		 		
		c.x[R1[i],CA2,XC2]

		c.x[R0[i],CA1,XC1].ccx[R0[i],XC1,CA1]		# if R0 and R1 are empty and CA1(XC1=0) is occupied, unloding from CA1			 
		c.ccx[XC1,CA1,R0[i]].x[R0[i],R1[i]]		# if CA1 is changed, R0 turns		
		q.cccx(c,R0[i],R1[i],CA1,XC1)			# if CA1 is changed, R0[i] are ocuippied and R1[i] is empty, XC1 turns		 		
		c.x[R1[i],CA1,XC1]

		c.x[R0[i],CA0,XC0].ccx[R0[i],XC0,CA0]		# if R0 and R1 are empty and CA0(XC0=0) is occupied, unloding from CA0			 
		c.ccx[XC0,CA0,R0[i]].x[R0[i],R1[i]]		# if CA0 is changed, R0 turns		
		q.cccx(c,R0[i],R1[i],CA0,XC0)			# if CA0 is changed, R0[i] are ocuippied and R1[i] is empty, XC0 turns		 		
		c.x[R1[i],CA0,XC0]

	#for i in range(N):					# second round (unloading only)
								# to prevnet a contradiction of periodic boundary conditions
	for i in range(4):

		c.x[R0[i],R2[i],CA2,XC2]
		c.ccx[R0[i],XC2,CA2]				# if R0 is empty and CA2(XC2=0)is occupied, unloading from CA2
		c.ccx[CA2,XC2,R2[i]].x[R2[i]]			# if CA2 is chnaged, R2 turns
		c.ccx[R2[i],CA2,XC2].x[R0[i],CA2,XC2]		# if CA2 is chnaged, XC2 turns
		
		c.x[R0[i],R2[i],CA1,XC1]
		q.ccccx(c,R0[i],R2[i],XC2,XC1,CA1)		# if R0 is empty, CA1(XC1=0)is occupied and R2 and XC2 is unchanged, unloading from CA1
		c.ccx[CA1,XC1,R2[i]].x[R2[i]]			# if CA1 is chnaged, R2 turns
		c.ccx[R2[i],CA1,XC1].x[R0[i],CA1,XC1]		# if CA1 is chnaged, XC1 turns

		c.x[R0[i],R2[i],CA0,XC0]
		q.cccccx(c,R0[i],R2[i],XC2,XC1,XC0,CA0)		# if R0 is empty, CA0(XC0=0)is occupied and R2,XC2 and XC1 is unchanged, unloading from CA0
		c.ccx[CA0,XC0,R2[i]].x[R2[i]]			# if CA0 is chnaged, R2 turns
		c.ccx[R2[i],CA0,XC0].x[R0[i],CA0,XC0]		# if CA0 is chnaged, XC0 turns
		
	for i in range(N):
		c.cx[R1[i],R0[i]]		
	#	c.cx[R2[i],R0[i]]

	for i in range(4):
		c.cx[R2[i],R0[i]]


	return

#---------- Main Body ----------------------------------------------------------------

pstep=0
ret='y'

while ret=='y':

	pstep+=1
	cst=time.time()

	#c=Circuit(R*N+6)
	c=Circuit(24)

	pinitial_a=probability_a	
	if pstep==1:
		pinitial_a=initial_a
			
	qcam.propinit(N,c,pinitial_a)

	t0=time.time()
	qproc()

	master_a=np.array(c.run())
	t1=time.time()

	num,vector_a,prob_a=qcam.qvextract(N,R,1,master_a)
	stdev,probability_a=qcam.qcalcd(N, num,vector_a,prob_a)
	qcam.qcresultout(N, pstep, num,pinitial_a,vector_a,prob_a,stdev,probability_a)

	cet = time.time()

	print(' Quantum Process Time =>',t1-t0)
	print(' Total Process Time   =>',cet-cst)

	stepdist_a[pstep,0:N]=probability_a[0:N]
	stdev_a[pstep]=stdev

	ret=input(' NEXT(Y/N)?')
	
qcam.qcfinal(N, pstep, initial_a ,stepdist_a, stdev_a )
qcam.qcamplot(N, pstep, stepdist_a)




