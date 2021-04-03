# q-cambb (q-cam periodical Ball & Box / carrior model ) rev.3 by Gigagulin  <Mar.2021>

from blueqat import Circuit
import numpy as np
import qcam
import qmcn as q
import time

# -------- Setting ------------------------------------------------------------------

N=7								# number of box (qubit) 
R=3								# number of registries

initial_a=np.array([0,0.9,0,0,0,0.9,0.9],dtype='float')	# initial probability distribution 

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
R2=list(range(2*N,R*N))
R2[N:N]=list(range(2*N,R*N))

CA0=R*N								# ball carrier
CA1=R*N+1
CA2=R*N+2
XC0=R*N+3
XC1=R*N+4
XC2=R*N+5
						
#-------  Q-Process Ball and Box (Soliton)  ----------------------------------

def Proc1():
	
	for i in range(N):
		c.cx[R0[i],R1[i]]				# copy to R1

	c.x[XC0,XC1,XC2]					# carrier flag sets when the carrier is empty, XCn=1

	for i in range(N):					# first round loading

		q.cccx(c,R1[i],CA1,CA2,XC2)			# XC2 error correction

		c.ccx[R1[i],XC0,CA0]				# if R1=1 and XC0=1,loading to CA0
		c.ccx[R1[i],CA0,XC0]				# if CA0 is loaded, XC0 turns
		
		q.cccx(c,R1[i],XC0,XC1,CA1)			# if R1=1 and XC1=XC0=1, loading to CA1
		c.ccx[R1[i],CA1,XC1]				# if CA1 is loaded, XC1 turns
	
		q.ccccx(c,R1[i],XC0,XC1,XC2,CA2)		# if R1=1 and XC2=XC1=XC0=1, loading to CA2
		c.ccx[R1[i],CA2,XC2]				# if CA2 is loaded, XC2 turns

		c.ccx[R1[i],CA2,XC1].ccx[R1[i],CA1,XC0]		# XC1 and XC0 error correction

								# first round unloading
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

	return


def Proc2():

	for i in range(N):					# second round (unloading only)
								# to prevnet a contradiction of periodic boundary conditions
		c.x[R0[i]]					
		c.ccx[R0[i],CA2,R2[i]].ccx[R0[i],R2[i],CA2]	# Unloading from CA2 to R2[i]
		c.ccx[R0[i],CA1,R2[i]].ccx[R0[i],R2[i],CA1]	# Unloading from CA1 to R2[i]
		q.cccx(c,R0[i],CA0,CA1,R2[i])			# error correction
		c.ccx[R0[i],CA0,R2[i]].ccx[R0[i],R2[i],CA0]	# Unloading from CA0 to R2[i]
		c.ccx[R0[i],CA0,R2[i]]				# error correction
		c.x[R0[i]]

	return


def Proc3():
		
	for i in range(N):
		c.cx[R2[i],R0[i]]
		c.cx[R1[i],R0[i]]
	return

#---------- Main Body ----------------------------------------------------------------

pstep=0
ret='y'


while ret=='y':

	pstep+=1
	cst=time.time()

	c=Circuit(R*N+6)

	pinitial_a=probability_a	
	if pstep==1:
		pinitial_a=initial_a
			
	qcam.propinit(N,c,pinitial_a)

	t0=time.time()
	Proc1()
	Proc2()
	Proc3()

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




