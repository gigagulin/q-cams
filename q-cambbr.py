# q-cambb (q-cam Ball & Box) by Gigagulin

from blueqat import Circuit
import numpy as np
import qcam
import qmcn as q
import time

# -------- Setting ------------------------------------------------------------------

N=9									# number of box (qubit) 
M=4									# number of ball
R=2									# number of registries

CI=[0,1,2,3,4,5,6,7,8,0,1,2,3]
flagA=2*N
flagB=2*N+1
flagC=2*N+2

initial_a=np.array([0,0,1,1,1,0,0,0,1],dtype='float')	# initial probability distribution 

prob=np.array([0]*100,dtype='float')
vector_a=np.array([[0]*N]*(2^N),dtype='int')
csum_a=np.array([[0]*N]*(2^N),dtype='float')
stepdist_a=np.array([[0]*N]*100,dtype='float')				# [step-number,cell-number]
probability_a=np.array([0]*N,dtype='float')
stepdist_a[0,0:N]=initial_a[0:N]
stdev_a=np.array([0]*100,dtype='float')


#-------  Q-Process Ball and Box (Soliton)  ----------------------------------

def Proc():

	for i in range(N):
		c.cx[i,i+N]

	c.x[flagA,flagB,flagC]

	for i in range(N):

		q.ccccx(c,i+N,CI[i+1],CI[i+2],CI[i+3],CI[i+4])
		c.ccx[CI[i+4],CI[i+3],flagA]
	
		q.ccccx(c,i+N,flagA,CI[i+1],CI[i+2],CI[i+3])
		c.ccx[CI[i+3],CI[i+2],flagB]
				
		q.cccx(c,i+N,flagB,CI[i+1],CI[i+2])
		c.ccx[CI[i+2],CI[i+1],flagC]
		
		c.ccx[i+N,flagC,CI[i+1]]

		c.ccx[CI[i+4],CI[i+3],flagA].ccx[CI[i+3],CI[i+2],flagB].ccx[CI[i+2],CI[i+1],flagC]

	for i in range(N):	
		c.cx[i+N,i]

	return

#---------- Main Body ----------------------------------------------------------------

pstep=0
ret='y'

while ret=='y':

	pstep+=1
	cst=time.time()

	c=Circuit(N*R+3)

	pinitial_a=probability_a	
	if pstep==1:
		pinitial_a=initial_a
			
	qcam.propinit(N,c,pinitial_a)

	t0=time.time()
	Proc()

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




