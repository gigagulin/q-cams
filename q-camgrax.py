# q-cam Grover revx by Gigagulin Dec.2020

from blueqat import Circuit
import numpy as np
import qcam
import qmcn

# -------- Setting -------------------------------------------------------------------------

N=6								# number of cells
R=1								# number of registries

initial_a=np.array([0.4,0.6,0.3,0.7,0.2,0.8],dtype='float')	# initial set
max_step=200							# maximum steps
stepdist_a=np.array([[0]*N]*max_step,dtype='float')		# [step-number,cell-number]
probability_a=np.array([0]*N,dtype='float')
stepdist_a[0,0:N]=initial_a[0:N]
stdev_a=np.array([0]*max_step,dtype='float')
marking_prob_a=np.array([0]*max_step,dtype='float')

Reg0sb=0							# start qubit number of Reg0
Reg1sb=N							# start qubit number of Reg1
Reg2sb=N*2							# start qubit number of Reg2
Reg3sb=N*3							# start qubit number of Reg3

# --------  Q-Process -----------------------------------------------------------------------

def qproc():
								
	c.x[2,3,4].h[N-1]	# qubit in x gate depends on vector you want to mark　
	qmcn.cccccx(c,0,1,2,3,4,5)				
	c.h[N-1].x[2,3,4]	# qubit in x gate depends on vector you want to mark

	c.h[:].x[:].h[N-1]					# amplifier
	qmcn.cccccx(c,0,1,2,3,4,5)
	c.h[N-1].x[:].h[:]

	return

# -------- final output  -------------------------------------------------------------------

def fout(fstep, finitial_a,fprob_a):

	sui=np.sum(finitial_a)
	ave=np.average(finitial_a)
	std=np.std(finitial_a)
	print(' ')
	print(' initial-set =',finitial_a)
	print(' sum=','{0:,.2f}'.format(sui),' average=','{0:,.2f}'.format(ave),' std=','{0:,.3f}'.format(std))
	for i in range(1,fstep+1):
		print(' Step',i,'   ','{0:,.1f}'.format(fprob_a[i])) 

	return
		
# -------- Main Body -----------------------------------------------------------------------

pstep=0
ret='y'
while ret=='y':

	pstep+=1

	c=Circuit(N)			
	qcam.propinit(N,c,initial_a)

	for i in range(pstep):
		qproc()

	master_a=np.array(c.run())
	
	num,vector_a,prob_a=qcam.qvextract(N,R,1,master_a)
	stdev,probability_a=qcam.qcalcd(N, num,vector_a,prob_a)
	qcam.qcresultout(N, pstep, num,initial_a,vector_a,prob_a,stdev,probability_a)
	
	stepdist_a[pstep,0:N]=probability_a[0:N]
	marking_prob_a[pstep]=prob_a[35]*100			#Decimal number of marked vector
	stdev_a[pstep]=stdev

	ret=input(' NEXT(Y/N)?')

qcam.qcfinal(N, pstep, initial_a ,stepdist_a, stdev_a )
fout(pstep, initial_a ,marking_prob_a )
qcam.qcamplot(N, pstep, stepdist_a)

