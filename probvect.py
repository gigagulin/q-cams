#probvect by Gigaulin

from blueqat import Circuit
import numpy as np
import qcam

# -------- Setting ------------------------------------------------------------------

N=8									# number of cells :You can change.
R=1									# number of registries		
initial_a=np.array([0.9,0,0,1,0,0.6,0,0],dtype='float')			# initial probability distribution :You can change.
vector_a=np.array([[0]*N]*(2^N),dtype='int')
csum_a=np.array([[0]*N]*(2^N),dtype='float')
final_a=np.array([0]*N,dtype='float')

#---------- Main Body ----------------------------------------------------------------

ret='y'
while ret=='y':

	c=Circuit(N)
	qcam.propinit(N,c,initial_a)

	master_a=np.array(c.run())

	num,vector_a,prob_a=qcam.qvextract(N,R,1,master_a)
	stdev,final_a=qcam.qcalcd(N, num,vector_a,prob_a)
	qcam.qcresultout(N, 0, num, initial_a, vector_a, prob_a,stdev, final_a)

	ret=input(' NEXT(Y/N)?')



