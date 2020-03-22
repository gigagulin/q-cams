#probshot by Gigaulin 

from blueqat import Circuit
import numpy as np
import qcam

# -------- Setting ------------------------------------------------------------------

N=8									# number of qubits :You can change.
nrun=1000								# number of runs :You can change.
initial_a=np.array([0.4,0,1,0.5,1,0,0,0.2],dtype='float')		# initial probability distribution :You can change.
vector_a=np.array([[0]*N]*(2^N),dtype='int')
csum_a=np.array([[0]*N]*(2^N),dtype='float')
final_a=np.array([0]*N,dtype='float')
prob_a=np.array([0]*(2^N),dtype='float')


#---------- Main Body ----------------------------------------------------------------

ret='y'
while ret=='y':

	c=Circuit(N)
	qcam.propinit(N,c,initial_a)
	
	ans=c.m[:].run(shots=nrun)
	num=len(ans)
	kk=list(ans.keys())
	vv=list(ans.values())

	for i in range(num):
		vector_a[i,0:N]=list(kk[i])				# Result Vectors
		prob_a[i]=vv[i]/nrun

	stdev,final_a=qcam.qcalcd(N,num,vector_a, prob_a)
	qcam.qcresultout(N, 0, num, initial_a, vector_a, prob_a,stdev, final_a)

	ret=input(' NEXT(Y/N)?')
	


