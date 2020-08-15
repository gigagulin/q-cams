# q-camgr (q-cam Grover) 15th Aug 2020 by Gigagulin

from blueqat import Circuit
import numpy as np
import qcam
import qmcn as q

# -------- Setting ------------------------------------------------------------------

N=6									# number of qubits :You can change.		
initial_a=np.array([0.075,0.205,0.156,0.14,0.204,0.22],dtype='float')	# initial probability distribution :You can change.
prob=np.array([0]*100,dtype='float')

#------- marking process -------------------------------------------------------------

def marking():								# marking |11111>

	c.h[N-1]
	q.cccccx(c,0,1,2,3,4,5)
	c.h[N-1]

	return

#------- Amplitude Amplification process ---------------------------------------------------------

def amp():

	c.h[:].x[:].h[N-1]
	q.cccccx(c,0,1,2,3,4,5)
	c.h[N-1].x[:].h[:]

	return

#---------- Main Body ----------------------------------------------------------------

pstep=0
ret='y'
while ret=='y':

	pstep+=1

	c=Circuit(N)
	qcam.propinit(N,c,initial_a)
	
	for i in range(pstep):
		marking()
		amp()
	
	master_a=np.array(c.run())

	for i in range(2**N):
		bstr = format(i, 'b')
		m=N-len(bstr)
		bstr='0'*m+bstr
		bstr = bstr[::-1]
		prob[pstep]=(master_a[i].real*master_a[i].real+master_a[i].imag*master_a[i].imag)*100

		cstr  = "{0.real: .2f} + {0.imag: .2f}j".format(master_a[i])
		print('  ',bstr,'  ', cstr, '  ','{0:,.1f}'.format(prob[pstep]),'%') 

	print(' ',pstep, end=' ')
	ret=input(' NEXT(Y/N)?')

# -------- Finalization  ----------------------------------------------------

sui=np.sum(initial_a)
ave=np.average(initial_a)
std=np.std(initial_a)


print(' initial-set =',initial_a)
print(' sum=','{0:,.2f}'.format(sui),' average=','{0:,.2f}'.format(ave),' std=','{0:,.3f}'.format(std))
for i in range(1,pstep+1):
	print(' TIME',i,'   ','{0:,.1f}'.format(prob[i])) 

