# q-cam rule184i (imaginary version) by Gigagulin

from blueqat import Circuit
import numpy as np
import qcam
import time

# -------- Setting -------------------------------------------------------------------------------

N=6								# number of cells
initial_a=np.array([0.4,0,1,0.5,1,0],dtype='complex')	# initial probability distribution
max_step=200							# maximum steps

stepdist_a=np.array([[0]*N]*max_step,dtype='complex')		# [step-number,cell-number]
probability_a=np.array([0]*N,dtype='complex')
stepdist_a[0,0:N]=initial_a[0:N]

Reg0sb=0							# start qubit number of Reg0
Reg1sb=N							# start qubit number of Reg1
Reg2sb=N*2							# start qubit number of Reg2
Reg3sb=N*3							# start qubit number of Reg3
Reg0lb=N-1							# last qubit number of Reg0
Reg1lb=N*2-1							# last qubit number of Reg1
Reg2lb=N*3-1							# last qubit number of Reg2
Reg3lb=N*4-1							# last qubit number of Reg3


# --------  Q-Process 1. indicating congestion cells -------------------------------------------

def Proc1():
	
	for i in range(N-1):
		c.ccx[i, 1+i, Reg1sb+i].x[Reg1sb+i]

	c.ccx[Reg0lb,Reg0sb,Reg1lb].x[Reg1lb]

	return


# --------  Q-Process 2.   indicating hopping correction cells  --------------------------------

def Proc2():

	for i in range(N):
		c.cx[i, Reg2sb+i].x[Reg2sb+i]

	c.ccx[Reg2lb,Reg2sb,Reg3sb]

	for i in range(1,N):
		c.ccx[Reg2sb+i, Reg2sb-1+i, Reg3sb+i]

	return


# -------- Q-Process 3.   hopping -------------------------------------------------------------

def Proc3():

	for i in range(N):
		c.cx[Reg1sb+i,Reg0sb+i]				# hopping with error
		c.cx[Reg3sb+i,Reg0sb+i]				# error correction
	
	return


# -------- Result extract -----------------------------------------------------------------------

def Extract(extract_a):

	num=0	
	reg0vectrs_list=['0']*500
	reg1vectrs_list=['0']*500
	reg0vector_a=np.array([['0']*N]*500)			# [result-number,cell-number]
	reg1vector_a=np.array([['0']*N]*500)			# [result-number,cell-number]
	eprobv_a=np.array([0]*500,dtype='complex')
	ereprov_a=np.array([0]*500,dtype='float')
	eimprov_a=np.array([0]*500,dtype='float')


	for i in range(len(extract_a)):				# Decimal 'i' is the very vector.  

				
		if abs(extract_a[i])>0.0001:				
			reg0vectrs_list[num]=bin(i)[-N:]	# # This is result vector. (Decimal to Binary)
			reg1vectrs_list[num]=bin(i)[-2*N:-N]	# Congestion Vectors
			eprobv_a[num]=extract_a[i]
			num+=1	

	for i in range(num):		
		tempo_list=list(reversed(reg0vectrs_list[i]))
		reg0vector_a[i,0:N]=tempo_list
		tempo_list=list(reversed(reg1vectrs_list[i]))
		if tempo_list==['b','0']:
			reg1vector_a[0,0:N]=['0']*N
		else:
			reg1vector_a[i,0:N]=tempo_list

		ereprov_a[i]=eprobv_a[i].real*eprobv_a[i].real
		eimprov_a[i]=eprobv_a[i].imag*eprobv_a[i].imag
		
	return num,reg0vector_a,reg1vector_a,ereprov_a,eimprov_a


# -------- Calculation of probability distribution ----------------------------------------------

def Calcprodistri(cnum, cvect_a, cprobare_a, cprobaim_a):
		
	cproreal_a=np.array([0]*N,dtype = 'float')	
	cproimag_a=np.array([0]*N,dtype = 'float')		
	csumre_a=np.array([[0]*N]*cnum,dtype = 'float')
	csumim_a=np.array([[0]*N]*cnum,dtype = 'float')

	for j in range(cnum):
		for i in range(N):
			csumre_a[j,i]=float(cvect_a[j,i])*cprobare_a[j]
			csumim_a[j,i]=float(cvect_a[j,i])*cprobaim_a[j]
			
	for i in range(N):
		qr=0
		qi=0
		for j in range(cnum):
			qr=qr+csumre_a[j,i]
			qi=qi+csumim_a[j,i]
		cproreal_a[i]=qr				#Probability of real
		cproimag_a[i]=qi				#Probability of imaginary

	return cproreal_a,cproimag_a


# -------- Calculation of flow rate -----------------------------------------------------------

def Calcflowr(fnum,fvect_a,fcong_a,fprore_a,fproim_a):

	fsumr_a=np.array([[0]*N]*fnum,dtype = 'float')
	fsumi_a=np.array([[0]*N]*fnum,dtype = 'float')
	fcflowre_a=np.array([0]*N,dtype = 'float')
	fcflowim_a=np.array([0]*N,dtype = 'float')

	for j in range(fnum):
		for i in range(N):
			fsumr_a[j,i]=float(fvect_a[j,i])*float(fcong_a[j,i])*fprore_a[j]
			fsumi_a[j,i]=float(fvect_a[j,i])*float(fcong_a[j,i])*fproim_a[j]

	for i in range(N):
		qfr=0
		qfi=0
		for j in range(fnum):
			qfr=qfr+fsumr_a[j,i]
			qfi=qfi+fsumi_a[j,i]
		fcflowre_a[i]=qfr
		fcflowim_a[i]=qfi		
	
	fr=0
	fi=0
	for i in range(N):
		fr=fr+fcflowre_a[i]
		fi=fi+fcflowim_a[i]
	return fr,fi,fcflowre_a,fcflowim_a

		
# -------- Result out -----------------------------------------------------------------------

def Resultout(rnum,rinitial_a,rresult_a,rpreal_a,rpimag_a,rpdreal_a,rpdimag_a,rfrre,rfrim,fcellre_a,fcellim_a):
	
	#Vector Results --------------------------------	
	for i in range(rnum):
		probre=round(100*rpreal_a[i],3)
		probim=round(100*rpimag_a[i],3)
		print(' >Step','{0:3g}'.format(pstep),' Result ','{0:3g}'.format(i+1), end='  =>  ')			
		for k in range(N):			
			print('{0:>5}'.format(rresult_a[i,k]), end=' ')
		print(' Real P.=','{0:,.3f}'.format(probre),'%', end=' ')
		print(' Imag P.=','{0:,.3f}'.format(probim),'%')
	print('')

	#Initial Complex represntation --------------------------------
	print(' >Step','{0:3g}'.format(pstep),' Initial Complex-rep.   ',end=' =>  ')
	sumc=0
	for i in range(N):
		print('{0:>5,.2f}'.format(rinitial_a[i]), end=' ')
		sumc=sumc+rinitial_a[i]
	print(' sum=  ','{0:>5,.2f}'.format(sumc))


	
	#Final Complex represntation --------------------------------
	print(' >Step','{0:3g}'.format(pstep),' Final Complex-rep.     ',end=' =>  ')	
	rsumx=0
	for i in range(N):
		finalc=rpdreal_a[i]+1j*rpdimag_a[i]
		rsumx=rsumx+finalc
		print('{0:>5,.2f}'.format(finalc), end=' ')		
	print(' sum=  ','{0:>5,.2f}'.format(rsumx))
	print('')

	#Cell Probability--------------------------------
	print(' >Step','{0:3g}'.format(pstep),' Real  Cell-Probability ',end=' => ')
	rsumr=0
	for i in range(N):
		print('{0:>5,.2f}'.format(rpdreal_a[i]), end=' ')
		rsumr=rsumr+rpdreal_a[i]
	print(' sum=  ','{0:>5,.2f}'.format(rsumr))

	
	print(' >Step','{0:3g}'.format(pstep),' Imag  Cell-Probability ',end=' => ')
	rsumi=0
	for i in range(N):
		print('{0:>5,.2f}'.format(rpdimag_a[i]), end=' ')
		rsumi=rsumi+rpdimag_a[i]
	print(' sum=  ','{0:>5,.2f}'.format(rsumi))


	#flow rate  ------------------------------------------
	print(' >Step','{0:3g}'.format(pstep),' Real Flow-Rate/Cell     ',end='=> ')
	for i in range(N):
		print('{0:>5,.2f}'.format(fcellre_a[i]), end=' ')

	print(' Total=','{0:>5,.2f}'.format(rfrre))

	print(' >Step','{0:3g}'.format(pstep),' Imag Flow-Rate/Cell     ',end='=> ')
	for i in range(N):
		print('{0:>5,.2f}'.format(fcellim_a[i]), end=' ')

	print(' Total=','{0:>5,.2f}'.format(rfrim))
	print(' ')

	return
	

# -------- Main Body ---------------------------------------------------------------------------

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

	Proc1()
	Proc2()
	Proc3()
	
	ans=c.run()
	t1=time.time()

	master_a=np.array(ans)
		
	print(' ')
	
	nov,reg0vectors_a,reg1vectors_a,reprob_a,improb_a=Extract(master_a)

	pdreal_a,pdimag_a=Calcprodistri(nov,reg0vectors_a,reprob_a,improb_a)

	fratere,frateim,flowreal_a,flowimag_a=Calcflowr(nov,reg0vectors_a,reg1vectors_a,reprob_a,improb_a)	

	Resultout(nov,pinitial_a,reg0vectors_a,reprob_a,improb_a,pdreal_a,pdimag_a,fratere,frateim,flowreal_a,flowimag_a)
	
	cet = time.time()

	print(' Quantum Process Time =>',t1-t0)
	print(' Total Process Time   =>',cet-cst)

	for i in range(N):
		probability_a[i]=pdreal_a[i]+1j*pdimag_a[i]

	stepdist_a[pstep,0:N]=probability_a[0:N]

	ret=input(' NEXT(Y/N)?')


# -------- Finalization  -----------------------------------------------------------------------------

for i in range(pstep+1):
	print(' Step ',i,end=' => ')
	
	for j in range(N):
		print('{0:>5,.2f}'.format(stepdist_a[i,j]),end=' ')

	print(' ')

