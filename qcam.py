
#------------------- qcam (qcam modules) ver.1.5  19th Mra.2020 by Gigagulin -------------------------------------------+
#															+
# Purrpose of qcam modules 	: Using for q-cam programs but regarding "propinit" and  "qvextract"			+								+
# 				  it's possible to use more general purpose for Blueqat.				+
#  															+
#-----------------------------------------------------------------------------------------------------------------------+

import numpy as np
import matplotlib.pyplot as plt 

# -------- qvextract (Result Extract)  ---------------------------------------------------------------------------------+
# Purpose	: Extraction for vectors and probability from Blueqat							+
# Arguments	: 1st: NN		: Number of cells								+
#		: 2nd: RN		: Number of registers								+
#		: 3rd: RW		: Registory number you want to extract ( 1 except rule184)			+
#		: 4th: extract_a	: not shos result of Blueqat  (single array or list of amplitude probability)	+		
# Return value	: nm : number of the etracted vectors									+
#		: reg0_a : the etracted vectors (2D array [result-number,cell-number] consits of  '0' or '1' string)	+
#		: ereprov_a : Probability of the extracted vector (Single array)					+
#-----------------------------------------------------------------------------------------------------------------------+

def qvextract(NN,RN,RW,extract_a):
	
	reg0_a=np.array([['0']*NN]*2**NN)					# [result-number,cell-number]
	tempo_list=['0']*NN
	eprobv_a=np.array([0]*2**NN,dtype='complex')
	ereprov_a=np.array([0]*2**NN,dtype='float')

	nm=0
	for i in range(len(extract_a)):									
		if abs(extract_a[i])>0.0001:					# Decimal 'i' is the very vector!  			
			tempo1_list=list(bin(i))				# vector extraction butreversed. (Decimal to Binary)
			del tempo1_list[0:2]					# Delete '0','b' in tempo list 		
			tempo1_list[0:0]=['0']*(NN*RN-len(tempo1_list))		# Zero-fill to adjust digits
			reg0_a[nm,0:NN]=tempo1_list[(1-RW)*NN-1:-RW*NN-1:-1]	# slice and reverse	
			eprobv_a[nm]=extract_a[i]				# Probability amplitude
			nm+=1


	for i in range(nm):				
		ereprov_a[i]=eprobv_a[i].real*eprobv_a[i].real+eprobv_a[i].imag*eprobv_a[i].imag
						
	return nm,reg0_a,ereprov_a


# ------- qcalcd ( Calculation of Probability)   -----------------------------------------------------------------------+
# Purpose	: Calculation for probability of each cells from the probability of the etracted vectors		+
# Arguments	: 1st: NN		: Number of cells								+
#		: 2nd: cnum		: number of the etracted vectors						+
#		: 3rd: cvect_a		: the etracted vectors								+
#		: 4th: cprobare_a	: probability of the etracted vectors						+	
# Return value	: cstd : standard deviation (population) of the cell probability sets					+					+
#		: cproreal_a : calcurated probability of each cells							+
#-----------------------------------------------------------------------------------------------------------------------+

def qcalcd(NN, cnum, cvect_a, cprobare_a):
		
	cproreal_a=np.array([0]*NN,dtype = 'float')		
	csum_a=np.array([[0]*NN]*cnum,dtype = 'float')

	for j in range(cnum):
		for i in range(NN):
			csum_a[j,i]=float(cvect_a[j,i])*cprobare_a[j]
			
	cproreal_a=np.sum(csum_a,axis=0)

	cstd=np.std(cproreal_a)

	return cstd,cproreal_a


# -------- propinit (Initial Cell)s (=Reg0 qubits) Q-Setting <Rotation around a Y-axis>  -------------------------------+
# Purpose	: Setting qubits based on each cell probability by rotation gate					+
# Arguments	: 1st: NN		: Number of cells								+
#		: 2nd: c		: Circuit object of blueqat							+
#		: 3rd: ppinit_a		: a single array of each cell probability					+
# Return value	: None													+
#-----------------------------------------------------------------------------------------------------------------------+

def propinit(NN, c, ppinit_a):
		
	for i in range(NN):					
		theta=2* np.arcsin(np.sqrt(ppinit_a[i]))		# Correlation between Probability and Rotation Angle
		c.ry(theta)[i]						# You can change from ry to rx.		
	return


# --------- rotangle (Correlation between Probability and Rotation Angle <for U transfer>) -----------------------------+
# Purpose	: This function is a subroutine of the below function. Calculation					+
#-----------------------------------------------------------------------------------------------------------------------+

def rotangle(probi):

	Areal=probi.real
	Bimag=probi.imag
	Zero=1-Areal-Bimag

	if Zero==1:
		phi=0
		lam=0
	else:
		phi=np.arccos(np.sqrt(Areal/(1-Zero)))
		lam=np.arcsin(np.sqrt(Bimag/(1-Zero)))

	theta=2*np.arccos(np.sqrt(Zero))

	return theta,phi,lam


# -------- impropinit( Initial Cells (=Reg0 qubits) Q-Setting including imaginary) <U3 transfer >  ---------------------+
# Purpose	: Setting qubits based on each cell probability by rotation gate					+
# Arguments	: 1st: NN		: Number of cells								+
#		: 2nd: c		: Circuit class of blueqat							+
#		: 3rd: ppinit_a		: a single array of each cell probability					+
# Return value	: None													+
#-----------------------------------------------------------------------------------------------------------------------+

def impropinit(NN, c, ppinit_a):
		
	for i in range(NN):					
		ptheta,pphi,plam=rotangle(ppinit_a[i])		
		c.u3(ptheta,pphi,plam)[i]
	return


# -------- qcfinal  Finalization  --------------------------------------------------------------------------------------+
# Purpose	: print out of proppagation result	(whole time)							+
#-----------------------------------------------------------------------------------------------------------------------+

def qcfinal(NN, laststep,finitial_a, fstepdist_a, fstdev_a):

	fstdev_a[0]=np.std(finitial_a)

	print(' Time   /    Probability Distribution /  Sum   /  STDEV ')  

	for i in range(laststep+1):
		msum=np.sum(fstepdist_a[i])
		print('{0:3g}'.format(i),end='   ')	
		for j in range(NN):		
			print('{0:>5,.2f}'.format(fstepdist_a[i,j]),end=' ')

		print('    ','{0:>5,.2f}'.format(msum), end=' ')
		print(' ','{0:>5,.2f}'.format(fstdev_a[i]))

	return


# -------- qcresultout (Result out) ------------------------------------------------------------------------------------+
# Purpose	: print out for single step time									+
#-----------------------------------------------------------------------------------------------------------------------+

def qcresultout(NN,rpstep,rnum,rinitial_a,rresult_a,rprobab_a,pstd,prodistri_a):	
		
	print(' >Time','{0:3g}'.format(rpstep-1),' Inital Cell-Probability ',end='=> ')
	pinitstd=np.std(rinitial_a)
	for i in range(NN):
		print('{0:>5,.2f}'.format(rinitial_a[i]), end=' ')
	print(' sum=  ','{0:>5,.2f}'.format(rinitial_a.sum()),end=' ')
	print(' stdev= ','{0:>5,.2f}'.format(pinitstd))
	
	for i in range(rnum):
		prob=round(100*rprobab_a[i],3)
		print(' >Step','{0:3g}'.format(rpstep),' Result ','{0:3g}'.format(i+1), end='             => ')			
		for k in range(NN):			
			print('{0:>5}'.format(rresult_a[i,k]), end=' ')
		print(' Probability=','{0:,.3f}'.format(prob),'%')

	print(' >Time','{0:3g}'.format(rpstep),' Final Cell-Probability ',end=' => ')

	for i in range(NN):
		print('{0:>5,.2f}'.format(prodistri_a[i]), end=' ')
	print(' sum=  ','{0:>5,.2f}'.format(prodistri_a.sum()), end=' ')
	print(' stdev= ','{0:>5,.2f}'.format(pstd))
	print(' ')

	return


# -------- qcamplot (Plot Result ) -------------------------------------------------------------------------------------+
# Purpose	: Plot Result in graph 	: gray scale represents the each cell probability				+
# Arguments	: 1st: NN		: Number of cells								+
#		: 2nd: pl		: Last step time								+								+
#		: 3rd: stepr_a		: 2D array : whole time cell probability set					+
# Return value	: None													+
#-----------------------------------------------------------------------------------------------------------------------+

def qcamplot(NN, pl, stepr_a):

	fig =plt.figure(figsize=(NN/2,pl/2.2))

	plt.xlim(0,NN)
	plt.ylim(0,pl+1)

	plt.xlabel('Cell',fontsize=10)
	plt.ylabel('Time',fontsize=10)
	
	for i in range(pl+1):
		for j in range(NN):
			grays=str(1.0000-stepr_a[i,j])
			plt.plot(j+0.5,i+0.5,marker="s", color = grays , markersize=20)

	plt.show()

	return


# -------- Calculation of flow raet for rule 184 -----------------------------------------------------------------------+
# Purpose	: 	Calculation of flow raet (for rule 184)								+
#-----------------------------------------------------------------------------------------------------------------------+

def calcflow(NN, fnum,fvect_a,fcong_a,fprob_a):
		
	fsum_a=np.array([[0]*NN]*fnum,dtype = 'float')
	fcflow_a=np.array([0]*NN,dtype = 'float')

	for j in range(fnum):
		for i in range(NN):
			fsum_a[j,i]=float(fvect_a[j,i])*float(fcong_a[j,i])*fprob_a[j]

	fcflow_a=np.sum(fsum_a,axis=0)
	fr=np.sum(fcflow_a)

	return fr,fcflow_a


# -------- jamout (Result out for rule 184) -------------------------------------------------------------------------------+
# Purpose	: 	print out for single step time (for rule 184)							+
#-----------------------------------------------------------------------------------------------------------------------+

def jamout(NN,rpstep,rnum,rinitial_a,rresult_a,rprobab_a,pstd,prodistri_a,fr,fcell_a):	
		
	print(' >Time','{0:3g}'.format(rpstep-1),' Inital Cell-Probability ',end='=> ')
	pinitstd=np.std(rinitial_a)
	for i in range(NN):
		print('{0:>5,.2f}'.format(rinitial_a[i]), end=' ')
	print(' sum=  ','{0:>5,.2f}'.format(rinitial_a.sum()),end=' ')
	print(' stdev= ','{0:>5,.2f}'.format(pinitstd))
	
	for i in range(rnum):
		prob=round(100*rprobab_a[i],3)
		print(' >Time','{0:3g}'.format(rpstep),' Result ','{0:3g}'.format(i+1), end='             => ')			
		for k in range(NN):			
			print('{0:>5}'.format(rresult_a[i,k]), end=' ')
		print(' Probability=','{0:,.3f}'.format(prob),'%')

	print(' >Time','{0:3g}'.format(rpstep),' Final Cell-Probability ',end=' => ')

	for i in range(NN):
		print('{0:>5,.2f}'.format(prodistri_a[i]), end=' ')
	print(' sum=  ','{0:>5,.2f}'.format(prodistri_a.sum()), end=' ')
	print(' stdev= ','{0:>5,.2f}'.format(pstd))

	print(' >Time','{0:3g}'.format(rpstep),' Flow-Rate   /Cell       ',end='=> ')
	for i in range(NN):
		print('{0:>5,.2f}'.format(fcell_a[i]), end=' ')

	print(' Total=','{0:>5,.2f}'.format(fr))
	print(' ')

	return


# -------- jamfinal  (Finalization for rule 184)  ----------------------------------------------------------------------+
# Purpose	: print out of proppagation result for rule 184	(whole time)						+
#-----------------------------------------------------------------------------------------------------------------------+

def jamfinal(NN, lastep,finitial_a, fstepdist_a, fstdev_a, ffr_a):

	fstdev_a[0]=np.std(finitial_a)

	print(' Time   /    Probability Distribution /  Sum   /  FR  /  STDEV ')  

	for i in range(lastep+1):
		msum=np.sum(fstepdist_a[i])
		print('{0:3g}'.format(i),end='   ')	
		for j in range(NN):		
			print('{0:>5,.2f}'.format(fstepdist_a[i,j]),end=' ')

		print('    ','{0:>5,.2f}'.format(msum), end=' ')
		print(' ','{0:>5,.2f}'.format(ffr_a[i]), end=' ')
		print(' ','{0:>5,.2f}'.format(fstdev_a[i]))

	return
