
#-------------- qmcn modules ( q- Multi Control Not Gates) Ver1.0   19th Mar.2020 by Gigagulin -------------------------+
#															+
# Purrpose of qmcn modules 	: Multi Control Not Gates for Blueqat							+				  
# Reqirments			: blueqat and numpy module are required.						+
# Arguments	: 1st   : c	: Circuit class of Blueqat								+
#		: 2nd-n : c0-cn : control qubits									+	 
#		: last		: target qubit										+		
# Return value	: None													+
#-----------------------------------------------------------------------------------------------------------------------+

import numpy as np

# -------- ccccccx -----------------------------------------------------------------------------------------------------+
# Purrpose	: 6 control qubits											+
#-----------------------------------------------------------------------------------------------------------------------+
	
def ccccccx(c,c0,c1,c2,c3,c4,c5,ta):
	
	ang=np.pi/32
	c.h[ta]								#Hadamard
									
	c.crz(ang)[c0,ta]						#Gray code 1
									#above this line is CX
	c.cx[c0,c1].crz(-ang)[c1,ta].cx[c0,c1].crz(ang)[c1,ta]		#Gray code 11 & 10
									#above this line is CCX(=Tofoli)
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 110 & 111
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 101 & 100
									#above this line is CCCX
	c.cx[c2,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1100 & 1101
	c.cx[c1,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1111 & 1110
	c.cx[c2,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1010 & 1011
	c.cx[c1,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1000 & 1000
									#above this line is CCCCX
	c.cx[c3,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11000 & 11001
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11011 & 11010
	c.cx[c2,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11110 & 11111
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11101 & 11100
	c.cx[c3,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10101 & 10111
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10110 & 10010
	c.cx[c2,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10010 & 10011
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10001 & 10000
									#above this line is CCCCCX
	c.cx[c4,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 110000 & 110001
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 110011 & 110010
	c.cx[c2,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 110110 & 110111
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 110101 & 110100
	c.cx[c3,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 111100 & 111101
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 111111 & 111110
	c.cx[c2,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 111010 & 111011
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 111001 & 111000
	c.cx[c4,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 101000 & 101001
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 101011 & 101010
	c.cx[c2,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 101110 & 101111
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 101101 & 101100
	c.cx[c3,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 100100 & 100101
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 100111 & 100110
	c.cx[c2,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 100010 & 100011
	c.cx[c1,c5].crz(-ang)[c5,ta].cx[c0,c5].crz(ang)[c5,ta]		#Gray Code 100001 & 100000
									#above this line is CCCCCCX
	c.h[ta]								#Hadamard

	return


# -------- cccccx ------------------------------------------------------------------------------------------------------+
# Purrpose	: 5 control qubits											+
#-----------------------------------------------------------------------------------------------------------------------+

def cccccx(c,c0,c1,c2,c3,c4,ta):
	
	ang=np.pi/16
	c.h[ta]								#Hadamard
									
	c.crz(ang)[c0,ta]						#Gray code 1
									#above this line is CX
	c.cx[c0,c1].crz(-ang)[c1,ta].cx[c0,c1].crz(ang)[c1,ta]		#Gray code 11 & 10
									#above this line is CCX(=Tofoli)
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 110 & 111
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 101 & 100
									#above this line is CCCX
	c.cx[c2,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1100 & 1101
	c.cx[c1,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1111 & 1110
	c.cx[c2,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1010 & 1011
	c.cx[c1,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1000 & 1000
									#above this line is CCCCX
	c.cx[c3,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11000 & 11001
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11011 & 11010
	c.cx[c2,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11110 & 11111
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 11101 & 11100
	c.cx[c3,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10101 & 10111
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10110 & 10010
	c.cx[c2,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10010 & 10011
	c.cx[c1,c4].crz(-ang)[c4,ta].cx[c0,c4].crz(ang)[c4,ta]		#Gray code 10001 & 10000
									#above this line is CCCCCX
	c.h[ta]								#Hadamard

	return


# -------- ccccx -------------------------------------------------------------------------------------------------------+
# Purrpose	: 4 control qubits											+
#-----------------------------------------------------------------------------------------------------------------------+

def ccccx(c,c0,c1,c2,c3,ta):								

	ang=np.pi/8
	c.h[ta]								#Hadamard
									
	c.crz(ang)[c0,ta]						#Gray code 1
									#above this line is CX
	c.cx[c0,c1].crz(-ang)[c1,ta].cx[c0,c1].crz(ang)[c1,ta]		#Gray code 11 & 10
									#above this line is CCX(=Tofoli)
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 110 & 111
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 101 & 100
									#above this line is CCCX
	c.cx[c2,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1100 & 1101
	c.cx[c1,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1111 & 1110
	c.cx[c2,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1010 & 1011
	c.cx[c1,c3].crz(-ang)[c3,ta].cx[c0,c3].crz(ang)[c3,ta]		#Gray code 1000 & 1000
									#above this line is CCCCX
	c.h[ta]								#Hadamard	

	return


# -------- cccx --------------------------------------------------------------------------------------------------------+
# Purrpose	: 3 control qubits											+
#-----------------------------------------------------------------------------------------------------------------------+

def cccx(c,c0,c1,c2,ta):								

	ang=np.pi/4
	c.h[ta]								#Hadamard
									
	c.crz(ang)[c0,ta]						#Gray code 1
									#above this line is CX
	c.cx[c0,c1].crz(-ang)[c1,ta].cx[c0,c1].crz(ang)[c1,ta]		#Gray code 11 & 10
									#above this line is CCX(=Tofoli)
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 110 & 111
	c.cx[c1,c2].crz(-ang)[c2,ta].cx[c0,c2].crz(ang)[c2,ta]		#Gray code 101 & 100
									#above this line is CCCX
	c.h[ta]								#Hadamard

	return
