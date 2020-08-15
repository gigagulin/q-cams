from blueqat import Circuit

dnum=input('input natural number less than 1024 > ')
temp=list(bin(int(dnum)))
del temp[0:2]
temp[0:0]=['0']*(10-len(temp))

c=Circuit(20)

for i in range(10):
	if temp[i]=='1':
		c.x[i]			# initial qubit setting

for i in range(9):
	c.cx[i,i+11]			# copy with shift
for i in range(10):
	c.cx[i,i+10]			# exclusive OR	

wgc=list(c.m[:].run(shots=1))	
xgc=wgc[0]

print(dnum,'=> Gray Code =', xgc[10:])