from datetime import datetime
import random
import numpy
from operator import add
from itertools import chain
print "Program starts at: ",datetime.now().time()

p=5;NV=15                                                       # p=Field Size,NV=number of variables 
print "No of variables=",NV,"FieldSize",p

#Defining Polynomial ring
R=PolynomialRing(GF(p), NV, var_array=['x'])
R.inject_variables()
x=(list)(R.gens())
xv=vector(x+[1])

Qx=[]
for i in range(NV):
          j=0
          while(j<NV):
               Qx.append(x[i]*x[j])
               j+=1
#Qx=Qx
Qxv=vector(Qx+x+[1])


# The function Exper2 generates a polynomial P of degree deg_of_poly defined by product of deg_of_poly linear polynomials. If deg_of_poly=2 then P=l1*l2. Further, it calculates the deg_of_poly+1 iterative derivatives of order deg_of_poly-1. And then proceed by checking that do linear polynomials lies in the span of derivatives or not.
def Exper2(deg_of_poly):
	print "degree of polynomial",deg_of_poly
	P=1
	l=[]	
	# Defining P=l1*l2......ldeg_of_poly
	# l contains the list form of linear polynomials l1,l2,....
	for i in range(deg_of_poly):
		l.append([random.randrange(p) for j in range(NV+1)])		
		P=P*(xv.dot_product(vector(l[i])))
	print l
	Der=[]
	# Calculate the iterative derivatives
	for i in range(deg_of_poly+1):
		tmp=P
		for j in range(deg_of_poly-1):
			r=[random.randrange(p) for i in range(NV)]
			tmp=tmp.subs(dict(zip(x,[x[i]+r[i] for i in range(NV)])))-tmp
		t=[0]*(NV+1)
		for j in range(NV):
			t[j]=tmp.monomial_coefficient(x[j])
		t[NV]=tmp.constant_coefficient()
		Der.append(t)
	arr=matrix(GF(p),Der)
	print arr
	# Check that does l[i] lies in the span of derivatives or not
	for i in range(len(l)):
		try:
			ret=arr.solve_left(vector(l[i]))
			print i,ret
		except ValueError:
			print i+"th linear poly. not in span"

# The code works for all degrees with suitable Field size p.
Exper2(3)




