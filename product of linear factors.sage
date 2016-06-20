# In this program, given a degree-d polynomial P known to be of the form l1*l2*....*ld  where l1,...,ld are randomly chosen polynomials of degree 1 and our goal is to recover back the polynomials l1,...,ld.

from datetime import datetime
import random
import numpy
from operator import add
from itertools import chain
print "Program starts at: ",datetime.now().time()

p=5;NV=15                                                       # Let p = Field Size, NV = number of variables 
print "No of variables=",NV,"FieldSize",p

#Define Polynomial ring 'R'
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


# The below defined function 'Exper2' generates a polynomial P of degree 'deg_of_poly' which is product of 'deg_of_poly' linear polynomials. Ex: if deg_of_poly=2 then P=l1*l2. Additionally, it calculates the deg_of_poly+1 iterative derivatives of order deg_of_poly-1. And then checks whether does linear polynomials lie in the span of derivatives?
def Exper2(deg_of_poly):
	print "degree of polynomial",deg_of_poly
	P=1
	l=[]	
	# Define P=l_1*l_2......l_deg_of_poly
	# l contains the list form of linear polynomials l_1,l_2,....,l_deg_of_poly
	for i in range(deg_of_poly):
		l.append([random.randrange(p) for j in range(NV+1)])		
		P=P*(xv.dot_product(vector(l[i])))
	print l
	Der=[]
	# Calculate the iterative derivatives and store the derivatives in the list 'Der'
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
	# Check that does linear polynomials in 'l' lies in the span of derivatives?
	for i in range(len(l)):
		try:
			ret=arr.solve_left(vector(l[i]))
			print i,"th linear is recoverable",ret
		except ValueError:
			print i+"th linear poly. not in span"

Exper2(3)




