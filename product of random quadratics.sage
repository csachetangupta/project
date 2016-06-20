# In this program, given a quartic polynomial P known to be of the form Q1*Q2  where Q1 and Q2 are randomly chosen quadratic polynomials and we try to recover back the polynomials Q1 and Q2. 

from datetime import datetime
import random
import numpy
from operator import add
from itertools import chain
print "Program starts at: ",datetime.now().time()

p=5;NV=5                                                         # Let p = Field Size, NV = number of variables 
SampleSize=p^2+200                                               # "SampleSize"=number of samples to be used for calculating Gowers norm
print "SampleSize",SampleSize
print "No of variables =",NV,"FieldSize",p

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
# Qx = Qx
Qxv=vector(Qx+x+[1])

def dectobasep(num,L):
          alpha=[0]*L;len=L-1
          while(num):
              alpha[len]=(int)(num%p)
              num=(int)(num/p)
              len-=1
          return alpha


def getPermutations(L,NV):
# getPermutations function returns sum of product of all the possible permutations of the lists in L.
# For example if L=[L1,L2,L3] then getPermutations returns L1*L2*L3+L1*L3*L2+ L2*L1*L3+L2*L3*L1+ L3*L1*L2+L3*L2*L1 where L1*L2=multiply(transpose(L1),L2)
         if(len(L)==1):
             return L[0]
         ret=matrix(GF(p),[[0]*(NV^(len(L)-1))]*NV)
         for i in range(len(L)):
             ret+=transpose(L[i])*getPermutations(L[:i]+L[i+1:],NV)
         return matrix(list(chain.from_iterable(ret)))


def CalculateGowersNorm(POLY,MaxDegree,NV):
# Below is the desciption of the variables:
# GN will be equal to the Gowers norm value of the polynomial POLY.
# Variable 'List' refers to the list of the random values of the variables of the derivative polynomial DP.
# 'L' contains the 'MaxDegree' number of equi-sized lists extracted from above random "List".
# SOL is a 1-dimensional matrix of size 'NV^MaxDegree' obtain by flattening the 'MaxDegree'-dimensional matrix 'M' which we got by appealing to getPermutation function defined above. The value of M[i][j][k] is equal to the L1[i].L2[j].L3[k]+L1[i].L2[k].L3[j]+L1[j].L2[i].L3[k]+L1[j].L2[k].L3[i]+L1[k].L2[i].L3[j]+L1[k].L2[j].L3[i]. In SOL, value of M[i][j][k] gets stored at position i*NV^2+j*NV+k in SOL.
# DPValueAtL is the value of the polynomial DP obtained by assigning variables with the values in list L.
                      pthroot = cos(2*pi/p)+I*sin(2*pi/p)
                      GN=0
                      for j in range(SampleSize):
                          List = [random.randrange(p) for i in range(NV*MaxDegree)]
                          L=[]
                          for i in range(MaxDegree):
                              L.append(matrix(List[i*NV:(i+1)*NV]))
                          SOL=matrix(GF(p),getPermutations(L,NV))
                          DPValueAtL=(POLY*transpose(SOL))[0][0]
                          GN+=pthroot^(int)(DPValueAtL)
                      GN=(float(abs(GN))*1/SampleSize)^(1/2^MaxDegree)
                      return GN

# generate polynomials P(=Q1*Q2), Q1 and Q2
def getPoly():
	Q1 = Qxv.dot_product(vector([random.randrange(p) for j in range(NV^2+NV+1)]))
	Q2 = Qxv.dot_product(vector([random.randrange(p) for j in range(NV^2+NV+1)]))
	return Q1*Q2,Q1,Q2
	
prod,Q1,Q2 = getPoly()
# print prod,Q1,Q2
h1=[]
h2=[]
derQ1withh1=[]
derQ1withh2=[]
doublederQ1=[]
derQ2withh1=[]
derQ2withh2=[]
doublederQ2=[]
derprod=[]
# calculate the derivatives
# derQ1withh1 will contain the derivatives of Q1 with first random direction h1 in all 3 iterations. Similarly, derQ2withh1, derQ1withh2, derQ2withh2.
# derprod will hold the iterated derivatives of P in directions h1, h2
for j in range(3):				#If P=Q1*Q2 then according to algorithm, 3 derivatives will be effective with high probability
	h1.append([random.randrange(p) for i in range(NV)])
	h2.append([random.randrange(p) for i in range(NV)]) 
	tmp=Q1
	# print tmp,dict(zip(x,[x[i]+h1[i] for i in range(NV)]))
	# print tmp.subs(dict(zip(x,[x[i]+h1[i] for i in range(NV)])))-tmp
	derQ1withh1.append(tmp.subs(dict(zip(x,[x[i]+h1[j][i] for i in range(NV)])))-tmp)
	derQ1withh2.append(tmp.subs(dict(zip(x,[x[i]+h2[j][i] for i in range(NV)])))-tmp)
	doublederQ1.append(derQ1withh1[j].subs(dict(zip(x,[x[i]+h2[j][i] for i in range(NV)])))-derQ1withh1[j])
	tmp=Q2
	derQ2withh1.append(tmp.subs(dict(zip(x,[x[i]+h1[j][i] for i in range(NV)])))-tmp)
	derQ2withh2.append(tmp.subs(dict(zip(x,[x[i]+h2[j][i] for i in range(NV)])))-tmp)
	doublederQ2.append(derQ2withh1[j].subs(dict(zip(x,[x[i]+h2[j][i] for i in range(NV)])))-derQ2withh1[j])
	tmp=prod	
	tmp=tmp.subs(dict(zip(x,[x[i]+h1[j][i] for i in range(NV)])))-tmp
	tmp=tmp.subs(dict(zip(x,[x[i]+h2[j][i] for i in range(NV)])))-tmp
	derprod.append(tmp)

def GenerateListForm(poly,degree):
# This function outputs the list form of the polynomial 'poly' with the only-one restriction that it consists only of the coefficients of monomials of degree="degree".
# Below ret_list represents the list form to be returned.
	Monom=poly.monomials()
	coeffs=poly.coefficients()
	ret_list=[0]*(NV^degree)
	for i in range(len(Monom)):
		if(Monom[i].total_degree()==degree):
			var=Monom[i].variables()
			count=0	
			for j in range(len(var)):
				deg_of_varj=Monom[i].degree(var[j])
				index=x.index(var[j])
				for m in range(deg_of_varj):
					count=count*NV+index
			ret_list[count]=coeffs[i]
	return ret_list

def recover_Quad():
# find the linear comb. which has high Gowers norm
	i=1
	mgn=0
	si=[0]*3
	while(i<p^3):
		i_inbasep=dectobasep(i,3)
		tmp=0
		poly=0
		for j in range(3):
			tmp=tmp+i_inbasep[j]*derprod[j]
			poly=poly+i_inbasep[j]*(derQ1withh1[j]*derQ2withh2[j]+derQ1withh2[j]*derQ2withh1[j]+doublederQ1[j]*(derQ2withh1[j]+derQ2withh2[j]) + doublederQ2[j]*(derQ1withh1[j]+derQ1withh2[j])+doublederQ1[j]*doublederQ2[j])
		ListForm=GenerateListForm(tmp,2)
		gn=CalculateGowersNorm(matrix(GF(p),ListForm),2,NV)
		R_poly=tmp-poly
		if(gn>mgn and R_poly==0):
			si=i_inbasep;mgn=gn
		i=i+1
	return si

si=recover_Quad()
print si


lowranklc=0	
# lowranklc refer to the linear combination of the derivatives of derprod which has the highest gowers norm
for j in range(3):
	lowranklc=lowranklc+si[j]*derprod[j]
print "lowranklc",lowranklc



def generate_list_form_of_linear(tmp):
# the below function generates the list form of linear polynomial tmp
		t=[0]*(NV+1)
		for j in range(NV):
			t[j]=tmp.monomial_coefficient(x[j])
		t[NV]=tmp.constant_coefficient()
		return t

Der=[]
tmp=lowranklc
# verify that does linear polynomials lie in the span of the derivative of tmp?
for i in range(17):
# we have evaluated 17 derivatives in random directons
		r=[random.randrange(p) for j in range(NV)]
		t=generate_list_form_of_linear(tmp.subs(dict(zip(x,[x[j]+r[j] for j in range(NV)])))-tmp)
		Der.append(t)
# arr is the matrix constructed using the derivatives in Der
arr=matrix(GF(p),Der)
for j in range(3):
		try:
			ret=arr.solve_left(vector(generate_list_form_of_linear(derQ1withh1[j])))
			print "derQ1withh1",j,ret
		except ValueError:
			print j+"th linear poly. not in span"

for j in range(3):
		try:
			ret=arr.solve_left(vector(generate_list_form_of_linear(derQ1withh2[j])))
			print "derQ1withh2",j,ret
		except ValueError:
			print j+"th linear poly. not in span"

for j in range(3):
		try:
			ret=arr.solve_left(vector(generate_list_form_of_linear(derQ2withh1[j])))
			print "derQ2withh1",j,ret
		except ValueError:
			print j+"th linear poly. not in span"

for j in range(3):
		try:
			ret=arr.solve_left(vector(generate_list_form_of_linear(derQ2withh2[j])))
			print "derQ2withh2",j,ret
		except ValueError:
			print j+"th linear poly. not in span"


# LET:-
# R1=doublederQ1[0]*(derQ2withh1[0]+derQ2withh2[0]) + doublederQ2[0]*(derQ1withh1[0]+derQ1withh2[0])+doublederQ1[0]*doublederQ2[0]
# R2=doublederQ1[1]*(derQ2withh1[1]+derQ2withh2[1]) + doublederQ2[1]*(derQ1withh1[1]+derQ1withh2[1])+doublederQ1[1]*doublederQ2[1]
# R3=doublederQ1[2]*(derQ2withh1[2]+derQ2withh2[2]) + doublederQ2[2]*(derQ1withh1[2]+derQ1withh2[2])+doublederQ1[2]*doublederQ2[2]

# WE SHOWED ABOVE THAT ALL OF THESE SINGLE DERIVATIVES can be recovered. Then, it is easy to recover R1,R2 and R3 since we explicitly know hi's.

R1=doublederQ1[0]*(derQ2withh1[0]+derQ2withh2[0]) + doublederQ2[0]*(derQ1withh1[0]+derQ1withh2[0])+doublederQ1[0]*doublederQ2[0]
R2=doublederQ1[1]*(derQ2withh1[1]+derQ2withh2[1]) + doublederQ2[1]*(derQ1withh1[1]+derQ1withh2[1])+doublederQ1[1]*doublederQ2[1]
R3=doublederQ1[2]*(derQ2withh1[2]+derQ2withh2[2]) + doublederQ2[2]*(derQ1withh1[2]+derQ1withh2[2])+doublederQ1[2]*doublederQ2[2]

print "R1",R1
print "R2",R2
print "R3",R3

R=[R1,R2,R3]
T=[0]*3
for j in range(3):
	T[j]=derprod[j]-(derQ1withh1[j]*derQ2withh2[j]+derQ1withh2[j]*derQ2withh1[j]+R[j])

# CAN WE RECOVER Q1 AND Q2
flag1=0;flag2=0
for i in range(p^3):
	i_inbasep=dectobasep(i,3)
	tmp=i_inbasep[0]*T[0]+i_inbasep[1]*T[1]+i_inbasep[2]*T[2]
	if(flag1==0 and tmp==Q1):
		flag1=1
		print "Q1 is recoverable"
	if(flag2==0 and tmp==Q2):
		flag2=1
		print "Q2 isrecoverable"


