# This code will run on polynomial P=L*Z  where L and Z are randomly chosen polynomials of degrees 1 and atmost 4 respectively and try to get back the polynomials L and Z. 

from datetime import datetime
import random
import numpy
from operator import add
from itertools import chain
p=5;NV=12                                                        # p=Field Size,NV=number of variables 
SampleSize=p^2+200                                               # "SampleSize"=number of samples to be used for calculating gower's norm
print "No of variables=",NV,"FieldSize",p
print "SampleSize",SampleSize

#Defining Polynomial ring
R=PolynomialRing(GF(p), NV, var_array=['x'])
R.inject_variables()
x=(list)(R.gens())

Qx=[]
for i in range(NV):
          j=0
          while(j<NV):
               Qx.append(x[i]*x[j])
               j+=1


def dectobasep(num,L):
          alpha=[0]*L;len=L-1
          while(num):
              alpha[len]=(int)(num%p)
              num=(int)(num/p)
              len-=1
          return alpha

def getPermutations(L,NV):
# getPermutations function returns sum of product of all the possible permutations of the lists in L.
# For example if L=[L1,L2,L3] then getPermutations returns L1*L2*L3+L1*L3*L2+ L2*L1*L3+L2*L3*L1+ L3*L1*L2+L3*L2*L1.
         if(len(L)==1):
             return L[0]
         ret=matrix(GF(p),[[0]*(NV^(len(L)-1))]*NV)
         for i in range(len(L)):
             ret+=transpose(L[i])*getPermutations(L[:i]+L[i+1:],NV)
         return matrix(list(chain.from_iterable(ret)))


def CalculateGowersNorm(POLY,MaxDegree,NV):
# GN equals to the gower's norm value of polynomial POLY.
# Variable "List" can be viewed as random values of variables of the polynomial DP as defined in the function BVLemma.
# L contains the MaxDegree number of lists extracted from above random "List".
# As an example if L=[L1,L2] then getPermutations returns L1*L2+L2*L1.
# SOL is a 1-dimensional matrix of size NV^MaxDegree obtained by flattening the MaxDegree-dimensional matrix M which equals to the sum of all possible permutations of the lists in L.
# The value of M[i][j][k] is equal to the L1[i].L2[j].L3[k]+L1[i].L2[k].L3[j]+L1[j].L2[i].L3[k]+L1[j].L2[k].L3[i]+L1[k].L2[i].L3[j]+L1[k].L2[j].L3[i]
# After flattening M, value of M[i][j][k] gets stored at position i*NV^2+j*NV+k in SOL.
# DPValueAtL is the value of the DP with list L as its variables values. Why DPValueAtL equals to (POLY*transpose(SOL))[0][0] (=dotproduct(POLY,SOL))? The reason is given in second point of the "Non-trivial tricks used in the program code for efficient computation" section of the paper.
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

def checkinspan(Der,l1):
# It checks that does linear polynomial l1 lies in the span of the polynomials in the list Der.
# Der contains the polynomials in algebraic form and l1 is the list form of L.
        List=[]
	for i in range(len(Der)):
		t=[0]*(NV)
		tmp=Der[i]
		for j in range(NV):
			t[j]=tmp.monomial_coefficient(x[j])
		List.append(t)
	try: 
		M=matrix(GF(p),List)						
		ret=M.solve_left(l1)
		print "l1 is in span",ret
		return true
	except ValueError:
		z=1
	return	false

def GenerateListForm(poly,degree):
# It will return the list form of the polynomial poly which consists of the coefficients of monomials of degree "degree" only.
# ret_list represents the list form to be returned.
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

def Factor(deg_of_sec_poly,l1,L,Z,Px):
# It generates the random derivatives of Px, finds the linear combination of the derivatives which has highest Gower's norm and recursively call Factor over new polynomial.

	# Px is the polynomial of degree deg_of_sec_poly+1 in the terms of x's.
	# Der will have collection of the random derivatives of Px.
	Der=[]
	# If deg_of_sec_poly==1 then generate the derivatives and check does linear polynomial l1 lies in the span of derivatives or not. 
	if(deg_of_sec_poly==1):
		for j in range(4):				# 4 represents the number of derivatives computed
                    tmp=Px
		    r=[random.randrange(p) for i in range(NV)]
                    tmp=tmp.subs(dict(zip(x,[x[i]+r[i] for i in range(NV)])))-tmp
                    Der.append(tmp)	
		return checkinspan(Der,l1)		
		

	# Else also generate the random derivatives
	# r in the random list of values of variables. In every iteration, we have randomly chosen r.
	# c[j] is equal to the derivative of L in random direction r chosen
	# dz[j] is equal to the derivative of Z in random direction r chosen
	c=[0]*random_deri
	dz=[0]*random_deri
	for j in range(random_deri):
		    tmp=Px
		    r=[random.randrange(p) for i in range(NV)]	
                    tmp=tmp.subs(dict(zip(x,[x[i]+r[i] for i in range(NV)])))-tmp
                    Der.append(tmp)
		    tmp=L
		    c[j]=tmp.subs(dict(zip(x,[x[i]+r[i] for i in range(NV)])))-tmp
		    tmp=Z
		    dz[j]=tmp.subs(dict(zip(x,[x[i]+r[i] for i in range(NV)])))-tmp
		    
		    #if(Der[j]==c[j]*Z+(L+c[j])*dz[j]):
		#	print "true"
		#    else:
		#	print "false"


	# the below nested loops find the linear combination of the derivatives which has the highest gowers norm
	# In the below loop, tmp be the linear combination of the derivatives.
	# R_Poly is equal to zero makes sure that tmp does not contain polynomial Z.
	# poly wil have the value equal to sum over j of i_inbasep[j]*(L+c[j])*dz[j] which is a part of polynomial Der[j].
	si=[0]*random_deri;mgn=0
	i=1
	while(i<p^random_deri):
		i_inbasep=dectobasep(i,random_deri)
		tmp=0;poly=0
		for j in range(random_deri):
			tmp=tmp+i_inbasep[j]*Der[j]
			poly=poly+i_inbasep[j]*(L+c[j])*dz[j]
		ListForm=GenerateListForm(tmp,deg_of_sec_poly)
		gn=CalculateGowersNorm(matrix(GF(p),ListForm),deg_of_sec_poly,NV)
		R_Poly=tmp-poly
		if(gn>mgn and R_Poly==0):
			si=i_inbasep;mgn=gn
		i=i+1
	
	#calculate new Z
	Z=0
	for j in range(random_deri):
		Z=Z+si[j]*dz[j]
	Monom=Z.monomials()
	Coeffs=Z.coefficients()
	Z=0
	for j in range(len(Monom)):
		if(Monom[j].total_degree()==deg_of_sec_poly-1):
			Z=Z+Coeffs[j]*Monom[j]
	

	# tmp refers the linear combination of the derivatives present in Der which has the highest gowers norm among all non-zero linear combinations
	tmp=0
	for j in range(random_deri):
		tmp=tmp+si[j]*Der[j]

	# Calculate new Px which will be equal to POLY
	Monom=tmp.monomials()
	Coeffs=tmp.coefficients()
	POLY=0
	for j in range(len(Monom)):
		if(Monom[j].total_degree()==deg_of_sec_poly):
			POLY=POLY+Coeffs[j]*Monom[j]
	Factor(deg_of_sec_poly-1,l1,L,Z,POLY)
	

cubic_terms=list(chain.from_iterable(transpose(matrix(x))*matrix(Qx)))
quartic_terms=list(chain.from_iterable(transpose(matrix(x))*matrix(cubic_terms)))
random_deri=2		# random_deri represents the number of random derivatives to be computed
deg_of_sec_poly=4	# deg_of_sec_poly is the degree of the second polynomial Z.
print "random_deri",random_deri
print "degree of second polynomial",deg_of_sec_poly

# The function Exper below generates the polynomial P=L.Z randomly where degrees of L and Z are 1 and deg_of_sec_poly repsectively and try to regain back the L and Z.
# The code works for degree atmost 4. It is easy to modify the code for higher degrees with suitable Field Size p by changing the code in function Exper.
def Exper(deg_of_sec_poly):		
	l1=vector([random.randrange(p) for i in range(NV)])
        l2=vector([random.randrange(p) for i in range(NV^deg_of_sec_poly)])
	L=vector(x).dot_product(l1)
	if(deg_of_sec_poly==1):
		Z=vector(x).dot_product(l2)
	if(deg_of_sec_poly==2):
		Z=vector(Qx).dot_product(l2)
	if(deg_of_sec_poly==3):
		Z=vector(cubic_terms).dot_product(l2)
	if(deg_of_sec_poly==4):
		Z=vector(quartic_terms).dot_product(l2)	
	#P=transpose(l1)*l2
	#P=numpy.array(P).flatten()%p
	Px=Z*L
	#print Px
	#print P,len(P)
	#print "l1,l2",l1,l2

	
	# Px will be the polynomial of the form L*Z.
	# l1 and l2 are list forms of polynomials L and Z respectively.
	# deg_of_sec_poly is the degree of the second polynomial Z.
	Factor(deg_of_sec_poly,l1,L,Z,Px)


for l in range(10):
	print l
	Exper(deg_of_sec_poly)

