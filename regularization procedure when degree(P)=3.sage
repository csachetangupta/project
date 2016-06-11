# Polynomials are represented in the list form which is defined below.
# Specifically, the term list form means the list of the coefficients of the monomials in lexical order with decreasing degree.
# Ex. If NV=2, then polynomial P of degree 3 is represented by the list of the coefficients of the monomials in the order [x0*x0*x0, x0*x0*x1, x0*x1*x0, x0*x1*x1, x1*x0*x0, x1*x0*x1, x1*x1*x0, x1*x1*x1, x0*x0,x0*x1,x1*x0,x1*x1,x0,x1,constant].
from datetime import datetime
import random
import numpy
from operator import add
from itertools import chain
print "Program starts at: ",datetime.now().time()

def dectobasep(num,L):
          alpha=[0]*L;len=L-1
          while(num):
              alpha[len]=(int)(num%p)
              num=(int)(num/p)
              len-=1
          return alpha


def BVLemma(Poly,k,MaxDegree,NV):
# Suppose DP is the polynomial obtained after taking iterative derivative MaxDegree times on the input polynomial "Poly".
# For instance if Poly= x[i]*x[j]*x[k], then DP=h1i*h2j*h3k+h1i*h3j*h2k+h2i*h1j*h3k+h2i*h3j*h1k+h3i*h2j*h1k+h3i*h1j*h2k.
# To identify variable "hai", we can use its index(hai) which is defined as (a-1)*NV+i.
# Below we have defined the variables which we will use:-
# CO will hold down the coefficients of the monomials of the DP.
# Ind is the collection of MaxDegree number of lists.
# Product of the variables defined by the Ind[i][j] over all i and for a particular j corresponds to a monomial of the DP.
# Ind[i] holds the indices of the i^th variable of the monomials of the DP.
# As instance for given above P, DP can also be represented by Ind and CO as defined below:-
# Ind[1]=[indices of {h1i,h1i,h2i,h2i,h3i,h3i}]=[i,i,NV+i,NV+i,2*NV+i,2*NV+i]
# Ind[2]=[indices of {h2j,h3j,h1j,h3j,h1j,h2j}]=[NV+j,2*NV+j,j,2*NV+j,j,NV+j]
# Ind[3]=[indices of {h3k,h2k,h3k,h1k,h2k,h1k}]=[2*NV+k,NV+k,2*NV+k,k,NV+k,k]
# and CO=[1,1,1,1,1,1].

#*******************PART 1******************
# Part1 generates the variables Ind and CO.
          #print "BVL starts",datetime.now().time()
          Ind=[[] for i in range(MaxDegree)]
          FactofMD=MaxDegree*(MaxDegree-1)
          CO=[]
          for i in range(len(Poly)):
              if(Poly[i]!=0):
                  tmp=i;L=[]
# L contains the indices of the variables of the monomial defined by index i.
# L=[i,j,k] for monomial=x[i]*x[j]*x[k] which was positioned at i*NV^2+j*NV+k in list Poly.
                  for j in range(MaxDegree):
                      L.append(int(tmp%NV))
                      tmp=int(tmp/NV)
                # CALCULATION OF MONOMIALS OF DP
                  if(MaxDegree==2):
                        if(L[0]!=L[1]):
                            Ind[0].append(L[0]);Ind[1].append(NV+L[1])
                            Ind[0].append(NV+L[0]);Ind[1].append(L[1])
                            CO+=[Poly[i]]*FactofMD
                        if(L[0]==L[1]):
                            Ind[0].append(L[0]);Ind[1].append(NV+L[0])
                            CO+=[2*Poly[i]]
                  if(MaxDegree==3):
                       A=L[0];B=L[1];C=L[2]
                       if(A==B and B==C):
                           Ind[0].append(A);Ind[1].append(NV+A);Ind[2].append(2*NV+A)
                           CO+=[6*Poly[i]]
                       if(A==B and B!=C):
                           Ind[0].append(A);Ind[1].append(NV+A);Ind[2].append(2*NV+C)
                           Ind[0].append(C);Ind[1].append(NV+A);Ind[2].append(2*NV+A)
                           Ind[0].append(A);Ind[1].append(NV+C);Ind[2].append(2*NV+A)
                           CO+=[2*Poly[i]]*MaxDegree
                       if(B==C and A!=C):
                           Ind[0].append(A);Ind[1].append(NV+B);Ind[2].append(2*NV+B)
                           Ind[0].append(B);Ind[1].append(NV+A);Ind[2].append(2*NV+B)
                           Ind[0].append(B);Ind[1].append(NV+B);Ind[2].append(2*NV+A)
                           CO+=[2*Poly[i]]*MaxDegree
                       if(A==C and B!=C):
                           Ind[0].append(B);Ind[1].append(NV+C);Ind[2].append(2*NV+C)
                           Ind[0].append(C);Ind[1].append(NV+B);Ind[2].append(2*NV+C)
                           Ind[0].append(C);Ind[1].append(NV+C);Ind[2].append(2*NV+B)
                           CO+=[2*Poly[i]]*MaxDegree
                       if(A!=C and B!=C and B!=A):
                           Ind[0].append(A);Ind[1].append(NV+B);Ind[2].append(2*NV+C)
                           Ind[0].append(A);Ind[1].append(2*NV+B);Ind[2].append(NV+C)
                           Ind[0].append(NV+A);Ind[1].append(B);Ind[2].append(2*NV+C)
                           Ind[0].append(NV+A);Ind[1].append(2*NV+B);Ind[2].append(C)
                           Ind[0].append(2*NV+A);Ind[1].append(B);Ind[2].append(NV+C)
                           Ind[0].append(2*NV+A);Ind[1].append(NV+B);Ind[2].append(C)
                           CO+=[Poly[i]]*FactofMD

#****************PART 2*********************
# Part2 uses the lists defined in part 1 and calculates the random derivatives of the DP.
# If P=2*x0*x1*x2, then from part 1 we have CO will be [2,2,2,2,2,2] with corresponding monomials MONOM=[h10*h21*h32, h10*h22*h31, h20*h31*h12, h20*h32*h12, h30*h11*h22, h30*h21*h12] which can also be defined as Ind1=[0,0,NV,NV,2*NV,2*NV], Ind2=[NV+1,2*NV+1,1,2*NV+1,1,NV+1], Ind3=[2*NV+2,NV+2,2*NV+2,2,NV+2,2].
# Below we have defined the variables which we have used:-
# Random list r specifies the random values of the variables of the DP.
# Random value of hi[j] will at position (i-1)*NV+j in list r. So, random_val(hi[j]) = r[(i-1)*NV+j].
# tmp is the collection of the random derivatives of the DP.
# constant corresponds to the constant term of the random derivative of the DP.
# a,b,c correponds to the values of the variables of a particular monomial of DP in random list r. Example if monomial=h10*h21*h32, then A=0,B=NV+1,C=2*NV+2 and a,b,c are their corresponding values in the random list r.
          #print "BVL middle",datetime.now().time()
          tmp=[]
          if(MaxDegree==2):
              for j in range(k):
                    constant=0
                    r=[random.randrange(p) for i in range(MaxDegree*NV)]
                    t=[0]*(MaxDegree*NV)
# t refers to the list form of the random derivative of the DP without constant term.
                    for i in range(len(Ind[0])):
                        A=Ind[0][i];B=Ind[1][i];COEF=CO[i]
                        a=r[A];b=r[B]
                        t[A]+=COEF*b
                        t[B]+=COEF*a
                        constant+=COEF*(a*b)
                    t.append(constant)
                    t=numpy.array(t)%p
                    tmp.append(t)
          if(MaxDegree==3):
              for j in range(k):
                    constant=0
                    r=[random.randrange(p) for i in range(MaxDegree*NV)]
                    t2=numpy.zeros(shape=(MaxDegree*NV,MaxDegree*NV),dtype=numpy.int)
                    t1=numpy.zeros(NV*MaxDegree,dtype=numpy.int)
# t2 and t1 hold degree 2 and 1 monomial coefficients of the random derivative of the DP respectively.
                    for i in range(len(Ind[0])):
                        A=Ind[0][i];B=Ind[1][i];C=Ind[2][i];COEF=CO[i]
                        a=r[A];b=r[B];c=r[C]
                        t2[A][B]+=COEF*c;t2[B][C]+=COEF*a;t2[A][C]+=COEF*b
                        t1[A]+=COEF*b*c
                        t1[B]+=COEF*a*c
                        t1[C]+=COEF*a*b
                        constant+=COEF*(a*b*c)
                    t2=numpy.triu(t2,0)+transpose(numpy.tril(t2,-1))
                    t2=t2.flatten()
                    t2=numpy.append(t2,numpy.append(t1,constant))
                    t2=t2%p
                    tmp.append(t2)
          #print "BVL ends",datetime.now().time()
          return tmp

def calculateDegrees(PL,NV):
# PL is the list of polynomials.
# DL represents the list of the degrees of the polynomials in the PL.
# Assumed min. degree is 1 for each polynomial.
         DL=[1]*len(PL)
         for j in range(len(PL)):
             if(len(PL[j])>=NV^3):
                  if(any(PL[j][NV^3:NV^3+NV^2])):
                      DL[j]=2
                  if(any(PL[j][:NV^3])):
                      DL[j]=3
             else:
		  if(len(PL[j])>=NV^2):
                  	if(any(PL[j][:NV^2])):
                      	    DL[j]=2
         return DL

def getPermutations(L,NV):
# getPermutations function returns sum of the product of the lists over all the possible permutations of the lists in L.
# For example if L=[L1,L2,L3] then getPermutations returns flattened value of M=L1*L2*L3+L1*L3*L2+ L2*L1*L3+L2*L3*L1+ L3*L1*L2+L3*L2*L1.
         if(len(L)==1):
             return L[0]
         ret=matrix(GF(p),[[0]*(NV^(len(L)-1))]*NV)
         for i in range(len(L)):
             ret+=transpose(L[i])*getPermutations(L[:i]+L[i+1:],NV)
         return matrix(list(chain.from_iterable(ret)))


def CalculateGowersNorm(POLY,MaxDegree,NV):
# GN equals to the gower's norm value of the polynomial POLY.
# Variable "List" can be viewed as random values of the variables of the polynomial DP as defined in the function BVLemma.
# L contains the MaxDegree number of lists having equal length extracted from above random "List".
# As an example if L=[L1,L2] then getPermutations returns L1*L2+L2*L1.
# SOL is a 1-dimensional matrix of size NV^MaxDegree obtained by flattening the MaxDegree-dimensional matrix M which equals to the sum of all possible permutations of the lists in L.
# The value of M[i][j][k] is equal to the L1[i].L2[j].L3[k]+L1[i].L2[k].L3[j]+L1[j].L2[i].L3[k]+L1[j].L2[k].L3[i]+L1[k].L2[i].L3[j]+L1[k].L2[j].L3[i]
# After flattening M, value of M[i][j][k] gets stored at position i*NV^2+j*NV+k in SOL.
# DPValueAtL is the value of the DP with list L as its variables values. Why DPValueAtL equals to (POLY*transpose(SOL))[0][0] (=dotproduct(POLY,SOL))? The reason is given in the second point of the "Non-trivial tricks used in the program code for efficient computation" section of the paper.
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


def CheckUniformityOfFunctions(factor,Epsilon,NV):
# This function ascertains the "Epsilon"-uniformity of the "factor".
# NV specifies the number of variables on which the polynomials of the "factor" are defined.
# Deglist is the list of the degrees of the polynomials present in the "factor".
# c (coefficient list of factor) is the list of coefficients of the polynomials in the "factor".
# MaxDegree holds the maximum degree among the polynomials having non zero coefficient value c[i].
# Variable "index" is equal to the index of the first polynomial in the "factor" having non zero c[i] value as well as deg=MaxDegree.

          #print "CUOF starts",datetime.now().time()
          Len=len(factor)
          Deglist=calculateDegrees(factor,NV)
          if(Len>0 and max(Deglist) > 1):
              K=1
              while(K<NoofCombinationsCheck):
                  c=dectobasep(random.randint(1,p^Len-1),Len)
                  K+=1
                  MaxDegree=0;index=0
                  for i in range(Len):
                      if(c[i]!=0):
                          if(MaxDegree<Deglist[i]):
                               MaxDegree=Deglist[i]
                               index=i
                  if(MaxDegree>1):
		      Poly=numpy.zeros(NV^MaxDegree,dtype=numpy.int)
                      PLD=numpy.zeros((int)((NV^(MaxDegree)-1)/(NV-1)),dtype=numpy.int)
# PLD holds the coefficients of the monomials of degree strictly less than the MaxDegree (i.e. beyond NV^MaxDegree).
# Poly gets assigned to the list of the coefficients of the monomials of degree exactly equal to MaxDegree (i.e. till NV^MaxDegree).
                      for i in range(Len):
                          if(c[i]!=0):
                              if(Deglist[i]==MaxDegree):
                                  tmp=factor[i][:NV^(MaxDegree)]
                                  Poly=Poly+c[i]*tmp
                                  tmp=factor[i][NV^(MaxDegree):]
                                  PLD=PLD+c[i]*tmp
                              else:
                                  tmp=numpy.append(numpy.zeros(len(PLD)-len(factor[i]),dtype=numpy.int),c[i]*factor[i])
                                  PLD=PLD+tmp
                      PLD=PLD%p
                      Poly=(Poly%p).tolist()
                      #print "Before GN",datetime.now().time()
                      GN=CalculateGowersNorm(matrix(GF(p),Poly),MaxDegree,NV)
                      #print "GN",GN
                      if(GN>Epsilon):
                               #print "Gowers Norm:",GN
                               #print "CUOF ends",datetime.now().time()
                               return c,Poly,PLD,index,MaxDegree
              #print "CUOF ends",datetime.now().time()
              return ([0]*Len),0,0,0,0
          else:
              return ([0]*Len),0,0,0,0


def MakeLinearIndependent(PL):
# PL is the collection of polynomials for which linear independent set has to be found.
# Used rank of the matrix for determining the linear dependence of the polynomial.
# The returning list PL will have the linear independent set of the PL.
                          maxLen=0
                          for i in range(len(PL)):
                              if(maxLen<len(PL[i])):
                                  maxLen=len(PL[i])
                          list=[]
                          i=len(PL)-1
                          while(i>=0):
                                   list.append([0]*(maxLen-len(PL[i]))+PL[i].tolist())
                                   rank=(matrix(GF(p),list)).rank()
                                   if(rank<len(list)):
                                        list.pop()
                                        del PL[i]
                                   i=i-1
                          return PL

def  UniformRefineOneStep(factor,Epsilon,k,NV):
# This function aims to make uniform the factor by a single step. The single step comprises of checking the "Epsilon" uniformity of the "factor". If it is already uniform, then return. Else estimate the polynomial Poly by random derivatives say D_1,....D_k by function BVLemma and again make uniform the factor defined by D_1,...D_k. It proceeds by lifting back the obtained uniform factor of the factor D_1,...D_k in terms of original variables.
# c,Poly and PLD are defined in function CheckUniformityOfFunctions.
# The variable "List" will hold the uniform factor output.
           #print "UROS starts",datetime.now().time()
           c,Poly,PLD,index,MaxDegree=CheckUniformityOfFunctions(factor,Epsilon,NV)
           List=[]
           if(any(c)):
                          P1=BVLemma(Poly,k,MaxDegree,NV)
                          P1=MakeLinearIndependent(P1)
                          del factor[index]
# In the below function call statement MaxDegree*NV will be our new NV.
                          Poly1=UniformRefine(P1,k,Epsilon,MaxDegree*NV)
                          for j in range(len(Poly1)):
                              Lin=Poly1[j][-(MaxDegree*NV+1):]
                              T1=numpy.zeros(NV,dtype=numpy.int)
# Below statement will bring back the derivatives whose variables are hi's to x's.
                              for l in range(MaxDegree):
                                  T1=T1+Lin[l*NV:(l+1)*NV]
                              T1=T1%p
                              P=numpy.append(T1,Lin[MaxDegree*NV])
                              temp_Poly=Poly1[j][:-(MaxDegree*NV+1)]
			       # lift back quadratic form if it exists
                              if((MaxDegree*NV)^2<=len(temp_Poly)):
                                  Quad=numpy.reshape(temp_Poly[-(MaxDegree*NV)^2:],(MaxDegree*NV,MaxDegree*NV))
                                  T2=numpy.zeros((NV,NV),dtype=numpy.int)
                                  for l in range(MaxDegree):
                                     for m in range(MaxDegree):
                                        T2=T2+Quad[l*NV:(l+1)*NV,m*NV:(m+1)*NV]
                                  T2=(T2.flatten())%p
                                  P=numpy.append(T2,P)
                                  temp_Poly=temp_Poly[:-(MaxDegree*NV)^2]
                              List.append(P)
                          List.append(PLD)
           #print "UROS ends",datetime.now().time()
           return c,factor,List


def UniformRefine(IPoly,k,Epsilon,NV):
# IPoly denotes the factor to be regularized.
# Below we have defined the variables, we used :-
# The role of the variable V is defined in the first point of the "Non-trivial tricks used in the program code for efficient computation" section of the paper.
# c is defined in function CheckUniformityOfFunctions.
# List represents the regularized version of the low degree factor defined by the derivatives D_1,....D_k (and D_1,....D_k are defined more properly in the function UniformRefineOneStep).
        V=[]
        while True:
          c,IPoly,List=UniformRefineOneStep(IPoly,Epsilon,k,NV)
          if(any(c)):
              V+=List
              continue
          else:
              if(V!=[]):
                  IPoly=IPoly+V
                  V=[]
                  continue
              return IPoly


def ConvertDeg3PolytoListform(Poly):
# MONOM holds all the monomials of Poly having non-zero coefficient value.
# COEFF holds the corresponding coefficient values.
# More clearly:- dotProduct(COEFF,MONOM)=Poly
# L will be the list form of the polynomial Poly.
# If Monomial=x1.x3.x4 then it's coefficient will be stored at position 1*(NV^2)+3*NV+4.
# If Monomial=x1.x6 then it's coefficient will be stored at position NV^3+1*NV+6.
# constant will be placed at last that is at NV^3+NV^2+NV.
          MONOM=Poly.monomials()
          COEFF=Poly.coefficients()
          L=numpy.zeros(NV^3+NV^2+NV+1,dtype=numpy.int)
          for i in range(len(MONOM)):
              t=MONOM[i].variables()
              if any(t):
                  ret=0
                  for k in range(len(t)):
                      index=x.index(t[k])
                      for j in range(MONOM[i].degree(t[k])):
                          ret=ret*NV+index
                  if(MONOM[i].total_degree()==3):
                      L[ret]=COEFF[i]
                  if(MONOM[i].total_degree()==2):
                      L[NV^3+ret]=COEFF[i]
                  if(MONOM[i].total_degree()==1):
                      L[NV^3+NV^2+ret]=COEFF[i]
              else:
                  L[NV^3+NV^2+NV]=COEFF[i]
          return L



def GenerateUniformFactor(IPolylist,k,Epsilon,NV):
# IPolyList has the given degree3 polynomial in the list form.
# UFL contains the polynomials of the uniform factor in list form.
# UF denotes the independent set of UFL.
           #print k,Epsilon
           UFL=UniformRefine([IPolylist],k,Epsilon,NV)
           UF=MakeLinearIndependent(UFL)
           return UF


p=5;NV=20                                                       # p=Field Size,NV=number of variables 
SampleSize=p^2+75                                               # "SampleSize"=number of samples to be used for calculating gower's norm
NoofCombinationsCheck=25				        # No_Of_Combinations_Check specifies how many coefficients lists to check to ascertain Epsilon uniformity of a factor.
print "No of variables=",NV,"FieldSize",p
print "SampleSize",SampleSize,"NoofCombinationsCheck",NoofCombinationsCheck

# Defining Polynomial ring
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
Qx=Qx+x+[1]
Qxv=vector(Qx)


def GenerateRandomCP():
# Generate a random cubic polynomial of the form pol1*pol2+pol3*pol4 
# where deg(pol1,pol3)=2 and deg(pol2,pol4)=1.
# IPolylist is the cubic polynomial in list(numpy.array) form.
# Let the term 'generator polynomials' specifies the polynomials pol1,pol2,pol3 and pol4.
        pol1=Qxv.dot_product(vector([random.randrange(p) for i in range(NV^2+NV+1)]))
        pol2=xv.dot_product(vector([random.randrange(p) for i in range(NV+1)]))
        pol3=Qxv.dot_product(vector([random.randrange(p) for i in range(NV^2+NV+1)]))
        pol4=xv.dot_product(vector([random.randrange(p) for i in range(NV+1)]))
        IPolylist=ConvertDeg3PolytoListform(pol1*pol2+pol3*pol4)
        return pol1, pol2, pol3, pol4, IPolylist

def CheckInSpan(M, Quad):
# This function checks that does "Quad" lies in the span of elements of M or not.
# Here, M will represent the uniform factor in list form and Quad is the polynomial in algebraic form.
	vec = [0]*(NV^2+NV+1)  
        for i in range(len(M)):
            Len=len(M[i])
	    M[i]=list(M[i])
            if(Len>len(vec)):		#if M contains a polynomial of degree greater than 2
                return 0
            else:
                if(Len!=len(vec)):
                    M[i]=[0]*(len(vec)-Len)+M[i]
        Monom = Quad.monomials()
        Coeffs = Quad.coefficients()
        for i in range(len(Monom)):
             t=Monom[i].variables()
             if any(t):
                 ret=0
                 for k in range(len(t)):
                     index=x.index(t[k])
                     for j in range(Monom[i].degree(t[k])):
                         ret=ret*NV+index
                 if(Monom[i].total_degree()==2):
                      vec[ret]=Coeffs[i]
                 if(Monom[i].total_degree()==1):
                      vec[NV^2+ret]=Coeffs[i]
             else:
                  vec[NV^2+NV] = Coeffs[i]
        return Matrix(GF(p), M).rank() == Matrix(GF(p), M + [vec]).rank()

def findEpsilon(newk,EpsilonS,EpsilonE,diff,total):
# It searches for Epsilon among the values defined by EpsilonS, EpsilonE and diff such that pol2 and pol4 can be written as a linear combination of the polynomials in the uniform factor.
# newk specifies the value of k(=#of derivatives to be evaluated in Bogdanov Viola Lemma).
# EpsilonS and EpsilonE specifies the starting and ending values of Epsilon.
# diff denotes the difference between two consecutive values of Epsilon.
# total specifies the number of times to execute the regularization procedure for any given Epsilon.
        Epsilon=EpsilonE
        k=newk
        while(Epsilon>=EpsilonS):
            bad=0
            for j in range(total):
                pol1, pol2, pol3, pol4, factor = GenerateRandomCP()
                M = GenerateUniformFactor(factor,k,Epsilon,NV)
                if not(CheckInSpan(M, pol2) and CheckInSpan(M,pol4)):
                    bad+=1
                DL=calculateDegrees(M,NV)
                print DL
            print bad,Epsilon
            Epsilon=Epsilon-diff

def RunChecksOverk(newk,Epsilon,total):
# This function finds the number of bad iterations out of total where bad iteration means that the generator polynomial cannot be written as a linear combination of the polynomials in the uniform factor.
# newk specifies the value of k(= # of derivatives to be evaluated in Bogdanov Viola Lemma).
# Epsilon specifies the value of epsilon to be used.
# total specifies the number of times to execute the regularization procedure.
        k = newk
        bad1=0
        bad2=0
        bad3=0
        bad4=0
        for i in range(total):
            pol1, pol2, pol3, pol4, factor = GenerateRandomCP()
            M = GenerateUniformFactor(factor,k,Epsilon,NV)
            if not(CheckInSpan(M, pol1)):
                bad1 += 1
            if not(CheckInSpan(M, pol2)):
                bad2 += 1
            if not(CheckInSpan(M, pol3)):
                bad3 += 1
            if not(CheckInSpan(M, pol4)):
                bad4 += 1
            DL=calculateDegrees(M,NV)
            print "Degree list",DL,len(DL)
        print "k",k
        print "bad1",bad1,"bad2",bad2,"bad3",bad3,"bad4",bad4

def PerfOverRandomPoly(newk,EpsilonS,EpsilonE,diff,total):
# This function generates the list of the lengths of the uniform factor for all Epsilon's, total number of times.
# newk specifies the value of k(=#of derivatives to be evaluated in Bogdanov Viola Lemma).
# EpsilonS and EpsilonE specifies the starting and ending values of Epsilon.
# diff denotes the difference between two consecutive values of Epsilon.
# Here, total specifies the number of times to execute the regularization procedure.
          for i in range(total):
             print "k",newk
             IPolylist=numpy.array([random.randrange(p) for j in range(NV^3+NV^2+NV+1)])
             Epsilon=EpsilonS;k=newk
             print "Input polynomial",IPolylist
             UFLlist=[]
             while(Epsilon<=EpsilonE):
                 UFL=GenerateUniformFactor(IPolylist,k,Epsilon,NV)
                 print Epsilon,"length of uniform factor",len(UFL)
                 Epsilon+=diff
                 if(len(UFL)>NV+1):
                     print calculateDegrees(UFL,NV)
                 UFLlist.append(len(UFL))
             print "Length of uniform factors",UFLlist

RunChecksOverk(8,0.68,5)
#PerfOverRandomPoly(6,0.55,0.91,0.05,3)
print "program ends at:",datetime.now().time()









