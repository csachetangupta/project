# Polynomials are represented in the list form.
# Specifically, the list form is the list of the coefficients of the monomials in lexical order with decreasing degree.
# Ex. If NV=2=number of variables, then polynomial P of degree 3 is represented as coefficients of monomials in the order [x0*x0*x0, x0*x0*x1, x0*x1*x0, x0*x1*x1, x1*x0*x0, x1*x0*x1, x1*x1*x0, x1*x1*x1, x0*x0,x0*x1,x1*x0,x1*x1,x0,x1,constant].
# The only difference between degree 3 and degree 4 codes is in the polynomial representation. If NV=2, then any polynomial of degree atmost 1 has form like [0,0,0,0,coeff(x0),coeff(x1),constant] in degree3 code. However, in degree 4 code it's form will be like [coeff(x0),coeff(x1),constant]. So, in degree 4 code polynomials have more precise represntation as compared to in degree 3.
from datetime import datetime
import random
import copy
import numpy
from operator import add
from itertools import chain
print "Program starts at: ",datetime.now().time()


def dectobasep(num,L):
#num is the decimal number to be converted into base p
#L corresponds to the maximum length of num in base p
          num_in_base_p=[0]*L;len=L-1
          while(num):
              num_in_base_p[len]=(int)(num%p)
              num=(int)(num/p)
              len-=1
          return num_in_base_p

def  CalculateFactorial(d):
        if(d==2):
            return 2
        return d*CalculateFactorial(d-1)


def B_VLemma(Poly,k,MaxDegree,NV):
# Suppose IP is the polynomial obtained after taking iterative derivative MaxDegree times on the input polynomial "Poly".
# For instance if P= x[i]*x[j]*x[k], then IP=h1i*h2j*h3k+h1i*h3j*h2k+h2i*h1j*h3k+h2i*h3j*h1k+h3i*h2j*h1k+h3i*h1j*h2k.
# NV used below is defined in UniformRefine and UniformRefineOneStep.
# To identify variable "hai", we can use its index(hai) which is =(a-1)*NV+i.
# Below we have defined the variables which we will use:-
# CO will hold down the coefficients of the monomials of the IP.
# NEW is the collection of MaxDegree number of lists.
# Product of the variables defined by the NEW[i][j] over all i and for a particular j corresponds to a monomial of the IP.
# NEW[i] holds the index of the ith variable of the monomial of the IP.
# As instance for given above Poly, IP can also be represented by NEW and CO as defined below:-
# NEW[1]=[indices of {h1i,h1i,h2i,h2i,h3i,h3i}]=[i,i,NV+i,NV+i,2*NV+i,2*NV+i]
# NEW[2]=[indices of {h2j,h3j,h1j,h3j,h1j,h2j}]=[NV+j,2*NV+j,j,2*NV+j,j,NV+j]
# NEW[3]=[indices of {h3k,h2k,h3k,h1k,h2k,h1k}]=[2*NV+k,NV+k,2*NV+k,k,NV+k,k]
# and CO=[1,1,1,1,1,1].

#*******************PART 1******************
# Part1 evaluates the variables NEW and CO.
          #print "BVL starts",datetime.now().time()
          NEW=[[] for i in range(MaxDegree)]
          fact_of_md=CalculateFactorial(MaxDegree)
          CO=[]
          for i in range(len(Poly)):
              if(Poly[i]!=0):
                  tmp=i;L=[]
# L will contain the indices of the variables of the corresponding monomial defined by index i.
# L=[i,j,k] for monomial=x[i]*x[j]*x[k] which was positioned at i*NV^2+j*NV+k in list Poly.
                  for j in range(MaxDegree):
                      L.append(int(tmp%NV))
                      tmp=int(tmp/NV)
                #CALCULATION OF MONOMIALS OF IP
                  if(MaxDegree==2):
                           NEW[0].append(L[0]);NEW[1].append(NV+L[1])
                           NEW[0].append(NV+L[0]);NEW[1].append(L[1])
                           CO=CO+[Poly[i]]*fact_of_md
                  if(MaxDegree==3):
                           A=L[0];B=L[1];C=L[2]
                           NEW[0].append(A);NEW[1].append(NV+B);NEW[2].append(2*NV+C)
                           NEW[0].append(A);NEW[1].append(NV+C);NEW[2].append(2*NV+B)
                           NEW[0].append(B);NEW[1].append(NV+A);NEW[2].append(2*NV+C)
                           NEW[0].append(B);NEW[1].append(NV+C);NEW[2].append(2*NV+A)
                           NEW[0].append(C);NEW[1].append(NV+A);NEW[2].append(2*NV+B)
                           NEW[0].append(C);NEW[1].append(NV+B);NEW[2].append(2*NV+A)
                           CO=CO+[Poly[i]]*fact_of_md
                  if(MaxDegree==4):
                           A=L[0];B=L[1];C=L[2];D=L[3]
                           NEW[0].append(A);NEW[1].append(NV+B);NEW[2].append(2*NV+C);NEW[3].append(3*NV+D)
                           NEW[0].append(A);NEW[1].append(NV+B);NEW[2].append(2*NV+D);NEW[3].append(3*NV+C)
                           NEW[0].append(A);NEW[1].append(NV+C);NEW[2].append(2*NV+B);NEW[3].append(3*NV+D)
                           NEW[0].append(A);NEW[1].append(NV+C);NEW[2].append(2*NV+D);NEW[3].append(3*NV+B)
                           NEW[0].append(A);NEW[1].append(NV+D);NEW[2].append(2*NV+B);NEW[3].append(3*NV+C)
                           NEW[0].append(A);NEW[1].append(NV+D);NEW[2].append(2*NV+C);NEW[3].append(3*NV+B)
                           
                           NEW[0].append(B);NEW[1].append(NV+A);NEW[2].append(2*NV+C);NEW[3].append(3*NV+D)
                           NEW[0].append(B);NEW[1].append(NV+A);NEW[2].append(2*NV+D);NEW[3].append(3*NV+C)
                           NEW[0].append(B);NEW[1].append(NV+C);NEW[2].append(2*NV+A);NEW[3].append(3*NV+D)
                           NEW[0].append(B);NEW[1].append(NV+C);NEW[2].append(2*NV+D);NEW[3].append(3*NV+A)
                           NEW[0].append(B);NEW[1].append(NV+D);NEW[2].append(2*NV+A);NEW[3].append(3*NV+C)
                           NEW[0].append(B);NEW[1].append(NV+D);NEW[2].append(2*NV+C);NEW[3].append(3*NV+A)
                            
                           NEW[0].append(C);NEW[1].append(NV+A);NEW[2].append(2*NV+B);NEW[3].append(3*NV+D)
                           NEW[0].append(C);NEW[1].append(NV+A);NEW[2].append(2*NV+D);NEW[3].append(3*NV+B)
                           NEW[0].append(C);NEW[1].append(NV+B);NEW[2].append(2*NV+A);NEW[3].append(3*NV+D)
                           NEW[0].append(C);NEW[1].append(NV+B);NEW[2].append(2*NV+D);NEW[3].append(3*NV+A)
                           NEW[0].append(C);NEW[1].append(NV+D);NEW[2].append(2*NV+A);NEW[3].append(3*NV+B)
                           NEW[0].append(C);NEW[1].append(NV+D);NEW[2].append(2*NV+B);NEW[3].append(3*NV+A)
                            
                           NEW[0].append(D);NEW[1].append(NV+A);NEW[2].append(2*NV+B);NEW[3].append(3*NV+C)
                           NEW[0].append(D);NEW[1].append(NV+A);NEW[2].append(2*NV+C);NEW[3].append(3*NV+B)
                           NEW[0].append(D);NEW[1].append(NV+B);NEW[2].append(2*NV+A);NEW[3].append(3*NV+C)
                           NEW[0].append(D);NEW[1].append(NV+B);NEW[2].append(2*NV+C);NEW[3].append(3*NV+A)
                           NEW[0].append(D);NEW[1].append(NV+C);NEW[2].append(2*NV+A);NEW[3].append(3*NV+B)
                           NEW[0].append(D);NEW[1].append(NV+C);NEW[2].append(2*NV+B);NEW[3].append(3*NV+A)
                           CO=CO+[Poly[i]]*fact_of_md
#****************PART 2*********************
# Part2 uses the lists defined in part 1 and calculates the random derivatives of the IP.
# If Poly=2*x0*x1*x2, then from part 1 we have CO will be [2,2,2,2,2,2] with corresponding monomials MONOM=[h10*h21*h32, h10*h22*h31, h20*h31*h12, h20*h32*h12, h30*h11*h22, h30*h21*h12] which can also be defined as NEW1=[0,0,NV,NV,2*NV,2*NV], NEW2=[NV+1,2*NV+1,1,2*NV+1,1,NV+1], NEW3=[2*NV+2,NV+2,2*NV+2,2,NV+2,2].
# Below we have defined the variables which we have used:-
# Random list r specifies the random values of the variables of the IP.
# Random value of hi[j] will at position (i-1)*NV+j in list r. So, random_val(hi[j]) = r[(i-1)*NV+j].
# tmp will hold the random derivatives of the IP.
# constant will correspond to the constant term of the derivative of the IP.
# a,b,c correponds to the values of the variables of a particular monomial in r. Example if monomial=h10*h21*h32, then A=0,B=NV+1,C=2*NV+2 as defined above and a,b,c are their corresponding values in the random list r.
          tmp=[]
          #print "BVL mid",datetime.now().time()
          if(MaxDegree==2):
              for j in range(k):
                    constant=0
                    r=[random.randrange(p) for i in range(MaxDegree*NV)]
                    t=[0]*(MaxDegree*NV)
# t will be the list form of the random derivative
                    for i in range(len(NEW[0])):
                        A=NEW[0][i];B=NEW[1][i];COEF=CO[i]
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
# t2 and t1 hold degree 2 and 1 monomial coefficients of the random derivative of IP respectively.
                    for i in range(len(NEW[0])):
                        A=NEW[0][i];B=NEW[1][i];C=NEW[2][i];COEF=CO[i]
                        a=r[A];b=r[B];c=r[C]
                        t2[A][B]+=COEF*c;t2[B][C]+=COEF*a;t2[A][C]+=COEF*b
                        t1[A]+=COEF*b*c
                        t1[B]+=COEF*a*c
                        t1[C]+=COEF*a*b
                        constant+=COEF*(a*b*c)
                    t2=t2.flatten()
                    t2=numpy.append(t2,numpy.append(t1,constant))
                    t2=t2%p
                    tmp.append(t2)
          if(MaxDegree==4):
              for j in range(k):
                    constant=0
                    r=[random.randrange(p) for i in range(MaxDegree*NV)]
                    t3=numpy.zeros(shape=(MaxDegree*NV,MaxDegree*NV,MaxDegree*NV),dtype=numpy.int)
                    t2=numpy.zeros(shape=(MaxDegree*NV,MaxDegree*NV),dtype=numpy.int)
                    t1=numpy.zeros(NV*MaxDegree,dtype=numpy.int)
# t3,t2 and t1 hold degree 3, 2 and 1 monomial coefficients of the random derivative of IP respectively.
                    for i in range(len(NEW[0])):
                        A=NEW[0][i];B=NEW[1][i];C=NEW[2][i];D=NEW[3][i];COEF=CO[i]
                        a=r[A];b=r[B];c=r[C];d=r[D]
                        t3[A][B][C]+=COEF*d
                        t3[A][B][D]+=COEF*c
                        t3[A][C][D]+=COEF*b
                        t3[B][C][D]+=COEF*a
                        t2[A][B]+=COEF*c*d;t2[A][C]+=COEF*b*d;t2[A][D]+=COEF*b*c
                        t2[B][C]+=COEF*a*d;t2[B][D]+=COEF*a*c
                        t2[C][D]+=COEF*a*b
                        t1[A]+=COEF*b*c*d
                        t1[B]+=COEF*a*c*d
                        t1[C]+=COEF*a*b*d
                        t1[D]+=COEF*a*b*c
                        constant+=COEF*(a*b*c*d)
                    #t2=numpy.triu(t2,0)+transpose(numpy.tril(t2,-1))
                    t3=t3.flatten()
                    t2=t2.flatten()
                    t3=numpy.append(t3,t2)
                    t1=numpy.append(t1,constant)
                    t3=numpy.append(t3,t1)
                    t3=t3%p
                    tmp.append(t3)
          #print "BVL ends",datetime.now().time()
          return tmp

def getPermutations(L,NV):
# getPermutations function returns sum of product of all the possible permutations of the lists in L.
# For example if L=[L1,L2,L3] then getPermutations returns L1*L2*L3+L1*L3*L2+ L2*L1*L3+L2*L3*L1+ L3*L1*L2+L3*L2*L1.
         if(len(L)==1):
             return L[0]
         ret=matrix(GF(p),[[0]*(NV^(len(L)-1))]*NV)
         for i in range(len(L)):
             ret=ret+transpose(L[i])*getPermutations(L[:i]+L[i+1:],NV)
         return matrix(list(chain.from_iterable(ret)))



def CalculateGowersNorm(Poly,MaxDegree,NV):
# GN equals to the gower's norm value of polynomial POLY.
# List can be viewed as random values of variables of the polynomial IP obtained by taking repeated derivatives deg(POLY)(=MaxDegree) times of the POLY.
# L will generate the MaxDegree number of lists by above random "List".
# Variable "List" can be viewed as random values of variables of the polynomial IP as defined in the function BVLemma.
# L contains the MaxDegree number of lists extracted from above random "List".
# If L=[L1,L2] then getPermutations returns L1*L2+L2*L1.
# SOL is a 1-dimensional matrix of size NV^MaxDegree obtained by flattening the MaxDegree-dimensional matrix M which equals to the sum of all possible permutations of the lists in L.
# The value of M[i][j][k] is equal to the L1[i].L2[j].L3[k]+L1[i].L2[k].L3[j]+L1[j].L2[i].L3[k]+L1[j].L2[k].L3[i]+L1[k].L2[i].L3[j]+L1[k].L2[j].L3[i]
# After flattening M, value of M[i][j][k] gets stored at position i*NV^2+j*NV+k in SOL.
# IPValueAtL is the value of the IP with list L as its variables values. Why IPValueAtL equals to (POLY*transpose(SOL))[0][0] (=dotproduct(POLY,SOL))? The reason is very easy to figure out.
                      pthroot = cos(2*pi/p)+I*sin(2*pi/p)
                      GN=0
                      for j in range(SampleSize):
                          List = [random.randrange(p) for i in range(NV*MaxDegree)]
                          L=[]
                          for i in range(MaxDegree):
                              L.append(matrix(List[i*NV:(i+1)*NV]))
                          SOL=getPermutations(L,NV)
                          SOL=matrix(GF(p),SOL)
                          IP_Value_At_L=(Poly*transpose(SOL))[0][0]
                          GN=GN+pthroot^(int)(IP_Value_At_L)
                      GN=(float(abs(GN))*1/SampleSize)^(1/2^MaxDegree)
                      return GN



def CheckUniformityOfFunctions(factor,Epsilon,NV):
# This function ascertains the "Epsilon"-uniformity of the "factor".
# NV specifies the number of variables on which polynomials of the "factor" is defined.
# Deglist is the list of the degrees of the polynomials present in the "factor".
# c (coefficient list of factor) is the list of coefficients of the polynomials in the "factor".
# MaxDegree holds the maximum degree among the polynomials having non zero coefficient value c[i].
# Variable "index" equals to the index of the first polynomial in the "factor" having non zero c[i] value as well as deg=MaxDegree.

          #print "CUOF starts",datetime.now().time()
          LenOfFact=len(factor)
          DegreeList=calculateDegrees(factor,NV)
          if(LenOfFact>0 and max(DegreeList) > 1):
              K=1
              while(K<No_Of_Combinations_Check):
                  c=dectobasep(random.randint(1,p^LenOfFact-1),LenOfFact)
                  K+=1
                  MaxDegree=0;index=0
                  for i in range(LenOfFact):
                      if(c[i]!=0):
                          if(MaxDegree<DegreeList[i]):
                               MaxDegree=DegreeList[i]
                               index=i
                  if(MaxDegree>1):
                      Poly=numpy.zeros(NV^MaxDegree,dtype=numpy.int)
                      PLD=numpy.zeros((int)((NV^(MaxDegree)-1)/(NV-1)),dtype=numpy.int)
# PLD holds the coefficients of all those monomials whose degree is less than MaxDegree (i.e. beyond NV^MaxDegree).
# Poly gets reassigned to the list of the coefficients of all those monomials whose degree equal to MaxDegree (i.e. till NV^MaxDegree).
                      for i in range(LenOfFact):
                          if(c[i]!=0):
                              if(DegreeList[i]==MaxDegree):
                                  tmp=factor[i][:NV^(MaxDegree)]
                                  Poly=Poly+c[i]*tmp
                                  tmp=factor[i][NV^(MaxDegree):]
                                  PLD=PLD+c[i]*tmp
                              else:
                                  tmp=numpy.append(numpy.zeros(len(PLD)-len(factor[i]),dtype=numpy.int),c[i]*factor[i])
                                  PLD=PLD+tmp
                      PLD=PLD%p
                      Poly=(Poly%p).tolist()
                      GN=CalculateGowersNorm(matrix(GF(p),Poly),MaxDegree,NV)
                      if(GN>Epsilon):
                               return c,Poly,PLD,index,MaxDegree
              #print "CUOF ends",datetime.now().time()
              return ([0]*LenOfFact),0,0,0,0
          else:
              return ([0]*LenOfFact),0,0,0,0


def MakeLinearIndependent(PL):
# PL is the collection of poly's of which linear independent set has to be found.
# Used rank of the matrix for determining the dependence of the polynomial.
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
# This function aims to make uniform the factor by a single step. The single step comprise of checking the "Epsilon" uniformity of the "factor". If it is already uniform, then return. Else estimate the polynomial Poly by random derivatives say D_1,....D_k by function BVLemma and again make uniform the factor defined by D_1,...D_k. It proceeds by lifting back the obtained uniform factor of the factor D_1,...D_k in terms of original variables.
# c,Poly and PLD are defined in function CheckUniformityOfFunctions.
# The variable "list" will hold the uniform factor.

           #print "UROS starts",datetime.now().time()
           c,Poly,PLD,index,MaxDegree=CheckUniformityOfFunctions(factor,Epsilon,NV)
           List=[]
           if(any(c)):
                          P1=B_VLemma(Poly,k,MaxDegree,NV)
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
			      # lift back cubic form if it exists
                              if((MaxDegree*NV)^3<=len(temp_Poly)):
                                  Cub=numpy.reshape(temp_Poly,(MaxDegree*NV,MaxDegree*NV,MaxDegree*NV))
                                  T3=numpy.zeros((NV,NV,NV),dtype=numpy.int)
                                  for l in range(MaxDegree):
                                      for m in range(MaxDegree):
                                         for k in range(MaxDegree):
                                            T3=T3+Cub[l*NV:(l+1)*NV,m*NV:(m+1)*NV,k*NV:(k+1)*NV]
                                  T3=(T3.flatten())%p
                                  P=numpy.append(T3,P)
                              List.append(P)
                          List.append(PLD)
           #print "UROS ends",datetime.now().time()
           return c,factor,List


def UniformRefine(IPoly,k,Epsilon,NV):
# IPoly is the factor to be regularized.
# Below we have defined the variables, we used :-
# The role of the variable V is defined in the first point of the "Non-trivial tricks used in the program code for efficient computation" section of the paper.
# c is defined in function CheckUniformityOfFunctions.
# List represents the regularized version of the factor defined by the derivatives of the polynomial obtained from the Bogdanov-Viola Lemma(function BVLemma).
        V=[]
        while True:
          c,IPoly,List=UniformRefineOneStep(IPoly,Epsilon,k,NV)
          if(any(c)):
              V=V+List
              continue
          else:
              if(V!=[]):
                  IPoly=IPoly+V
                  V=[]
                  continue
              return IPoly


def Convert_Deg4Poly_to_List_form(Poly):
# MONOM holds all the monomials of Poly having non-zero coefficient value
# COEFF holds the corresponding coefficient values. More clearly:- dotProduct(COEFF,MONOM)=Poly
# Li holds the coefficient of monomials of degree i.
          MONOM=Poly.monomials()
          COEFF=Poly.coefficients()
          L4=[0]*(NV^4)
          L3=[0]*(NV^3)
          L2=[0]*(NV^2)
          L1=[0]*(NV+1)
          for i in range(len(MONOM)):
              t=MONOM[i].variables()
              if any(t):
                  ret=0
                  for k in range(len(t)):
                      index=x.index(t[k])
                      for j in range(MONOM[i].degree(t[k])):
                          ret=ret*NV+index
                  if(MONOM[i].total_degree()==4):
                      L4[ret]=COEFF[i]
                  if(MONOM[i].total_degree()==3):
                      L3[ret]=COEFF[i]
                  if(MONOM[i].total_degree()==2):
                      L2[ret]=COEFF[i]
                  if(MONOM[i].total_degree()==1):
                      L1[ret]=COEFF[i]
              else:
                  L1[NV]=COEFF[i]
          return numpy.array(L4+L3+L2+L1)


def GenerateUniformFactor(IPolylist,k,Epsilon,NV):
# IPolyList holds the input degree 4 polynomial in the list form.
# UFL contains the polynomials of uniform factor in list form.
# UF contains the independent set of UFL
           #print k,Epsilon
           UFL=UniformRefine([IPolylist],k,Epsilon,NV)
           UF=MakeLinearIndependent(UFL)
           return UF

def calculateDegrees(PL,NV):
# PL is the list of polynomials.
# DL is the degree of PL's.
# Assumed min. degree is 1 for each polynomial.
         DL=[1]*len(PL)
         for j in range(len(PL)):
             if(len(PL[j])>=NV^4):
                  if(any(PL[j][NV^4+NV^3:NV^4+NV^3+NV^2])):
                      DL[j]=2
                  if(any(PL[j][NV^4:NV^4+NV^3])):
                      DL[j]=3
                  if(any(PL[j][:NV^4])):
                      DL[j]=4
             else:
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

p=5;NV=30                                                      # p=Field Size,k=no. of polynomials to be generated in BVLemma
SampleSize=p^2+75                                             # "SampleSize"=number of samples to be used for calculating gower's norm
No_Of_Combinations_Check=20				      # No_Of_Combinations_Check specifies how many coefficients lists to generate to ascertain Epsilon uniformity of a factor.
print "No of variables=",NV,"FieldSize",p
print "SampleSize",SampleSize,"No_Of_Combinations_Check",No_Of_Combinations_Check

#Defining Polynomial ring
R=PolynomialRing(GF(p), NV, var_array=['x'])
R.inject_variables()
x=(list)(R.gens())
#xv=vector(x+[1])

Qx=[]
for i in range(NV):
          j=0
          while(j<NV):
               Qx.append(x[i]*x[j])
               j+=1
Qx=Qx+x+[1]
Qxv=vector(Qx)

def Generate_Random_QP():
# Generate quartic polynomial of form pol1*pol2+pol3*pol4 where
# deg(pol1,pol2,pol3,pol4)=2
# IPolylist is equal to the quartic polynomial in list(numpy.array) form
# Let the term 'generator polynomials' specifies the polynomials pol1,pol2,pol3 and pol4.
        pol1=Qxv.dot_product(vector([random.randrange(p) for i in range(NV^2+NV+1)]))
        pol2=Qxv.dot_product(vector([random.randrange(p) for i in range(NV^2+NV+1)]))
        pol3=Qxv.dot_product(vector([random.randrange(p) for i in range(NV^2+NV+1)]))
        pol4=Qxv.dot_product(vector([random.randrange(p) for i in range(NV^2+NV+1)]))
        IPolylist=Convert_Deg4Poly_to_List_form(pol1*pol2+pol3*pol4)
        return pol1, pol2, pol3, pol4, IPolylist

def CheckInSpan(M, Quad):
# This function checks that does "Quad" lies in the span of elements of M or not.
# Here, M will represent the uniform factor in list form and Quad is the polynomial in algebraic form.
        vec = [0]*(NV^2+NV+1)
        for i in range(len(M)):
            Len=len(M[i])
            M[i]=list(M[i])
            if(Len>len(vec)):    #if M contains a polynomial of degree greater than 2
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
# It searches for Epsilon among the values defined by EpsilonE, EpsilonS and diff such that pol2 and pol4 can be written as a linear combination of the polynomials in the uniform factor.
# newk specifies the value of k(=#of derivatives to be evaluated in Bogdanov Viola Lemma)
# EpsilonS and EpsilonE specifies the starting and ending values of Epsilon
# diff denotes the difference between two consecutive values of Epsilon
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
# This function finds the number of bad iterations out of total where bad means that particular generator polynomial cannot be written as a linear combination of the polynomials in the uniform factor.
# newk specifies the value of k(=#of derivatives to be evaluated in Bogdanov Viola Lemma)
# Epsilon specifies the value of epsilon to be used
# total specifies the number of times to execute the regularization procedure.
        k = newk
        bad1=0
        bad2=0
        bad3=0
        bad4=0
        for i in range(total):
            pol1, pol2, pol3, pol4, factor = Generate_Random_QP()
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
            print "Degree list",DL
        print "k",k
        print "bad1",bad1,"bad2",bad2,"bad3",bad3,"bad4",bad4

def PerfOverRandomPoly(newk,EpsilonS,EpsilonE,diff,total):
# This function generates the list of lengths of the uniform factor for Epsilon's, total number of times.
# newk specifies the value of k(=#of derivatives to be evaluated in Bogdanov Viola Lemma)
# EpsilonS and EpsilonE specifies the starting and ending values of Epsilon
# diff denotes the difference between two consecutive values of Epsilon
# Here, total specifies the number of times to execute the regularization procedure.
          for i in range(total):
             print "k",newk
             IPolylist=numpy.array([random.randrange(p) for j in range(NV^4+NV^3+NV^2+NV+1)])
             Epsilon=EpsilonS;k=newk
	     print "Input Polynomial",IPolylist
             UFL_list=[]
             while(Epsilon<=EpsilonE):
                 UFL=GenerateUniformFactor(IPolylist,k,Epsilon,NV)
                 print Epsilon,"length of uniform factor",len(UFL)
                 Epsilon+=diff
                 if(len(UFL)>NV+1):
                     print calculateDegrees(UFL,NV)
                 UFL_list.append(len(UFL))
             print "Length of factors",UFL_list

RunChecksOverk(3,0.8,4)
#PerfOverRandomPoly(2,0.55,0.91,0.05,2)
print "program ends at:",datetime.now().time()











