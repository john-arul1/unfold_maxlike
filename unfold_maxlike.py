# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 12:48:00 2019
John Arul, IGCAR, India
Program to unfold spectrum from measured counts
"""

###main driver
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from numpy import zeros
import numpy as np
import math
import sys
import re

## functions

def constr(A,x,b):
    b1=np.zeros(len(b))
    for i in range(len(b)):
        for j in range(len(x)):
            b1[i]+=A[i,j]*x[j]

    return b1-b

def constr_norm(x):

    return sum(x)-1.0

def dist(x,phi):
    res=0.0
    for i in range(len(phi)):
#       if phi[i] < 1.0E-15:
#          phi[i]=1.0E-15
       res+=(x[i]-phi[i])**2
    return res


def Normal(x,R,a):
    c1=R.dot(x)
    res = (c1-a)**2
    res =sum(res)
    return res

def Poisson(x,R,a):

    res=0.0
    for i in range(len(a)):
       ci=0.0
       for j in range(len(x)):
          ci=ci+R[i,j]*x[j]
       if ci > 1.0E-15:
          res += ci-a[i]*math.log(ci)
       else:
           # to be checked
           res+=-15*a[i]

    return res


def L(x,R,a,phi):
    res=0
    #     Likelyhood +      dist to calculated phi+ dist to SVD phi
    res+=Poisson(x,R,a)+dist(x,phi)
    return res


## read files
def read_data(fname):

    try:
        fi=open(fname,'r')    
    except OSError:
        print('cannot open', fname)
        sys.exit(1)
        
    fi.readline() #skip  title
    # two column or one column file
    dim_txt=fi.readline()
    dim_txt=dim_txt.rstrip()
    m,n=dim_txt.split(',')
    m=int(m)
    n=int(n)
    data=fi.readlines()
    fi.close()

    m1=len(data)
    if m != m1:
        print("possible error in  input data")

    B=zeros([m,1])
    i=0
    for line in data:
        line=line.rstrip()
        if n == 2:

            sn,B[i,0]=re.split('\s+',line)  # any white space
            #sn,B[i,0]=line.split("\t",)
        else:
           B[i,0]=line
        i+=1
        if i == m:  # in case files have blank lines
           break

    return B
#------

###main driver

# assumes  minimum 2 input files are given, response and counts
# further files are theoretical flux, guess flux
# there can be standard files of guess flux and theoretical flux
# this is standard flux files have to be implemented

flag_theory=0
flag_guess=0


arguments = sys.argv[:]
narg = len(arguments)
if narg < 2:
    print ("pl. give input file name\n")
    sys.exit(1)

fstr=sys.argv[1]

try:
   file_list=open(fstr,'r')
except OSError:
   print('cannot open', fstr)
   sys.exit(1)

file_list.readline() # ignore first title line
lines=file_list.readlines()

dict_fname={}
for line in lines:
    line=line.rstrip()
    if line=='end':
        break
    print (line)
    key,fname=line.split(',')
    dict_fname.update({key:fname})

if 'response' not in dict_fname:
    print ("Error response matrix file not specified")
    sys.exit(1)

if "counts" not in  dict_fname:
    print ("Measured counts file not specified")
    sys.exit(1)

    
try:
   f1=open(dict_fname["response"],"r")
except OSError:
   print('cannot open', dict_fname["response"])
   sys.exit(1)


f1.readline() #title
dim_txt=f1.readline()
dim_txt=dim_txt.rstrip()
#assumes flux x measurement matrix, col number of energy vector
# E is  first or last coloumn
n,m,Ecol,unit=dim_txt.split(',')
n,m,Ecol=map(int,[n,m,Ecol])

resp_data=f1.readlines()
f1.close()

n1=len(resp_data)

if n1 != n:
    print("response data input error")
    #n=n1   # stll can proceed


Rinp=zeros([n,m])
i=0
for line in resp_data:
    line=line.rstrip()
    values=line.split()
    Rinp[i,:]=values
    i+=1

if Ecol ==1:
    R=Rinp[:,1:m].transpose()  #response matrix
elif Ecol==m:
    R=Rinp[:,:m-1].transpose()  #response matrix

E=np.copy(Rinp[:,0])
E=list(map(float,E))

if 'counts' in dict_fname:
    # read teoretical flux
    CT=read_data('counts')
    B=np.copy(CT[:,0])

if 'theory' in dict_fname:
    # read teoretical flux
    TH=read_data('theory')
    phic=np.copy(TH[:,0])
    flag_theory=1

if 'guess' in dict_fname:
    # read teoretical flux
    GE=read_data('guess')
    phi0=np.copy(GE[:,0])
    phi0n=phi0/sum(phi0)
    flag_guess=1


'''
tikhonov regularization
minimize (Rx-c)+lam *x
'''
###--------------------------------
rp=0.1
X=R.T.dot(R)
X=X+np.identity(n)*rp
XI=np.linalg.inv(X)
x1=XI.dot(R.T.dot(B))
x1s=x1
##---------------------------------


if not flag_theory:

    plt.semilogx(E,x1,label="Regularized")
    plt.xlabel("Energy bin ("+unit+")")
    plt.ylabel(" flux (n/sq.cm/s)")
    plt.grid()
    bn1=R.dot(x1)
    error=(bn1-B)/B*100.0
    print ("Error percent=",error)

else:


 # R . phics = B
    cp=R.dot(phic)
    scale=B.dot(B)/(B.dot(cp))

    phics = phic*scale   #scaled phi to match measured counts
    flnorm=sum(phics)
    phicn=phics/flnorm
    Bn=B/flnorm

# xinit options  i) guess phi0 ii) from regularisation iii) theory
# if guess is given use it, else chack for theoretical value, else default is phi0


    x1n=x1/flnorm


    if flag_guess:
        xinit=phi0n
    elif flag_theory:
        xinit=phicn
    else:
        xinit=phicn#x1n


    ##------------------------------------------------- use SLSQP
    bnds=tuple([(0.0, None) for l in range(len(phic))])
    #cons=({'type': 'eq', 'fun': lambda x:am.dot(x)-m})
    cons=({'type': 'eq', 'fun': lambda x:constr_norm(x)})

    opt={'eps': 1.0E-2, 'maxiter': 1000, 'ftol': 1E-12}
    #opt={'maxiter': 1000, 'ftol': 1e-10}
    #res = minimize(Entro,x0=phicn, args=(phicn),method='SLSQP',bounds=bnds,constraints=cons,options=opt)

    res = minimize(L,x0=xinit,args=(R,Bn,phicn),method='SLSQP',bounds=bnds,constraints=cons,options=opt)
    x2s=res.x*flnorm
    ###-------------------------------------------------


    '''
    Plot and save
    '''

    handle=plt.semilogx(E,x1s,'+',E,x2s,'-',E,phics,'--')
    ax=plt.gca()
    ax.legend(handle,["Regularized","MaxLikelihood","Theory"])
    #E=[1,2,3,4,5]
    plt.xlabel("Energy bin ("+unit+")")
    plt.ylabel(" flux (n/sq.cm/s)")
    #plt.legend()


    plt.grid()
    plt.savefig("unfold_.png")
    #ax.legend(handle1,["ML algorithm"])

    bn1=R.dot(x2s)
    bnreg=R.dot(x1s)
    error=(bn1-B)/B*100.0
    err_reg=(bnreg-B)/B*100.0

    print ("error%=",error,"\n",err_reg)

    fl=open('result_flux.txt','w')

    fl.write("Energy bin ("+unit+")\t Flux in n/cm2/s\n")
    for i in range(len(E)):
        fl.write("{:f},{:.3E}\n".format(E[i],x2s[i]))
    fl.close()

    erf=open('result_error.txt','w')
    erf.write("Percent Error between measured and predicted count")
    erf.write("CPS - Measured,predicted,regularised,error(pre-m)/m,error(reg-m)/m")
    for i in range(len(B)):
        l1=[B[i],bn1[i],bnreg[i],error[i],err_reg[i]]
        for j in l1:
            erf.write("{:.3f}\t".format(j))
        erf.write("\n")
    erf.close()

'''
    hdle=plt.plot(x,amber,'o',x,ambers,'x',x,ambersvd,'+')
    ax=plt.gca()
    ax.legend(hdle,["without entropy","with entropy","svd result"])
    plt.xlabel("data point")
    plt.ylabel("percent error in predcited counts")
    plt.grid()
    plt.savefig("error_ambe2.png")

'''

