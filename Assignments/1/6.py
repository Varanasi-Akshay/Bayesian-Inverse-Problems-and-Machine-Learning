import numpy as np
from scipy.sparse.linalg import cg
#import tensorflow as tf
import time


def conjugate_grad(A, b, sigma ):
    """
    Description
    -----------
    Solve a linear equation Ax = b with conjugate gradient method.

    Parameters
    ----------
    A: 2d numpy.array of positive semi-definite (symmetric) matrix
    b: 1d numpy.array
    x: 1d numpy.array of initial point

    Returns
    -------
    1d numpy.array x such that Ax = b
    """
    n = len(b)

    # if not x:
    x = np.ones(n)
    r = np.dot(A, x) - b
    p = - r
    r_k_norm = np.dot(r, r)
    for i in range(2*n):
        Ap = np.dot(A, p)
        alpha = r_k_norm / np.dot(p, Ap)
        x += alpha * p
        r += alpha * Ap
        r_kplus1_norm = np.dot(r, r)
        beta = r_kplus1_norm / r_k_norm
        r_k_norm = r_kplus1_norm
        if pow(r_kplus1_norm,2) < (n+1)*pow(sigma,2):
            print ('Itr:', i)
            break
        p = beta * p - r
    return x

# # x_ph = tf.placeholder('float32', [None, None])
# # r = tf.matmul(A, x_ph) - b

# if __name__ == '__main__':
#     n = 1000
#     P = np.random.normal(size=[n, n])
#     A = np.dot(P.T, P)
#     b = np.ones(n)

#     t1 = time.time()
#     print ('start')
#     x = conjugate_grad(A, b)
#     t2 = time.time()
#     print (t2 - t1)
#     x2 = np.linalg.solve(A, b)
#     t3 = time.time()
#     print (t3 - t2)
#     x3 = cg(A, b)
#     t4 = time.time()
#     print (t4 - t3)

#     # print np.dot(A, x) - b


import numpy as np
import math 
import matplotlib.pyplot as plt
from numpy.linalg import inv
from numpy.linalg import norm

## 1.6 a

Nobs=100
beta = 0.2
Const = 1/math.sqrt(2*math.pi*pow(beta,2))

Yobs=np.zeros(Nobs+1)
m=np.zeros(Nobs+1)
sj=np.zeros(Nobs+1)
#si=np.zeros(Nobs+1)
A=np.zeros((Nobs+1,Nobs+1))
fs_exact=np.zeros(Nobs+1)
fs_synthetic=np.zeros(Nobs+1)

for j in range(Nobs+1):
    sj[j]=j/Nobs

    # True solution
    fs_exact[j]=math.sin(2*math.pi*sj[j])

    # Doing the integration using rectangle rule
    for x in range(Nobs+1):
        Yobs[j]=Yobs[j]+Const*math.exp(-1*pow((x/Nobs-sj[j]),2)/(2*pow(beta,2)))*math.sin(2*math.pi*x/Nobs)/Nobs
    
    # Finding the kernel matrix
    for i in range(Nobs+1):
        si=i/Nobs
        A[i,j]= Const*math.exp(-1/(2*pow(beta,2))*pow((si-sj[j]),2))


k=1

gamma=np.zeros((Nobs+1,Nobs+1))
for i in range(Nobs):
    gamma[i,i]=1
    gamma[i,i+1]=-1
gamma[Nobs,Nobs]=1    
fs_synthetic=inv(np.transpose(A)@A+k*np.transpose(gamma)@gamma)@(np.transpose(A)@Yobs)

# Error
mean=0
sigma = 0.05*np.amax(fs_exact)

Error=np.random.normal(mean,sigma,Nobs+1)
fs_synthetic_noise=inv(np.transpose(A)@A+k*np.transpose(gamma)@gamma)@(np.transpose(A)@(Yobs-Error))

plt.figure()
#plt.plot(fs_exact,'g',fs_synthetic,'r--')
plt.title(r'$\kappa =1 $ ')
without,=plt.plot(fs_synthetic,'b')
withnoise,=plt.plot(fs_synthetic_noise,'g')
gt,=plt.plot(fs_exact,'y')
plt.legend([without,withnoise,gt],['Without noise','With noise','Actual'])
#plt.show()


## 1.6 b Should plot the misfit with different k values

Min=0.01 #min value for k 
Max=5 #max value for k
K=20 # No.of values between min and maz

misfit=np.zeros(K)
k=np.zeros(K)
for i in range(K):
    k[i]=(i+1)*(Max-Min)/K+Min
    m=inv(np.transpose(A)@A+k[i]*np.transpose(gamma)@gamma)@(np.transpose(A)@(Yobs-Error))
    misfit[i]=pow((norm((A@m-fs_exact),2)),2)
plt.figure()
plt.plot(k,misfit)
plt.ylabel('Misfit')
plt.xlabel('k')
# x=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
plt.title('Misfit vs k')



## 1.6 c using Conjugate Gradient 
kappa=0.01
m_cg=conjugate_grad(np.transpose(A)@A,np.transpose(A)@(Yobs-Error),sigma)
m_tik=inv(np.transpose(A)@A+kappa*np.transpose(gamma)@gamma)@(np.transpose(A)@(Yobs-Error))
plt.figure()
plt.title('Comparision of Tikhonov and Conjugate Gradient')
CG,=plt.plot(m_cg,'b')
Tik,=plt.plot(m_tik,'g')
gt,=plt.plot(fs_exact,'y')
plt.legend([CG,Tik,gt],['Conjugate Gradient','Tikhonov k=0.01','Actual'])

plt.figure()
plt.title('Comparision of Tikhonov and Conjugate Gradient')
CG,=plt.plot(m_cg,'b')
Tik,=plt.plot(m_tik,'g')
#gt,=plt.plot(fs_exact,'y')
plt.legend([CG,Tik],['Conjugate Gradient','Tikhonov k=0.01'])


plt.show()

# -*- coding: utf-8 -*-

