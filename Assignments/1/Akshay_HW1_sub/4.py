
import numpy as np
import math 
import matplotlib.pyplot as plt
from numpy.linalg import inv

## 1.4 a

Nobs=100
beta = pow(10,-10)#0.1#pow(10,-10)#10 # Vary it 10 # pow(10,-10)

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



fs_synthetic=np.matmul(inv(A),Yobs)


## 1.4 b
# Error
mean=0
sigma = 0.05*np.amax(fs_exact)

Error=np.random.normal(mean,sigma,Nobs+1)
fs_synthetic_noise=np.matmul(inv(A),Yobs-Error)  

plt.figure()
#plt.plot(fs_exact,'g',fs_synthetic,'r--')
plt.title(r'$\beta = 10^{-10} $ ')
without,=plt.plot(fs_synthetic,'b')
withnoise,=plt.plot(fs_synthetic_noise,'g')
plt.legend([without,withnoise],['Without noise','With noise'])
plt.show()

 

# print(Error.shape)

# plt.plot(fs_exact,'g',fs_synthetic,'r--')
# plt.figure()
# plt.plot(Yobs)


# plt.figure()
# plt.plot(fs_exact)
# plt.show()

