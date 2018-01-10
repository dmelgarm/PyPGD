from numpy.random import multivariate_normal
from numpy import array
import ndtest
from matplotlib import pyplot as plt


# Draw two samples from the SAME distribution
mean=[4,7]
Nsamples=10000
cov=array([[4,1],[2,1]])
s1=multivariate_normal(mean,cov,Nsamples)
s2=multivariate_normal(mean,cov,Nsamples)

#Run the KS test
result=ndtest.ks2d2s(s1[:,0], s1[:,1], s2[:,0], s2[:,1])

#plot the two data sets
bins=50
xl=[0,10]
yl=[4,10]
vmin=0
vmax=10

plt.figure(figsize=(10,4))
plt.subplot(121)
plt.hist2d(s1[:,0],s1[:,1],bins=50,cmap=plt.cm.magma_r,vmin=vmin,vmax=vmax)
plt.colorbar()
plt.ylim(yl)
plt.xlim(xl)

plt.subplot(122)
plt.hist2d(s2[:,0],s2[:,1],bins=50,cmap=plt.cm.magma_r,vmin=vmin,vmax=vmax)
plt.colorbar()
plt.ylim(yl)
plt.xlim(xl)

plt.show()