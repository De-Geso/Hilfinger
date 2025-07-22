# Plot the tests for the mrna-gene system. Here there is no feedback.
# The production rate of x1 is a constant.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fname = "mrna_hill_tests.dat"
fpath = "data/"

data = np.loadtxt(fpath+fname, skiprows=1)

alpha = data[:,0]
beta = data[:,1]
tau1 = data[:,2]
tau2 = data[:,3]
mean1 = data[:,4]
mean2 = data[:,5]
Tmean1 = data[:,6]
Tmean2 = data[:,7]
covar11 = data[:,8]
Tcovar11 = data[:,9]
covar12 = data[:,10]
Tcovar12 = data[:,11]
covar22 = data[:,12]
Tcovar22 = data[:,13]



# Plot mean of x1
fig1, ax1 = plt.subplots()
ax1.set_title(r'Mean mRNA Check. Hill function (k=n=1).')
ax1.set_ylabel(r'$\left< x_1 \right>$')
ax1.set_xlabel(r'$\left< R(x_2) \right> \tau_1$')
ax1.plot(Tmean1, mean1, '.')
lims = [
    np.min([ax1.get_xlim(), ax1.get_ylim()]),  # min of both axes
    np.max([ax1.get_xlim(), ax1.get_ylim()]),  # max of both axes
]
# now plot both limits against eachother
ax1.plot(lims, lims, 'k--', alpha=0.8)
ax1.set_xlim(lims)
ax1.set_ylim(lims)
ax1.set_aspect('equal')
fig1.tight_layout()


# Plot mean of x2
fig2, ax2 = plt.subplots()
ax2.set_title(r'Mean Gene Check. Hill function (k=n=1).')
ax2.set_ylabel(r'$\left< x_2 \right>$')
ax2.set_xlabel(r'$\left< x_1 \right> \beta \tau_2$')
ax2.plot(Tmean2, mean2, '.')
lims = [
    np.min([ax2.get_xlim(), ax2.get_ylim()]),  # min of both axes
    np.max([ax2.get_xlim(), ax2.get_ylim()]),  # max of both axes
]
# now plot both limits against eachother
ax2.plot(lims, lims, 'k--', alpha=0.8)
ax2.set_xlim(lims)
ax2.set_ylim(lims)
ax2.set_aspect('equal')
fig2.tight_layout()


# Plot covariance between x1 and x2
fig12, ax12 = plt.subplots()
ax12.set_title(r'Offdiagonal Covariance Check. Hill function (k=n=1).')
ax12.set_ylabel(r'$\eta_{12(21)}$')
ax12.set_xlabel(r'$\frac{1}{1+\tau_2/\tau_1}$' \
#	r'$\left( \frac{1}{\left<x_1\right>} + \frac{cov(x_1,R(x_2))}{\left<x_1\right>\left<R\right>}\right) + $' \
	r'$\eta_{11}$'
	r'$+ \frac{\tau_2/\tau_1}{1+\tau_2/\tau_1} \frac{cov(x_2,R(x_2))}{\left<x_2\right>\left<R\right>}$')
ax12.plot(Tcovar12 , covar12, '.')
lims = [
    np.min([ax12.get_xlim(), ax12.get_ylim()]),  # min of both axes
    np.max([ax12.get_xlim(), ax12.get_ylim()]),  # max of both axes
]
# now plot both limits against eachother
ax12.plot(lims, lims, 'k--', alpha=0.8)
ax12.set_xlim(lims)
ax12.set_ylim(lims)
ax12.set_aspect('equal')
fig12.tight_layout()


# Plot autocovariance of x1
fig11, ax11 = plt.subplots()
ax11.set_title(r'mRNA Variance Check. Hill function (k=n=1).')
ax11.set_ylabel(r'$\eta_{11}$')
ax11.set_xlabel(r'$\frac{1}{\left< x_1 \right>} + \frac{cov(x_1,R(x_2))}{\left<x_1\right>\left<R\right>}$')
ax11.plot(Tcovar11, covar11, '.')
lims = [
    np.min([ax11.get_xlim(), ax11.get_ylim()]),  # min of both axes
    np.max([ax11.get_xlim(), ax11.get_ylim()]),  # max of both axes
]
# now plot both limits against eachother
ax11.plot(lims, lims, 'k--', alpha=0.8)
ax11.set_xlim(lims)
ax11.set_ylim(lims)
ax11.set_aspect('equal')
fig11.tight_layout()


# Plot autocovariance of x2
fig22, ax22 = plt.subplots()
ax22.set_title(r'Gene Variance Check. Hill function (k=n=1).')
ax22.set_ylabel(r'$\eta_{22}$')
ax22.set_xlabel(r'$\frac{1}{\left< x_2 \right>} + \eta_{12}$')
ax22.plot(Tcovar22, covar22, '.')
lims = [
    np.min([ax22.get_xlim(), ax22.get_ylim()]),  # min of both axes
    np.max([ax22.get_xlim(), ax22.get_ylim()]),  # max of both axes
]
# now plot both limits against eachother
ax22.plot(lims, lims, 'k--', alpha=0.8)
ax22.set_xlim(lims)
ax22.set_ylim(lims)
ax22.set_aspect('equal')
fig22.tight_layout()


plt.show()
