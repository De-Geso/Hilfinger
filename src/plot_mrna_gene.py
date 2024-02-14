# Plot the tests for the mrna-gene system. Here there is no feedback.
# The production rate of x1 is a constant.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fname = "mrna_data.csv"
fpath = "data/"

data = pd.read_csv(fpath+fname, header="infer")

for col in data.columns:
    print(col)

    
# Plot mean of x1
fig1, ax1 = plt.subplots()
ax1.set_title(r'Mean mRNA Check. Constant production rate $\alpha.$')
ax1.set_ylabel(r'$\left< x_1 \right>$')
ax1.set_xlabel(r'$\alpha / d_1$')
ax1.plot(data.loc[:,'alpha']/data.loc[:,'decay1'], data.loc[:,'mean1'], '.')
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
ax2.set_title(r'Mean Gene Check.')
ax2.set_ylabel(r'$\left< x_2 \right>$')
ax2.set_xlabel(r'$\left< x_1 \right> \beta / d_2$')
ax2.plot(data.loc[:,'mean1']*data.loc[:,'beta']/data.loc[:,'decay2'], data.loc[:,'mean2'], '.')
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
ax12.set_title(r'Offdiagonal Covariance Check.')
ax12.set_ylabel(r'$\eta_{12(21)}$')
ax12.set_xlabel(r'$\frac{1}{\left< x_1 \right>}\frac{d_2}{d_1+d_2}$')
ax12.plot(1/data.loc[:,'mean1']*data.loc[:,'decay2']/(data.loc[:,'decay2']+data.loc[:,'decay1']), data.loc[:,'covar12'], '.')
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



# Plot autocovariance of x2
fig22, ax22 = plt.subplots()
ax22.set_title(r'Gene Variance Check.')
ax22.set_ylabel(r'$\eta_{22}$')
ax22.set_xlabel(r'$\frac{1}{\left< x_2 \right>} + \eta_{12(21)}$')
ax22.plot(1/data.loc[:,'mean2']+data.loc[:,'covar12'], data.loc[:,'covar22'], '.')
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
