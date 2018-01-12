import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from scipy.optimize import curve_fit
import scipy.stats as stats
import seaborn as sns
import xlrd

df = pd.read_excel(r'C:\Burkina\Copy of NFP PhD human based data_Andre Lin.xls')
new_infection_df = pd.read_csv(r'C:\Uganda\new_infection_Garki_with_error.csv')

infectiousness = [float(df['Num infected mosquito'][i])/float(df['Num dissected mosquitoes'][i]) for i in range(len(df['Num infected mosquito']))]
indices_of_density =  df['gf'][df['gf']>0].index.values
infectiousness_subset = [infectiousness[x] for x in indices_of_density]
gam_micro_density = np.log10(df['gf'][df['gf']>0])

# plt.scatter(gam_micro_density,infectiousness_subset)
# plt.show()

def sigmoid(x, x0, k):
    y = 1 / (1 + np.exp(-k * (x - x0)))

    return y

xdata = gam_micro_density
ydata = infectiousness_subset

popt, pcov = curve_fit(sigmoid, xdata, ydata,maxfev=10000)
print popt

x = np.linspace(0, 5, 100)
y = sigmoid(x, *popt)

pylab.plot(xdata, ydata, 'o', label='data')
pylab.plot(x, y, label='fit',linewidth = 2)
pylab.ylim(0, 1)
pylab.legend(loc='best')
pylab.show()

short_mfi_ni = new_infection_df[new_infection_df.interval <=80]
long_mfi_ni = new_infection_df[new_infection_df.interval >=365]

short_mfi_gametocytemia = short_mfi_ni['gametocyte_density']
long_mfi_gametocytemia = long_mfi_ni['gametocyte_density']



short_output = sigmoid(short_mfi_gametocytemia,*popt)
long_output = sigmoid(long_mfi_gametocytemia,*popt)

short_mean = np.mean(short_output)
long_mean = np.mean(long_output)

p_value = stats.ks_2samp(short_output, long_output)

print(short_mean,long_mean,long_mean/short_mean)
print(p_value)

_ = sns.kdeplot(short_output, alpha=0.5, color='#224099', shade=True,bw = 0.075)
_ = sns.kdeplot(long_output, alpha=0.5, color='#a41d2e', shade=True, bw = 0.075)

plt.show()