import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import LogLocator, AutoLocator
import scipy.stats.mstats as stats
import numpy.ma as ma


#load the banks of simmatched individuals as dictionaries
sim_match_post = np.load('sim_match_MFIrange.npy').item()
sim_match_pre = np.load('sim_match_noMFI.npy').item()

new_infection_df = pd.read_csv(r'C:\Uganda\new_infection_Garki_with_error.csv')

post_list_intervals = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
post_list_ages = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
pre_list_intervals = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
pre_list_ages = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
pre_list_of_asexual_measurements = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
pre_list_of_gametocyte_measurements = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
pre_integrated_infectiousness = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
post_list_of_asexual_measurements = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
post_list_of_gametocyte_measurements = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]
post_integrated_infectiousness = [ [ [],[],[] ],[ [],[],[] ],[ [],[],[] ] ]

interval_bins = [0,80,365,10000]
age_bins = [0,5*365,11*365,10000*365]

counter_try = 0
counter_except = 0

for j in sim_match_post.keys():
    interval = sim_match_post[j]['interval']
    for m in range(len(age_bins)-1):
        if (sim_match_post[j]['age'] > age_bins[m]) & (sim_match_post[j]['age'] <= age_bins[m+1]): #& (sim_match_post[j]['age'] <11*365):
            for k in range(len(interval_bins)-1):
                if (interval > interval_bins[k]) & (interval <= interval_bins[k+1]):
                    try:

                        asexual_to_match = sim_match_post[j]['asexual_to_match']
                        post_list_of_asexual_measurements[m][k].append(asexual_to_match)
                        post_list_intervals[m][k].append(interval)
                        post_list_ages[m][k].append(sim_match_post[j]['age'])
                        index_of_maximum_asexual_density = sim_match_post[j]['pos_asexual_slides'].index(max(sim_match_post[j]['pos_asexual_slides']))

                        gametocytemia_after_maximum_asexual_density = np.mean(sim_match_post[j]['true_gametocytes'][index_of_maximum_asexual_density+7:index_of_maximum_asexual_density+27])

                        post_list_of_gametocyte_measurements[m][k].append(gametocytemia_after_maximum_asexual_density)

                        integration_of_infectiousness_after_maximum = sum(sim_match_post[j]['infectiousness'])
                        post_integrated_infectiousness[m][k].append(float(integration_of_infectiousness_after_maximum))
                        counter_try += 1
                    except:
                        print(str(j) + '_' + str(k) + '_' + str(index_of_maximum_asexual_density))
                        counter_except += 1
                    break

for j in sim_match_pre.keys():
    interval = sim_match_pre[j]['interval']
    for m in range(len(age_bins)-1):
        if (sim_match_pre[j]['age'] > age_bins[m]) & (sim_match_pre[j]['age'] <= age_bins[m+1]): #& (sim_match_pre[j]['age'] <11*365):
            for k in range(len(interval_bins)-1):
                if (interval > interval_bins[k]) & (interval <= interval_bins[k+1]):
                    try:

                        asexual_to_match = sim_match_pre[j]['asexual_to_match']
                        pre_list_of_asexual_measurements[m][k].append(asexual_to_match)
                        pre_list_intervals[m][k].append(interval)
                        pre_list_ages[m][k].append(sim_match_pre[j]['age'])
                        index_of_maximum_asexual_density = sim_match_pre[j]['pos_asexual_slides'].index(max(sim_match_pre[j]['pos_asexual_slides']))

                        gametocytemia_after_maximum_asexual_density = np.mean(sim_match_pre[j]['true_gametocytes'][index_of_maximum_asexual_density+7:index_of_maximum_asexual_density+27])
                        pre_list_of_gametocyte_measurements[m][k].append(gametocytemia_after_maximum_asexual_density)

                        integration_of_infectiousness_after_maximum = sum(sim_match_pre[j]['infectiousness'])
                        pre_integrated_infectiousness[m][k].append(float(integration_of_infectiousness_after_maximum))
                        counter_try += 1
                    except:
                        print(str(j) + '_' + str(k) + '_' + str(index_of_maximum_asexual_density))
                        counter_except += 1
                    break

#for plotting asexuals
a_post = post_list_of_gametocyte_measurements[0][0]
b_post = post_list_of_gametocyte_measurements[0][1]
c_post = post_list_of_gametocyte_measurements[0][2]
d_post = post_list_of_gametocyte_measurements[1][0]
e_post = post_list_of_gametocyte_measurements[1][1]
f_post = post_list_of_gametocyte_measurements[1][2]
g_post = post_list_of_gametocyte_measurements[2][0]
h_post = post_list_of_gametocyte_measurements[2][1]
j_post = post_list_of_gametocyte_measurements[2][2]
j_post[305] = 0

# a_post_logtransform = [-416* np.log(1-float(a_post[i])/200) if a_post[i] < 200 else 2205 for i in range(len(a_post))]
# b_post_logtransform = [-416* np.log(1-float(b_post[i])/200) if b_post[i] < 200 else 2205 for i in range(len(b_post))]
# c_post_logtransform = [-416* np.log(1-float(c_post[i])/200) if c_post[i] < 200 else 2205 for i in range(len(c_post))]
# d_post_logtransform = [-416* np.log(1-float(d_post[i])/200) if d_post[i] < 200 else 2205 for i in range(len(d_post))]
# e_post_logtransform = [-416* np.log(1-float(e_post[i])/200) if e_post[i] < 200 else 2205 for i in range(len(e_post))]
# f_post_logtransform = [-416* np.log(1-float(f_post[i])/200) if f_post[i] < 200 else 2205 for i in range(len(f_post))]
# g_post_logtransform = [-416* np.log(1-float(g_post[i])/200) if g_post[i] < 200 else 2205 for i in range(len(g_post))]
# h_post_logtransform = [-416* np.log(1-float(h_post[i])/200) if h_post[i] < 200 else 2205 for i in range(len(h_post))]
# j_post_logtransform = [-416* np.log(1-float(j_post[i])/200) if j_post[i] < 200 else 2205 for i in range(len(j_post))]

a_pre = pre_list_of_gametocyte_measurements[0][0]
b_pre = pre_list_of_gametocyte_measurements[0][1]
c_pre = pre_list_of_gametocyte_measurements[0][2]
d_pre = pre_list_of_gametocyte_measurements[1][0]
e_pre = pre_list_of_gametocyte_measurements[1][1]
f_pre = pre_list_of_gametocyte_measurements[1][2]
g_pre = pre_list_of_gametocyte_measurements[2][0]
h_pre = pre_list_of_gametocyte_measurements[2][1]
j_pre = pre_list_of_gametocyte_measurements[2][2]
#
# a_pre_logtransform = [-416* np.log(1-float(a_pre[i])/200) if a_pre[i] < 200 else 2205 for i in range(len(a_pre))]
# b_pre_logtransform = [-416* np.log(1-float(b_pre[i])/200) if b_pre[i] < 200 else 2205 for i in range(len(b_pre))]
# c_pre_logtransform = [-416* np.log(1-float(c_pre[i])/200) if c_pre[i] < 200 else 2205 for i in range(len(c_pre))]
# d_pre_logtransform = [-416* np.log(1-float(d_pre[i])/200) if d_pre[i] < 200 else 2205 for i in range(len(d_pre))]
# e_pre_logtransform = [-416* np.log(1-float(e_pre[i])/200) if e_pre[i] < 200 else 2205 for i in range(len(e_pre))]
# f_pre_logtransform = [-416* np.log(1-float(f_pre[i])/200) if f_pre[i] < 200 else 2205 for i in range(len(f_pre))]
# g_pre_logtransform = [-416* np.log(1-float(g_pre[i])/200) if g_pre[i] < 200 else 2205 for i in range(len(g_pre))]
# h_pre_logtransform = [-416* np.log(1-float(h_pre[i])/200) if h_pre[i] < 200 else 2205 for i in range(len(h_pre))]
# j_pre_logtransform = [-416* np.log(1-float(j_pre[i])/200) if j_pre[i] < 200 else 2205 for i in range(len(j_pre))]

fig, axarr = plt.subplots(nrows=2, ncols=3,sharex = True)

axarr[0, 0].hist(a_pre, bins=np.logspace(0, 4, 100), alpha=1, color='#224099', linewidth=0)
axarr[1, 0].hist(c_post, bins=np.logspace(0, 4, 100), alpha=1, color='#a41d2e', linewidth=0)

axarr[0, 1].hist(d_pre, bins=np.logspace(0, 4, 100), alpha=1, color='#224099', linewidth=0)
axarr[1, 1].hist(f_post, bins=np.logspace(0, 4, 100), alpha=1, color='#a41d2e', linewidth=0)

axarr[0, 2].hist(g_pre, bins=np.logspace(0, 4, 100), alpha=1, color='#224099', linewidth=0)


axarr[1, 2].hist(j_post, bins=np.logspace(0, 4, 100), alpha=1, color='#a41d2e', linewidth=0)

plt.gca().set_xscale("symlog")
axarr[0,0].set_yscale("symlog", nonposy='clip')
axarr[1,0].set_yscale("symlog", nonposy='clip')
axarr[0,1].set_yscale("symlog", nonposy='clip')
axarr[1,1].set_yscale("symlog", nonposy='clip')
axarr[0,2].set_yscale("symlog", nonposy='clip')
axarr[1,2].set_yscale("symlog", nonposy='clip')
plt.show()


fig, axarr = plt.subplots(1, 3)
sns.kdeplot(np.log10(a_pre), alpha=0.5, color='#224099', shade=True, ax=axarr[0])
sns.kdeplot(np.log10(c_post), alpha=0.5, color='#a41d2e', shade=True, ax=axarr[0])
sns.kdeplot(np.log10(d_pre), alpha=0.5, color='#224099', shade=True, ax=axarr[1])
sns.kdeplot(np.log10(f_post), alpha=0.5, color='#a41d2e', shade=True, ax=axarr[1])
sns.kdeplot(np.log10(g_pre), alpha=0.5, color='#224099', shade=True, ax=axarr[2])
sns.kdeplot(np.log10(j_post), alpha=0.5, color='#a41d2e', shade=True, ax=axarr[2])

plt.show()

a_inf = post_integrated_infectiousness[0][0]
b_inf = post_integrated_infectiousness[0][1]
c_inf = post_integrated_infectiousness[0][2]
d_inf = post_integrated_infectiousness[1][0]
e_inf = post_integrated_infectiousness[1][1]
f_inf = post_integrated_infectiousness[1][2]
g_inf = post_integrated_infectiousness[2][1]
h_inf = post_integrated_infectiousness[2][1]
j_inf = post_integrated_infectiousness[2][2]


print(np.mean(a_inf))
print(np.mean(c_inf))
print(np.mean(d_inf))
print(np.mean(f_inf))
print(np.mean(g_inf))
print(np.mean(j_inf))

fig, axarr = plt.subplots(1, 3)
sns.kdeplot(np.log10(a_inf), alpha=0.5, color='#224099', shade=True, ax=axarr[0])
sns.kdeplot(np.log10(c_inf), alpha=0.5, color='#a41d2e', shade=True, ax=axarr[0])
sns.kdeplot(np.log10(d_inf), alpha=0.5, color='#224099', shade=True, ax=axarr[1])
sns.kdeplot(np.log10(f_inf), alpha=0.5, color='#a41d2e', shade=True, ax=axarr[1])
sns.kdeplot(np.log10(g_inf), alpha=0.5, color='#224099', shade=True, ax=axarr[2])
sns.kdeplot(np.log10(j_inf), alpha=0.5, color='#a41d2e', shade=True, ax=axarr[2])

plt.gca().set_xscale("linear")
plt.show()
print('a')
#
#
# c_gam = [x for x in pre_list_of_gametocyte_measurements[0] if str(x) != 'nan']#[-416* np.log(1-float(post_list_of_gametocyte_measurements[0][i])/200) if post_list_of_gametocyte_measurements[0][i] < 200 else 2205 for i in range(len(post_list_of_gametocyte_measurements[0]))]
# d_gam = [x for x in post_list_of_gametocyte_measurements[2] if str(x) != 'nan']#[-416* np.log(1-float(post_list_of_gametocyte_measurements[2][i])/200) if post_list_of_gametocyte_measurements[2][i] < 200 else 2205 for i in range(len(post_list_of_gametocyte_measurements[2]))]
#
# e_inf = post_integrated_infectiousness[0]
# f_inf = post_integrated_infectiousness[2]
#
# e_mean = stats.gmean([x for x in k if x >0])
# f_mean = stats.gmean([x for x in f if x >0])
# print(e_mean)
# print(f_mean)

#
# fig, axarr = plt.subplots(nrows=2, ncols=2,sharex = True)
# #
# axarr[0,0].hist(g, bins=np.logspace(0, 4, 100), alpha=1, color='#c6c6c6', linewidth=0) #np.logspace(0, 6, 100)
# axarr[1,0].hist(b, bins=np.logspace(0, 4, 100), alpha=1, color='#7c7c7c', linewidth=0)
# axarr[0,1].hist(i, bins=np.logspace(0, 4, 100), alpha=1, color='#c6c6c6', linewidth=0)
# axarr[1,1].hist(d, bins=np.logspace(0, 4, 100), alpha=1, color='#7c7c7c', linewidth=0)
# # axarr[0, 1].hist(c, bins=np.logspace(0, 6, 100), alpha=1, color='#c6c6c6', linewidth=0)
# # axarr[1, 1].hist(d, bins=np.logspace(0, 6, 100), alpha=1, color='#7c7c7c', linewidth=0)
# # axarr[0, 2].hist(e, bins=np.logspace(0, 6, 100), alpha=1, color='#c6c6c6', linewidth=0)
# # axarr[1, 2].hist(f, bins=np.logspace(0, 6, 100), alpha=1, color='#7c7c7c', linewidth=0)
# #
# plt.gca().set_xscale("symlog")
# axarr[0,0].set_yscale("symlog", nonposy='clip')
# axarr[1,0].set_yscale("symlog", nonposy='clip')
# axarr[0,1].set_yscale("symlog", nonposy='clip')
# axarr[1,1].set_yscale("symlog", nonposy='clip')
# plt.show()
#
# fig_,axarr_ = plt.subplots(nrows= 2, ncols= 1, sharex=True)
#
# axarr_[0].hist(k, bins=100, alpha=1, color='#c6c6c6', linewidth=0)
# axarr_[1].hist(f, bins=100, alpha=1, color='#7c7c7c', linewidth=0)
# axarr_[0].set_yscale("log", nonposy='clip')
# axarr_[1].set_yscale("log", nonposy='clip')
# plt.show()

# axarr[0, 1].set_yscale("log", nonposy='clip')
# axarr[1, 1].set_yscale("log", nonposy='clip')
# axarr[0, 2].set_yscale("log", nonposy='clip')
# axarr[1, 2].set_yscale("log", nonposy='clip')
#
#
#
#
#
#     try:
#         max_asexual_density = max(sim_match_post[j]['true_asexual_parasites'])
#         max_gametocyte_density = max(sim_match_post[j]['true_gametocytes'])
#         infectiousness = sum(sim_match_post[j]['infectiousness'])
#
#         maxes_of_indiviudal_asexuals.append(max_asexual_density)
#         maxes_of_indiviudal_gametocytes.append(max_gametocyte_density)
#         infectiousness_list.append(infectiousness)
#
#     except:
#         print(str(j)+'_'+str(i))
#
#
# # -416*np.log(1-float(df_Garki.loc[row, 'pfa'])/(df_Garki.loc[row,'exam']))
# pre_maxes_of_indiviudal_asexuals =[]
# pre_maxes_of_indiviudal_gametocytes =[]
# pre_infectiousness_list = []
#
# for j in sim_match_pre.keys():
#     try:
#         max_asexual_density = max(sim_match_pre[j]['true_asexual_parasites'])
#         max_gametocyte_density = max(sim_match_pre[j]['true_gametocytes'])
#         infectiousness = sum(sim_match_pre[j]['infectiousness'])
#
#         pre_maxes_of_indiviudal_asexuals.append(max_asexual_density)
#         pre_maxes_of_indiviudal_gametocytes.append(max_gametocyte_density)
#         pre_infectiousness_list.append(infectiousness)
#     except:
#         print(str(j)+'_'+str(i))
#
#
# ind = np.arange(2)
#
# fig, axarr = plt.subplots(1,5)
# axarr[0].boxplot([np.log10(pre_maxes_of_indiviudal_asexuals)])
#
# axarr[1].boxplot([np.log10(maxes_of_indiviudal_asexuals)])
#
# axarr[2].boxplot(np.log10([i for i in pre_maxes_of_indiviudal_gametocytes if i>0]))
# axarr[3].boxplot(np.log10([i for i in maxes_of_indiviudal_gametocytes if i>0]))
#
# # axarr[0].bar(ind, [stats.gmean(pre_maxes_of_indiviudal_asexuals),stats.gmean(maxes_of_indiviudal_asexuals)],width = 0.5)
# # axarr[1].bar(ind, [stats.gmean([i for i in pre_maxes_of_indiviudal_gametocytes if i>0]),stats.gmean([i for i in maxes_of_indiviudal_gametocytes if i>0])],width = 0.5)
# axarr[4].bar(ind, [np.mean(pre_infectiousness_list),np.mean(infectiousness_list)],width = 0.5)
#
# # axarr[0].hist(maxes_of_indiviudal_asexuals,alpha = 0.5, color = 'r')
# # axarr[0].hist(pre_maxes_of_indiviudal_asexuals,alpha = 0.5, color = 'b')
# # axarr[1].hist(maxes_of_indiviudal_gametocytes,alpha = 0.5, color = 'r')
# # axarr[1].hist(pre_maxes_of_indiviudal_gametocytes,alpha = 0.5, color = 'b')
# # axarr[2].hist(infectiousness_list,alpha = 0.5, color = 'r')
# # axarr[2].hist(pre_infectiousness_list,alpha = 0.5, color = 'b')
#
# plt.show()
