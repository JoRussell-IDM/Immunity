import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import json
from pprint import pprint
from malaria.createSimDirectoryMap import createSimDirectoryMap
from operator import itemgetter
import itertools
import scipy.stats as stats
import math

expname = 'Challenge_bite_noMFI_infectiousness'
exp_dir = createSimDirectoryMap(expname)

exp_dir.sort_values('start_day', inplace=True)
exp_dir.reset_index(drop= True, inplace=True)
sim_bank = {}

generate_sim_bank = True

if generate_sim_bank == True:
    for i in range(exp_dir.shape[0]):
        interval = exp_dir.start_day[i]
        MPR_path = os.path.join(exp_dir.outpath[i],'output')
        files = [x for x in os.listdir(MPR_path) if 'Survey' in x]
        sim_bank[interval] = {}
        try:
            with open(os.path.join(MPR_path,files[0])) as MPR:

                data = json.load(MPR)

                #print(data['patient_array'][0]['id'])
                number_patients = len(data['patient_array'])
                sim_bank[interval]['id'] = [data['patient_array'][x]['id'] for x in range(number_patients)]
                sim_bank[interval]['interval'] = interval
                sim_bank[interval]['age'] = [data['patient_array'][x]['initial_age'] for x in range(number_patients)]
                sim_bank[interval]['infectiousness'] = [data['patient_array'][x]['infectiousness'] for x in range(number_patients)]
                sim_bank[interval]['pos_asexual_slides'] = [data['patient_array'][x]['pos_asexual_fields'] for x in range(number_patients)]
                sim_bank[interval]['true_asexuals'] = [data['patient_array'][x]['true_asexual_parasites'] for x in
                                                            range(number_patients)]

                sim_bank[interval]['pos_gametocyte_slides'] = [data['patient_array'][x]['pos_gametocyte_fields'] for x in range(number_patients)]
                sim_bank[interval]['true_gametocytes'] = [data['patient_array'][x]['true_gametocytes'] for x in
                                                       range(number_patients)]
                sim_bank[interval]['temps'] = [data['patient_array'][x]['temps'] for x in range(number_patients)]
        except:
            print('interval: '+str(interval)+" Malaria Patient Report could not be found")

fever_proportions = []

age_bins = [0, 1, 5, 11, 2000]
for age in range(len(age_bins) - 1):
    for j in range(0):#20, 350, 10):
        try:
            sim_bank_subset_indices = [i for i, x in enumerate(sim_bank[j]['age']) if
                                       ((x > 365 * age_bins[age]) & (x < 365 * age_bins[age + 1]))]

            population = len(sim_bank_subset_indices)
            fever_sums = []
            for k in sim_bank_subset_indices:

                fever_duration = np.nansum(
                [np.nan if x == (-1) else 1 for x in sim_bank[j]['temps'][k]])
                fever_sums.append(fever_duration)
            proportion = (float(np.count_nonzero(fever_sums)) / float(population))
            fever_proportions.append([j, proportion])
        except:
            print(j)
            pass
fig, axarr = plt.subplots(nrows=1, ncols=4, sharex=True, sharey = True)

#make a dataframe of Garki reference first infections after a MFI
df = pd.read_csv('C:\Uganda\DF_Garki_08_07.csv', index_col=0, parse_dates=['datecoll'])
df.asexual_density[df.asexual_density == 'maxed'] = float(2205)
df.asexual_density = pd.to_numeric(df.asexual_density, errors='coerce')
df.gametocyte_density = pd.to_numeric(df.gametocyte_density, errors='coerce')

#
def new_infection_events_Garki(ind_df):
    """
    Identify based on microscopy and LAMP the time interval since last malaria-positive measurement,
    and return re-infection events with information on id, date, age, malaria-free-interval, etc.
    """
    first_positive = []
    lapsed_positive = []
    first_malaria_free_date = None  # first time point is left-censored on malaria-free interval
    for ind, measurements in ind_df.iterrows():
        # initialize positive flag as False

        positive = False
        if measurements.asexual_density > 0:
            positive = True
        elif measurements.gametocyte_density > 0:
            positive = True
        if first_malaria_free_date is None:
            if positive:
                first_malaria_free_date = None
            else:
                first_malaria_free_date = measurements.datecoll
        else:
            interval = (measurements.datecoll - first_malaria_free_date) / np.timedelta64(1, 'D')
            if positive:
                time_since_last = (measurements.datecoll - ind_df.loc[ind-1,:].datecoll) / np.timedelta64(1, 'D')
                first_malaria_free_date = None

                if time_since_last < 100:

                    first_positive.append({'date': measurements.datecoll, 'area': measurements.area, 'vname':measurements.vname,
                                        'age': measurements.age, 'age_at_enroll':measurements.age_at_enroll,'IRS_status':measurements.IRS_status,'id': measurements.id, 'fever': measurements.fever,
                                        'asexual_density': measurements.pfa, 'gametocyte_density': measurements.pfg, 'examined': measurements.exam, 'interval': interval, 'time_since_last': time_since_last,'frequency': measurements.infection_frequency})
                else:
                    lapsed_positive.append({'date': measurements.datecoll, 'area': measurements.area,'vname':measurements.vname,
                                        'age': measurements.age, 'age_at_enroll':measurements.age_at_enroll,'IRS_status':measurements.IRS_status,'id': measurements.id, 'fever': measurements.fever,
                                        'asexual_density': measurements.pfa, 'gametocyte_density': measurements.pfg, 'examined': measurements.exam, 'interval': interval,'time_since_last': time_since_last,'frequency': measurements.infection_frequency})
                if measurements.IRS_status == 'post_IRS':
                    break
            else:
                pass

    return first_positive,lapsed_positive

first_positives = []
lapsed_positives = []

recreate_new_infections_df = False

if recreate_new_infections_df == True:
    for ind_id, ind_df in df.groupby('id'):
        fp, lp = new_infection_events_Garki(ind_df)

        first_positives += fp
        lapsed_positives += lp

    new_infection_df = pd.DataFrame.from_dict(first_positives)  # .drop_duplicates(subset='id', keep='last').set_index('id')
    lapsed_new_infection_df = pd.DataFrame.from_dict(lapsed_positives)

    new_infection_df = new_infection_df[new_infection_df.area == 'A1']
    new_infection_df.fever = new_infection_df.fever == 1
    new_infection_df.reset_index(inplace=True)
    # new_infection_df = new_infection_df[new_infection_df.IRS_status == 'post_IRS']


    new_infection_df['fraction_asexual_fields'] = new_infection_df['asexual_density'].astype(float)/new_infection_df['examined'].astype(float)
    new_infection_df['error'] = [float(new_infection_df['examined'][x])*math.sqrt((new_infection_df['fraction_asexual_fields'][x]*(1-new_infection_df['fraction_asexual_fields'][x]))/float(new_infection_df['examined'][x])) for x in new_infection_df.index]

else:
    new_infection_df = pd.read_csv(r'C:\Uganda\new_infection_Garki_with_error.csv')

sim_match_set = []
def round_down(num, divisor):
    return num - (num%divisor)

#new_infection_df_test = new_infection_df[(new_infection_df.interval<110) &(new_infection_df.interval>60)]
# new_infection_df = new_infection_df[new_infection_df.IRS_status == 'post_IRS']

unmatched_infections = []
sim_match = {}
for i,infection in new_infection_df.iterrows():

    #compare in order age,interval,asexual density, gam +/-
    #compare to ref interval rounded down to multiples of 10s
    interval_to_match = 0#round_down(int(infection.interval),10)
    age_to_match = int(365*infection.age)
    gametocytes_to_match = bool(infection.gametocyte_density)
    asexual_to_match = int(infection.asexual_density)
    error_on_density = int(infection.error)

    #first produce a subset of list for all the individuals within 1 years of ref ind in this interval bin
    sim_bank_indices_age_matched = [k for k, x in enumerate(sim_bank[interval_to_match]['age']) if (x < age_to_match+190) & (x > age_to_match-190)]
    match_list = []
    interval_bank_list = []
    for person in sim_bank_indices_age_matched:
        try:
            start_compare_range = min([l for l, x in enumerate(sim_bank[interval_to_match]['pos_asexual_slides'][person]) if x > 0])
            dt = 90
            #among individuals of like ages, who has a density measurement equal to our ref ind (on any day)
            if any(range(asexual_to_match - error_on_density, asexual_to_match + error_on_density)) in sim_bank[interval_to_match]['pos_asexual_slides'][person]:
                match_list.append(person)
                interval_bank_list.append(interval_to_match)
            #what is that day of measurement?
            #index_day = sim_bank[interval_to_match]['pos_asexual_slides'][person].index(asexual_to_match)
            #now if that day fo measrement has a match of +/- gametocytes
            #if bool(sim_bank[interval_to_match]['pos_gametocyte_slides'][person][index_day]) == gametocytes_to_match:
                #add this matched sim individual to the matched dictionary with all associated fields
        except:
            print('positive slides empty for '+str(person)+'_'+str(i))
    #if there are any matches
    if match_list:
        #select one at random
        random_selection = np.random.choice(match_list)
        #and append the relevnat fields to the dictionary matched individuals sim_match
        sim_match[int(infection.id)] = {'age':sim_bank[interval_to_match]['age'][random_selection],
                                           'id':sim_bank[interval_to_match]['id'][random_selection],
                                           'interval': interval_to_match,
                                           'asexual_to_match': asexual_to_match,
                                           'infectiousness': sim_bank[interval_to_match]['infectiousness'][random_selection],
                                           'pos_gametocyte_slides': sim_bank[interval_to_match]['pos_gametocyte_slides'][random_selection],
                                           'true_gametocytes': sim_bank[interval_to_match]['true_gametocytes'][
                                            random_selection],
                                           'pos_asexual_slides': sim_bank[interval_to_match]['pos_asexual_slides'][random_selection],
                                           'true_asexual_parasites': sim_bank[interval_to_match]['true_asexuals'][
                                            random_selection],
                                           'temps':sim_bank[interval_to_match]['temps'][random_selection]
                                           }


    if not match_list:
        unmatched_infections.append(infection)


np.save('sim_match_noMFI.npy', sim_match)

print(len(unmatched_infections))
print(sim_match.keys())

