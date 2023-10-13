import pandas as pd
import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from datetime import datetime as dt
from calc_endsplit import wtd_mean


sites = ['ARIK', 'BIGC', 'BLDE', 'BLUE', 'BLWA', 'COMO', 'FLNT', 'HOPB', 'KING', 'LECO', 'LEWI',
            'MART', 'MAYF', 'MCDI', 'POSE', 'PRIN', 'REDB', 'SYCA', 'TECR', 'TOMB', 'WALK', 'WLOU']

# corresponding sites with precipitation data in the same order
precip = ['ARIK', 'SJER', 'YELL', 'BLUE', 'DELA', 'NIWO', 'JERC', 'HARV', 'KONZ', 'GRSM', 'BLAN',
            'WREF', 'TALL', 'KONZ', 'SCBI', 'PRIN', 'REDB', 'SYCA', 'SJER', 'LENO', 'ORNL', 'NIWO']

# Tells if primary is available for precip sites listed above. If primary precipitation is not available, use secondary
primary = [True, True, True, True, False, True, False, True, True, False, False, True, False, True, True, True, True,
           False, True, False, True, True]

def time(dt_str):
    t = dt.strptime(dt_str, '%Y-%m-%dT%H:%MZ')
    return t

def time_d(dt_str):
    t = dt.strptime(dt_str, '%Y-%m-%dT%H:%M:%SZ')
    return t

# Combine all precipitation isotope data and save as pickle file

# directory where I've saved all precipitation isotope data downloaded from NEON
maindir = r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\NEON\NEON_isotope-ppt\NEON_isotope-ppt'
"""

summary = {}
cols = ['start_date', 'end_date', 'd2h', 'd2h_sd', 'wt']
for i in range(len(sites)):
    data_summary = pd.DataFrame(columns=['start_date', 'end_date', 'd2h', 'd2h_sd', 'wt'])
    for dirname in os.listdir(maindir):
        if precip[i] in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if 'isoPerSample' in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname) 
                    for x in range(len(data['collectDate'])):
                        # slice dates to remove time information at the end
                        new_row = {'start_date': [data['setDate'][x]], 'end_date': [data['collectDate'][x]], 'd2h': [data['d2HWater'][x]], 'd2h_sd': [data['d2HsdWater'][x]], 'wt': [np.nan]}
                        data_summary = pd.concat([data_summary, pd.DataFrame(new_row)], keys=cols, ignore_index=True)
    summary[precip[i]] = data_summary

# df with cols listed above, but empty 'wt' column
with open('NEON_precip.pickle', 'wb') as f:
    pickle.dump(summary, f)

"""

# Copy primary or secondary precipitation data with date information and save to pickle file
maindir = r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\NEON\NEON_precipitation\NEON_precipitation'

"""
precip_30m = {}
for i in range(len(sites)):
    data_summary = pd.DataFrame(columns=['date', 'ppt_mm'])
    for dirname in os.listdir(maindir):
        if precip[i] in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if primary[i]:
                    f = 'PRIPRE_30min'
                    c = 'priPrecipBulk'
                    qf = 'priPrecipFinalQF'
                else:
                    f = 'SECPRE_30min'
                    c = 'secPrecipBulk'
                    qf = 'secPrecipSciRvwQF'
                if f in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname)
                    data = data[data[qf] != 1].reset_index()   # remove records if quality flag is 1
                    data = data.rename(columns={'endDateTime': 'date', c: 'ppt_mm'})
                    data_summary = pd.concat([data_summary, data.loc[:, ('date', 'ppt_mm')]], keys=['date', 'ppt_mm'], ignore_index=True)
    precip_30m[precip[i]] = data_summary


with open('NEON_precip_30min.pickle', 'wb') as f:
    pickle.dump(precip_30m, f)

"""

with open('NEON_precip.pickle', 'rb') as f:
    summary = pickle.load(f)

with open('NEON_precip_30min.pickle', 'rb') as f:
    precip_30m = pickle.load(f)

#calculate mean 30 min precipitation amount for each month from daily
        #create a new column for month, then group by month and calculate mean

month_30m_means = {}
for i in precip:
    precip_rates = precip_30m[i]
    month = []
    for j in precip_rates['date']:
        month.append(j.split('-')[1])
    precip_rates['month'] = month
    month_means = precip_rates.drop('date', axis=1).groupby(['month']).mean()
    month_30m_means[i] = month_means


#calclate mean 30 min precipitation amount for each sampling period
        #for any that are zero, make 'wt'

for i in precip:
    data = summary[i]
    daily = precip_30m[i]
    wt = []
    month = []
    for j in range(len(data['start_date'])):
        month.append(time(data['start_date'][j]).strftime('%m'))
    data['month'] = month
    for j in range(len(data['start_date'])):
        n = len(data[data['month'] == data['month'][j]])
        month = time(data['start_date'][j]).strftime('%m')
        month_ppt = month_30m_means[i].loc[month, 'ppt_mm']
        ppt = []
        for d in range(len(daily['date'])):
            if time_d(daily['date'][d]) >= time(data['start_date'][j]) and time_d(daily['date'][d]) <= time(data['end_date'][j]):
                ppt.append(daily['ppt_mm'][d])
        if np.isnan(sum(ppt)) or sum(ppt) == 0:
            mean_ppt = month_ppt    # replace nans and zeros with the average monthly value because zero values are
                                    # impossible and there are a lot of gaps in the precipitation data
        else:
            mean_ppt = sum(ppt)/len(ppt)
        weight = (1 + (mean_ppt - month_ppt) / month_ppt) * (month_ppt / n)
        wt.append(weight)  #W = ( 1 + (Pi - Pm) / Pm) (Pm / N)
    data['wt'] = wt
    summary[i] = data
    print(i)


#with open('NEON_precip_wt.pickle', 'wb') as f:
#    pickle.dump(summary, f)


# calculate seasonal (April – September v. October – March) isotope values as weighted means
# apply elevation correction
# start date defines which month the sample is for
# if the sum of precipitation over the sample time period is zero (precip must fall to get a sample) or nan,
        # then the average 30-min interval value associate with that month is used
# weights are calculated from the equation W = ( 1 + (Qi - Qm) / Qm) (Qm / N)

#with open('NEON_precip_wt.pickle', 'rb') as f:
#    summary = pickle.load(f)

elevation = pd.read_csv(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\elevations.csv')

pw = []
pw_se = []
ps = []
ps_se = []
for i in precip:
    data = summary[i]
    index = list(elevation['precip_site']).index(i)
    elev_diff = elevation['runoff_avg_elev_m'][index] - elevation['precip_elev_m'][index]
    d2h_correction = elev_diff * -2.24 / 100
    iso_val_sw = [[] for n in range(2)]
    wt_sw = [[] for n in range(2)]
    for j in range(len(data['start_date'])):
        month = time(data['start_date'][j]).strftime('%m')
        if int(month) > 3 and int(month) < 10:
            season = 0
        elif int(month) < 4 or int(month) > 9:
            season = 1
        iso_val_sw[season].append(data['d2h'][j])
        wt_sw[season].append(data['wt'][j])
    summer_pdel = wtd_mean(iso_val_sw[0], wt_sw[0])
    winter_pdel = wtd_mean(iso_val_sw[1], wt_sw[1])
    pw.append(winter_pdel[0] + d2h_correction)
    pw_se.append(winter_pdel[1])
    ps.append(summer_pdel[0] + d2h_correction)
    ps_se.append(summer_pdel[1])
    print(i, winter_pdel[0] + d2h_correction, winter_pdel[1], summer_pdel[0] + d2h_correction, summer_pdel[1])

df = pd.DataFrame(data={'site': precip, 'pw': pw, 'pw_se': pw_se, 'ps': ps, 'ps_se': ps_se}, columns=['site', 'pw', 'pw_se', 'ps', 'ps_se'])
#df.to_csv(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\p_d2h_neon.csv')
quit()


sizes = {'ARIK': 2, 'SJER': 20, 'YELL': 4, 'BLUE': 2, 'DELA': 20, 'NIWO': 3, 'JERC': 30, 'HARV': 2, 'KONZ': 3,
         'GRSM': 30, 'BLAN': 30, 'WREF': 1, 'TALL': 2, 'SCBI': 4, 'PRIN': 3, 'REDB': 30, 'SYCA': 30, 'SJER': 20,
         'LENO': 30, 'ORNL': 2}



    #fig, ax = plt.subplots()
    #dates = []
    #for d in data['end_date']:
    #    dates.append(dt.strptime(d, '%Y-%m-%dT%H:%MZ'))
    #plt.errorbar(x=dates, y=data['d2h'], yerr=data['d2h_sd']*10, linewidth=0.5)
    #plt.scatter(x=list(dates), y=list(data['d2h']), s=list(data['wt']*sizes[i]))
    #plt.title(i)
    #plt.ylabel("δ$^{2}$H of Precipitation")
    #plt.xlabel("Date")
    #for x in ax.xaxis.get_major_ticks():
    #    label = x.label
    #    label.set_fontsize(10)
    #    label.set_rotation('vertical')
    #fig.subplots_adjust(bottom=0.2)
    #plt.show()


print(summary['LENO'])
# None should have 0 weights

for i in precip:
    data = summary[i]
    daily = precip_daily[i]
    #calculate mean 30 min precipitation amount for each month from daily
        #create a new column for month, then group by month and calculate mean
    #calclate mean 30 min precipitation amount for each sampling period (wt)
        #for any that are zero, make 'wt' the same as the mean for that month
    #then use the equation W = ( 1+ (Qi - Qm) / Qm) (Qm / N)

