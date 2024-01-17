import pandas as pd
import numpy as np
import statistics as stats
import os
import pickle
from datetime import datetime as dt
from datetime import date, timedelta
from NEON_runoff import wtd_mean


days_in_month = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# Tells if primary is available for precip sites listed above. If primary precipitation is not available, use secondary
primary_precip = {'ARIK': True, 'SJER': True, 'YELL': True, 'BLUE': True, 'DELA': False, 'NIWO': True, 'JERC': False, 
                  'HARV': True, 'KONZ': True, 'GRSM': False, 'BLAN': False, 'WREF': True, 'TALL': False, 
                  'SCBI': True, 'PRIN': True, 'REDB': True, 'SYCA': False, 'LENO': False, 'ORNL': True}
p_sites = list(primary_precip.keys())

maindir = r'C:\Users\User\Documents\UNR\NEON Project\Data\NEON\NEON_precipitation'
'''
precip = {}
for i in p_sites:
    data_summary = pd.DataFrame(columns=['start_date', 'end_date', 'ppt_mm'])
    for dirname in os.listdir(maindir):
        if i in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if primary_precip[i]:
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
                    for x in range(len(data['startDateTime'])):
                        data_summary.loc[len(data_summary.index)] = \
                            [data['startDateTime'][x], data['endDateTime'][x], data[c][x]]
    print(i)
    precip[i] = data_summary

with open('NEON_precip.pickle', 'wb') as f:
    pickle.dump(precip, f)
'''

with open('NEON_precip.pickle', 'rb') as f:
    precip = pickle.load(f)
    
maindir = r'C:\Users\User\Documents\UNR\NEON Project\Data\NEON\NEON_isotope-ppt'
precip_d2h = {}
for i in p_sites:
    data_summary = pd.DataFrame(columns=['start_date', 'end_date', 'd2h', 'd2h_sd'])
    for dirname in os.listdir(maindir):
        if i in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if 'isoPerSample' in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname) 
                    for x in range(len(data['collectDate'])):
                        data_summary.loc[len(data_summary.index)] = \
                        [data['setDate'][x], data['collectDate'][x], data['d2HWater'][x], data['d2HsdWater'][x]]
    precip_d2h[i] = data_summary


def extract_month(row):
    start = row['start_date']
    if len(start) == 17:
        start_date = dt.strptime(row['start_date'], '%Y-%m-%dT%H:%MZ') 
    elif len(start) == 20:
        start_date = dt.strptime(row['start_date'], '%Y-%m-%dT%H:%M:%SZ') 
    else:
        raise Exception('Date format is not covered in extract_month function')
    if start_date.day == 1 and start_date.hour == 0 and start_date.minute == 0 and start_date.month > 1:
        return start_date.month - 1
    elif start_date.day == 1 and start_date.hour == 0 and start_date.minute == 0 and start_date.month == 1:
        return 12
    else:
        return start_date.month

precip_monthly = {}
for i in p_sites:
    precip_d2h[i]['month'] = precip_d2h[i].apply(extract_month, axis=1)
    precip[i]['month'] = precip[i].apply(extract_month, axis=1)
    precip_monthly[i] = pd.DataFrame(precip[i].groupby(by=['month'], as_index=False)['ppt_mm'].mean())
    for m in range(1, 13):
        if not m in list(precip_monthly[i]['month']):
            if m > 1 and m < 12:
                p = (precip_monthly[i]['ppt_mm'][m - 1] + precip_monthly[i]['ppt_mm'][m - 2]) / 2
            elif m == 1 or m == 12:
                p = (precip_monthly[i]['ppt_mm'][0] + \
                     precip_monthly[i]['ppt_mm'][len(precip_monthly[i].index) - 1]) / 2
            precip_monthly[i].loc[len(precip_monthly[i])] = {'month':m, 'ppt_mm':p}
    precip_monthly[i] = precip_monthly[i].sort_values('month').reset_index(drop=True)

# Round date to the nearest 30 minutes 
def round_dt_30min(time):
    tm = dt.strptime(time, '%Y-%m-%dT%H:%MZ')
    if tm.minute <= 15:
        return dt(tm.year, tm.month, tm.day, tm.hour, 0)
        
    elif tm.minute > 15 and tm.minute < 45:
        return dt(tm.year, tm.month, tm.day, tm.hour, 30)
    
    elif tm.minute >= 45:
        if tm.hour == 23:
            day = date(tm.year, tm.month, tm.day) + timedelta(days=1)
            return dt(day.year, day.month, day.day, 0, 0)
        else:
            return dt(tm.year, tm.month, tm.day, tm.hour + 1, 0)
    else:
        raise Exception('round_dt_30min found invalid minute value')

# Calculate precipitation rate for each isotope sample
for i in p_sites:
    Pi = []
    for r in range(len(precip_d2h[i]['start_date'])):
        start = round_dt_30min(precip_d2h[i]['start_date'][r])
        end = round_dt_30min(precip_d2h[i]['end_date'][r])
        precip_mean = []
        for p in range(len(precip[i]['start_date'])):
            if dt.strptime(precip[i]['start_date'][p], '%Y-%m-%dT%H:%M:%SZ') >= start and \
                dt.strptime(precip[i]['end_date'][p], '%Y-%m-%dT%H:%M:%SZ') <= end:
                precip_mean.append(precip[i]['ppt_mm'][p])
        if len([x for x in precip_mean if not np.isnan(x)]) >= 3 and sum([x for x in precip_mean if not np.isnan(x)]) < 0:
            p = stats.mean([x for x in precip_mean if not np.isnan(x)])
        else:
            p = precip_monthly[i]['ppt_mm'][start.month - 1]
        Pi.append(p)
    precip_d2h[i]['Pi'] = Pi

# Calculate N for each month
for i in p_sites:
    monthly_n = precip_d2h[i]['month'].value_counts()
    n_filled = [0]*12
    for j in monthly_n.index:
        n_filled[j - 1] = monthly_n[j]
    precip_monthly[i]['N'] = n_filled

    
for i in p_sites:
    D = []
    N = []
    Pm = []
    for m in precip_d2h[i]['month']:
        D.append(days_in_month[m - 1])
        N.append(precip_monthly[i]['N'][m - 1])
        Pm.append(precip_monthly[i]['ppt_mm'][m - 1])
    precip_d2h[i]['Pm'] = Pm
    precip_d2h[i]['N'] = N
    precip_d2h[i]['D'] = D
    

def calc_Pdel_wts(row):
    wt = (1 + (row['Pi'] - row['Pm']) / row['Pm']) * (row['Pm'] * row ['D'] / row['N'])
    return wt

winter_months = [1, 2, 3, 10, 11, 12]
summer_months = [4, 5, 6, 7, 8, 9]
Pdel_s = {}
Pdel_w = {}
for i in p_sites:
    precip_d2h[i]['Pdel_wt'] = precip_d2h[i].apply(calc_Pdel_wts, axis=1)
    Pdel_s[i] = wtd_mean(list(precip_d2h[i].loc[precip_d2h[i]['month'].isin(summer_months), 'd2h']), \
                         list(precip_d2h[i].loc[precip_d2h[i]['month'].isin(summer_months),'Pdel_wt']))
    Pdel_w[i] = wtd_mean(list(precip_d2h[i].loc[precip_d2h[i]['month'].isin(winter_months), 'd2h']), \
                         list(precip_d2h[i].loc[precip_d2h[i]['month'].isin(winter_months),'Pdel_wt']))

with open('NEON_Pdel_s.pickle', 'wb') as f:
    pickle.dump(Pdel_s, f)
    
with open('NEON_Pdel_w.pickle', 'wb') as f:
    pickle.dump(Pdel_w, f)
    
for i in p_sites:
    precip_d2h[i].to_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\NEON_precipitation\precip_d2h' + i + '.csv.')
    precip_monthly[i].to_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\NEON_precipitation\precip_monthly' + i + '.csv.')
   
        
   
        

    
    