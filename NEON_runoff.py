import pandas as pd
import numpy as np
import math
from scipy.stats import sem
import statistics as stats
import os
import copy
import pickle


# Calculate annual runoff from field data and output annual runoff (mm), standard error, and 
# the count of months missing stream data for each site

# Areas calculated from NEON Aquatic Watershed boundary
site_area_sq_km = {'ARIK': 2632.3, 'BIGC': 10.9009, 'BLDE': 37.839, 
                   'BLUE': 322.15, 'BLWA': 16171.6, 'COMO': 3.56629, 
                   'FLNT': 14986.7, 'HOPB': 11.882, 'KING': 13.0413, 
                   'LECO': 9.12099, 'LEWI': 11.9086, 'MART': 6.34773, 
                   'MAYF': 14.3621, 'MCDI': 22.5796, 'MCRA': 3.93222, 
                   'POSE': 2.02731, 'PRIN': 48.9114, 'REDB': 16.7148, 
                   'SYCA': 280.475, 'TECR': 2.97479, 'TOMB': 47109.9, 
                   'WALK': 1.08996, 'WLOU': 4.90348}
sites = list(site_area_sq_km.keys())
days_in_month = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

Qtot = {}
Qtot_se = {}
Qtot_missing_mo = {}

maindir = r'C:\Users\User\Documents\UNR\NEON Project\Data\NEON\NEON_discharge-field'
def calc_Qtot(site):
    runoff_sample = [[] for _ in range(12)]
    runoff_monthly = [0] * 12
    runoff_monthly_se_sq = [0] * 12
    # Compile all data for each month for each site and convert to mm per day over the watershed area for each month
    for dirname in os.listdir(maindir):
        if site in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if 'dsc_fieldData' in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname)
                    for j in range(len(data['totalDischarge'])):
                        month = data['startDate'][j].split('-')[1]
                        if data['totalDischargeUnits'][j] == 'litersPerSecond':
                            runoff_sample[int(month) - 1].append((data['totalDischarge'][j] / 1000) * 86400 / (1000000 * site_area_sq_km[i]) * 1000)
                        elif data['totalDischargeUnits'][j] == 'cubicMetersPerSecond':
                            runoff_sample[int(month) - 1].append(data['totalDischarge'][j] * 86400 / (1000000 * site_area_sq_km[i]) * 1000)

    # Calculate monthly runoff and standard error from the list of samples, count months with <2 values
    count_missing = 0
    count_one_sample = 0
    for m in range(12):
        if len(runoff_sample[m]) > 1:
            runoff_monthly[m] = stats.mean(runoff_sample[m]) * days_in_month[m]
            runoff_monthly_se_sq[m] = (sem([x for x in runoff_sample[m] if np.isnan(x) == False]) ** 2) * days_in_month[m]
        else:
            if len(runoff_sample[m]) == 1:
                runoff_monthly[m] = runoff_sample[m][0]
                count_one_sample += 1
            if len(runoff_sample[m]) == 0:
                count_missing += 1
            
    # Calculate annual runoff and standard error by summing monthly runoff and error propagation
    runoff_annual_mm = sum(runoff_monthly)
    runoff_annual_se = math.sqrt(sum(runoff_monthly_se_sq))
    return [(runoff_annual_mm, runoff_annual_se), count_missing, count_one_sample]

for i in sites:
    Qtot[i], Qtot_missing_mo[i], Qtot_mo_one_measurement = calc_Qtot(i)
    
with open('NEON_Qtot.pickle', 'wb') as f:
    pickle.dump(Qtot, f)
    
# Compile NEON isotope surface water data into a table with date, sample isotope value, month,
# days_in_month and include sample runoff and average monthly runoff from the Rhea dataset

Qdel = {}
runoff_isotope_samples_df = {}

maindir = r'C:\Users\User\Documents\UNR\NEON Project\Data\NEON\NEON_isotope-surfacewater'
def calc_Qdeli(site):
    runoff_isotopes = []
    runoff_isotope_dates = []
    months = []
    for dirname in os.listdir(maindir):
        if site in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if 'externalLabH2OIsotopes' in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname)
                    for j in range(len(data['d2HWater'])):
                        runoff_isotopes.append(data['d2HWater'][j])
                        runoff_isotope_dates.append(data['collectDate'][j])  # date format 2016-01-01
                        months.append(data['collectDate'][j].split('-')[1])
    df = pd.DataFrame({'Qdel_date': runoff_isotope_dates, 'Qdeli': runoff_isotopes, 'month': months})
    return df

for i in sites:
    runoff_isotope_samples_df[i] = calc_Qdeli(i)

monthly_runoff = {}
def calc_Q(daily_q, isotope_df, column):    
    isotope_sample_q = [np.nan]*len(isotope_df)
    month = []
    for i in range(len(daily_q)):   # date format 2016-01-01
        date = daily_q.index[i]
        month.append(date.split('-')[1])
        if date in isotope_df['Qdel_date'].unique():
            index = list(isotope_df['Qdel_date']).index(date)
            isotope_sample_q[index] = daily_q.loc[date].iloc[0]
    daily_q['month'] = month
    monthly_q = pd.DataFrame(daily_q.groupby(by=["month"], as_index=False)[column].mean())
    month_counts = []      # count the number of isotope samples in each month
    for m in range(1,13):    
        month_counts.append([int(i.split('-')[1]) for i in isotope_df['Qdel_date']].count(m))
    monthly_q['sample_count'] = month_counts
    monthly_q['days_in_month'] = days_in_month
    isotope_df['Qi'] = isotope_sample_q
    for i in range(len(isotope_df['Qi'])):
        if np.isnan(isotope_df['Qi'].iloc[i]):
            isotope_df.loc[i, 'Qi'] = monthly_q.iloc[monthly_q.index == isotope_df.loc[i, 'month']]
    return monthly_q, isotope_df
  
for site in [x for x in sites if x not in ['TECR', 'BIGC', 'WLOU']]:
    daily_runoff_predictions = pd.read_csv(
        r'C:\Users\User\Documents\UNR\NEON Project\Data\Rhea\modeled_discharge_2022-10-13\q_lm_outdata\predictions\\' 
        + site + '.csv').replace(-9999.0, np.nan)
    if daily_runoff_predictions.columns[0] == 'datetime':
        daily_runoff_predictions['date'] = [d[:10] for d in daily_runoff_predictions['datetime']]
    daily_runoff_predictions = pd.DataFrame(daily_runoff_predictions.groupby(by=["date"])['Q_predicted'].mean())
    monthly_runoff[site], runoff_isotope_samples_df[site] = calc_Q(daily_runoff_predictions, runoff_isotope_samples_df[site], 'Q_predicted')

for site in ['TECR']:
    daily_runoff_predictions = pd.read_csv(
        r'C:\Users\User\Documents\UNR\NEON Project\Data\Rhea\modeled_discharge_2022-10-13\lstm_outdata\predictions\\'
        + site +'.csv').replace(-9999.0, np.nan)
    daily_runoff_predictions = pd.DataFrame(daily_runoff_predictions.groupby(by=["date"])['Q_predicted'].mean())
    monthly_runoff[site], runoff_isotope_samples_df[site] = calc_Q(daily_runoff_predictions, runoff_isotope_samples_df[site], 'Q_predicted')  

for site in ['BIGC', 'WLOU']:
    daily_runoff = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Rhea\StreamData.csv').replace(-9999.0, np.nan).loc[:,['Date', site]] 
    daily_runoff = daily_runoff.dropna(subset=daily_runoff.columns.values).rename(columns={site: 'Q'})
    daily_runoff['date'] = [d[:10] for d in daily_runoff['Date']] # date format 2015-01-01
    daily_runoff = pd.DataFrame(daily_runoff.groupby(by=["date"])['Q'].mean())
    monthly_runoff[site], runoff_isotope_samples_df[site] = calc_Q(daily_runoff, runoff_isotope_samples_df[site], 'Q')    


def remove_spikes(df, column):
    df = copy.deepcopy(df)
    while True:   # repeat until there are no outliers, then break
        length = int(len(df))
        mean = np.nanmean(df[column])
        sd = np.nanstd(df[column])
        # if any values are further than 3 standard deviations
        # away from the mean, remove the row
        df = df[(df[column] >= (mean - 3 * sd)) & (df[column] <= (mean + 3 * sd))]
        if len(df) == length:
            break
    df = df.reset_index()
    return df

for site in sites:
    runoff_isotope_samples_df[site] = remove_spikes(runoff_isotope_samples_df[site], 'Qdeli')


# Calculate Qdel_wts from equation (1+(Qi-Qm)/Qm)(Qm*D/N), the calculate Qdel weighted mean
Qdel = {}

def calc_Qdel_wts(row):
    wt = (1 + (row['Qi'] - row['Qm']) / row['Qm']) * (row['Qm'] * row ['D'] / row['N'])
    return wt

for site in sites:
    if site in ['BIGC', 'WLOU']:
        column = 'Q'
    else:
        column = 'Q_predicted'
    mo = runoff_isotope_samples_df[site]['month'].values
    runoff_isotope_samples_df[site]['Qm'] = [monthly_runoff[site].loc[monthly_runoff[site]['month'] == x, column].iloc[0] for x in mo]
    runoff_isotope_samples_df[site]['N'] = [monthly_runoff[site].loc[monthly_runoff[site]['month'] == x, 'sample_count'].iloc[0] for x in mo]
    runoff_isotope_samples_df[site]['D'] = [monthly_runoff[site].loc[monthly_runoff[site]['month'] == x, 'days_in_month'].iloc[0] for x in mo]
    runoff_isotope_samples_df[site]['Qdel_wt'] = runoff_isotope_samples_df[site].apply(calc_Qdel_wts, axis=1)
    
def wtd_mean(x, wt=None):
    if wt is None:
        wt = [1] * len(x)
    remove = []
    for item in range(len(x)):
        if math.isnan(x[item]) or math.isnan(wt[item]):
            remove.append(item)
    for item in sorted(remove, reverse=True):
        del x[item]
        del wt[item]
    if len(x) != len(wt):
        raise Exception("error in wtd_mean: x and wt have different lengths")
    for item in range(len(x)):
        if wt[item] < 0:
            raise Exception("error in wtd_mean: negative weights")
    sumwt = sum(wt)
    sq_list = []
    for item in range(len(x)):
        sq_list.append(wt[item]**2)
    sumsq = sum(sq_list)
    n_eff = (sumwt*sumwt)/sumsq
    xbar_list = []
    for item in range(len(x)):
        xbar_list.append(x[item]*wt[item])
    xbar = sum(xbar_list)/sumwt
    varx_list = []
    for item in range(len(x)):
        varx_list.append(wt[item] * ((x[item] - xbar)**2))
    varx = (sum(varx_list)/sumwt) * n_eff/(n_eff - 1.0)
    return xbar, math.sqrt(varx/n_eff)

for site in sites:
    Qdel[site] = wtd_mean(list(runoff_isotope_samples_df[site]['Qdeli']), list(runoff_isotope_samples_df[site]['Qdel_wt']))

with open('NEON_Qdel.pickle', 'wb') as f:
    pickle.dump(Qdel, f)
   
    
        

