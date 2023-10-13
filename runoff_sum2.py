import pandas as pd
import numpy as np
import math
from scipy.stats import sem
import statistics as stats
import os

site_data = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\site_data.csv')
stream_d2h = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\Stream_d2H.csv')
stream_d2h = stream_d2h.replace(-9999.0, np.nan)

sites = ['ARIK', 'BIGC', 'BLDE', 'BLUE', 'BLWA', 'COMO', 'FLNT', 'HOPB', 'KING', 'LECO', 'LEWI',
            'MART', 'MAYF', 'MCDI', 'MCRA', 'POSE', 'PRIN', 'REDB', 'SYCA', 'TECR', 'TOMB', 'WALK', 'WLOU']

# if all columns of the subset contain NaN for a given row, delete that row
stream_d2h = stream_d2h.dropna(subset=sites, how='all')

discharge_sum = pd.DataFrame(data=[], columns=['site', 'Qtot_ann_mm', 'Qtot_ann_se', 'missing_months'])

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

maindir = r'C:\Users\User\Downloads\NEON_discharge-field\NEON_discharge-field'
for i in sites:
    index = list(site_data['SITEID']).index(i)
    area = site_data['area_km_sq'][index]
    for dirname in os.listdir(maindir):
        if i in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if 'dsc_fieldData' in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname)   # monthly file
                    data = data[data.totalDischarge >= 0]
                    # slice dates to remove time information at the end
                    data['date_only'] = [x[0:10] for x in data['startDate']]
                    month = data['startDate'][0].split('-')[1]
                    # calculate mean discharge rate for each day sampled
                    data.groupby(by=["date_only"])['totalDischarge'].mean()
                    for j in data.index:
                        if data['date_only'][j] in list(stream_d2h['Date']):
                            iso_index = list(stream_d2h['Date']).index(data['startDate'][j][0:10])
                            weights[iso_index] = data['totalDischarge'][j]
                        if data['totalDischargeUnits'][j] == 'litersPerSecond':
                            discharge_sample[int(month) - 1].append((data['totalDischarge'][j] / 1000) * 86400 / (1000000 * area) * 1000)    #should be using mean of daily discharge if discharge is not 0
                        elif data['totalDischargeUnits'][j] == 'cubicMetersPerSecond':
                            discharge_sample[int(month) - 1].append(data['totalDischarge'][j] * 86400 / (1000000 * area) * 1000)

    count = 0
    for m in range(12):
        if len(discharge_sample[m]) != 0:
            discharge_monthly[m] = round(stats.mean(discharge_sample[m]) * days_in_month[m], 3)
            discharge_monthly_se_sq[m] = (sem([x for x in discharge_sample[m] if np.isnan(x) == False])) ** 2 * days_in_month[m]
        else:
            count += 1
    discharge_annual_mm = sum(discharge_monthly)
    discharge_annual_se = math.sqrt(sum([x for x in discharge_monthly_se_sq if np.isnan(x) == False]))
    discharge_sum.loc[len(discharge_sum.index)] = [i, discharge_annual_mm, discharge_annual_se, count]
    stream_d2h_wts[i] = stream_d2h[i]
    stream_d2h_wts[i + '_wt'] = weights


#with pd.option_context('display.max_rows', None,
#                       'display.max_columns', None,
#                       'display.precision', 3,
#                       ):
    #print(discharge_sum)
    #print(stream_d2h_wts)

#discharge_sum.to_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\discharge_sum.csv')




