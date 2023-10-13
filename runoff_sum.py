import pandas as pd
import numpy as np
import math
from scipy.stats import sem
import statistics as stats
import os

site_data = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\site_data.csv')

discharge_sum = pd.DataFrame(data=[], columns=['site', 'Qtot (mm)', 'Qtot_se', 'missing_months'])

sites = ['ARIK', 'BIGC', 'BLDE', 'BLUE', 'BLWA', 'COMO', 'FLNT', 'HOPB', 'KING', 'LECO', 'LEWI',
            'MART', 'MAYF', 'MCDI', 'MCRA', 'POSE', 'PRIN', 'REDB', 'SYCA', 'TECR', 'TOMB', 'WALK', 'WLOU']

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

maindir = r'C:\Users\User\Downloads\NEON_discharge-field\NEON_discharge-field'

for i in sites:
    index = list(site_data['SITEID']).index(i)
    area = site_data['area_km_sq'][index]
    month = []
    discharge_sample = [[] for _ in range(12)]
    discharge_monthly = [0] * 12
    discharge_monthly_se_sq = [0] * 12
    for dirname in os.listdir(maindir):
        if i in dirname:
            for fname in os.listdir(maindir + '\\' + dirname):
                if 'dsc_fieldData' in fname:
                    data = pd.read_csv(maindir + '\\' + dirname + '\\' + fname)
                    for j in range(len(data['totalDischarge'])):
                        month = data['startDate'][j].split('-')[1]
                        if data['totalDischargeUnits'][j] == 'litersPerSecond':
                            discharge_sample[int(month) - 1].append((data['totalDischarge'][j] / 1000) * 86400 / (1000000 * area) * 1000)
                        elif data['totalDischargeUnits'][j] == 'cubicMetersPerSecond':
                            discharge_sample[int(month) - 1].append(data['totalDischarge'][j] * 86400 / (1000000 * area) * 1000)

    count = 0
    for m in range(12):
        if len(discharge_sample[m]) != 0:
            discharge_monthly[m] = stats.mean(discharge_sample[m]) * days_in_month[m]
            discharge_monthly_se_sq[m] = (sem([x for x in discharge_sample[m] if np.isnan(x) == False])) ** 2 * days_in_month[m]
        else:
            count += 1
    discharge_annual_mm = sum(discharge_monthly)
    discharge_annual_se = math.sqrt(sum([x for x in discharge_monthly_se_sq if np.isnan(x) == False]))
    discharge_sum.loc[len(discharge_sum.index)] = [i, discharge_annual_mm, discharge_annual_se, count]

with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
    print(discharge_sum)

#discharge_sum.to_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\discharge_sum.csv')




