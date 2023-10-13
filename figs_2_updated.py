import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
import seaborn as sn
from sklearn.ensemble import RandomForestRegressor
import math
from textwrap import wrap
import statsmodels.api as sm
import scipy.stats as stats
from sklearn.feature_selection import SequentialFeatureSelector as sfs

# read data and print column names
df = pd.read_csv(r'C:\Users\User\Documents\UNR\Summer2023\IsotopeMapping\model_data.csv')
df_corr = df.drop(labels='f_ET_from_Ps_se', axis=1)
df_corr = df_corr.drop(labels='f_Ps_to_ET_se', axis=1)

df_corr = df_corr.drop(['WS_area_km_sq', 'elevation_m', 'GEDI_pct_above_5m', 'GEDI_pct_above_10m', 'GEDI_pct_above_20m',
                     'EVIamp', 'EVIarea', 'ratio_EVI_area_min', 'latitude'], axis=1)
df_corr = df_corr.rename(columns={'f_ET_from_Ps': 'Fraction of ET from Ps', 'f_Ps_to_ET': 'Fraction of Ps to ET',
           'Ptot': 'Mean Annual Precipitation (MAP)',
           'ratio_Ps_to_P': 'Ratio of Ps to MAP',
           'ratio_ET_to_P': 'Ratio of ET to MAP', 'meanT_C': 'Mean Annual Temperature',
           'GEDI_mean': 'Mean Canopy Height', 'GEDI_stdev': 'Standard Deviation of Canopy Height',
           'EVImin': 'Minimum EVI', 'seasonLength': 'Growing Season Length'})
print(df_corr.columns)

#'''
# use seaborn to plot a correlation matrix for all but the first column

plt.rcParams['figure.facecolor'] = 'white'
plt.figure(figsize=(14, 12))
sn.set_theme(style='white')
sn.set(font_scale=1.6)
sn.heatmap(df_corr[1:].corr(), cmap='seismic', annot=True, fmt='.2f', annot_kws={"size":14})
plt.subplots_adjust(bottom=0.4, left=0.4)
plt.savefig(r'C:\Users\User\Documents\UNR\Summer2023\IsotopeMapping\poster\corrplot_simple.svg', dpi=500)
plt.show()

df_supp = pd.read_csv(r'C:\Users\User\Documents\UNR\Summer2023\IsotopeMapping\code_input_f_ET_f_Ps.csv')
df['ratio_Ps_to_P_se'] = df_supp['ratio_se']
df = df.sort_values('f_ET_from_Ps')
df['rank'] = list(range(23, 0, -1))


plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['axes.axisbelow'] = 'line'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.spines.left'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False


plt.figure(figsize=(5, 11))
plt.axvline(0, c='black', linestyle='dashed')
plt.axvline(1, c='black', linestyle='dashed')
plt.scatter(data=df, x='ratio_Ps_to_P', y='rank', facecolors='none', edgecolors='dodgerblue', s=40, marker='o', label='\n'.join(wrap('Ratio of Summer Precipitation to Total Precipitation', 20)))
plt.scatter(data=df, x='f_ET_from_Ps', y='rank', color='green', s=40, label='\n'.join(wrap('Fraction of ET from Summer Precipitation', 20)))
for i in range(len(df['f_ET_from_Ps'])):
    plt.plot((df['ratio_Ps_to_P'][i] - df['ratio_Ps_to_P_se'][i], df['ratio_Ps_to_P'][i] +
              df['ratio_Ps_to_P_se'][i]), (df['rank'][i], df['rank'][i]), color='dodgerblue',
             linewidth=1)
    plt.plot((df['f_ET_from_Ps'][i] - df['f_ET_from_Ps_se'][i], df['f_ET_from_Ps'][i] +
              df['f_ET_from_Ps_se'][i]), (df['rank'][i], df['rank'][i]), color='green', linewidth=3)


plt.xlim([-1, 1.2])
plt.ylim([0.5, 23.2])
plt.margins(y=3)
plt.legend(bbox_to_anchor=[0.62, -0.14], loc='center', frameon=False)
plt.tick_params(axis='both', labelsize=14, labelleft=False)
plt.xticks([-1, -0.5, 0, 0.5, 1], fontsize=20)

plt.subplots_adjust(bottom=0.25)
plt.tight_layout()
for i in range(len(df)):
    print(df['rank'][i], df['site'][i], df['f_ET_from_Ps'][i])

#plt.savefig(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\f_ET_sites.svg', dpi=500)
plt.show()


