import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
import seaborn as sn
from sklearn.ensemble import RandomForestRegressor
import math
from textwrap import wrap
import matplotlib
import matplotlib.cm as cm
import statsmodels.api as sm
import scipy.stats as stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.feature_selection import SequentialFeatureSelector as sfs

# read data and print column names
df = pd.read_csv(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\model_data.csv')
#df = df.drop(labels='f_ET_from_Ps_se', axis=1)
#df = df.drop(labels='f_Ps_to_ET_se', axis=1)

#df = df.rename(columns={'f_ET_from_Ps': 'Fraction of ET from Ps', 'f_Ps_to_ET': 'Fraction of Ps to ET',
#           'Ptot': 'Mean Annual Precipitation', 'WS_area_km_sq': 'Watershed Area',
#           'ratio_Ps_to_P': 'Ratio of Ps to Mean Annual Precipitation',
#           'ratio_ET_to_P': 'Ratio of ET to Mean Annual Precipitation', 'meanT_C': 'Mean Annual Temperature',
#           'elevation_m': 'Elevation', 'GEDI_pct_above_5m': 'Percent of Canopy Above 5m',
#           'GEDI_pct_above_10m': 'Percent of Canopy Above 10m', 'GEDI_pct_above_20m': 'Percent of Canopy Above 20m',
#           'GEDI_mean': 'Mean Canopy Height', 'GEDI_stdev': 'Standard Deviation of Canopy Height',
#           'EVImin': 'Minimum EVI', 'EVIamp': 'EVI Amplitude', 'EVIarea': 'EVI Area Under the Curve',
#           'seasonLength': 'Growing Season Length', 'ratio_EVI_area_min': 'Ratio of EVI Area to EVI Minimum'})


#'''
# use seaborn to plot a correlation matrix for all but the first column

plt.figure(figsize=(16, 16))
sn.set(font_scale=1.6)
sn.heatmap(df[1:].corr(), cmap='seismic', annot=True, fmt='.2f', annot_kws={"size":12})
plt.subplots_adjust(bottom=0.32, left=0.32)
#plt.savefig(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\corrplot.svg', dpi=500)
plt.show()

#'''
#print(df.iloc[:, :5].sort_values(by=['f_ET_from_Ps']))
cmap='viridis'
vars = [5, 15, 16]


predictions = []
for i in range(len(df)):
    df_cv = df.drop(i, axis=0)
    wls = sm.WLS(df_cv['f_ET_from_Ps'], sm.add_constant(df_cv.iloc[:, vars]), weights=1/df_cv['f_ET_from_Ps_se']).fit()
    pars = wls.params.tolist()
    loocv = df.iloc[i, vars].values.tolist()
    predictions.append(pars[0] + pars[1] * loocv[0] + pars[2] * loocv[1] + pars[3] * loocv[2])

plt.figure(figsize=(8, 6))
plt.plot([0, 1], [0, 1])
plt.scatter(df['f_ET_from_Ps'], predictions, s=80, c=df['latitude'], cmap=cmap)
plt.errorbar(x=df['f_ET_from_Ps'], y=predictions, xerr=df['f_ET_from_Ps_se'], ecolor='black', fmt='none')
plt.title('Fraction of ET from Summer Precipitation')
plt.ylabel('Predicted', fontsize=20)
plt.xlabel('Observed', fontsize=20)
plt.xlim([-1.5, 1.5])
plt.text(2.15, 0.4, 'Latitude', fontsize=20, rotation=270)
plt.colorbar()
#plt.savefig(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\f_ET_from_Ps_loocv.svg', dpi=500)
plt.show()


diff = 0
diff_sq = 0
for i in range(len(predictions)):
    diff += abs(df['f_ET_from_Ps'][i] - predictions[i])
    diff_sq += (df['f_ET_from_Ps'][i] - predictions[i])**2
mae1 = 1/len(predictions) * diff
rmse1 = math.sqrt(diff_sq / len(predictions))
print(mae1, rmse1)

predictions = []
for i in range(len(df)):
    df_cv = df.drop(i, axis=0)
    wls = sm.WLS(df_cv['f_Ps_to_ET'], sm.add_constant(df_cv.iloc[:, vars]), weights=1/df_cv['f_Ps_to_ET_se']).fit()
    pars = wls.params.tolist()
    loocv = df.iloc[i, vars].values.tolist()
    predictions.append(pars[0] + pars[1] * loocv[0] + pars[2] * loocv[1] + pars[3] * loocv[2])

plt.figure(figsize=(8, 6))
ax = plt.axes()
plt.plot([0, 1], [0, 1])
plt.errorbar(x=df['f_Ps_to_ET'], y=predictions, xerr=df['f_Ps_to_ET_se'], ecolor='black', fmt='none')
plt.scatter(df['f_Ps_to_ET'], predictions, s=80, c=df['latitude'], cmap=cmap)
plt.title('Fraction of Summer Precipitation to ET')
plt.ylabel('Predicted', fontsize=20)
plt.xlabel('Observed', fontsize=20)
plt.xlim([-1.5, 1.5])
plt.text(2.15, 0.455, 'Latitude', fontsize=20, rotation=270)
plt.colorbar()

#plt.savefig(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\f_Ps_to_ET_loocv.svg', dpi=500)
plt.show()

diff = 0
diff_sq = 0
for i in range(len(predictions)):
    diff += abs(df['f_Ps_to_ET'][i] - predictions[i])
    diff_sq += (df['f_Ps_to_ET'][i] - predictions[i]) ** 2
mae2 = 1/len(predictions) * diff
rmse2 = math.sqrt(diff_sq / len(predictions))
print(mae2, rmse2)


# Build step forward feature selection
#sfs1 = sfs(estimator=clf, n_features_to_select=3, direction='forward', scoring='neg_mean_absolute_error', cv=23)

fig, ax = plt.subplots(2, len(vars), figsize=(12, 10))

ax[0, 0].set_ylabel('\n'.join(wrap('Fraction of Summer Precipitation to ET (unitless)', 30)))
ax[1, 0].set_ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitless)', 30)))

for i in range(len(vars)):
    ols = sm.OLS(df['f_Ps_to_ET'], sm.add_constant(df.iloc[:, vars[i]])).fit()
    ax[0, i].axhline(y=1, color="grey", linestyle="dashed", zorder=1)
    ax[0, i].axhline(y=0, color="grey", linestyle="dashed", zorder=1)
    ax[0, i].plot(df.iloc[:, vars[i]], ols.fittedvalues)
    ax[0, i].errorbar(x=df.iloc[:, vars[i]], y=df['f_Ps_to_ET'], yerr=df['f_Ps_to_ET_se'], fmt='none',
                      ecolor='black', elinewidth=1, zorder=3)
    ax[0, i].scatter(x=df.iloc[:, vars[i]], y=df['f_Ps_to_ET'], s=80, c=df['latitude'], cmap=cmap, zorder=2)
    ax[0, i].set_ylim([-1.5, 1.5])

    ols = sm.OLS(df['f_ET_from_Ps'], sm.add_constant(df.iloc[:, vars[i]])).fit()
    ax[1, i].axhline(y=1, color="grey", linestyle="dashed", zorder=1)
    ax[1, i].axhline(y=0, color="grey", linestyle="dashed", zorder=1)
    ax[1, i].plot(df.iloc[:, vars[i]], ols.fittedvalues)
    ax[1, i].errorbar(x=df.iloc[:, vars[i]], y=df['f_ET_from_Ps'], yerr=df['f_ET_from_Ps_se'], fmt='none',
                      ecolor='black', elinewidth=1, zorder=3)
    ax[1, i].scatter(x=df.iloc[:, vars[i]], y=df['f_ET_from_Ps'], s=80, c=df['latitude'], cmap=cmap, zorder=2)
    ax[1, i].set_ylim([-1.5, 1.5])

ax[1, 0].set_xlabel('\n'.join(wrap('Mean Annual Precipitation (mm)', 18)))
ax[1, 1].set_xlabel('\n'.join(wrap('Standard Deviation of Canopy Height (m)', 20)))
ax[1, 2].set_xlabel('\n'.join(wrap('Minimum EVI (unitless)', 12)))

ax[0, 0].text(2500, 1.2, 'A')
ax[0, 1].text(14, 1.2, 'B')
ax[0, 2].text(2700, 1.2, 'C')

ax[1, 0].text(2500, 1.2, 'D')
ax[1, 1].text(14, 1.2, 'E')
ax[1, 2].text(2700, 1.2, 'F')

fig.subplots_adjust(bottom=0.2, right=0.87, wspace=0.3)
cb_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
norm = plt.Normalize(min(df['latitude']), max(df['latitude']))
cb = fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), cax=cb_ax, shrink=0.8)
cb.set_label("Latitude", rotation=270, labelpad=20)
yax = cb.ax.get_yaxis()

yax.set_ticks([32, 34, 36, 38, 40, 42, 44])
yax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
cb.ax.minorticks_off()

#plt.savefig(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\top_predictors.svg', dpi=500)
plt.show()
quit()
print(round(stats.linregress(df['f_P_from_Ps'] * df['ptot_mm'], df['EVI_min'])[2]**2, 4))
print(round(stats.linregress(df['f_P_from_Ps'], df['EVI_min'])[2]**2, 4))
print(round(stats.linregress(df['ptot_mm'], df['EVI_min'])[2]**2, 4))
print('column   f_Ps_to_ET   f_ET_from_Ps')
for i in range(5, len(df.columns)):
    print(df.columns[i], round(stats.linregress(df.iloc[:, i], df['f_Ps_to_ET'])[2]**2, 4), round(stats.linregress(df.iloc[:, i], df['f_ET_from_Ps'])[2]**2, 4))


'''
def get_stats(x_cols, y):
    x = df[x_cols]
    results = sm.OLS(y, x).fit()
    print(results.summary())

x_cols.remove('EVI_area_mean')
x_cols.remove('EVI_amp_mean')
x_cols.remove('Peak_late_CCI')
x_cols.remove('st_dev_med_CCI')
x_cols.remove('peak_CCI')
x_cols.remove('min_CCI')
x_cols.remove('EVI_min')
x_cols.remove('forest_height_mean')
get_stats(x_cols, y)

x = df[x_cols]
linear_model = LinearRegression()

linear_model.fit(x, y)
predictions = []
for i in range(len(df)):
    df_cv = df.drop(i, axis=0)
    linear_model.fit(df_cv.iloc[:, [2, 3, 4, 9]], df_cv['f_ET_from_Ps'])
    loocv = df.iloc[i, [2, 3, 4, 9]].values.reshape(1, -1)
    print(i, loocv)
    predictions.append(linear_model.predict(loocv))

plt.scatter(x=predictions, y=df['f_ET_from_Ps'])
for i, txt in enumerate(df['siteID']):
    plt.annotate(txt, (predictions[i], df['f_ET_from_Ps'][i]))
plt.plot([0, 0.8], [0, 0.8])
plt.xlabel('Predicted')
plt.ylabel('Observed')
plt.show()






# separate training and testing data based on 7 randomly selected indices
random.seed(1)
test = random.choice(np.arange(len(df)), size=7, replace=False)
df_train = df.drop(test, axis=0)
df_test = df.iloc[test]

# set hyperparameters for random forest algorithm
rf = RandomForestRegressor(n_estimators=100, random_state=2)

# conduct leave-one-out cross validation for model selection
print(df_train.iloc[:, [2, 3, 4, 9]])
df_train.reset_index(drop=True, inplace=True)
print(df_train)
predictions = []
for i in range(len(df_train)):
    df_cv = df_train.drop(i, axis=0)
    rf.fit(df_cv.iloc[:, [2, 3, 4, 9]], df_cv['f_ET_from_Ps'])
    loocv = df_train.iloc[i, [2, 3, 4, 9]].values.reshape(1, -1)
    print(i, loocv)
    predictions.append(rf.predict(loocv))

plt.scatter(x=predictions, y=df_train['f_ET_from_Ps'])
for i, txt in enumerate(df_train['siteID']):
    plt.annotate(txt, (predictions[i], df_train['f_ET_from_Ps'][i]))
plt.plot([0, 0.8], [0, 0.8])
plt.xlabel('Predicted')
plt.ylabel('Observed')
plt.show()


rf.fit(df_train.iloc[:, [2, 3, 4, 9]], df_train['f_ET_from_Ps'])
predictions = rf.predict(df_test.iloc[:, [2, 3, 4, 9]])
errors = abs(predictions - df_test['f_ET_from_Ps'])
print(errors)
'''


