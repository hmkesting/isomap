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
df = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\model_data.csv')
plt.hist(df['EVI_min'])
df = df.drop(labels='f_ET_from_Ps_se', axis=1)
df = df.drop(labels='f_Ps_to_ET_se', axis=1)
df = df.drop(labels='latitude', axis=1)
df.columns = ['SiteID', 'Fraction of ET from Ps', 'Fraction of Ps to ET', 'Coefficient of Variation of CCI', 'Standard Deviation of CCI', 'Peak - Late Season CCI', 'Peak CCI', 'Minimum CCI', '\n'.join(wrap('Ratio of Summer Precip. to Annual Precip.', 30)), 'Mean Annual Precip.', 'Forest Height', 'EVI Amplitude', 'EVI Area-Under-the-Curve', 'EVI Minimum']
#'''
# use seaborn to plot a correlation matrix for all but the first column

plt.figure(figsize=(11, 11))
sn.set(font_scale=1.6)
sn.heatmap(df[1:].corr(), cmap='seismic', annot=True, fmt='.2f', annot_kws={"size":12})
plt.subplots_adjust(bottom=0.32, left=0.32)
#plt.savefig(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\corrplot.svg', dpi=500)
plt.show()
quit()
#'''
print(df.iloc[:, :5].sort_values(by=['f_ET_from_Ps']))
cmap='viridis'

predictions = []
for i in range(len(df)):
    df_cv = df.drop(i, axis=0)
    wls = sm.WLS(df_cv['f_ET_from_Ps'], sm.add_constant(df_cv.iloc[:, [10, 15]]), weights=1/df_cv['f_ET_from_Ps_se']).fit()
    pars = wls.params.tolist()
    loocv = df.iloc[i, [10, 15]].values.tolist()
    predictions.append(pars[0] + pars[1] * loocv[0] + pars[2] * loocv[1])

plt.figure(figsize=(8, 6))
plt.plot([0, 1], [0, 1])
plt.scatter(df['f_ET_from_Ps'], predictions, c=df['latitude'], cmap=cmap)
plt.errorbar(x=df['f_ET_from_Ps'], y=predictions, xerr=df['f_ET_from_Ps_se'], ecolor='black', fmt='none')
plt.title('Fraction of ET from Summer Precipitation')
plt.ylabel('Predicted', fontsize=11)
plt.xlabel('Observed', fontsize=11)
plt.text(1.85, 0.4, 'Latitude', fontsize=11, rotation=270)
plt.colorbar()
#plt.savefig(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\f_ET_from_Ps_loocv.svg', dpi=500)
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
    wls = sm.WLS(df_cv['f_Ps_to_ET'], sm.add_constant(df_cv.iloc[:, [5, 11, 14, 15]]), weights=1/df_cv['f_Ps_to_ET_se']).fit()
    pars = wls.params.tolist()
    loocv = df.iloc[i, [5, 11, 14, 15]].values.tolist()
    predictions.append(pars[0] + pars[1] * loocv[0] + pars[2] * loocv[1] + pars[3] * loocv[2] + pars[4] * loocv[3])

plt.figure(figsize=(8, 6))
plt.plot([0, 1], [0, 1])
plt.errorbar(x=df['f_Ps_to_ET'], y=predictions, xerr=df['f_Ps_to_ET_se'], ecolor='black', fmt='none')
plt.scatter(df['f_Ps_to_ET'], predictions, c=df['latitude'], cmap=cmap)
plt.title('Fraction of Summer Precipitation to ET')
plt.ylabel('Predicted', fontsize=11)
plt.xlabel('Observed', fontsize=11)
plt.text(1.89, -0.02, 'Latitude', fontsize=11, rotation=270)
plt.colorbar()
#plt.savefig(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\f_Ps_to_ET_loocv.svg', dpi=500)
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


fig, ax = plt.subplots(2, 2)
ax[0, 0].set_ylabel('\n'.join(wrap('Fraction of Summer Precipitation to ET (unitless)', 30)))
ols = sm.OLS(df['f_Ps_to_ET'], sm.add_constant(df.iloc[:, 11])).fit()
ax[0, 0].plot(df.iloc[:, 11], ols.fittedvalues)
ax[0, 0].scatter(x=df.iloc[:, 11], y=df['f_Ps_to_ET'], s=3/df['f_Ps_to_ET_se'] + 1, c=df['latitude'], cmap=cmap)
ax[0, 0].set_xlabel('Mean Annual Precipitation (mm)')
ols = sm.OLS(df['f_Ps_to_ET'], sm.add_constant(df.iloc[:, 15])).fit()
ax[0, 1].plot(df.iloc[:, 15], ols.fittedvalues)
ax[0, 1].scatter(x=df.iloc[:, 15], y=df['f_Ps_to_ET'], s=3/df['f_Ps_to_ET_se'] + 1, c=df['latitude'], cmap=cmap)
ax[0, 1].set_xlabel('Minimum EVI (unitless)')

ax[1, 0].set_ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitless)', 30)))
ols = sm.OLS(df['f_ET_from_Ps'], sm.add_constant(df.iloc[:, 10])).fit()
ax[1, 0].plot(df.iloc[:, 10], ols.fittedvalues)
ax[1, 0].scatter(x=df.iloc[:, 10], y=df['f_ET_from_Ps'], s=1/df['f_ET_from_Ps_se'] + 2, c=df['latitude'], cmap=cmap)
ax[1, 0].set_xlabel('\n'.join(wrap('Ratio of Summer Precipitation to Annual Precipitation', 30)))
ols = sm.OLS(df['f_ET_from_Ps'], sm.add_constant(df.iloc[:, 15])).fit()
ax[1, 1].plot(df.iloc[:, 15], ols.fittedvalues)
ax[1, 1].scatter(x=df.iloc[:, 15], y=df['f_ET_from_Ps'], s=1/df['f_ET_from_Ps_se'] + 2, c=df['latitude'], cmap=cmap)
ax[1, 1].set_xlabel('Minimum EVI (unitless)')
ax[0, 0].text(2550, 1, 'A')
ax[0, 1].text(2600, 0.98, 'B')
ax[1, 0].text(0.73, 1, 'C')
ax[1, 1].text(2600, 0.98, 'D')
plt.tight_layout()
#plt.savefig(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\top_predictors.svg', dpi=500)
plt.show()
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
'''




'''
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


