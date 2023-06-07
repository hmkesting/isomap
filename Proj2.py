import copy
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy
from textwrap import wrap

streamflow = {}
sites = ['ARIK', 'BLDE', 'BLUE', 'BLWA', 'COMO', 'FLNT', 'HOPB', 'KING', 'LECO', 'LEWI', 'MART', 'MAYF', 'MCDI',
         'MCRA', 'POSE', 'PRIN', 'REDB', 'SYCA', 'TOMB', 'WALK']

for i in sites:
    streamflow[i] = pd.read_csv(
        'C:\\Users\\User\\Downloads\\modeled_discharge_2022-10-13\\modeled_discharge_2022-10-13\\q_lm_outdata\\predictions\\' + i + '.csv')
    streamflow[i] = streamflow[i].replace(-9999.0, np.nan)
streamflow['TECR'] = pd.read_csv(
    'C:\\Users\\User\\Downloads\\modeled_discharge_2022-10-13\\modeled_discharge_2022-10-13\\lstm_outdata\\predictions\\TECR.csv')
# Missing BIGC, WLOU, CRAM, LIRO, PRLA, PRPO


streamflow_alt = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\StreamData.csv')
streamflow_alt = streamflow_alt.loc[:, ['Date', 'BIGC', 'WLOU']]
streamflow_alt = streamflow_alt.replace(-9999.0, np.nan)

# CRAM, LIRO, PRLA, and PRPO are lakes


all_cols = ['Date', 'ARIK', 'BIGC', 'BLDE', 'BLUE', 'BLWA', 'COMO', 'FLNT', 'HOPB', 'KING', 'LECO', 'LEWI',
            'MART', 'MAYF', 'MCDI', 'MCRA', 'POSE', 'PRIN', 'REDB', 'SYCA', 'TECR', 'TOMB', 'WALK', 'WLOU']
stream_d2h = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\Stream_d2H.csv')
stream_d2h = stream_d2h.loc[:, all_cols]
stream_d2h = stream_d2h.replace(-9999.0, np.nan)

# Mean of instantaneous predictions for each day
for i in streamflow:
    if streamflow[i].columns[0] == 'datetime':
        date = []
        for j in streamflow[i]['datetime']:
            date.append(j[0:10])
        streamflow[i]['date'] = date
streamflow_daily = {}
for i in streamflow:
    streamflow_daily[i] = streamflow[i].groupby(by=["date"])['Q_predicted'].mean()
'''
for i in streamflow_daily:
    plt.scatter(streamflow_daily[i].index, streamflow_daily[i])
    plt.title(i)
    plt.xlabel('Date')
    plt.ylabel('Runoff')
    plt.show()
'''
# The 4 lowest records in MART should be manually removed because streamflow was not unusually high
# Analyses should be repeated with and without these 4 points
def remove_spikes(stream):
    stream = copy.deepcopy(stream)
    count = 1
    while count > 0:
        mean = np.nanmean(stream)
        sd = np.nanstd(stream)   # standard error
        new_count = 0
        for x in range(len(stream)):
            if stream[x] < mean - 3 * sd or stream[x] > mean + 3 * sd:
                stream[x] = np.nan
                new_count += 1
        if new_count == 0:
            count = 0
    return stream
print(np.count_nonzero(~np.isnan([np.nan, 0, 0.5, np.nan, 0.5])))
for i in all_cols[1:]:
    stream_d2h[i] = remove_spikes(stream_d2h[i])
    print(i, np.count_nonzero(~np.isnan(list(stream_d2h[i]))))
quit()


months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
monthly_q = {}
for i in all_cols[1:]:  # error is from COMO nans
    for j in months:
        q = []
        if i in streamflow_daily.keys():
            indices = streamflow_daily[i].index
            for x in indices:
                if x.split('-')[1] == j:
                    q.append(streamflow_daily[i][x])
            monthly_q[i + j] = np.nanmean(q)
            monthly_q[i + j + 'se'] = np.nanstd(q)
        elif i in streamflow_alt.columns:
            for x in range(len(streamflow_alt[i])):
                if streamflow_alt['Date'][x].split('-')[1] == j:
                    q.append(streamflow_alt[i][x])
            monthly_q[i + j] = np.nanmean(q)
            monthly_q[i + j + 'se'] = np.nanstd(q)
        else:
            monthly_q[i + j] = 1
            monthly_q[i + j + 'se'] = 0

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

def weight_isotopes(stream, i):
    qdel_sample = []
    q_sample = []
    month_indices = []
    mean_daily_q = [0] * 12
    days_in_month = [0] * 12
    for j in range(len(months)):
        mean_daily_q[j] = monthly_q[i + months[j]]
    for x in range(len(stream[i])):
        if i in ['BIGC', 'WLOU']:
            index = list(streamflow_alt['Date']).index(stream['Date'][x])
            if not np.isnan(stream[i][x]) and not np.isnan(streamflow_alt[i][index]):
                qdel_sample.append(stream[i][x])
                q_sample.append(streamflow_alt[i][index])
                month_index = int(stream['Date'][x].split('-')[1]) - 1
                month_indices.append(month_index)
                days_in_month[month_index] += 1
        elif i in ['CRAM', 'LIRO', 'PRLA', 'PRPO']:
            if not np.isnan(stream[i][x]):
                qdel_sample.append(stream[i][x])
                q_sample.append(1)
                month_index = int(stream['Date'][x].split('-')[1]) - 1
                month_indices.append(month_index)
                days_in_month[month_index] += 1
        else:
            if stream['Date'][x] in list(streamflow_daily[i].index):
                index = list(streamflow_daily[i].index).index(stream['Date'][x])
                if not np.isnan(stream[i][x]) and not np.isnan(streamflow_daily[i][index]):
                    qdel_sample.append(stream[i][x])
                    q_sample.append(streamflow_daily[i][index])
                    month_index = int(stream['Date'][x].split('-')[1]) - 1
                    month_indices.append(month_index)
                    days_in_month[month_index] += 1
    wts = []
    for x in range(len(qdel_sample)):
        mean_q = mean_daily_q[month_indices[x]]
        wts.append((1 + (q_sample[x] - mean_q) / mean_q) * (mean_q / days_in_month[month_indices[x]]))
    return wtd_mean(qdel_sample, wts)


annual_Q = {}
for i in all_cols[1:]:
    annual_Q[i] = weight_isotopes(stream_d2h, i)

ppt_iso = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\monthly_data2.csv')
site_data = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\site_data.csv')
q_sum = pd.read_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\discharge_sum.csv')

def err_prop_add_sub(std_errs):
    std_errs_sq = []
    for i in std_errs:
        if not np.isnan(i):
            std_errs_sq.append(i **2)
    std_err_out = math.sqrt(sum(std_errs_sq))
    return std_err_out

def err_prop_mult_div(value, vars, errs):
    errs_sq_scaled = []
    for i in range(len(vars)):
        errs_sq_scaled.append((errs[i] / vars[i]) ** 2)
    std_err_out = abs(value) * math.sqrt(sum(errs_sq_scaled))
    return std_err_out

final_data = pd.DataFrame(data=[], columns=['site', 'representation of Ps in Q', 'SE of representation of Ps in Q',\
                                            'f_ET_from_Ps', 'f_ET_from_Ps_se', 'f_Ps_to_ET', 'f_Ps_to_ET_se', 'Ptot (mm)', 'Ptot_se', 'WS_area_km_sq', \
                                            'Qtot (mm)', 'Qtot_se', 'ratio_Ps_to_P', 'ratio_se', 'f_Q_from_Ps', \
                                            'f_Q_from_Ps_se', 'Pwdel (d2h)', 'Pwdel_se', 'Psdel (d2h)', 'Psdel_se', \
                                            'Qdel (d2h)', 'Qdel_se', 'ratio_ET_to_Ptot', 'meanT (C)', 'elevation (m)'])


for i in all_cols[1:]:
    Pw_del_wtd = []
    Pw_wts = []
    Pw_wts_se = []
    Pw_del_se = []
    Pw_del_wtd_se = []
    Ps_del_wtd = []
    Ps_wts = []
    Ps_wts_se = []
    Ps_del_se = []
    Ps_del_wtd_se = []
    for j in range(0, 4):
        index = list(ppt_iso['SITEID']).index(i)
        Pdel_wtd = ppt_iso.loc[index, 'wtd_d2h_' + months[j] + '_SUM']
        Pw_del_wtd.append(Pdel_wtd)
        Pwts = ppt_iso.loc[index, 'ppt_' + months[j] + '_SUM']
        Pw_wts.append(Pwts)
        Pwts_se = ppt_iso.loc[index, 'ppt_' + months[j] + '_SUM'] * 0.055  # monthly P se is 5.5%
        Pw_wts_se.append(Pwts_se)
        Pdel_se = math.sqrt(ppt_iso.loc[index, 'd2h_' + months[j] + '_mean_se_sq_MEAN'])
        Pw_del_se.append(Pdel_se)
    for j in range(10, 12):
        index = list(ppt_iso['SITEID']).index(i)
        Pdel_wtd = ppt_iso.loc[index, 'wtd_d2h_' + months[j] + '_SUM']
        Pw_del_wtd.append(Pdel_wtd)
        Pwts = ppt_iso.loc[index, 'ppt_' + months[j] + '_SUM']
        Pw_wts.append(Pwts)
        Pwts_se = ppt_iso.loc[index, 'ppt_' + months[j] + '_SUM'] * 0.055  # monthly P se is 5.5%
        Pw_wts_se.append(Pwts_se)
        Pdel_se = math.sqrt(ppt_iso.loc[index, 'd2h_' + months[j] + '_mean_se_sq_MEAN'])
        Pw_del_se.append(Pdel_se)
    for j in range(4, 10):
        index = list(ppt_iso['SITEID']).index(i)
        Pdel_wtd = ppt_iso.loc[index, 'wtd_d2h_' + months[j] + '_SUM']
        Ps_del_wtd.append(Pdel_wtd)
        Pwts = ppt_iso.loc[index, 'ppt_' + months[j] + '_SUM']
        Ps_wts.append(Pwts)
        Pwts_se = ppt_iso.loc[index, 'ppt_' + months[j] + '_SUM'] * 0.055  # monthly P se is 5.5%
        Ps_wts_se.append(Pwts_se)
        Pdel_se = math.sqrt(ppt_iso.loc[index, 'd2h_' + months[j] + '_mean_se_sq_MEAN'])
        Ps_del_se.append(Pdel_se)

    Psdel = sum(Ps_del_wtd) / sum(Ps_wts)
    Psdel_se = err_prop_mult_div(Psdel, [sum(Ps_del_wtd), sum(Ps_wts)], [err_prop_add_sub(Ps_del_se),
                                                                           err_prop_add_sub(Ps_wts_se)])
    Pwdel = sum(Pw_del_wtd) / sum(Pw_wts)
    Pwdel_se = err_prop_mult_div(Pwdel, [sum(Pw_del_wtd), sum(Pw_wts)], [err_prop_add_sub(Pw_del_se),
                                                                           err_prop_add_sub(Pw_wts_se)])

    index = list(site_data['SITEID']).index(i)
    elev_diff = site_data['DEM_10m'][index] - site_data['DEM_1km'][index]
    d2h_correction = elev_diff * -2.24 / 100
    area = site_data['area_km_sq'][index]
    meanT = site_data['Tmean_C'][index]
    elev = site_data['DEM_10m'][index]

    Psdel += d2h_correction
    Pwdel += d2h_correction

    Qdel, Qdel_se = annual_Q[i]

    s_day = 86400
    cubicm_l = 0.001
    sqm_sqkm = 1000000
    m_to_mm = 1000

    Ps = ppt_iso.loc[index, 'ppt_05_MEAN'] + ppt_iso.loc[index, 'ppt_06_MEAN'] + ppt_iso.loc[index, 'ppt_07_MEAN'] \
        + ppt_iso.loc[index, 'ppt_08_MEAN'] + ppt_iso.loc[index, 'ppt_09_MEAN'] + ppt_iso.loc[index, 'ppt_10_MEAN']
    Ps_se = Ps * 0.055
    Pw = ppt_iso.loc[index, 'ppt_01_MEAN'] + ppt_iso.loc[index, 'ppt_02_MEAN'] + ppt_iso.loc[index, 'ppt_03_MEAN'] \
        + ppt_iso.loc[index, 'ppt_04_MEAN'] + ppt_iso.loc[index, 'ppt_11_MEAN'] + ppt_iso.loc[index, 'ppt_12_MEAN']
    Pw_se = Ps * 0.055
    Ptot = Ps + Pw
    Ptot_se = Ptot * 0.055
    f_Q_from_Ps = (Qdel - Pwdel) / (Psdel - Pwdel)
    f_Q_from_Ps_se = err_prop_mult_div(f_Q_from_Ps, [Qdel - Pwdel, Psdel - Pwdel], [err_prop_add_sub([Qdel_se, Pwdel_se]),
                                                                            err_prop_add_sub([Psdel_se, Pwdel_se])])

    if i in ['CRAM', 'LIRO', 'PRLA', 'PRPO']:
        Qtot = np.nan
        Qtot_se = np.nan
    else:
        index = list(q_sum['site']).index(i)
        Qtot = q_sum['Qtot (mm)'][index]
        Qtot_se = q_sum['Qtot_se'][index]
    ratio_Ps_to_P = Ps / (Ps + Pw)
    ratio_se = err_prop_mult_div(ratio_Ps_to_P, [Ps, Ps + Pw], [Ps_se, err_prop_add_sub([Ps_se, Pw_se])])
    soi = f_Q_from_Ps / ratio_Ps_to_P
    soi_se = err_prop_mult_div(soi, [f_Q_from_Ps, ratio_Ps_to_P], [f_Q_from_Ps_se, ratio_se])

    f_ET_from_Ps = (Ps - Qtot * f_Q_from_Ps) / (Ptot - Qtot)
    f_ET_from_Ps_se = err_prop_mult_div(f_ET_from_Ps, [Ps - Qtot * f_Q_from_Ps, Ptot - Qtot],
        [err_prop_add_sub([Ps_se,
         err_prop_mult_div(Qtot * f_Q_from_Ps, [Qtot, f_Q_from_Ps], [Qtot_se, f_Q_from_Ps_se])]),
         err_prop_add_sub([Ptot_se, Qtot_se])])

    f_Ps_to_ET = f_ET_from_Ps * ((Ptot - Qtot) / Ps)
    f_Ps_to_ET_se = err_prop_mult_div(f_ET_from_Ps * ((Ptot - Qtot) / Ps), [f_ET_from_Ps, ((Ptot - Qtot) / Ps)],
                                      [f_ET_from_Ps_se, err_prop_mult_div((Ptot - Qtot) / Ps, [Ptot - Qtot, Ps],
                                                                          [err_prop_add_sub([Ptot_se, Qtot_se]), Ps_se])])

    final_data.loc[len(final_data.index)] = [i, soi, soi_se, f_ET_from_Ps, f_ET_from_Ps_se, f_Ps_to_ET, f_Ps_to_ET_se, Ptot, Ptot_se, area, Qtot, \
                                             Qtot_se, ratio_Ps_to_P, ratio_se, f_Q_from_Ps, f_Q_from_Ps_se, Pwdel, \
                                             Pwdel_se, Psdel, Psdel_se, Qdel, Qdel_se, (Ptot - Qtot) / Ptot, meanT, elev]


with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
    print(final_data)

#final_data.to_csv(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\code_output4.csv')

def calc_ci(t, s_err, n, x, x2, y2):
    ci = t * s_err * np.sqrt(1 / n + (x2 - np.mean(x)) ** 2 / np.sum((x - np.mean(x)) ** 2))
    ci_upp = y2 + ci
    ci_low = y2 - ci
    return ci_low, ci_upp

def ols(x_vals, y_vals, x_val_ci=None):
    x, y = zip(*sorted(zip(x_vals, y_vals)))
    X = sm.add_constant(x)
    model_ols = sm.OLS(y, X).fit()
    slope_pval = model_ols.pvalues[1]
    if slope_pval < 0.1:
        print(slope_pval)
        y_pred = model_ols.fittedvalues.tolist()
        df = model_ols.df_resid
        t_crit = abs(scipy.stats.t.ppf(q=0.025, df=df))
        y_resid_sq = []
        for i in range(len(y)):
            y_resid_sq.append((y[i] - y_pred[i]) ** 2)
        s_err = np.sqrt(np.sum(y_resid_sq) / df)
        y_fitted = [0] * len(x_val_ci)
        for i in range(len(x_val_ci)):
            y_fitted[i] = model_ols.params[1] * x_val_ci[i] + model_ols.params[0]
        ci_low, ci_upp = calc_ci(t_crit, s_err, len(x), x, x_val_ci, y_fitted)
        plt.fill_between(x_val_ci, ci_upp, ci_low, alpha=.4)
        plt.plot(x, y_pred, linewidth=2)
    plt.scatter(x, y, c='tab:blue')


ols(final_data['Ptot (mm)'], final_data['f_ET_from_Ps'], range(430, 2700, 25))
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitess)', 30)), fontsize=20)
plt.xlabel('Mean Annual Precipitation (mm)', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()


ols(final_data['ratio_ET_to_Ptot'], final_data['f_ET_from_Ps'], [0.27, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitess)', 30)), fontsize=20)
plt.xlabel('Ratio of ET to Precipitation (unitless)', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()

ols(final_data['elevation (m)'], final_data['f_ET_from_Ps'], [0.27, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitess)', 30)), fontsize=20)
plt.xlabel('Elevation (m)', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()

'''
plt.scatter(data=final_data, y='f_ET_from_Ps', x='Ptot (mm)')
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitess)', 30)), fontsize=20)
plt.xlabel('Mean Annual Precipitation (mm)', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()
'''
plt.scatter(data=final_data, y='f_ET_from_Ps', x='WS_area_km_sq')
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitless)', 30)), fontsize=20)
plt.xlabel('Log(Watershed Area (km$^{2}$))', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()
'''
plt.scatter(data=final_data, y='f_ET_from_Ps', x='meanT (C)')
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitless)', 30)), fontsize=20)
plt.xlabel('Mean Annual Temperature (\N{DEGREE SIGN}C)', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()


plt.scatter(data=final_data, y='f_ET_from_Ps', x='ratio_ET_to_Ptot')
plt.ylabel('\n'.join(wrap('Fraction of ET from Summer Precipitation (unitless)', 30)), fontsize=20)
plt.xlabel('Ratio of ET to Precipitation (unitless)', fontsize=20)
plt.tick_params(axis='both', labelsize=16)
plt.tight_layout()
plt.show()

'''

rank = [2, 20, 6, 10, 19, 11, 18, 12, 3, 5,  7, 22, 17,  4, 23,  9,  8, 14, 15, 21, 16, 13,  1]
final_data['rank'] = rank
for i in range(len(final_data)):
    print(final_data['site'][i], final_data['rank'][i], round(final_data['f_ET_from_Ps'][i], 1), round(final_data['ratio_Ps_to_P'][i], 1), round(final_data['WS_area_km_sq'][i], 1), round(final_data['Ptot (mm)'][i], 1), round(final_data['meanT (C)'][i], 1))


plt.figure(figsize=(5, 9))
#plt.axvline(0, c='black', linestyle='dashed')
#plt.axvline(1, c='black', linestyle='dashed')
plt.scatter(data=final_data, x='ratio_Ps_to_P', y='rank', facecolors='none', edgecolors='dodgerblue', s=40, marker='o', label='\n'.join(wrap('Ratio of Summer Precipitation to Total Precipitation', 20)))
plt.scatter(data=final_data, x='f_ET_from_Ps', y='rank', color='green', s=40, label='\n'.join(wrap('Fraction of ET from Summer Precipitation', 20)))
for i in range(len(final_data['f_ET_from_Ps'])):
    plt.plot((final_data['ratio_Ps_to_P'][i] - final_data['ratio_se'][i], final_data['ratio_Ps_to_P'][i] +
              final_data['ratio_se'][i]), (final_data['rank'][i], final_data['rank'][i]), color='dodgerblue',
             linewidth=1)
    plt.plot((final_data['f_ET_from_Ps'][i] - final_data['f_ET_from_Ps_se'][i], final_data['f_ET_from_Ps'][i] +
              final_data['f_ET_from_Ps_se'][i]), (final_data['rank'][i], final_data['rank'][i]), color='green', linewidth=3)

plt.margins(y=3)
plt.legend()
plt.tick_params(axis='both', labelsize=14)
plt.tight_layout()
#plt.savefig(r'C:\Users\User\Documents\UNR\LASTSEMESTER!!\Project2\f_ET_sites.svg', dpi=500)
plt.show()
