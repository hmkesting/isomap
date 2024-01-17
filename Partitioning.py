import pickle
import pandas as pd
import math
import numpy as np
from operator import add
import seaborn as sn
import matplotlib.pyplot as plt

# removed BIGC: SJER, TECR: SJER, BLDE: YELL, and REDB: REDB because they were missing several months of precipitation isotope data
# write this into precipitation script later
                
precip_sites = {'ARIK': 'ARIK', 'BLUE': 'BLUE', 'BLWA': 'DELA', 'COMO': 'NIWO', 
                'FLNT': 'JERC', 'HOPB': 'HARV', 'KING': 'KONZ', 'LECO': 'GRSM', 'LEWI': 'BLAN', 'MART': 'WREF', 
                'MAYF': 'TALL', 'MCDI': 'KONZ', 'POSE': 'SCBI', 'PRIN': 'PRIN', 'SYCA': 'SYCA',
                'TOMB': 'LENO', 'WALK': 'ORNL', 'WLOU': 'NIWO'}
sites = list(precip_sites.keys())

other_sites = ['BIGC', 'TECR', 'BLDE', 'REDB', 'MCRA']

# Elevations copied from NEON site descriptions
precip_elev_m = {'ARIK': 1194.89, 'SJER': 391.66, 'YELL': 2133.54, 'BLUE': 293.43, 'DELA': 67.55, 'NIWO': 3488.10, 'JERC': 89.21, 
                  'HARV': 340.12, 'KONZ': 416.61, 'GRSM': 621.42, 'BLAN': 191.42, 'WREF': 357.44, 'TALL': 202.02, 
                  'SCBI': 354.67, 'PRIN': 262.08, 'REDB': 1717.96, 'SYCA': 655.01, 'LENO': 60.45, 'ORNL': 245.53}

# Elevation corrections to precipitation values
correction_factor = -2.24 / 100 
with open('NEON_Pdel_s.pickle', 'rb') as f:
    NEON_Pdel_s = pickle.load(f)
    
with open('NEON_Pdel_w.pickle', 'rb') as f:
    NEON_Pdel_w = pickle.load(f)
    
wtd_elev = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\wtd_elev.csv')
NEON_pdel = {}
NEON_pdel_se = {}
for i in sites:
    idx = list(wtd_elev['site']).index(i)
    summer_elev = wtd_elev['summer_wtd_elev'][idx]
    winter_elev = wtd_elev['winter_wtd_elev'][idx]
    idx = precip_sites[i]
    NEON_pdel[i] = (list(NEON_Pdel_s[idx])[0] + (summer_elev - precip_elev_m[idx]) * correction_factor, \
                    list(NEON_Pdel_w[idx])[0] + (winter_elev - precip_elev_m[idx]) * correction_factor)                   
    NEON_pdel_se[i] = (list(NEON_Pdel_s[idx])[1], list(NEON_Pdel_w[idx])[1])

precip_grid = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\precip_d2h_grid.csv')
for i in range(len(precip_grid['site'])):
    precip_grid.loc[i, 'p_d2h_summer'] += (precip_grid['elev_30m'][i] - precip_grid['elev_1km'][i]) * correction_factor
    precip_grid.loc[i, 'p_d2h_winter'] += (precip_grid['elev_30m'][i] - precip_grid['elev_1km'][i]) * correction_factor

columns = ['site', 'Pdels_grid', 'Pdels_neon', 'Pdels_avg', 'Pdels_avg_se', \
           'Pdelw_grid', 'Pdelw_neon', 'Pdelw_avg', 'Pdelw_avg_se', 'Ps', 'Ps_se', 'Pw', 'Pw_se', \
            'Qdel', 'Qdel_se', 'Q', 'Q_se']

with open('NEON_Qtot.pickle', 'rb') as f:
    Qtot = pickle.load(f)
    
with open('NEON_Qdel.pickle', 'rb') as f:
    Qdel = pickle.load(f)
   
partitioning_fractions = pd.DataFrame(data=[], columns=['site', 'f_ET_from_Ps', 'f_ET_from_Ps_se', 'f_Ps_to_ET', 'f_Ps_to_ET_se'])

def end_member_mixing(site, pdel_bar, qdel_bar, pdel_se, qdel_se, ptot, qtot, ptot_se, qtot_se):
    et = sum(ptot) - qtot
    f = [[0, 0], [0, 0]]
    f_se = [[0, 0], [0, 0]]

    denom = pdel_bar[0] - pdel_bar[1]

    f[0][0] = (qdel_bar - pdel_bar[1]) / denom
    f[0][1] = 1 - f[0][0]
    f_se[0][0] = math.sqrt((qdel_se / denom) ** 2 + (pdel_se[0] * (-f[0][0] / denom)) ** 2 + (
                    pdel_se[1] * (pdel_bar[0] - qdel_bar) / denom ** 2) ** 2)
    f_se[0][1] = f_se[0][0]

    f[1][0] = (ptot[0] - qtot * (qdel_bar - pdel_bar[1]) / denom) / (ptot[0] + ptot[1] - qtot)
    f_se[1][0] = math.sqrt((ptot_se[0] * (1 - f[1][0]) / et) ** 2
                           + (ptot_se[1] * f[1][0] / et) ** 2
                           + (qtot_se * (f[1][0] - f[0][0]) / et) ** 2
                           + (qdel_se * qtot / (pdel_bar[0] - pdel_bar[1]) / et) ** 2
                           + (pdel_se[0] * qtot * f[0][0] / (pdel_bar[0] - pdel_bar[1]) / et) ** 2
                           + (pdel_se[1] * qtot * f[0][1] / (pdel_bar[1] - pdel_bar[0]) / et) ** 2)
    f[1][1] = 1 - f[1][0]
    f_se[1][1] = f_se[1][0]

    eta = [[0, 0], [0, 0]]
    eta_se = [[0, 0], [0, 0]]

    for i in range(2):
        eta[0][i] = f[0][i] * qtot / ptot[i]
        eta_se[0][i] = abs(eta[0][i]) * math.sqrt((f_se[0][i] / f[0][i]) ** 2
                                                  + (qtot_se / qtot) ** 2 + (ptot_se[i] / ptot[i]) ** 2)

    for i in range(2):
        eta[1][i] = 1 - eta[0][i]
        eta_se[1][i] = eta_se[0][i]
    partitioning_fractions.loc[len(partitioning_fractions.index)] = [site, f[1][0], f_se[1][0], eta[1][0], eta_se[1][0]]
    return f, f_se, eta, eta_se


Pdels_neon = []
Pdels_grid = []
Pdelw_neon = []
Pdelw_grid = []
Q_del = []
Ps = []
Pw = []
Q = []
for i in sites:
    idx = list(precip_grid['site']).index(i)
    Pdel = [(NEON_pdel[i][0] + precip_grid.loc[idx, 'p_d2h_summer'])/2, \
            (NEON_pdel[i][1] + precip_grid.loc[idx, 'p_d2h_winter'])/2]
    Pdel_se = [math.sqrt(NEON_pdel_se[i][0] ** 2 + precip_grid.loc[idx, 'p_d2h_summer_se'] ** 2)/2, \
              math.sqrt(NEON_pdel_se[i][1] ** 2 + precip_grid.loc[idx, 'p_d2h_winter_se'] ** 2)/2]
    Ptot = [precip_grid['p_summer'][idx], precip_grid['p_winter'][idx]]
    Ptot_se = [precip_grid['p_summer'][idx] * 0.055, precip_grid['p_winter'][idx] * 0.055]

    end_member_mixing(i, Pdel, Qdel[i][0], Pdel_se, Qdel[i][1], Ptot, Qtot[i][0], Ptot_se, Qtot[i][1])

    Pdels_neon.append(NEON_pdel[i][0])
    Pdels_grid.append(precip_grid.loc[precip_grid['site'] == i, 'p_d2h_summer'].iloc[0])
    Pdelw_neon.append(NEON_pdel[i][1])
    Pdelw_grid.append(precip_grid.loc[precip_grid['site'] == i, 'p_d2h_winter'].iloc[0])
    Q_del.append(Qdel[i][0])
    Ps.append(Ptot[0])
    Pw.append(Ptot[1])
    Q.append(Qtot[i][0])

for i in other_sites:
    idx = list(precip_grid['site']).index(i)
    Pdel = [precip_grid.loc[idx, 'p_d2h_summer'], precip_grid.loc[idx, 'p_d2h_winter']]
    Pdel_se = [precip_grid.loc[idx, 'p_d2h_summer_se'], precip_grid.loc[idx, 'p_d2h_winter_se']]
    Ptot = pdel = [precip_grid.loc[idx, 'p_summer'], precip_grid.loc[idx, 'p_winter']]
    Ptot_se = [precip_grid.loc[idx, 'p_summer'] * 0.055, precip_grid.loc[idx, 'p_winter'] * 0.055]
    
    end_member_mixing(i, Pdel, Qdel[i][0], Pdel_se, Qdel[i][1], Ptot, Qtot[i][0], Ptot_se, Qtot[i][1])
    
    Pdels_neon.append(np.nan)
    Pdels_grid.append(Pdel[0])
    Pdelw_neon.append(np.nan)
    Pdelw_grid.append(Pdel[1])
    Q_del.append(Qdel[i][0])
    Ps.append(Ptot[0])
    Pw.append(Ptot[1])
    Q.append(Qtot[i][0])

partitioning_fractions['pdels_neon'] = Pdels_neon
partitioning_fractions['pdels_grid'] = Pdels_grid
partitioning_fractions['pdelw_neon'] = Pdelw_neon
partitioning_fractions['pdelw_grid'] = Pdelw_grid
partitioning_fractions['qdel'] = Q_del
partitioning_fractions['ps'] = Ps
partitioning_fractions['pw'] = Pw
partitioning_fractions['ptot'] = list(map(add, Ps, Pw))
partitioning_fractions['q'] = Q
partitioning_fractions['ratio_ps_to_ptot'] = np.divide(Ps, list(partitioning_fractions['ptot']))
ratio_q_to_ptot = np.divide(Q, list(partitioning_fractions['ptot']))
partitioning_fractions['ratio_et_to_ptot'] = [1 - x for x in ratio_q_to_ptot]

site_area_sq_km = {'ARIK': 2632.3, 'BIGC': 10.9009, 'BLDE': 37.839, 
                   'BLUE': 322.15, 'BLWA': 16171.6, 'COMO': 3.56629, 
                   'FLNT': 14986.7, 'HOPB': 11.882, 'KING': 13.0413, 
                   'LECO': 9.12099, 'LEWI': 11.9086, 'MART': 6.34773, 
                   'MAYF': 14.3621, 'MCDI': 22.5796, 'MCRA': 3.93222, 
                   'POSE': 2.02731, 'PRIN': 48.9114, 'REDB': 16.7148, 
                   'SYCA': 280.475, 'TECR': 2.97479, 'TOMB': 47109.9, 
                   'WALK': 1.08996, 'WLOU': 4.90348}

# Add remote sensing variables
min_evi = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\min_EVI.csv')
amp_evi = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\amp_EVI.csv')
area_evi = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\area_EVI.csv')
greenup = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\greenup_EVI.csv')
dormancy = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\dormancy_EVI.csv')
mean_gedi = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\means_GEDI.csv')
stdev_gedi = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\stdev_GEDI.csv')
hist_gedi = pd.read_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\hist_GEDI.csv')

partitioning_fractions['temperature'] = [wtd_elev.loc[wtd_elev['site'] == x, 'temperature'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['elevation'] = [precip_grid.loc[precip_grid['site'] == x, 'elev_30m'].iloc[0] for x in partitioning_fractions['site']]

hist_gedi['pct_abv_5m'] = [sum(map(float, x[1:-1].split(',')[1:]))/sum(map(float, x[1:-1].split(','))) for x in hist_gedi['means']]
hist_gedi['pct_abv_10m'] = [sum(map(float, x[1:-1].split(',')[2:]))/sum(map(float, x[1:-1].split(','))) for x in hist_gedi['means']]
hist_gedi['pct_abv_20m'] = [sum(map(float, x[1:-1].split(',')[4:]))/sum(map(float, x[1:-1].split(','))) for x in hist_gedi['means']]

partitioning_fractions['area_sq_km'] = [site_area_sq_km[x] for x in partitioning_fractions['site']]
partitioning_fractions['min_evi'] = [min_evi.loc[min_evi['SiteID'] == x, 'mean'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['amp_evi'] = [amp_evi.loc[amp_evi['SiteID'] == x, 'mean'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['area_evi'] = [area_evi.loc[area_evi['SiteID'] == x, 'mean'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['len_grow_season'] = [dormancy.loc[dormancy['SiteID'] == x, 'mean'].iloc[0] - \
                                             greenup.loc[greenup['SiteID'] == x, 'mean'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['ratio_evi_area_to_min'] = np.divide(list(partitioning_fractions['area_evi']), list(partitioning_fractions['min_evi']))

    
partitioning_fractions['mean_gedi'] = [mean_gedi.loc[mean_gedi['SiteID'] == x, 'mean'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['stdev_gedi'] = [stdev_gedi.loc[stdev_gedi['SiteID'] == x, 'stdDev'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['gedi_abv_5m'] = [hist_gedi.loc[hist_gedi['SiteID'] == x, 'pct_abv_5m'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['gedi_abv_10m'] = [hist_gedi.loc[hist_gedi['SiteID'] == x, 'pct_abv_10m'].iloc[0] for x in partitioning_fractions['site']]
partitioning_fractions['gedi_abv_20m'] = [hist_gedi.loc[hist_gedi['SiteID'] == x, 'pct_abv_20m'].iloc[0] for x in partitioning_fractions['site']]

partitioning_fractions.to_csv(r'C:\Users\User\Documents\UNR\NEON Project\Data\Compiled\partitioning.csv.')


plot_df = partitioning_fractions.drop(labels=['site', 'f_ET_from_Ps_se', 'f_Ps_to_ET_se', 'pdels_neon', 'pdels_grid', 'pdelw_neon', 'pdelw_grid', \
                                              'qdel', 'ps', 'pw', 'q'],  axis=1)

plot_df = plot_df.rename(columns={'f_ET_from_Ps': 'Fraction of ET from Ps', 'f_Ps_to_ET': 'Fraction of Ps to ET',
           'ptot': 'Mean Annual Precipitation', 'area_sq_km': 'Watershed Area',
           'ratio_ps_to_ptot': 'Ratio of Ps to Mean Annual Precipitation',
           'ratio_et_to_ptot': 'Ratio of ET to Mean Annual Precipitation', 'temperature': 'Mean Annual Temperature',
           'elevation': 'Elevation', 'gedi_abv_5m': 'Percent of Canopy Above 5m',
           'gedi_abv_10m': 'Percent of Canopy Above 10m', 'gedi_abv_20m': 'Percent of Canopy Above 20m',
           'mean_gedi': 'Mean Canopy Height', 'stdev_gedi': 'Standard Deviation of Canopy Height',
           'min_evi': 'Minimum EVI', 'amp_evi': 'EVI Amplitude', 'area_evi': 'EVI Area Under the Curve',
           'len_grow_season': 'Growing Season Length', 'ratio_evi_area_to_min': 'Ratio of EVI Area to EVI Minimum'})

plt.figure(figsize=(16, 16))
sn.set(font_scale=1.6)
sn.heatmap(plot_df.corr(), cmap='seismic', annot=True, fmt='.2f', annot_kws={"size":12})
plt.subplots_adjust(bottom=0.32, left=0.32)
plt.savefig(r'C:\Users\User\Documents\UNR\NEON Project\corrplot.svg', dpi=500)
plt.show()


