import pandas as pd
import numpy as np
import math

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

data = pd.read_csv(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\code_input_f_ET_f_Ps.csv', encoding= 'unicode_escape')

final_data = pd.DataFrame(data=[], columns=['site', 'f_ET_from_Ps', 'f_ET_from_Ps_se', 'f_Ps_to_ET', 'f_Ps_to_ET_se'])


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
    final_data.loc[len(final_data.index)] = [site, f[1][0], f_se[1][0], eta[1][0], eta_se[1][0]]
    return f, f_se, eta, eta_se

for i in range(len(data['Psdel'])):
    end_member_mixing(data['site'][i], [data['Psdel'][i], data['Pwdel'][i]], data['Qdel'][i],
                      [data['Psdel_se'][i], data['Pwdel_se'][i]], data['Qdel_se'][i],
                      [data['Ps'][i], data['Pw'][i]], data['Qtot'][i],
                      [data['Ps_se'][i], data['Pw_se'][i]], data['Qtot_se'][i])

# Add columns for additional variables
final_data['Ptot'] = data['Ptot']
final_data['WS_area_km_sq'] = data['WS_area_km_sq']
final_data['ratio_Ps_to_P'] = data['ratio_Ps_to_P']
final_data['ratio_ET_to_P'] = data['ratio_ET_to_Ptot']
final_data['meanT_C'] = data['meanT_C']
final_data['elevation_m'] = data['elevation (m)']

# going to manually add EVI metrics and GEDI histogram metrics to the output
final_data.to_csv(r'C:\Users\User\Documents\UNR\Summer 2023\IsotopeMapping\code_output_f_ET_f_Ps.csv')