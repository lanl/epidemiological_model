"""
This script fits time-varying logistic growth parameters to the Total Mosquito population of Toronto.

"data" is obtained from the Mosquito PBM output

@author: Marina Mancuso (mmancuso@lanl.gov)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import leastsq

# note: I'm not sure where this data is on Darwin, but I can send over the CSV if needed- Marina
# import data
df = pd.read_csv("/Users/mbarnard/Documents/MosqPBM_out.csv")

# find column names
col_names = list(df.columns)
print(col_names)
    # col 1 = Unamed: 0
    # col 2 = times
    # col 3 = date
    # col 4=  location
    # col 5 = observed.mosq, i.e., from trap data
    # col 6 = modeled.mosq, i.e., active mosquitoes
    # col 7 = TotMosq, i.e., total mosquitoes


# See locations contained in data
df['location'].unique()

# Group by location
byloc = df.groupby('location')
print(byloc.groups)

# subset by location
Toronto = byloc.get_group('Toronto')

# sort the time series
Toronto_sort = Toronto.sort_values(by = 'date')    

# find length of time series
T_length = len(Toronto_sort)
print(T_length)

# find date range of time series
T_start = pd.to_datetime(Toronto_sort['date'].min())
print(T_start)

T_end = pd.to_datetime(Toronto_sort['date'].max())
print(T_end)

T_range = pd.date_range(start = T_start, end = T_end)
print(T_range)

# check if time series length equals length of mosquito observations
len(T_range) == T_length

# find mean of mosquitoes by location
mean_mosq = byloc['TotMosq'].agg([np.mean])
print(mean_mosq)

mean_mosq_array = np.array(mean_mosq)   # transform to np.array

Kv_mean_T = mean_mosq_array[2]
Kv_mean_BC = mean_mosq_array[0]
Kv_mean_C = mean_mosq_array[1]


# Convert mosquito time series to np.array and view mosquito time series
T_mosq = np.array(Toronto_sort['TotMosq'])
plt.plot(T_mosq)

# Check to make sure there are no NA values for mosquitoes or dates
for x in Toronto_sort['TotMosq'].isnull():
    if x == True:
        print('Mosquito has no value')
        
for x in Toronto_sort['date'].isnull():
    if x == True:
        print('Date has no value')

#%% Functions
        
def logistic(x, t, min_rv, max_rv, min_K, max_K):
    scale_rv = (max_rv - min_rv)/2
    base_rv = min_rv + scale_rv
    
    scale_K = (max_K - min_K)/2
    base_K = min_K + scale_K
 
    p = 365/np.pi
    
    r_v = base_rv + scale_rv*(np.sin(2*t/p - np.pi/2))
    K_v = base_K + scale_K*np.sin(2*t/p - np.pi/2)

    dxdt = (1- (x/K_v))*r_v*x 
    return dxdt

def cost(init_params, mosq): # pass in initial paramter values, mosquito data
    init_params = tuple(init_params)
    t = np.arange(0,len(mosq))
    x0 = list(mosq)[0]
    model = odeint(logistic, x0, t, args = init_params).flatten()
    res = model - mosq
    return res.flatten()

def fit_data(init_params, mosq): # pass in initial paramter values, mosquito data

    fitted_params = leastsq(cost, init_params, args = mosq)[0]
    print('fitted min r_v = ', fitted_params[0])
    print('fitted max r_v = ', fitted_params[1])
    print('fitted min K_v = ', fitted_params[2])
    print('fitted max K_v = ', fitted_params[3])
    
    x0 = list(mosq)[0]
    t = np.arange(0,len(mosq))
    fit_x = odeint(logistic, x0, t, args = tuple(fitted_params))

    plt.plot(t, mosq, 'r', linewidth = 3)
    plt.plot(t, fit_x, 'b-', linewidth = 3)
    plt.legend(['Data','Model'], loc = 'best')
    plt.xlabel('t (days)')
    plt.ylabel('(Susceptible) Mosquitoes')
    plt.show()
    
    return (fit_x, fitted_params)


#%% fit by year

season_start = [25, 370, 720, 1101, 1460, 1828, 2190, 2560, 2920, 3286, 3655, 4009, 4390, 4760] 
season_end = [200, 559, 925, 1274, 1650, 2033, 2390, 2748, 3122, 3469, 3843, 4238, 4535, 4871]

start_year = 2004

n = np.arange(0,len(season_start)-1)

def fit_subset(i):  # for running one at a time
    # mosquito season
    T_mosq_sub = T_mosq[season_start[n[i]]:season_end[n[i]]]

    # fit
    min_rv = -0.07
    max_rv = 0.07
    min_K = np.mean(T_mosq_sub) - 100
    max_K = np.mean(T_mosq_sub) + 100

    T_params = [min_rv, max_rv, min_K, max_K]
    T_params = np.array(T_params)

    print('Year', i, ':', i + start_year)
    print('-----------------------------')
  
    T_fit = fit_data(T_params, T_mosq_sub)

    # find r_b, r_s, K_b, K_s
    rv_min = T_fit[1][0]      # rv min
    rv_max = T_fit[1][1]      # rv max
    Kv_min = T_fit[1][2]      # Kv min
    Kv_max = T_fit[1][3]      # Kv max

    r_s = (rv_max - rv_min)/2
    r_b = rv_min + r_s
    K_s = (Kv_max - Kv_min)/2
    K_b = Kv_min + K_s

    print('r_b = ', r_b)
    print('r_s = ', r_s)
    print('K_b = ', K_b)
    print('K_s = ', K_s)
    
    # Mosq season begin 
    td_season_start = pd.Timedelta(days = season_start[i])
    print('Mosq season begins: ', T_start + td_season_start)

    # Mosq season end
    td_season_end = pd.Timedelta(days = season_end[i])
    print('Mosq season ends: ', T_start + td_season_end)

    # Mosq season length
    print('Mosq season length = ', td_season_end - td_season_start)
    return T_fit[0]

def plot_together():
    
    # determine each year's fitting interval
    fit_length = np.array(season_end) - np.array(season_start)
    
    # create array to collect all years' fits together
    max_length = max(fit_length)
    fit_array = np.zeros((len(n), max_length))

    # fit each year
    for i in n:
        fit_array[i][0:fit_length[i]] = (fit_subset(i)).T
 
    # plot mosquito PBM
    plt.plot(T_mosq, 'r.')  
    
    # plot each years' fit togther
    for i in n:
        plt.plot(np.arange(season_start[i],season_end[i]), fit_array[i][0:fit_length[i]], 'b', linewidth = 3)

#plot_together()
