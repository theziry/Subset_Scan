import numpy as np
import LoadData as ld
import pandas as pd
import LuminosityOptimization as lo
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import time as t
#from lmfit import Model
from scipy import integrate
import scipy.optimize
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


#_____________________________________________
plt.close("all")

#Font setting
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica",
  "font.size": 12
})

#defining the start time of the program
start=t.time()

#Selecting Current Year
year=16

#plotting
plot=True

#model parameters 
if year==16:
   n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2016()
elif year==17:
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2017()
elif year==18:
    n_i, k_b, B_s, E_s, B_r, G_r, S_int, n_c, N_i, T_hc, T_ph, S_z, S_s, Fe, f_rev, Xi, Eps= lo.Parameters2018()


#loading fill number
FillNumber16, FillNumber17, FillNumber18 = ld.FillNumber()


if year==16:
    FillNumber=FillNumber16
elif year==17:
    FillNumber=FillNumber17
    FillNumber_Prev=FillNumber16
    previous_year=16
elif year==18:
    FillNumber=FillNumber18
    FillNumber_Prev=FillNumber17
    previous_year=17
#print(FillNumber)

# Fill to delete:
skip16=[5111,5112,5116,5117,5257,5258,5264,5406,5427,5439,5451]
skip17=[5837,5840,5856,6019,6055,6060,6082,6084,6089,6116,6152,6156,6158,6160,6167,6168,6169,6170,6193,6258,6268,6271,6272,6275,6283,6285,6287,6288,6291,6298,6303,6304,6305,6308,6311,6312,6314,6317,6324,6325,6337,6346,6356,6360,6362,6364,6371]
skip18=[6640,6654,6659,6688,6690,6693,6694,6696,6700,6706,6729,6741,6747,6752,6782,6778,6882,6890,6892,6901,6929,6939,6940,6961,7026,7033,7036,7039,7056,7069,7105,7118,7120,7135,7145,7217,7245,7259,7271,7324]
skip = skip16 + skip17 + skip18
for fill in skip:
    FillNumber=np.delete(FillNumber, (np.where(FillNumber==(fill))[0]))
for i in range(len(FillNumber)):
    if FillNumber[i] in skip:
        continue

with open('20{}/best_param/20{}_1param_last_best_FitCoeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close()
    
with open('20{}/best_param/20{}_1param_first_best_FitCoeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close()

with open('20{}/best_param/20{}_2param_last_best_FitCoeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close()
    
with open('20{}/best_param/20{}_2param_first_best_FitCoeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close()
    
with open('20{}/best_param/20{}_3param_Coeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close() 

with open('20{}/best_param/20{}_delta_2param_Coeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close() 

with open('20{}/best_param/20{}_delta_1param_Coeff.txt'.format(str(year),str(year)), 'w') as f:
    f.write('')
    f.close() 


#defining the double-exponential decay fit function
def fit(x, a, b, c, d):
    return (a*np.exp((-b)*x))+(c*np.exp((-d)*x))


### $\chi^2$ function
def chi_2(original, estimate):
    """Chi squared function
    Parameters
    ----------
    original : ndarray
        ideal values
    estimate : ndarray
        estimated values
    Returns
    -------
    float
        chi squared value
    """    
    return np.sum(np.power(estimate - original, 2) / original)


def L_model_3param(x, eps_ni, B, k):
        x = np.where(x>1, x, 1.1 ) 
        np.nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None) #if np.isnan(x): x = 0.0 
        #D = B / (np.log(x) ** k)   
        D = B * np.power(k / (2 * np.exp(1)), k) / (np.power(np.log(x), k))
        L = (1/(1 + eps_ni * (x-1))**2) - ((1 + D**2)*np.exp(-(D**2))) * ((2 - (1 + D**2)*np.exp(-(D**2)/2)))
        return L 


#defining the DA MODEL-4 fit function
def L_model_2param(x, eps_ni, B):
        x = np.where(x>1, x, 1.1)
        np.nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None) 
        #D = B / (np.log(x) ** k)
        D = B * np.power(k / (2 * np.exp(1)), k) / (np.power(np.log(x), k))
        L = (1/(1 + eps_ni * (x-1))**2) - ((1 + D**2)*np.exp(-(D**2))) * ((2 - (1 + D**2)*np.exp(-(D**2)/2)))
        return L

def L_model_1param(x, eps_ni):
        x = np.where(x>1, x, 2 )
        np.nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None) 
        #D = B / (np.log(x) ** k)
        D = B * np.power(k / (2 * np.exp(1)), k) / (np.power(np.log(x), k))
        L = (1/(1 + eps_ni * (x-1))**2) - ((1 + D**2)*np.exp(-(D**2))) * ((2 - (1 + D**2)*np.exp(-(D**2)/2)))
        return L

def Cut_Fit(year, text):
    """Function that performs the necessary cut on the current fill

    Args:
        year (int): current year
        text (str): current fill

    Returns:
        L_fit: cutted data
        T_fit_real: times in second for the cutted data fit
        Y: Luminosity evolution form the fit
        a: fitting parameter
        b: fitting parameter
        c: fitting parameter
        d: fitting parameter
        chi: reduced chi square of the fit
        L_evol: raw luminosity data
        Times: raw Unix time
    """
    year=str(year)
    f=open('ATLAS/ATLAS_fill_20{}/{}_lumi_ATLAS.txt'.format(year, text),"r")
    lines=f.readlines()
    L_evolx=[]
    times=[]
    for x in lines:
        times.append(int(x.split(' ')[0]))  
        L_evolx.append(float(x.split(' ')[2]))
        
    f.close()
    Times = np.array(times)
    L_evol = np.array(L_evolx)

    #deleting the null values of the luminosity
    zero=np.where(L_evol<100)
    L_zero=np.delete(L_evol, zero)
    T_zero=np.delete(Times, zero)
        
    #check for enough points
    if len(L_zero)<10:
        zero=np.where(L_evol<5)
        L_zero=np.delete(L_evol, zero)
        T_zero=np.delete(Times, zero)

    #defining the derivative 
    dy = np.zeros(L_zero.shape)
    dy[0:-1] = np.diff(L_zero)/np.diff(T_zero)


    #start to slim down the fit interval       
    L_tofit=[]
    T_tofit=[]
    for idx in range(len(L_zero)):
        #cancelling too strong derivative points
        if dy[idx]<0 and dy[idx]>-1.5:
            L_tofit.append(L_zero[idx])
            T_tofit.append(T_zero[idx])
        if dy[idx]>0 or dy[idx]<-1.5:
            continue     
        
    #evaluating the differences between two subsequent points
    diff=np.diff(L_tofit)
        
    #deleting the discrepancies
    thr=np.max(abs(diff))*0.05
    idx_diff= np.where(abs(diff)>thr)[0]+1
        
    #new slim down of data
    L_tofit2=np.delete(L_tofit, idx_diff)
    T_tofit2=np.delete(T_tofit, idx_diff)
        
    #check for enough points
    if len(L_tofit2) < 30:
        L_tofit2=L_tofit
        T_tofit2=T_tofit
        
    L_fit=L_tofit2
    T_fit=T_tofit2     

    L_fit=np.array(L_fit)
    T_fit=np.array(T_fit)

    #transforming the times from unix in seconds
    T_fit_real=T_fit-np.amin(T_fit)

    #Normlise L_fit 
    L_min = min(L_fit)
    L_max = max(L_fit)

    L_norm = L_fit/L_max 

    L_ma = np.array(L_max)
    L_Tau = np.array(L_norm)

    return T_fit_real,L_fit, L_ma, L_Tau, L_evol, Times


for i in range(len(FillNumber)):

    text = str(int(FillNumber[i])) #number of current fill
    T_fit_real,L_fit, L_ma, L_Tau, L_evol, Times = Cut_Fit(year, text)

    Turn=[] 
    Turn=np.array(Turn)
    for el in T_fit_real:
       tau=(f_rev*el+1)
       Turn=np.append(Turn, tau) 

    if year==16:
        p0 = (4.5e-10, 155, 0.9)
    elif year==17:
        p0 = (7.5e-10, 155, 0.9)
    elif year==18:
        p0 = (9.5e-10, 155, 0.9)
        
    #p0 = (2.5e-10, 155, 0.9)

    params, pcov = curve_fit(L_model_3param, Turn,L_Tau,p0,bounds=([2e-10, 50., 0.7], [2e-9, 200., 2]),maxfev=5000)
    eps_all, B_all, k_all = params


    # The fit output
    L = L_model_3param(Turn, eps_all, B_all, k_all)

    # RESIDUALS 
    AbsErr = L_Tau - L
    #print('absErr', AbsErr)
    SE = np.square(AbsErr) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE

    # calculate the R-squared value
    r_squared = r2_score(L_Tau, L)

    # number of observations
    n = len(L_Tau)
    # number of degree of freedom
    p = len(params) - 1

    # calculate the adjusted R-squared value
    R_adj = 1 - (1 - r_squared) * (n - 1) / (n - p - 1)

    # calculate the chi_square value
    chi_square = sum(((L_Tau - L)**2/L_Tau))
    # calculate the Reduced_chi_square
    R_chi = chi_square/(len(L_Tau) - 3)
    
    fig, ax = plt.subplots() 
    ax.plot(T_fit_real/3600, L_Tau, "b.", markersize=4, label='Data')
    ax.plot(T_fit_real/3600, L, "r-", markersize=4, label='3 param $Model$-$2$')
    ax.set_title('FillNumber: {}'.format(text))
    ax.set_xlabel(' Time [h]')
    ax.set_ylabel('${L}/{L_i}$' ) 
    plt.legend(loc='best')
    plt.savefig('20{}/best_mod/{}_3param_Fit.pdf'.format(str(year),text))

    with open('20{}/best_param/20{}_3param_Coeff.txt'.format(str(year),str(year)), 'a') as f:
        f.write(text)
        f.write(' ')
        f.write(str(eps_all))
        f.write(' ')
        f.write(str(B_all))
        f.write(' ')
        f.write(str(k_all))
        f.write(' ')
        f.write(str(R_adj))
        f.write(' ')
        f.write(str(chi_square))
        f.write('\n')
       
       


    with open('20{}/Coeff/{}_2param_last_FitCoeff.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 
        
    with open('20{}/Coeff/{}_2param_first_FitCoeff.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/Coeff/{}_1param_last_FitCoeff.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 
        
    with open('20{}/Coeff/{}_1param_first_FitCoeff.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 
        
        
    with open('20{}/data/{}_1param_Fitdata.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/data/{}_2param_Fitdata.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/data/{}_3param_Fitdata.txt'.format((str(year)),text), 'w') as f:
        f.write('')
        f.close() 
         

    p0 = (2.5e-10, 155)
    
    subsets = []
    subsets1 = []
    subsets2 = []
    subsets3 = []
    
   

    k_min = 0.1
    k_max = 2
    dk = 0.01
    #k_values = np.arange(k_min, k_max+dk, dk)
    k_values = np.linspace(k_min, k_max, int((k_max-k_min)/dk)+1)

    #k_arr = np.linspace(k_min, k_max, int((k_max-k_min)/dk)+1)
    L_Tau_fit = np.zeros((len(k_values), len(Turn)))
    

    eps_ni = np.zeros((len(k_values)))
    B = np.zeros((len(k_values)))
    R_squared = np.zeros((len(k_values)))
    chi_squared= np.zeros((len(k_values)))
    
    
    min_chi_squared = np.inf
    min_chi_squared_index = None
    
    for i, k in enumerate(k_values):
        #p0 = (2.5e-10, 155)
        
                
        popt, pcov = curve_fit(L_model_2param, Turn, L_Tau, p0=p0, bounds=([2e-10, 50], [2.e-9, 200.]), maxfev=5000)
        L_Tau_fit[i] = L_model_2param(Turn, popt[0], popt[1])
        eps_ni[i] = popt[0]
        B[i] = popt[1]

        residuals = L_Tau - L_model_2param(Turn, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((L_Tau - np.mean(L_Tau))**2)

        subsets.append((T_fit_real, L_Tau_fit[i]))     
   
        if ss_tot == 0:
            r_2 = 1
        else:
            r_2 = 1 - (ss_res / ss_tot)     

        # number of observations
        n = len(L_Tau)
        # number of predictors (not including the constant term)0.008*
        p = len(popt) - 1
        # calculate the adjusted R-squared value
        R_squared[i] = 1 - (1 - r_2) * (n - 1) / (n - p - 1)
        chi_squared[i] = np.sum(((residuals)**2)/L_Tau)
        
        
        print('FillNumber', text,'R_squared', R_squared)
        print('popt is',popt) 
        print('iteration i is',i) 
        #p0 = popt
                
        if chi_squared[i] < min_chi_squared:
           min_chi_squared = chi_squared[i]
           min_chi_squared_index = i
           
        print('min_chi_squared is',min_chi_squared) 
        
        
    best_k0 = k_values[min_chi_squared_index]
    best_B0 = B[min_chi_squared_index]
    best_eps_ni0 = eps_ni[min_chi_squared_index]
    best_R_squared0 = R_squared[min_chi_squared_index]
    best_chi_squared0 = chi_squared[min_chi_squared_index]
    k = best_k0
    best_L_Tau_fit0 = L_model_2param(Turn, best_eps_ni0, best_B0)

    print('i ---- best_k0 -- best_kB0 -- best_eps_n0 is',i, best_k0, best_B0, best_eps_ni0 )
    print('min_chi_squared and min_chi_squared_index is',min_chi_squared_index, min_chi_squared)

    fig, ax = plt.subplots() 
    ax.plot(T_fit_real/3600, L_Tau, "b.", markersize=4, label='Data')
    ax.plot(T_fit_real/3600, best_L_Tau_fit0, "r-", markersize=4, label='2 param full interval')
    ax.set_title('FillNumber: {}'.format(text))
    ax.set_xlabel(' Time [h]')
    ax.set_ylabel('${L}/{L_i}$' ) 
    plt.legend(loc='best')
    plt.savefig('20{}/best_mod/first_check/{}_2param_first_interval_Fit.pdf'.format(str(year),text))
        


    with open('20{}/best_param/20{}_2param_first_best_FitCoeff.txt'.format(str(year),str(year)), 'a') as f:
        f.write(text)
        f.write(' ')
        f.write(str(best_eps_ni0))
        f.write(' ')
        f.write(str(best_B0))
        f.write(' ')
        f.write(str(best_k0)) 
        f.write(' ')
        f.write(str(best_R_squared0))
        f.write(' ')
        f.write(str(best_chi_squared0))
        f.write('\n')      

        
    k_min_last = best_k0 - 0.009
    k_max_last = best_k0 + 0.009
    dk_last = 0.001
    #k_values_last = np.arange(k_min_last, k_max_last+dk_last, dk_last)
    k_values_last = np.linspace(k_min_last, k_max_last, int((k_max_last-k_min_last)/dk_last)+1)
    
    L_Tau_fit_last = np.zeros((len(k_values_last), len(Turn)))

    eps_ni_last = np.zeros((len(k_values_last)))
    B_last = np.zeros((len(k_values_last)))
    R_squared_last = np.zeros((len(k_values_last)))
    chi_squared_last = np.zeros((len(k_values_last)))

    min_chi_squared_last = np.inf
    min_chi_squared_index_last = None

    p0 = (best_eps_ni0, best_B0)
    
        
    for i, k in enumerate(k_values_last): 
        #p0 = (best_eps_ni0, best_B0) 
        

        popt2, pcov2 = curve_fit(L_model_2param, Turn, L_Tau, p0=p0, bounds=([2e-10, 50], [2.e-9, 200.]), maxfev=5000)
        L_Tau_fit_last[i] = L_model_2param(Turn, popt2[0], popt2[1])
        eps_ni_last[i] = popt2[0]
        B_last[i] = popt2[1]

        residuals2 = L_Tau - L_model_2param(Turn, *popt2)
        ss_res2 = np.sum(residuals2**2)
        ss_tot2 = np.sum((L_Tau - np.mean(L_Tau))**2)

        subsets2.append((T_fit_real, L_Tau_fit_last[i]))     

        if ss_tot2 == 0:
            r2 = 1
        else:
            r2 = 1 - (ss_res2 / ss_tot2)     

        R_squared_last[i] = 1 - (1 - r2) * (n - 1) / (n - p - 1)
        chi_squared_last[i] = np.sum(((residuals2)**2)/L_Tau)

        print('popt_last is',popt2) 
        print('iteration i is',i)
        
        if chi_squared_last[i] < min_chi_squared_last:
            min_chi_squared_last = chi_squared_last[i]
            min_chi_squared_index_last = i

        print('min_chi_squared_last is',min_chi_squared_last)       
        

    best_k = k_values_last[min_chi_squared_index_last]
    k = best_k
    best_B = B_last[min_chi_squared_index_last]
    best_eps_ni = eps_ni_last[min_chi_squared_index_last]
    best_R_squared = R_squared_last[min_chi_squared_index_last]
    best_chi_squared = chi_squared_last[min_chi_squared_index_last]
    best_L_Tau_fit = L_model_2param(Turn, best_eps_ni, best_B)

    print('i_last ---- best_k -- best_kB -- best_eps_n is', i, best_k, best_B, best_eps_ni )
    print('min_chi_squared_last and min_chi_squared_index_last is', min_chi_squared_index_last, min_chi_squared_last)

    with open('20{}/best_param/20{}_2param_last_best_FitCoeff.txt'.format(str(year),str(year)), 'a') as f:
        f.write(text)
        f.write(' ')
        f.write(str(best_eps_ni))
        f.write(' ')
        f.write(str(best_B))
        f.write(' ')
        f.write(str(best_k)) 
        f.write(' ')
        f.write(str(best_R_squared))
        f.write(' ')
        f.write(str(best_chi_squared))
        f.write('\n')

    fig, ax = plt.subplots() 
    ax.plot(T_fit_real/3600, L_Tau, "b.", markersize=4, label='Data')
    ax.plot(T_fit_real/3600, best_L_Tau_fit, "r-", markersize=4, label='2 param Model-2')
    ax.set_title('FillNumber: {}'.format(text))
    ax.set_xlabel(' Time [h]')
    ax.set_ylabel('${L}/{L_i}$' ) 
    plt.legend(loc='best')
    #plt.savefig('20{}/best_mod/{}_reduced_interval_Fit_mod2_best.pdf'.format(str(year),text))
    plt.savefig('20{}/best_mod/{}_2param_Fit.pdf'.format(str(year),text))
        
        
    fig, ax = plt.subplots() 
    ax.plot(T_fit_real/3600, L_Tau, "b.", markersize=4, label='Data')
    ax.plot(T_fit_real/3600, best_L_Tau_fit, "r-", markersize=4, label='2 param Model-2')
    ax.plot(T_fit_real/3600, L, "g-", markersize=4, label='3 param Model-2')
    ax.set_title('FillNumber: {}'.format(text))
    ax.set_xlabel(' Time [h]')
    ax.set_ylabel('${L}/{L_i}$' ) 
    plt.legend(loc='best')
    plt.savefig('20{}/best_mod/{}_2param_3param_Fits.pdf'.format(str(year),text))
        
        
    for x in range(len(k_values_last)):
        with open('20{}/Coeff/{}_2param_last_FitCoeff.txt'.format((str(year)),text), 'a') as f:
            f.write(str(k_values_last[x]))
            f.write(' ')
            f.write(str(eps_ni_last[x]))
            f.write(' ')
            f.write(str(B_last[x]))
            f.write(' ')
            f.write(str(R_squared_last[x]))
            f.write(' ')
            f.write(str(chi_squared_last[x]))
            f.write('\n')
            
    for x in range(len(k_values)):
        with open('20{}/Coeff/{}_2param_first_FitCoeff.txt'.format((str(year)),text), 'a') as f:
            f.write(str(k_values[x]))
            f.write(' ')
            f.write(str(eps_ni[x]))
            f.write(' ')
            f.write(str(B[x]))
            f.write(' ')
            f.write(str(R_squared[x]))
            f.write(' ')
            f.write(str(chi_squared[x]))
            f.write('\n')
            
            
    for x in range(len(T_fit_real)):
        with open('20{}/data/{}_2param_Fitdata.txt'.format((str(year)),text), 'a') as f:  
            f.write(str(T_fit_real[x]))
            f.write(' ')
            f.write(str(L_Tau[x]))
            f.write(' ')
            f.write(str(best_L_Tau_fit[x]))
            f.write('\n')
               
        
        
    fig, ax = plt.subplots()
    ax.plot(T_fit_real/3600, L_Tau, "b.", markersize=4, label='Data')  
                  
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, (T_fit_real, L_Tau_fit_last) in enumerate(subsets2):
        color = colors[i % len(colors)]
        #label = f'Subset [i % len(colors)]'
        #label = f'$\kappa$={k_values_last[i]}'
        label = f'k_values {i+1}'
        
        ax.plot(T_fit_real/3600, L_Tau_fit_last, color=color, label=label)
            
    ax.set_title('FillNumber: {}'.format(text))
    ax.set_xlabel(' Time [h]')
    ax.set_ylabel('${L}/{L_i}$' ) 
    #plt.legend(loc='best')
    plt.savefig('20{}/mod/{}_Fit_2param_scan.pdf'.format(str(year),text))
    #plt.show()


    delta_k = (best_k - k_all).tolist()
    delta_B = (best_B - B_all).tolist()
    delta_eps_ni = (best_eps_ni - eps_all).tolist()
    delta_R_squared = (best_R_squared - R_adj).tolist()
    delta_chi_squared = (best_chi_squared- chi_square).tolist()

    with open('20{}/best_param/20{}_delta_2param_Coeff.txt'.format(str(year),str(year)), 'a') as f:
        f.write(text)
        f.write(' ')
        f.write(str(delta_eps_ni))
        f.write(' ')
        f.write(str(delta_B))
        f.write(' ')
        f.write(str(delta_k))
        f.write(' ')
        f.write(str(delta_R_squared))
        f.write(' ')
        f.write(str(delta_chi_squared))
        f.write('\n')
    
   
#defining the stop time of the program      
stop=t.time()
#evaluating the run time 
runtime_seq=stop-start
print('The runtime is:', runtime_seq, '[s]')      

        


