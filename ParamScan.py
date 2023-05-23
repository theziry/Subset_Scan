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
from tqdm.autonotebook import tqdm


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

def R_2(original, Fit):
    """R-squared function
    Parameters
    ----------
    original : ndarray
        ideal values
    estimate : ndarray
        estimated values
    Returns
    -------
    float
        R-squared value
    -------
    This equivalent to calculate the R-squared value
    As fellow: R_squared = r2_score(L_Tau, L)
    """
    residuals = original - Fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((original - np.mean(original))**2)
    if ss_tot == 0:
        r_2 = 1
    else:
        r_2 = 1 - (ss_res / ss_tot)  
    return r_2


def L_model_3param(x, eps_ni, B, k):
        x = np.where(x>1, x, 1.1 ) 
        np.nan_to_num(x, copy=True, nan=0.0, posinf=None, neginf=None) #if np.isnan(x): x = 0.0 
        #D = B / (np.log(x) ** k)   
        D = B * np.power(k / (2 * np.exp(1)), k) / (np.power(np.log(x), k))
        L = (1/(1 + eps_ni * (x-1))**2) - ((1 + D**2)*np.exp(-(D**2))) * ((2 - (1 + D**2)*np.exp(-(D**2)/2)))
        return L 



#EXPLORE THE KAPPA AND FIT ONLY TWO PARAMETER (Rho and EPSILON*Ni): from tqdm.autonotebook import tqdm
def fit_2param(Turn,L_Tau,k_min, k_max, k_samples):
    k_values = np.linspace(k_min, k_max, k_samples)
    pars = []
    errs = []
    co_pars = []
    Radj_pars = []
    L_Tau_fit = []

    for k in tqdm(k_values):
        par, co_par = curve_fit(
            lambda x, eps_ni, B: L_model_3param(x,  eps_ni, B, k),
            Turn,
            L_Tau,
            p0=(9.5e-10, 155),
            bounds=([2e-10, 150], [2.e-9, 170.]),
            maxfev=5000
        )
        n = len(L_Tau)
        p = len(par) - 1
        r2 = R_2(L_Tau, L_model_3param(Turn, par[0], par[1], k))
        L_Tau_fit.append(L_model_3param(Turn, par[0], par[1], k))
        #print('L_Tau --- ',L_Tau, '----- L_Tau_fit ------', L_Tau_fit)
        pars.append(par)
        co_pars.append(co_par)
        errs.append(chi_2(L_Tau, L_model_3param(Turn, par[0], par[1], k)))
        Radj = (1 - (1 - r2) * (n - 1) / (n - p - 1))
        Radj_pars.append(Radj)  

    return np.asarray(L_Tau_fit), np.asarray(pars), np.asarray(errs), np.asarray(co_pars), np.array(Radj_pars)


def fit_1param(Turn,L_Tau,k_min, k_max, k_samples, B_min, B_max, B_samples):
    k_values = np.linspace(k_min, k_max, k_samples)
    B_values = np.linspace(B_min, B_max, B_samples)
    pars = []
    errs = []
    co_pars = []
    Radj_pars = []
    L_Tau_fit = []
    k_val = []
    B_val = []

    for k in tqdm(k_values):
        for B in tqdm(B_values):
            par, co_par = curve_fit(
                lambda x, eps_ni: L_model_3param(x,  eps_ni, B, k),
                Turn,
                L_Tau,
                p0=(9.5e-10),
                bounds=([2e-10], [2.e-9]),
                maxfev=5000
            )
            n = len(L_Tau)
            p = len(par) - 1
            r2 = R_2(L_Tau, L_model_3param(Turn, par[0], B, k))
            L_Tau_fit.append(L_model_3param(Turn, par[0],B, k))
            #print('L_Tau --- ',L_Tau, '----- L_Tau_fit ------', L_Tau_fit)
            k_val.append(k)
            B_val.append(B)
            pars.append(par)
            co_pars.append(co_par)
            errs.append(chi_2(L_Tau, L_model_3param(Turn, par[0],B, k)))
            Radj = (1 - (1 - r2) * (n - 1) / (n - p - 1))
            Radj_pars.append(Radj)  

    return np.array(L_Tau_fit),np.array(k_val),np.array(B_val), np.array(pars), np.array(errs), np.array(co_pars), np.array(Radj_pars)
   



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

    
    with open('20{}/best_param/20{}_1param_last_best_FitCoeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close()
        
    with open('20{}/best_param/20{}_1param_first_best_FitCoeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close()

    with open('20{}/best_param/20{}_2param_last_best_FitCoeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close()
        
    with open('20{}/best_param/20{}_2param_first_best_FitCoeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close()
        
    with open('20{}/best_param/20{}_3param_Coeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/best_param/20{}_delta_2param_Coeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close() 

    with open('20{}/best_param/20{}_delta_1param_Coeff.txt'.format(str(year),text), 'w') as f:
        f.write('')
        f.close()

    print('FillNumber',text)

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

    start_time = T_fit_real[0]
    #start_time = 3600
    print('T_fit_real[0] start_time', T_fit_real[0])
    end_time = T_fit_real[-1]
    print('T_fit_real[-1] end_time', T_fit_real[-1])

    step_size = 3600  # Fixed step size

    sub = []
    sub1 = []
    sub2 = []
    subsets = []
    subsets1 = []
    subsets2 = []
    subsets_1param = []
    subsets_2param = []

    #for start_time in range(int(start_time), int(end_time) + 1, int(step_size)):
    while start_time <= (end_time):
        # Select a subset of the data
        start_idx = int(start_time)
        end_idx = int(start_time + step_size)


        Time_subset = T_fit_real[T_fit_real <=end_idx]
        Time_subset[-1] = end_idx

        if end_idx >= T_fit_real[-1]:
            Time_subset[-1] = T_fit_real[-1]



        #L_Tau_subset = L_Tau[:len(Time_subset)]
        L_Tau_subset = L_Tau[T_fit_real <=end_idx]

        

        Turn_subset=[] 
        Turn_subset=np.array(Turn_subset)
        for el in Time_subset:
           tau=(f_rev*el+1)
           Turn_subset=np.append(Turn_subset, tau)
           
        
        
        #print('length Turn', len(Turn))
        print('length Time', len(T_fit_real))
        print('length Time_subset', len(Time_subset))
        print('length Turn_subset', len(Turn_subset))
        
        tim = int(end_idx)
        tmp = int(Turn_subset[-1])

        HR = tim/3600

        print('time time time ',(Time_subset[-1]))

        params, pcov = curve_fit(L_model_3param, Turn_subset,L_Tau_subset,p0,bounds=([2e-10, 100., 0.7], [2e-9, 200., 1.2]),maxfev=5000)
        eps_all, B_all, k_all = params


        # The fit output
        L = L_model_3param(Turn_subset, eps_all, B_all, k_all)

        

        # RESIDUALS 
        AbsErr = L_Tau_subset - L
        #print('absErr', AbsErr)
        SE = np.square(AbsErr) # squared errors
        MSE = np.mean(SE) # mean squared errors
        RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE

        # calculate the R-squared value
        r_squared = r2_score(L_Tau_subset, L)

        # number of observations
        n = len(L_Tau_subset)
        # number of degree of freedom
        p = len(params) - 1

        # calculate the adjusted R-squared value
        R_adj = 1 - (1 - r_squared) * (n - 1) / (n - p - 1)

        # calculate the chi_square value
        chi_square = sum(((L_Tau_subset - L)**2/L_Tau_subset))
        # calculate the Reduced_chi_square
        R_chi = chi_square/(len(L_Tau_subset) - 3)

        print('FillNumber', text,'R_squared', R_adj)




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

        #_____________________________________________________________________ NOW THE SCANNING OF KAPPA AND WE FIT THE TWO OTHER PARAMETER ____________________________________________________
                
        k_min = 0.7
        k_max = 1.2
        dk = 0.01
        k_samples = (int((k_max-k_min)/dk)+1)
        k_values = np.linspace(k_min, k_max, k_samples)

        
        
        real_L_Tau_fit,real_pars, real_errs, real_co_pars, real_Radj_pars  = fit_2param(
            Turn,L_Tau,k_min, k_max, k_samples
        )


        L_Tau_fit=real_L_Tau_fit
        eps_ni=real_pars[:, 0]
        B=real_pars[:,1]
        R_squared=real_Radj_pars
        chi_squared=real_errs

        sub2.append((T_fit_real, L_Tau_fit)) 

        min_err_index = np.argmin(real_errs)

        
            
        print('min_chi_squared is',min_err_index) 

        best_k0 = k_values[np.argmin(real_errs)]
        best_B0 = B[np.argmin(real_errs)]
        best_eps_ni0 = eps_ni[np.argmin(real_errs)]
        best_R_squared0 = R_squared[np.argmin(real_errs)]
        best_chi_squared0 = chi_squared[np.argmin(real_errs)]
        best_L_Tau_fit0 = L_model_3param(Turn_subset, best_eps_ni0, best_B0,best_k0)
        

        print('---- min_err_index ---- best_k0 -- best_kB0 -- best_eps_n0 is', min_err_index, best_k0, best_B0, best_eps_ni0 )
        print('min_chi_squared and min_chi_squared_index is',min_err_index, best_chi_squared0, 'best_R_squared0', min_err_index, best_R_squared0)

            
        k_min_last = best_k0 - 0.009
        k_max_last = best_k0 + 0.009
        dk_last = 0.001
        k_samples_last = (int((k_max_last-k_min_last)/dk_last)+1)
        k_values_last = np.linspace(k_min_last, k_max_last, k_samples_last)

        min_chi_squared_last = np.inf
        min_chi_squared_index_last = None
        
        real_L_Tau_fit_last, real_pars_last, real_errs_last, real_co_pars_last, real_Radj_pars_last  = fit_2param(
            Turn,L_Tau,k_min_last, k_max_last, k_samples_last
        )

        L_Tau_fit_last=(real_L_Tau_fit_last)
        eps_ni_last=(real_pars_last[:, 0])
        B_last=(real_pars_last[:, 1])
        R_squared_last=(real_Radj_pars_last)
        chi_squared_last=(real_errs_last)

                    

        last_min_err_index = np.argmin(chi_squared_last)

        print('last_min_chi_squared is',last_min_err_index)
    

        best_k = k_values_last[np.argmin(real_errs_last)]
        best_B = B_last[last_min_err_index]
        best_eps_ni = eps_ni_last[last_min_err_index]
        best_R_squared = R_squared_last[last_min_err_index]
        best_chi_squared = chi_squared_last[last_min_err_index]
        best_L_Tau_fit = L_model_3param(Turn_subset, best_eps_ni, best_B,best_k)

        print('last_index ---- best_k -- best_kB -- best_eps_n is', i, last_min_err_index, best_k, best_B, best_eps_ni )
        print('min_chi_squared_last and min_chi_squared_index_last is', last_min_err_index, min_chi_squared_last)


        #_____________________________________________________________________ NOW THE SCANNING OF KAPPA Rho AND WE FIT ONLY ONE PARAMETER EPSILON*Ni ____________________________________________________

        k1_min = 0.5
        k1_max = 1.5
        dk1 = 0.01
        k1_samples = (int((k1_max-k1_min)/dk1)+1)
        k1_values = np.linspace(k1_min, k1_max, k1_samples)

        B1_min = 100
        B1_max = 300
        dB1 = 1
        B1_samples = (int((B1_max-B1_min)/dB1)+1)
        B1_values = np.linspace(B1_min, B1_max, B1_samples)

            
        #min_chi_squared = np.inf
        #min_chi_squared_index = None
        
        real_L_Tau_fit1,real_k1_values,real_B1_values,real_pars1,real_errs1, real_co_pars1, real_Radj_pars1  = fit_1param(
            Turn,L_Tau,k1_min, k1_max, k1_samples, B1_min, B1_max, B1_samples
        )

    
        #print('real_pars --- ',real_pars[0], 'kkkkkkkkkk', k_values)

        k1_val=real_k1_values
        B1_val=real_B1_values

        L_Tau_fit1=real_L_Tau_fit1
        eps_ni1=real_pars1[:, 0]
        R_squared1=real_Radj_pars1
        chi_squared1=real_errs1

        sub1.append((T_fit_real, L_Tau_fit1)) 

        min_err_index1 = np.argmin(real_errs1)

        #min_index = np.unravel_index(np.argmin(real_errs),real_errs.shape)

        #print('Eps --- ',eps_ni, '----- B ------', B, 'kappa', k_values)

            
        print('min_chi_squared is',min_err_index1)


        best1_k0 = k1_val[min_err_index1]
        best1_B0 = B1_val[min_err_index1]
        best1_eps_ni0 = eps_ni1[min_err_index1]
        best1_R_squared0 = R_squared1[min_err_index1]
        best1_chi_squared0 = chi_squared1[min_err_index1]
        best1_L_Tau_fit0 = L_model_3param(Turn, best1_eps_ni0, best1_B0,best1_k0)

        #subsets1.append((T_fit_real, L_Tau_fit[i])) 
        
        

        print('---- min_err_index ---- best_k0 -- best_kB0 -- best_eps_n0 is', min_err_index1, best1_k0, best1_B0, best1_eps_ni0 )
        print('min_chi_squared and min_chi_squared_index is',min_err_index1, best1_chi_squared0, 'best_R_squared0', min_err_index1, best1_R_squared0)

             

            
        k1_min_last = best1_k0 - 0.009
        k1_max_last = best1_k0 + 0.009
        dk1_last = 0.001
        k1_samples_last = (int((k1_max_last-k1_min_last)/dk1_last)+1)
        k1_values_last = np.linspace(k1_min_last, k1_max_last, k1_samples_last)

        min_chi_squared_last = np.inf
        min_chi_squared_index_last = None
        
        B1_min_last = best1_B0 - 0.9
        B1_max_last = best1_B0 + 0.9
        dB1_last = 0.1
        B1_samples_last = (int((B1_max_last-B1_min_last)/dB1_last)+1)
        B1_values_last = np.linspace(B1_min_last, B1_max_last, B1_samples_last)

        
        real_L_Tau_fit1_last,real_k1_values_last,real_B1_values_last, real_pars1_last, real_errs1_last, real_co_pars1_last, real_Radj_pars1_last  = fit_1param(
            Turn,L_Tau,k1_min_last, k1_max_last, k1_samples_last, B1_min_last, B1_max_last, B1_samples_last
        )

        k1_val_last=real_k1_values_last
        B1_val_last=real_B1_values_last

        L_Tau_fit1_last=(real_L_Tau_fit1_last)
        eps_ni1_last=(real_pars1_last[:, 0])
        R_squared1_last=(real_Radj_pars1_last)
        chi_squared1_last=(real_errs1_last)
            

        last_min_err_index1 = np.argmin(chi_squared1_last)

        print('last_min_chi_squared is',last_min_err_index1)
    

        best_k1 = k1_val_last[np.argmin(real_errs1_last)]
        best_B1 = B1_val_last[last_min_err_index1]
        best_eps_ni1 = eps_ni1_last[last_min_err_index1]
        best_R_squared1 = R_squared1_last[last_min_err_index1]
        best_chi_squared1 = chi_squared1_last[last_min_err_index1]
        best_L_Tau_fit1 = L_model_3param(Turn, best_eps_ni1, best_B1,best_k1)

        print('last_index ---- best_k -- best_kB -- best_eps_n is', i, last_min_err_index1, best_k1, best_B1, best_eps_ni1)
        print('min_chi_squared_last and min_chi_squared_index_last is', last_min_err_index1, '---np.argmin(chi_squared_last--',np.argmin(chi_squared1_last))


        subsets.append((Time_subset, L))  
        subsets1.append((Time_subset, best1_L_Tau_fit0)) 
        subsets_1param.append((Time_subset, best_L_Tau_fit1))
        subsets2.append((Time_subset, best_L_Tau_fit0)) 
        subsets_2param.append((Time_subset, best_L_Tau_fit))


        fig, ax = plt.subplots()
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                    
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, L) in enumerate(subsets):
            color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f'subset {i+1}'
            
            ax.plot(Time_subset/3600, L, color=color, label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/{}_3param_Fit.pdf'.format(str(year),text))
        plt.close('all') 
    
          

        with open('20{}/best_param/{}_3param_Coeff.txt'.format(str(year),text), 'a') as f:
            f.write(str(tim))
            f.write(' ')
            f.write(str(tmp))
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


        fig, ax = plt.subplots()
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, best_L_Tau_fit0) in enumerate(subsets2):
            color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f' {i+1}'
            
            ax.plot(Time_subset/3600, best_L_Tau_fit0, color=color, label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/first_check/{}_2param_Fit.pdf'.format(str(year),text))
        plt.close('all')



        with open('20{}/best_param/20{}_2param_first_best_FitCoeff.txt'.format(str(year),text), 'a') as f:
            f.write(str(tim))
            f.write(' ')
            f.write(str(tmp))
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

        

        fig, ax = plt.subplots()
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, best_L_Tau_fit) in enumerate(subsets_2param):
            color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f' {i+1}'
            
            ax.plot(Time_subset/3600, best_L_Tau_fit, color=color, label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/{}_2param_Fit.pdf'.format(str(year),text))
        plt.close('all')


        fig, ax = plt.subplots()
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        #for i, (Time_subset, best_L_Tau_fit) in enumerate(subsets_2param):
        for i, (subset, subset_2param) in enumerate(zip(subsets, subsets_2param)):
            Time_subset, L = subset
            Time_subset_2param, best_L_Tau_fit = subset_2param
    
            #color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f' {i+1}'
            
            ax.plot(Time_subset/3600, best_L_Tau_fit, color='red', label=label)
            ax.plot(Time_subset/3600, L, color='lime', label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/{}_2param_3param_Fits.pdf'.format(str(year),text))
        plt.close('all')


        

        with open('20{}/best_param/20{}_2param_last_best_FitCoeff.txt'.format(str(year),text), 'a') as f:
            f.write(str(tim))
            f.write(' ')
            f.write(str(tmp))
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
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, (Time_subset, best1_L_Tau_fit0) in enumerate(subsets_2param):
            color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f' {i+1}'
            
            ax.plot(Time_subset/3600, best1_L_Tau_fit0, color=color, label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/first_check/{}_1param_Fit.pdf'.format(str(year),text))
        plt.close('all')

        fig, ax = plt.subplots()
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        #for i, (Time_subset, best_L_Tau_fit) in enumerate(subsets_2param):
        for i, (subset, subset_1param) in enumerate(zip(subsets, subsets_1param)):
            Time_subset, L = subset
            Time_subset_1param, best_L_Tau_fit1 = subset_1param
    
            #color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f' {i+1}'
            
            ax.plot(Time_subset/3600, best_L_Tau_fit1, color='red', label=label)
            ax.plot(Time_subset/3600, L, color='lime', label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/{}_1param_3param_Fits.pdf'.format(str(year),text))
        plt.close('all')


        fig, ax = plt.subplots()
        ax.plot(Time_subset/3600, L_Tau_subset, "b.", markersize=4, label='Data')  
                  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        #for i, (Time_subset, best_L_Tau_fit) in enumerate(subsets_2param):
        for i, (subset, subset_1param) in enumerate(zip(subsets, subsets_1param)):
            Time_subset, L = subset
            Time_subset_1param, best_L_Tau_fit1 = subset_1param
            Time_subset_2param, best_L_Tau_fit2 = subset_2param
    
            #color = colors[i % len(colors)]
            #label = f'Subset {i+1}'
            label = f' {i+1}'
            
            ax.plot(Time_subset/3600, best_L_Tau_fit1, color='red', label=label)
            ax.plot(Time_subset/3600, best_L_Tau_fit2, color='yellow', label=label)
            ax.plot(Time_subset/3600, L, color='lime', label=label)
            
        ax.set_title('FillNumber: {}'.format(text))
        ax.set_xlabel(' Time [h]')
        ax.set_ylabel('${L}/{L_i}$' ) 
        #plt.legend(loc='best')
        plt.savefig('20{}/best_mod/{}_1param_2param_3param_Fits.pdf'.format(str(year),text))
        plt.close('all')
            


        with open('20{}/best_param/20{}_1param_first_best_FitCoeff.txt'.format(str(year),text), 'a') as f:
            f.write(str(tim))
            f.write(' ')
            f.write(str(tmp))
            f.write(' ')
            f.write(str(best1_eps_ni0))
            f.write(' ')
            f.write(str(best1_B0))
            f.write(' ')
            f.write(str(best1_k0)) 
            f.write(' ')
            f.write(str(best1_R_squared0))
            f.write(' ')
            f.write(str(best1_chi_squared0))
            f.write('\n') 

        with open('20{}/best_param/20{}_1param_last_best_FitCoeff.txt'.format(str(year),text), 'a') as f:
            f.write(str(tim))
            f.write(' ')
            f.write(str(tmp))
            f.write(' ')
            f.write(str(best_eps_ni1))
            f.write(' ')
            f.write(str(best_B1))
            f.write(' ')
            f.write(str(best_k1)) 
            f.write(' ')
            f.write(str(best_R_squared1))
            f.write(' ')
            f.write(str(best_chi_squared1))
            f.write('\n') 


        

        

        start_time += step_size
        
        

   
   
#defining the stop time of the program      
stop=t.time()
#evaluating the run time 
runtime_seq=stop-start
print('The runtime is:', runtime_seq, '[s]')  