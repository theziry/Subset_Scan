import numpy as np
import LoadData as ld
import pandas as pd
import LuminosityOptimization as lo
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import matplotlib.ticker as ticker
import time as t
#from lmfit import Model
from scipy import integrate
import scipy.optimize
from scipy.optimize import curve_fit
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
# create an empty list to store the L_fit data for each fill number

year = 18

#Plotting Luminosity Evolution Graphs
plot=True

#model parameters 

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

#load turnaround times and fill times 
data_ta16, data_tf16, data_ta17, data_tf17, data_ta18, data_tf18 = ld.loadFill()
data_ta16_sec = data_ta16*3600 
data_tf16_sec = data_tf16*3600  
data_ta17_sec = data_ta17*3600 
data_tf17_sec = data_tf17*3600
data_ta18_sec = data_ta18*3600 
data_tf18_sec = data_tf18*3600

skip16=[5111,5112,5116,5117,5257,5258,5264,5406,5427,5439,5451]
skip17=[5837,5840,5856,6019,6055,6060,6082,6084,6089,6116,6152,6156,6158,6160,6167,6168,6169,6170,6193,6258,6268,6271,6272,6275,6283,6285,6287,6288,6291,6298,6303,6304,6305,6308,6311,6312,6314,6317,6324,6325,6337,6346,6356,6360,6362,6364,6371]
skip18=[6640,6654,6659,6688,6690,6693,6694,6696,6700,6706,6729,6741,6747,6752,6782,6778,6882,6890,6892,6901,6929,6939,6940,6961,7026,7033,7036,7039,7056,7069,7105,7118,7120,7135,7145,7217,7245,7259,7271,7324]
skip = skip16 + skip17 + skip18
for fill in skip:
    FillNumber=np.delete(FillNumber, (np.where(FillNumber==(fill))[0]))
for i in range(len(FillNumber)):
    if FillNumber[i] in skip:
        continue


for i in range(len(FillNumber)):
    text = str(int(FillNumber[i])) #number of current fill

    with open('20{}/best_param/20{}_1param_last_best_FitCoeff.txt'.format(str(year),str(year)), 'r') as f:
        lines=f.readlines()
        fill_scan=[]
        k_scan=[] #fill number
        eps_ni_scan=[]
        B_scan=[]
        R_scan=[]
        chi_scan=[]

        


        for x in lines:
            fill_scan.append(float(x.split(' ')[0]))
            eps_ni_scan.append(float(x.split(' ')[1]))
            B_scan.append(float(x.split(' ')[2]))
            k_scan.append(float(x.split(' ')[3])) 
            R_scan.append(float(x.split(' ')[4]))
            chi_scan.append(float(x.split(' ')[5])) 

    fill_list_scan = np.array(fill_scan)     
    Eps_list_scan = np.array(eps_ni_scan)
    B_list_scan = np.array(B_scan)
    k_list_scan = np.array(k_scan)
    R_list_scan = np.array(R_scan)
    chi_list_scan = np.array(chi_scan)


    with open('20{}/best_param/20{}_2param_last_best_FitCoeff.txt'.format(str(year),str(year)), 'r') as f:
        lines=f.readlines()
        fill_scan2=[]
        k_scan2=[] #fill number
        eps_ni_scan2=[]
        B_scan2=[]
        R_scan2=[]
        chi_scan2=[]

        


        for x in lines:
            fill_scan2.append(float(x.split(' ')[0]))
            eps_ni_scan2.append(float(x.split(' ')[1]))
            B_scan2.append(float(x.split(' ')[2]))
            k_scan2.append(float(x.split(' ')[3])) 
            R_scan2.append(float(x.split(' ')[4]))
            chi_scan2.append(float(x.split(' ')[5])) 

    fill_list_scan2 = np.array(fill_scan2)     
    Eps_list_scan2 = np.array(eps_ni_scan2)
    B_list_scan2 = np.array(B_scan2)
    k_list_scan2 = np.array(k_scan2)
    R_list_scan2 = np.array(R_scan2)
    chi_list_scan2 = np.array(chi_scan2)


    with open('20{}/best_param/20{}_3param_Coeff.txt'.format(str(year),str(year)), 'r') as f:
        lines=f.readlines()
        fill=[]
        k=[] #fill number
        eps_ni=[]
        B=[]
        R=[]
        chi=[]

        


        for x in lines:
            fill.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[1]))
            B.append(float(x.split(' ')[2]))
            k.append(float(x.split(' ')[3])) 
            R.append(float(x.split(' ')[4]))
            chi.append(float(x.split(' ')[5])) 

    fill_list = np.array(fill)     
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    k_list = np.array(k)
    R_list = np.array(R)
    chi_list = np.array(chi)

    with open('20{}/best_param/20{}_delta_1param_Coeff.txt'.format(str(year),str(year)), 'r') as f:
        lines=f.readlines()
        fillnbr=[]
        dk=[] #fill number
        deps_ni=[]
        dB=[]
        dR=[]
        dchi=[]

         
        for x in lines:
            fillnbr.append(float(x.split(' ')[0]))
            deps_ni.append(float(x.split(' ')[1]))
            dB.append(float(x.split(' ')[2]))
            dk.append(float(x.split(' ')[3])) 
            dR.append(float(x.split(' ')[4]))
            dchi.append(float(x.split(' ')[5])) 

    fillNbr_list = np.array(fillnbr)     
    delta_eps_ni = np.array(deps_ni)
    delta_B = np.array(dB)
    delta_k = np.array(dk)
    delta_R = np.array(dR)*100
    delta_chi = np.array(dchi)


    with open('20{}/best_param/20{}_delta_2param_Coeff.txt'.format(str(year),str(year)), 'r') as f:
        lines=f.readlines()
        fillnbr2=[]
        dk2=[] #fill number
        deps_ni2=[]
        dB2=[]
        dR2=[]
        dchi2=[]

         
        for x in lines:
            fillnbr2.append(float(x.split(' ')[0]))
            deps_ni2.append(float(x.split(' ')[1]))
            dB2.append(float(x.split(' ')[2]))
            dk2.append(float(x.split(' ')[3])) 
            dR2.append(float(x.split(' ')[4]))
            dchi2.append(float(x.split(' ')[5])) 

    fillNbr_list2 = np.array(fillnbr2)     
    delta_eps_ni2 = np.array(deps_ni2)
    delta_B2 = np.array(dB2)
    delta_k2 = np.array(dk2)
    delta_R2 = np.array(dR2)*100
    delta_chi2 = np.array(dchi2)

#__________________________________________________ 1 PARAMETER FIT PLOTS - 1 PARAMETER FIT PLOTS - 3 PARAMETER FIT PLOTS ______________________________________________

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, R_list_scan*100, color='red', marker='d', s=25, label='$1$-$Fit$ $Parameter$')
ax.scatter(fill_list_scan, R_list_scan2*100, color='green', marker='*', s=35, label='$2$-$Fit$ $Parameters$')
ax.scatter(fill_list_scan, R_list*100, color='blue', marker='o', s=20, label='$3$-$Fit$ $Parameters$')
ax.set_title('$20%s$: ${R_{adj}}$$^2$ for $Model$-$2$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}' ) 
plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_R_123Fit_Param.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_R_123Fit_Param.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, chi_list_scan, color='red', marker='d', s=25, label='$1$-$Fit$ $Parameter$')
ax.scatter(fill_list_scan, chi_list_scan2, color='green', marker='*', s=35, label='$2$-$Fit$ $Parameters$')
ax.scatter(fill_list_scan, chi_list, color='blue', marker='o', s=20, label='$3$-$Fit$ $Parameters$')
ax.set_title('$20%s$: ${\chi}^2$ for $Model$-$2$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${\chi}^2$') 
plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_chi_123Fit_Param.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_chi_123Fit_Param.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, k_list_scan, color='red', marker='d', s=25, label='$1$-$Fit$ $Parameter$')
ax.scatter(fill_list_scan, k_list_scan2, color='green', marker='*', s=35, label='$2$-$Fit$ $Parameters$')
ax.scatter(fill_list_scan, k_list, color='blue', marker='o', s=20, label='$3$-$Fit$ $Parameters$')
ax.set_title('$20%s$: $\kappa$ for $Model$-$2$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('$\kappa$') 
plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_k_123Fit_Param.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_k_123Fit_Param.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, B_list_scan, color='red', marker='d', s=25, label='$1$-$Fit$ $Parameter$')
ax.scatter(fill_list_scan, B_list_scan2, color='green', marker='*', s=35, label='$2$-$Fit$ $Parameters$')
ax.scatter(fill_list_scan, B_list, color='blue', marker='o', s=20, label='$3$-$Fit$ $Parameters$')
ax.set_title(''r'$20%s$: $\rho_\ast$ for $Model$-$2$' % str(year))    
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'$\rho_\ast$') 
plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_Rho_123Fit_Param.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Rho_123Fit_Param.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, Eps_list_scan*1e9, color='red', marker='d', s=25, label='$1$-$Fit$ $Parameter$')
ax.scatter(fill_list_scan, Eps_list_scan2*1e9, color='green', marker='*', s=35, label='$2$-$Fit$ $Parameters$')
ax.scatter(fill_list_scan, Eps_list*1e9, color='blue', marker='o', s=20, label='$3$-$Fit$ $Parameters$')
ax.set_title(''r'$20%s$: ${\epsilon \cdot N_i}$ for $Model$-$2$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r' ${\epsilon \cdot N_i}${(10$^{-9}$)}') 
plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_Eps_123Fit_Param.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Eps_123Fit_Param.png'.format(str(year),str(year)))
plt.close('all')

##________________ COMPARE DELTA ____________________##


fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_k, color='red', marker='.', s=50, label=r'${\Delta \kappa}_{\mathrm{(Difference\ between\ 1\ and\ 3\ Fit\ Parameters)}}$')
ax.scatter(fill_list, delta_k2, color='green', marker='.', s=50, label=r'${\Delta \kappa}_{\mathrm{(Difference\ between\ 2\ and\ 3\ Fit\ Parameters)}}$')
#ax.set_title(r'$20%s$: $\Delta \kappa$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \kappa$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_delta_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_delta_k.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_B, color='red', marker='.', s=50, label=r'${\Delta \rho_\ast}_{\mathrm{(Difference\ between\ 1\ and\ 3\ Fit\ Parameters)}}$')
ax.scatter(fill_list, delta_B2, color='green', marker='.', s=50,label=r'${\Delta \rho_\ast}_{\mathrm{(Difference\ between\ 2\ and\ 3\ Fit\ Parameters)}}$')
#ax.set_title(r'$20%s$: $\Delta \rho_\ast$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \rho_\ast$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_delta_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_delta_Rho.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_eps_ni*1e9, color='red', marker='.', s=50, label=r'${\Delta {\epsilon \cdot N_i}}_{\mathrm{(Difference\ between\ 1\ and\ 3\ Fit\ Parameters)}}$')
ax.scatter(fill_list, delta_eps_ni2*1e9, color='green', marker='.', s=50, label=r'${\Delta {\epsilon \cdot N_i}}_{\mathrm{(Difference\ between\ 2\ and\ 3\ Fit\ Parameters)}}$')
#ax.set_title(r'$20%s$: $\Delta {\epsilon}*{N_i}$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta {\epsilon \cdot N_i}${(10$^{-9}$)}')  
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_delta_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_delta_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_R, color='red', marker='.', s=50, label=r'${\Delta {R_{adj}}^2}_{\mathrm{(Difference\ between\ 1\ and\ 3\ Fit\ Parameters)}}$')
ax.scatter(fill_list, delta_R2, color='green', marker='.', s=50, label=r'${\Delta {R_{adj}}^2}_{\mathrm{(Difference\ between\ 2\ and\ 3\ Fit\ Parameters)}}$')
#ax.set_title(r'$20%s$: $\Delta {R_{adj}}^2$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.165), loc='upper center')
ax.set_xlabel(' FillNumber') 
ax.set_ylabel(r'$\Delta {R_{adj}}^2{({\%})}$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_delta_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_delta_R.png'.format(str(year),str(year)))

plt.close('all')
##________________ COMPARE RELATIVE DIFFERENCES ____________________## ax.set_ylabel(''r'${({\epsilon \cdot N_i})_{\mathrm{ (1 , 2)\ Fit\ Param}}}/{({\epsilon \cdot N_i})_{\mathrm{ 3 Fit\ Param}}}$' ) 


fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_eps_ni/Eps_list, color='red', marker='.', s=50, label=r'${\Delta ({\epsilon \cdot N_i})_{\mathrm{Difference\ between\ 1\ and\ 3\ Fit\ Parameters}} / ( {\epsilon \cdot N_i})}_{\mathrm{ 3\ Fit\ Parameters}}$')
ax.scatter(fill_list, delta_eps_ni2/Eps_list, color='green', marker='.', s=50, label=r'${\Delta ({\epsilon \cdot N_i})_{\mathrm{ Difference\ between\ 2\ and\ 3\ Fit\ Parameters}} / ( {\epsilon \cdot N_i})}_{\mathrm{ 3\ Fit\ Parameters}}$')
#ax.set_title(r'$20%s$: $(\Delta {\epsilon \cdot N_i})/({\epsilon \cdot N_i})$' % str(year))
# Set the legend outside the graph
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$(\Delta {\epsilon \cdot N_i})/({\epsilon \cdot N_i})$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_delta_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_delta_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_B/B_list, color='red', marker='.', s=50, label=r'${{\Delta \rho_\ast}_{\mathrm{ Difference\ between\ 1\ and\ 3\ Fit\ Parameters}} / \rho_\ast}_{\mathrm{ 3\ Fit\ Parameters}}$')
ax.scatter(fill_list, delta_B2/B_list, color='green', marker='.', s=50, label=r'${{\Delta \rho_\ast}_{\mathrm{ Difference\ between\ 2\ and\ 3\ Fit\ Parameters}} / \rho_\ast}_{\mathrm{ 3\ Fit\ Parameters}}$')
#ax.set_title(r'$20%s$: $\Delta \rho_\ast / \rho_\ast$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$\Delta \rho_\ast / \rho_\ast$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_delta_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_delta_Rho.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_k/k_list, color='red', marker='.', s=50, label=r'${{\Delta \kappa}_{\mathrm{ Difference\ between\ 1\ and\ 3\ Fit\ Parameters}} / \kappa}_{\mathrm{3\ Fit\ Parameters}}$')
ax.scatter(fill_list, delta_k2/k_list, color='green', marker='.', s=50, label=r'${{\Delta \kappa}_{\mathrm{ Difference\ between\ 1\ and\ 3\ Fit\ Parameters}} / \kappa}_{\mathrm{3\ Fit\ Parameters}}$')
#ax.set_title(r'$20%s$: $\Delta \kappa / \kappa$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$\Delta \kappa / \kappa $' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_delta_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_delta_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, delta_R/(R_list*100), color='red', marker='.', s=50, label=r'${{\Delta {R_{adj}}^2}_{\mathrm{ Difference\ between\ 1\ and\ 3\ Fit\ Parameters}} / {R_{adj}}^2}_{\mathrm{3\ Fit\ Parameters}}$')
ax.scatter(fill_list, delta_R2/(R_list*100), color='green', marker='.', s=50, label=r'${{\Delta {R_{adj}}^2}_{\mathrm{ Difference\ between\ 2\ and\ 3\ Fit\ Parameters}} / {R_{adj}}^2}_{\mathrm{3\ Fit\ Parameters}}$')
#ax.set_title(r'$20%s$: $\Delta {R_{adj}}^2 / {R_{adj}}^2$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.165), loc='upper center')
ax.set_xlabel('FillNumber')
ax.set_ylabel(r'$\Delta {R_{adj}}^2 / {R_{adj}}^2$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_delta_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_delta_R.png'.format(str(year),str(year)))

plt.close('all')



fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, Eps_list_scan/Eps_list, color='red',marker='.', s=50, label=r'${{(\epsilon \cdot N_i)}_{\mathrm{ 1\ Fit\ Parameter}}}/{{(\epsilon \cdot N_i)}_{\mathrm{ 3\ Fit\ Parameters}}}$')
ax.scatter(fill_list, Eps_list_scan2/Eps_list, color='green',marker='.', s=50, label=r'${{(\epsilon \cdot N_i)}_{\mathrm{ 2\ Fit\ Parameters}}}/{{(\epsilon \cdot N_i)}_{\mathrm{ 3\ Fit\ Parameters}}}$')
#ax.set_title(r'$20%s$: ${({\epsilon \cdot N_i})_{\mathrm{ (1 , 2)\ Fit\ Param}}}/{({\epsilon \cdot N_i})_{\mathrm{ 3\ Fit\ Param}}}$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel('FillNumber' )
ax.set_ylabel(''r'${({\epsilon \cdot N_i})_{\mathrm{ (1 , 2)\ Fit\ Param}}}/{({\epsilon \cdot N_i})_{\mathrm{ 3 Fit\ Param}}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, B_list_scan/B_list, color='red',marker='.', s=50, label=r'${{\rho_\ast}_{\mathrm{ (1\ Fit\ Parameter)}}}/{{\rho_\ast}_{\mathrm{ (3\ Fit\ Parameters)}}}$')
ax.scatter(fill_list, B_list_scan2/B_list, color='green',marker='.', s=50, label=r'${{\rho_\ast}_{\mathrm{ (2\ Fit\ Parameters)}}}/{{\rho_\ast}_{\mathrm{ (3\ Fit\ Parameters)}}}$')
#ax.set_title(r'$20%s$: ${\rho_\ast_{\mathrm{ (1 , 2)\ Fit\ Param}}}/{\rho_\ast_{\mathrm{ 3\ Fit\ Param}}}$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel('FillNumber')
#ax.set_ylabel(''r'$d \rho_\ast$' ) 
ax.set_ylabel(r'${{\rho_\ast}_{\mathrm{ ((1 , 2)\ Fit\ Param)}}}/{{\rho_\ast}_{\mathrm{ (3\ Fit\ Param)}}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_Rho.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, k_list_scan/k_list, color='red', marker='.', s=50, label=r'${{\kappa}_{\mathrm{ (1\ Fit\ Parameter)}}}/{{\kappa}_{\mathrm{ (3\ Fit\ Parameters)}}}$')
ax.scatter(fill_list, k_list_scan2/k_list, color='green', marker='.', s=50, label=r'${{\kappa}_{\mathrm{ (2\ Fit\ Parameters)}}}/{{\kappa}_{\mathrm{ (3\ Fit\ Parameters)}}}$')
#ax.set_title(r'$20%s$: ${{\kappa}_{\mathrm{ (1 , 2)\ Fit\ Param}}/{{\kappa}_{\mathrm{ 3\ Fit\ Param}}$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.15), loc='upper center')
ax.set_xlabel('FillNumber')
#ax.set_ylabel(''r'$d \kappa$' )
ax.set_ylabel(r'${{\kappa}_{\mathrm{ ((1 , 2)\ Fit\ Param)}}}/{{\kappa}_{\mathrm{ (3\ Fit\ Param)}}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots(figsize=(10, 8)) 
ax.scatter(fill_list, R_list_scan/R_list, color='red', marker='.', s=50, label=r'${{{R_{adj}}^2}_{\mathrm{ (1\ Fit\ Parameter)}}}/{{{R_{adj}}^2}_{\mathrm{ (3\ Fit\ Parameters)}}}$')
ax.scatter(fill_list, R_list_scan2/R_list, color='green', marker='.', s=50, label=r'${{{R_{adj}}^2}_{\mathrm{ (2\ Fit\ Parameters)}}}/{{{R_{adj}}^2}_{\mathrm{ (3\ Fit\ Parameters)}}}$')
#ax.set_title(r'$20%s$: ${{{R_{adj}}^2}_{\mathrm{ (1 , 2)\ Fit\ Param}}/{{{R_{adj}}^2}_{\mathrm{ 3\ Fit\ Param}}$' % str(year))
legend = ax.legend(title=r'20%s' % str(year), bbox_to_anchor=(0.5, 1.165), loc='upper center')
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'${{{R_{adj}}^2}_{\mathrm{ ((1 , 2)\ Fit\ Param)}}}/{{{R_{adj}}^2}_{\mathrm{ (3\ Fit\ Param)}}}$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/compare/{}_123Fit_Param_relative_diff_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_123Fit_Param_relative_diff_R.png'.format(str(year),str(year)))

plt.close('all')



#_______________________________________________________________________________________ 1 PARAMETER FIT PLOTS ______________________________________________

fig, ax = plt.subplots() 
#ax.plot(fill_list, R_list*100, "g.", markersize=6, label='FillsYear')
ax.scatter(fill_list_scan, R_list*100, color='blue', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${R_{adj}}$$^2$ for $Model$-$2$ with $3$-$Fit Parameters$'% str(year))
ax.set_xlabel('FillNumbers')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_R_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_R_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, R_list_scan*100, color='red', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${R_{adj}}$$^2$ for $Model$-$2$ with $1$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_R_1parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_R_1parm.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
#ax.plot(fill_list, chi_list, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list_scan, chi_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${\chi}^2$ for $Model$-$2$ with $3$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${\chi}^2$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_chi_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_chi_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, chi_list_scan,color='red', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${\chi}^2$ for $Model$-$2$ with $1$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${\chi}^2$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_chi_1parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_chi_1parm.png'.format(str(year),str(year)))
plt.close('all')



#### Parameters
fig, ax = plt.subplots() 
#ax.plot(fill_list, k_list, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list_scan, k_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: $\kappa$ for $Model$-$2$ with $3$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('$\kappa$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_k_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_k_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, k_list_scan,color='red', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: $\kappa$ for $Model$-$2$ with $1$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('$\kappa$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_k_1parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_k_1parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
#ax.plot(fill_list, B_list, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list_scan, B_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: $\rho_\ast$ for $Model$-$2$ with $3$-$Fit Parameters$' % str(year))    
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'$\rho_\ast$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Rho_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Rho_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, B_list_scan, color='red', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: $\rho_\ast$ for $Model$-$2$ with $1$-$Fit Parameters$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'$\rho_\ast$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Rho_1parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Rho_1parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
#ax.plot(fill_list, Eps_list*1e9, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list_scan, Eps_list*1e9, color='blue', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: ${\epsilon}*{N_i}$ for $Model$-$2$ with $3$-$Fit Parameter$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r' ${\epsilon}*{N_i}${(10$^{-9}$)}') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Eps_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Eps_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan, Eps_list_scan*1e9, color='red', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: ${\epsilon}*{N_i}$ for $Model$-$2$ with $1$-$Fit Parameter$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r' ${\epsilon}*{N_i}${(10$^{-9}$)}') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Eps_1parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Eps_1parm.png'.format(str(year),str(year)))
plt.close('all')



##________________ DELTA ____________________##


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_k, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \kappa$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \kappa$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_delta_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_delta_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_k, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \kappa$' % str(year))
ax.set_ylim([-0.101, 0.101])
ax.set_yticks(np.linspace(-0.1, 0.1, 11))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \kappa$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_delta_k_zoom.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_delta_k_zoom.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_B, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \rho_\ast$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \rho_\ast$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_delta_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_delta_Rho.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_eps_ni*1e9, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta {\epsilon}*{N_i}$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta {\epsilon}*{N_i}${(10$^{-9}$)}' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_delta_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_delta_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_R, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta {R_{adj}}^2$' % str(year))
ax.set_xlabel(' FillNumber')
#ax.set_ylabel(''r'$\Delta {R_{adj}}$^{2}$${({\%})}' ) 
ax.set_ylabel(r'$\Delta {R_{adj}}^2{({\%})}$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_delta_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_delta_R.png'.format(str(year),str(year)))

plt.close('all')

##________________ RELATIVE DIFFERENCES ____________________##


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_eps_ni/Eps_list, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $(\Delta {\epsilon}*{N_i})/({\epsilon}*{N_i})$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$(\Delta {\epsilon}*{N_i})/({\epsilon}*{N_i})$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_delta_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_delta_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_B/B_list, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \rho_\ast / \rho_\ast$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$\Delta \rho_\ast / \rho_\ast$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_delta_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_delta_Rho.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_k/k_list, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \kappa / \kappa$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$\Delta \kappa / \kappa $' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_delta_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_delta_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_R/(R_list*100), color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta {R_{adj}}^2 / {R_{adj}}^2$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(r'$\Delta {R_{adj}}^2 / {R_{adj}}^2$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_delta_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_delta_R.png'.format(str(year),str(year)))

plt.close('all')




fig, ax = plt.subplots() 
ax.scatter(fill_list, Eps_list_scan/Eps_list, color='red',marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${({\epsilon}*{N_i})_{1 param}}/{({\epsilon}*{N_i})_{3 param}}$' % str(year))
ax.set_xlabel('FillNumber' )
#ax.set_ylabel(''r'$d {\epsilon}*{N_i}$' ) 
ax.set_ylabel(''r'${({\epsilon}*{N_i})_{1 param}}/{({\epsilon}*{N_i})_{3 param}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, B_list_scan/B_list, color='red',marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${{\rho_\ast}_{1 param}}/{{\rho_\ast}_{3 param}}$' % str(year))
ax.set_xlabel('FillNumber')
#ax.set_ylabel(''r'$d \rho_\ast$' ) 
ax.set_ylabel(''r'${{\rho_\ast}_{1 param}}/{{\rho_\ast}_{3 param}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_Rho.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, k_list_scan/k_list, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${{\kappa}_{1 param}}/{{\kappa}_{3 param}}$' % str(year))
ax.set_xlabel('FillNumber')
#ax.set_ylabel(''r'$d \kappa$' )
ax.set_ylabel(''r'${{\kappa}_{1 param}}/{{\kappa}_{3 param}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_1param_relative_diff_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, R_list_scan/R_list, color='red', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${{{R_{adj}}^2}_{1 param}}/{{{R_{adj}}^2}_{3 param}}$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'${{{R_{adj}}^2}_{1 param}}/{{{R_{adj}}^2}_{3 param}}$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_1param_relative_diff_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}1param_relative_diff_R.png'.format(str(year),str(year)))

plt.close('all')

#_______________________________________________________________________________________ 2 PARAMETER FIT PLOTS ______________________________________________

fig, ax = plt.subplots() 
#ax.plot(fill_list, R_list*100, "g.", markersize=6, label='FillsYear')
ax.scatter(fill_list, R_list*100, color='blue', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${R_{adj}}$$^2$ for $Model$-$2$ with $3$-$Fit Parameters$'% str(year))
ax.set_xlabel('FillNumbers')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_R_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_R_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan2, R_list_scan2*100, color='green', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${R_{adj}}$$^2$ for $Model$-$2$ with $2$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${R_{adj}}$$^2${({\%})}' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_R_2parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_R_2parm.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
#ax.plot(fill_list, chi_list, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list, chi_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${\chi}^2$ for $Model$-$2$ with $3$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${\chi}^2$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_chi_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_chi_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan2, chi_list_scan2, color='green', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: ${\chi}^2$ for $Model$-$2$ with $2$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('${\chi}^2$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_chi_2parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_chi_2parm.png'.format(str(year),str(year)))
plt.close('all')



#### Parameters
fig, ax = plt.subplots() 
#ax.plot(fill_list, k_list, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list, k_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: $\kappa$ for $Model$-$2$ with $3$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('$\kappa$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_k_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_k_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan2, k_list_scan2, color='green', marker='.', s=50, label='delta_k')
ax.set_title('$20%s$: $\kappa$ for $Model$-$2$ with $2$-$Fit Parameters$'% str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel('$\kappa$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_k_2parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_k_2parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
#ax.plot(fill_list, B_list, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list, B_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: $\rho_\ast$ for $Model$-$2$ with $3$-$Fit Parameters$' % str(year))    
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'$\rho_\ast$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Rho_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Rho_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan2, B_list_scan2, color='green', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: $\rho_\ast$ for $Model$-$2$ with $2$-$Fit Parameters$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'$\rho_\ast$') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Rho_2parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Rho_2parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
#ax.plot(fill_list, Eps_list*1e9, "g.", markersize=6, label='delta_k')
ax.scatter(fill_list, Eps_list, color='blue', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: ${\epsilon}*{N_i}$ for $Model$-$2$ with $3$-$Fit Parameter$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r' ${\epsilon}*{N_i}${(10$^{-9}$)}') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Eps_3parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Eps_3parm.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list_scan2, Eps_list_scan2*1e9,color='green', marker='.', s=50, label='delta_k')
ax.set_title(''r'$20%s$: ${\epsilon}*{N_i}$ for $Model$-$2$ with $2$-$Fit Parameter$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r' ${\epsilon}*{N_i}${(10$^{-9}$)}') 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_Eps_2parm.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_Eps_2parm.png'.format(str(year),str(year)))
plt.close('all')

##________________ DELTA ____________________##

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_k2,color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \kappa$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \kappa$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_delta_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_delta_k.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_k2, color='green',marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \kappa$' % str(year))
ax.set_ylim([-0.101, 0.101])
ax.set_yticks(np.linspace(-0.1, 0.1, 11))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \kappa$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_delta_k_zoom.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_delta_k_zoom.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_B2, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \rho_\ast$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta \rho_\ast$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_delta_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_delta_Rho.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_eps_ni2*1e9, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta {\epsilon}*{N_i}$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(''r'$\Delta {\epsilon}*{N_i}${(10$^{-9}$)}' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_delta_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_delta_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_R2,color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta {R_{adj}}^2$' % str(year))
ax.set_xlabel(' FillNumber')
#ax.set_ylabel(''r'$\Delta {R_{adj}}$^{2}$${({\%})}' ) 
ax.set_ylabel(r'$\Delta {R_{adj}}^2{({\%})}$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_delta_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_delta_R.png'.format(str(year),str(year)))

plt.close('all')

##________________ RELATIVE DIFFERENCES ____________________##


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_eps_ni2/Eps_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $(\Delta {\epsilon}*{N_i})/({\epsilon}*{N_i})$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$(\Delta {\epsilon}*{N_i})/({\epsilon}*{N_i})$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_delta_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_delta_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_B2/B_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \rho_\ast / \rho_\ast$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$\Delta \rho_\ast / \rho_\ast$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_delta_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_delta_Rho.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_k2/k_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta \kappa / \kappa$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(''r'$\Delta \kappa / \kappa $' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_delta_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_delta_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, delta_R2/(R_list*100), color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: $\Delta {R_{adj}}^2 / {R_{adj}}^2$' % str(year))
ax.set_xlabel('FillNumber')
ax.set_ylabel(r'$\Delta {R_{adj}}^2 / {R_{adj}}^2$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_delta_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_delta_R.png'.format(str(year),str(year)))

plt.close('all')




fig, ax = plt.subplots() 
ax.scatter(fill_list, Eps_list_scan2/Eps_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${({\epsilon}*{N_i})_{2 param}}/{({\epsilon}*{N_i})_{3 param}}$' % str(year))
ax.set_xlabel('FillNumber' )
#ax.set_ylabel(''r'$d {\epsilon}*{N_i}$' ) 
ax.set_ylabel(''r'${({\epsilon}*{N_i})_{2 param}}/{({\epsilon}*{N_i})_{3 param}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_Eps.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_Eps.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, B_list_scan2/B_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${{\rho_\ast}_{2 param}}/{{\rho_\ast}_{3 param}}$' % str(year))
ax.set_xlabel('FillNumber')
#ax.set_ylabel(''r'$d \rho_\ast$' ) 
ax.set_ylabel(''r'${{\rho_\ast}_{2 param}}/{{\rho_\ast}_{3 param}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_Rho.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_Rho.png'.format(str(year),str(year)))
plt.close('all')

fig, ax = plt.subplots() 
ax.scatter(fill_list, k_list_scan2/k_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${{\kappa}_{2 param}}/{{\kappa}_{3 param}}$' % str(year))
ax.set_xlabel('FillNumber')
#ax.set_ylabel(''r'$d \kappa$' )
ax.set_ylabel(''r'${{\kappa}_{2 param}}/{{\kappa}_{3 param}}$' ) 
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_k.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}_2param_relative_diff_k.png'.format(str(year),str(year)))
plt.close('all')


fig, ax = plt.subplots() 
ax.scatter(fill_list, R_list_scan2/R_list, color='green', marker='.', s=50, label='delta_k')
ax.set_title(r'$20%s$: ${{{R_{adj}}^2}_{2 param}}/{{{R_{adj}}^2}_{3 param}}$' % str(year))
ax.set_xlabel(' FillNumber')
ax.set_ylabel(r'${{{R_{adj}}^2}_{2 param}}/{{{R_{adj}}^2}_{3 param}}$')
#plt.legend(loc='best')
plt.savefig('20{}/best_param/{}_2param_relative_diff_R.pdf'.format(str(year),str(year)))
plt.savefig('20{}/best_param/img/{}2param_relative_diff_R.png'.format(str(year),str(year)))

plt.close('all')



        
#defining the stop time of the program      
stop=t.time()
#evaluating the run time 
runtime_seq=stop-start
print('The runtime is:', runtime_seq, '[s]')      

        
