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

year = 16

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
    
    with open('20{}/Coeff/{}_2param_first_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        k=[] #fill number
        eps_ni=[]
        B=[]
        R_squared=[]
        chi_squared=[]

        


        for x in lines:
            k.append(float(x.split(' ')[0]))
            eps_ni.append(float(x.split(' ')[1]))
            B.append(float(x.split(' ')[2])) 
            R_squared.append(float(x.split(' ')[3]))
            chi_squared.append(float(x.split(' ')[4])) 

          
    k_list = np.array(k)
    Eps_list = np.array(eps_ni)
    B_list = np.array(B)
    R_squared_list = np.array(R_squared)
    chi_squared_list = np.array(chi_squared)

    with open('20{}/Coeff/{}_2param_last_FitCoeff.txt'.format((str(year)),text), 'r') as f:
        lines=f.readlines()
        k_last=[] #fill number
        eps_ni_last=[]
        B_last=[]
        R_squared_last=[]
        chi_squared_last=[]

        


        for x in lines:
            k_last.append(float(x.split(' ')[0]))
            eps_ni_last.append(float(x.split(' ')[1]))
            B_last.append(float(x.split(' ')[2])) 
            R_squared_last.append(float(x.split(' ')[3]))
            chi_squared_last.append(float(x.split(' ')[4])) 

          
    k_list_last = np.array(k_last)
    Eps_list_last = np.array(eps_ni_last)
    B_list_last = np.array(B_last)
    R_squared_list_last = np.array(R_squared_last)
    chi_squared_list_last = np.array(chi_squared_last)
    
           
            

    fig, ax = plt.subplots()
    ax.plot(k_list, chi_squared_list, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))
    #ax.set_ylim([-0.4, 10])
    #ax.set_yticks(np.linspace(0.5, 1.5, 16))           
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel('${\chi}^2$') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/quality/{}_k_chi_full.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/quality/img/{}_k_chi_full.png'.format(str(year),text))

    fig, ax = plt.subplots()
    ax.plot(k_list, chi_squared_list, "b.", markersize=4, label='k')  
    #ax.set_xlim([0., 2.1])
    #ax.set_xticks(np.linspace(0.1, 2, 20))
    ax.set_ylim([0., 1])
    #ax.set_yticks(np.linspace(0., 1, 16))           
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel('${\chi}^2$') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/quality/{}_CHI_k_full_zoom.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/quality/img/{}_CHI_k_full_zoom.png'.format(str(year),text))


    fig, ax = plt.subplots()
    ax.plot(k_list, R_squared_list*100, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))          
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel('${R_{adj}}$$^2${({\%})}' )
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/quality/{}_k_R_full.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/quality/img/{}_k_R_full.png'.format(str(year),text))

    fig, ax = plt.subplots()
    ax.plot(k_list, B_list, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))         
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel(r'$\rho_\ast$') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/param/{}_k_Rho_full.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/param/img/{}_k_Rho_full.png'.format(str(year),text))

    fig, ax = plt.subplots()
    ax.plot(k_list, Eps_list*1e9, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))         
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel(''r' ${\epsilon}*{N_i}${(10$^{-9}$)}') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/param/{}_k_Eps_full.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/param/img/{}_k_Eps_full.png'.format(str(year),text))


        
        
        


    fig, ax = plt.subplots()
    ax.plot(k_list_last, chi_squared_list_last, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))
    #ax.set_ylim([-0.4, 10])
    #ax.set_yticks(np.linspace(0.5, 1.5, 16))           
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel('${\chi}^2$') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/quality/{}_k_chi_last.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/param/img/{}_k_chi_last.png'.format(str(year),text))

    fig, ax = plt.subplots()
    ax.plot(k_list_last, chi_squared_list_last, "b.", markersize=4, label='k')  
    #ax.set_xlim([0., 2.1])
    #ax.set_xticks(np.linspace(0.1, 2, 20))
    ax.set_ylim([0., 0.1])
    #ax.set_yticks(np.linspace(0., 1, 16))           
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel('${\chi}^2$') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/quality/{}_CHI_k_last_zoom.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/param/img/{}_CHI_k_last_zoom.png'.format(str(year),text))


    fig, ax = plt.subplots()
    ax.plot(k_list_last, R_squared_list_last*100, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))          
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel('${R_{adj}}$$^2${({\%})}' )
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/quality/{}_k_R_last.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/quality/img/{}_k_R_last.png'.format(str(year), text))

    fig, ax = plt.subplots()
    ax.plot(k_list_last, B_list_last, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))         
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel(r'$\rho_\ast$') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/param/{}_k_Rho_last.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/param/img/{}_k_Rho_last.png'.format(str(year),text))


    fig, ax = plt.subplots()
    ax.plot(k_list_last, Eps_list_last*1e9, "b.", markersize=4, label='k')  
    #ax.set_xlim([0.4, 1.6])
    #ax.set_xticks(np.linspace(0.5, 1.5, 16))         
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel(''r' ${\epsilon}*{N_i}${(10$^{-9}$)}') 
    ax.set_title('FillNumber: {}'.format(text))
    plt.savefig('20{}/scan/param/{}_k_Eps_last.pdf'.format(str(year),text))
    plt.savefig('20{}/scan/param/img/{}_k_Eps_last.png'.format(str(year),text))




    

        
#defining the stop time of the program      
stop=t.time()
#evaluating the run time 
runtime_seq=stop-start
print('The runtime is:', runtime_seq, '[s]')      

        
