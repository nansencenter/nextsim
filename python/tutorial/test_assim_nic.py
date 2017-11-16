# test_assim_nic.py
# - to see what FiniteElement::assimilate_topazForecastAmsr2OsisafNicIce does
import numpy as np
M_nic_conc        = np.linspace(0,.9,3)
M_nic_weekly_conc = np.array([0.,.05,.1,.25,.5,.7,.9,.95,1.])
TEST_DAILY        = 1
TEST_WEEKLY       = 1

if TEST_DAILY:
   print('\nprinting conc brackets for NIC daily')
   for i,ci in enumerate(M_nic_conc):
      if(M_nic_conc[i]<=0.45):
         # CT18
         # .1 - .8
         alpha_up=(0.45-M_nic_conc[i])/0.45;
         thin_conc_obs_min=0.*alpha_up+0.1*(1-alpha_up);
         thin_conc_obs_max=0.1*alpha_up+0.8*(1-alpha_up);
      elif(M_nic_conc[i]<=0.9):
         # CT81
         # .8 - 1.
         # alpha_up=(0.9-M_nic_conc[i])/0.45; # original code
         alpha_up=(0.9-M_nic_conc[i])/0.9;
         thin_conc_obs_min=0.1*alpha_up+0.8*(1-alpha_up);
         thin_conc_obs_max=0.8*alpha_up+1.0*(1-alpha_up);
      elif(M_nic_conc[i]<=1.):
         thin_conc_obs_min=0.8;
         thin_conc_obs_max=1.;
      else:
         # should not happen
         thin_conc_obs_min=1.;
         thin_conc_obs_max=1.;

      print(ci,thin_conc_obs_min,thin_conc_obs_max)

if TEST_WEEKLY:
   print('\nprinting conc brackets for NIC weekly')
   for i,ci in enumerate(M_nic_weekly_conc):
      if(M_nic_weekly_conc[i]<=0.05):
          # CT01
          # 0 - .1
          alpha_up=(0.05-M_nic_weekly_conc[i])/(0.05-0.);
          thin_conc_obs_min=0.0*alpha_up+0.0*(1-alpha_up);
          thin_conc_obs_max=0.0*alpha_up+0.1*(1-alpha_up);
      elif(M_nic_weekly_conc[i]<=0.10):
          # CT02
          # 0 - .2
          alpha_up=(0.10-M_nic_weekly_conc[i])/(0.10-0.);
          thin_conc_obs_min=0.0*alpha_up+0.0*(1-alpha_up);
          thin_conc_obs_max=0.0*alpha_up+0.2*(1-alpha_up);
      elif(M_nic_weekly_conc[i]<=0.25):
          # CT14, CT13, CT24
          # 0.1 - .4
          alpha_up=(0.25-M_nic_weekly_conc[i])/(0.25-0.);
          thin_conc_obs_min=0.0*alpha_up+0.1*(1-alpha_up);
          thin_conc_obs_max=0.1*alpha_up+0.4*(1-alpha_up);
      elif(M_nic_weekly_conc[i]<=0.50):
          # CT46
          # 0.4 - .6
          alpha_up=(0.50-M_nic_weekly_conc[i])/(0.50-0.25);
          thin_conc_obs_min=0.1*alpha_up+0.40*(1-alpha_up);
          thin_conc_obs_max=0.4*alpha_up+0.60*(1-alpha_up);
      elif(M_nic_weekly_conc[i]<=0.70):
          # T68
          # 0.6 - .8
          alpha_up=(0.70-M_nic_weekly_conc[i])/(0.70-0.50);
          thin_conc_obs_min=0.40*alpha_up+0.60*(1-alpha_up);
          thin_conc_obs_max=0.60*alpha_up+0.80*(1-alpha_up);
      elif(M_nic_weekly_conc[i]<=0.90):
          # CT81
          # .8 - 1.
          alpha_up=(0.90-M_nic_weekly_conc[i])/(0.90-0.70);
          thin_conc_obs_min=0.60*alpha_up+0.80*(1-alpha_up);
          thin_conc_obs_max=0.80*alpha_up+1.0*(1-alpha_up);
      elif(M_nic_weekly_conc[i]<=1.):
          # CT92
          # .9 - 1.
          alpha_up=(1.0-M_nic_weekly_conc[i])/(1.0-0.9);
          thin_conc_obs_min=0.8*alpha_up+1.0*(1-alpha_up);
          thin_conc_obs_max=1.0*alpha_up+1.0*(1-alpha_up);
      else:
          # should not happen
          thin_conc_obs_min=1.;
          thin_conc_obs_max=1.;
   
      print(ci,thin_conc_obs_min,thin_conc_obs_max)
