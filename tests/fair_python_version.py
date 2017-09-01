# Original Python implementation of FAIR.
# From Richard Millar: January 23, 2017

import numpy as np
from scipy.optimize import root
from scipy.ndimage.filters import gaussian_filter1d

def iirf_interp_funct(alp_b,a,b,t_iirf,iirf,x,e_dims,p_dims):

    alp_b = alp_b.reshape((e_dims,p_dims))
    iirf_arr = alp_b*(np.sum(a*b*(1.0 - np.exp(-t_iirf[...,np.newaxis]/(b*alp_b[...,np.newaxis]))),axis=-1))
    return np.ndarray.flatten(iirf_arr   -  iirf[...,x])

def fair_scm(in_driver,other_rf=0.0,input_params=np.array([1.75,2.5,4.1,239.0,0.2173,0.2240,0.2824,0.2763,1000000,394.4,36.54,4.304,100.0,35.0,0.02,4.5,3.74,278.0,2.123,97.0]),restart_in=False,restart_out=False,mode='emissions_driven'):

    #Broadcast everything to the same shape
    #Check that all parameter shapes are the same or 1
    if in_driver.ndim==1:
      e_dims = 1
      in_driver = in_driver.reshape((1,in_driver.shape[0]))
    else:
      e_dims = in_driver.shape[0]
    t_dims = in_driver.shape[-1]
    if input_params.ndim==1:
      p_dims = 1
      input_params = input_params.reshape((1,input_params.shape[0]))
    else:
      p_dims = input_params.shape[0]



    in_params = input_params[np.newaxis,...]

    ecstcr = in_params[...,0:2]
    d = in_params[...,2:4]
    a = in_params[...,4:8]
    b = in_params[...,8:12]
    t_iirf = in_params[...,12]
    r0 = in_params[...,13]
    rc = in_params[...,14]
    rt = in_params[...,15]
    F_2x = in_params[...,16]
    pre_indust_co2 = in_params[...,17]
    c = in_params[...,18]
    iirf_max = in_params[...,19]



    #specific TCR and ECS
    k = 1.0 - (d/70.0)*(1.0 - np.exp(-70.0/d))
    q = (np.transpose((1.0 / F_2x) * (1.0/(k[...,0]-k[...,1])) * np.array([ecstcr[...,0]-ecstcr[...,1]*k[...,1],ecstcr[...,1]*k[...,0]-ecstcr[...,0]]))).reshape((1,p_dims,2))


    #Form the output timeseries variables
    radiative_forcing = np.zeros(shape=(e_dims,p_dims,t_dims))
    cumuptake = np.zeros(shape=(e_dims,p_dims,t_dims))
    iirf = np.zeros(shape=(e_dims,p_dims,t_dims))
    R = np.zeros(shape=(e_dims,p_dims,t_dims,b.shape[-1]))
    T = np.zeros(shape=(e_dims,p_dims,t_dims,q.shape[-1]))

    #Reshape the emissions arrays
    if type(other_rf) != float:
      if other_rf.ndim == 1:
        other_rf = other_rf[np.newaxis,...]
      other_rf = other_rf[:,np.newaxis,...]
    if mode == 'emissions_driven':
      emissions = in_driver[:,np.newaxis,...]
      out_concs = np.zeros(shape=(e_dims,p_dims,t_dims))
      out_temp = np.zeros(shape=(e_dims,p_dims,t_dims))
    if mode == 'emissions_back':
      out_concs = in_driver[:,np.newaxis,...] - pre_indust_co2[...,np.newaxis]
      emissions = np.zeros(shape=(e_dims,p_dims,t_dims))
      out_temp = np.zeros(shape=(e_dims,p_dims,t_dims))
    if mode == 'forcing_driven':
      radiative_forcing = in_driver[:,np.newaxis,...]
      out_temp = np.zeros(shape=(e_dims,p_dims,t_dims))
      out_concs = np.zeros(shape=(e_dims,p_dims,t_dims))
    if mode == 'forcing_back':
      out_temp = in_driver[:,np.newaxis,...]
      radiative_forcing = np.zeros(shape=(e_dims,p_dims,t_dims))
      out_concs = np.zeros(shape=(e_dims,p_dims,t_dims))


    if restart_in:
      R[...,0,:]=restart_in[0]
      T[...,0,:]=restart_in[1]
      cumuptake[...,0] = restart_in[2]
      if mode != 'forcing_back':
        out_temp[...,0] = np.sum(T[...,0,:],axis=-1)
      if mode != 'emissions_back':
        out_concs[...,0] = np.sum(R[...,0,:],axis=-1)
      if mode != 'emissions_driven':
        emissions[...,0] = restart_in[3]
      if mode != 'forcing_driven':
        if type(other_rf) == float:
          radiative_forcing[...,0] = (F_2x/np.log(2.)) * np.log((out_concs[...,0] + pre_indust_co2) /pre_indust_co2) + other_rf
        else:
          radiative_forcing[...,0] = (F_2x/np.log(2.)) * np.log((out_concs[...,0] + pre_indust_co2) /pre_indust_co2) + other_rf[...,0]

    elif (mode in ['emissions_driven','emissions_back']):

      #Initialise the carbon pools to be correct for first timestep in numerical method
      R[...,0,:] = a * emissions[...,0,np.newaxis] / (c[...,np.newaxis]) * 0.5
      cumuptake[...,0] = emissions[...,0]



    for x in range(1,in_driver.shape[-1]):

      if (mode in ['emissions_driven','emissions_back']):
        #Calculate the parametrised iIRF and check if it is over the maximum allowed value
        iirf[...,x] = (rc * cumuptake[...,x-1]  + rt* out_temp[...,x-1]  + r0   )
        iirf[(iirf[...,x]-iirf_max)>0,x] = np.broadcast_to(iirf_max,(e_dims,p_dims))[(iirf[...,x]-iirf_max)>0]

        #Linearly interpolate a solution for alpha
        if x == 1:
          time_scale_sf = (root(iirf_interp_funct,0.16*np.ones((e_dims,p_dims)),args=(a,b,t_iirf,iirf,x,e_dims,p_dims)))['x'].reshape(e_dims,p_dims)
        else:
          time_scale_sf = (root(iirf_interp_funct,time_scale_sf,args=(a,b,t_iirf,iirf,x,e_dims,p_dims)))['x'].reshape(e_dims,p_dims)


        #Multiply default timescales by scale factor
        b_new = b * time_scale_sf[...,np.newaxis]

      if mode == 'emissions_driven':
        #Do forward numerical method, carbon first followed by forcing followed by temperature
        R[...,x,:] = R[...,x-1,:]*(np.exp((-1.0)/b_new)) + 0.5*(a)*(emissions[...,x-1,np.newaxis]+emissions[...,x,np.newaxis])/ c[...,np.newaxis]
        out_concs[...,x] = np.sum(R[...,x,:],axis=-1)
        cumuptake[...,x] =  cumuptake[...,x-1] + emissions[...,x] - (out_concs[...,x] - out_concs[...,x-1])*c

      if mode == 'emissions_back':
        R[...,x,:] = R[...,x-1,:]*(np.exp((-1.0)/b_new))
        emissions[...,x]  = (c*(out_concs[...,x] - np.sum(R[...,x,:],axis=-1)) - 0.5*np.sum(a,axis=-1)*emissions[...,x-1]) / (0.5*np.sum(a,axis=-1))
        R[...,x,:] =  R[...,x,:] + 0.5*(a)*(emissions[...,x-1,np.newaxis]+emissions[...,x,np.newaxis])/ c[...,np.newaxis]
        cumuptake[...,x] =  cumuptake[...,x-1] + emissions[...,x] - (out_concs[...,x] - out_concs[...,x-1])*c

      if (mode in ['emissions_driven','emissions_back']):
        if type(other_rf) == float:

          #NB/ CALC OTHER RFS COULD BE DONE IN HERE INTERACTIVLY BY CALLING BOX MODEL FUNCTIONS
          radiative_forcing[...,x] = (F_2x/np.log(2.)) * np.log((out_concs[...,x] + pre_indust_co2) /pre_indust_co2) + other_rf
        else:
          radiative_forcing[...,x] = (F_2x/np.log(2.)) * np.log((out_concs[...,x] + pre_indust_co2) /pre_indust_co2) + other_rf[...,x]

      if (mode in ['emissions_driven','emissions_back','forcing_driven']):
        T[...,x,:] = T[...,x-1,:]*(np.exp((-1.0)/d)) + 0.5*(q)*(radiative_forcing[...,x-1,np.newaxis]+radiative_forcing[...,x,np.newaxis])*(1-np.exp((-1.0)/d))
        out_temp[...,x]=np.sum(T[...,x,:],axis=-1)

      if mode == 'forcing_back':
        T[...,x,:] = T[...,x-1,:]*(np.exp((-1.0)/d))
        radiative_forcing[...,x] = (out_temp[...,x] - np.sum(T[...,x,:],axis=-1) - 0.5*(radiative_forcing[...,x-1])*np.sum(q*(1-np.exp((-1.0)/d)),axis=-1) ) / (0.5*np.sum((q*(1-np.exp((-1.0)/d)) ),axis=-1) )
        T[...,x,:] = T[...,x,:] + 0.5*(q)*(radiative_forcing[...,x-1,np.newaxis]+radiative_forcing[...,x,np.newaxis])*(1-np.exp((-1.0)/d))

    return_concs = out_concs + pre_indust_co2[...,np.newaxis]
    return_temps = out_temp
    if mode in ['emissions_driven','emissions_back']:
      return_emissions = emissions
    return_forcing = radiative_forcing
    if e_dims==1 and p_dims!=1:
      return_concs = return_concs[0]
      return_temps = return_temps[0]
      if mode in ['emissions_driven','emissions_back']:
        return_emissions = emissions[0]
      return_forcing = radiative_forcing[0]
    if e_dims!=1 and p_dims==1:
      return_concs = return_concs[:,0]
      return_temps = return_temps[:,0]
      if mode in ['emissions_driven','emissions_back']:
        return_emissions = emissions[:,0]
      return_forcing = radiative_forcing[:,0]
    if e_dims==1 and p_dims==1:
      return_concs = return_concs[0,0]
      return_temps = return_temps[0,0]
      if mode in ['emissions_driven','emissions_back']:
        return_emissions = emissions[0,0]
      return_forcing = radiative_forcing[0,0]

    #Get the required output
    if mode == 'emissions_driven':
      if restart_out:
        restart_out_val=(R[...,-1,:],T[...,-1,:],cumuptake[...,-1])
        return [return_concs, return_temps], restart_out_val
      else:
        return [return_concs, return_temps]
    if mode == 'emissions_back':
      #Smooth the backed out emissions
      emissions = gaussian_filter1d(emissions,sigma=2,axis=-1)
      if restart_out:
        restart_out_val=(R[...,-1,:],T[...,-1,:],cumuptake[...,-1])
        return [return_emissions, return_temps], restart_out_val
      else:
        return [return_emissions, return_temps]
    if mode == 'forcing_driven':
      if restart_out:
        restart_out_val=(R[...,-1,:],T[...,-1,:],cumuptake[...,-1])
        return return_temps, restart_out_val
      else:
        return return_temps
    if mode == 'forcing_back':
      #Smooth the backed out emissions
      radiative_forcing = gaussian_filter1d(radiative_forcing,sigma=2,axis=-1)
      if restart_out:
        restart_out_val=(R[...,-1,:],T[...,-1,:],cumuptake[...,-1],emissions[...,-1])
        return return_forcing, restart_out_val
      else:
        return return_forcing
