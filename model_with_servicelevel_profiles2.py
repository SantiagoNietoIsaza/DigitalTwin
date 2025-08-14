#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:07:54 2021

@author: santi
"""

import math
import random
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import simpy
import pandas as pd
from scipy import stats
from scipy.stats import rv_discrete
import os
import csv
import statistics
from collections import defaultdict
import time

# simulation:

class ParcelArrival(object):
    "this class represents arrivals of parcels"
    def __init__(self,env,key,old_key):
        self.od=key
        self.old_od=old_key
        self.env=env
        self.store=simpy.Store(env)
        self.candidates=None
        self.num_routings=0
    
    def put(self,timewindow):
        global box_caps
        global BackupCost_cop
        global BackupCost
        if timewindow is None: # new parcel
            self.candidates = [i for i in factor_columns[0] if i != self.od[0]] # saves the possible SPs for transitioning
            timewindow =(self.env.now, self.env.now+service_level, Distances[self.od], self.candidates, self.od[0], self.env.now, self.num_routings) # time window tuple from arrival time at origin to latest delivery time at destination, now I include the od distance
            parcels[self.od[0], self.od[1], self.env.now, self.env.now+service_level] = 0
            if origin_wait_nocap:
                Parcel_arrivals[self.od].append(timewindow)
                Parcel_arrivals_cop[self.od].append(timewindow) #copy of parcel arrivals
                #box_caps[self.od[0]] -= 1  # Revisar si tengo que quitar esto para modelar no uso de capacidad en boxes en el origen (creo que no)
            else:
                if capacitated and box_caps[self.od[0]] <= 0:
                    BackupCost_cop += pc
                    BackupCost += pc
                    parcel_key = (self.od[0], self.od[1], self.env.now, self.env.now+service_level)
                    if parcel_key not in parcel_legs:
                        parcel_legs[parcel_key] = [0,1] # start counter for number of legs
                    else:
                        #parcel_legs[parcel_key][0] += 1 # updates counter for number of legs
                        parcel_legs[parcel_key][1] += 1
                    if compute_shadows:
                        if policy == 3:
                            if discrete_MDP:
                                transport_time = math.ceil((Distances[self.od]/avg_speed) * tau_p)
                                #tau_1 = int(env.now / tau_p)
                                #tau_2 = tau_1 + int(service_level / tau_p) - transport_time
                                nper = int(service_level/tau_p)
                                deadline = timewindow[1]  # Deadline in absolute time (clock time)
                                t_left = deadline - env.now
                                t_elapsed = service_level - t_left
                                tau_2 = int((t_elapsed + transport_time) / tau_p)
                                tau_3 = int(t_elapsed / tau_p)
                                curr_state = (self.od[0],tau_3,self.od[1],nper)
                            else:
                                transport_time = (((Distances[self.od]/avg_speed)/3600)/tau_p) * tau_p
                                tau_1 = round(int(env.now / tau_p) * tau_p, 4)
                                tau_2 = round(tau_1 + int(service_level/tau_p) * tau_p - transport_time, 4)
                                curr_state = (self.od[0],self.od[1],tau_1,tau_2)
                            shadow_prices_comp[self.od[0],tau_3] += pc - V_st[curr_state] # chequear esto, si es continuo deberían estar por unidad de tiempo
                            shadow_counts[self.od[0],tau_3] += 1
                        elif policy == 4:
                            if discrete_MDP:
                                tau_1 = int(env.now / tau_p)
                            else:
                                tau_1 = round(int(env.now / tau_p) * tau_p, 4)
                            transport_time = ((Distances[self.od]/avg_speed)/3600)
                            remaining_time = service_level - transport_time
                            #shadow_prices[self.od[0],tau_1] += stats.expon.cdf(remaining_time, loc=0, scale=OC_lambdas[self.od]) * (pc -r1) # Revise
                            shadow_prices_comp[self.od[0],tau_1] += stats.expon.cdf(remaining_time, loc=0, scale=OC_lambdas[self.od]*profiles[int(env.now)+1]) * (pc + r1 * my_distances[self.od])
                            shadow_counts[self.od[0],tau_1] += 1
                else:
                    Parcel_arrivals[self.od].append(timewindow)
                    Parcel_arrivals_cop[self.od].append(timewindow) #copy of parcel arrivals
                    box_caps[self.od[0]] -= 1
                    if compute_shadows:
                        if discrete_MDP:
                            tau_1 = int(env.now / tau_p)
                        else:
                            tau_1 = round(int(env.now / tau_p) * tau_p, 4)
                        shadow_counts[self.od[0],tau_1] += 1
        else: # transferred parcel
            Parcel_arrivals[self.od].append(timewindow)
            Parcel_arrivals[self.od].sort(key=lambda tup: tup[0]) # sort list of parcels in ascending arrival order
            Parcel_arrivals_cop[self.od].append(timewindow)
            Parcel_arrivals_cop[self.od].sort(key=lambda tup: tup[0]) # sort list of parcels in ascending arrival order in the copy
            if origin_wait_nocap: # Check is infinite capacity is allowed at origin locations
                if self.old_od[0] != timewindow[4]:  # Only increment capacity if the transfer os not from the parcel's origin location
                    box_caps[self.old_od[0]] += 1
                else:
                    pass
            else:
                box_caps[self.old_od[0]] += 1
            #box_caps[self.od[0]] += 1
            box_caps[self.od[0]] -= 1
    
    def pop(self):
        # delete assigned parcel (= first in list) from list
        Parcel_arrivals[self.od].pop(0)
        Parcel_arrivals_cop[self.od].pop(0) # delete parcel in the copy

    
    
    
    
def arrival_parcels(env,parcel):
    "simulates new parcel arrivals"
    global totalParcels
    global box_caps
    while env.now <timeframe_P: 
        #wait for next arrival
        yield env.timeout(random.expovariate(Parcel_lambdas[parcel.od]))
        if env.now <=timeframe_P: 
            parcel.put(None)
            totalParcels +=1

def usage_atDestination(env, assigned_parcel, time_to_deadline):
    """
    Adds the estimated usage cost of box at destination.
    
    env: Simulation environment
    time_to_deadline: Time until delivery deadline
    delta: Time length of discrete periods
    """
    global box_caps
    #time_to_deadline = Parcel_arrivals[assigned_parcel.od][0][1] - (env.now+transport_time)
    destination = assigned_parcel.od[1]
    if env.now <=timeframe_P: 
        yield env.timeout(time_to_deadline)
        box_caps[destination] += 1

class CourierArrival(object):
    "this class represents arrivals of couriers"
    def __init__(self,env,key):
        self.od=key
        self.env=env
        self.store=simpy.Store(env)
       

    def Assignment(self, prob_dist):
         global CourierCost1
         global CourierCost2
         global assigned
         global transshipment
         global used_OCs
         global box_caps
         
         Options={} # dict of parcels available and their expected transport time to destination
         Options2={} # dict of parcels available and their expected cost to destination
         transport_time = (Distances[self.od]/avg_speed)/3600 # travel time of OC in hours
         for key in Parcel_arrivals:
             if capacitated and box_caps[self.od[1]] <= 0: # Breaks the cycle if there is no capacity in the destination box
                 # if discrete_MDP:
                 #     tau_1 = int(env.now / tau_p)
                 #     tau_2 = int(Parcel_arrivals[key][i][1] / tau_p)
                 #     curr_state = (self.od[0],self.od[1],tau_1,tau_2)
                 #     wait_state = (self.od[0],self.od[1],tau_1 + 1,tau_2)
                 # else:
                 #     tau_1 = int(env.now / tau_p) * tau_p
                 #     tau_2 = int(Parcel_arrivals[key][i][1] / tau_p) * tau_p
                 #     curr_state = (self.od[0],self.od[1],tau_1,tau_2)
                 #     wait_state = (self.od[0],self.od[1],tau_1 + tau_p,tau_2)
                 # break
                 pass
             if key[0] == self.od[0]: # parcel at same origin location as courier
                 Options[key] = [] # create lists of parcels with the same key that might be assigned
                 Options2[key] = []
                 if Parcel_arrivals[key]: # if list not empty = parcels already arrived at origin
                     for i in range(len(Parcel_arrivals[key])):
                        if Parcel_arrivals[key][i][0] <= env.now: # arrival of parcel in queue
                             exp_time_wait = (Distances[key]/avg_speed)/3600 + (1/(OC_lambdas[key]*profiles[int(env.now+1)])) # expected delivery time if waiting for direct transport to destination
                             if key != self.od:
                                 rest = (self.od[1],key[1])
                                 t_prime = env.now + (Distances[self.od]/avg_speed)/3600
                                 if policy == 0:
                                     exp_remaining_time = (Distances[rest]/avg_speed)/3600 + (1/(OC_lambdas[rest]*profiles[int(t_prime)+1])) # travel time + time until next OC arrives
                                 elif policy == 1:
                                     exp_remaining_time = (Distances[rest]/avg_speed)/3600 # considers only travel time for myopic policy
                                 elif policy == 3:
                                     exp_remaining_time = (Distances[rest]/avg_speed)/3600 + (1/(OC_lambdas[rest]*profiles[int(t_prime)+1])) # I will not use this variable for the optimal policy, so I just assign any value to avoid errors
                                 elif policy == 4:
                                     exp_remaining_time = (Distances[rest]/avg_speed)/3600 + (1/(OC_lambdas[rest]*profiles[int(t_prime)+1]))
                                 remaining_distance = Distances[rest] # computes de remaining distance to final destination
                             else:
                                 exp_remaining_time = 0
                                 remaining_distance = 0
                             exp_time = transport_time + exp_remaining_time 
                             if policy == 0: # executes type of policy time look-ahead
                                 if exp_time < exp_time_wait: # expected time until delivery gets reduced by assignment to OC
                                     if env.now + exp_time <= Parcel_arrivals[key][i][1]: # if exp. arrival time at final destination is within timewindow -> possible assignment
                                         Options[key].append(exp_time) # the parcel is candidate for courier assignment
                             elif policy == 1: # executes policy type distance myopic
                                 if remaining_distance < Parcel_arrivals[key][i][2]: # if the crowd courier takes me closer
                                     if env.now + exp_time <= Parcel_arrivals[key][i][1]: # if exp. arrival time at final destination is within timewindow -> possible assignment
                                         Options[key].append(exp_time) # the parcel is candidate for courier assignment
                             elif policy == 3: # executes optimal policy
                                 if discrete_MDP:
                                     #tau_0 = int(env.now / tau_p)
                                     #tau_1 = int((env.now + transport_time) / tau_p)
                                     #tau_2 = round((tau_1/nper - int(tau_1/nper)) * nper)
                                     #tau_3 = round((tau_0/nper - int(tau_1/nper)) * nper)
                                     deadline = Parcel_arrivals[key][i][1]  # Deadline in absolute time (clock time)
                                     t_left = deadline - env.now
                                     t_elapsed = service_level - t_left
                                     if non_stationary:
                                         nper = int(deadline/tau_p) # aplica para demanda no estacionaria
                                         tau_2 = int((env.now + transport_time) / tau_p) # aplica para el caso con demanda no estacionaria
                                         tau_3 = int(env.now / tau_p) # aplica para el caso con demanda no estacionaria
                                     else:
                                         nper = int(service_level/tau_p) # aplica para demanda estacionaria
                                         tau_2 = int((t_elapsed + transport_time) / tau_p) # aplica para el caso con demanda estacionaria
                                         tau_3 = int(t_elapsed / tau_p) # aplica para el caso con demanda estacionaria
                                     
                                     #curr_state = (self.od[0],self.od[1],tau_1,tau_2)
                                     next_state = (self.od[1], tau_2, key[1], nper)
                                     #wait_state = (self.od[0],self.od[1],tau_1 + 1,tau_2)
                                     wait_state = (self.od[0], tau_3 + 1, key[1], nper)
                                     curr_state = (self.od[0], tau_3, key[1], nper)
                                 else: # Aún no implemento la actualizacion de esto. No lo necesito (deprecated)
                                     tau_1 = round(int(env.now / tau_p) * tau_p, 4)
                                     tau_2 = round(int(Parcel_arrivals[key][i][1] / tau_p) * tau_p, 4)
                                     curr_state = (self.od[0],self.od[1],tau_1,tau_2)
                                     wait_state = (self.od[0],self.od[1],tau_1 + tau_p,tau_2)
                                 #if curr_state not in V_st:
                                 #    print(f"Estado actual: {curr_state}")
                                 #if next_state not in V_st:
                                 #    print(f"Siguiente estado: {next_state}")
                                 #if wait_state not in V_st:
                                 #    print(f"Estado espera: {wait_state}")
                                 if wait_state in V_st:
                                     if next_state in V_st:
                                         if V_st[next_state] + r1 * my_distances[self.od] < V_st[wait_state] + shadow_prices[(self.od[0], tau_3 + 1)]:
                                         #if V_st[next_state] + r1 < V_st[wait_state]:
                                             Options[key].append(exp_time)
                                             Options2[key].append(V_st[next_state])
                             elif policy == 4: # Deprecated
                                 origin_lambda = sum(OC_lambdas[self.od[1],k] for k in Parcel_arrivals[key][i][3] if k!= self.od[1])
                                 prob_wait = 1 - prob_dist.cdf(tau_p, loc=0, scale=origin_lambda)
                                 if env.now + exp_time <= Parcel_arrivals[key][i][1]: # if exp. arrival time at final destination is within timewindow -> possible assignment
                                     if key != self.od:
                                         if discrete_MDP:
                                             tau_1 = int(env.now / tau_p)
                                         else:
                                             tau_1 = round(int(env.now / tau_p) * tau_p, 4) # The time of courier arrival in standard units
                                             tau_2 = {k: round(int((env.now + transport_time + tau_p + (Distances[(self.od[1], k)]/avg_speed)/3600) / tau_p) * tau_p, 4) for k in Parcel_arrivals[key][i][3] if (k != self.od[1] and k != key[1])}
                                         # Revise from here
                                         #time_left = Parcel_arrivals[key][i][1] - (env.now + exp_time) # revisar: creo que restar el tiempo promedio (1/OC_lambdas[rest]) está de más
                                         #exp_wait_direct = prob_wait * prob_dist.cdf(time_left - tau_p, loc=0, scale=OC_lambdas[rest]) * (r1 + shadow_prices[self.od[1],tau_1]) + r1 # (me faltan los shadows)
                                         #exp_wait_prof = prob_wait * (1 - prob_dist.cdf(time_left - tau_p, loc=0, scale=OC_lambdas[rest])) * (pc + shadow_prices[self.od[1],tau_1]) + r1 # (me faltan los shadows que tienen que estar por unidad de tiempo)
                                         times_left = {k: Parcel_arrivals[key][i][1] - (env.now + (transport_time + tau_p + (Distances[(self.od[1], k)]/avg_speed)/3600 + (1/OC_lambdas[(k,key[1])]) + (Distances[(k,key[1])]/avg_speed)/3600)) for k in Parcel_arrivals[key][i][3] if (k != self.od[1] and k != key[1])} # revisar: creo que restar (1/OC_lambdas[(k,key[1])]) está de más Revisado Sep2024, ok
                                         times_left[key[1]] = Parcel_arrivals[key][i][1] - (env.now + (transport_time + tau_p + (Distances[(self.od[1], key[1])]/avg_speed)/3600))
                                         sum_left_lambdas = sum(OC_lambdas[self.od[1],k] for k in times_left if times_left[k] > 0)
                                         #exp_nowait_direct = (1 - prob_wait) * sum((OC_lambdas[self.od[1],k]/sum_left_lambdas)*r1 + (prob_dist.cdf(times_left[k], loc=0, scale=OC_lambdas[k,key[1]])) * r1 for k in times_left if (times_left[k] > 0 and k != key[1])) # revisar esto
                                         #exp_nowait_direct += (1 - prob_wait) * (OC_lambdas[self.od[1],key[1]]/sum_left_lambdas)*r1
                                         #exp_nowait_prof = (1 - prob_wait) * sum((OC_lambdas[self.od[1],k]/sum_left_lambdas)*r1 + (1 - prob_dist.cdf(times_left[k], loc=0, scale=OC_lambdas[k,key[1]])) * pc for k in times_left if (times_left[k] > 0 and k != key[1])) # revisar esto
                                         #Options2[key].append(exp_wait_direct + exp_wait_prof + exp_nowait_direct + exp_nowait_prof)
                                         # till here
                                         TR_dc = Parcel_arrivals[key][i][1] - (env.now + transport_time + tau_p + (Distances[rest]/avg_speed)/3600)
                                         gamma2_p = r1 + shadow_prices[self.od[1],tau_1] + prob_wait * (prob_dist.cdf(TR_dc, loc=0, scale=OC_lambdas[rest]) * r1 + (1- prob_dist.cdf(TR_dc, loc=0, scale=OC_lambdas[rest])) * pc)
                                         
                                         tree_b1 = (1 - prob_wait) * sum((OC_lambdas[self.od[1],k]/sum_left_lambdas) * 
                                                                         prob_dist.cdf(times_left[k], loc=0, scale=OC_lambdas[k,key[1]]) *
                                                                         (2 * r1 + shadow_prices[k,tau_2[k]]) # No es tau_1!!! Corregido.
                                                                         for k in times_left if (times_left[k] > 0 and k != key[1])) # First branches in the CFA tree as specified in the paper
                                         
                                         tree_b2 = (1 - prob_wait) * sum((OC_lambdas[self.od[1],k]/sum_left_lambdas) * 
                                                                         prob_dist.cdf(times_left[k], loc=0, scale=OC_lambdas[k,key[1]]) *
                                                                         (pc + r1 + shadow_prices[k,tau_2[k]]) 
                                                                         for k in times_left if (times_left[k] > 0 and k != key[1])) # Second branches
                                         
                                         tree_b3 = (1 - prob_wait) * (OC_lambdas[self.od[1],key[1]]/sum_left_lambdas) * r1 # Third branch
                                         
                                         gamma1_p = r1 + shadow_prices[self.od[1],tau_1] + tree_b1 + tree_b2 + tree_b3
                                         
                                         Options2[key].append(gamma1_p + gamma2_p)
                                         
                                     else:
                                         Options2[key].append(r1)
                                         
                                     Options[key].append(exp_time)
                                     
                 if not Options[key]: # delete empty lists
                     Options.pop(key)
                     Options2.pop(key)
         assigned_to_OC = 0
         for i in range(crowd_courier_cap): # transport max range parcels
            #if capacitated and box_caps[self.od[1]] <= 0: # Breaks the cycle if there is no capacity in the destination box
            #    break
            if Options:
                if policy >= 3:
                    a = min(Options2, key=Options2.get) # get the key with the minimum expected delivery cost to final destination
                else:
                    a = min(Options, key=Options.get) # get the key with the minimum expected delivery time to final destination
                if a != self.od:
                    transport_time = ((Distances[self.od[0], a[1]]/avg_speed)/3600) # tau_ocdp in paper
                    tau_ocdc = ((Distances[self.od[0], self.od[1]]/avg_speed)/3600) # same name in paper
                    tau_dcdp = ((Distances[self.od[1], a[1]]/avg_speed)/3600) # same name in paper
                    remaining_time = Parcel_arrivals[a][0][1] - (env.now + transport_time) # t_Roc in paper
                    t_Rdc = Parcel_arrivals[a][0][1] - (env.now + tau_ocdc + tau_dcdp)
                    t_prime = int(env.now + tau_ocdc)
                    lambda_ocdp = OC_lambdas[self.od[0], a[1]] * profiles[int(env.now)+1]
                    lambda_dcdp = OC_lambdas[self.od[1], a[1]] * profiles[t_prime+1]
                    exp_togo_waitDirect = pc * (1 - stats.expon.cdf(remaining_time, loc=0, scale=lambda_ocdp)) + r1 * my_distances[self.od[0], a[1]] * (stats.expon.cdf(remaining_time, loc=0, scale=lambda_ocdp))
                    exp_togo_route = r1 * my_distances[self.od] + r1 * my_distances[self.od[1], a[1]] * (stats.expon.cdf(t_Rdc, loc=0, scale=lambda_dcdp)) + pc * (1 - stats.expon.cdf(t_Rdc, loc=0, scale=lambda_dcdp))
                else:
                    transport_time = ((Distances[self.od[0], a[1]]/avg_speed)/3600) # tau_ocdp in paper
                    remaining_time = Parcel_arrivals[a][0][1] - (env.now + transport_time) # t_Roc in paper
                    lambda_ocdp = OC_lambdas[self.od[0], a[1]] * profiles[int(env.now)+1]
                    exp_togo_waitDirect = pc * (1 - stats.expon.cdf(remaining_time, loc=0, scale=lambda_ocdp)) + r1 * my_distances[self.od[0], a[1]] * (stats.expon.cdf(remaining_time, loc=0, scale=lambda_ocdp))
                    exp_togo_route = r1 * my_distances[self.od]
                
                if env.now + Options[a][0] <= Parcel_arrivals[a][0][1]: # check this statement, when policy==3 does not work. I created the options2 dictionary
                    if discrete_MDP: # de pronto esta parte es innecesaria
                        #tau_1 = int(env.now / tau_p)
                        #t_wait = tau_1 + 1
                        deadline = Parcel_arrivals[a][0][1]  # Deadline in absolute time (clock time)
                        t_left = deadline - env.now
                        t_elapsed = service_level - t_left
                        if non_stationary:
                            nper = int(deadline/tau_p) # aplica para demanda no estacionaria
                            tau_2 = int((env.now + transport_time) / tau_p) # aplica para demanda no estacionaria
                            tau_1 = int(env.now / tau_p) # aplica para demanda no estacionaria
                        else:
                            nper = int(service_level/tau_p) # aplica para demanda estacionaria
                            tau_2 = int((t_elapsed + transport_time) / tau_p) # aplica para demanda estacionaria
                            tau_1 = int(t_elapsed / tau_p) # aplica para demanda estacionaria
                        
                    else:
                        tau_1 = round(int(env.now / tau_p) * tau_p, 4)
                        t_wait = tau_1 + tau_p
                    
                    if capacitated and box_caps[self.od[1]] <= 0:
                        shadow_counts[self.od[1],tau_1] += 1
                        #curr_state = (self.od[0], a[1], tau_1, Parcel_arrivals[a][0][1], )
                        #wait_state = (self.od[0], a[1], t_wait, Parcel_arrivals[a][0][1], )
                        if compute_shadows:
                            if policy == 3:
                                #nper = int(service_level/tau_p)
                                #deadline = Parcel_arrivals[key][i][1]  # Deadline in absolute time (clock time)
                                #t_left = deadline - env.now
                                #t_elapsed = service_level - t_left
                                #tau_2 = int((t_elapsed + transport_time) / tau_p)
                                #tau_3 = int(t_elapsed / tau_p)
                                next_state = (self.od[1], tau_2, a[1], nper)
                                wait_state = (self.od[0], tau_1 + 1, a[1], nper)
                                shadow_prices_comp[self.od[1],tau_1] += V_st.get(wait_state, 1.1*pc) - (V_st.get(next_state, 1.1*pc) + r1 * my_distances[self.od])
                                #print((wait_state, next_state), V_st[wait_state] - (V_st[next_state] + r1))
                            elif policy == 4:
                                #transport_time = ((Distances[self.od[0], a[1]]/avg_speed)/3600) # tau_ocdp in paper
                                #tau_ocdc = ((Distances[self.od[0], self.od[1]]/avg_speed)/3600) # same name in paper
                                #tau_dcdp = ((Distances[self.od[1], a[1]]/avg_speed)/3600) # same name in paper
                                #remaining_time = Parcel_arrivals[a][0][1] - (env.now + transport_time) # t_Roc in paper
                                #t_Rdc = Parcel_arrivals[a][0][1] - (env.now + tau_ocdc + tau_dcdp)
                                #shadow_prices[self.od[1],tau_1] += stats.expon.cdf(remaining_time, loc=0, scale=OC_lambdas[self.od[0], a[1]]) * (pc -r1) # revisar esto
                                #exp_togo_waitDirect = pc * (1 - stats.expon.cdf(remaining_time, loc=0, scale=OC_lambdas[self.od[0], a[1]])) + r1 * (stats.expon.cdf(remaining_time, loc=0, scale=OC_lambdas[self.od[0], a[1]]))
                                #exp_togo_route = r1 + r1 * (stats.expon.cdf(t_Rdc, loc=0, scale=OC_lambdas[self.od[1], a[1]])) + pc * (1 - stats.expon.cdf(t_Rdc, loc=0, scale=OC_lambdas[self.od[1], a[1]]))
                                shadow_prices_comp[self.od[1],tau_1] += max(0, exp_togo_waitDirect - exp_togo_route)
                                #shadow_counts[self.od[0],tau_1] += 1 # Esto estaba mal aquí, ya estoy contando al inicio del condicional
                                #Options.pop(a) # Eliminates job from the list. Nota: es esto correcto? eliminar toda la llave?
                        if destination_nowait:  # If parcel do not wait at destination ask if it's destination node to let the cycle finish
                            if self.od[1] != a[1]:
                                Options[a].pop(0) # Corregido. Ojo, esto no deberìa estar indentado así.
                                if policy >= 3:
                                    #Options2.pop(a) # Eliminates job from list2 
                                    Options2[a].pop(0) # Elimina 1 trabajo a la vez
                                if not Options[a]:
                                    Options.pop(a)
                                    if policy >= 3:
                                        Options2.pop(a)
                                continue # Continues to compute shadows for each feasible job
                        else: # if it waits at destiantion then it cannot be routed in any case
                            Options[a].pop(0) # Corregido. Ojo, esto no deberìa estar indentado así.
                            if policy >= 3:
                                #Options2.pop(a) # Eliminates job from list2 
                                Options2[a].pop(0) # Elimina 1 trabajo a la vez
                            if not Options[a]:
                                Options.pop(a)
                                if policy >= 3:
                                    Options2.pop(a)
                            continue # Continues to compute shadows for each feasible job
                    else:
                        shadow_counts[self.od[1],tau_1] += 1
                    # Revisar desde aquí para incluir el la condición de costo esperado positivo
                    curr_num_routings = Parcel_arrivals[a][0][6]
                    if policy != 3:
                        if a != self.od: # Check whether is profitable to route parcels according to expected cost to route and max routings
                            if (exp_togo_waitDirect < exp_togo_route) or curr_num_routings >= max_num_routings:
                                Options[a].pop(0) 
                                if policy >= 3:
                                    Options2[a].pop(0) # Elimina 1 trabajo a la vez
                                if not Options[a]:
                                    Options.pop(a)
                                    if policy >= 3:
                                        Options2.pop(a)
                                continue
                        else: # Revisar si esto es necesario
                            if (exp_togo_waitDirect < exp_togo_route) or curr_num_routings > max_num_routings:
                                Options[a].pop(0) 
                                if policy >= 3:
                                    Options2[a].pop(0) # Elimina 1 trabajo a la vez
                                if not Options[a]:
                                    Options.pop(a)
                                    if policy >= 3:
                                        Options2.pop(a)
                                continue
                    assigned_key = a
                    assigned_parcel = ParcelArrival(env, assigned_key, None)
                    if assigned_key != self.od:
                        new_key = (self.od[1],assigned_key[1])
                        update_parcel = ParcelArrival(env, new_key, assigned_key)
                        curr_candidates = Parcel_arrivals[assigned_key][0][3].copy()
                        try:
                            curr_candidates.pop(curr_candidates.index(self.od[1]))
                        except:
                            #pass # correcto? que significa que trate y no encuentre?
                            #continue  # será más correcto un continue? o tal vez un break?
                            #break Creo que esto es incorrecto, si rompe el ciclo no asigna nada y lo unico que pasa es que seguramente
                            pass  # Creo que esto es más correcto
                        curr_num_routings += 1
                        timewindow = (env.now+transport_time, Parcel_arrivals[assigned_key][0][1], Distances[new_key], curr_candidates, Parcel_arrivals[assigned_key][0][4], Parcel_arrivals[assigned_key][0][5], curr_num_routings)
                        update_parcel.put(timewindow) # put parcel in queue of new key
                        transshipment +=1 # count number of transshipments
                    else: # parcel reaches final destination
                        dtime = (env.now+transport_time-Parcel_arrivals[assigned_key][0][0]) # time from start at origin to final destination
                        Deliverytime[assigned_key].append(dtime)
                        if origin_wait_nocap: # Check is infinite capacity is allowed at origin locations
                            if self.od[0] != Parcel_arrivals[assigned_key][0][4]:  # Only increment capacity if the transfer os not from the parcel's origin location
                                box_caps[self.od[0]] += 1 # fataba esto, tengo que actualizar la capacidad del box de origen cuando ya va a ser entregado. Estoy asumiendo que la recogida en el destino es instantánea (verificar)
                        else:
                            box_caps[self.od[0]] += 1
                        if destination_nowait:
                            pass # not finishied
                        else:
                            time_to_deadline = Parcel_arrivals[assigned_key][0][1] - (env.now+transport_time)
                            box_caps[self.od[1]] -= 1
                            env.process(usage_atDestination(env, assigned_parcel, time_to_deadline))
                    parcel_key = (Parcel_arrivals[assigned_key][0][4], assigned_key[1], Parcel_arrivals[assigned_key][0][5], Parcel_arrivals[assigned_key][0][1])
                    if parcel_key not in parcel_legs:
                        parcel_legs[parcel_key] = [1,0] # start counter for number of legs
                    else:
                        parcel_legs[parcel_key][0] += 1 # updates counter for number of legs
                        #parcel_legs[parcel_key][1] += 1
                    assigned_parcel.pop()
                    CourierCost1 = CourierCost1 + r1 * my_distances[self.od] # update OC costs reward scheme 1 paid per parcel
                    assigned +=1 # count number of assigned parcels
                    assigned_to_OC +=1 
                    Options[assigned_key].pop(0)
                    if policy >= 3:
                        Options2[assigned_key].pop(0)
                    if not Options[assigned_key]:
                        Options.pop(assigned_key)
                        if policy >= 3:
                            Options2.pop(assigned_key)
                    if not Options:
                        break
                else:
                    Options.pop(a)
                    if policy >= 3:
                        Options2.pop(a)
         if assigned_to_OC >=1: # check whether OC was used
             used_OCs +=1
             CourierCost2 = CourierCost2 + r2*Distances[self.od] # update OC costs reward scheme 2 paid once per distance
         Options.clear()
         Options2.clear()
         
         #dispatch_prof(env, self.od)
        

def arrival_couriers(env,courier):
    "simulates courier arrivals"
    global totalCouriers
    global box_caps
    while True: 
        #wait for next arrival
        yield env.timeout(random.expovariate(OC_lambdas[courier.od]*profiles[int(env.now)+1]))
        #yield env.timeout(random.expovariate(OC_lambdas[courier.od]))
        courier.Assignment(prob_dist=stats.expon)
        dispatch_prof(env, courier.od)
        totalCouriers +=1

def dispatch_prof(env, od):
    "dispatches parcels due soon with professional couriers"
    global box_caps
    global BackupCost_cop
    #for key in Parcel_arrivals:
    indexes = []
    if Parcel_arrivals_cop[od]: # if list not empty = parcels already arrived at origin
        for i in range(len(Parcel_arrivals_cop[od])):
            if Parcel_arrivals_cop[od][i][0] <= env.now: # arrival of parcel in queue
                transport_time = (Distances[od]/avg_speed)/3600 # travel time of OC in hours
                if env.now + transport_time >= Parcel_arrivals_cop[od][i][1]:
                    if origin_wait_nocap: # Check is infinite capacity is allowed at origin locations
                        if od[0] != Parcel_arrivals_cop[od][i][4]:  # Only increment capacity if the transfer os not from the parcel's origin location
                            box_caps[od[0]] += 1
                    else:
                        box_caps[od[0]] += 1
                    BackupCost_cop += pc
                    indexes.append(i)
                    parcel_key = (Parcel_arrivals_cop[od][i][4], od[1], Parcel_arrivals_cop[od][i][5], Parcel_arrivals_cop[od][i][1])
                    if parcel_key not in parcel_legs:
                        parcel_legs[parcel_key] = [0,1] # start counter for number of legs
                    else:
                        #parcel_legs[parcel_key][0] += 1 # updates counter for number of legs
                        parcel_legs[parcel_key][1] += 1
                    
        for i in range(len(indexes)):
            Parcel_arrivals_cop[od].pop(indexes[i])
            indexes = [x - 1 if index >= i else x for index, x in enumerate(indexes)]
        
def dynamic_program(f, prob_dist, tau):
    # factor_columns[0]: Service points
    # time_per_day: Horizon
    # pc: professional courier cost
    # r1 = 1.50 # fixed reward 
    # r2 = 0.0005 # reward paid per m
    # service_level: time to delivery
    global states, exp_states
    states = []
    #exp_states = {}
    for o in factor_columns[0]:
        for d in factor_columns[0]:
            #if o != d:
            tau_1 = 0
            #tau = 1/_lambda_avg_OC # Check this, data transformation should be done. I'll do data transformation based on tau and tau_od. 1/_lambda_avg_OC is the minimum value for the decision epoch
            try:
                tau_od = 1/OC_lambdas[(o,d)] # if this depends on the action as well, then I can easily model pricing
                transport_time = (Distances[(o,d)]/avg_speed)/3600
            except:
                tau_od = 0
                transport_time = 0
            #while tau_1 <= int(time_per_day/tau) - service_level/tau: # tengo que actualizar service_level a unidades tau. No era necesario, me devolví a solamente incrementar con tau, luego en la asignación cambio el tiempo a de acuerdo con los incrementos discretos
            while tau_1 <= time_per_day - service_level:
                #tau_2 = tau_1 + tau
                tau_2 = tau_1
                #while tau_2 <= int(time_per_day/tau): # no era necesario llevar a periodos
                while tau_2 <= time_per_day:
                    curr_state = (o,d,tau_1,tau_2)
                    #states.append(curr_state)
                    #exp_states[curr_state] = pc/tau_od
                    #print(curr_state)
                    if o == d:
                        exp_states[curr_state] = 0
                    elif tau_1 + transport_time > tau_2: # no falta el condicional para cuando o=d? Sí, hay que incluirlo. Ya está hecho
                        exp_states[curr_state] = 999999
                    else:
                        exp_states[curr_state] = pc/tau_od
                        eps = 0.0
                        P_transf = 0.0
                        wait_state = (o,d,tau_1+tau, tau_2) # check this seems wrong. Corrected.
                        if wait_state not in exp_states:
                        #    pass
                            exp_states[wait_state] = (pc/tau_od)
                        #probs = [1 - prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)]) for k in factor_columns[0] if k!= o]
                        origin_lambda = sum(OC_lambdas[o,k] for k in factor_columns[0] if k!= o)
                        prob_wait = 1 - prob_dist.cdf(tau, loc=0, scale=origin_lambda)
                        #for k in filter(lambda x: x != o, factor_columns[0]): # ok semms that this was very inefficient
                        for k in factor_columns[0]:
                            if k != o:
                                transport_time1 = int(((Distances[(o,k)]/avg_speed)/3600)/tau) * tau # revisar esto, la funcion ceil lo lleva a un entero. Ya quité la funcion ceil
                                try:
                                    transport_time2 = int(((Distances[(k,d)]/avg_speed)/3600)/tau) * tau # revisar esto, la funcion ceil lo lleva a un entero. Ya quite la funcion ceil
                                except:
                                    transport_time2 = 0
                                eta = tau_1 + transport_time1 + transport_time2
                                if eta <= tau_2:
                                    try:
                                        target_dep_time = tau_1 + int(((Distances[(o,d)]/avg_speed)/3600)/tau) * tau # cambié kd por od, esto ya está bien
                                    except:
                                        target_dep_time = tau_1
                                    target_arr_time = tau_2
                                    target_state = (k,d,target_dep_time, target_arr_time) # chequear target state, particularmente el tiempo actual del estado (Ok chequeado y corregido sí es t_od)
                                    if target_state in exp_states:
                                        if exp_states[target_state] + f/tau_od < exp_states[wait_state]:
                                            #prob_cur_targ = prob_dist.cdf(target_dep_time, loc=0, scale=1/tau) # aqui hay un error chequear el lambda, chequear también el valor de x de la probabilidad
                                            #try:
                                                #prob_cur_targ = prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,d)]) # Lambda od?, no debería ser kd? chequear conjuncion de eventos para definir la probabilidad
                                            #    prob_cur_targ = prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)])
                                            #except:
                                                #prob_cur_targ = 0
                                            #myfilter = list(filter(lambda x: x!= o and x!=k, factor_columns[0])) #este filtro se puede sacar del ciclo?
                                            #probs = [1 - prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,l)]) for l in myfilter]
                                            #prob_cur_targ = prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)])
                                            #prob_cur_targ = math.prod(probs, start=prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)]))
                                            #prob_cur_targ = math.prod(probs) * (prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)]) / (1 - prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)])))
                                            prob_cur_targ = (1 - prob_wait) * (OC_lambdas[(o,k)] / origin_lambda)
                                            eps += (tau/tau_od) * prob_cur_targ * (exp_states[target_state] + f/tau_od) # corrected with data transformation (per decision-epoch time unit)
                                            P_transf += prob_cur_targ
                        #eps += (tau/tau_od) * (1-P_transf) * (exp_states[wait_state] + shadow_prices[(o,tau_1)]) # corrected with data transformation (per decision-epoch time unit)
                        eps += (tau/tau_od) * (prob_wait) * (exp_states[wait_state] + shadow_prices[(o,tau_1)]) # corrected with data transformation (per decision-epoch time unit)
                        exp_states[curr_state] = min(pc/tau_od, eps) # yo cambié la indentación de esta instrucción porque debe estar dentro del ciclo del else. # corrected with data transformation (per decision-epoch time unit)
                    tau_2 += tau
                    tau_2 = round(tau_2,4)
                tau_1 += tau
                tau_1 = round(tau_1,4)
                    

def dynamic_program_disc(f, prob_dist, tau): # Hay que verificar cómo hace la recursividad, esto no está funcionando correctamente (lo mismo para el caso continuo)
    # factor_columns[0]: Service points
    # time_per_day: Horizon
    # pc: professional courier cost
    # r1 = 1.50 # fixed reward 
    # r2 = 0.0005 # reward paid per m
    # service_level: time to delivery
    # number of periods are in minutes
    
    global states, exp_states
    states = []
    #exp_states = {}
    num_per = int(round(time_per_day / tau, 0))
    for o in factor_columns[0]:
        for t in range(num_per + 1):
            for t_bar in range(t, t + int(service_level*60) + 1):
                exp_states[o,o,t,t_bar] = 0
    
    factor_columns_pairs = [(j,k) for j in factor_columns[0] for k in factor_columns[0] if j!=k]
    for pair in factor_columns_pairs:
        o = pair[0]
        d = pair[1]
    #for o in factor_columns[0]:
        #for d in factor_columns[0]:
        try:
            #tau_od = 1/OC_lambdas[(o,d)]
            transport_time = math.ceil((Distances[(o,d)]/avg_speed)/60)
        except:
            #tau = 0
            transport_time = 0
        for t in range(num_per + 1): # periods of 1 minute
        #while tau_1 <= int(time_per_day/tau_od) - service_level:
            #tau_2 = tau_1 + tau_od
            #tau_2 = tau_1
            #while tau_2 <= int(time_per_day/tau_od):
            for t_bar in range(t, t + int(service_level*60) + 1): # aqui hay error, hay que pasar el service level a minutos. Yalo hice al multiplicar por 60
                curr_state = (o,d,t,t_bar)
                #if curr_state == ('Trudering', 'Neuperlach Zentrum', 0, 24): 
                #    print("aqui esta")
                #states.append(curr_state)
                #exp_states[curr_state] = pc
                #print(curr_state)
                if o == d:
                    exp_states[curr_state] = 0
                elif t + transport_time > t_bar:
                    exp_states[curr_state] = 999999
                else:
                    #exp_states[curr_state] = pc
                    eps = 0
                    P_transf = 0
                    wait_state = (o,d,t+1,t_bar)
                    if wait_state not in exp_states:
                        pass
                        #exp_states[wait_state] = 999999
                        #if t + 1 + transport_time > t_bar:
                        #    exp_states[wait_state] = 999999
                        #pass
                        #exp_states[wait_state] = pc
                    #for k in filter(lambda x: x != o, factor_columns[0]): # ok semms that this was very inefficient
                    #probs = [1 - prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)]) for k in factor_columns[0] if k!= o]
                    origin_lambda = sum(OC_lambdas[o,k] for k in factor_columns[0] if k!= o)
                    prob_wait = 1 - prob_dist.cdf(tau, loc=0, scale=origin_lambda)
                    for k in factor_columns[0]:
                        if k != o:
                            transport_time1 = math.ceil((Distances[(o,k)]/avg_speed)/60)
                            try:
                                transport_time2 = math.ceil((Distances[(k,d)]/avg_speed)/60)
                            except:
                                transport_time2 = 0
                            eta = t + transport_time1 + transport_time2
                            if eta <= t_bar:
                                try:
                                    target_dep_time = t + math.ceil((Distances[(o,k)]/avg_speed)/60) # también cambié a od
                                except:
                                    target_dep_time = t
                                target_arr_time = t_bar
                                target_state = (k,d,target_dep_time, target_arr_time) # chequear target state, particularmente el tiempo actual del estado
                                if (target_state in exp_states):
                                #if (target_state in exp_states) and (wait_state in exp_states):
                                    if wait_state in exp_states:
                                        RHS = exp_states[wait_state]
                                    else:
                                        RHS = 999999
                                    if exp_states[target_state] + f < RHS:
                                        #prob_cur_targ = prob_dist.cdf(1/60, loc=0, scale=OC_lambdas[(o,k)]) #cambié el lambda por el lambda correspondiente al par od y el x a 1 minuto (es decir, la probabilidad de que pasa menos de 1 minuto hasta la llegada del courier correspondiente)
                                        #prob_cur_targ = math.prod(probs) * (prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)]) / (1 - prob_dist.cdf(tau, loc=0, scale=OC_lambdas[(o,k)])))
                                        prob_cur_targ = (1 - prob_wait) * (OC_lambdas[(o,k)] / origin_lambda)
                                        eps += prob_cur_targ * (exp_states[target_state] + f)
                                        P_transf += prob_cur_targ
                    #eps += (1-P_transf) * (exp_states[wait_state] + shadow_prices[(o,t)])
                    if wait_state in exp_states:
                        eps += (prob_wait) * (exp_states[wait_state] + shadow_prices[(o,t)])
                    exp_states[curr_state] = min(pc, eps) # yo cambié la indentación de esta instrucción porque debe estar dentro del ciclo del else
                    #if eps > 0:
                    #    eps += (prob_wait) * (exp_states[wait_state] + shadow_prices[(o,t)])
                    #    exp_states[curr_state] = min(pc, eps) # yo cambié la indentación de esta instrucción porque debe estar dentro del ciclo del else
                        #pass
                    #else:
                    #    pass
                        #exp_states[curr_state] = pc
                #tau_2 += tau_od
            #tau_1 += tau_od

# def BackupDelivery(env):
#     global BackupCost
#     global backup
#     while True:
#         for key in Parcel_arrivals:
#             direct_transport = (Distances[key]/avg_speed)/3600
#             for i in range(len(Parcel_arrivals[key])):
#                 if not env.now + direct_transport < Parcel_arrivals[key][0][1]:
#                     backup_parcel = ParcelArrival(env,key)
#                     backup_parcel.pop()
#                     BackupCost = BackupCost + pc
#                     backup +=1
#                     if not Parcel_arrivals[key]:
#                         break

# New functions for solving the stochastic dynamic program
def build_space_time_network(physical_nodes, physical_connections, distances, avg_speed, time_horizon, tau):
    """
    Constructs a space-time network and computes incoming and outgoing nodes.

    Parameters:
        physical_nodes (list): List of physical node identifiers.
        physical_connections (list of tuples): List of (origin, destination) tuples for physical connections.
        distances (dict): Dictionary with keys as (origin, destination) and values as distances (meters).
        avg_speed (float): Average speed (meters/second).
        time_horizon (int): Total length of the time horizon (in seconds).
        tau (int): Length of each time period (in seconds).

    Returns:
        space_time_nodes (list): List of space-time nodes represented as (node, time).
        space_time_arcs (list): List of space-time arcs represented as (origin, time_origin, destination, time_destination).
        incoming_nodes (dict): Dictionary where keys are space-time nodes and values are lists of incoming nodes.
        outgoing_nodes (dict): Dictionary where keys are space-time nodes and values are lists of outgoing nodes.
    """
    # Calculate the number of periods
    num_periods = int(round(time_horizon / tau, 0))

    # Generate space-time nodes (location, time)
    space_time_nodes = [(node, t) for node in physical_nodes for t in range(num_periods + 1)]

    # Initialize space-time arcs and dictionaries for incoming and outgoing nodes
    space_time_arcs = []
    incoming_nodes = {node: [] for node in space_time_nodes}
    outgoing_nodes = {node: [] for node in space_time_nodes}

    # Calculate travel times between physical nodes
    travel_times = {}
    for (o, d), distance in distances.items():
        if distance > 0:  # Only consider valid connections
            travel_time = math.ceil((distance / avg_speed) / tau)  # Travel time in periods
            travel_times[(o, d)] = travel_time

    # Generate space-time arcs
    for (o, d) in physical_connections:
        if distances.get((o, d), 0) > 0:  # Valid physical arc
            travel_time = travel_times[(o, d)]
            for t in range(num_periods + 1):
                t_destination = t + travel_time
                if t_destination <= num_periods:  # Ensure within time horizon
                    arc = (o, t, d, t_destination)
                    space_time_arcs.append(arc)
                    origin_node = (o, t)
                    destination_node = (d, t_destination)
                    outgoing_nodes[origin_node].append(destination_node)
                    incoming_nodes[destination_node].append(origin_node)

    # Add waiting arcs
    for node in physical_nodes:
        for t in range(num_periods):
            arc = (node, t, node, t + 1)
            space_time_arcs.append(arc)
            origin_node = (node, t)
            destination_node = (node, t + 1)
            outgoing_nodes[origin_node].append(destination_node)
            incoming_nodes[destination_node].append(origin_node)

    return space_time_nodes, space_time_arcs, incoming_nodes, outgoing_nodes

def subset_space_time_network(space_time_nodes, space_time_arcs, incoming_nodes, outgoing_nodes, origin, destination, OC_lambdas, pc, r1, mydistances):
    """
    Generates a subset of the space-time network by removing arcs and nodes outside the time horizon determined
    by the origin and destination. Also computes arc weights for transfer, waiting, and P-arcs.

    Parameters:
        space_time_nodes (list): List of space-time nodes.
        space_time_arcs (list): List of space-time arcs.
        incoming_nodes (dict): Dictionary of incoming nodes for each space-time node.
        outgoing_nodes (dict): Dictionary of outgoing nodes for each space-time node.
        origin (tuple): Origin space-time node as (node, time).
        destination (tuple): Destination space-time node as (node, time).
        OC_lambdas (dict): Dictionary with keys as (i', j') and values as average rates (lambdas).
        pc (float): Constant weight for P-arcs.
        r1 (float): Base cost parameter for transfer arcs.

    Returns:
        subset_arcs (list): Subset of space-time arcs.
        p_arcs (list): List of P-arcs connecting to the final destination.
        arc_weights (dict): Dictionary with three keys ('regular', 'waiting', 'prof') and values as dictionaries of weights.
        arc_probabilities (dict): Dictionary with three keys ('regular', 'waiting', 'prof') and values as dictionaries of probabilities.
        updated_incoming_nodes (dict): Updated dictionary of incoming nodes.
        updated_outgoing_nodes (dict): Updated dictionary of outgoing nodes.
    """
    final_destination = destination[0]
    origin_time = origin[1]
    destination_time = destination[1]

    def F(tau, lambda_sum):
        """Computes the probability that the event occurs within the current time period."""
        if lambda_sum == 0:
            return 0
        return 1 - math.exp(-lambda_sum * tau)

    # Filter arcs within the time horizon
    subset_arcs = [
        arc for arc in space_time_arcs
        if origin_time <= arc[1] <= destination_time and origin_time <= arc[3] <= destination_time
    ]

    # Filter P-arcs connecting to the final destination, excluding waiting arcs
    p_arcs = [
        arc for arc in subset_arcs
        if arc[2] == final_destination and arc[0] != arc[2]  # Exclude waiting arcs
    ]

    # Filter nodes and update dictionaries
    valid_nodes = {node for arc in subset_arcs for node in [(arc[0], arc[1]), (arc[2], arc[3])]}  # Nodes in valid arcs

    updated_incoming_nodes = {
        node: [n for n in neighbors if n in valid_nodes]
        for node, neighbors in incoming_nodes.items() if node in valid_nodes
    }
    updated_outgoing_nodes = {
        node: [n for n in neighbors if n in valid_nodes]
        for node, neighbors in outgoing_nodes.items() if node in valid_nodes
    }

    # Initialize weights and probabilities
    arc_weights = {'regular': {}, 'waiting': {}, 'prof': {}}
    arc_probabilities = {'regular': {}, 'waiting': {}, 'prof': {}}

    # Assign weights and probabilities for arcs
    for arc in subset_arcs:
        origin_node = (arc[0], arc[1])
        destination_node = (arc[2], arc[3])
        if arc[0] != arc[2]:  # Regular transfer arc
            outgoing = updated_outgoing_nodes[origin_node]
            #lambda_sum = sum(OC_lambdas.get((arc[0], neighbor[0]), 0) for neighbor in outgoing)
            lambda_sum = _lambda_avg_OC * profiles[int(arc[1]*tau_p)+1]
            if lambda_sum > 0:
                #F_tau = F(1, lambda_sum)  # F(tau) for a single time unit. I think there's a mistake here, it fhould be F(tau_p, lambda_sum)
                F_tau = F(tau_p, lambda_sum)
                p_ij = F_tau * ((OC_lambdas.get((arc[0], arc[2]), 0)*profiles[int(arc[1]*tau_p)+1]) / lambda_sum)
                arc_weights['regular'][arc] = r1 * mydistances[(arc[0], arc[2])]
                arc_probabilities['regular'][arc] = p_ij
            else:
                arc_weights['regular'][arc] = r1 * mydistances[(arc[0], arc[2])]
                arc_probabilities['regular'][arc] = 0
        elif arc[0] == arc[2]:  # Waiting arc
            outgoing = updated_outgoing_nodes[origin_node]
            lambda_sum = sum(OC_lambdas.get((arc[0], neighbor[0])*profiles[int(arc[1]*tau_p)+1], 0) for neighbor in outgoing)
            #F_tau = F(1, lambda_sum) if lambda_sum > 0 else 0
            F_tau = F(tau_p, lambda_sum) if lambda_sum > 0 else 0
            p_ij = 1 - F_tau
            arc_weights['waiting'][arc] = 0
            arc_probabilities['waiting'][arc] = p_ij

    # Assign weights and probabilities for P-arcs
    for arc in p_arcs:
        arc_weights['prof'][arc] = pc
        arc_probabilities['prof'][arc] = 1.0

    return subset_arcs, p_arcs, arc_weights, arc_probabilities, updated_incoming_nodes, updated_outgoing_nodes

def solve_stochastic_dynamic_program(subset_arcs, p_arcs, arc_weights, arc_probabilities, updated_incoming_nodes, updated_outgoing_nodes, destination, pc, shadow_prices, destination_nowait):
    """
    Solves the stochastic dynamic program using backward recursion for the last-mile shipping problem.

    Parameters:
        subset_arcs (list): Subset of space-time arcs.
        p_arcs (list): List of P-arcs connecting to the final destination.
        arc_weights (dict): Dictionary with weights for regular, waiting, and professional arcs.
        arc_probabilities (dict): Dictionary with probabilities for regular, waiting, and professional arcs.
        updated_incoming_nodes (dict): Updated dictionary of incoming nodes.
        updated_outgoing_nodes (dict): Updated dictionary of outgoing nodes.
        destination (tuple): Destination space-time node.

    Returns:
        V (dict): Dictionary of expected cost-to-go for each space-time node.
        policy (dict): Dictionary of optimal actions for each space-time node.
    """
    nodes = list(updated_incoming_nodes.keys())[::-1]  # Reverse order of nodes
    V = {destination: 0.0}  # Initialize value function
    policy = {}  # Initialize policy
    if destination_nowait:
        indicator_var = {node: 0 if node[0]==destination[0] else 1 for node in shadow_prices.keys()} # since I am not using capacity in the origin, Ishould assign also 0 to the indicatior variable in that case
    else:
        for i in nodes:
            if i[0] == destination[0]:
                V[i] = sum(shadow_prices[(i[0],t)] for t in range(i[1],destination[1]+1))
        indicator_var = {node: 1 for node in shadow_prices.keys()}

    # Perform backward recursion
    for i in nodes:
        for j in updated_incoming_nodes[i]:
            value_wait_action = (float('inf'), "wait")
            value_route_prof = (float('inf'), f"route-to-{i}")
            value_route_crowd = (float('inf'), "route-crowd")
            arc = j + i  # Define the arc
            wait_arc = (j[0], j[1], j[0], j[1]+1)
            value_wait_action = (arc_weights['waiting'][wait_arc] + V.get((j[0], j[1]+1), 1.1*pc), "wait")
            if arc in arc_weights['prof']:  # Professional courier action
                value_route_prof = (arc_weights['prof'][arc] + V[i], f"prof-route-to-{i}")
            if arc in arc_weights['regular']:  # Crowd courier action
                value_next_crowd = 0.0
                prob_transf = 0.0
                for k in updated_outgoing_nodes[j]:
                    value_next_node = V.get(k, 1.1*pc) + arc_weights['regular'][arc]
                    if value_next_node < value_wait_action[0]:
                        value_next_crowd += arc_probabilities['regular'].get((j[0], j[1], k[0], k[1]), 0) * value_next_node
                        prob_transf += arc_probabilities['regular'].get((j[0], j[1], k[0], k[1]), 0)
                value_next_crowd += (1 - prob_transf) * (V.get((j[0], j[1]+1), 1.1*pc) + shadow_prices[(j[0], j[1]+1)]*indicator_var[(j[0], j[1]+1)])
                value_route_crowd = (value_next_crowd, "route-crowd")
            # Update the best value and action if this is better (revisar)
            value = min(value_wait_action , value_route_prof, 
                        value_route_crowd, key=lambda x: x[0])
            if j in V:
                if value[0] < V[j]:
                    V[j], policy[j] = value
            else:
                V[j], policy[j] = value

    return V, policy

def runSimulation():
    global backup
    global BackupCost
    global TotalCost1
    global TotalCost2
    global deliveries
    global Avg_DT
    global Deliverytime
    global Avg_DT_OC
    global Avg_DT_total
    global OC_ratio
    global avg_parcels_per_OC
    global avg_compensation1_per_OC
    global avg_compensation2_per_OC

    env.run(until=time_per_day) # simulate 1 whole day 
    #while env.peek() < time_per_day:
    #    print(env.peek())
    #    env.step()
    
    # count the total number of parcels not assigned -> backup delivery
    for key in Parcel_arrivals: # remaining parcels after all allocations are done
        backup = backup + len(Parcel_arrivals[key])
    BackupCost += backup*pc
    
    TotalCost1 = CourierCost1 + BackupCost # costs with reward scheme 1
    TotalCost2 = CourierCost2 + BackupCost # costs with reward scheme 2
    
    
    # deliveries solely performed by OCs
    deliveries = totalParcels-backup
    
    
    for key in Deliverytime:
        if Deliverytime[key]:
            Avg_DT[key]= sum(Deliverytime[key])/len(Deliverytime[key])
    
    try:
        Avg_DT_OC = sum(Avg_DT.values())/len(Avg_DT) # average deliverytime of a parcel delivered solely by OCs
    except:
        print("Check problem")
    
    # assumption: PF delivery = service level
    try:
        Avg_DT_total = (deliveries*Avg_DT_OC+backup*service_level)/totalParcels
    except:
        print("Check this problem too")
    
    try:
        OC_ratio = used_OCs/totalCouriers # share of matched OCs 
        avg_parcels_per_OC= assigned/used_OCs # average amount of parcels assigned to a matched OC
    except:
        print("Check this third problem")
    
    # average amount of reward paid per used OC
    try:
        avg_compensation1_per_OC = CourierCost1/used_OCs
        avg_compensation2_per_OC = CourierCost2/used_OCs
    except:
        print("Check this fourth problem")

# Load data and run experiments

MasterFile=open(os.path.abspath('')+"/masterfile.csv",'r')
MasterMatrix=np.loadtxt(MasterFile, dtype=str, delimiter=',', skiprows=1)
MasterFile.close()

#for exp in range(1,21):
for exp in range(1,2):
    num_rep = int(MasterMatrix[exp-1][11]) # Number of replications for the simulation
    TotalCost1_rep = []
    TotalCost2_rep = []
    CourierCost1_rep = []
    CourierCost2_rep = []
    totalParcels_rep = []
    totalCouriers_rep = []
    assigned_rep = []
    backup_rep = []
    transshipment_rep = []
    used_OCs_rep = []
    BackupCost_cop_rep = [] # alternative backup cost calculation
    BackupCost_rep = []
    Parcel_arrivals_rep = {} 
    Parcel_arrivals_cop_rep = {} 
    Deliverytime_rep = {} 
    Avg_DT_rep = {}
    parcel_legs_rep = {}
    
    # Fixing the seed
    fix_seed = bool(MasterMatrix[exp-1][1])
    if fix_seed == True:
        seed_value = 42
        random.seed(seed_value)
    
    time_per_day = float(MasterMatrix[exp-1][2]) #1.667 # in hours
    service_level = float(MasterMatrix[exp-1][3]) #0.4 # in hours
    timeframe_P = time_per_day-service_level # time in which new parcels can arrive
    num_parcels = int(MasterMatrix[exp-1][4])#20000 # Total number of parcels per day
    _lambda_avg_P = num_parcels/timeframe_P # number of parcels per hour in whole system
    sup_dem_ratio = float(MasterMatrix[exp-1][5])
    num_couriers = num_parcels * sup_dem_ratio#10000
    _lambda_avg_OC = num_couriers/time_per_day # number of OCs per hour in whole system
    avg_speed = 30.5 # in m/s ca. 37.8 km/h
    #r1 = float(MasterMatrix[exp-1][6]) # fixed reward 
    r2 = 0.0005 # reward paid per m
    variable_compensation = False
    non_stationary = False
    
    if variable_compensation:
        r1 = r2
    else:
        r1 = float(MasterMatrix[exp-1][6]) # fixed reward 
        
    pc = float(MasterMatrix[exp-1][7])
    # TotalCost = 0
    # CourierCost1 = 0
    # CourierCost2 = 0
    # totalParcels = 0
    # totalCouriers = 0
    # assigned = 0
    # backup = 0
    # transshipment = 0
    # used_OCs = 0
    policy = int(MasterMatrix[exp-1][8])   # variable to select the type of policy (0 time look-ahead, 1 distance myopic, 3 optimal)
    capacitated = bool(MasterMatrix[exp-1][9]) # capacicites in boxes considered
    origin_wait_nocap = True # allows that the parcel waits for a crowd courier at its origin, even if there is no capacity at the origin box
    destination_nowait = False
    crowd_courier_cap = 1 # Crowd courier capacity
    box_cap = int(MasterMatrix[exp-1][12])
    # BackupCost_cop = 0 # alternative backup cost calculation
    # BackupCost = 0
    exp_states = {} # expected cost under the optimal policy
    discrete_MDP = True
    compute_shadows = bool(MasterMatrix[exp-1][10]) # true if shadow prices are estimated
    max_num_routings = int(pc/r1) # Determines the maximum number of routings
    time_granularity = float(MasterMatrix[exp-1][14])  # Determines the lenght of time periods
    
    # environment with 25 SPs
    #Source files
    myPath=os.path.abspath('./Input_Data')
    #distanceMatrix_file = open(myPath+"\DMatrixTuple_25.csv",'r')
    distanceMatrix_file = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/DMatrixTuple_25.csv",'r')
    #DistrictsFactors_file=open(myPath+"\Districts_factors_25.csv",'r', newline='')    
    DistrictsFactors_file = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/Districts_factors_25.csv",'r')
    DistrictsFactors_file_demand = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/Districts_factors_25_demand.csv",'r')
    #ProfilesFile = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/crowdProfile.csv",'r')
    #BoxCaps_file = open("./Experiment"+str(exp+1)+"/Input_data/25N_C"+str(box_cap)+".csv",'r', newline='')
    
    
    # DistrictsFactors:
    # each value in each column is appended to a list in order to afterwards create the dictionary
    # with factor_columns[0] being the list of station names and factor_columns[1] being the list of generation factors
    factor_columns = defaultdict(list) 
    with DistrictsFactors_file:
        #reader = csv.reader(DistrictsFactors_file,delimiter=";")
        reader = csv.reader(DistrictsFactors_file,delimiter=",")
        next(reader)
        for row in reader:
            for (i,v) in enumerate(row):
                factor_columns[i].append(v)
    
    # Same for demand factors
    factor_columns_demand = defaultdict(list) 
    with DistrictsFactors_file_demand:
        #reader = csv.reader(DistrictsFactors_file,delimiter=";")
        reader = csv.reader(DistrictsFactors_file_demand,delimiter=",")
        next(reader)
        for row in reader:
            for (i,v) in enumerate(row):
                factor_columns_demand[i].append(v)
    
    df = pd.read_csv("./25_nodes/Experiment"+str(exp)+"/Input_Data/crowdProfile.csv", header=None)
    profiles = {int(row[0]): int(row[1]) for row in df.values}
    
    box_caps = {}
    MaxBox_caps = {}
    #with BoxCaps_file:
    #    #reader = csv.reader(DistrictsFactors_file,delimiter=";")
    #    reader = csv.reader(BoxCaps_file,delimiter=",")
    #    next(reader)
    #    for row in reader:
    #        #for (i,v) in enumerate(row):
    #        box_caps[row[0]] = int(row[1])
    #box_caps = {factor_columns[0][i]: box_cap for i in range(len(factor_columns[0]))} # fill the box capacities dictionary
    
    # I convert the generation/attraction factor strings into floats
    generation_factor = [float(item) for item in factor_columns[1]]
    attraction_factor = [float(item) for item in factor_columns[2]]
    generation_factor_demand = [float(item) for item in factor_columns_demand[1]]
    attraction_factor_demand = [float(item) for item in factor_columns_demand[2]]
    
    
    # dictionary of od-tuples and associated parcel lambdas_ij
    Parcel_lambdas = {} 
    OD_tuple = ()
    for i in range(len(factor_columns_demand[0])):
        for j in range(len(factor_columns_demand[0])):
           if i != j: # if the station names of origin and destination are not the same, a tuple gets created
            OD_tuple = (factor_columns_demand[0][i],factor_columns_demand[0][j])
            P_lambda_ij = generation_factor_demand[i]*attraction_factor_demand[j]*_lambda_avg_P #conditional lambda_ij as value for the tuple_ij
            Parcel_lambdas[OD_tuple] = P_lambda_ij
    
    
    # Define time granularity
    tau_p = 1/time_granularity
    num_per = int(round(time_per_day / tau_p, 0))
    
    # dictionary of od-tuples and associated OC lambdas_ij        
    OC_lambdas = {} 
    OD_tuple = ()
    for i in range(len(factor_columns[0])):
        for j in range(len(factor_columns[0])):
           if i != j:
            OD_tuple = (factor_columns[0][i],factor_columns[0][j])
            OC_lambda_ij = generation_factor[i]*attraction_factor[j]*_lambda_avg_OC
            OC_lambdas[OD_tuple] = OC_lambda_ij
    
    tot_OC_lambda = sum(OC_lambdas[i] for i in OC_lambdas)
    
    
    
    # Distances
    #each value in each column is appended to a list
    distance_columns = defaultdict(list) 
    with distanceMatrix_file:
        reader = csv.reader(distanceMatrix_file,delimiter=";")
        next(reader)
        for row in reader:
            for (i,v) in enumerate(row):
                distance_columns[i].append(v)
    
    distance_list = [float(item) for item in distance_columns[2]]
    
    # dictionary of od-tuples and associated distances
    Distances = {} 
    distance_tuple = ()
    for i in range(len(distance_columns[0])):
        if distance_columns[0][i] != distance_columns[1][i]:
            distance_tuple = (distance_columns[0][i],distance_columns[1][i])
            distance = distance_list[i]
            Distances[distance_tuple] = distance
    
    
        
    
    # run dynamic program
    if discrete_MDP:
        #tau_p = 1/60 # in hours (1/60 means 1 minute periods)
        shadow_prices = {(s,t): 0 for s in factor_columns[0] for t in range(num_per + 1)} # stores the shadow prices per SP and epoch
        shadow_prices_comp = {(s,t): 0 for s in factor_columns[0] for t in range(num_per + 1)} # stores a copy of the shadow prices for computational purposes
        shadow_counts = {(s,t): 0 for s in factor_columns[0] for t in range(num_per + 1)} # stores the counts to compute average shadow prices per SP and epoch
        V_st = {}
        if policy == 3:
            # runs first without shadow prices
            #dynamic_program_disc(f=r1, prob_dist=stats.expon, tau=tau_p)
            #V_st = {}
            space_time_nodes, space_time_arcs, incoming_nodes, outgoing_nodes = build_space_time_network(
                physical_nodes=factor_columns[0], physical_connections=list(Distances.keys()), distances=Distances, avg_speed=avg_speed*3600, time_horizon=time_per_day, tau=tau_p)
            if variable_compensation:
                my_distances = Distances
            else:
                my_distances = {key:1.0 for key in Distances.keys()}
            if non_stationary:
                dynamic_range = num_per - int(service_level/tau_p)
            else:
                dynamic_range = 1
            for t in range(dynamic_range):
            #for t in range(3):
                origin = ('Laim', t)
                for dest in factor_columns[0]:
                    destination = (dest, t + int(service_level/tau_p))
                    subset_arcs, p_arcs, arc_weights, arc_probabilities, updated_incoming_nodes, updated_outgoing_nodes = subset_space_time_network(
                        space_time_nodes, space_time_arcs, incoming_nodes, outgoing_nodes, origin, destination, OC_lambdas, pc, r1, my_distances)
                    V, policies = solve_stochastic_dynamic_program(
                        subset_arcs, p_arcs, arc_weights, arc_probabilities, updated_incoming_nodes, updated_outgoing_nodes, destination, pc, shadow_prices, destination_nowait)
                    for key in V:
                        V_st[key + destination] = V[key]
                
            V_st_before = V_st.copy()
    else:
        tau_p = round(1/_lambda_avg_OC, 4) # this is the minimum decision period
        #shadow_prices = {(s,t): 0 for s in factor_columns[0] for t in np.arange(0,time_per_day - service_level,tau_p)}
        #shadow_counts = {(s,t): 0 for s in factor_columns[0] for t in np.arange(0,time_per_day - service_level,tau_p)}
        shadow_prices = {}
        shadow_prices_comp = {}
        shadow_counts = {}
        for s in factor_columns[0]:
            #for t in np.arange(0,time_per_day - service_level,tau_p):
            for t in np.arange(0,time_per_day,tau_p):
                mytuple = (s,round(t,4))
                shadow_prices[mytuple] = 0
                shadow_prices_comp[mytuple] = 0
                shadow_counts[mytuple] = 0
        if policy == 3:
            # runs first without shadow prices
            dynamic_program(f=r1, prob_dist=stats.expon, tau=tau_p)
    
    init_time = time.time()
    if compute_shadows and capacitated:
        # run simulations to estimate shadow prices
        for rep in range(num_rep):
            TotalCost1 = 0
            TotalCost2 = 0
            CourierCost1 = 0
            CourierCost2 = 0
            totalParcels = 0
            totalCouriers = 0
            assigned = 0
            backup = 0
            transshipment = 0
            used_OCs = 0
            BackupCost_cop = 0 # alternative backup cost calculation
            BackupCost = 0
            
            # Setup and start the simulation
    
            Parcel_arrivals = {} 
            Parcel_arrivals_cop = {} # copy to pop parcels delivered by profesional couriers
            Deliverytime = {} # dictionary of parcel's total delivery time if delivered by OC
            Avg_DT = {} # avergage deliverytimes of parcels delivered by OCs for every key
            #shadow_prices = {} # stores the shadow prices per SP and epoch
            parcel_legs = {}
            parcels = {}
            # fix_seed=1100011
            # random.seed(fix_seed)
            
            # Initialize box capacities
            BoxCaps_file = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/25N_C"+str(box_cap)+".csv",'r')
            with BoxCaps_file:
                #reader = csv.reader(DistrictsFactors_file,delimiter=";")
                reader = csv.reader(BoxCaps_file,delimiter=",")
                next(reader)
                for row in reader:
                    #for (i,v) in enumerate(row):
                    box_caps[row[0]] = int(row[1])
    
            env = simpy.Environment()
    
            for key in Parcel_lambdas:    
                Parcel_arrivals[key]=[] # queue of arriving parcels 
                Parcel_arrivals_cop[key]=[]
                Deliverytime[key]=[] # list of deliverytimes for every key if delivered by OCs
                if Parcel_lambdas[key] > 0:
                    myParcel=ParcelArrival(env,key,None) # generate new parcels for every key
                myCourier=CourierArrival(env,key) # generate new OCs for every key
                if Parcel_lambdas[key] > 0:
                    env.process(arrival_parcels(env,myParcel))
                env.process(arrival_couriers(env,myCourier))
            
            runSimulation()
            # Compute average show prices dividing by the counts
            for key in shadow_prices_comp:
                #shadow_prices = {key: shadow_prices[key]/shadow_counts[key] for key in shadow_prices}
                if shadow_counts[key] > 0:
                    shadow_prices[key] = shadow_prices_comp[key] / shadow_counts[key]
                    
        if discrete_MDP: # Runs the dynamic programs again with updated sahdow prices
            if policy == 3:
                # runs first without shadow prices
                #dynamic_program_disc(f=r1, prob_dist=stats.expon, tau=tau_p)
                V_st = {}
                space_time_nodes, space_time_arcs, incoming_nodes, outgoing_nodes = build_space_time_network(
                    physical_nodes=factor_columns[0], physical_connections=list(Distances.keys()), distances=Distances, avg_speed=avg_speed*3600, time_horizon=time_per_day, tau=tau_p)
                for t in range(dynamic_range):
                #for t in range(1):
                    origin = ('Laim', t)
                    for dest in factor_columns[0]:
                        destination = (dest, t + int(service_level/tau_p))
                        subset_arcs, p_arcs, arc_weights, arc_probabilities, updated_incoming_nodes, updated_outgoing_nodes = subset_space_time_network(
                            space_time_nodes, space_time_arcs, incoming_nodes, outgoing_nodes, origin, destination, OC_lambdas, pc, r1, my_distances)
                        V, policies = solve_stochastic_dynamic_program(
                            subset_arcs, p_arcs, arc_weights, arc_probabilities, updated_incoming_nodes, updated_outgoing_nodes, destination, pc, shadow_prices, destination_nowait)
                        for key in V:
                            V_st[key + destination] = V[key]
                V_st_after = V_st.copy()
    
    # run simulations to assess performance
    compute_shadows = False # go back to False to avoid double computation of shadow prices
    for rep in range(num_rep):
        
        TotalCost1 = 0
        TotalCost2 = 0
        CourierCost1 = 0
        CourierCost2 = 0
        totalParcels = 0
        totalCouriers = 0
        assigned = 0
        backup = 0
        transshipment = 0
        used_OCs = 0
        BackupCost_cop = 0 # alternative backup cost calculation
        BackupCost = 0
        
        # Setup and start the simulation
    
        Parcel_arrivals = {} 
        Parcel_arrivals_cop = {} # copy to pop parcels delivered by profesional couriers
        Deliverytime = {} # dictionary of parcel's total delivery time if delivered by OC
        Avg_DT = {} # avergage deliverytimes of parcels delivered by OCs for every key
        #shadow_prices = {} # stores the shadow prices per SP and epoch
        parcel_legs = {}
        parcels = {}
        
        # Initialize box capacities
        BoxCaps_file = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/25N_C"+str(box_cap)+".csv",'r')
        with BoxCaps_file:
            #reader = csv.reader(DistrictsFactors_file,delimiter=";")
            reader = csv.reader(BoxCaps_file,delimiter=",")
            next(reader)
            for row in reader:
                #for (i,v) in enumerate(row):
                box_caps[row[0]] = int(row[1])
        
        # fix_seed=1100011
        # random.seed(fix_seed)
    
        env = simpy.Environment()
    
        for key in Parcel_lambdas:    
            Parcel_arrivals[key]=[] # queue of arriving parcels 
            Parcel_arrivals_cop[key]=[]
            Deliverytime[key]=[] # list of deliverytimes for every key if delivered by OCs
            if Parcel_lambdas[key] > 0:
                myParcel=ParcelArrival(env,key,None) # generate new parcels for every key
            myCourier=CourierArrival(env,key) # generate new OCs for every key
            if Parcel_lambdas[key] > 0:
                env.process(arrival_parcels(env,myParcel))
            env.process(arrival_couriers(env,myCourier))
        
        runSimulation()
        
        TotalCost1_rep.append(TotalCost1)
        TotalCost2_rep.append(TotalCost2)
        CourierCost1_rep.append(CourierCost1)
        CourierCost2_rep.append(CourierCost2)
        totalParcels_rep.append(totalParcels)
        totalCouriers_rep.append(totalCouriers)
        assigned_rep.append(assigned)
        backup_rep.append(backup)
        transshipment_rep.append(transshipment)
        used_OCs_rep.append(used_OCs)
        BackupCost_cop_rep.append(BackupCost_cop)
        BackupCost_rep.append(BackupCost)
        
        Parcel_arrivals_rep[rep] = Parcel_arrivals.copy()
        Parcel_arrivals_cop_rep[rep] = Parcel_arrivals_cop.copy()
        Deliverytime_rep[rep] = Deliverytime.copy()
        Avg_DT_rep[rep] = Avg_DT.copy()
        parcel_legs_rep[rep] = parcel_legs.copy()
    
    # Compute run time
    run_time = time.time() - init_time
    
    # Write results
    outputs_file = open("./25_nodes/Experiment"+str(exp)+"/Output_data/outputs.csv",'w', newline='')
    outputs_DTs = open("./25_nodes/Experiment"+str(exp)+"/Output_data/outputs_DTs.csv",'w', newline='')
    outputs_Legs = open("./25_nodes/Experiment"+str(exp)+"/Output_data/outputs_Legs.csv",'w', newline='')
    outputs_Values = open("./25_nodes/Experiment"+str(exp)+"/Output_data/outputs_Vals.csv",'w', newline='')
    outputs_LMult = open("./25_nodes/Experiment"+str(exp)+"/Output_data/outputs_LMults.csv",'w', newline='')
    
    V_st0={key:V_st[key] for key in V_st.keys() if key[1]==0}
    
    BoxCaps_file = open("./25_nodes/Experiment"+str(exp)+"/Input_Data/25N_C"+str(box_cap)+".csv",'r')
    with BoxCaps_file:
        #reader = csv.reader(DistrictsFactors_file,delimiter=";")
        reader = csv.reader(BoxCaps_file,delimiter=",")
        next(reader)
        for row in reader:
            #for (i,v) in enumerate(row):
            MaxBox_caps[row[0]] = int(row[1])
    # Theoretical bounds
    if policy==3:
        exp_capPenalty = sum(shadow_prices[key]*MaxBox_caps[key[0]] for key in shadow_prices.keys())
        upper_bound = sum(V_st0[key]*Parcel_lambdas[(key[0],key[2])]*14 for key in V_st0.keys() if key[0]!=key[2]) * (num_parcels/sum(Parcel_lambdas[(key[0],key[2])]*14 for key in V_st0.keys() if key[0]!=key[2]))
        lower_bound = upper_bound - exp_capPenalty
    else:
        exp_capPenalty = "NA"
        upper_bound = "NA"
        lower_bound = "NA"
    
    with outputs_file:
        writer = csv.writer(outputs_file, delimiter=',')
        #reader(DistrictsFactors_file,delimiter=",")
        row1 = TotalCost1_rep
        row1.insert(0, "Total Cost by Replication")
        writer.writerow(row1)
        
        row2 = CourierCost1_rep
        row2.insert(0, "Courier Cost by Replication")
        writer.writerow(row2)
        
        row3 = totalParcels_rep
        row3.insert(0, "Total Parcels by Replication")
        writer.writerow(row3)
        
        row4 = totalCouriers_rep
        row4.insert(0, "Total OCs by Replication")
        writer.writerow(row4)
        
        row5 = assigned_rep
        row5.insert(0, "Assigned OCs by Replication")
        writer.writerow(row5)
        
        row6 = backup_rep
        row6.insert(0, "Professional Shipments by Replication")
        writer.writerow(row6)
        
        row7 = BackupCost_rep
        row7.insert(0, "Professional Cost by Replication")
        writer.writerow(row7)
        
        writer.writerow(["Run time", run_time])
        
        #Write upper and lower bounds
        writer.writerow(["Upper Bound", upper_bound])
        writer.writerow(["Lower Bound", lower_bound])
    
    with outputs_DTs:
        writer = csv.writer(outputs_DTs, delimiter=',')
        #reader(DistrictsFactors_file,delimiter=",")
        #row1 = Avg_DT_rep
        #row1.insert(0, "Total Cost by Replication")
        for key in Avg_DT_rep:
            #row = "Replication" + str(key)
            #writer.writerow(row)
            #row = []
            for pair in Avg_DT_rep[key]:
                row = []
                row.append(pair)
                row.append(Avg_DT_rep[key][pair])
                writer.writerow(row)
    
    with outputs_Legs:
        writer = csv.writer(outputs_Legs, delimiter=',')
        #reader(DistrictsFactors_file,delimiter=",")
        #row1 = Avg_DT_rep
        #row1.insert(0, "Total Cost by Replication")
        for key in parcel_legs_rep:
            #row = "Replication" + str(key)
            #writer.writerow(row)
            #row = []
            for pair in parcel_legs_rep[key]:
                row = []
                row.append(pair)
                row.extend(parcel_legs_rep[key][pair])
                writer.writerow(row)
                
    with outputs_Values:
        writer = csv.writer(outputs_Values, delimiter=',')
        #reader(DistrictsFactors_file,delimiter=",")
        #row1 = Avg_DT_rep
        #row1.insert(0, "Total Cost by Replication")
        for key in V_st:
            #row = "Replication" + str(key)
            #writer.writerow(row)
            #row = []
            #for pair in parcel_legs_rep[key]:
            row = []
            #row.append(key, V_st[key])
            row.extend([key, V_st[key]])
            writer.writerow(row)
            
    with outputs_LMult:
        writer = csv.writer(outputs_LMult, delimiter=',')
        #reader(DistrictsFactors_file,delimiter=",")
        #row1 = Avg_DT_rep
        #row1.insert(0, "Total Cost by Replication")
        for key in shadow_prices:
            #row = "Replication" + str(key)
            #writer.writerow(row)
            #row = []
            #for pair in parcel_legs_rep[key]:
            row = []
            #row.append(key, V_st[key])
            row.extend([key, shadow_prices[key]])
            writer.writerow(row)


