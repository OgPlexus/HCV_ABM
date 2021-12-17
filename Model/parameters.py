import numpy as np

Params = {#'identity': 'locSIS_rndTreatA0',
          'identity': 'reproducibility_test',
          'p0': 1,
          'pclr0': 10*[0.26],
          #'pclr0': np.linspace(0, 0.5, 10),
          'pclr1': np.zeros(10),
          'adhere': np.linspace(0, 1, 10)[9],
          'treatment start': 1825, # 4000
          'sis start': 4000, # 1825
          'rand treat': True,
          'local sis': True,
          'use': 0.5,
          'decay': 0.96,
          'share': 0.10,
          'relocate': 1/90,
          'nep': [0] * 10,
          'treatment': np.linspace(0.05, 0.5, 10),
          'connect': 1 / 9,
          'nep daily rate': 1 / 7,
          'death': 1/365}
          #'clr nep': np.linspace(0, 1, 10),
          #'clr treat': np.linspace(0, 1, 10)}


"""
Description of 'identity' labels:
=================================

'popularX0':       beta = [[c, x], [y, c]], w/ x+y=c and p0 = 0.X0
'connectX0':       beta = [[x, y], [y, x]], w/ x**2 + y**2 = c**2 and p0 = 0.X0 
'deg_sweep00':     beta = [[c, 0], [0, 0]], w/ c in [0.3 * L, L], L = 1/sqrt(1000)
'const_beta00':    beta = [[c, 0], [0, 0]], w/ c = L/2, pclr0 ranges in [0.2, 0.3]
'const_beta50':    "" "", pclr0 ranges in [0.00001, 0.5]
'clr_degXY':       beta = [[c, 0], [0, 0]], w/ c in [0.3 * L, L], X in [A, Z] & Y in [0, 9] |-> [0.00001, 0.5] (clr frac.) 
'treat_adhereXY':  beta = const. w/ c = L/2, run_index specifies treatment frac., Y in [0, 9] |-> not-sharing frac. (0, 1)
'inf_time_locXY':  beta = const. w/ c = L/2, X = A, B, ..., Y in [0, 9] |-> repetitions
'rnd_treatXY':     beta = const. XY for repetitions, random selection for treatment, run_index runs over treat. frac. (0.05 to 0.5), sis start = 4000
'deg_treatXY':     beta = const. XY for repetitions, treatment selection by degree, run_index runs over treat. frac. (0.05 to 0.5), sis start = 4000
'local_SISXY':     beta = const. XY for reps., run_index runs over 10 identical runs (pclr0 = 0.26), local means 1 SIS
'nlocal_SISXY':    beta = const. XY for reps., run_index runs over 10 identical runs (pclr0 = 0.26), nlocal means 25 SIS locations
'sis_randtreatXY:  beta = const. XY for reps., run_index runs over (random) treatment frac: (0.05 to 0.5)

"""



"""
Description of 'trial' labels:
=============================

'T(N)CLRS':     testing treatment among (non-)clearers only, with p_share = 0.45
'T14(N)CLRS':   testing treatment among (non-)clearers only, with p_share = 0.14
'32T(N)CLR:     test of the 0 to 32 clearance percent range (p_share = 0.2)
'32NOINT':      no interventions for the (0, 32) range
'33NOINT':      no interventions for the (1, 33) range
'NEP(N)CLR':    all in NEP are (non-)clearers, sweep over p_clear
'SWP':          sweeping over p_clear w/o interventions
'rSWP':         a re-do of SWP where the R0 and graph components can be accounted for
'NEP':          fixed p_clear, sweep over NEP clearer fraction
'NEP10'         fixed p_clear (26), sweep over NEP clearer fraction (w/ 10% of pop. in NEPs)
'NEPswp'        fixed p_clear (= 26) sweep across NEP % in total pop. (from 0 to 100%)
'TRTswp'        same as NEPswp for treatment
"""



"""
Variable name:      Description
==============================================================================================================

'nep':            # fraction of people in Needle Exchange Program  
'trial':          # label to distiguish parameter classes: 'Y' => 'exchange' = 0.02, 'Z' => 'exchange' = 0.01          
'use':            # probability a node will inject drugs per day   
'lend':           # depreciated                                    
'decay':          # fraction of infected needles that clear HCV daily  
'radius':         # distance people will connect within                
'share':          # probability a node will share drugs that day      
'max nep':        # maximum number of individuals in NEP (can set to 0 to turn off NEPs)   
'break':          # probability of an edge break (per edge, per time step)         
'max treat':      # maximum number of individuals in treatment program (can set to 0 to turn off treatment)
'exchange':       # depreciated (now we use p_share**2 or 0, if one user in edge is in treatment)            
'connect':        # probability that, given nodes are with 'radius', they will form an edge   
'nep type':       # indicates which group the program focuses on 
'treat type':     # indicates which group the program focuses on
'nep daily rate': # probability per day that an individual in an NEP will swap out needles
'recover':        # probability of recovery if user can clear the infection                   
'death':          # probability of leaving the network                                          
'clearances':     # range of clearance percentages to consider

=============================================================================================================
"""

""" 
Description of infection-avenue abbreviations:
==============================================

'nss':    needle becomes infected due to a ss-interaction (sick-sick)
'nhs':    needle becomes infected due to a hs-interaction (healthy person is the one to share the needle)
'nsh':    needle becomes infected due to a sh-interaction (sick person shares the needle)
'ns':     needle becomes infected due to a self-injection event (where person is sick)
'hh':     two healthy users become infected when one shares an infected needle with the other
'hs':     healthy user becomes infected when healthy user shares an infected needle (using it first)
'sh':     healthy user becomes infected when an infected user shares a needle
'h':      healthy user becomes infected in self-injection event (using an infected needle already in their collection)

"""
