import networkx as nx
import random as rnd
import numpy as np
from random import choice
from parameters import Params
from functions import dist, clearance, Beta, GBM2, HCV_ABM
import pickle
import sys

# these vars interface with the cluster
JOB_ID = int(sys.argv[1])
TASK_ID = int(sys.argv[2])

# save_obj and load_obj save and load resp. '.pkl' files:
def save_obj(obj, name):
    with open('Python/obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('Python/Python_Home/obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

# parameters:
run_index = 0
N = 1000
t_steps = 3650

p0 = Params['p0']
pclr0 = Params['pclr0'][run_index]
pclr1 = Params['pclr1'][run_index]
identity = Params['identity']
beta = Beta(p0, identity[:-2])[run_index]

g = GBM2(p0, pclr0, pclr1, beta)

# randomly selects patient 0 and assigns attributes:
patient_0 = rnd.choice(range(N))
g.nodes[patient_0]['health'] = 1      # <-- patient 0 is infected

""" running the sim. """
p_decay = Params['decay']
p_death = Params['death']

(AB, Needles, deg_dist, relocations) = HCV_ABM(g, beta, p_decay, p_death, t_steps, run_index)
sim_data = [AB, Needles, deg_dist, relocations]

# save_obj(multi_pops, 'GBM pops beta%s%s_%s%s %s' %(identity, run_index, int(beta[0][0] * 100), int(beta[0][1] * 100), TASK_ID - 1))
# save_obj(inf_gen, 'GBM inf gen beta%s%s_%s%s %s' %(identity, run_index, int(beta[0][0] * 100), int(beta[0][1] * 100), TASK_ID - 1))

save_obj(sim_data, 'sim data %s beta%s run%s' %(identity, run_index, TASK_ID-1))

