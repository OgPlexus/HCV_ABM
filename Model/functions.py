import networkx as nx
import random as rnd
import numpy as np
from parameters import Params

""" parameters from Params """
p_use = Params['use']  # 0.5
p_share = Params['share']  # 0.45
p_decay = Params['decay']  # 0.96
p_death = Params['death']  # 3 * 10 ** (-4)
p_connect = Params['connect']  # 1 / 9
p_adhere = Params['adhere']
t_start = Params['treatment start']
sis_start = Params['sis start']
rand_treat = Params['rand treat']
local_SIS = Params['local sis']

# --------------------------------------------------------------------------------------------------------- #

""" function 1. recovery """
def recovery(Node):
    if Node['health'] == 1 and rnd.random() < Node['immunity']:
        Node['health'] = 0  # the infected user recovers with a probability per time step
        return True
    else:
        return False
# --------------------------------------------------------------------------------------------------------- #

""" function 2. Euclidean 2-distance """
def dist(a, b):
    ax, ay = a
    bx, by = b
    return np.sqrt((bx - ax) ** 2 + (by - ay) ** 2)

# --------------------------------------------------------------------------------------------------------- #

""" function 3. returns the beta matrices for the right occasion ('connect', 'popular', etc.) """
def Beta(p0, ID_type):
    beta = []
    if ID_type == 'connect':
        alpha = 2 * p0 * (1 - p0)
        a = np.sqrt(1 - alpha)
        b = np.sqrt(alpha)
        for i in range(10):
            theta = np.pi * i/18
            r = 1/np.sqrt(1000)  # <-- the average distance between nodes
            u, v = r * np.cos(theta), r * np.sin(theta)
            beta.append([[u/a, v/b], [v/b, u/a]])
    elif ID_type == 'popular':
        for i in range(10):
            c = 1/np.sqrt(1000)
            x = c * i/9
            y = c - x
            beta.append([[c, x], [y, c]])
    elif ID_type == 'deg_sweep':
        avg_len = 1/np.sqrt(1000)
        C = np.linspace(0.3 * avg_len, avg_len, 10)
        for c in C:
            beta.append([[c, 0], [0, 0]])
    elif ID_type == 'clr_deg':
        avg_len = 1 / np.sqrt(1000)
        C = np.linspace(0.3 * avg_len, avg_len, 10)
        for c in C:
            beta.append([[c, 0], [0, 0]])
    else:  # <-- we leave the constant beta case as the default!
        c = 1/(2 * np.sqrt(1000))
        for i in range(10):
            beta.append([[c, 0], [0, 0]])
    return beta

# --------------------------------------------------------------------------------------------------------- #

""" function 4. counts the numbers of infected people """
def pop_count(nodes):
    a_clr, a_nclr, b_clr, b_nclr = 0, 0, 0, 0
    for x in nodes:
        if nodes[x]['health'] == 1:
            if nodes[x]['group'] == 0 and nodes[x]['immunity'] > 0:  # <-- in group 0 and can clear
                a_clr += 1
            elif nodes[x]['group'] == 0 and nodes[x]['immunity'] == 0:  # <-- in group 0 and cannot clear
                a_nclr += 1
            elif nodes[x]['group'] == 1 and nodes[x]['immunity'] > 0:  # <-- in group 1 and can clear
                b_clr += 1
            elif nodes[x]['group'] == 1 and nodes[x]['immunity'] == 0:  # <-- in group 1 and cannot clear
                b_nclr += 1
    return (a_clr, a_nclr, b_clr, b_nclr)

# --------------------------------------------------------------------------------------------------------- #

""" function 5. returns the fraction of area in the overlap of two circles """
def overlap_fraction(d, r1, r2):
    # if circles don't overlap
    if r1 + r2 <= d:
        return 0

    # if edges of circles touch
    elif abs(r1 - r2) <= d:
        s1 = r1**2 * np.arccos((r1**2 - r2**2 + d**2)/(2 * r1 * d))  # area of sector 1
        s2 = r2**2 * np.arccos((r2**2 - r1**2 + d**2)/(2 * r2 * d))  # area of sector 2
        q = np.sqrt(4 * r1**2 * r2**2 - (r1**2 + r2**2 - d**2)**2)/2
        area = s1 + s2 - q
        return area/(np.pi * r1**2 + np.pi * r2**2 - area)

    # if one circle contains the other
    else:
        return (min(r1, r2)/max(r1, r2))**2

# --------------------------------------------------------------------------------------------------------- #

""" function 6. returns a list of tuples indicating who within each group can clear """
def clearance(alpha0, num, p0, p1, recover_rate=0.03333):

    if num > 1:
        n_0 = int(alpha0 * num)  # number of type 0
        n_1 = num - int(alpha0 * num)  # number of type 1

        n_00 = n_0 - int(p0 * n_0)  # num of type 0 w/o clearance allele
        n_01 = int(p0 * n_0)  # num of type 0 w/ clearance allele
        n_10 = n_1 - int(p1 * n_1)  # num of type 1 w/o clearance allele
        n_11 = int(p1 * n_1)  # num of type 1 w/ clearance allele

        """ the first 0 or 1 indicates the group type """
        """ the second 0 or 1 indicates whether the clearance allele is present """
        """ if value in array is nonzero, clearance allele is present """

        a_00 = list(zip(np.zeros(n_00), [0] * n_00))
        a_01 = list(zip([recover_rate] * n_01, [0] * n_01))
        a_10 = list(zip(np.zeros(n_10), [1] * n_10))
        a_11 = list(zip([recover_rate] * n_11, [1] * n_11))

        c = a_00 + a_01 + a_10 + a_11
        rnd.shuffle(c)
        return c

    else:  # the case when num = 1
        if rnd.random() < alpha0:  # randomly decides if type 0
            if rnd.random() < p0:  # if agent has clearance capability
                return recover_rate, 0  # type 0 w/ clearance capability
            else:
                return 0, 0  # type 0 w/o clearance capability
        elif rnd.random() < p1:  # if the agent is type 1 then if has clearance capability
            return recover_rate, 1  # type 1 w/ clearance capability
        else:
            return 0, 1  # type 1 w/o clearance capability

# --------------------------------------------------------------------------------------------------------- #

""" function 7. gives a list of the n highest degree nodes in g (in decreasing order of degree) '"""
def high_degrees(g, n):
    nodes = list(g.nodes)
    degs = [len(list(g.neighbors(x))) for x in g]
    max_degs = []
    for i in range(n):
        max_deg = degs.index(max(degs))
        degs[max_deg] = -1
        max_degs.append(nodes[max_deg])
    return max_degs

# --------------------------------------------------------------------------------------------------------- #

""" function 8. puts a given frac. of the pop. in a treatment group """
def treatment(g, pTreat, random=True):
    treat_num = int(round(pTreat * 1000))
    nodes = list(g.nodes)
    to_treat = []
    if random:  # <-- randomly chosen
        for i in range(treat_num):
            new_treat = rnd.choice([x for x in nodes if not x in to_treat])
            to_treat.append(new_treat)
            g.nodes[new_treat]['treat'] = 1
    else:  # <-- chosen based off of degree
        high_degs = high_degrees(g, treat_num)
        for i in high_degs:
            g.nodes[i]['treat'] = 1

# --------------------------------------------------------------------------------------------------------- #

""" function 9a. initiates safe injection sites (SISs) """
def safe_injection(g, localized=True):

    if localized:  # <-- localized SIS version

        centre = (0.5, 0.5)
        search_radius = np.sqrt(300 / (np.pi * 1000))

        candidates = []
        for each in g:
            pos = g.nodes[each]['pos']
            if dist(centre, pos) < search_radius:
                candidates.append(each)

        rnd.shuffle(candidates)
        for i in range(250):
            new_participant = candidates[i]
            g.nodes[new_participant]['sis'] = ('0.5', '0.5')
            g.nodes[new_participant]['p_share'] = 0

    else:  # <-- distributed SIS version

        mixed_nodes = list(g.nodes)
        rnd.shuffle(mixed_nodes)

        loc_map = {0: 0, 1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 4}

        candidates = {(i, j): [] for i in range(5) for j in range(5)}
        for each in mixed_nodes:
            pos = g.nodes[each]['pos']
            x, y = int(round(6 * pos[0])), int(round(6 * pos[1]))
            i, j = loc_map[x], loc_map[y]

            # assign SIS locations:
            if len(candidates[(i, j)]) < 10:
                candidates[(i, j)].append(each)
                g.nodes[each]['sis'] = (i, j)
                g.nodes[each]['p_share'] = 0

    return candidates

# --------------------------------------------------------------------------------------------------------- #

""" function 9b. finds a new SIS participant """
def find_SIS_participant(g, location):
    mixed_nodes = list(g.nodes)
    rnd.shuffle(mixed_nodes)
    if location == ('0.5', '0.5'):
        deficit = 250 - len([x for x in g if g.nodes[x]['sis'] == location])
        centre = (0.5, 0.5)
        search_radius = np.sqrt(300 / (np.pi * 1000))
        for search in range(deficit):  # <-- attempts to make up the deficit in that SIS location
            for new in mixed_nodes:
                pos = g.nodes[new]['pos']
                if dist(centre, pos) < search_radius and type(g.nodes[new]['sis']) == str:
                    g.nodes[new]['sis'] = location
                    g.nodes[new]['p_share'] = 0
                    break
    else:
        loc_map = {0: 0, 1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 4}
        deficit = 10 - len([x for x in g if g.nodes[x]['sis'] == location])
        for search in range(deficit):  # <-- attempts to make up the deficit in that SIS location
            for new in mixed_nodes:
                pos = g.nodes[new]['pos']
                x, y = int(round(6 * pos[0])), int(round(6 * pos[1]))
                i, j = loc_map[x], loc_map[y]
                if (i, j) == location and type(g.nodes[new]['sis']) == str:
                    g.nodes[new]['sis'] = location
                    g.nodes[new]['p_share'] = 0
                    break

# --------------------------------------------------------------------------------------------------------- #

""" function 10. finds a new patient for treatment """
def find_new_patient(g, node, ptreat, random=True, find_someone_else=False):
    g.nodes[node]['treat'] = 0  # <-- no longer in treatment
    g.nodes[node]['p_share'] = p_share  # <-- resumes sharing w/ others
    g.nodes[node]['treat_time'] = 0  # <-- treat time starts over
    if random:
        if find_someone_else:
            new_patient = rnd.choice([x for x in g if g.nodes[x]['treat'] == 0 and not x == node])
            g.nodes[new_patient]['treat'] = 1
        else:
            new_patient = rnd.choice([x for x in g if g.nodes[x]['treat'] == 0])
            g.nodes[new_patient]['treat'] = 1
    else:
        if find_someone_else:
            high_degs = high_degrees(g, 1+int(round(ptreat * 1000)))  # <-- we add 1 in case 'node' is in high_degs
            for new_patient in high_degs:
                if g.nodes[new_patient]['treat'] == 0 and not new_patient == node:
                    g.nodes[new_patient]['treat'] = 1
                    break
        else:
            high_degs = high_degrees(g, int(round(ptreat * 1000)))
            for new_patient in high_degs:
                if g.nodes[new_patient]['treat'] == 0:
                    g.nodes[new_patient]['treat'] = 1
                    break
    return new_patient

# --------------------------------------------------------------------------------------------------------- #

""" function 11a. builds a random geometric block model graph, v.1 """
def GBM1(p0, pclr0, pclr1, b):

    """ this version forms edges with probability given by <overlap_fraction> (no edge-weights)"""

    """
    p0:    fraction of graph in group 0
    pclr0: fraction of clearers in group 0
    pclr1: fraction of clearers in group 1
    bXY:   connection distances (in matrix form) for folks in group X to group Y

    """

    N = 1000
    p_use = Params['use']
    p_share = Params['share']
    immunity = clearance(p0, N, pclr0, pclr1)
    radius = max(b[0] + b[1])  # <-- will be used to record neighbors

    base = nx.empty_graph()
    nhbrs = {x: [] for x in range(N)}  # <-- will store nearby nodes for each node
    loc = {}

    for node in range(N):

        new_pos = (rnd.random(), rnd.random())
        loc[node] = new_pos

        base.add_node(
               node,
               group=immunity[node][1],
               pos=new_pos,
               needles=[rnd.randint(1, 3), 0],
               health=0,
               immunity=immunity[node][0],
               p_use=p_use,
               p_share=p_share,
               treat_time=0,
               treat=0,
               nep=0)

        for nhbr in loc:
            # we use 2*radius bc that is the maximum distance that folks can communicate within
            if dist(loc[nhbr], new_pos) < 2*radius and not nhbr == node:
                nhbrs[node].append(nhbr)

    seen_pairs = []  # <-- to store exhausted node-pairs
    for node in nhbrs:
        pos = loc[node]
        proximity = nhbrs[node]
        node_grp = immunity[node][1]
        for each in proximity:
            nhbr_grp = immunity[each][1]
            d = dist(pos, loc[each])
            r1 = b[node_grp][nhbr_grp]
            r2 = b[nhbr_grp][node_grp]
            pair = (min(node, each), max(node, each))
            if not pair in seen_pairs:
                seen_pairs.append(pair)
                if rnd.random() < overlap_fraction(d, r1, r2):
                    base.add_edge(node, each)

    return base

# --------------------------------------------------------------------------------------------------------- #

""" function 11b. builds a random geometric block model graph, v.2 """
def GBM2(p0, pclr0, pclr1, b_matrix):

    """ this version forms all edges with <overlap_fraction> > 0 and weighs edges by that fraction """

    """
    p0:    fraction of graph in group 0
    pclr0: fraction of clearers in group 0
    pclr1: fraction of clearers in group 1
    bXY:   connection distances (in matrix form) for folks in group X to group Y

    """

    N = 1000
    p_use = Params['use']
    p_share = Params['share']
    immunity = clearance(p0, N, pclr0, pclr1)
    max_radius = max(b_matrix[0] + b_matrix[1])  # <-- max. radius that nodes could connect within

    base = nx.empty_graph()
    loc = {}

    for node in range(N):

        new_pos = (rnd.random(), rnd.random())

        base.add_node(
            node,
            group=immunity[node][1],
            pos=new_pos,
            needles=[rnd.randint(1, 3), 0],
            health=0,
            immunity=immunity[node][0],
            p_use=p_use,
            p_share=p_share,
            treat_time=0,
            treat=0,
            sis='not in SIS',
            nep=0)

        for nhbr in loc:
            d = dist(loc[nhbr], new_pos)  # <-- distance between new node and old nodes
            if d < 2*max_radius:  # <-- checks if nodes are close enough to connect
                grp1 = base.nodes[node]['group']  # <-- 'node' group
                grp2 = base.nodes[nhbr]['group']  # <-- 'nhbr' group
                r1 = b_matrix[grp1][grp2]  # <-- distance grp1 folks will go to link up with grp2 folks
                r2 = b_matrix[grp2][grp1]  # <-- distance grp2 folks will go to link up with grp1 folks
                overlap = overlap_fraction(d, r1, r2)
                if overlap > 0:  # <-- if radii overlap
                    base.add_edge(node, nhbr, weight=overlap)

        # lastly, we add the new location to the locations-dictionary, 'loc'
        loc[node] = new_pos

    return base


# --------------------------------------------------------------------------------------------------------- #

""" function 12. runs through the interaction cases """
def interaction_cases(node, g):
    nodes = g.nodes
    Node = nodes[node]

    nhbrs = list(g.neighbors(node))
    rnd.shuffle(nhbrs)

    """ to record who becomes infected """
    infections = []

    for nhbr in nhbrs:  # iterates through neighbors of the node

        p_exchange = g.get_edge_data(nhbr, node)['weight']  # <-- probability of an exchange of needles

        if rnd.random() < p_exchange:
            # ------------------------------------------------------------------------------------ #
            # ------ We run through the interaction cases ---------------------------------------- #
            # ------------------------------------------------------------------------------------ #

            # like meets like case:
            if Node['health'] == g.nodes[nhbr]['health']:

                # this is the sick meets sick case
                if Node['health'] == 1:

                    # host, x, with the higher p_lend will lend
                    x = nhbr if rnd.random() < 0.5 else node  # coin flip to see who shares

                    ax = nodes[x]['needles'][0]  # num of clean needles for user x
                    bx = nodes[x]['needles'][1]  # num of infected needles for user x

                    # modifying the num. of clean needles in an infected user's collection
                    if rnd.random() < ax / (ax + bx):  # chances that the needle is clean

                        nodes[x]['needles'][0] -= 1  # num of clean needles for user x
                        nodes[x]['needles'][1] += 1  # num of infected needles for user x

                # this is the healthy meets healthy case
                else:

                    # coin flip to see who shares
                    x = nhbr if rnd.random() < 0.5 else node

                    ax = nodes[x]['needles'][0]  # num of clean needles for user x
                    bx = nodes[x]['needles'][1]  # num of infected needles for user x

                    # user injects first then shares (there's a chance the user has an infected needle)
                    if rnd.random() < bx / (ax + bx):  # chances that an infected needle is selected

                        nodes[nhbr]['health'] = 1  # both users become
                        Node['health'] = 1         # infected

                        # record the infections (user x first)
                        infections.extend([x, [y for y in [node, nhbr] if not y == x][0]])

            else:  # healthy meets sick case

                # picks out the healthy one in the pair
                healthy_user = [user for user in [node, nhbr] if nodes[user]['health'] == 0][0]

                # picks out the sick one in the pair
                sick_user = [user for user in [node, nhbr] if nodes[user]['health'] == 1][0]

                # host with the higher p_lend gives the needle
                # in this case, if the healthy_user selects the needle
                if rnd.random() < 0.5:  # coin flip to see if healthy_user shares

                    a = nodes[healthy_user]['needles'][0]  # num of clean needles for user x
                    b = nodes[healthy_user]['needles'][1]  # num of infected needles for user x

                    # if a clean needle was selected by the uninfected user
                    if rnd.random() < a / (a + b):

                        nodes[healthy_user]['needles'][0] -= 1  # num of clean needles of healthy_user
                        nodes[healthy_user]['needles'][1] += 1  # num of infctd needles of healthy_user

                    else:
                        # if selection of an infected needle occurs, the healthy user is
                        # assumed to use it before sharing, and therefore become infected
                        nodes[healthy_user]['health'] = 1

                        # record the infection
                        infections.append(healthy_user)

                # infected user gives needle to uninfected user (again, the lender is assumed
                # to inject first then lend, in which case the healthy user becomes infected)
                else:

                    a = nodes[sick_user]['needles'][0]  # num of clean needles for user x
                    b = nodes[sick_user]['needles'][1]  # num of infected needles for user x

                    # healthy user becomes infected regardless of the type of needle selected
                    nodes[healthy_user]['health'] = 1

                    # record the infection
                    infections.append(healthy_user)

                    # if a clean needle is selected by the infected user:
                    if rnd.random() < a / (a + b):
                        nodes[sick_user]['needles'][0] -= 1  # num of clean needles for user x
                        nodes[sick_user]['needles'][1] += 1  # num of infected needles for user x

            # ------------------------------------------------------------------------------------ #
            # ------------------------------------------------------------------------------------ #
            # ------------------------------------------------------------------------------------ #

    return infections

# --------------------------------------------------------------------------------------------------------- #

""" function 13. relocates the nodes """
def relocate(node, g, b_matrix):

    nodes = g.nodes

    # identify prior edge connections
    nhbrs = list(g.neighbors(node))
    old_edges = [(node, x) for x in nhbrs]

    # remove prior edge connections
    for edge in old_edges:
        g.remove_edge(*edge)

    # generate a new node position
    new_pos = (rnd.random(), rnd.random())

    # change to new position
    g.nodes[node]['pos'] = new_pos

    # max. distance nodes can connect within
    max_radius = max(b_matrix[0] + b_matrix[1])

    # dictionary of all node positions
    loc = {x: nodes[x]['pos'] for x in g}

    # forms new edge connections if radii overlap
    for nhbr in loc:
        d = dist(new_pos, loc[nhbr])  # <-- distance between new node and old nodes
        if 0 < d < 2 * max_radius:  # <-- checks if nodes are close enough to connect (and distinct)
            grp1 = g.nodes[node]['group']  # <-- 'node' group
            grp2 = g.nodes[nhbr]['group']  # <-- 'nhbr' group
            r1 = b_matrix[grp1][grp2]  # <-- distance grp1 folks will go to link up with grp2 folks
            r2 = b_matrix[grp2][grp1]  # <-- distance grp2 folks will go to link up with grp1 folks
            overlap = overlap_fraction(d, r1, r2)
            if overlap > 0:  # <-- if radii overlap
                g.add_edge(node, nhbr, weight=overlap)

    return new_pos

# --------------------------------------------------------------------------------------------------------- #

""" function 14. the simulation """
def HCV_ABM(g, b_matrix, p_decay, p_death, t_steps, run_index):

    p_treat = Params['treatment'][run_index]

    total_in_treat, total_in_sis, total_in_both = [], [], []

#    num_in_SIS = {('0.5', '0.5'): np.zeros(t_steps)}
#    num_in_SIS = {(i, j): np.zeros(t_steps) for i in range(5) for j in range(5)}

    components = []
    sis_dicts = []

    nodes = g.nodes

    """ to record the time and location of an infection/recovery (starting w/ patient 0) """
    for node in g:
        if g.nodes[node]['health'] == 1:
            pos_0 = g.nodes[node]['pos']
            break
    time_loc = [(0, 'i', pos_0)]

    """ to store the degree distribution over time """
    deg_dist = []

    """ storing initial infected populations """
    a_clr, a_nclr, b_clr, b_nclr = pop_count(nodes)

    A_clr, A_nclr = np.zeros(t_steps+1), np.zeros(t_steps+1)
    B_clr, B_nclr = np.zeros(t_steps+1), np.zeros(t_steps+1)

    A_clr[0], A_nclr[0] = a_clr, a_nclr
    B_clr[0], B_nclr[0] = b_clr, b_nclr

    clean_ndls = np.array([sum([nodes[nds]['needles'][0] for nds in nodes])])  # initial num of clean needles
    infct_ndls = np.array([sum([nodes[nds]['needles'][1] for nds in nodes])])  # initial num of infected needles

    """ a time step corresponds to approximately one day """
    for t in range(t_steps):

        NODES = list(g.nodes())
        rnd.shuffle(NODES)

        """ initiates therapy program at predetermined start """
        if t == t_start:
            p_treat = Params['treatment'][run_index]
            treatment(g, p_treat, random=rand_treat)

        """ initiates SIS start """
        if t == sis_start:
            safe_injection(g, localized=local_SIS)
            cut_ties = []
            for each in g:
                if type(g.nodes[each]['sis']) == tuple:
                    cut_ties.extend(list(g.neighbors(each)))
#        num_in_treat.append(len([x for x in g if g.nodes[x]['treat'] == 1]))

        # iterate through all nodes
        for node in NODES:

            # store attributes of the node
            Node = g.nodes[node]

            # possibility of recovery by clearers
            recovered = recovery(Node)

            if recovered:
                time_loc.append((t, 'r', Node['pos']))

            """ if user decides to inject """
            if rnd.random() < p_use:
                
                if type(Node['sis']) == tuple:
                    Node['p_share'] = 0

                # if user decides to share needles
                if rnd.random() < Node['p_share']:

                    # run through the cases:
                    infections = interaction_cases(node, g)

                    # record the infections:
                    for each in infections:
                        pos = g.nodes[each]['pos']
                        time_loc.append((t, 'i', pos))

                # if user decides not to share
                else:
                    a = Node['needles'][0]  # num of clean needles for node
                    b = Node['needles'][1]  # num of infected needles for node

                    # if user is infected and a clean needle is selected
                    if Node['health'] == 1 and rnd.random() < a / (a + b):

                        Node['needles'][0] -= 1  # num of clean needles for user x
                        Node['needles'][1] += 1  # num of infected needles for user x

                    # if user is uninfected and infected needle is selected
                    elif Node['health'] == 0 and rnd.random() < b / (a + b):

                        Node['health'] = 1  # user becomes infected

                        # record the infection
                        time_loc.append((t, 'i', Node['pos']))

            """ infection on needles can decay away with probability p_decay """
            if rnd.random() < p_decay and Node['needles'][1] > 0:

                Node['needles'][0] += 1
                Node['needles'][1] -= 1

            """ Needle-exchange program """
            if Node['nep'] == 1:
                if rnd.random() < nep_per_day:
                    num_inf_needles = Node['needles'][1]
                    Node['needles'][1] -= num_inf_needles
                    Node['needles'][0] += num_inf_needles

            """ Treatment """
            if Node['treat'] == 1 and Node['health'] == 1:
                Node['treat_time'] += 1
                if rnd.random() < p_adhere:
                    Node['p_share'] = 0  # <-- assumes person isn't sharing needles while in treatment
                if Node['treat_time'] > 3.5 * 30:  # <-- takes about 3-4 months to treat HCV
                    Node['health'] = 0             # <-- cured after 3.5 months
                    Node['p_share'] = p_share      # <-- resumes sharing w/ others
                    Node['treat_time'] = 0         # <-- treat time starts over

                    time_loc.append((t, 'r', Node['pos']))  # <-- record the recovery

            """ nodes relocate """
            if rnd.random() < Params['relocate']:
                old_pos = Node['pos']
                new_pos = relocate(node, g, b_matrix)
                if Node['health'] == 1:  # <-- an infection relocates
                    time_loc.append((t, 'r', old_pos))
                    time_loc.append((t, 'i', new_pos))
                if g.nodes[node]['treat'] == 1:   # <-- find new patient
                    new_treat = find_new_patient(g, node, p_treat, random=rand_treat)
                if type(g.nodes[node]['sis']) == tuple:  # <-- if in SIS
                    SIS_loc = g.nodes[node]['sis']
                    g.nodes[node]['sis'] = 'not in SIS'
                    g.nodes[node]['p_share'] = p_share  # after relocating, node resumes sharing w/ others
                    find_SIS_participant(g, SIS_loc)

                """ death and birth: """
            if rnd.random() < p_death:

                """ if a death occurs, exactly one birth will occur to keep the population fixed """

                """
                we generate a new node number, 'new_node', which can also be
                used to identify new nodes in the order of their appearance
                """
                new_node = 1 + max(g)

                # new_node will take on several of the old node's attributes (used g.nodes to get up-to-date info)
                immunity = g.nodes[node]['immunity']
                group = g.nodes[node]['group']
                in_nep = g.nodes[node]['nep']
                in_treat = g.nodes[node]['treat']
                in_SIS = g.nodes[node]['sis']

                # if old node was infected, need to record it as a recovery in time_loc
                if Node['health'] == 1:
                    time_loc.append((t, 'r', Node['pos']))

                # if old node was in treatment program, will find a new patient
                if in_treat == 1:
                    new_treat = find_new_patient(g, node, p_treat, random=rand_treat, find_someone_else=True)

                # if old node participated in an SIS, find a new participant
                if type(in_SIS) == tuple:
                    find_SIS_participant(g, in_SIS)

                # old node leaves network
                g.remove_node(node)

                # we randomly select a new location:
                new_pos = (rnd.random(), rnd.random())

                # adding new_node to the graph
                g.add_node(new_node,
                           pos=new_pos,
                           group=group,
                           needles=[rnd.randint(1, 3), 0],
                           health=0,  # <-- have chosen 0 to mean not infected
                           immunity=immunity,
                           p_use=p_use,
                           p_share=p_share,
                           treat_time=0,
                           nep=in_nep,
                           sis='not in SIS',
                           treat=0)

                # --- this forms new edges based on the modified geometric block model------------------------ #

                # max. distance nodes can connect within
                max_radius = max(b_matrix[0] + b_matrix[1])

                # dictionary of all node positions
                loc = {x: nodes[x]['pos'] for x in g}

                # forms connections if radii overlap
                for nhbr in loc:
                    d = dist(new_pos, loc[nhbr])  # <-- distance between new node and old nodes
                    if 0 < d < 2 * max_radius:    # <-- checks if nodes are close enough to connect (and distinct)
                        grp1 = g.nodes[new_node]['group']  # <-- 'node' group
                        grp2 = g.nodes[nhbr]['group']      # <-- 'nhbr' group
                        r1 = b_matrix[grp1][grp2]  # <-- distance grp1 folks will go to link up with grp2 folks
                        r2 = b_matrix[grp2][grp1]  # <-- distance grp2 folks will go to link up with grp1 folks
                        overlap = overlap_fraction(d, r1, r2)
                        if overlap > 0:  # <-- if radii overlap
                            g.add_edge(new_node, nhbr, weight=overlap)

                # ------------------------------------------------------------------------------------------- #

        # updates the num of clean needles
        clean_ndls = np.append(clean_ndls, sum([nodes[nds]['needles'][0] for nds in nodes]))

        # updates the num of infected needles
        infct_ndls = np.append(infct_ndls, sum([nodes[nds]['needles'][1] for nds in nodes]))

        """ storing updated infected populations """
        a_clr, a_nclr, b_clr, b_nclr = pop_count(nodes)

        A_clr[t+1], A_nclr[t+1] = a_clr, a_nclr
        B_clr[t+1], B_nclr[t+1] = b_clr, b_nclr

        # records the degree distribution after each time step:
        deg_dist.append([len(list(g.neighbors(x))) for x in g])

        # records SIS numbers:
#        for each in g:
#            sis_tup = g.nodes[each]['sis']
#            if type(sis_tup) == tuple:
#                num_in_SIS[sis_tup][t] += 1

        # records the components and whether node is in SIS:
        # in_SIS = [x for x in g if type(g.nodes[x]['sis']) == tuple]
        # sis_subgraph = g.subgraph(in_SIS)
        # components.append(list(nx.connected_components(sis_subgraph)))
        # temp_sis_dict = {(i, j): [] for i in range(5) for j in range(5)}
        # for each in in_SIS:
        #     sis_loc = g.nodes[each]['sis']
        #     temp_sis_dict[sis_loc].append(each)
        #
        # sis_dicts.append(temp_sis_dict)
            
#        num_in_treat = 0
#        num_in_sis = 0
#        num_in_both = 0
#        for x in g:
#            a, b = 0, 0
#            if g.nodes[x]['treat'] == 1:
#                a = 1
#            if type(g.nodes[x]['sis']) == tuple:
#                b = 1
#            if a + b > 1:
#                num_in_both += 1
#            num_in_treat += a
#            num_in_sis += b

#        total_in_treat.append(num_in_treat)
#        total_in_sis.append(num_in_sis)
#        total_in_both.append(num_in_both)

    As = [A_clr, A_nclr]
    Bs = [B_clr, B_nclr]

    AB = [As, Bs]

    Needles = [clean_ndls, infct_ndls]

#    aux = [total_in_treat, total_in_sis, total_in_both]  # looking at how much overlap between SIS and treatment
#    aux = cut_ties  # for when SIS is used
    aux = 'n/a'

    return (AB, Needles, deg_dist, aux)

# --------------------------------------------------------------------------------------------------------- #

