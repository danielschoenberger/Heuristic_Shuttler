from Circles_v3 import *
import random


def create_starting_config(perc, graph, seed=None):
    if seed != None:
        random.seed(seed)
        random_starting_traps = random.sample(range(n_of_traps), (math.ceil(perc*n_of_traps)))
        starting_traps = []
        for trap in random_starting_traps:
            starting_traps.append([edges for edges in graph.edges() if graph.get_edge_data(edges[0], edges[1])['edge_type'] == 'trap'][trap])
    else:
        starting_traps = [edges for edges in graph.edges() if graph.get_edge_data(edges[0], edges[1])['edge_type'] == 'trap'][:math.ceil(perc*n_of_traps)+1]
    number_of_registers = len(starting_traps)

    # place ions onto traps (ion0 on starting_trap0)
    ion_chains = {}
    for ion, idc in enumerate(starting_traps):
        ion_chains[ion] = idc

    return ion_chains, number_of_registers


archs = [
[2, 2, 1, 5],
[2, 2, 1, 11],
[2, 2, 1, 19],
[2, 2, 1, 29],

[3, 3, 1, 1],
[4, 4, 1, 1],
[5, 5, 1, 1],
[6, 6, 1, 1]
]
seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
perc = 0.5
results = {}

for j, arch in enumerate(archs):
    timestep_arr = []
    for k, seed in enumerate(seeds):
        m, n, v, h = arch
        # create dummy graph
        graph = GraphCreator(m, n, v, h).get_graph()
        n_of_traps = len([trap for trap in graph.edges() if graph.get_edge_data(trap[0], trap[1])['edge_type'] == 'trap'])

        ion_chains, number_of_registers = create_starting_config(perc, graph, seed=seed)

        sequence = list(range(len(ion_chains)))

        max_timesteps = 500

        Mem1 = MemoryZone(m, n, v, h, ion_chains, sequence, max_timesteps)

        timestep = 0
        while timestep < max_timesteps:
            
            ### move all chains until they need to rotate because they are at a junction (move if path is free and not at a junction)
            need_rotate = [False] * len(sequence)
            while sum(need_rotate) < len(sequence):
                for i, rotate_chain in enumerate(sequence):
                    edge_idc = Mem1.ion_chains[rotate_chain]
                    next_edge = Mem1.find_next_edge(edge_idc)
                    
                    state_edges_idx = Mem1.get_state_idxs()
                    if Mem1.have_common_junction_node(edge_idc, next_edge) == False and get_idx_from_idc(Mem1.idc_dict, next_edge) not in state_edges_idx:
                        # update ion chains
                        Mem1.ion_chains[rotate_chain] = next_edge
                    else:
                        need_rotate[i] = True        

            ### calc distance to entry for all chains and determine which chains can rotate
            path_length_sequence = {}
            move_sequence = []
            for i, rotate_chain in enumerate(sequence):
                edge_idc = Mem1.ion_chains[rotate_chain]
                path_to_go = nx.shortest_path(Mem1.graph, edge_idc[0], Mem1.graph_creator.processing_zone, lambda edge0, edge1, edge_attr_dict : (edge_attr_dict['edge_type'] == 'entry')*1e8+1)
                path_length_sequence[rotate_chain] = len(path_to_go)
                
                # first chain always moves (also always move if in entry or exit edge)
                if i == 0:
                    move_sequence.append(rotate_chain)
                # wenn path von chain größer als alle vorherigen ist -> sum = länge move_sequence -> lassen, sonst -> remove
                elif sum(np.array([path_length_sequence[rotate_chain]]*len(move_sequence))>np.array([path_length_sequence[chain] for chain in move_sequence])) == len(move_sequence):
                    move_sequence.append(rotate_chain)

            all_circles = {}
            ### create circles for all chains in move_sequence (dictionary with chain as key and circle_idcs as value)
            for rotate_chain in move_sequence:
                edge_idc = Mem1.ion_chains[rotate_chain]
                next_edge = Mem1.find_next_edge(edge_idc)

                # make edge_idc and next_edge consistent
                edge_idc, next_edge = Mem1.find_ordered_edges(edge_idc, next_edge)

                # if next_edge (now over junction) is free -> circle not needed -> only edge_idc + next_edge
                # TODO checken ob in nächstem Schritt frei? -> Müsste hier schon alles nacheinander fahren (also neue Logik einbauen)
                if Mem1.check_if_edge_is_filled(next_edge) == False:
                    all_circles[rotate_chain] = [edge_idc] + [next_edge]
                elif get_idx_from_idc(Mem1.idc_dict, next_edge) == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.exit_edge) or get_idx_from_idc(Mem1.idc_dict, next_edge) == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.entry_edge):
                    all_circles[rotate_chain] = [edge_idc] + [next_edge] + Mem1.create_outer_circle(edge_idc, next_edge, path_over_pz=True, include_path_to_exit_edge=False) + [edge_idc]
                else:
                    all_circles[rotate_chain] = [edge_idc] + [next_edge] + Mem1.create_outer_circle(edge_idc, next_edge) + [edge_idc]   # edge idc is added twice to close circle
            
            # find circles that can move while first seq ion is moving
            nonfree_circles, free_circle_combs = Mem1.find_circle_idxs_sharing_nodes(all_circles)
            free_circle_seq_idxs = [move_sequence[0]]
            for seq_circ in move_sequence[1:]:
                nonfree = False
                for mov_circ in free_circle_seq_idxs:
                    if (seq_circ, mov_circ) in nonfree_circles or (mov_circ, seq_circ) in nonfree_circles:
                        nonfree = True
                        break
                if nonfree == False:
                    free_circle_seq_idxs.append(seq_circ)

            # need circles given in indexes for rotate function
            free_circle_idxs = {}
            for seq_idx in free_circle_seq_idxs:
                free_circle_idxs[seq_idx] = [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in all_circles[seq_idx]]
                # rotate chains
                Mem1.rotate(free_circle_idxs[seq_idx])

            timestep += 1
            if get_idx_from_idc(Mem1.idc_dict, Mem1.ion_chains[sequence[0]]) == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.entry_edge):
                print('\ntime step: %s, chain %s is at entry edge' %(timestep, sequence[0]))
                if len(sequence) == 1:
                    timestep_arr.append(timestep)
                    print('full circuit executed, resulting time steps: %s' %timestep)
                    break
                sequence = sequence[1:]
            
    timestep_mean = np.mean(timestep_arr)
    results[j] = timestep_mean

print(results)

    
