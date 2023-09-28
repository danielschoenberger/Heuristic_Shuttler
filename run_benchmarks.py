import contextlib
import math
import random
import time

import networkx as nx
import numpy as np

from Cycles import GraphCreator, MemoryZone, get_idx_from_idc, get_path_to_node


def create_starting_config(perc, graph, seed=None):
    if seed is not None:
        random.seed(seed)
        random_starting_traps = random.sample(range(n_of_traps), (math.ceil(perc * n_of_traps)))
        starting_traps = []
        traps = [edges for edges in graph.edges() if graph.get_edge_data(edges[0], edges[1])["edge_type"] == "trap"]
        for trap in random_starting_traps:
            starting_traps.append(traps[trap])
    else:
        starting_traps = [
            edges for edges in graph.edges() if graph.get_edge_data(edges[0], edges[1])["edge_type"] == "trap"
        ][: math.ceil(perc * n_of_traps) + 1]
    number_of_registers = len(starting_traps)

    # place ions onto traps (ion0 on starting_trap0)
    ion_chains = {}
    for ion, idc in enumerate(starting_traps):
        ion_chains[ion] = idc

    return ion_chains, number_of_registers


archs = [[3, 3, 1, 1]]
seeds = [1]  # , 2, 3, 4, 5, 6, 7, 8, 9, 10]
perc = 0.5
results = {}
cpu_time_results = {}
start_time_all = time.time()


for j, arch in enumerate(archs):
    timestep_arr = []
    cpu_time_arr = []
    for _k, seed in enumerate(seeds):
        start_time = time.time()
        m, n, v, h = arch
        # create dummy graph
        graph = GraphCreator(m, n, v, h).get_graph()
        n_of_traps = len(
            [trap for trap in graph.edges() if graph.get_edge_data(trap[0], trap[1])["edge_type"] == "trap"]
        )

        ion_chains, number_of_registers = create_starting_config(perc, graph, seed=seed)

        sequence = [0, 1, 3, 2]  # TODO no immediately repeating seq elements are possible, e.g. [0, 1, 1, 2]
        # sequence = list(range(number_of_registers))
        # 2-qubit sequence -> [2] means, that the 3rd ion in the sequence (ion with id 2) is a 2-qubit gate
        two_qubit_sequence = [0, 2]
        assert (
            two_qubit_sequence[-1] != len(sequence) - 1
        ), "2-qubit sequence is not valid (last element can not be 2-qubit gate)"
        for i, elem in enumerate(two_qubit_sequence[:-1]):
            assert (
                two_qubit_sequence[i + 1] > elem + 1
            ), "2-qubit sequence is not valid (can only be 2-qubit gate if next element is at least 2 steps away -> can not do two 2-qubit gate on same ion)"

        seq_element_counter = 0

        print(f"arch: {arch}, seed: {seed}, registers: {number_of_registers}\n")
        max_timesteps = 50000

        Mem1 = MemoryZone(m, n, v, h, ion_chains, sequence, max_timesteps)
        timestep = 0
        print("time step: %s" % timestep)

        seq_ion_was_at_entry = 0
        next_seq_ion_in_exit = 0

        seq_ion_was_at_entry = False
        timestep = 1
        while timestep < max_timesteps:
            ########### PREPROCESSING ###########
            ### move all chains until they need to rotate because they are at a junction (move if path is free and not at a junction)
            need_rotate = [False] * len(sequence)
            while sum(need_rotate) < len(sequence):
                for i, rotate_chain in enumerate(sequence):
                    edge_idc = Mem1.ion_chains[rotate_chain]
                    next_edge = Mem1.find_next_edge(edge_idc)

                    state_edges_idx = Mem1.get_state_idxs()
                    if (
                        Mem1.have_common_junction_node(edge_idc, next_edge) is False
                        and get_idx_from_idc(Mem1.idc_dict, next_edge) not in state_edges_idx
                    ):
                        # update ion chains
                        Mem1.ion_chains[rotate_chain] = next_edge
                    else:
                        need_rotate[i] = True
            print("time step: %s, moved all chains as far as possible" % timestep)
            # Mem1.graph_creator.plot_state([get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()], show_plot=False)
            # plt.show()

            ######### MOVE SEQUENCE #########
            ### calc distance to entry for all chains and determine which chains can rotate
            path_length_sequence = {}
            move_sequence = []
            for i, rotate_chain in enumerate(sequence):
                edge_idc = Mem1.ion_chains[rotate_chain]
                path_to_go = nx.shortest_path(
                    Mem1.graph,
                    edge_idc[0],
                    Mem1.graph_creator.processing_zone,
                    lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "entry") * 1e8 + 1,
                )
                path_length_sequence[rotate_chain] = len(path_to_go)

                # first chain always moves (also always move if in entry or exit edge)
                # or: wenn path von chain größer als alle vorherigen ist -> sum = länge move_sequence -> lassen, sonst -> remove
                if i == 0 or sum(
                    np.array([path_length_sequence[rotate_chain]] * len(move_sequence))
                    > np.array([path_length_sequence[chain] for chain in move_sequence])
                ) == len(move_sequence):
                    move_sequence.append(rotate_chain)

                # change path of chain if it is in the entry edge (to 0, so it doesn't block others when it is the first one)
                if get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ):
                    path_length_sequence[rotate_chain] = 0
                    # if this chain is already present in move_sequence, remove it and add it again at the beginning
                    with contextlib.suppress(Exception):
                        move_sequence.remove(rotate_chain)

                    move_sequence = [rotate_chain, *move_sequence]

            ######### ROTATE CHAINS #########
            all_circles = {}
            ### create circles for all chains in move_sequence (dictionary with chain as key and circle_idcs as value)
            for rotate_chain in move_sequence:
                edge_idc = Mem1.ion_chains[rotate_chain]
                next_edge = Mem1.find_next_edge(edge_idc)

                # make edge_idc and next_edge consistent
                edge_idc, next_edge = Mem1.find_ordered_edges(edge_idc, next_edge)

                ### if next_edge (now over junction) is free -> circle not needed -> only edge_idc + next_edge
                # real TODO checken ob in nächstem Schritt frei? -> Müsste hier schon alles nacheinander fahren (also neue Logik einbauen)
                if not Mem1.check_if_edge_is_filled(next_edge):
                    all_circles[rotate_chain] = [edge_idc, next_edge]

                ### move into exit (if next edge is exit)
                elif get_idx_from_idc(Mem1.idc_dict, next_edge) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.exit_edge
                ):
                    # now here same logic as for entry edge below
                    top_left_free_edge_idc = Mem1.bfs_free_edge((0, 0))
                    path0 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[0],
                        exclude_exit=True,
                        exclude_entry=False,
                    )
                    path1 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[1],
                        exclude_exit=True,
                        exclude_entry=False,
                    )
                    # find "circle" (is only a path here) to free edge top left
                    # include free edge in path -> find node with longer path -> is outer node -> this path also includes the free edge itself
                    if len(path1) > len(path0):
                        path = [*path0, (top_left_free_edge_idc[0], top_left_free_edge_idc[1])]
                    else:
                        path = [*path1, (top_left_free_edge_idc[1], top_left_free_edge_idc[0])]

                    all_circles[rotate_chain] = [edge_idc, next_edge, *path]
                    assert next_edge == (Mem1.graph_creator.exit, Mem1.graph_creator.processing_zone)

                ### move into entry (if next edge is entry)
                elif get_idx_from_idc(Mem1.idc_dict, next_edge) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ):
                    top_left_free_edge_idc = Mem1.bfs_free_edge((0, 0))
                    path0 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[0],
                        exclude_exit=True,
                        exclude_entry=False,
                    )
                    path1 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[1],
                        exclude_exit=True,
                        exclude_entry=False,
                    )
                    # find "circle" (is only a path here) to free edge top left
                    # include free edge in path -> find node with longer path -> is outer node -> this path also includes the free edge itself
                    if len(path1) > len(path0):
                        path = [*path0, (top_left_free_edge_idc[0], top_left_free_edge_idc[1])]
                    else:
                        path = [*path1, (top_left_free_edge_idc[1], top_left_free_edge_idc[0])]

                    all_circles[rotate_chain] = [edge_idc, *path]

                ### move out of entry (if edge is entry)
                elif get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ):
                    # same logic as for entry edge above
                    top_left_free_edge_idc = Mem1.bfs_free_edge((0, 0))
                    path0 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[0],
                        exclude_exit=True,
                        exclude_entry=False,
                    )
                    path1 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[1],
                        exclude_exit=True,
                        exclude_entry=False,
                    )
                    # find "circle" (is only a path here) to free edge top left
                    # include free edge in path -> find node with longer path -> is outer node -> this path also includes the free edge itself
                    if len(path1) > len(path0):
                        path = [*path0, (top_left_free_edge_idc[0], top_left_free_edge_idc[1])]
                    else:
                        path = [*path1, (top_left_free_edge_idc[1], top_left_free_edge_idc[0])]

                    all_circles[rotate_chain] = path

                ### move within memory zone (circles)
                else:
                    all_circles[rotate_chain] = [
                        edge_idc,
                        next_edge,
                        *Mem1.create_outer_circle(edge_idc, next_edge),
                        edge_idc,
                    ]  # edge idc is added twice to close circle

            # if pz node is in circles -> combine them to one circle
            first_circle = -1
            # need to loop over keys instead of dict, because dict changes size during loop (keys should be rotate chains)
            all_circles_keys = list(all_circles.keys())

            for circle in all_circles_keys:
                inner_if_broken = False  # flag to track inner if clause break
                for edge in all_circles[circle]:
                    if Mem1.graph_creator.processing_zone in edge:
                        if first_circle == -1:
                            first_circle = circle
                            inner_if_broken = True
                            break
                        else:
                            # here new function needed
                            if all(edge in all_circles[circle] for edge in all_circles[first_circle]):
                                all_circles[first_circle] = all_circles[circle].copy()
                                del all_circles[circle]
                                move_sequence.remove(circle)
                            elif all(edge in all_circles[first_circle] for edge in all_circles[circle]):
                                del all_circles[circle]
                                move_sequence.remove(circle)
                            else:
                                all_circles[first_circle] = Mem1.combine_paths_over_pz(
                                    all_circles[first_circle], all_circles[circle]
                                )
                                del all_circles[circle]
                                move_sequence.remove(circle)
                            inner_if_broken = True
                            break
                    if inner_if_broken:
                        break

            # find circles that can move while first seq ion is moving
            nonfree_circles, free_circle_combs = Mem1.find_nonfree_and_free_circle_idxs(all_circles)
            free_circle_seq_idxs = [move_sequence[0]]
            for seq_circ in move_sequence[1:]:
                nonfree = False
                for mov_circ in free_circle_seq_idxs:
                    if (seq_circ, mov_circ) in nonfree_circles or (mov_circ, seq_circ) in nonfree_circles:
                        nonfree = True
                        break
                if nonfree is False:
                    free_circle_seq_idxs.append(seq_circ)

            # need circles given in indexes for rotate function
            free_circle_idxs = {}
            for seq_idx in free_circle_seq_idxs:
                free_circle_idxs[seq_idx] = [
                    get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in all_circles[seq_idx]
                ]
                # rotate chains
                Mem1.rotate(free_circle_idxs[seq_idx])
            # Mem1.graph_creator.plot_state([get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()], show_plot=True)
            timestep += 1

            # if seq ion was at entry last timestep -> is now back in memory -> remove from sequence
            if seq_ion_was_at_entry is True:
                # seq_ion_was_at_entry is only True if seq ion is now not in entry anymore (for 2-qubit gates it has to wait in entry):
                if seq_element_counter == two_qubit_sequence[0] and next_seq_ion_in_exit < 2:
                    pass
                else:
                    if len(sequence) == 1:
                        cpu_time_arr.append(time.time() - start_time)
                        timestep_arr.append(timestep)
                        print("full circuit executed, resulting time steps: %s" % timestep)
                        break
                    last_seq_element = sequence[0]
                    sequence = sequence[1:]
                    seq_ion_was_at_entry = False
                    seq_element_counter += 1

            # if seq ion is at entry -> has to move out of entry -> then remove from sequence with if clause above in next iteration
            if get_idx_from_idc(Mem1.idc_dict, Mem1.ion_chains[sequence[0]]) == get_idx_from_idc(
                Mem1.idc_dict, Mem1.graph_creator.entry_edge
            ):
                seq_ion_was_at_entry = True
                print(f"\ntime step: {timestep}, chain {sequence[0]} is at entry edge")

    timestep_mean = np.mean(timestep_arr)
    cpu_time_mean = np.mean(cpu_time_arr)
    results[j] = timestep_mean
    cpu_time_results[j] = cpu_time_mean

print("\n archs: \n", dict(enumerate(archs)))
print("results: \n", results)
print("cpu time results: \n", cpu_time_results)
print("time all: \n", time.time() - start_time_all)
