import contextlib
import math
import random
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import networkx as nx
import numpy as np

from Cycles import GraphCreator, MemoryZone, get_idx_from_idc, get_path_to_node
from get_sequence import parse_qasm, write_qft_to_file


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


archs = [[3, 3, 1, 1]]  # , [5, 5, 1, 1], [6, 6, 1, 1]]#, [20, 20, 1, 1], [5, 5, 10, 10]]#[5, 5, 1, 1],
seeds = [2]  # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
perc = 0.5
results = {}
cpu_time_results = {}
start_time_all = time.time()
show_plot = False
save_plot = True
# assert either show_plot or save_plot is True, "Either show_plot or save_plot must be True"
assert not (show_plot and save_plot), "Only one of show_plot or save_plot can be True"

for j, arch in enumerate(archs):
    timestep_arr = []
    cpu_time_arr = []
    for _k, seed in enumerate(seeds):
        start_time = time.time()

        if save_plot:
            # Create a folder for each run with a timestamp (plot widget)
            run_folder = Path(f'plots/run_{datetime.now().strftime("%Y%m%d_%H%M%S")}')
            run_folder.mkdir(parents=True, exist_ok=True)
        else:
            run_folder = ""

        m, n, v, h = arch
        # create dummy graph
        graph = GraphCreator(m, n, v, h).get_graph()
        n_of_traps = len(
            [trap for trap in graph.edges() if graph.get_edge_data(trap[0], trap[1])["edge_type"] == "trap"]
        )

        ion_chains, number_of_registers = create_starting_config(perc, graph, seed=seed)

        # generate sequence and two-qubit sequence
        N = number_of_registers
        filename = "QASM_files/qft_%squbits.qasm" % N
        # /Users/danielschonberger/Desktop/Heuristic_Github/
        write_qft_to_file(N, filename)
        print(f"QFT for {N} qubits written to {filename}\n")

        seq = parse_qasm(filename)
        flat_seq = [item for sublist in seq for item in sublist]
        # change last 3 seq elements
        seq = seq[:-3] + [(4, 5), (5,)]
        print(seq)

        two_qubit_seq = []
        i = 0
        s = 0
        while i < len(seq):
            if len(seq[i]) == 2:
                two_qubit_seq.append(i + s)
                s += 1
            i += 1

        # TODO no immediately repeating seq elements are possible, e.g. [0, 1, 1, 2]
        sequence = flat_seq  # [0, 1, 0, 2, 0, 3, 0, 4]
        # sequence = list(range(number_of_registers))

        # chain in entry is first sequence element at the start
        chain_in_entry = sequence[0]

        # # 2-qubit sequence -> [2] means, that the 3rd ion in the sequence (ion with id 2) is a 2-qubit gate
        # two_qubit_sequence = [0, 2, 4, 6]
        # # two_qubit_sequence = [1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 16, 17, 19]
        two_qubit_sequence = two_qubit_seq  # [1, 3, 5, 7, 9, 12, 14, 16, 18, 21, 23, 25, 28, 30, 33]#[]
        two_qubit_sequence.append(
            -1
        )  # -1 at the end, so two-qubit sequence is not empty after all 2-qubit gates have been processed

        # assert (
        #    two_qubit_sequence[-2] != len(sequence) - 1
        # ), "2-qubit sequence is not valid (last element can not be 2-qubit gate)"

        # unique sequence is sequence without repeating elements (for move_sequence and 2-qubit gates)
        unique_sequence = []
        for seq_elem in sequence:
            if seq_elem not in unique_sequence:
                unique_sequence.append(seq_elem)

        seq_element_counter = 0

        print(f"arch: {arch}, seed: {seed}, registers: {number_of_registers}\n")
        max_timesteps = 50000
        # time_2qubit_gate only works for 2-qubit gates -> need also for 1-qubit gates
        time_2qubit_gate = 3

        Mem1 = MemoryZone(m, n, v, h, ion_chains, sequence, max_timesteps, time_2qubit_gate=time_2qubit_gate)
        timestep = 0
        print("time step: %s" % timestep)

        # Save the current plot with a meaningful filename (plot widget)
        plot_filename = Path(run_folder) / f"plot_{1:03d}.png"

        Mem1.graph_creator.plot_state(
            [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
            labels=[
                "time step %s" % timestep,
                "next 3 seq elem: %s" % seq[seq_element_counter : seq_element_counter + 3],
            ],
            show_plot=show_plot,
            save_plot=save_plot,
            filename=[plot_filename if save_plot else None][0],
        )

        seq_ion_was_at_entry = 0
        next_seq_ion_in_exit = 0

        seq_ion_was_at_entry = False
        # timestep = 1
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
            print(state_edges_idx, "s", Mem1.ion_chains)
            print("time step: %s, moved all chains as far as possible" % timestep)
            print("sequence: ", sequence)
            # Mem1.graph_creator.plot_state([get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()], show_plot=False)
            # plt.show()

            ######### MOVE SEQUENCE #########
            ### calc distance to entry for all chains and determine which chains can rotate
            path_length_sequence = {}
            move_sequence = []
            for i, rotate_chain in enumerate(unique_sequence):
                # if 2-qubit gate -> check if second ion is in exit or in parking edge
                # maybe not in parking edge -> move to exit from parking edge to do 2 qubit gate?
                if (
                    seq_element_counter == two_qubit_sequence[0]
                    and rotate_chain == sequence[1]
                    and (
                        get_idx_from_idc(Mem1.idc_dict, Mem1.ion_chains[rotate_chain])
                        == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.exit_edge)
                        or get_idx_from_idc(Mem1.idc_dict, Mem1.ion_chains[rotate_chain])
                        == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.parking_edge)
                    )
                ):
                    next_seq_ion_in_exit += 1

                edge_idc = Mem1.ion_chains[rotate_chain]
                path_to_go = nx.shortest_path(
                    Mem1.graph,
                    edge_idc[0],
                    Mem1.graph_creator.processing_zone,
                    lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "entry") * 1e8 + 1,
                )
                path_length_sequence[rotate_chain] = len(path_to_go)

                print(
                    "here 0",
                    rotate_chain,
                    path_to_go,
                    path_length_sequence[rotate_chain],
                    sum(
                        np.array([path_length_sequence[rotate_chain]] * len(move_sequence))
                        > np.array([path_length_sequence[chain] for chain in move_sequence])
                    )
                    == len(move_sequence),
                    (
                        np.array([path_length_sequence[rotate_chain]] * len(move_sequence))
                        > np.array([path_length_sequence[chain] for chain in move_sequence]),
                        [path_length_sequence[rotate_chain]],
                        [path_length_sequence[chain] for chain in move_sequence],
                    ),
                )

                # first chain always moves (also always move if in entry or exit edge)
                # or: wenn path von chain größer als alle vorherigen ist -> sum = länge move_sequence -> lassen, sonst -> remove
                if i == 0 or sum(
                    np.array([path_length_sequence[rotate_chain]] * len(move_sequence))
                    > np.array([path_length_sequence[chain] for chain in move_sequence])
                ) == len(move_sequence):
                    move_sequence.append(rotate_chain)
                print("move_sequence", move_sequence)
                # TODO funktioniert immer noch nicht + erstes mal 0 in entry reinfahren -> 1 könnte eig instant nach aber wartet?

                # change path of chain if it is in entry edge (to 0, so it doesn't block others when it is the first one)
                # but only if it is not the next or next after that in the sequence (needed again)
                if get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ) and (
                    rotate_chain not in sequence[1:2]
                    or (
                        rotate_chain in sequence[1:2]
                        and seq_element_counter == two_qubit_sequence[0]
                        and next_seq_ion_in_exit < Mem1.time_2qubit_gate
                    )
                ):
                    # problem fixed with extra if clause at the end of "move into exit" below (QFT specific)
                    path_length_sequence[rotate_chain] = 0
                    # if this chain is already present in move_sequence, remove it and add it again at the beginning
                    with contextlib.suppress(Exception):
                        move_sequence.remove(rotate_chain)

                    move_sequence = [rotate_chain, *move_sequence]

            ######### ROTATE CHAINS #########
            # first find chain in entry edge (if present -> then also present in move_sequence)
            for rotate_chain in move_sequence:
                if get_idx_from_idc(Mem1.idc_dict, Mem1.ion_chains[rotate_chain]) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ):
                    chain_in_entry = rotate_chain

            ### create circles for all chains in move_sequence (dictionary with chain as key and circle_idcs as value)
            all_circles = {}
            for rotate_chain in move_sequence:
                edge_idc = Mem1.ion_chains[rotate_chain]
                if (
                    rotate_chain in sequence[1:]
                ):  # if chain is needed again (is present in rest of sequence) -> move (only out of entry) towards exit instead of top left
                    next_edge = Mem1.find_next_edge(edge_idc, towards="exit")
                    # print('next edge towards exit', rotate_chain, sequence[1:])
                else:
                    next_edge = Mem1.find_next_edge(edge_idc, towards=(0, 0))
                    # print('next edge towards top left', rotate_chain, sequence[1:])

                # make edge_idc and next_edge consistent
                edge_idc, next_edge = Mem1.find_ordered_edges(edge_idc, next_edge)

                ### if next_edge (now over junction) is free -> circle not needed -> only edge_idc + next_edge
                # real TODO checken ob in nächstem Schritt frei? -> Müsste hier schon alles nacheinander fahren (also neue Logik einbauen)
                if not Mem1.check_if_edge_is_filled(next_edge):
                    all_circles[rotate_chain] = [edge_idc, next_edge]

                    # new try
                    # if rotate_chain == sequence[0] and get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                    # Mem1.idc_dict, Mem1.graph_creator.entry_edge
                    # ):
                    #     all_circles[rotate_chain] = [edge_idc, Mem1.graph_creator.parking_edge]
                    #     continue

                    # new park logic
                    if rotate_chain in sequence[1:3] and get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                        Mem1.idc_dict, Mem1.graph_creator.entry_edge
                    ):
                        all_circles[rotate_chain] = [edge_idc, Mem1.graph_creator.parking_edge]

                    # for two-qubit gates: leave ion in entry until gate is processed (other ion was in exit long enough)
                    # entry
                    # if 2-qubit gate and in entry -> need to wait even if next edge is free
                    print("\nnext edge free")
                    print("seq_element_counter: %s" % seq_element_counter)
                    print("next_seq_ion_in_exit: %s" % next_seq_ion_in_exit, "\n")
                    if (
                        get_idx_from_idc(Mem1.idc_dict, edge_idc)
                        == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.entry_edge)
                        and seq_element_counter == two_qubit_sequence[0]
                    ):
                        if next_seq_ion_in_exit < Mem1.time_2qubit_gate:
                            all_circles[rotate_chain] = [edge_idc, edge_idc]
                        else:
                            next_seq_ion_in_exit = 0
                            two_qubit_sequence.pop(0)

                ### move into exit (if next edge is exit)
                elif get_idx_from_idc(Mem1.idc_dict, next_edge) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.exit_edge
                ):
                    if (
                        chain_in_entry in sequence[1:]
                    ):  # if chain in entry is needed again (is present in rest of sequence) -> move (only out of entry) towards exit instead of top left
                        top_left_free_edge_idc = Mem1.bfs_free_edge(Mem1.graph_creator.exit)
                    else:
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

                    # extra clause needed, because if first ion of 2 qubit gate moves out of pz -> path length is 0 in move seq -> doesn't block third ion to enter, even though first ion of 2 qubit gate is needed immediately again (QFT)
                    # summary: only enters pz if it is the next or next next in sequence (only two ions in pz)
                    # if rotate_chain not in sequence[:2]:
                    #    all_circles[rotate_chain] = [edge_idc, edge_idc]

                ### move into entry (if next edge is entry)
                elif get_idx_from_idc(Mem1.idc_dict, next_edge) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ):
                    if (
                        chain_in_entry in sequence[1:]
                    ):  # if chain in entry is needed again (is present in rest of sequence) -> move (only out of entry) towards exit instead of top left
                        top_left_free_edge_idc = Mem1.bfs_free_edge(Mem1.graph_creator.exit)
                    else:
                        top_left_free_edge_idc = Mem1.bfs_free_edge((0, 0))
                    path0 = get_path_to_node(
                        Mem1.graph,
                        Mem1.graph_creator.processing_zone,
                        top_left_free_edge_idc[0],
                        exclude_exit=True,
                        exclude_entry=False
                        # exclude=Mem1.ion_chains[
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

                    # new move out of parking logic
                    if get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                        Mem1.idc_dict, Mem1.graph_creator.parking_edge
                    ):
                        all_circles[rotate_chain] = [edge_idc, next_edge]

                ### move out of entry (if edge is entry)
                elif get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                    Mem1.idc_dict, Mem1.graph_creator.entry_edge
                ):
                    # new try
                    # if rotate_chain == sequence[0]:
                    #     all_circles[rotate_chain] = [edge_idc, Mem1.graph_creator.parking_edge]
                    #     continue

                    # same logic as for entry edge above
                    if (
                        rotate_chain in sequence[1:]
                    ):  # if chain in entry is needed again (is present in rest of sequence) -> move (only out of entry) towards exit instead of top left
                        top_left_free_edge_idc = Mem1.bfs_free_edge(Mem1.graph_creator.exit)
                    else:
                        top_left_free_edge_idc = Mem1.bfs_free_edge((0, 0))
                    # assert chain_in_entry == get_idx_from_idc(Mem1.idc_dict, edge_idc), "chain_in_entry is not equal to edge_idc, which is the chain in entry in this case"
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

                    # new park logic
                    if rotate_chain in sequence[1:3] and get_idx_from_idc(Mem1.idc_dict, edge_idc) == get_idx_from_idc(
                        Mem1.idc_dict, Mem1.graph_creator.entry_edge
                    ):
                        all_circles[rotate_chain] = [edge_idc, Mem1.graph_creator.parking_edge]

                    # if 2-qubit gate and in entry -> need to wait
                    print("\npath through pz")
                    print("seq_element_counter: %s" % seq_element_counter)
                    print("two_qubit_sequence: %s" % two_qubit_sequence)
                    print("next_seq_ion_in_exit: %s" % next_seq_ion_in_exit, "\n")
                    if seq_element_counter == two_qubit_sequence[0]:
                        if next_seq_ion_in_exit < Mem1.time_2qubit_gate:
                            all_circles[rotate_chain] = [edge_idc, edge_idc]
                        else:
                            next_seq_ion_in_exit = 0
                            two_qubit_sequence.pop(0)

                ### move within memory zone (circles)
                else:
                    all_circles[rotate_chain] = [
                        edge_idc,
                        next_edge,
                        *Mem1.create_outer_circle(edge_idc, next_edge),
                        edge_idc,
                    ]  # edge idc is added twice to close circle

            print("all circles before pz combination: %s" % all_circles)

            # if pz node is in circles -> combine them to one circle
            first_circle = -1
            # need to loop over keys instead of dict, because dict changes size during loop (keys should be rotate chains)
            all_circles_keys = list(all_circles.keys())

            for circle in all_circles_keys:
                inner_if_broken = False  # flag to track inner if clause break
                for edge in all_circles[circle]:
                    # if processing zone is in circle -> combine and delete circles (if length <=2 -> was a free edge <- Update: still need it for first circle)
                    if Mem1.graph_creator.processing_zone in edge:  # and len(all_circles[circle]) > 2:
                        # check if it is the first circle in pz -> else continue with combining
                        if first_circle == -1:
                            first_circle = circle
                            inner_if_broken = True
                            break
                        # combine and delete
                        else:
                            # if first circle is move from parking to entry + second circle is move from exit to entry
                            # -> delete second circle, so that first circle can move to entry and out of parking
                            if get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.parking_edge) in [
                                get_idx_from_idc(Mem1.idc_dict, edge) for edge in all_circles[first_circle]
                            ] and get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.entry_edge) in [
                                get_idx_from_idc(Mem1.idc_dict, edge) for edge in all_circles[circle]
                            ]:
                                del all_circles[circle]
                                move_sequence.remove(circle)

                            # already combines that correctly?
                            # if Mem1.graph_creator.parking_node in edge:
                            #     print('parking: ', all_circles[first_circle], all_circles[circle])
                            #     assert get_idx_from_idc(Mem1.idc_dict, all_circles[first_circle][0]) == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.entry_edge), "first circle is not starting from entry edge"
                            #     print(Mem1.combine_paths_over_pz(
                            #         all_circles[first_circle], all_circles[circle]
                            #     ))

                            # check if all edges of first circle are in second circle -> then first circle is now second circle
                            elif all(edge in all_circles[circle] for edge in all_circles[first_circle]):
                                # check if a "wait" move (edge, edge) is in first circle -> happens for 2-qubit gates in entry (ion in entry waits) -> don't combine circles
                                counter_first_circle = Counter(all_circles[first_circle])
                                if 2 in counter_first_circle.values():
                                    break
                                all_circles[first_circle] = all_circles[circle].copy()
                                del all_circles[circle]
                                move_sequence.remove(circle)
                            # check if all edges of second circle are in first circle -> delete second circle
                            elif all(edge in all_circles[first_circle] for edge in all_circles[circle]):
                                counter_first_circle = Counter(all_circles[first_circle])
                                del all_circles[circle]  # always delete circle and keep first circle (priority queue)
                                move_sequence.remove(circle)
                            else:
                                # if already parking edge in circle -> don't combine? and delete third circle?
                                all_circles[first_circle] = Mem1.combine_paths_over_pz(
                                    all_circles[first_circle], all_circles[circle]
                                )
                                del all_circles[circle]
                                move_sequence.remove(circle)
                            inner_if_broken = True
                            break
                    if inner_if_broken:
                        break

            print("all circles after pz combination: %s" % all_circles)

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
                print("rotate seq_idx", seq_idx)
                Mem1.rotate(free_circle_idxs[seq_idx])

            timestep += 1
            print("new timestep: %s" % timestep)

            # Save the current plot with a meaningful filename (plot widget)
            plot_filename = Path(run_folder) / f"plot_{timestep:03d}.png"
            print(sequence, "seq_element_counter", seq_element_counter)
            Mem1.graph_creator.plot_state(
                [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
                labels=["time step %s" % timestep, "next seq elem: %s" % sequence[0]],
                show_plot=show_plot,
                save_plot=save_plot,
                filename=[plot_filename if save_plot else None][0],
            )

            # if seq ion was at entry last timestep -> is now back in memory -> remove from sequence
            if seq_ion_was_at_entry is True:
                # seq_ion_was_at_entry is only True if seq ion is now not in entry anymore (for 2-qubit gates it has to wait in entry):
                if seq_element_counter == two_qubit_sequence[0] and next_seq_ion_in_exit < Mem1.time_2qubit_gate:
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
                    unique_sequence = []
                    for seq_elem in sequence:
                        if seq_elem not in unique_sequence:
                            unique_sequence.append(seq_elem)

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


# TODO repeating sequence elements

# 17.11.23
# two chains park at the end + 5 is exiting for last gate even though it could just stay

# e.g. time step 50 -> 0 goes out of entry towards top left (but only because free edge is in second row from bottom -> was correctly searched from exit, but lowest row is full on the left side -> should search at bottom be bfs?)
