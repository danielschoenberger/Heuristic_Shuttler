import contextlib
import math
import random
import time
from collections import Counter

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


archs = [[10, 10, 1, 1]]
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

        # sequence = [
        #     0,
        #     1,
        #     3,
        #     2,
        #     1,
        #     3,
        #     2,
        #     0,
        #     1,
        # ]  # TODO no immediately repeating seq elements are possible, e.g. [0, 1, 1, 2]
        sequence = [
            0,
            1,
            0,
            2,
            0,
            3,
            0,
            4,
            0,
            5,
            0,
            1,
            2,
            1,
            3,
            1,
            4,
            1,
            5,
            1,
            2,
            3,
            2,
            4,
            2,
            5,
            2,
            3,
            4,
            3,
            5,
            3,
            4,
            5,
            4,
            5,
        ]
        sequence = [
            0,
            1,
            0,
            2,
            0,
            3,
            0,
            4,
            0,
            5,
            0,
            6,
            0,
            7,
            0,
            8,
            0,
            9,
            0,
            10,
            0,
            11,
            0,
            12,
            0,
            13,
            0,
            14,
            0,
            15,
            0,
            16,
            0,
            17,
            0,
            18,
            0,
            19,
            0,
            1,
            2,
            1,
            3,
            1,
            4,
            1,
            5,
            1,
            6,
            1,
            7,
            1,
            8,
            1,
            9,
            1,
            10,
            1,
            11,
            1,
            12,
            1,
            13,
            1,
            14,
            1,
            15,
            1,
            16,
            1,
            17,
            1,
            18,
            1,
            19,
            1,
            2,
            3,
            2,
            4,
            2,
            5,
            2,
            6,
            2,
            7,
            2,
            8,
            2,
            9,
            2,
            10,
            2,
            11,
            2,
            12,
            2,
            13,
            2,
            14,
            2,
            15,
            2,
            16,
            2,
            17,
            2,
            18,
            2,
            19,
            2,
            3,
            4,
            3,
            5,
            3,
            6,
            3,
            7,
            3,
            8,
            3,
            9,
            3,
            10,
            3,
            11,
            3,
            12,
            3,
            13,
            3,
            14,
            3,
            15,
            3,
            16,
            3,
            17,
            3,
            18,
            3,
            19,
            3,
            4,
            5,
            4,
            6,
            4,
            7,
            4,
            8,
            4,
            9,
            4,
            10,
            4,
            11,
            4,
            12,
            4,
            13,
            4,
            14,
            4,
            15,
            4,
            16,
            4,
            17,
            4,
            18,
            4,
            19,
            4,
            5,
            6,
            5,
            7,
            5,
            8,
            5,
            9,
            5,
            10,
            5,
            11,
            5,
            12,
            5,
            13,
            5,
            14,
            5,
            15,
            5,
            16,
            5,
            17,
            5,
            18,
            5,
            19,
            5,
            6,
            7,
            6,
            8,
            6,
            9,
            6,
            10,
            6,
            11,
            6,
            12,
            6,
            13,
            6,
            14,
            6,
            15,
            6,
            16,
            6,
            17,
            6,
            18,
            6,
            19,
            6,
            7,
            8,
            7,
            9,
            7,
            10,
            7,
            11,
            7,
            12,
            7,
            13,
            7,
            14,
            7,
            15,
            7,
            16,
            7,
            17,
            7,
            18,
            7,
            19,
            7,
            8,
            9,
            8,
            10,
            8,
            11,
            8,
            12,
            8,
            13,
            8,
            14,
            8,
            15,
            8,
            16,
            8,
            17,
            8,
            18,
            8,
            19,
            8,
            9,
            10,
            9,
            11,
            9,
            12,
            9,
            13,
            9,
            14,
            9,
            15,
            9,
            16,
            9,
            17,
            9,
            18,
            9,
            19,
            9,
            10,
            11,
            10,
            12,
            10,
            13,
            10,
            14,
            10,
            15,
            10,
            16,
            10,
            17,
            10,
            18,
            10,
            19,
            10,
            11,
            12,
            11,
            13,
            11,
            14,
            11,
            15,
            11,
            16,
            11,
            17,
            11,
            18,
            11,
            19,
            11,
            12,
            13,
            12,
            14,
            12,
            15,
            12,
            16,
            12,
            17,
            12,
            18,
            12,
            19,
            12,
            13,
            14,
            13,
            15,
            13,
            16,
            13,
            17,
            13,
            18,
            13,
            19,
            13,
            14,
            15,
            14,
            16,
            14,
            17,
            14,
            18,
            14,
            19,
            14,
            15,
            16,
            15,
            17,
            15,
            18,
            15,
            19,
            15,
            16,
            17,
            16,
            18,
            16,
            19,
            16,
            17,
            18,
            17,
            19,
            17,
            18,
            19,
            18,
            19,
        ]

        # chain in entry is first sequence element at the start
        chain_in_entry = sequence[0]

        # 2-qubit sequence -> [2] means, that the 3rd ion in the sequence (ion with id 2) is a 2-qubit gate
        two_qubit_sequence = [
            1,
            3,
            5,
            7,
            9,
            11,
            13,
            15,
            17,
            19,
            21,
            23,
            25,
            27,
            29,
            31,
            33,
            35,
            37,
            40,
            42,
            44,
            46,
            48,
            50,
            52,
            54,
            56,
            58,
            60,
            62,
            64,
            66,
            68,
            70,
            72,
            74,
            77,
            79,
            81,
            83,
            85,
            87,
            89,
            91,
            93,
            95,
            97,
            99,
            101,
            103,
            105,
            107,
            109,
            112,
            114,
            116,
            118,
            120,
            122,
            124,
            126,
            128,
            130,
            132,
            134,
            136,
            138,
            140,
            142,
            145,
            147,
            149,
            151,
            153,
            155,
            157,
            159,
            161,
            163,
            165,
            167,
            169,
            171,
            173,
            176,
            178,
            180,
            182,
            184,
            186,
            188,
            190,
            192,
            194,
            196,
            198,
            200,
            202,
            205,
            207,
            209,
            211,
            213,
            215,
            217,
            219,
            221,
            223,
            225,
            227,
            229,
            232,
            234,
            236,
            238,
            240,
            242,
            244,
            246,
            248,
            250,
            252,
            254,
            257,
            259,
            261,
            263,
            265,
            267,
            269,
            271,
            273,
            275,
            277,
            280,
            282,
            284,
            286,
            288,
            290,
            292,
            294,
            296,
            298,
            301,
            303,
            305,
            307,
            309,
            311,
            313,
            315,
            317,
            320,
            322,
            324,
            326,
            328,
            330,
            332,
            334,
            337,
            339,
            341,
            343,
            345,
            347,
            349,
            352,
            354,
            356,
            358,
            360,
            362,
            365,
            367,
            369,
            371,
            373,
            376,
            378,
            380,
            382,
            385,
            387,
            389,
            392,
            394,
            397,
        ]  # [0, 2, 4, 6]
        # two_qubit_sequence = [1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 16, 17, 19]
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
        time_pz = 2

        Mem1 = MemoryZone(m, n, v, h, ion_chains, sequence, max_timesteps, time_pz=time_pz)
        timestep = 0
        print("time step: %s" % timestep)
        Mem1.graph_creator.plot_state(
            [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
            labels=timestep,
            show_plot=True,
        )

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
            for i, rotate_chain in enumerate(unique_sequence):
                # if 2-qubit gate -> check if second ion is in exit
                if (
                    seq_element_counter == two_qubit_sequence[0]
                    and rotate_chain == sequence[1]
                    and get_idx_from_idc(Mem1.idc_dict, Mem1.ion_chains[rotate_chain])
                    == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.exit_edge)
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
                        and next_seq_ion_in_exit < Mem1.time_pz
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

                    # for two-qubit gates: leave ion in entry until gate is processed (other ion was in exit long enough)

                    # entry
                    # if 2-qubit gate and in entry -> need to wait even if next edge is free
                    print("\nnext edge free")
                    print("seq_element_counter: %s" % seq_element_counter)
                    print("two_qubit_sequence: %s" % two_qubit_sequence)
                    print("next_seq_ion_in_exit: %s" % next_seq_ion_in_exit, "\n")
                    if (
                        get_idx_from_idc(Mem1.idc_dict, edge_idc)
                        == get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.entry_edge)
                        and seq_element_counter == two_qubit_sequence[0]
                    ):
                        if next_seq_ion_in_exit < Mem1.time_pz:
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
                    if rotate_chain not in sequence[:2]:
                        all_circles[rotate_chain] = [edge_idc, edge_idc]

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

                    # if 2-qubit gate and in entry -> need to wait
                    print("\npath through pz")
                    print("seq_element_counter: %s" % seq_element_counter)
                    print("two_qubit_sequence: %s" % two_qubit_sequence)
                    print("next_seq_ion_in_exit: %s" % next_seq_ion_in_exit, "\n")
                    if seq_element_counter == two_qubit_sequence[0]:
                        if next_seq_ion_in_exit < Mem1.time_pz:
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
                    if Mem1.graph_creator.processing_zone in edge:
                        if first_circle == -1:
                            first_circle = circle
                            inner_if_broken = True
                            break
                        else:  # Bricht hier, weil sollte warten und circles in pz nicht kombinieren, aber beide circles sind jetzt trotzdem in move sequence -> schiebt trotzdem durch
                            # wait move is special, nur in pz, könnte man daher auch speziell behandeln, zb wait move zählt nur erste node mit statt gar keine
                            if all(edge in all_circles[circle] for edge in all_circles[first_circle]):
                                # check if a "wait" move (edge, edge) is in first circle -> happens for 2-qubit gates in entry (ion in entry waits) -> don't combine circles
                                counter_first_circle = Counter(all_circles[first_circle])
                                if 2 in counter_first_circle.values():
                                    break
                                all_circles[first_circle] = all_circles[circle].copy()
                                del all_circles[circle]
                                move_sequence.remove(circle)
                            elif all(edge in all_circles[first_circle] for edge in all_circles[circle]):
                                counter_first_circle = Counter(all_circles[first_circle])
                                del all_circles[circle]  # always delete circle and keep first circle (priority queue)
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
                print("seq_idx", seq_idx)
                Mem1.rotate(free_circle_idxs[seq_idx])

            timestep += 1
            print("new timestep: %s" % timestep)
            Mem1.graph_creator.plot_state(
                [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
                str(timestep),
                show_plot=True,
            )

            # if seq ion was at entry last timestep -> is now back in memory -> remove from sequence
            if seq_ion_was_at_entry is True:
                # seq_ion_was_at_entry is only True if seq ion is now not in entry anymore (for 2-qubit gates it has to wait in entry):
                if seq_element_counter == two_qubit_sequence[0] and next_seq_ion_in_exit < Mem1.time_pz:
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

# TODO Schalter
# TODO cut sequence? Unique sequence reasonable?
# TODO repeating sequence elements

# e.g. time step 50 -> 0 goes out of entry towards top left (but only because free edge is in second row from bottom -> was correctly searched from exit, but lowest row is full on the left side -> should search at bottom be bfs?)
