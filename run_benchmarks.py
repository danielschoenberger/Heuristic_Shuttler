import contextlib
import math
import random
import time
from datetime import datetime
from pathlib import Path

import networkx as nx
import numpy as np
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dagdependency

from compilation_new import manual_copy_dag, remove_node, update_sequence
from Cycles import GraphCreator, MemoryZone, get_idx_from_idc


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


def preprocess(memorygrid, sequence):
    need_rotate = [False] * len(sequence)
    while sum(need_rotate) < len(sequence):
        for i, rotate_chain in enumerate(sequence):
            edge_idc = memorygrid.ion_chains[rotate_chain]
            next_edge = memorygrid.find_next_edge(edge_idc)
            state_edges_idx = memorygrid.get_state_idxs()

            if (
                memorygrid.have_common_junction_node(edge_idc, next_edge) is False
                and get_idx_from_idc(memorygrid.idc_dict, next_edge) not in state_edges_idx
            ):
                memorygrid.ion_chains[rotate_chain] = next_edge
            else:
                need_rotate[i] = True
    return memorygrid


def create_move_list(memorygrid, sequence, max_length=5):
    """
    max_length: max length of move_list (if sequence is longer than max_length, only first max_length elements are considered)
    """
    # unique sequence is sequence without repeating elements (for move_list and 2-qubit gates)
    unique_sequence = []
    for seq_elem in sequence:
        if seq_elem not in unique_sequence:
            unique_sequence.append(seq_elem)
            if len(unique_sequence) == max_length:
                break

    path_length_sequence = {}
    move_list = []
    for i, rotate_chain in enumerate(unique_sequence):
        edge_idc = memorygrid.ion_chains[rotate_chain]
        # TODO shortest path here maybe not optimal?
        path_to_go = nx.shortest_path(
            memorygrid.graph,
            edge_idc[0],
            memorygrid.graph_creator.processing_zone,
            lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "entry") * 1e8 + 1,
        )
        path_length_sequence[rotate_chain] = len(path_to_go)

        if i == 0 or sum(
            np.array([path_length_sequence[rotate_chain]] * len(move_list))
            > np.array([path_length_sequence[chain] for chain in move_list])
        ) == len(move_list):
            move_list.append(rotate_chain)

        # if ion in entry edge -> move back to memory zone (highest priority)
        if get_idx_from_idc(memorygrid.idc_dict, edge_idc) == get_idx_from_idc(
            memorygrid.idc_dict, memorygrid.graph_creator.entry_edge
        ):
            path_length_sequence[rotate_chain] = 0
            with contextlib.suppress(Exception):
                move_list.remove(rotate_chain)
            move_list = [rotate_chain, *move_list]

    # get chains in all entry edges and place in front
    # chain in entry must move out
    chains_in_entry_connections = []
    for ion, chain_idx in enumerate(memorygrid.get_state_idxs()):
        if chain_idx in memorygrid.graph_creator.path_from_pz_idxs:
            if chain_idx == memorygrid.graph_creator.entry_edge:
                # place chain in entry at the end of move_list -> so later looping over list leads to chain in entry being first
                chains_in_entry_connections.append(ion)
            else:
                chains_in_entry_connections.insert(0, ion)

    if len(chains_in_entry_connections) > 0:
        for ion in chains_in_entry_connections:
            with contextlib.suppress(Exception):
                move_list.remove(ion)
            move_list = [ion, *move_list]

    return move_list


archs = [[3, 4, 2, 2], [5, 5, 1, 1], [6, 6, 1, 1]]  # , [20, 20, 1, 1], [5, 5, 10, 10]]#[5, 5, 1, 1],
seeds = [1]  # , 2, 3, 4, 5, 6, 7, 8, 9, 10]
perc = 0.5
results = {}
cpu_time_results = {}
start_time_all = time.time()
show_plot = False
save_plot = False  # not show_plot
use_gate_execution = True

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
        # max_chains_in_pz = 4 # TODO needed?
        max_chains_in_parking = 3

        # two_qubit_sequence = []
        # i = 0
        # while i < len(seq):
        #     if len(seq[i]) == 2:
        #         two_qubit_sequence.append(i)
        #     i += 1
        # two_qubit_sequence.append(
        #     -1
        # )  # -1 at the end, so two-qubit sequence is not empty after all 2-qubit gates have been processed
        # assert (
        #     two_qubit_sequence[-2] != len(flat_seq) - 1
        # ), "2-qubit sequence is not valid (last element can not be 2-qubit gate)"

        # seq_element_counter = 0

        print(f"arch: {arch}, seed: {seed}, registers: {number_of_registers}\n")
        max_timesteps = 50000
        # time_2qubit_gate only works for 2-qubit gates -> need also for 1-qubit gates
        time_2qubit_gate = 3
        time_1qubit_gate = 1
        time_in_pz_counter = 0

        Mem1 = MemoryZone(
            m,
            n,
            v,
            h,
            ion_chains,
            max_timesteps,
            max_chains_in_parking,
            time_2qubit_gate=time_2qubit_gate,
            time_1qubit_gate=time_1qubit_gate,
        )
        timestep = 0
        print("time step: %s" % timestep)

        seq_ion_was_at_entry = 0
        next_seq_ion_in_exit = 0

        seq_ion_was_at_entry = False
        gate_execution_finished = True

        # generate sequence and two-qubit sequence
        N = number_of_registers
        filename = "QASM_files/qft_%squbits.qasm" % N
        # write_qft_to_file(N, filename)
        # print(f"QFT for {N} qubits written to {filename}\n")
        # seq = parse_qasm(filename)

        qc = QuantumCircuit.from_qasm_file(filename)
        dag_dep = circuit_to_dagdependency(qc)

        Mem1.update_distance_map()
        gate_ids, next_node = update_sequence(dag_dep, Mem1.distance_map)
        seq = [tuple(gate) for gate in gate_ids]
        flat_seq = [item for sublist in seq for item in sublist]
        next_gate_is_two_qubit_gate = len(seq[0]) == 2

        # Save the current plot (plot widget)
        plot_filename = Path(run_folder) / f"plot_{1:03d}.png"

        Mem1.graph_creator.plot_state(
            [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
            labels=[
                "time step %s" % timestep,
                f"next seq elem: {seq[0]}",
            ],
            show_plot=show_plot,
            save_plot=save_plot,
            filename=[plot_filename if save_plot else None][0],
        )

        # timestep = 1
        while timestep < max_timesteps:
            rotate_entry = False

            # update state_idxs
            Mem1.get_state_idxs()

            ######### PREPROCESS #########
            Mem1 = preprocess(Mem1, flat_seq)

            Mem1.update_distance_map()

            # TODO update here or only update sequence after gate execution?
            # new_sequence = update_sequence(dag_dep, Mem1.distance_map)
            # print(new_sequence, '\n', flat_seq)

            ######### CREATE MOVE SEQUENCE #########
            move_list = create_move_list(Mem1, flat_seq)

            ######### CREATE CIRCLES #########
            ### create circles for all chains in move_list (dictionary with chain as key and circle_idcs as value)
            chain_to_park = Mem1.find_chain_in_edge(Mem1.graph_creator.path_to_pz[-1])
            if Mem1.count_chains_in_parking() < Mem1.max_num_parking or gate_execution_finished:
                parking_open = True
            else:
                parking_open = False

            all_circles = {}
            stop_exit_edges = []
            # need to find next_edges before for bfs search of "out of entry move"
            next_edges = {}
            for rotate_chain in move_list:
                edge_idc = Mem1.ion_chains[rotate_chain]
                # if chain is needed again (is present in rest of sequence) -> move (only out of entry) towards exit instead of top left
                towards = "exit" if rotate_chain in flat_seq[1:] else (0, 0)
                next_edges[rotate_chain] = Mem1.find_next_edge(edge_idc, towards=towards)

            for rotate_chain in move_list:
                edge_idc = Mem1.ion_chains[rotate_chain]
                next_edge = next_edges[rotate_chain]

                # make edge_idc and next_edge consistent
                edge_idc, next_edge = Mem1.find_ordered_edges(edge_idc, next_edge)

                # moves in pz
                if get_idx_from_idc(Mem1.idc_dict, next_edge) in [
                    *Mem1.graph_creator.path_to_pz_idxs,
                    get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.parking_edge),
                ]:
                    all_circles[rotate_chain] = [edge_idc, next_edge]
                    # block moves to pz if parking is full
                    if (
                        get_idx_from_idc(Mem1.idc_dict, next_edge)
                        in [
                            *Mem1.graph_creator.path_to_pz_idxs,
                            get_idx_from_idc(Mem1.idc_dict, Mem1.graph_creator.parking_edge),
                        ]
                        and parking_open is False
                    ) and (get_idx_from_idc(Mem1.idc_dict, next_edge) in stop_exit_edges or stop_exit_edges == []):
                        all_circles[rotate_chain] = [edge_idc, edge_idc]
                        # needed later for blocking moves to parking
                        stop_exit_edges.append(get_idx_from_idc(Mem1.idc_dict, edge_idc))
                # moves without circle
                elif not Mem1.check_if_edge_is_filled(next_edge):
                    all_circles[rotate_chain] = [edge_idc, next_edge]
                # moves with circle
                else:
                    # create circle (deleted in create_outer_circle: in parking circle is a "stop move")
                    all_circles[rotate_chain] = Mem1.create_outer_circle(edge_idc, next_edge, next_edges.values())

            # move chain out of parking edge if needed
            chains_in_parking = Mem1.find_chains_in_parking()
            # if pz full and no chain is moving out (not in state_idxs entry edge) but chain is moving in
            if (
                Mem1.count_chains_in_parking() >= Mem1.max_num_parking
                and gate_execution_finished
                and chain_to_park is not None
            ):
                # find least important chain in parking edge
                chain_to_move_out_of_pz = Mem1.find_least_import_chain_in_parking(
                    flat_seq, [*chains_in_parking, chain_to_park]
                )
                if chain_to_move_out_of_pz != chain_to_park:
                    # move it to entry
                    rotate_entry = True
                    # change its path/circle to a stop move
                    all_circles[chain_to_move_out_of_pz] = [
                        Mem1.graph_creator.path_from_pz[0],
                        Mem1.graph_creator.path_from_pz[0],
                    ]

            ######### FIND CIRCLES THAT CAN MOVE #########
            # find circles that can move while first seq ion is moving
            nonfree_circles, free_circle_combs = Mem1.find_nonfree_and_free_circle_idxs(all_circles)
            free_circle_seq_idxs = [move_list[0]]
            for seq_circ in move_list[1:]:
                nonfree = False
                for mov_circ in free_circle_seq_idxs:
                    if (seq_circ, mov_circ) in nonfree_circles or (mov_circ, seq_circ) in nonfree_circles:
                        nonfree = True
                        break
                if nonfree is False:
                    free_circle_seq_idxs.append(seq_circ)

            ######### ROTATE CIRCLES #########
            # need circles given in idxs for rotate function
            free_circle_idxs = {}
            for seq_idx in free_circle_seq_idxs:
                free_circle_idxs[seq_idx] = [
                    get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in all_circles[seq_idx]
                ]
                # rotate chains
                new_state_dict = Mem1.rotate(free_circle_idxs[seq_idx])
                if rotate_entry:
                    Mem1.ion_chains[chain_to_move_out_of_pz] = Mem1.graph_creator.path_from_pz[0]

            ######### PLOT #########
            # Save the current plot (plot widget)
            plot_filename = Path(run_folder) / f"plot_{timestep:03d}.png"

            ######### UPDATE SEQUENCE / PROCESS GATE #########
            gate = seq[0]
            chains_in_parking = Mem1.find_chains_in_parking()
            if sum((gate_element in chains_in_parking) for gate_element in gate) == len(gate):
                # new TODO use gate_execution_finished or not?
                if use_gate_execution is True:
                    gate_execution_finished = False

                # track start of gate execution
                start_execution = time_in_pz_counter == 0

                time_in_pz_counter += 1
                Mem1.graph_creator.plot_state(
                    [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
                    labels=[
                        "time step %s" % timestep,
                        f"seq elem {seq[0]} execution",
                    ],
                    show_plot=show_plot,
                    save_plot=save_plot,
                    filename=[plot_filename if save_plot else None][0],
                )

                print(f"\ntime step: {timestep}, gate {seq[0]} is executed,")
                time_gate = time_2qubit_gate if next_gate_is_two_qubit_gate is True else time_1qubit_gate

                if time_in_pz_counter == time_gate:
                    ######### END IF SEQUENCE IS FINISHED #########
                    if len(seq) == 1:
                        print("\nFull Sequence executed")
                        cpu_time_arr.append(time.time() - start_time)
                        timestep_arr.append(timestep)
                        break

                    for _ in gate:
                        flat_seq.pop(0)
                    # if seq_element_counter == two_qubit_sequence[0]:
                    #     two_qubit_sequence.pop(0)
                    # seq_element_counter += 1
                    time_in_pz_counter = 0
                    gate_execution_finished = True

                    # update dag
                    remove_node(dag_dep, next_node)
                    dag_dep = manual_copy_dag(dag_dep)

                    gate_ids, next_node = update_sequence(dag_dep, Mem1.distance_map)
                    seq = [tuple(gate) for gate in gate_ids]
                    flat_seq = [item for sublist in seq for item in sublist]
                    next_gate_is_two_qubit_gate = len(seq[0]) == 2

            else:
                Mem1.graph_creator.plot_state(
                    [get_idx_from_idc(Mem1.idc_dict, edge_idc) for edge_idc in Mem1.ion_chains.values()],
                    labels=["time step %s" % timestep, f"next seq elem: {seq[0]}"],
                    show_plot=show_plot,
                    save_plot=save_plot,
                    filename=[plot_filename if save_plot else None][0],
                )

            ######### SETUP NEW TIME STEP #########
            timestep += 1

    timestep_mean = np.mean(timestep_arr)
    cpu_time_mean = np.mean(cpu_time_arr)
    results[j] = timestep_mean
    cpu_time_results[j] = cpu_time_mean

    # Create a Path object for the file
    file_path = Path("test_benchmark_results.txt")

    # Open the file using the Path object
    with file_path.open("a") as file:
        line = f"& {arch[0]} {arch[1]} {arch[2]} {arch[3]} & {number_of_registers}/{n_of_traps} ({int(100*number_of_registers/n_of_traps)}%) & {timestep_mean} & {cpu_time_mean} {'s'} & {timestep_mean} \\\\"
        file.write(line + "\n")

print("\n archs: \n", dict(enumerate(archs)))
print("results: \n", results)
print("cpu time results: \n", cpu_time_results)
print("time all: \n", time.time() - start_time_all)


# TODO
# - move out of parking -> double-check if still works correctly [done]
# - move from parking to entry [done]
# - add time for 1-qubit gates [done]
# - check if correct chain moves out of parking [x]
# - check if chain moves to correct edge after moving out of entry [x]
# - rewrite find_nonfree_and_free_circle_idxs(self, circles_dict) -> some things not needed anymore [x]
# - clean-up [x]
